use glob::glob;
use noodles::vcf::variant::{
    record::samples::{keys::key::GENOTYPE, series::value::genotype::Phasing},
    record_buf::samples::sample::value::Genotype,
};

use clap::Parser;
use serde::Deserialize;
use std::{
    error::Error,
    fs::File,
    io::{BufReader, BufWriter},
    iter::zip,
    path::{Path, PathBuf},
    str::FromStr,
};

/**
* Argument serialization
*/
#[derive(Parser, Debug)]
#[command(name = "Barcode VCF", version = env!("CARGO_PKG_VERSION"), about = "Sample barcoding", long_about = "Barcodes a sample using SNP file and a bgzip-compressed VCF.")]
struct Args {
    #[arg(short, long)]
    dir_in: String,

    #[arg(short, long)]
    ext_in: String,

    #[arg(short, long)]
    out: String,
}

/**
* File serialization
*/

#[derive(Deserialize)]
struct SNP {
    // #[serde(rename = "CHROM")]
    // chr: String,
    //
    // #[serde(rename = "POS")]
    // pos: usize,
    //
    // #[serde(rename = "REF")]
    // var_ref: char,
    //
    // #[serde(rename = "ALT")]
    // var_alt: char,
    //
    // #[serde(rename = "AF")]
    // af: f32,
    #[serde(rename = "BETA")]
    beta: f32,
}

fn read_snps(path: &str) -> Result<Vec<SNP>, Box<dyn Error>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_reader(BufReader::new(File::open(Path::new(path))?));

    let snps: Vec<SNP> = reader
        .deserialize()
        .into_iter()
        .map(|r| r.unwrap())
        .collect();

    Ok(snps)
}

fn read_calls(path: &str, gt_name: &str) -> Result<Vec<Genotype>, Box<dyn Error>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(path)?;

    let headers = reader.byte_headers()?;
    let gt_index = headers
        .iter()
        .position(|name| name == gt_name.as_bytes())
        .expect(format!("Could not find column {} in header", gt_name).as_str());

    let calls: Vec<Genotype> = reader
        .records()
        .map(|r| {
            let record = r.unwrap();

            if let Some(str) = record.get(gt_index) {
                return Genotype::from_str(str).unwrap_or_default();
            }

            return Genotype::default();
        })
        .collect();

    Ok(calls)
}

fn get_barcode_paths(dir_path: &str, extension: &str) -> Result<Vec<PathBuf>, Box<dyn Error>> {
    let dir = Path::new(dir_path);

    if !Path::is_dir(dir) {
        return Err(format!("Path is not a directory: {}", dir_path).into());
    }

    let paths: Vec<_> = glob(format!("{}/*{}", dir_path, extension).as_str())
        .expect("Glob pattern failure")
        .filter_map(|entry| match entry {
            Ok(path) => Some(path),
            Err(_) => None,
        })
        .collect();

    Ok(paths)
}

fn sum_beta_per_haplotype(
    gts: &Vec<Genotype>,
    beta: &Vec<f32>,
) -> (f32, f32, usize, usize, usize, usize, usize) {
    let mut paternal_sum = 0.0;
    let mut maternal_sum = 0.0;

    let mut unphased_het_snps = 0;
    let mut paternal_het_snps = 0;
    let mut maternal_het_snps = 0;
    let mut hom_snps_alt = 0;
    let mut hom_snps_ref = 0;
    for (gt, beta) in zip(gts, beta) {
        let alleles = gt.as_ref();
        let paternal = alleles.first();
        let maternal = alleles.last();

        if paternal.is_none() || maternal.is_none() {
            continue;
        }

        let paternal_allele = paternal.unwrap();
        let maternal_allele = maternal.unwrap();

        if paternal_allele.phasing() != maternal_allele.phasing() {
            continue;
        }

        if paternal_allele.position().is_some()
            && paternal_allele.position().unwrap() == 0
            && maternal_allele.position().is_some()
            && maternal_allele.position().unwrap() == 0
        {
            hom_snps_ref += 1;
            continue;
        }

        if paternal_allele.phasing() == Phasing::Phased {
            if (paternal_allele.position().is_some() && paternal_allele.position().unwrap() > 0)
                && (maternal_allele.position().is_none()
                    || maternal_allele.position().unwrap() == 0)
            {
                paternal_sum += beta;
                paternal_het_snps += 1;
                continue;
            } else if (paternal_allele.position().is_none()
                || paternal_allele.position().unwrap() == 0)
                && (maternal_allele.position().is_some() && maternal_allele.position().unwrap() > 0)
            {
                maternal_sum += beta;
                maternal_het_snps += 1;
                continue;
            }
        } else {
            if (paternal_allele.position().is_some() && paternal_allele.position().unwrap() > 0)
                && (maternal_allele.position().is_some() && maternal_allele.position().unwrap() > 0)
            {
                paternal_sum += beta;
                maternal_sum += beta;
                hom_snps_alt += 1;
            } else if (paternal_allele.position().is_some()
                && paternal_allele.position().unwrap() > 0)
                || (maternal_allele.position().is_some() && maternal_allele.position().unwrap() > 0)
            {
                paternal_sum += beta;
                maternal_sum += beta;
                unphased_het_snps += 1;
            }
        }
    }

    return (
        paternal_sum,
        maternal_sum,
        unphased_het_snps,
        paternal_het_snps,
        maternal_het_snps,
        hom_snps_alt,
        hom_snps_ref,
    );
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();

    let paths: Vec<PathBuf> = get_barcode_paths(&args.dir_in, &args.ext_in)?;

    let gts_matrix: Vec<Vec<Genotype>> = paths
        .iter()
        .map(|path| path.clone().into_os_string().into_string().unwrap())
        .filter_map(|path| match read_calls(path.as_str(), GENOTYPE) {
            Ok(gt) => Some(gt),
            Err(_) => None,
        })
        .collect();

    let snps = read_snps(
        paths
            .first()
            .unwrap()
            .clone()
            .into_os_string()
            .into_string()
            .unwrap()
            .as_str(),
    )?;

    let betas: Vec<f32> = snps.into_iter().map(|snp| snp.beta).collect();

    let allele_betas: Vec<(f32, f32, usize, usize, usize, usize, usize)> = gts_matrix
        .iter()
        .map(|gts| sum_beta_per_haplotype(gts, &betas))
        .collect();

    let file_names: Vec<&str> = paths
        .iter()
        .map(|path| path.file_name().unwrap().to_str().unwrap())
        .collect();

    // Init writer
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_writer(BufWriter::new(File::create(Path::new(&args.out))?));

    writer.write_record(&[
        "SAMPLE",
        "PATERNAL",
        "MATERNAL",
        "N_HET_UNPHASED",
        "N_HET_PATERNAL",
        "N_HET_MATERNAL",
        "N_HOM_ALT",
        "N_HOM_REF",
    ])?;
    for (sample, beta_sums) in file_names.into_iter().zip(allele_betas) {
        writer.write_record(&[
            sample,
            beta_sums.0.to_string().as_str(),
            beta_sums.1.to_string().as_str(),
            beta_sums.2.to_string().as_str(),
            beta_sums.3.to_string().as_str(),
            beta_sums.4.to_string().as_str(),
            beta_sums.5.to_string().as_str(),
            beta_sums.6.to_string().as_str(),
        ])?;
    }

    Ok(())
}
