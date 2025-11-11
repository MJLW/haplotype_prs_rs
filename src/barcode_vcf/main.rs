use std::{
    collections::HashMap,
    error::Error,
    fs::File,
    io::{self, BufReader, BufWriter},
    iter::{Peekable, zip},
    path::Path,
};

use clap::Parser;
use serde::Deserialize;

use noodles::{
    bgzf,
    core::Position,
    cram::{self, crai},
    fasta::{self, fai},
    sam,
    vcf::{
        self, Record,
        variant::{
            record::{
                self, AlternateBases,
                samples::{
                    Series,
                    keys::key,
                    series::value::{Genotype, genotype::Phasing},
                },
            },
            record_buf::samples::sample::{self, value::genotype::Allele},
        },
    },
};

/**
* Argument serialization
*/
#[derive(Parser, Debug)]
#[command(name = "Barcode VCF", version = env!("CARGO_PKG_VERSION"), about = "Sample barcoding", long_about = "Barcodes a sample using SNP file and a bgzip-compressed VCF.")]
struct Args {
    #[arg(short, long)]
    snp_in: String,

    #[arg(short, long)]
    vcf_in: String,

    #[arg(short, long)]
    cram_in: String,

    #[arg(short, long)]
    fasta_in: String,

    #[arg(short, long)]
    out: String,
}

/**
* File serialization
*/

#[derive(Deserialize)]
struct SNP {
    #[serde(rename = "CHROM")]
    chr: String,

    #[serde(rename = "POS")]
    pos: usize,

    #[serde(rename = "REF")]
    var_ref: char,

    #[serde(rename = "ALT")]
    var_alt: char,

    #[serde(rename = "AF")]
    af: f32,

    #[serde(rename = "BETA")]
    beta: f32,
}

/**
* Genotype extensions
*/
trait GenotypeToString {
    fn to_string(&self) -> String;
}

impl GenotypeToString for sample::value::Genotype {
    fn to_string(&self) -> String {
        let is_phased = self
            .iter()
            .map(|r| r.unwrap())
            .any(|(_, phasing)| phasing == Phasing::Phased);

        let separator = if is_phased { "|" } else { "/" };

        let string_alleles: Vec<String> = self
            .iter()
            .map(|r| r.unwrap())
            .map(|(allele, _)| match allele {
                Some(al) => al.to_string(),
                None => ".".to_string(),
            })
            .collect();

        return string_alleles.join(separator);
    }
}

trait DefaultMissing {
    fn missing() -> Self;
}

impl DefaultMissing for sample::value::Genotype {
    fn missing() -> Self {
        return sample::value::Genotype::from_iter([
            Allele::new(Some(0), Phasing::Unphased),
            Allele::new(Some(0), Phasing::Unphased),
        ]);
    }
}

/**
* SNPs
*/

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

/**
* VCFs
*/

fn open_vcf(path: &str) -> Result<vcf::io::Reader<bgzf::io::Reader<File>>, Box<dyn Error>> {
    let reader = File::open(path)
        .map(bgzf::io::Reader::new)
        .map(vcf::io::Reader::new)?;

    Ok(reader)
}

fn setup_contig_hashmap(header: &vcf::Header) -> HashMap<&str, usize> {
    return HashMap::from_iter(
        header
            .contigs()
            .keys()
            .into_iter()
            .enumerate()
            .map(|(i, key)| (key.as_str(), i)),
    );
}

fn get_contig_id(
    contig: &str,
    contig_hashmap: &HashMap<&str, usize>,
) -> Result<usize, Box<dyn Error>> {
    return match contig_hashmap.get(contig) {
        Some(id) => Ok(*id),
        None => Err(format!("Contig '{}' does not exist in VCF header.", contig).into()),
    };
}

// fn get_af(
//     record: &vcf::Record,
//     header: &vcf::Header,
//     tag: &str,
//     allele_index: usize,
// ) -> Result<f32, Box<dyn Error>> {
//     let info = record.info();
//
//     let info_af = info
//         .get(header, tag)
//         .expect(format!("Could not find tag INFO/{}", tag).as_str())
//         .expect("")
//         .expect("");
//
//     match info_af {
//         record::info::field::Value::Float(af_value) => return Ok(af_value),
//         record::info::field::Value::Array(af_array) => {
//             match af_array {
//                 record::info::field::value::Array::Float(af_array_float) => {
//                     let target_af = af_array_float.iter().enumerate().find(|(i, _)| (i+1) == allele_index);
//
//                     if target_af.is_none() { return Err(format!("Could not find target allele AF in INFO/{}", tag).into()); }
//
//                     match target_af.unwrap().1? {
//                         Some(af) => return Ok(af)
//                         ,
//                         None => return Err("".into())
//                     }
//                 },
//                 _ => return Err(format!(
//                     "INFO/{} is not neither a float nor an Array of Float's. Cannot parse it as an AF tag.",
//                     tag
//                 ).into()),
//             };
//         }
//         _ => Err(format!(
//             "INFO/{} is neither a Float nor an Array of Float's. Cannot parse it as an AF tag.",
//             tag
//         )
//         .into()),
//     }
// }

fn get_gt(
    record: &vcf::Record,
    header: &vcf::Header,
) -> Result<sample::value::Genotype, Box<dyn Error>> {
    let samples = record.samples();
    let gt_series = samples.select(key::GENOTYPE).expect("Missing GT field");

    let gt_value = gt_series
        .get(header, 0)
        .expect("No GT field.") // Result
        .expect("No GT value.")?; // Option

    // record::samples::series::Value::Genotype(genotyp)
    let record::samples::series::Value::Genotype(genotype) = gt_value else {
        panic!("Failed to cast genotype series element, underlying data is invalid.");
    };

    Ok(genotype.as_ref().try_into().unwrap())
}

fn filter_gt_to_allele(
    gt: sample::value::Genotype,
    allele_index: usize,
) -> sample::value::Genotype {
    return sample::value::Genotype::from_iter(gt.as_ref().iter().map(|allele: &Allele| {
        match allele.position() {
            Some(position) => {
                if position == 0 || position == allele_index {
                    allele.clone()
                } else {
                    Allele::new(None, allele.phasing())
                }
            }
            None => allele.clone(),
        }
    }));
}

fn find_matching_record<I>(
    snp: &SNP,
    contig_hashmap: &HashMap<&str, usize>,
    records_iter: &mut Peekable<I>,
    header: &vcf::Header,
    // info_af_tag: &Option<String>,
) -> Result<sample::value::Genotype, Box<dyn Error>>
where
    I: Iterator<Item = io::Result<Record>>,
{
    let snp_contig_id = get_contig_id(&snp.chr, contig_hashmap)?;
    let snp_pos = Position::try_from(snp.pos).unwrap();

    while let Some(peeked) = records_iter.peek() {
        let peeked_record = peeked.as_ref().unwrap();
        let record_contig_id =
            get_contig_id(peeked_record.reference_sequence_name(), &contig_hashmap)?;
        let record_pos = peeked_record.variant_start().unwrap().unwrap();

        // If peeked record is still upstream of the SNP, consume and attempt next
        if record_contig_id < snp_contig_id
            || (record_contig_id == snp_contig_id && record_pos < snp_pos)
        {
            records_iter.next();
            continue;
        }

        // If peeked record is downstream of the SNP, the SNP must not be present
        if record_contig_id > snp_contig_id
            || (record_contig_id == snp_contig_id && record_pos > snp_pos)
        {
            return Ok(sample::value::Genotype::missing());
        }

        // SNP and VCF record must have same coordinates, now check if alleles match
        let alts = peeked_record.alternate_bases();

        let matched_alt = alts.iter().enumerate().find(|(_, alt)| {
            peeked_record.reference_bases().len() == 1
                && alt.as_ref().unwrap().len() == 1
                && alt.as_ref().unwrap().chars().last().unwrap() == snp.var_alt
        });

        if matched_alt.is_none() {
            // If multiallelics are split, the SNP might be on the next line instead, so just consume and skip.
            // If they are not split, the next peek will be downstream and will return missing anyway
            records_iter.next();
            continue;
        }

        // Record has been found, huzzaah!
        let alt_allele = matched_alt.unwrap().0 + 1;
        let record = records_iter.next().unwrap()?;

        // // Fetch and match AF
        // if info_af_tag.is_some() {
        //     let tag = info_af_tag.as_ref().unwrap();
        //     let variant_af = get_af(&record, header, tag, alt_allele)?;
        //
        //     // TODO: Is there any point when this will not meet HWE for SNP AF anyway?
        //     const THRESHOLD: f32 = 0.3;
        //     if f32::abs(snp.af - variant_af) > THRESHOLD {
        //         return Ok(sample::value::Genotype::missing());
        //     }
        // }

        // Fetch and match GT
        let gt = get_gt(&record, header)?;
        let snp_gt = filter_gt_to_allele(gt, alt_allele);

        return Ok(snp_gt);
    }

    Ok(sample::value::Genotype::missing())
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();

    let snps: Vec<SNP> = read_snps(&args.snp_in)?;
    let mut vcf_reader = open_vcf(&args.vcf_in)?;
    let vcf_header: vcf::Header = vcf_reader.read_header()?;
    let contig_hashmap = setup_contig_hashmap(&vcf_header);

    // Define peekable iterator to look at downstream variants without consuming them.
    // They might still be needed to match other SNPs.
    let mut records_iter = vcf_reader.records().peekable();
    // let af_tag = &args.info_af_tag;

    let gts: Vec<_> = snps
        .iter()
        .map(|snp| {
            find_matching_record(
                &snp,
                &contig_hashmap,
                &mut records_iter,
                &vcf_header,
                // af_tag,
            )
        })
        .collect();

    // Shouldn't ever happen, but very bad things could happen in downstream analysis if this ever evaluated to true
    if gts.len() != snps.len() {
        return Err(format!(
            "Mismatch between number of SNPs ({}) and retrieved calls ({}). Report this.",
            gts.len(),
            snps.len()
        )
        .into());
    }

    // TODO: Add option for checking coverage in cram file
    // let mut cram_reader = cram::io::Reader::new(BufReader::new(File::open(&args.cram_in)?));
    // let cram_header: sam::Header = cram_reader.read_header()?;
    // let cram_index = crai::fs::read(format!("{}.crai", &args.cram_in))?;
    //
    // // Preload fasta for CRAM access
    // let fasta_reader =
    //     fasta::io::indexed_reader::Builder::default().build_from_path(&args.fasta_in)?;
    // let fasta_adapter = fasta::repository::adapters::IndexedReader::new(fasta_reader);

    // Write output to file
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_writer(BufWriter::new(File::create(Path::new(&args.out))?));

    // Write out header
    writer.write_record(&["CHROM", "POS", "REF", "ALT", "AF", "BETA", "GT"])?;

    for (snp, gt) in zip(snps, gts) {
        let call = gt.unwrap().to_string();
        writer.serialize((
            snp.chr,
            snp.pos,
            snp.var_ref,
            snp.var_alt,
            snp.af,
            snp.beta,
            call,
        ))?;
    }

    Ok(())
}
