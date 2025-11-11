use std::{
    char,
    error::Error,
    fs::File,
    io::{self, BufReader, BufWriter},
    path::Path,
};

use noodles::{
    core::Position,
    fasta::{self, Repository, repository::adapters::IndexedReader},
};

use clap::Parser;

use serde::{Deserialize, Serialize};

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    #[arg(short, long)]
    fasta: String,

    #[arg(short, long)]
    snp_in: String,

    #[arg(short, long)]
    out: String,
}

#[derive(Deserialize, Serialize)]
struct SNP {
    #[serde(rename = "CHROM")]
    chr: String,

    #[serde(rename = "POS")]
    pos: usize,

    #[serde(rename = "EFFECT")]
    effect: char,

    #[serde(rename = "OTHER")]
    other: char,

    #[serde(rename = "AF")]
    af: f32,

    #[serde(rename = "BETA")]
    beta: f32,
}

enum SNPOrder {
    Reference,
    Alternate,
}

fn open_fasta(path: &str) -> Result<Repository, Box<dyn Error>> {
    let reader = fasta::io::indexed_reader::Builder::default().build_from_path(path)?;
    let adapter = IndexedReader::new(reader);
    let repository = fasta::Repository::new(adapter);

    Ok(repository)
}

fn get_position_from_reference(
    reference: &Repository,
    region: &str,
    pos: usize,
) -> Result<u8, Box<dyn Error>> {
    let sequence = reference
        .get(region.as_bytes())
        .transpose()?
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Invalid region: {}", region),
            )
        })?;

    let start = Position::try_from(pos)?;
    let base = sequence.get(start..=start).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::NotFound,
            format!("Position {} does not exist in region '{}'", pos, region),
        )
    })?;

    Ok(base[0])
}

fn snp_as_complementary(snp: SNP) -> SNP {
    return SNP {
        chr: snp.chr,
        pos: snp.pos,
        effect: snp.other,
        other: snp.effect,
        af: 1.0 - snp.af,
        beta: -snp.beta,
    };
}

fn check_snp_order(snp: &SNP, reference: &Repository) -> Result<SNPOrder, Box<dyn Error>> {
    let ref_base = char::from(get_position_from_reference(reference, &snp.chr, snp.pos)?);

    // Effect should be ALT, Other should be REF
    if ref_base == snp.effect {
        return Ok(SNPOrder::Alternate); // Effect == REF
    } else if ref_base == snp.other {
        return Ok(SNPOrder::Reference); // Other == REF
    } else {
        return Err(format!(
            "SNP at {}:{} does not match reference base '{}'. Found effect='{}', other='{}'.",
            snp.chr, snp.pos, ref_base, snp.effect, snp.other
        )
        .into());
    }
}

fn read_snps(path: &str, reference: &Repository) -> Result<Vec<SNP>, Box<dyn Error>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_reader(BufReader::new(File::open(Path::new(path))?));

    let snps: Vec<SNP> = reader
        .deserialize()
        .into_iter()
        .filter_map(|r| {
            let snp: SNP = r.unwrap();
            match check_snp_order(&snp, reference) {
                Ok(SNPOrder::Reference) => Some(snp),
                Ok(SNPOrder::Alternate) => Some(snp_as_complementary(snp)),
                Err(_) => None,
            }
        })
        .collect();

    Ok(snps)
}

fn write_snps(path: &str, snps: Vec<SNP>) -> Result<(), Box<dyn Error>> {
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_writer(BufWriter::new(File::create(Path::new(path))?));

    // Write out header
    writer.write_record(&["CHROM", "POS", "REF", "ALT", "AF", "BETA"])?;

    for snp in snps {
        // Convert from EFFECT, OTHER (ALT, REF) to REF, ALT to make downstream comparison against
        // VCFs easier
        writer.serialize((snp.chr, snp.pos, snp.other, snp.effect, snp.af, snp.beta))?;
    }

    writer.flush()?;

    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();

    // const FASTA: &str = "data/hs_ref_GRCh37.p5_all_contigs.fa";
    // const SNP_IN: &str = "data/EA4_additive.sorted.tsv";
    // const SNP_OUT: &str = "out.tsv";

    let repository = open_fasta(&args.fasta)?;
    let snps = read_snps(&args.snp_in, &repository)?;
    write_snps(&args.out, snps)?;

    Ok(())
}
