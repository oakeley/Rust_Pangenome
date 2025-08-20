use anyhow::*;
use clap::Parser;
use std::path::PathBuf;
use fractal_pangenome_core::index::HierIndex;
use fractal_pangenome_core::io::{save_index};
use needletail::{parse_fastx_file, Sequence};
use indicatif::{ProgressBar, ProgressStyle};

#[derive(Parser, Debug)]
struct Args {
    /// Path to database (bincode-serialized index)
    #[arg(long)]
    db: PathBuf,
    /// Add genome FASTA
    #[arg(long)]
    add_genome: PathBuf,
    /// Genome ID
    #[arg(long)]
    genome_id: String,
    /// k-mer length
    #[arg(long, default_value_t = 31)]
    k: u8,
    /// Threads
    #[arg(long, default_value_t = 0)]
    threads: usize,
}

fn main() -> Result<()> {
    let args = Args::parse();
    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global().ok();
    }

    let idx = HierIndex::new(args.k);
    let pb = ProgressBar::new_spinner().with_style(ProgressStyle::default_spinner().template("{spinner} {msg}").unwrap());
    pb.enable_steady_tick(std::time::Duration::from_millis(100));
    pb.set_message("Reading FASTA and ingesting kmers...");

    let mut reader = parse_fastx_file(&args.add_genome)?;
    let mut kmers: Vec<Vec<u8>> = Vec::new();
    while let Some(record) = reader.next() {
        let rec = record?;
        let seq = rec.seq().to_vec();
        // naive kmerization (can be SIMD-optimized later)
        for i in 0..=seq.len().saturating_sub(args.k as usize) {
            kmers.push(seq[i..i+args.k as usize].to_vec());
        }
        if kmers.len() > 200_000 {
            idx.ingest_parallel(&args.genome_id, std::mem::take(&mut kmers))?;
        }
    }
    if !kmers.is_empty() {
        idx.ingest_parallel(&args.genome_id, kmers)?;
    }
    pb.finish_with_message("Ingestion complete. Saving index...");
    save_index(&args.db, &idx)?;
    eprintln!("Saved index to {}", args.db.display());
    Ok(())
}
