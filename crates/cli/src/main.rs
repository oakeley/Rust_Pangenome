use clap::{Parser, Subcommand};
use anyhow::Result;
use core::{hilbert::HilbertMapper, index::HierarchicalIndex, io::load_index_hdf5};

#[derive(Parser)]
#[command(name="pango-cli", version, about="Fractal Pangenome CLI")]
struct Args { #[arg(long, default_value="./data/pangenome.h5")] db: String, #[arg(long, default_value_t=48)] threads: usize, #[command(subcommand)] cmd: Cmd }

#[derive(Subcommand)] enum Cmd { Stats }

fn main() -> Result<()> {
    let a = Args::parse();
    rayon::ThreadPoolBuilder::new().num_threads(a.threads).build_global().ok();
    let idx = HierarchicalIndex::new(HilbertMapper::new(256), &[1,4,16,64]);
    load_index_hdf5(&a.db, &idx)?;
    match a.cmd { Cmd::Stats => { println!("exact={}, levels={}, roads={}", idx.exact.len(), idx.levels.len(), idx.road_graph.edges.len()); } }
    Ok(())
}
