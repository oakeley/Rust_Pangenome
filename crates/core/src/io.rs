use anyhow::Result;
use hdf5::File;
use crate::{index::HierarchicalIndex, model::SegmentId};
use std::collections::HashMap;
use log::info;

pub fn load_index_hdf5(db: &str, idx: &HierarchicalIndex) -> Result<usize> {
    let f = File::open(db)?;
    let root = f.group("/")?;
    if !root.member_names()?.contains(&"path_segments".to_string()) { return Ok(0); }
    let segs = root.group("path_segments")?;
    let mut loaded = 0usize;
    // Load segments
    for name in segs.member_names()? {
        let g = segs.group(&name)?;
        // Remove the ? operator from link_exists
        if g.link_exists("kmer_sequence") {
            let dset = g.dataset("kmer_sequence")?;
            let vec: Vec<u32> = dset.read_raw()?.to_vec();
            let freq: u32 = g.attr("frequency").ok().and_then(|a| a.read_scalar().ok()).unwrap_or(1);
            
            // Fix the string reading issue
            let genome_id: String = match g.attr("genome_id") {
                Ok(attr) => {
                    // Try reading as variable-length string first
                    if let Ok(s) = attr.read_scalar::<hdf5::types::VarLenUnicode>() {
                        s.to_string()
                    } else if let Ok(s) = attr.read_scalar::<hdf5::types::VarLenAscii>() {
                        s.to_string()
                    } else if let Ok(s) = attr.read_scalar::<hdf5::types::FixedAscii<64>>() {
                        s.to_string()
                    } else {
                        "unknown".to_string()
                    }
                },
                Err(_) => "unknown".to_string(),
            };
            
            let k: u8 = g.attr("k").ok().and_then(|a| a.read_scalar().ok()).unwrap_or(31u8);
            let id: SegmentId = name.parse().unwrap_or(0);
            idx.add_segment(id, &genome_id, &vec, freq, k);
            loaded += 1;
        }
    }
    info!("Loaded {loaded} segments");
    // Load genome paths if present
    if root.member_names()?.contains(&"paths".to_string()) {
        let paths = root.group("paths")?;
        for gname in paths.member_names()? {
            let g = paths.group(&gname)?;
            // Remove the ? operator from link_exists
            if g.link_exists("segments") {
                let segs: Vec<u64> = g.dataset("segments")?.read_raw()?.to_vec();
                let genome_id = gname.clone();
                let segs_u64: Vec<SegmentId> = segs.iter().copied().collect();
                idx.add_genome_path_from_segments(&genome_id, segs_u64);
            }
        }
    } else {
        // Try to build simple paths per-genome by ordering segments by their hilbert projection
        let mut by_genome: HashMap<String, Vec<(u64,(u16,u16))>> = HashMap::new();
        for r in idx.meta.iter() {
            let (sid, m) = (r.key().clone(), r.value().clone());
            by_genome.entry(m.genome_id.clone()).or_default().push((sid, m.hilbert_xy));
        }
        for (genome, mut list) in by_genome.into_iter() {
            // sort by hilbert xy (row-major) as a heuristic path
            list.sort_by_key(|(_,xy)| (xy.0 as u32)<<16 | (xy.1 as u32));
            let segs: Vec<SegmentId> = list.into_iter().map(|(s,_)| s).collect();
            if segs.len() > 1 {
                idx.add_genome_path_from_segments(&genome, segs);
            }
        }
    }
    Ok(loaded)
}
