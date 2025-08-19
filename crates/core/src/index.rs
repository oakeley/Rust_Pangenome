use std::sync::Arc;
use dashmap::DashMap;
use hashbrown::{HashMap, HashSet};
use parking_lot::RwLock;
// Remove unused rayon import
use ahash::RandomState;
use crate::hilbert::HilbertMapper;
use crate::model::{SegmentId, SegmentMeta};
use crate::graph::RoadGraph;

#[derive(Debug)]
pub struct HierarchicalIndex {
    pub hilbert: HilbertMapper,
    pub exact: DashMap<u64, SegmentId>,
    pub levels: Vec<RwLock<HashMap<(u16,u16), Vec<SegmentId>, RandomState>>>,
    pub meta: DashMap<SegmentId, SegmentMeta>,
    pub seg_seqs: DashMap<SegmentId, Vec<u32>>,
    pub road_graph: Arc<RoadGraph>,
}

impl HierarchicalIndex {
    pub fn new(hilbert: HilbertMapper, zooms: &[u16]) -> Self {
        let levels = zooms.iter().map(|_| RwLock::new(HashMap::with_hasher(RandomState::new()))).collect();
        Self { hilbert, exact: DashMap::new(), levels, meta: DashMap::new(), seg_seqs: DashMap::new(), road_graph: Arc::new(RoadGraph::new()) }
    }
    #[inline]
    fn coords_to_grid(&self, xy: (u16,u16), zoom: u16) -> (u16,u16) {
        let grid = self.hilbert.size / zoom; (xy.0 / grid, xy.1 / grid)
    }
    pub fn add_segment(&self, seg_id: SegmentId, genome_id: &str, kmer_seq: &[u32], frequency: u32, k: u8) {
        if kmer_seq.is_empty() { return; }
        use xxhash_rust::xxh3::Xxh3;
        let mut h = Xxh3::default();
        let bytes = unsafe { std::slice::from_raw_parts(kmer_seq.as_ptr() as *const u8, kmer_seq.len()*4) };
        h.update(bytes);
        let key = h.digest();
        self.exact.insert(key, seg_id);
        let rep = kmer_seq[0] % self.hilbert.total_positions;
        let xy = self.hilbert.index_to_xy(rep);
        let zooms: [u16;4] = [1,4,16,64];
        for (i,z) in zooms.iter().enumerate() {
            let grid = self.coords_to_grid(xy, *z);
            let mut m = self.levels[i].write();
            m.entry(grid).or_default().push(seg_id);
        }
        self.meta.insert(seg_id, SegmentMeta { genome_id: genome_id.to_string(), length: kmer_seq.len() as u32, frequency, hilbert_xy: xy, k });
        self.seg_seqs.insert(seg_id, kmer_seq.to_vec());
    }

    pub fn add_genome_path_from_segments(&self, genome_id: &str, segs: Vec<SegmentId>) {
        // add to road graph, and ensure segments exist
        self.road_graph.add_genome_path(genome_id, segs);
    }

    pub fn exact_match(&self, kmer_seq: &[u32]) -> Option<SegmentId> {
        use xxhash_rust::xxh3::Xxh3;
        let mut h = Xxh3::default();
        let bytes = unsafe { std::slice::from_raw_parts(kmer_seq.as_ptr() as *const u8, kmer_seq.len()*4) };
        h.update(bytes);
        let key = h.digest();
        self.exact.get(&key).map(|e| *e.value())
    }
    pub fn candidates_by_tile(&self, kmer_seq: &[u32], max_n: usize) -> Vec<SegmentId> {
        if kmer_seq.is_empty() { return vec![]; }
        let rep = kmer_seq[0] % self.hilbert.total_positions;
        let xy = self.hilbert.index_to_xy(rep);
        let zooms: [u16;4] = [1,4,16,64];
        let mut seen = HashSet::new();
        let mut out = Vec::new();
        'outer: for (i,z) in zooms.iter().enumerate() {
            let grid = self.coords_to_grid(xy, *z);
            let map = self.levels[i].read();
            for dx in -1i16..=1 { for dy in -1i16..=1 {
                let nx = grid.0 as i16 + dx; let ny = grid.1 as i16 + dy;
                if nx<0 || ny<0 { continue; }
                if let Some(v)=map.get(&(nx as u16, ny as u16)) {
                    for &sid in v {
                        if seen.insert(sid) { out.push(sid); if out.len()>=max_n { break 'outer; } }
                    }
                }
            }}
        }
        out
    }
}
