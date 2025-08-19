use dashmap::DashMap;
use petgraph::graphmap::DiGraphMap;
use crate::model::SegmentId; // Remove unused RoadEdge import
use parking_lot::RwLock;

#[derive(Default, Debug)] // Add Debug here
pub struct RoadGraph {
    // adjacency: directed edges -> importance counter
    pub edges: DashMap<(SegmentId, SegmentId), u64>,
    // per-genome paths (store vector of SegmentId per genome to reconstruct roads)
    pub genome_paths: DashMap<String, Vec<SegmentId>>,
    // fast lookup adjacency map for traversal (constructed on demand)
    pub adj: RwLock<Option<DiGraphMap<SegmentId, u64>>>,
}

impl RoadGraph {
    pub fn new() -> Self { Self::default() }

    pub fn add_genome_path(&self, genome_id: &str, path: Vec<SegmentId>) {
        // store the path
        self.genome_paths.insert(genome_id.to_string(), path.clone());
        // increase edge importance for successive pairs
        for w in path.windows(2) {
            let from = w[0]; let to = w[1];
            self.edges.entry((from,to)).and_modify(|e| { *e += 1 }).or_insert(1);
        }
        // invalidate adj cache
        let mut a = self.adj.write();
        *a = None;
    }

    pub fn downgrade_or_add_edge(&self, from: SegmentId, to: SegmentId) {
        self.edges.entry((from,to)).and_modify(|e| { *e += 1 }).or_insert(1);
        let mut a = self.adj.write();
        *a = None;
    }

    pub fn build_adj_if_needed(&self) {
        let mut a = self.adj.write();
        if a.is_some() { return; }
        let mut g = DiGraphMap::new();
        for r in self.edges.iter() {
            let ((f,t), imp) = r.pair();
            g.add_edge(*f, *t, *imp);
        }
        *a = Some(g);
    }

    pub fn neighbors(&self, seg: SegmentId) -> Vec<SegmentId> {
        self.build_adj_if_needed();
        if let Some(ref g) = *self.adj.read() {
            return g.neighbors(seg).collect();
        }
        vec![]
    }

    pub fn edge_importance(&self, from: SegmentId, to: SegmentId) -> u64 {
        self.edges.get(&(from,to)).map(|v| *v.value()).unwrap_or(0)
    }
}
