use serde::{Serialize, Deserialize};
pub type SegmentId = u64;
pub type RoadId = u64;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SegmentMeta {
    pub genome_id: String,
    pub length: u32,
    pub frequency: u32,
    pub hilbert_xy: (u16, u16),
    pub k: u8,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RoadEdge {
    pub from: SegmentId,
    pub to: SegmentId,
    pub importance: u64, // how many genomes reuse this edge
}
