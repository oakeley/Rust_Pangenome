use crate::index::HierarchicalIndex;
use crate::encoding::{encode_kmers, reconstruct_seq_from_kmers};
use crate::model::SegmentId;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use serde::{Serialize, Deserialize};

// Keep original AlignHit for backward compatibility
#[derive(Debug, serde::Serialize)]
pub struct AlignHit { 
    pub segment_id: SegmentId, 
    pub genome_id: String, 
    pub score: i32, 
    pub cigar: String 
}

// Enhanced alignment result with pangenome features
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EnhancedAlignHit {
    pub segment_id: SegmentId,
    pub genome_id: String,
    pub score: i32,
    pub cigar: String,
    pub alignment_type: AlignmentType,
    pub splice_sites: Vec<SpliceSite>,
    pub strand: Option<char>,
    pub mapping_quality: u8,
    pub hilbert_coordinates: (u16, u16),
    pub segment_path: Vec<SegmentId>,  // For multi-segment alignments
    pub road_importance: Vec<u64>,     // Importance of each segment
    pub variant_calls: Vec<VariantCall>, // Detected variants
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum AlignmentType {
    Genomic,
    Spliced,
    SplitRead,
    ChimericRead,
    RepeatRegion,
    StructuralVariant,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpliceSite {
    pub donor_segment: SegmentId,
    pub acceptor_segment: SegmentId,
    pub confidence: f64,
    pub canonical: bool,
    pub junction_sequence: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VariantCall {
    pub position: u32,
    pub reference_allele: String,
    pub alternate_allele: String,
    pub variant_type: VariantType,
    pub allele_frequency: f64,
    pub supporting_segments: Vec<SegmentId>,
    pub novel: bool,
    pub quality_score: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum VariantType {
    SNV,
    Insertion,
    Deletion,
    Complex,
    StructuralVariant,
}

// Pangenome-aware SAM record
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PangenomeSamRecord {
    pub qname: String,          // Query name
    pub flag: u16,              // Bitwise flags
    pub rname: String,          // Reference sequence name (segment_id)
    pub pos: u32,               // 1-based leftmost mapping position
    pub mapq: u8,               // Mapping quality
    pub cigar: String,          // CIGAR string
    pub rnext: String,          // Reference name of the mate/next read
    pub pnext: u32,             // Position of the mate/next read
    pub tlen: i32,              // Template length
    pub seq: String,            // Segment sequence
    pub qual: String,           // ASCII of Phred-scaled base qualities
    
    // Custom pangenome tags
    pub pangenome_path: Vec<SegmentId>,     // PG:Z: - Path through pangenome
    pub hilbert_coords: Vec<(u16, u16)>,    // HC:Z: - Hilbert coordinates
    pub road_importance: Vec<u64>,          // RI:Z: - Road importance values
    pub genome_matches: Vec<String>,        // GM:Z: - Matching genomes
    pub alignment_type: AlignmentType,      // AT:Z: - Alignment type
    pub variant_calls: Vec<VariantCall>,    // VC:Z: - Detected variants
}

// Keep original function signature for compatibility
pub fn align_reads_with_cigar_splice(index: &HierarchicalIndex, reads_text: &str, k: usize, top_n: usize) -> Vec<AlignHit> {
    let enhanced_hits = align_reads_enhanced(index, reads_text, k, top_n);
    enhanced_hits.into_iter().map(|hit| AlignHit {
        segment_id: hit.segment_id,
        genome_id: hit.genome_id,
        score: hit.score,
        cigar: hit.cigar,
    }).collect()
}

// Enhanced alignment function with full pangenome features
pub fn align_reads_enhanced(index: &HierarchicalIndex, reads_text: &str, k: usize, top_n: usize) -> Vec<EnhancedAlignHit> {
    let seqs: Vec<(String, Vec<u8>)> = parse_reads(reads_text);

    seqs.par_iter().flat_map_iter(|(read_id, seq)| {
        align_single_read_enhanced(index, read_id, seq, k, top_n)
    }).collect()
}

// Enhanced alignment with pangenome SAM output
pub fn align_reads_to_pangenome_sam(index: &HierarchicalIndex, reads_text: &str, k: usize, top_n: usize) -> Vec<PangenomeSamRecord> {
    let enhanced_hits = align_reads_enhanced(index, reads_text, k, top_n);
    
    enhanced_hits.into_iter().map(|hit| {
        convert_to_sam_record(hit, index)
    }).collect()
}

// Generate VCF from pangenome analysis
pub fn generate_pangen
