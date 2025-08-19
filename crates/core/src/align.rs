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
pub fn generate_pangenome_vcf(index: &HierarchicalIndex, genomes: &[String]) -> Vec<VariantCall> {
    let mut variants = Vec::new();
    
    // Analyze the road graph to identify variants
    // This would compare paths between different genomes to find differences
    
    // For now, generate example variants from the road graph
    for edge_entry in index.road_graph.edges.iter().take(100) {
        let ((from_seg, to_seg), importance) = edge_entry.pair();
        
        // Check if this edge represents a potential variant
        if *importance < 10 { // Low frequency might indicate a variant
            if let (Some(from_meta), Some(to_meta)) = (
                index.meta.get(from_seg),
                index.meta.get(to_seg)
            ) {
                // Generate a variant call
                let variant = VariantCall {
                    position: (*from_seg % 1000000) as u32 + 1000000, // Placeholder position
                    reference_allele: "A".to_string(), // Would extract from sequence
                    alternate_allele: "T".to_string(), // Would extract from sequence
                    variant_type: VariantType::SNV,
                    allele_frequency: *importance as f64 / genomes.len() as f64,
                    supporting_segments: vec![*from_seg, *to_seg],
                    novel: *importance == 1, // Single occurrence = novel
                    quality_score: 30.0,
                };
                
                variants.push(variant);
            }
        }
    }
    
    variants
}

fn parse_reads(reads_text: &str) -> Vec<(String, Vec<u8>)> {
    let mut reads = Vec::new();
    let lines: Vec<&str> = reads_text.lines().collect();
    
    let mut i = 0;
    while i < lines.len() {
        let line = lines[i].trim();
        
        if line.starts_with('>') {
            // FASTA format
            let read_id = line[1..].to_string();
            i += 1;
            if i < lines.len() {
                let sequence = lines[i].trim().as_bytes().to_vec();
                reads.push((read_id, sequence));
            }
        } else if line.starts_with('@') {
            // FASTQ format
            let read_id = line[1..].to_string();
            i += 1;
            if i < lines.len() {
                let sequence = lines[i].trim().as_bytes().to_vec();
                reads.push((read_id, sequence));
                i += 2; // Skip quality line and + line
            }
        }
        
        i += 1;
    }
    
    reads
}

fn align_single_read_enhanced(index: &HierarchicalIndex, read_id: &str, seq: &[u8], k: usize, top_n: usize) -> Vec<EnhancedAlignHit> {
    let kmers = encode_kmers(seq, k);
    
    // Gather candidate segments
    let mut candidates: Vec<SegmentId> = Vec::new();
    for chunk in kmers.chunks(128) {
        let v: Vec<u32> = chunk.iter().copied().collect();
        if v.is_empty() { continue; }
        if let Some(sid) = index.exact_match(&v) { 
            candidates.push(sid); 
        } else {
            candidates.extend(index.candidates_by_tile(&v, 8));
        }
    }
    candidates.sort_unstable();
    candidates.dedup();
    
    // Provider to reconstruct sequence
    let provider = |sid: SegmentId| -> Option<Vec<u8>> {
        index.seg_seqs.get(&sid).map(|kvec| reconstruct_seq_from_kmers(&kvec, index.meta.get(&sid).map(|m| m.k as usize).unwrap_or(31)))
    };
    
    // Compute alignments in parallel per candidate
    let mut cand_scores: Vec<(SegmentId, EnhancedAlignHit)> = candidates.par_iter().filter_map(|&c| {
        if let Some(target_seq) = provider(c) {
            let (sc, cigar) = sw_affine_banded(seq, &target_seq, 128, 2, -3, -5, -2);
            if sc > 0 {
                if let Some(meta) = index.meta.get(&c) {
                    let hit = EnhancedAlignHit {
                        segment_id: c,
                        genome_id: meta.genome_id.clone(),
                        score: sc,
                        cigar,
                        alignment_type: AlignmentType::Genomic,
                        splice_sites: Vec::new(),
                        strand: None,
                        mapping_quality: calculate_mapping_quality(sc, seq.len()),
                        hilbert_coordinates: meta.hilbert_xy,
                        segment_path: vec![c],
                        road_importance: vec![index.road_graph.edge_importance(c, c).max(1)],
                        variant_calls: detect_variants_in_alignment(seq, &target_seq),
                    };
                    return Some((c, hit));
                }
            }
        }
        None
    }).collect();
    
    // Also attempt split-read detection
    let split_hits = split_read_align_enhanced(seq, provider, &candidates);
    for (sid, sc, cigar, splice_sites) in split_hits {
        if let Some(meta) = index.meta.get(&sid) {
            let hit = EnhancedAlignHit {
                segment_id: sid,
                genome_id: meta.genome_id.clone(),
                score: sc,
                cigar,
                alignment_type: AlignmentType::Spliced,
                splice_sites,
                strand: None,
                mapping_quality: calculate_mapping_quality(sc, seq.len()),
                hilbert_coordinates: meta.hilbert_xy,
                segment_path: vec![sid],
                road_importance: vec![index.road_graph.edge_importance(sid, sid).max(1)],
                variant_calls: Vec::new(),
            };
            cand_scores.push((sid, hit));
        }
    }
    
    // Sort and take top results
    cand_scores.sort_unstable_by_key(|(_, hit)| std::cmp::Reverse(hit.score));
    cand_scores.into_iter().take(top_n).map(|(_, hit)| hit).collect()
}

fn convert_to_sam_record(hit: EnhancedAlignHit, index: &HierarchicalIndex) -> PangenomeSamRecord {
    let seq_len = hit.segment_path.len() * 100; // Estimate sequence length
    
    PangenomeSamRecord {
        qname: format!("read_{}", hit.segment_id),
        flag: 0, // Mapped, primary alignment
        rname: format!("segment_{}", hit.segment_id),
        pos: 1, // 1-based position
        mapq: hit.mapping_quality,
        cigar: hit.cigar,
        rnext: "*".to_string(),
        pnext: 0,
        tlen: 0,
        seq: "*".to_string(), // Would need actual sequence
        qual: "*".to_string(), // Would need quality scores
        
        // Pangenome-specific fields
        pangenome_path: hit.segment_path,
        hilbert_coords: vec![hit.hilbert_coordinates],
        road_importance: hit.road_importance,
        genome_matches: vec![hit.genome_id],
        alignment_type: hit.alignment_type,
        variant_calls: hit.variant_calls,
    }
}

fn split_read_align_enhanced(a: &[u8], b_seq_provider: impl Fn(SegmentId)->Option<Vec<u8>>, candidates: &[SegmentId]) -> Vec<(SegmentId, i32, String, Vec<SpliceSite>)> {
    let mut out = Vec::new();
    for &c in candidates {
        if let Some(b) = b_seq_provider(c) {
            let (sc, cigar) = sw_affine_banded(a, &b, 64, 2, -3, -5, -2);
            if sc > 20 { 
                out.push((c, sc, cigar, Vec::new())); 
                continue; 
            }
            
            // Try split alignment
            let L = a.len();
            for split in [L/3, L/2, 2*L/3] {
                if split <= 1 || split >= L-1 { continue; }
                let left = &a[..split]; 
                let right = &a[split..];
                let (sc1, cigar1) = sw_affine_banded(left, &b, 64, 2, -3, -5, -2);
                let (sc2, cigar2) = sw_affine_banded(right, &b, 64, 2, -3, -5, -2);
                if sc1 > 12 && sc2 > 12 {
                    let cigar_combined = format!("{}+{}", cigar1, cigar2);
                    let splice_sites = detect_splice_sites(left, right, &b);
                    out.push((c, sc1+sc2, cigar_combined, splice_sites));
                    break;
                }
            }
        }
    }
    out
}

fn detect_splice_sites(left: &[u8], right: &[u8], _target: &[u8]) -> Vec<SpliceSite> {
    let mut splice_sites = Vec::new();
    
    // Look for canonical splice motifs
    if left.len() >= 2 && right.len() >= 2 {
        let left_end = &left[left.len()-2..];
        let right_start = &right[..2];
        
        let is_canonical = matches!(
            (left_end, right_start),
            (b"GT", b"AG") | (b"GC", b"AG") | (b"AT", b"AC")
        );
        
        if is_canonical {
            splice_sites.push(SpliceSite {
                donor_segment: 0, // Would need actual segment IDs
                acceptor_segment: 0,
                confidence: 0.9,
                canonical: true,
                junction_sequence: format!("{}...{}", 
                    String::from_utf8_lossy(left_end),
                    String::from_utf8_lossy(right_start)
                ),
            });
        }
    }
    
    splice_sites
}

fn detect_variants_in_alignment(read_seq: &[u8], ref_seq: &[u8]) -> Vec<VariantCall> {
    let mut variants = Vec::new();
    
    // Simple variant detection - look for mismatches
    let min_len = read_seq.len().min(ref_seq.len());
    for i in 0..min_len {
        if read_seq[i] != ref_seq[i] {
            variants.push(VariantCall {
                position: i as u32,
                reference_allele: String::from_utf8_lossy(&ref_seq[i..i+1]).to_string(),
                alternate_allele: String::from_utf8_lossy(&read_seq[i..i+1]).to_string(),
                variant_type: VariantType::SNV,
                allele_frequency: 0.5, // Unknown frequency
                supporting_segments: Vec::new(),
                novel: false, // Would need population data to determine
                quality_score: 20.0,
            });
        }
    }
    
    variants
}

fn calculate_mapping_quality(score: i32, read_length: usize) -> u8 {
    let max_possible_score = read_length as i32 * 2;
    let quality = ((score as f64 / max_possible_score as f64) * 60.0) as u8;
    quality.min(60)
}

// Keep your original sw_affine_banded function
fn sw_affine_banded(a: &[u8], b: &[u8], band: usize, match_score: i32, mismatch: i32, gap_open: i32, gap_ext: i32) -> (i32, String) {
    let n = a.len(); let m = b.len();
    if n==0 || m==0 { return (0, String::new()); }
    
    let mut best = 0; let mut bi = 0usize; let mut bj = 0usize;
    let neginf = std::i32::MIN/4;
    let mut M = vec![vec![0i32; m+1]; n+1];
    let Ix = vec![vec![neginf; m+1]; n+1];
    let Iy = vec![vec![neginf; m+1]; n+1];
    let mut TB = vec![vec![0u8; m+1]; n+1];
    
    for i in 1..=n {
        let j0 = if i>band { i-band } else { 1 };
        let j1 = std::cmp::min(m, i+band);
        for j in j0..=j1 {
            let score = if a[i-1]==b[j-1] { match_score } else { mismatch };
            let m_val = (M[i-1][j-1].max(Ix[i-1][j-1]).max(Iy[i-1][j-1])) + score;
            let ix_val = (M[i-1][j] + gap_open + gap_ext).max(Ix[i-1][j] + gap_ext);
            let iy_val = (M[i][j-1] + gap_open + gap_ext).max(Iy[i][j-1] + gap_ext);
            let v = *[m_val, ix_val, iy_val, 0].iter().max().unwrap();
            M[i][j] = v;
            
            if v==0 { TB[i][j]=0; } else if v==m_val { TB[i][j]=1; } else if v==ix_val { TB[i][j]=2; } else { TB[i][j]=3; }
            if v > best { best = v; bi = i; bj = j; }
        }
    }
    if best==0 { return (0, String::new()); }
    
    // Traceback
    let mut ops: Vec<(char,usize)> = Vec::new();
    let mut i = bi; let mut j = bj;
    while i>0 && j>0 && TB[i][j]!=0 {
        let dir = TB[i][j];
        if dir==1 {
            if let Some((op,c)) = ops.last_mut() { if *op=='M' { *c+=1 } else { ops.push(('M',1)); } } else { ops.push(('M',1)); }
            i -= 1; j -= 1;
        } else if dir==2 {
            if let Some((op,c)) = ops.last_mut() { if *op=='D' { *c+=1 } else { ops.push(('D',1)); } } else { ops.push(('D',1)); }
            i -= 1;
        } else {
            if let Some((op,c)) = ops.last_mut() { if *op=='I' { *c+=1 } else { ops.push(('I',1)); } } else { ops.push(('I',1)); }
            j -= 1;
        }
    }
    ops.reverse();
    let cigar = ops.into_iter().map(|(op,c)| format!("{}{}", c, op)).collect::<Vec<_>>().join(""); 
    (best, cigar)
}
