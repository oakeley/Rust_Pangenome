use serde::{Serialize, Deserialize};
use std::collections::HashMap;
use dashmap::DashMap;
use crate::model::SegmentId;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EnhancedAnnotationTag {
    pub feature_type: FeatureType,
    pub gene_id: Option<String>,
    pub gene_name: Option<String>,
    pub transcript_id: Option<String>,
    pub protein_id: Option<String>,
    pub functional_class: Option<String>,
    pub regulatory_type: Option<RegulatoryType>,
    pub clinical_significance: Option<ClinicalSignificance>,
    pub conservation_score: Option<f64>,
    pub expression_level: Option<f64>,
    pub tissue_specificity: Vec<String>,
    pub pathway_associations: Vec<String>,
    pub attributes: HashMap<String, String>, // Flexible key-value for custom annotations
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum FeatureType {
    Gene,
    Exon,
    Intron,
    CDS,
    UTR5,
    UTR3,
    Promoter,
    Enhancer,
    Silencer,
    Insulator,
    TFBS, // Transcription Factor Binding Site
    CpGIsland,
    RepeatElement,
    SNV,
    Indel,
    StructuralVariant,
    CopyNumberVariant,
    Karyotype(KaryotypeInfo),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum RegulatoryType {
    Promoter { strength: f64 },
    Enhancer { target_genes: Vec<String> },
    Silencer { repressed_genes: Vec<String> },
    Insulator,
    CTCF,
    PolII,
    H3K4me3,
    H3K27ac,
    H3K9me3,
    DNaseHypersensitive,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ClinicalSignificance {
    Pathogenic,
    LikelyPathogenic,
    VUS, // Variant of Uncertain Significance
    LikelyBenign,
    Benign,
    DrugResponse,
    RiskFactor,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KaryotypeInfo {
    pub chromosome: String,
    pub arm: ChromosomeArm,
    pub band: String,
    pub sub_band: Option<String>,
    pub cytogenetic_location: String, // e.g., "17q21.31"
    pub genomic_size: u64,
    pub staining_pattern: StainingPattern,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ChromosomeArm {
    P, // Short arm
    Q, // Long arm
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum StainingPattern {
    Gpos, // G-positive (dark bands)
    Gneg, // G-negative (light bands)
    Gvar, // G-variable
    Acen, // Centromere
    Stalk,
}

#[derive(Default)]
pub struct EnhancedAnnotationStore {
    // Hierarchical organization by scale
    pub by_segment: DashMap<String, DashMap<SegmentId, Vec<EnhancedAnnotationTag>>>,
    
    // Spatial indices for different feature types
    pub gene_spatial_index: DashMap<(u16, u16), Vec<SegmentId>>, // Hilbert coordinates -> gene segments
    pub regulatory_spatial_index: DashMap<(u16, u16), Vec<SegmentId>>,
    pub variant_spatial_index: DashMap<(u16, u16), Vec<SegmentId>>,
    pub karyotype_spatial_index: DashMap<(u16, u16), Vec<SegmentId>>,
    
    // Feature type indices for rapid querying
    pub by_feature_type: DashMap<String, Vec<(String, SegmentId)>>, // feature_type -> [(genome_id, segment_id)]
    pub by_gene_name: DashMap<String, Vec<(String, SegmentId)>>,
    pub by_pathway: DashMap<String, Vec<(String, SegmentId)>>,
    pub by_clinical_significance: DashMap<String, Vec<(String, SegmentId)>>,
}

impl EnhancedAnnotationStore {
    pub fn new() -> Self {
        Self::default()
    }
    
    pub fn attach_enhanced(&self, genome_id: &str, seg: SegmentId, tag: EnhancedAnnotationTag, hilbert_coord: (u16, u16)) {
        // Store in main index
        let inner = self.by_segment.entry(genome_id.to_string()).or_insert_with(DashMap::new);
        inner.entry(seg).or_default().push(tag.clone());
        
        // Update spatial indices based on feature type
        match &tag.feature_type {
            FeatureType::Gene | FeatureType::Exon | FeatureType::CDS => {
                self.gene_spatial_index.entry(hilbert_coord).or_default().push(seg);
            },
            FeatureType::Promoter | FeatureType::Enhancer | FeatureType::Silencer => {
                self.regulatory_spatial_index.entry(hilbert_coord).or_default().push(seg);
            },
            FeatureType::SNV | FeatureType::Indel | FeatureType::StructuralVariant => {
                self.variant_spatial_index.entry(hilbert_coord).or_default().push(seg);
            },
            FeatureType::Karyotype(_) => {
                self.karyotype_spatial_index.entry(hilbert_coord).or_default().push(seg);
            },
            _ => {}
        }
        
        // Update feature type indices
        let feature_key = format!("{:?}", tag.feature_type);
        self.by_feature_type.entry(feature_key).or_default().push((genome_id.to_string(), seg));
        
        if let Some(gene_name) = &tag.gene_name {
            self.by_gene_name.entry(gene_name.clone()).or_default().push((genome_id.to_string(), seg));
        }
        
        for pathway in &tag.pathway_associations {
            self.by_pathway.entry(pathway.clone()).or_default().push((genome_id.to_string(), seg));
        }
        
        if let Some(clinical_sig) = &tag.clinical_significance {
            let sig_key = format!("{:?}", clinical_sig);
            self.by_clinical_significance.entry(sig_key).or_default().push((genome_id.to_string(), seg));
        }
    }
    
    pub fn query_spatial_region(&self, center: (u16, u16), radius: u16, feature_types: &[FeatureType]) -> Vec<SegmentId> {
        let mut results = Vec::new();
        
        for feature_type in feature_types {
            let spatial_index = match feature_type {
                FeatureType::Gene | FeatureType::Exon | FeatureType::CDS => &self.gene_spatial_index,
                FeatureType::Promoter | FeatureType::Enhancer | FeatureType::Silencer => &self.regulatory_spatial_index,
                FeatureType::SNV | FeatureType::Indel | FeatureType::StructuralVariant => &self.variant_spatial_index,
                FeatureType::Karyotype(_) => &self.karyotype_spatial_index,
                _ => continue,
            };
            
            // Search in spatial radius
            for dx in -(radius as i32)..=(radius as i32) {
                for dy in -(radius as i32)..=(radius as i32) {
                    let x = (center.0 as i32 + dx).max(0) as u16;
                    let y = (center.1 as i32 + dy).max(0) as u16;
                    
                    if let Some(segments) = spatial_index.get(&(x, y)) {
                        results.extend(segments.iter().copied());
                    }
                }
            }
        }
        
        results.sort_unstable();
        results.dedup();
        results
    }
    
    pub fn find_genes_in_pathway(&self, pathway: &str) -> Vec<(String, SegmentId)> {
        self.by_pathway.get(pathway).map(|v| v.clone()).unwrap_or_default()
    }
    
    pub fn find_clinical_variants(&self, significance: ClinicalSignificance) -> Vec<(String, SegmentId)> {
        let key = format!("{:?}", significance);
        self.by_clinical_significance.get(&key).map(|v| v.clone()).unwrap_or_default()
    }
}
