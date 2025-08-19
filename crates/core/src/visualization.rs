use crate::hilbert::HilbertMapper;
use crate::enhanced_annotations::{EnhancedAnnotationStore, FeatureType, KaryotypeInfo, StainingPattern, ChromosomeArm};
use crate::model::SegmentId;
use serde::{Serialize, Deserialize};
use std::collections::HashMap;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VisualizationLayer {
    pub layer_type: LayerType,
    pub zoom_level: u8, // 0-255, where 0 is most zoomed out
    pub elements: Vec<VisualElement>,
    pub style: LayerStyle,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum LayerType {
    Karyotype,
    Genes,
    RegulatoryElements,
    Variants,
    Expression,
    Conservation,
    Roads, // The k-mer segment connections
    Density,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VisualElement {
    pub id: String,
    pub hilbert_coord: (u16, u16),
    pub size: f64,
    pub color: [u8; 4], // RGBA
    pub shape: ElementShape,
    pub metadata: HashMap<String, String>,
    pub connected_to: Vec<String>, // For showing connections
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ElementShape {
    Circle,
    Rectangle,
    Arrow { direction: f64 }, // Angle in radians for gene direction
    Line { to: (u16, u16) },
    Polygon { points: Vec<(f64, f64)> },
    Chromosome { arm_ratio: f64 }, // For karyotype visualization
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LayerStyle {
    pub opacity: f64,
    pub stroke_width: f64,
    pub show_labels: bool,
    pub label_size: f64,
    pub color_scheme: ColorScheme,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ColorScheme {
    Default,
    Heatmap { min_color: [u8; 4], max_color: [u8; 4] },
    Categorical { colors: Vec<[u8; 4]> },
    Conservation { high: [u8; 4], medium: [u8; 4], low: [u8; 4] },
    Expression { low: [u8; 4], high: [u8; 4] },
}

pub struct MultiScaleVisualizer {
    pub hilbert: HilbertMapper,
    pub annotation_store: EnhancedAnnotationStore,
    pub layers: Vec<VisualizationLayer>,
    pub current_zoom: u8,
    pub viewport: ViewPort,
}

#[derive(Debug, Clone)]
pub struct ViewPort {
    pub center: (u16, u16),
    pub radius: u16,
    pub rotation: f64,
    pub aspect_ratio: f64,
}

impl MultiScaleVisualizer {
    pub fn new(hilbert: HilbertMapper, annotation_store: EnhancedAnnotationStore) -> Self {
        Self {
            hilbert,
            annotation_store,
            layers: Vec::new(),
            current_zoom: 128, // Start at medium zoom
            viewport: ViewPort {
                center: (128, 128), // Center of 256x256 Hilbert space
                radius: 64,
                rotation: 0.0,
                aspect_ratio: 1.0,
            },
        }
    }
    
    pub fn generate_karyotype_layer(&mut self, genome_id: &str) -> VisualizationLayer {
        let mut elements = Vec::new();
        
        // Get all karyotype annotations for this genome
        if let Some(genome_annotations) = self.annotation_store.by_segment.get(genome_id) {
            for segment_entry in genome_annotations.iter() {
                let (segment_id, annotations) = segment_entry.pair();
                
                for annotation in annotations.iter() {
                    if let FeatureType::Karyotype(karyotype_info) = &annotation.feature_type {
                        // Map genomic position to Hilbert coordinate
                        let hilbert_coord = self.estimate_hilbert_coord_for_segment(*segment_id);
                        
                        elements.push(VisualElement {
                            id: format!("karyotype_{}_{}", karyotype_info.chromosome, karyotype_info.cytogenetic_location),
                            hilbert_coord,
                            size: (karyotype_info.genomic_size as f64).log10() * 2.0, // Size based on genomic size
                            color: self.karyotype_color(&karyotype_info),
                            shape: ElementShape::Chromosome { 
                                arm_ratio: self.calculate_arm_ratio(&karyotype_info)
                            },
                            metadata: {
                                let mut meta = HashMap::new();
                                meta.insert("chromosome".to_string(), karyotype_info.chromosome.clone());
                                meta.insert("band".to_string(), karyotype_info.band.clone());
                                meta.insert("cytogenetic_location".to_string(), karyotype_info.cytogenetic_location.clone());
                                meta.insert("size".to_string(), karyotype_info.genomic_size.to_string());
                                meta
                            },
                            connected_to: Vec::new(),
                        });
                    }
                }
            }
        }
        
        VisualizationLayer {
            layer_type: LayerType::Karyotype,
            zoom_level: 0, // Visible at all zoom levels
            elements,
            style: LayerStyle {
                opacity: 0.8,
                stroke_width: 2.0,
                show_labels: true,
                label_size: 12.0,
                color_scheme: ColorScheme::Categorical { 
                    colors: self.generate_chromosome_colors() 
                },
            },
        }
    }
    
    pub fn generate_gene_layer(&mut self, genome_id: &str, zoom_level: u8) -> VisualizationLayer {
        let mut elements = Vec::new();
        
        // Only show genes at appropriate zoom levels
        if zoom_level < 50 {
            return VisualizationLayer {
                layer_type: LayerType::Genes,
                zoom_level,
                elements,
                style: LayerStyle::default(),
            };
        }
        
        // Query genes in current viewport
        let gene_segments = self.annotation_store.query_spatial_region(
            self.viewport.center,
            self.viewport.radius,
            &[FeatureType::Gene, FeatureType::Exon, FeatureType::CDS]
        );
        
        for segment_id in gene_segments {
            if let Some(genome_annotations) = self.annotation_store.by_segment.get(genome_id) {
                if let Some(annotations) = genome_annotations.get(&segment_id) {
                    for annotation in annotations.iter() {
                        match &annotation.feature_type {
                            FeatureType::Gene => {
                                let hilbert_coord = self.estimate_hilbert_coord_for_segment(segment_id);
                                
                                elements.push(VisualElement {
                                    id: format!("gene_{}", annotation.gene_id.as_ref().unwrap_or(&"unknown".to_string())),
                                    hilbert_coord,
                                    size: self.calculate_gene_size(annotation),
                                    color: self.gene_color(annotation),
                                    shape: ElementShape::Arrow { 
                                        direction: self.gene_direction(annotation) 
                                    },
                                    metadata: {
                                        let mut meta = HashMap::new();
                                        if let Some(gene_id) = &annotation.gene_id {
                                            meta.insert("gene_id".to_string(), gene_id.clone());
                                        }
                                        if let Some(gene_name) = &annotation.gene_name {
                                            meta.insert("gene_name".to_string(), gene_name.clone());
                                        }
                                        if let Some(expression) = annotation.expression_level {
                                            meta.insert("expression".to_string(), expression.to_string());
                                        }
                                        meta
                                    },
                                    connected_to: Vec::new(),
                                });
                            },
                            _ => {}
                        }
                    }
                }
            }
        }
        
        VisualizationLayer {
            layer_type: LayerType::Genes,
            zoom_level,
            elements,
            style: LayerStyle {
                opacity: 0.9,
                stroke_width: 1.0,
                show_labels: zoom_level > 100,
                label_size: 10.0,
                color_scheme: ColorScheme::Expression { 
                    low: [100, 100, 255, 255], 
                    high: [255, 100, 100, 255] 
                },
            },
        }
    }
    
    pub fn generate_road_network_layer(&mut self, _importance_threshold: u64) -> VisualizationLayer {
        let elements = Vec::new();
        
        // This would interface with your RoadGraph to show k-mer connections
        // For each edge with importance > threshold, create a line element
        
        VisualizationLayer {
            layer_type: LayerType::Roads,
            zoom_level: self.current_zoom,
            elements,
            style: LayerStyle {
                opacity: 0.6,
                stroke_width: 1.0,
                show_labels: false,
                label_size: 8.0,
                color_scheme: ColorScheme::Heatmap { 
                    min_color: [200, 200, 200, 100], 
                    max_color: [255, 0, 0, 255] 
                },
            },
        }
    }
    
    pub fn zoom_to_region(&mut self, center: (u16, u16), zoom_level: u8) {
        self.viewport.center = center;
        self.current_zoom = zoom_level;
        self.viewport.radius = std::cmp::max(1, 256 / (zoom_level as u16 + 1));
        
        // Regenerate visible layers at new zoom level
        self.update_layers();
    }
    
    pub fn zoom_to_gene(&mut self, gene_name: &str, genome_id: &str) {
        // First collect the gene segments to avoid borrowing conflicts
        let gene_segments: Vec<(String, SegmentId)> = self.annotation_store.by_gene_name
            .get(gene_name)
            .map(|segments| segments.clone())
            .unwrap_or_default();
            
        for (gid, segment_id) in gene_segments.iter() {
            if gid == genome_id {
                let hilbert_coord = self.estimate_hilbert_coord_for_segment(*segment_id);
                self.zoom_to_region(hilbert_coord, 200); // High zoom for gene view
                break;
            }
        }
    }
    
    pub fn zoom_to_chromosome(&mut self, chromosome: &str, genome_id: &str) {
        // Find all karyotype elements for this chromosome
        let mut chromosome_coords = Vec::new();
        
        if let Some(genome_annotations) = self.annotation_store.by_segment.get(genome_id) {
            for segment_entry in genome_annotations.iter() {
                let (_segment_id, annotations) = segment_entry.pair();
                
                for annotation in annotations.iter() {
                    if let FeatureType::Karyotype(karyotype_info) = &annotation.feature_type {
                        if karyotype_info.chromosome == chromosome {
                            let hilbert_coord = self.estimate_hilbert_coord_for_segment(*_segment_id);
                            chromosome_coords.push(hilbert_coord);
                        }
                    }
                }
            }
        }
        
        if !chromosome_coords.is_empty() {
            // Calculate centroid
            let center_x = chromosome_coords.iter().map(|(x, _)| *x as u32).sum::<u32>() / chromosome_coords.len() as u32;
            let center_y = chromosome_coords.iter().map(|(_, y)| *y as u32).sum::<u32>() / chromosome_coords.len() as u32;
            
            self.zoom_to_region((center_x as u16, center_y as u16), 100); // Medium zoom for chromosome view
        }
    }
    
    fn estimate_hilbert_coord_for_segment(&self, segment_id: SegmentId) -> (u16, u16) {
        // This should interface with your HierarchicalIndex to get the actual Hilbert coordinate
        // For now, returning a placeholder
        let idx = (segment_id % self.hilbert.total_positions as u64) as u32;
        self.hilbert.index_to_xy(idx)
    }
    
    fn karyotype_color(&self, karyotype_info: &KaryotypeInfo) -> [u8; 4] {
        // Color based on staining pattern
        match &karyotype_info.staining_pattern {
            StainingPattern::Gpos => [50, 50, 50, 255],     // Dark
            StainingPattern::Gneg => [200, 200, 200, 255],  // Light
            StainingPattern::Gvar => [100, 100, 100, 255],  // Variable
            StainingPattern::Acen => [255, 100, 100, 255],  // Centromere - red
            StainingPattern::Stalk => [100, 255, 100, 255], // Stalk - green
        }
    }
    
    fn calculate_arm_ratio(&self, karyotype_info: &KaryotypeInfo) -> f64 {
        // This would calculate p-arm vs q-arm ratio for chromosome visualization
        match &karyotype_info.arm {
            ChromosomeArm::P => 0.4, // p-arm is typically shorter
            ChromosomeArm::Q => 0.6, // q-arm is typically longer
        }
    }
    
    fn generate_chromosome_colors(&self) -> Vec<[u8; 4]> {
        // Generate distinct colors for each chromosome
        vec![
            [255, 100, 100, 255], [100, 255, 100, 255], [100, 100, 255, 255],
            [255, 255, 100, 255], [255, 100, 255, 255], [100, 255, 255, 255],
            // Add more colors for all chromosomes...
        ]
    }
    
    fn calculate_gene_size(&self, annotation: &crate::enhanced_annotations::EnhancedAnnotationTag) -> f64 {
        // Size based on expression level or conservation
        if let Some(expression) = annotation.expression_level {
            2.0 + expression * 8.0 // Size between 2-10 based on expression
        } else {
            4.0 // Default size
        }
    }
    
    fn gene_color(&self, annotation: &crate::enhanced_annotations::EnhancedAnnotationTag) -> [u8; 4] {
        // Color based on expression level
        if let Some(expression) = annotation.expression_level {
            let intensity = (expression * 255.0) as u8;
            [intensity, 100, 255 - intensity, 255]
        } else {
            [150, 150, 150, 255] // Gray for unknown expression
        }
    }
    
    fn gene_direction(&self, _annotation: &crate::enhanced_annotations::EnhancedAnnotationTag) -> f64 {
        // This would determine gene direction from strand information
        // For now, return a placeholder
        0.0 // 0 radians = pointing right (5' to 3')
    }
    
    fn update_layers(&mut self) {
        // Regenerate layers based on current zoom and viewport
        self.layers.clear();
        
        // Always show karyotype at low zoom
        if self.current_zoom < 150 {
            // Add karyotype layer generation here
        }
        
        // Show genes at higher zoom
        if self.current_zoom > 50 {
            // Add gene layer generation here
        }
        
        // Show regulatory elements at highest zoom
        if self.current_zoom > 150 {
            // Add regulatory element layer generation here
        }
    }
}

impl Default for LayerStyle {
    fn default() -> Self {
        Self {
            opacity: 1.0,
            stroke_width: 1.0,
            show_labels: true,
            label_size: 10.0,
            color_scheme: ColorScheme::Default,
        }
    }
}
