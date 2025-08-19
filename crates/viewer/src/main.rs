use eframe::egui::{self, Color32, Pos2, Stroke, Rect, Vec2};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

// Import from core - you'll need to make sure these are available
use core::{
    HilbertMapper, 
    EnhancedAnnotationStore, 
    MultiScaleVisualizer,
    HierarchicalIndex,
    align_reads_enhanced,
    AlignmentType,
};

#[derive(Default)]
struct App { 
    // Server connection
    server: String, 
    connected: bool, 
    
    // Legacy road query
    segment_id: u64, 
    radius: usize, 
    edges: Vec<Edge>,
    
    // Visualization system
    visualizer: Option<MultiScaleVisualizer>,
    current_zoom: f32,
    pan_offset: Vec2,
    selected_layers: HashMap<String, bool>,
    show_genes: bool,
    show_karyotype: bool,
    show_roads: bool,
    show_variants: bool,
    show_annotations: bool,
    
    // Navigation and search
    gene_search: String,
    current_genome: String,
    chromosome_view: String,
    
    // Database management
    database_path: String,
    database_loaded: bool,
    index: Option<HierarchicalIndex>,
    
    // File browser state
    show_file_browser: bool,
    selected_fasta_path: String,
    selected_gff_path: String,
    selected_annotation_path: String,
    new_genome_id: String,
    
    // Alignment and analysis
    reads_input: String,
    alignment_results: Vec<PangenomeAlignment>,
    vcf_results: Vec<VariantCall>,
    show_alignment_panel: bool,
    show_vcf_panel: bool,
    
    // Creation mode
    create_new_database: bool,
    new_database_path: String,
}

#[derive(Deserialize, Clone)]
struct Edge { 
    from: u64, 
    to: u64, 
    importance: u64 
}

// Enhanced alignment result with pangenome coordinates
#[derive(Debug, Clone, Serialize, Deserialize)]
struct PangenomeAlignment {
    pub read_id: String,
    pub segment_path: Vec<u64>,  // Path through pangenome segments
    pub hilbert_coordinates: Vec<(u16, u16)>,  // Hilbert space coordinates
    pub genome_matches: Vec<String>,  // Which genomes this path appears in
    pub alignment_score: i32,
    pub mapping_quality: u8,
    pub cigar: String,
    pub alignment_type: AlignmentType,
    pub road_importance: Vec<u64>,  // Importance of each segment in the path
    pub strand: Option<char>,
    pub start_pos: Option<u64>,  // Genomic start position if available
    pub end_pos: Option<u64>,    // Genomic end position if available
}

// VCF-like variant calls from pangenome
#[derive(Debug, Clone, Serialize, Deserialize)]
struct VariantCall {
    pub chromosome: String,
    pub position: u64,
    pub reference_allele: String,
    pub alternate_alleles: Vec<String>,
    pub allele_frequencies: Vec<f64>,  // Frequency in pangenome
    pub allele_counts: Vec<u32>,       // Count in pangenome
    pub total_genomes: u32,
    pub variant_type: VariantType,
    pub novel: bool,  // True if this variant is novel
    pub supporting_segments: Vec<u64>,  // Pangenome segments supporting this variant
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

impl App { 
    fn new() -> Self { 
        let mut app = Self { 
            server: "http://127.0.0.1:8787".into(), 
            connected: false, 
            segment_id: 0, 
            radius: 1, 
            edges: vec![],
            
            visualizer: None,
            current_zoom: 1.0,
            pan_offset: Vec2::ZERO,
            selected_layers: HashMap::new(),
            show_genes: true,
            show_karyotype: true,
            show_roads: false,
            show_variants: false,
            show_annotations: true,
            
            gene_search: String::new(),
            current_genome: "Human_T2T_CHM13".to_string(),
            chromosome_view: "chr1".to_string(),
            
            database_path: "./data/pangenome.h5".to_string(),
            database_loaded: false,
            index: None,
            
            show_file_browser: false,
            selected_fasta_path: String::new(),
            selected_gff_path: String::new(),
            selected_annotation_path: String::new(),
            new_genome_id: String::new(),
            
            reads_input: String::new(),
            alignment_results: Vec::new(),
            vcf_results: Vec::new(),
            show_alignment_panel: false,
            show_vcf_panel: false,
            
            create_new_database: false,
            new_database_path: String::new(),
        };
        
        // Initialize layer selections
        app.selected_layers.insert("Karyotype".to_string(), true);
        app.selected_layers.insert("Genes".to_string(), true);
        app.selected_layers.insert("Roads".to_string(), false);
        app.selected_layers.insert("Variants".to_string(), false);
        app.selected_layers.insert("Annotations".to_string(), true);
        
        app
    }
    
    fn initialize_visualizer(&mut self) {
        if self.index.is_some() {
            let hilbert = HilbertMapper::new(256);
            let annotation_store = EnhancedAnnotationStore::new();
            self.visualizer = Some(MultiScaleVisualizer::new(hilbert, annotation_store));
            println!("Visualizer initialized successfully");
        }
    }
    
    fn load_database(&mut self) {
        // Load the HDF5 database
        let hilbert = HilbertMapper::new(256);
        let index = HierarchicalIndex::new(hilbert, &[1,4,16,64]);
        
        if let Ok(_loaded_count) = core::io::load_index_hdf5(&self.database_path, &index) {
            self.index = Some(index);
            self.database_loaded = true;
            self.connected = true;
            self.initialize_visualizer();
            println!("Database loaded successfully from {}", self.database_path);
        } else {
            eprintln!("Failed to load database from {}", self.database_path);
        }
    }
    
    fn create_new_database_from_files(&mut self) {
        if !self.selected_fasta_path.is_empty() && !self.new_genome_id.is_empty() {
            // Create new database
            let hilbert = HilbertMapper::new(256);
            let index = HierarchicalIndex::new(hilbert, &[1,4,16,64]);
            
            // TODO: Implement FASTA parsing and k-mer generation
            // This would:
            // 1. Parse FASTA file
            // 2. Generate k-mers
            // 3. Add to hierarchical index
            // 4. If GFF3 provided, parse annotations
            // 5. Save to HDF5 format
            
            self.index = Some(index);
            self.database_loaded = true;
            self.initialize_visualizer();
            
            println!("Created new database with genome: {}", self.new_genome_id);
            println!("FASTA: {}", self.selected_fasta_path);
            if !self.selected_gff_path.is_empty() {
                println!("GFF3: {}", self.selected_gff_path);
            }
        }
    }
    
    fn align_reads_to_pangenome(&mut self) {
        if let Some(ref index) = self.index {
            if !self.reads_input.is_empty() {
                // Perform enhanced alignment
                let enhanced_hits = align_reads_enhanced(index, &self.reads_input, 31, 10);
                
                // Convert to pangenome alignments
                self.alignment_results = enhanced_hits.into_iter().enumerate().map(|(i, hit)| {
                    let hilbert_coord = if let Some(meta) = index.meta.get(&hit.segment_id) {
                        vec![meta.hilbert_xy]
                    } else {
                        vec![(0, 0)]
                    };
                    
                    PangenomeAlignment {
                        read_id: format!("read_{}", i),
                        segment_path: vec![hit.segment_id],
                        hilbert_coordinates: hilbert_coord,
                        genome_matches: vec![hit.genome_id],
                        alignment_score: hit.score,
                        mapping_quality: hit.mapping_quality,
                        cigar: hit.cigar,
                        alignment_type: hit.alignment_type,
                        road_importance: vec![1], // Would get from road graph
                        strand: hit.strand,
                        start_pos: None, // Would calculate from segment position
                        end_pos: None,
                    }
                }).collect();
                
                self.show_alignment_panel = true;
                println!("Aligned {} reads to pangenome", self.alignment_results.len());
            }
        }
    }
    
    fn generate_vcf_from_pangenome(&mut self) {
        if let Some(ref _index) = self.index {
            // Generate VCF-like variant calls from pangenome
            // This would analyze the road graph to find variants
            self.vcf_results = Vec::new();
            
            // Example variant generation (placeholder)
            for i in 0..10 {
                self.vcf_results.push(VariantCall {
                    chromosome: format!("chr{}", (i % 22) + 1),
                    position: 1000000 + i * 1000,
                    reference_allele: "A".to_string(),
                    alternate_alleles: vec!["T".to_string()],
                    allele_frequencies: vec![0.1],
                    allele_counts: vec![5],
                    total_genomes: 50,
                    variant_type: VariantType::SNV,
                    novel: i % 5 == 0, // Every 5th variant is novel
                    supporting_segments: vec![i as u64],
                    quality_score: 30.0,
                });
            }
            
            self.show_vcf_panel = true;
            println!("Generated {} variant calls from pangenome", self.vcf_results.len());
        }
    }

    fn handle_viewport_interactions(&mut self, response: &egui::Response, ctx: &egui::Context) {
        // Handle pan
        if response.dragged() {
            self.pan_offset += response.drag_delta() / self.current_zoom;
        }
        
        // Handle zoom
        if let Some(_hover_pos) = response.hover_pos() {
            let scroll_delta = ctx.input(|i| i.raw_scroll_delta.y);
            if scroll_delta != 0.0 {
                let zoom_factor = 1.0 + scroll_delta * 0.001;
                self.current_zoom *= zoom_factor;
                self.current_zoom = self.current_zoom.clamp(0.1, 100.0);
            }
        }
    }
    
    fn draw_all_layers(&self, painter: &egui::Painter, rect: &Rect, transform: &dyn Fn(Pos2) -> Pos2) {
        // Draw background grid
        self.draw_background_grid(painter, rect);
        
        // Draw layers in order
        if self.show_karyotype {
            self.draw_karyotype_layer(painter, transform);
        }
        
        if self.show_roads {
            self.draw_roads_layer(painter, transform);
        }
        
        if self.show_genes {
            self.draw_genes_layer(painter, transform);
        }
        
        if self.show_variants {
            self.draw_variants_layer(painter, transform);
        }
        
        if self.show_annotations {
            self.draw_annotations_layer(painter, transform);
        }
    }
    
    fn draw_background_grid(&self, painter: &egui::Painter, rect: &Rect) {
        let grid_size = 20.0 * self.current_zoom;
        if grid_size < 5.0 { return; } // Don't draw grid when too small
        
        let start_x = (rect.min.x / grid_size).floor() * grid_size;
        let start_y = (rect.min.y / grid_size).floor() * grid_size;
        
        let mut x = start_x;
        while x <= rect.max.x {
            painter.line_segment(
                [Pos2::new(x, rect.min.y), Pos2::new(x, rect.max.y)],
                Stroke::new(0.5, Color32::from_gray(30))
            );
            x += grid_size;
        }
        
        let mut y = start_y;
        while y <= rect.max.y {
            painter.line_segment(
                [Pos2::new(rect.min.x, y), Pos2::new(rect.max.x, y)],
                Stroke::new(0.5, Color32::from_gray(30))
            );
            y += grid_size;
        }
    }
    
    fn draw_karyotype_layer(&self, painter: &egui::Painter, transform: &dyn Fn(Pos2) -> Pos2) {
        // Draw chromosome ideograms
        for i in 1..=22 {
            let x = (i as f32 - 1.0) * 40.0;
            let rect = Rect::from_min_size(
                transform(Pos2::new(x, 20.0)),
                Vec2::new(30.0 * self.current_zoom, 120.0 * self.current_zoom)
            );
            
            let color = self.chromosome_color(i);
            painter.rect_filled(rect, 8.0, color);
            painter.rect_stroke(rect, 8.0, Stroke::new(2.0, Color32::BLACK));
            
            // Centromere
            let centromere_y = rect.center().y;
            painter.circle_filled(Pos2::new(rect.center().x, centromere_y), 5.0 * self.current_zoom, Color32::RED);
            
            // Labels
            if self.current_zoom > 0.5 {
                painter.text(
                    rect.center() + Vec2::new(0.0, 80.0 * self.current_zoom),
                    egui::Align2::CENTER_CENTER,
                    format!("{}", i),
                    egui::FontId::proportional(12.0),
                    Color32::WHITE,
                );
            }
        }
        
        // Sex chromosomes
        for (i, chr) in ["X", "Y"].iter().enumerate() {
            let x = 22.0 * 40.0 + (i as f32) * 40.0;
            let height = if chr == &"X" { 100.0 } else { 70.0 };
            let rect = Rect::from_min_size(
                transform(Pos2::new(x, 20.0)),
                Vec2::new(30.0 * self.current_zoom, height * self.current_zoom)
            );
            
            let color = if chr == &"X" { Color32::from_rgb(255, 100, 150) } else { Color32::from_rgb(100, 150, 255) };
            painter.rect_filled(rect, 8.0, color);
            painter.rect_stroke(rect, 8.0, Stroke::new(2.0, Color32::BLACK));
            
            if self.current_zoom > 0.5 {
                painter.text(
                    rect.center() + Vec2::new(0.0, (height + 20.0) * self.current_zoom),
                    egui::Align2::CENTER_CENTER,
                    chr.to_string(),
                    egui::FontId::proportional(12.0),
                    Color32::WHITE,
                );
            }
        }
    }
    
    fn draw_roads_layer(&self, painter: &egui::Painter, transform: &dyn Fn(Pos2) -> Pos2) {
        // Draw k-mer roads with importance-based styling
        for edge in &self.edges {
            let from_pos = self.segment_to_hilbert_position(edge.from);
            let to_pos = self.segment_to_hilbert_position(edge.to);
            
            let start = transform(from_pos);
            let end = transform(to_pos);
            
            let importance = edge.importance as f32;
            let width = (importance.ln().max(1.0) * self.current_zoom * 0.5).min(8.0);
            
            // Color encoding: frequency -> color intensity
            let color = if importance > 50.0 {
                Color32::from_rgb(255, 50, 50)   // Very high frequency - bright red
            } else if importance > 20.0 {
                Color32::from_rgb(255, 100, 50)  // High frequency - orange-red
            } else if importance > 10.0 {
                Color32::from_rgb(255, 150, 50)  // Medium frequency - orange
            } else if importance > 5.0 {
                Color32::from_rgb(200, 200, 100) // Low frequency - yellow
            } else {
                Color32::from_rgb(150, 150, 150) // Very low frequency - gray
            };
            
            painter.line_segment([start, end], Stroke::new(width, color));
            
            // Draw direction arrows for high-importance roads
            if importance > 20.0 && self.current_zoom > 2.0 {
                let direction = (end - start).normalized();
                let mid_point = start + (end - start) * 0.5;
                let arrow_size = 6.0 * self.current_zoom;
                
                let perpendicular = Vec2::new(-direction.y, direction.x);
                let arrow_points = [
                    mid_point + direction * arrow_size,
                    mid_point - direction * arrow_size/2.0 + perpendicular * arrow_size/2.0,
                    mid_point - direction * arrow_size/2.0 - perpendicular * arrow_size/2.0,
                ];
                
                painter.add(egui::Shape::convex_polygon(
                    arrow_points.to_vec(),
                    Color32::from_rgba_unmultiplied(color.r(), color.g(), color.b(), 180),
                    Stroke::NONE
                ));
            }
        }
    }
    
    fn draw_genes_layer(&self, painter: &egui::Painter, transform: &dyn Fn(Pos2) -> Pos2) {
        // Only show genes at sufficient zoom
        if self.current_zoom < 3.0 { return; }
        
        // Draw gene representations
        for i in 0..30 {
            let x = (i % 10) as f32 * 80.0 + 50.0;
            let y = (i / 10) as f32 * 60.0 + 200.0;
            
            let start = transform(Pos2::new(x, y));
            let end = transform(Pos2::new(x + 60.0, y));
            
            // Gene body - thicker line
            painter.line_segment(
                [start, end],
                Stroke::new(4.0 * self.current_zoom, Color32::from_rgb(100, 220, 100))
            );
            
            // Exons as rectangles
            for j in 0..3 {
                let exon_start = start + Vec2::new(j as f32 * 20.0 * self.current_zoom, 0.0);
                let exon_rect = Rect::from_min_size(
                    exon_start - Vec2::new(0.0, 3.0 * self.current_zoom),
                    Vec2::new(15.0 * self.current_zoom, 6.0 * self.current_zoom)
                );
                painter.rect_filled(exon_rect, 2.0, Color32::from_rgb(80, 180, 80));
            }
            
            // Direction arrow
            let arrow_size = 8.0 * self.current_zoom;
            let arrow_points = [
                end,
                end + Vec2::new(-arrow_size, -arrow_size/3.0),
                end + Vec2::new(-arrow_size, arrow_size/3.0),
            ];
            painter.add(egui::Shape::convex_polygon(
                arrow_points.to_vec(),
                Color32::from_rgb(100, 220, 100),
                Stroke::NONE
            ));
            
            // Gene labels at high zoom
            if self.current_zoom > 6.0 {
                painter.text(
                    start + Vec2::new(0.0, -20.0 * self.current_zoom),
                    egui::Align2::LEFT_BOTTOM,
                    format!("GENE{}", i + 1),
                    egui::FontId::proportional(10.0),
                    Color32::WHITE,
                );
            }
        }
    }
    
    fn draw_variants_layer(&self, painter: &egui::Painter, transform: &dyn Fn(Pos2) -> Pos2) {
        if self.current_zoom < 5.0 { return; } // Only show variants at high zoom
        
        for (i, variant) in self.vcf_results.iter().enumerate() {
            let x = (i % 20) as f32 * 30.0 + 100.0;
            let y = (i / 20) as f32 * 20.0 + 350.0;
            
            let pos = transform(Pos2::new(x, y));
            let radius = 4.0 * self.current_zoom;
            
            let color = match variant.variant_type {
                VariantType::SNV => if variant.novel { Color32::YELLOW } else { Color32::BLUE },
                VariantType::Insertion => Color32::GREEN,
                VariantType::Deletion => Color32::RED,
                VariantType::StructuralVariant => Color32::from_rgb(128, 0, 128), // Purple
                VariantType::Complex => Color32::from_rgb(255, 165, 0), // Orange
            };
            
            painter.circle_filled(pos, radius, color);
            
            if variant.novel {
                painter.circle_stroke(pos, radius + 2.0, Stroke::new(2.0, Color32::YELLOW));
            }
            
            // Frequency indicator
            let freq_height = (variant.allele_frequencies[0] * 20.0 * self.current_zoom as f64) as f32;
            let freq_rect = Rect::from_min_size(
                pos + Vec2::new(radius + 2.0, -freq_height/2.0),
                Vec2::new(3.0, freq_height)
            );
            painter.rect_filled(freq_rect, 1.0, Color32::from_rgba_unmultiplied(255, 255, 255, 150));
        }
    }
    
    fn draw_annotations_layer(&self, painter: &egui::Painter, transform: &dyn Fn(Pos2) -> Pos2) {
        if self.current_zoom < 4.0 { return; }
        
        // Draw regulatory elements, promoters, enhancers, etc.
        for i in 0..15 {
            let x = i as f32 * 60.0 + 200.0;
            let y = 450.0;
            
            let pos = transform(Pos2::new(x, y));
            
            // Different shapes for different annotation types
            match i % 4 {
                0 => { // Promoter
                    let rect = Rect::from_center_size(pos, Vec2::splat(12.0 * self.current_zoom));
                    painter.rect_filled(rect, 2.0, Color32::from_rgb(255, 200, 100));
                    painter.rect_stroke(rect, 2.0, Stroke::new(1.0, Color32::BLACK));
                },
                1 => { // Enhancer
                    painter.circle_filled(pos, 6.0 * self.current_zoom, Color32::from_rgb(100, 255, 200));
                    painter.circle_stroke(pos, 6.0 * self.current_zoom, Stroke::new(1.0, Color32::BLACK));
                },
                2 => { // TFBS
                    let points = [
                        pos + Vec2::new(0.0, -8.0 * self.current_zoom),
                        pos + Vec2::new(6.0 * self.current_zoom, 4.0 * self.current_zoom),
                        pos + Vec2::new(-6.0 * self.current_zoom, 4.0 * self.current_zoom),
                    ];
                    painter.add(egui::Shape::convex_polygon(
                        points.to_vec(),
                        Color32::from_rgb(255, 100, 255),
                        Stroke::new(1.0, Color32::BLACK)
                    ));
                },
                _ => { // CpG Island
                    painter.circle_filled(pos, 4.0 * self.current_zoom, Color32::from_rgb(200, 200, 255));
                }
            }
        }
    }
    
    fn draw_info_overlay(&self, ui: &mut egui::Ui, response: &egui::Response) {
        if let Some(hover_pos) = response.hover_pos() {
            let center = response.rect.center();
            let world_pos = self.screen_to_world(hover_pos, center);
            
            // Show coordinates and any overlapping features
            let overlay_text = format!(
                "Hilbert: ({:.0}, {:.0})\nZoom: {:.2}x\nSegments: {}",
                world_pos.x, world_pos.y, self.current_zoom, self.edges.len()
            );
            
            ui.painter().text(
                hover_pos + Vec2::new(15.0, -40.0),
                egui::Align2::LEFT_TOP,
                overlay_text,
                egui::FontId::monospace(11.0),
                Color32::WHITE,
            );
        }
    }
    
    fn query_roads(&mut self) {
        let url = format!("{}/roads/query", self.server);
        let body = serde_json::json!({
            "segment_id": self.segment_id, 
            "max_depth": self.radius
        });
        
        match reqwest::blocking::Client::new().post(&url).json(&body).send() {
            Ok(res) => {
                match res.json::<serde_json::Value>() {
                    Ok(json) => {
                        if let Some(edges) = json["edges"].as_array() {
                            self.edges = edges
                                .iter()
                                .filter_map(|v| serde_json::from_value(v.clone()).ok())
                                .collect();
                        } else if let Some(neighbors) = json["neighbors"].as_array() {
                            self.edges = neighbors
                                .iter()
                                .filter_map(|v| v.as_u64())
                                .map(|neighbor| Edge {
                                    from: self.segment_id,
                                    to: neighbor,
                                    importance: 1,
                                })
                                .collect();
                        }
                    },
                    Err(e) => eprintln!("Failed to parse JSON: {}", e),
                }
            },
            Err(e) => eprintln!("Failed to query roads: {}", e),
        }
    }
    
    fn segment_to_hilbert_position(&self, segment_id: u64) -> Pos2 {
        if let Some(ref index) = self.index {
            if let Some(meta) = index.meta.get(&segment_id) {
                return Pos2::new(meta.hilbert_xy.0 as f32 * 3.0, meta.hilbert_xy.1 as f32 * 3.0);
            }
        }
        // Fallback calculation
        let x = (segment_id % 256) as f32 * 3.0;
        let y = ((segment_id / 256) % 256) as f32 * 3.0;
        Pos2::new(x, y)
    }
    
    fn chromosome_color(&self, chromosome: i32) -> Color32 {
        let colors = [
            Color32::from_rgb(255, 100, 100), Color32::from_rgb(100, 255, 100), Color32::from_rgb(100, 100, 255),
            Color32::from_rgb(255, 255, 100), Color32::from_rgb(255, 100, 255), Color32::from_rgb(100, 255, 255),
            Color32::from_rgb(255, 150, 100), Color32::from_rgb(150, 255, 100), Color32::from_rgb(100, 150, 255),
            Color32::from_rgb(255, 100, 150), Color32::from_rgb(150, 100, 255), Color32::from_rgb(100, 255, 150),
            Color32::from_rgb(200, 200, 100), Color32::from_rgb(200, 100, 200), Color32::from_rgb(100, 200, 200),
            Color32::from_rgb(255, 200, 150), Color32::from_rgb(200, 255, 150), Color32::from_rgb(150, 200, 255),
            Color32::from_rgb(255, 150, 200), Color32::from_rgb(200, 150, 255), Color32::from_rgb(150, 255, 200),
            Color32::from_rgb(180, 180, 180),
        ];
        colors.get((chromosome - 1) as usize).copied().unwrap_or(Color32::GRAY)
    }
    
    fn screen_to_world(&self, screen_pos: Pos2, center: Pos2) -> Pos2 {
        let relative = screen_pos - center - self.pan_offset;
        Pos2::new(relative.x / self.current_zoom, relative.y / self.current_zoom)
    }
    
    fn zoom_to_gene(&mut self) {
        if let Some(ref mut visualizer) = self.visualizer {
            visualizer.zoom_to_gene(&self.gene_search, &self.current_genome);
            self.current_zoom = 8.0;
        }
        println!("Zooming to gene: {}", self.gene_search);
    }
    
    fn zoom_to_chromosome(&mut self, chromosome: &str) {
        if let Some(ref mut visualizer) = self.visualizer {
            visualizer.zoom_to_chromosome(chromosome, &self.current_genome);
            self.current_zoom = 2.0;
        }
        
        let chr_num = chromosome.strip_prefix("chr").unwrap_or("1");
        if let Ok(num) = chr_num.parse::<i32>() {
            self.pan_offset = Vec2::new(-(num as f32) * 40.0, -80.0);
        } else if chr_num == "X" {
            self.pan_offset = Vec2::new(-22.0 * 40.0, -80.0);
        } else if chr_num == "Y" {
            self.pan_offset = Vec2::new(-23.0 * 40.0, -80.0);
        }
        
        println!("Zooming to chromosome: {}", chromosome);
    }
    
    fn zoom_to_genome(&mut self) {
        self.pan_offset = Vec2::ZERO;
        self.current_zoom = 1.0;
        println!("Zooming to whole genome view");
    }
    
    fn reset_view(&mut self) {
        self.pan_offset = Vec2::ZERO;
        self.current_zoom = 1.0;
        self.gene_search.clear();
    }
    
    fn export_sam_format(&self) {
        let mut sam_content = String::new();
        sam_content.push_str("@HD\tVN:1.6\tSO:coordinate\n");
        sam_content.push_str("@PG\tID:fractal_pangenome\tPN:FractalPangenome\tVN:1.0\n");
        
        if let Some(ref index) = self.index {
            for meta_entry in index.meta.iter() {
                let (segment_id, meta) = meta_entry.pair();
                sam_content.push_str(&format!(
                    "@SQ\tSN:segment_{}\tLN:{}\tGI:{}\tHC:{},{}\n",
                    segment_id, meta.length, meta.genome_id, meta.hilbert_xy.0, meta.hilbert_xy.1
                ));
            }
        }
        
        for alignment in &self.alignment_results {
            let flag = 0;
            let rname = format!("segment_{}", alignment.segment_path[0]);
            let pos = alignment.start_pos.unwrap_or(1);
            let mapq = alignment.mapping_quality;
            let cigar = &alignment.cigar;
            let rnext = "*";
            let pnext = 0;
            let tlen = 0;
            let seq = "*";
            let qual = "*";
            
            let tags = format!(
                "PG:Z:{}\tHC:Z:{:?}\tRI:Z:{}\tAT:Z:{:?}\tSP:Z:{:?}",
                alignment.genome_matches.join(","),
                alignment.hilbert_coordinates,
                alignment.road_importance.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(","),
                alignment.alignment_type,
                alignment.segment_path
            );
            
            sam_content.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                alignment.read_id, flag, rname, pos, mapq, cigar, 
                rnext, pnext, tlen, seq, qual, tags
            ));
        }
        
        println!("SAM format generated ({} alignments)", self.alignment_results.len());
        println!("Content preview:\n{}", sam_content.lines().take(10).collect::<Vec<_>>().join("\n"));
    }
    
    fn export_vcf_format(&self) {
        let mut vcf_content = String::new();
        vcf_content.push_str("##fileformat=VCFv4.3\n");
        vcf_content.push_str("##source=FractalPangenome\n");
        vcf_content.push_str("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency in Pangenome\">\n");
        vcf_content.push_str("##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele Count in Pangenome\">\n");
        vcf_content.push_str("##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total Number of Alleles in Pangenome\">\n");
        vcf_content.push_str("##INFO=<ID=NOVEL,Number=0,Type=Flag,Description=\"Novel variant not seen in reference\">\n");
        vcf_content.push_str("##INFO=<ID=SEGMENTS,Number=.,Type=String,Description=\"Supporting pangenome segments\">\n");
        vcf_content.push_str("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n");
        vcf_content.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
        
        for (i, variant) in self.vcf_results.iter().enumerate() {
            let id = format!("var_{}", i);
            let alt = variant.alternate_alleles.join(",");
            let qual = format!("{:.1}", variant.quality_score);
            let filter = "PASS";
            
            let mut info_fields = vec![
                format!("AF={:.6}", variant.allele_frequencies[0]),
                format!("AC={}", variant.allele_counts[0]),
                format!("AN={}", variant.total_genomes * 2),
            ];
            
            if variant.novel {
                info_fields.push("NOVEL".to_string());
            }
            
            info_fields.push(format!("SEGMENTS={}", 
                variant.supporting_segments.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")));
            
            info_fields.push(format!("SVTYPE={:?}", variant.variant_type));
            
            let info = info_fields.join(";");
            
            vcf_content.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                variant.chromosome, variant.position, id, variant.reference_allele,
                alt, qual, filter, info
            ));
        }
        
        println!("VCF format generated ({} variants)", self.vcf_results.len());
        println!("Content preview:\n{}", vcf_content.lines().take(15).collect::<Vec<_>>().join("\n"));
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> { 
    let opts = eframe::NativeOptions::default(); 
    eframe::run_native(
        "Fractal Pangenome Viewer - Complete", 
        opts, 
        Box::new(|_| Ok(Box::new(App::new())))
    )?;
    Ok(()) 
}

impl eframe::App for App {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // Menu bar
        egui::TopBottomPanel::top("menu_bar").show(ctx, |ui| {
            egui::menu::bar(ui, |ui| {
                ui.menu_button("Database", |ui| {
                    if ui.button("📂 Load Existing Database").clicked() {
                        self.load_database();
                        ui.close_menu();
                    }
                    
                    if ui.button("🆕 Create New Database").clicked() {
                        self.create_new_database = true;
                        self.show_file_browser = true;
                        ui.close_menu();
                    }
                    
                    ui.separator();
                    
                    if ui.button("💾 Save Database").clicked() {
                        ui.close_menu();
                    }
                    
                    if ui.button("📊 Database Statistics").clicked() {
                        ui.close_menu();
                    }
                });
                
                ui.menu_button("Analysis", |ui| {
                    if ui.button("🧬 Align Reads").clicked() {
                        self.show_alignment_panel = true;
                        ui.close_menu();
                    }
                    
                    if ui.button("🔬 Generate VCF").clicked() {
                        self.generate_vcf_from_pangenome();
                        ui.close_menu();
                    }
                    
                    if ui.button("📈 Population Analysis").clicked() {
                        ui.close_menu();
                    }
                });
                
                ui.menu_button("View", |ui| {
                    ui.checkbox(&mut self.show_genes, "Show Genes");
                    ui.checkbox(&mut self.show_karyotype, "Show Karyotype");
                    ui.checkbox(&mut self.show_roads, "Show K-mer Roads");
                    ui.checkbox(&mut self.show_variants, "Show Variants");
                    ui.checkbox(&mut self.show_annotations, "Show Annotations");
                });
            });
        });
        
        // Top control panel
        egui::TopBottomPanel::top("controls").show(ctx, |ui| {
            ui.horizontal(|ui| {
                if self.database_loaded {
                    ui.colored_label(Color32::GREEN, "✓ Database Loaded");
                } else {
                    ui.colored_label(Color32::RED, "⚠ No Database");
                    ui.text_edit_singleline(&mut self.database_path);
                    if ui.button("Load").clicked() {
                        self.load_database();
                    }
                }
                
                ui.separator();
                
                ui.label("Server:"); 
                ui.text_edit_singleline(&mut self.server); 
                if ui.button("Connect").clicked() { 
                    self.connected = true;
                }
                
                ui.separator();
                
                ui.label("Segment:"); 
                ui.add(egui::DragValue::new(&mut self.segment_id)); 
                if ui.button("Query").clicked() && self.connected {
                    self.query_roads();
                }
            });
            
            ui.horizontal(|ui| {
                ui.label("Genome:");
                ui.text_edit_singleline(&mut self.current_genome);
                
                ui.separator();
                
                ui.label("Chromosome:");
                egui::ComboBox::from_label("")
                    .selected_text(&self.chromosome_view)
                    .show_ui(ui, |ui| {
                        for i in 1..=22 {
                            ui.selectable_value(&mut self.chromosome_view, format!("chr{}", i), format!("chr{}", i));
                        }
                        ui.selectable_value(&mut self.chromosome_view, "chrX".to_string(), "chrX");
                        ui.selectable_value(&mut self.chromosome_view, "chrY".to_string(), "chrY");
                    });
                
                if ui.button("🔍 View").clicked() {
                    let chromosome_view = self.chromosome_view.clone();
                    self.zoom_to_chromosome(&chromosome_view);
                }
                
                ui.separator();
                
                ui.label("Zoom:");
                ui.add(egui::Slider::new(&mut self.current_zoom, 0.1..=10.0).logarithmic(true));
                
                ui.separator();
                
                ui.label("Gene:");
                ui.text_edit_singleline(&mut self.gene_search);
                if ui.button("🔍").clicked() && !self.gene_search.is_empty() {
                    self.zoom_to_gene();
                }
                
                if ui.button("Reset View").clicked() {
                    self.reset_view();
                }
            });
        });
        
        // File browser dialog
        if self.show_file_browser {
            egui::Window::new("Create New Database")
                .collapsible(false)
                .resizable(true)
                .show(ctx, |ui| {
                    ui.heading("Select Files for New Pangenome");
                    
                    ui.horizontal(|ui| {
                        ui.label("Genome ID:");
                        ui.text_edit_singleline(&mut self.new_genome_id);
                    });
                    
                    ui.separator();
                    
                    ui.horizontal(|ui| {
                        ui.label("FASTA file:");
                        ui.text_edit_singleline(&mut self.selected_fasta_path);
                        if ui.button("Browse").clicked() {
                            self.selected_fasta_path = "/path/to/genome.fasta".to_string();
                        }
                    });
                    
                    ui.horizontal(|ui| {
                        ui.label("GFF3 file (optional):");
                        ui.text_edit_singleline(&mut self.selected_gff_path);
                        if ui.button("Browse").clicked() {
                            self.selected_gff_path = "/path/to/annotations.gff3".to_string();
                        }
                    });
                    
                    ui.horizontal(|ui| {
                        ui.label("Additional annotations (optional):");
                        ui.text_edit_singleline(&mut self.selected_annotation_path);
                        if ui.button("Browse").clicked() {
                        }
                    });
                    
                    ui.separator();
                    
                    ui.horizontal(|ui| {
                        if ui.button("Create Database").clicked() {
                            self.create_new_database_from_files();
                            self.show_file_browser = false;
                        }
                        
                        if ui.button("Cancel").clicked() {
                            self.show_file_browser = false;
                        }
                    });
                });
        }
        
        // Alignment panel
        if self.show_alignment_panel {
            egui::Window::new("Read Alignment")
                .collapsible(true)
                .resizable(true)
                .show(ctx, |ui| {
                    ui.heading("Align Reads to Pangenome");
                    
                    ui.label("Paste FASTQ/FASTA reads:");
                    ui.text_edit_multiline(&mut self.reads_input);
                    
                    ui.horizontal(|ui| {
                        if ui.button("🧬 Align to Pangenome").clicked() {
                            self.align_reads_to_pangenome();
                        }
                        
                        if ui.button("💾 Export SAM").clicked() {
                            self.export_sam_format();
                        }
                        
                        if ui.button("❌ Close").clicked() {
                            self.show_alignment_panel = false;
                        }
                    });
                    
                    ui.separator();
                    
                    if !self.alignment_results.is_empty() {
                        ui.heading("Alignment Results");
                        
                        egui::ScrollArea::vertical().show(ui, |ui| {
                            for (i, alignment) in self.alignment_results.iter().enumerate() {
                                ui.group(|ui| {
                                    ui.horizontal(|ui| {
                                        ui.label(&alignment.read_id);
                                        ui.label(format!("Score: {}", alignment.alignment_score));
                                        ui.label(format!("Quality: {}", alignment.mapping_quality));
                                        ui.label(format!("Type: {:?}", alignment.alignment_type));
                                    });
                                    
                                    ui.horizontal(|ui| {
                                        ui.label("Segment Path:");
                                        ui.label(format!("{:?}", alignment.segment_path));
                                    });
                                    
                                    ui.horizontal(|ui| {
                                        ui.label("Hilbert Coords:");
                                        ui.label(format!("{:?}", alignment.hilbert_coordinates));
                                    });
                                    
                                    ui.horizontal(|ui| {
                                        ui.label("Genomes:");
                                        ui.label(alignment.genome_matches.join(", "));
                                    });
                                });
                                
                                if i < self.alignment_results.len() - 1 {
                                    ui.separator();
                                }
                            }
                        });
                    }
                });
        }
        
        // VCF panel
        if self.show_vcf_panel {
            egui::Window::new("Variant Calls")
                .collapsible(true)
                .resizable(true)
                .show(ctx, |ui| {
                    ui.heading("Pangenome Variant Analysis");
                    
                    ui.horizontal(|ui| {
                        if ui.button("🔄 Refresh").clicked() {
                            self.generate_vcf_from_pangenome();
                        }
                        
                        if ui.button("💾 Export VCF").clicked() {
                            self.export_vcf_format();
                        }
                        
                        if ui.button("❌ Close").clicked() {
                            self.show_vcf_panel = false;
                        }
                    });
                    
                    ui.separator();
                    
                    if !self.vcf_results.is_empty() {
                        egui::ScrollArea::vertical().show(ui, |ui| {
                            for variant in &self.vcf_results {
                                ui.group(|ui| {
                                    ui.horizontal(|ui| {
                                        ui.label(&variant.chromosome);
                                        ui.label(format!("Pos: {}", variant.position));
                                        ui.label(format!("{}→{}", variant.reference_allele, variant.alternate_alleles.join(",")));
                                        ui.label(format!("AF: {:.3}", variant.allele_frequencies[0]));
                                        
                                        if variant.novel {
                                            ui.colored_label(Color32::YELLOW, "NOVEL");
                                        }
                                        
                                        ui.label(format!("Type: {:?}", variant.variant_type));
                                    });
                                    
                                    ui.horizontal(|ui| {
                                        ui.label(format!("Count: {}/{}", variant.allele_counts[0], variant.total_genomes));
                                        ui.label(format!("Quality: {:.1}", variant.quality_score));
                                        ui.label(format!("Segments: {:?}", variant.supporting_segments));
                                    });
                                });
                            }
                        });
                    }
                });
        }
        
        // Left panel with layers and navigation
        egui::SidePanel::left("layers").default_width(250.0).show(ctx, |ui| {
            ui.heading("Visualization Control");
            
            ui.separator();
            
            ui.group(|ui| {
                ui.label("Display Layers:");
                ui.checkbox(&mut self.show_karyotype, "🧬 Karyotype");
                ui.checkbox(&mut self.show_genes, "🧪 Genes");
                ui.checkbox(&mut self.show_roads, "🛣️ K-mer Roads");
                ui.checkbox(&mut self.show_variants, "🔬 Variants");
                ui.checkbox(&mut self.show_annotations, "📝 Annotations");
            });
            
            ui.separator();
            
            ui.group(|ui| {
                ui.label("Quick Navigation:");
                
                if ui.button("🌍 Whole Genome View").clicked() {
                    self.zoom_to_genome();
                }
                
                ui.horizontal(|ui| {
                    if ui.button("Chr 1").clicked() { self.zoom_to_chromosome("chr1"); }
                    if ui.button("Chr 2").clicked() { self.zoom_to_chromosome("chr2"); }
                    if ui.button("Chr X").clicked() { self.zoom_to_chromosome("chrX"); }
                });
                
                if ui.button("🎯 Find BRCA1").clicked() {
                    self.gene_search = "BRCA1".to_string();
                    self.zoom_to_gene();
                }
                
                if ui.button("🎯 Find TP53").clicked() {
                    self.gene_search = "TP53".to_string();
                    self.zoom_to_gene();
                }
            });
            
            ui.separator();
            
            ui.group(|ui| {
                ui.label("System Status:");
                
                if self.database_loaded {
                    ui.colored_label(Color32::GREEN, "✓ Database Loaded");
                } else {
                    ui.colored_label(Color32::RED, "⚠ No Database");
                }
                
                if self.connected {
                    ui.colored_label(Color32::GREEN, "✓ Server Connected");
                } else {
                    ui.colored_label(Color32::YELLOW, "⚠ Server Disconnected");
                }
                
                ui.label(format!("Zoom: {:.2}x", self.current_zoom));
                ui.label(format!("Pan: ({:.0}, {:.0})", self.pan_offset.x, self.pan_offset.y));
                ui.label(format!("K-mer Roads: {}", self.edges.len()));
                ui.label(format!("Alignments: {}", self.alignment_results.len()));
                ui.label(format!("Variants: {}", self.vcf_results.len()));
            });
        });
        
        // Main visualization panel
        egui::CentralPanel::default().show(ctx, |ui| {
            let (response, painter) = ui.allocate_painter(ui.available_size(), egui::Sense::click_and_drag());
            
            self.handle_viewport_interactions(&response, ctx);
            
            let center = response.rect.center();
            let transform = |pos: Pos2| -> Pos2 {
                let scaled = pos.to_vec2() * self.current_zoom;
                center + scaled + self.pan_offset
            };
            
            self.draw_all_layers(&painter, &response.rect, &transform);
            
            self.draw_info_overlay(ui, &response);
        });
    }
}
