use eframe::egui::{self, Color32, Pos2, Stroke, Rect, Vec2};
use anyhow::Result;
use serde::Deserialize;
use std::collections::HashMap;

#[derive(Default)]
struct App { 
    server: String, 
    connected: bool, 
    segment_id: u64, 
    radius: usize, 
    edges: Vec<Edge>,
    
    // New visualization fields
    visualizer: Option<MultiScaleVisualizer>,
    current_zoom: f32,
    pan_offset: Vec2,
    selected_layers: HashMap<String, bool>,
    show_genes: bool,
    show_karyotype: bool,
    show_roads: bool,
    gene_search: String,
    current_genome: String,
}

#[derive(Deserialize, Clone)]
struct Edge { 
    from: u64, 
    to: u64, 
    importance: u64 
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
            gene_search: String::new(),
            current_genome: "Human_T2T_CHM13".to_string(),
        };
        
        // Initialize layer selections
        app.selected_layers.insert("Karyotype".to_string(), true);
        app.selected_layers.insert("Genes".to_string(), true);
        app.selected_layers.insert("Roads".to_string(), false);
        app.selected_layers.insert("Variants".to_string(), false);
        
        app
    }
    
    fn initialize_visualizer(&mut self) {
        // This would be called when you load a database
        // For now, we'll create a placeholder
        
        /* Uncomment when you have the modules available:
        let hilbert = HilbertMapper::new(256);
        let annotation_store = EnhancedAnnotationStore::default();
        self.visualizer = Some(MultiScaleVisualizer::new(hilbert, annotation_store));
        */
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> { 
    let opts = eframe::NativeOptions::default(); 
    eframe::run_native(
        "Fractal Pangenome Viewer", 
        opts, 
        Box::new(|_| Ok(Box::new(App::new())))
    )?;
    Ok(()) 
}

impl eframe::App for App {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // Top panel with controls
        egui::TopBottomPanel::top("top").show(ctx, |ui| {
            ui.horizontal(|ui| {
                // Server connection
                ui.label("Server:"); 
                ui.text_edit_singleline(&mut self.server); 
                if ui.button("Connect").clicked() { 
                    self.connected = true;
                    self.initialize_visualizer();
                }
                
                ui.separator();
                
                // Legacy road query controls
                ui.label("Segment ID:"); 
                ui.add(egui::DragValue::new(&mut self.segment_id)); 
                ui.label("Radius:"); 
                ui.add(egui::DragValue::new(&mut self.radius));
                
                if ui.button("Query roads").clicked() && self.connected {
                    self.query_roads();
                }
            });
            
            // Second row with visualization controls
            ui.horizontal(|ui| {
                ui.label("Genome:");
                ui.text_edit_singleline(&mut self.current_genome);
                
                ui.separator();
                
                ui.label("Zoom:");
                ui.add(egui::Slider::new(&mut self.current_zoom, 0.1..=10.0).logarithmic(true));
                
                ui.separator();
                
                ui.label("Search Gene:");
                ui.text_edit_singleline(&mut self.gene_search);
                if ui.button("🔍").clicked() && !self.gene_search.is_empty() {
                    self.zoom_to_gene();
                }
                
                ui.separator();
                
                if ui.button("Reset View").clicked() {
                    self.reset_view();
                }
            });
        });
        
        // Left panel with layer controls
        egui::SidePanel::left("layers").default_width(200.0).show(ctx, |ui| {
            ui.heading("Visualization Layers");
            
            ui.separator();
            
            // Layer toggles
            ui.checkbox(&mut self.show_karyotype, "🧬 Karyotype");
            ui.checkbox(&mut self.show_genes, "🧪 Genes");
            ui.checkbox(&mut self.show_roads, "🛣️ K-mer Roads");
            
            ui.separator();
            
            ui.heading("Navigation");
            
            if ui.button("🌍 Whole Genome").clicked() {
                self.zoom_to_genome();
            }
            
            if ui.button("🧬 Chromosome 1").clicked() {
                self.zoom_to_chromosome("chr1");
            }
            
            if ui.button("🧬 Chromosome X").clicked() {
                self.zoom_to_chromosome("chrX");
            }
            
            ui.separator();
            
            ui.heading("Information");
            ui.label(format!("Zoom: {:.2}x", self.current_zoom));
            ui.label(format!("Pan: ({:.0}, {:.0})", self.pan_offset.x, self.pan_offset.y));
            ui.label(format!("Edges: {}", self.edges.len()));
            
            if self.connected {
                ui.colored_label(Color32::GREEN, "✓ Connected to server");
            } else {
                ui.colored_label(Color32::RED, "⚠ Not connected to server");
            }
        });
        
        // Main visualization panel
        egui::CentralPanel::default().show(ctx, |ui| {
            let (response, painter) = ui.allocate_painter(ui.available_size(), egui::Sense::click_and_drag());
            
            // Handle pan and zoom
            if response.dragged() {
                self.pan_offset += response.drag_delta() / self.current_zoom;
            }
            
            if let Some(hover_pos) = response.hover_pos() {
                let scroll_delta = ctx.input(|i| i.scroll_delta.y);
                if scroll_delta != 0.0 {
                    let zoom_factor = 1.0 + scroll_delta * 0.001;
                    self.current_zoom *= zoom_factor;
                    self.current_zoom = self.current_zoom.clamp(0.1, 100.0);
                }
            }
            
            // Calculate view transformation
            let center = response.rect.center();
            let transform = |pos: Pos2| -> Pos2 {
                let scaled = pos.to_vec2() * self.current_zoom;
                center + scaled + self.pan_offset
            };
            
            // Draw background grid
            self.draw_background_grid(&painter, &response.rect);
            
            // Draw karyotype if enabled
            if self.show_karyotype {
                self.draw_karyotype(&painter, &transform);
            }
            
            // Draw genes if enabled  
            if self.show_genes {
                self.draw_genes(&painter, &transform);
            }
            
            // Draw k-mer roads if enabled
            if self.show_roads {
                self.draw_roads(&painter, &transform);
            }
            
            // Draw legacy edges for compatibility
            self.draw_legacy_edges(&painter, &transform);
            
            // Draw coordinate info
            if let Some(hover_pos) = response.hover_pos() {
                let world_pos = self.screen_to_world(hover_pos, center);
                ui.painter().text(
                    hover_pos + Vec2::new(10.0, -20.0),
                    egui::Align2::LEFT_TOP,
                    format!("Hilbert: ({:.0}, {:.0})", world_pos.x, world_pos.y),
                    egui::FontId::monospace(12.0),
                    Color32::WHITE,
                );
            }
        });
    }
}

impl App {
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
    
    fn draw_background_grid(&self, painter: &egui::Painter, rect: &Rect) {
        let grid_size = 20.0 * self.current_zoom;
        let start_x = (rect.min.x / grid_size).floor() * grid_size;
        let start_y = (rect.min.y / grid_size).floor() * grid_size;
        
        let mut x = start_x;
        while x <= rect.max.x {
            painter.line_segment(
                [Pos2::new(x, rect.min.y), Pos2::new(x, rect.max.y)],
                Stroke::new(0.5, Color32::from_gray(40))
            );
            x += grid_size;
        }
        
        let mut y = start_y;
        while y <= rect.max.y {
            painter.line_segment(
                [Pos2::new(rect.min.x, y), Pos2::new(rect.max.x, y)],
                Stroke::new(0.5, Color32::from_gray(40))
            );
            y += grid_size;
        }
    }
    
    fn draw_karyotype(&self, painter: &egui::Painter, transform: &dyn Fn(Pos2) -> Pos2) {
        // Draw chromosome representations as colored rectangles
        for i in 1..=22 {
            let x = (i as f32 - 1.0) * 30.0;
            let rect = Rect::from_min_size(
                transform(Pos2::new(x, 10.0)),
                Vec2::new(20.0 * self.current_zoom, 100.0 * self.current_zoom)
            );
            
            let color = self.chromosome_color(i);
            painter.rect_filled(rect, 5.0, color);
            painter.rect_stroke(rect, 5.0, Stroke::new(1.0, Color32::BLACK));
            
            // Label
            painter.text(
                rect.center() + Vec2::new(0.0, 60.0 * self.current_zoom),
                egui::Align2::CENTER_CENTER,
                format!("{}", i),
                egui::FontId::proportional(10.0),
                Color32::WHITE,
            );
        }
        
        // X and Y chromosomes
        for (i, chr) in ["X", "Y"].iter().enumerate() {
            let x = 22.0 * 30.0 + (i as f32) * 30.0;
            let rect = Rect::from_min_size(
                transform(Pos2::new(x, 10.0)),
                Vec2::new(20.0 * self.current_zoom, 80.0 * self.current_zoom)
            );
            
            let color = if chr == &"X" { Color32::from_rgb(255, 100, 150) } else { Color32::from_rgb(100, 150, 255) };
            painter.rect_filled(rect, 5.0, color);
            painter.rect_stroke(rect, 5.0, Stroke::new(1.0, Color32::BLACK));
            
            painter.text(
                rect.center() + Vec2::new(0.0, 50.0 * self.current_zoom),
                egui::Align2::CENTER_CENTER,
                chr.to_string(),
                egui::FontId::proportional(10.0),
                Color32::WHITE,
            );
        }
    }
    
    fn draw_genes(&self, painter: &egui::Painter, transform: &dyn Fn(Pos2) -> Pos2) {
        // Draw gene representations as arrows
        // This is a placeholder - in a real implementation, you'd query your annotation store
        
        if self.current_zoom > 2.0 { // Only show genes at sufficient zoom
            for i in 0..20 {
                let x = i as f32 * 50.0 + 100.0;
                let y = 200.0 + (i % 3) as f32 * 30.0;
                
                let start = transform(Pos2::new(x, y));
                let end = transform(Pos2::new(x + 40.0, y));
                
                // Gene body
                painter.line_segment(
                    [start, end],
                    Stroke::new(3.0, Color32::from_rgb(100, 200, 100))
                );
                
                // Arrow head for direction
                let arrow_size = 5.0 * self.current_zoom;
                let arrow_points = [
                    end,
                    end + Vec2::new(-arrow_size, -arrow_size/2.0),
                    end + Vec2::new(-arrow_size, arrow_size/2.0),
                ];
                painter.add(egui::Shape::convex_polygon(
                    arrow_points.to_vec(),
                    Color32::from_rgb(100, 200, 100),
                    Stroke::NONE
                ));
                
                // Gene label
                if self.current_zoom > 5.0 {
                    painter.text(
                        start + Vec2::new(0.0, -15.0),
                        egui::Align2::LEFT_BOTTOM,
                        format!("GENE{}", i + 1),
                        egui::FontId::proportional(9.0),
                        Color32::WHITE,
                    );
                }
            }
        }
    }
    
    fn draw_roads(&self, painter: &egui::Painter, transform: &dyn Fn(Pos2) -> Pos2) {
        // Draw k-mer roads (connections between segments)
        for edge in &self.edges {
            let from_pos = self.segment_to_position(edge.from);
            let to_pos = self.segment_to_position(edge.to);
            
            let start = transform(from_pos);
            let end = transform(to_pos);
            
            // Line width based on importance
            let importance = edge.importance as f32;
            let width = (importance.ln().max(1.0) * self.current_zoom).min(10.0);
            
            // Color based on importance
            let color = if importance > 10 {
                Color32::from_rgb(255, 100, 100) // High importance = red
            } else if importance > 5 {
                Color32::from_rgb(255, 200, 100) // Medium importance = orange
            } else {
                Color32::from_rgb(150, 150, 150) // Low importance = gray
            };
            
            painter.line_segment([start, end], Stroke::new(width, color));
            
            // Draw directional arrow for high importance roads
            if importance > 10 && self.current_zoom > 3.0 {
                let direction = (end - start).normalized();
                let mid_point = start + (end - start) * 0.5;
                let arrow_size = 8.0 * self.current_zoom;
                
                let perpendicular = Vec2::new(-direction.y, direction.x);
                let arrow_points = [
                    mid_point + direction * arrow_size,
                    mid_point - direction * arrow_size/2.0 + perpendicular * arrow_size/2.0,
                    mid_point - direction * arrow_size/2.0 - perpendicular * arrow_size/2.0,
                ];
                
                painter.add(egui::Shape::convex_polygon(
                    arrow_points.to_vec(),
                    color,
                    Stroke::NONE
                ));
            }
        }
    }
    
    fn draw_legacy_edges(&self, painter: &egui::Painter, transform: &dyn Fn(Pos2) -> Pos2) {
        // Keep your original edge drawing for compatibility
        for edge in &self.edges {
            let fx = (edge.from % 256) as f32 * 2.0; 
            let fy = ((edge.from / 256) % 256) as f32 * 2.0;
            let tx = (edge.to % 256) as f32 * 2.0; 
            let ty = ((edge.to / 256) % 256) as f32 * 2.0;
            
            let start = transform(Pos2::new(fx, fy));
            let end = transform(Pos2::new(tx, ty));
            
            let w = (edge.importance as f32).ln().max(1.0) * self.current_zoom;
            
            painter.line_segment([start, end], Stroke { 
                width: w, 
                color: Color32::from_rgb(200, 100, 100) 
            });
        }
    }
    
    fn segment_to_position(&self, segment_id: u64) -> Pos2 {
        // Convert segment ID to Hilbert coordinates
        // This is a placeholder - you'd use your HilbertMapper here
        let x = (segment_id % 256) as f32 * 3.0;
        let y = ((segment_id / 256) % 256) as f32 * 3.0;
        Pos2::new(x, y)
    }
    
    fn chromosome_color(&self, chromosome: i32) -> Color32 {
        // Generate distinct colors for chromosomes
        let colors = [
            Color32::from_rgb(255, 100, 100), Color32::from_rgb(100, 255, 100), Color32::from_rgb(100, 100, 255),
            Color32::from_rgb(255, 255, 100), Color32::from_rgb(255, 100, 255), Color32::from_rgb(100, 255, 255),
            Color32::from_rgb(255, 150, 100), Color32::from_rgb(150, 255, 100), Color32::from_rgb(100, 150, 255),
            Color32::from_rgb(255, 100, 150), Color32::from_rgb(150, 100, 255), Color32::from_rgb(100, 255, 150),
            Color32::from_rgb(200, 200, 100), Color32::from_rgb(200, 100, 200), Color32::from_rgb(100, 200, 200),
            Color32::from_rgb(255, 200, 150), Color32::from_rgb(200, 255, 150), Color32::from_rgb(150, 200, 255),
            Color32::from_rgb(255, 150, 200), Color32::from_rgb(200, 150, 255), Color32::from_rgb(150, 255, 200),
            Color32::from_rgb(180, 180, 180), // Chromosome 22
        ];
        
        colors.get((chromosome - 1) as usize).copied()
            .unwrap_or(Color32::from_rgb(128, 128, 128))
    }
    
    fn screen_to_world(&self, screen_pos: Pos2, center: Pos2) -> Pos2 {
        let relative = screen_pos - center - self.pan_offset;
        Pos2::new(relative.x / self.current_zoom, relative.y / self.current_zoom)
    }
    
    fn zoom_to_gene(&mut self) {
        // In a real implementation, this would query your annotation store
        // For now, just zoom to a placeholder position
        self.pan_offset = Vec2::new(-200.0, -200.0);
        self.current_zoom = 8.0;
        println!("Zooming to gene: {}", self.gene_search);
    }
    
    fn zoom_to_chromosome(&mut self, chromosome: &str) {
        // Zoom to show entire chromosome
        let chr_num = chromosome.strip_prefix("chr").unwrap_or("1");
        if let Ok(num) = chr_num.parse::<i32>() {
            self.pan_offset = Vec2::new(-(num as f32) * 30.0, -50.0);
            self.current_zoom = 3.0;
        } else if chr_num == "X" {
            self.pan_offset = Vec2::new(-22.0 * 30.0, -50.0);
            self.current_zoom = 3.0;
        } else if chr_num == "Y" {
            self.pan_offset = Vec2::new(-23.0 * 30.0, -50.0);
            self.current_zoom = 3.0;
        }
        println!("Zooming to chromosome: {}", chromosome);
    }
    
    fn zoom_to_genome(&mut self) {
        // Reset to show entire genome
        self.pan_offset = Vec2::ZERO;
        self.current_zoom = 1.0;
        println!("Zooming to whole genome view");
    }
    
    fn reset_view(&mut self) {
        self.pan_offset = Vec2::ZERO;
        self.current_zoom = 1.0;
        self.gene_search.clear();
    }
}
