use axum::{
    extract::{Query, State},  // Add State here
    response::Json,
    routing::{get, post},     // Add post here
    Router,
};
use serde::Deserialize;
use std::collections::HashMap;
use std::sync::Arc;           // Add Arc import
use std::path::PathBuf;       // Add PathBuf import
use tower_http::cors::{CorsLayer, Any};  // Add Any import
use core::encoding::encode_kmers;  // Keep only one import of encode_kmers

// Remove duplicate imports and unused imports
// use core::encoding::encode_kmers;  // Remove this duplicate line

// Add your other necessary imports here based on what you're using in the file
// For example, if you're using these types, add them:
use core::index::HierarchicalIndex;
use core::hilbert::HilbertMapper;
// use core::annotations::AnnotationStore;  // Add if this exists
// use core::align::align_reads_with_cigar_splice;  // Add if this exists

#[derive(Deserialize)]
struct LoadReq {
    gtf_path: String,
    genome_id: String,
}

#[derive(Deserialize)]
struct TileQuery {
    z: u16,
    x: u16,
    y: u16,
    genome_id: String,
}

#[derive(Deserialize)]
struct PortReq {
    // Add fields as needed
}

#[derive(Deserialize)]
struct GeneQuery {
    gene_name: String,
    genome_id: String,
}

#[derive(Deserialize)]
struct AddGenomeReq {
    fasta_path: String,
    genome_id: String,
    k: u8,
}

#[derive(Deserialize)]
struct AlignReq {
    reads: String,
    k: usize,
    top_n: usize,
}

#[derive(Deserialize)]
struct RoadsQuery {
    segment_id: u64,
    max_depth: usize,
}

// You'll need to define or import AnnotationStore
#[derive(Default)]
struct AnnotationStore {
    // Add fields as needed
}

#[derive(Clone)]
struct AppState { 
    index: Arc<HierarchicalIndex>, 
    ann: Arc<AnnotationStore> 
}

#[tokio::main]
async fn main() {
    env_logger::init();
    let index = Arc::new(HierarchicalIndex::new(HilbertMapper::new(256), &[1,4,16,64]));
    
    let ann = Arc::new(AnnotationStore::default());
    
    let app = Router::new()
        .route("/", get(root))
        .route("/stats", get(stats))
        .route("/annotations/load", post(load_gtf))
        .route("/annotations/by_tile", get(ann_by_tile))
        .route("/annotations/port", post(port_annotations))
        .route("/annotations/search_gene", post(search_gene))
        .route("/genomes/add", post(add_genome))
        .route("/align/reads", post(align_reads))
        .route("/roads/query", post(query_roads))
        .with_state(AppState { index, ann })
        .layer(CorsLayer::new().allow_methods(Any).allow_origin(Any).allow_headers(Any));

    let listener = tokio::net::TcpListener::bind("0.0.0.0:8787").await.unwrap();
    println!("Server running on http://0.0.0.0:8787");
    axum::serve(listener, app).await.unwrap();
}

async fn root() -> &'static str {
    "Fractal Pangenome Server"
}

async fn stats(State(st): State<AppState>) -> Json<serde_json::Value> {
    let total_segs = st.index.meta.len();
    let total_genomes = st.index.road_graph.genome_paths.len();
    Json(serde_json::json!({
        "segments": total_segs,
        "genomes": total_genomes
    }))
}

async fn load_gtf(State(st): State<AppState>, Json(req): Json<LoadReq>) -> Json<serde_json::Value> {
    // Implementation depends on your GTF loading logic
    // This is a placeholder - you'll need to implement the actual GTF loading
    Json(serde_json::json!({
        "status": "success",
        "message": format!("Would load GTF from {} for genome {}", req.gtf_path, req.genome_id)
    }))
}

async fn ann_by_tile(State(st): State<AppState>, Query(q): Query<TileQuery>) -> Json<serde_json::Value> {
    // Implementation depends on your annotation logic
    // This is a placeholder
    Json(serde_json::json!({
        "items": [],
        "tile": format!("{}/{}/{}", q.z, q.x, q.y),
        "genome_id": q.genome_id
    }))
}

async fn port_annotations(State(st): State<AppState>, Json(req): Json<PortReq>) -> Json<serde_json::Value> {
    // Implementation depends on your annotation porting logic
    Json(serde_json::json!({
        "status": "success",
        "message": "Annotations ported"
    }))
}

async fn search_gene(State(st): State<AppState>, Json(q): Json<GeneQuery>) -> Json<serde_json::Value> {
    // Implementation depends on your gene search logic
    Json(serde_json::json!({
        "gene": q.gene_name,
        "genome_id": q.genome_id,
        "results": []
    }))
}

async fn add_genome(State(st): State<AppState>, Json(req): Json<AddGenomeReq>) -> Json<serde_json::Value> {
    let p = PathBuf::from(&req.fasta_path);
    if !p.exists() {
        return Json(serde_json::json!({
            "error": "FASTA file not found"
        }));
    }
    
    // Read FASTA and process
    match std::fs::read_to_string(&p) {
        Ok(content) => {
            let mut seq = String::new();
            for line in content.lines() {
                if !line.starts_with('>') {
                    seq.push_str(line.trim());
                }
            }
            
            if !seq.is_empty() {
                let k = req.k as usize;
                let kmers = encode_kmers(seq.as_bytes(), k);
                // Process kmers and add to index
                // This is a simplified version - you may need more logic here
                
                Json(serde_json::json!({
                    "status": "success",
                    "genome_id": req.genome_id,
                    "sequence_length": seq.len(),
                    "kmers_count": kmers.len()
                }))
            } else {
                Json(serde_json::json!({
                    "error": "Empty sequence"
                }))
            }
        },
        Err(e) => Json(serde_json::json!({
            "error": format!("Failed to read file: {}", e)
        }))
    }
}

async fn align_reads(State(st): State<AppState>, Json(req): Json<AlignReq>) -> Json<serde_json::Value> {
    // Implementation depends on your alignment logic
    // You may need to import and use align_reads_with_cigar_splice
    Json(serde_json::json!({
        "alignments": [],
        "reads_processed": req.reads.lines().count()
    }))
}

async fn query_roads(State(st): State<AppState>, Json(q): Json<RoadsQuery>) -> Json<serde_json::Value> {
    // Implementation depends on your road graph query logic
    let neighbors = st.index.road_graph.neighbors(q.segment_id);
    Json(serde_json::json!({
        "segment_id": q.segment_id,
        "neighbors": neighbors
    }))
}
