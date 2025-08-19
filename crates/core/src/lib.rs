pub mod hilbert;
pub mod encoding;
pub mod index;
pub mod io;
pub mod model;
pub mod gtf;
pub mod annotations; 
pub mod align; 
pub mod graph;

pub mod enhanced_annotations; // The enhanced annotation system
pub mod visualization; // The visualization system

pub use align::{AlignHit, align_reads_with_cigar_splice, align_reads_enhanced, EnhancedAlignHit};
pub use enhanced_annotations::{EnhancedAnnotationStore, EnhancedAnnotationTag, FeatureType, KaryotypeInfo};
pub use visualization::{MultiScaleVisualizer, VisualizationLayer, LayerType};
