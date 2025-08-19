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

pub use hilbert::HilbertMapper;
pub use enhanced_annotations::{EnhancedAnnotationStore, EnhancedAnnotationTag, FeatureType};
pub use visualization::{MultiScaleVisualizer, VisualizationLayer, LayerType};
pub use index::HierarchicalIndex;
pub use align::{align_reads_enhanced, EnhancedAlignHit, AlignmentType};
pub use io::load_index_hdf5;
