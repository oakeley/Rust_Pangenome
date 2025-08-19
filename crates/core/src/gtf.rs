use serde::{Serialize, Deserialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GtfRecord {
    pub seqname: String,
    pub source: String,
    pub feature: String,
    pub start: u64,
    pub end: u64,
    pub score: Option<f32>,
    pub strand: Option<char>,
    pub frame: Option<u8>,
    pub attributes: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Exon {
    pub gene_id: String,
    pub transcript_id: Option<String>,
    pub start: u64,
    pub end: u64,
    pub strand: Option<char>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Gene {
    pub gene_id: String,
    pub gene_name: Option<String>,
    pub seqname: String,
    pub exons: Vec<Exon>,
}

pub fn parse_gtf(path: &str) -> anyhow::Result<Vec<GtfRecord>> {
    use std::io::BufRead;
    let f = std::fs::File::open(path)?;
    let r = std::io::BufReader::new(f);
    let mut out = Vec::new();
    for line in r.lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') { continue; }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 9 { continue; }
        let score = if parts[5] == "." { None } else { parts[5].parse().ok() };
        let strand = parts[6].chars().next();
        let frame = if parts[7] == "." { None } else { parts[7].parse().ok() };
        out.push(GtfRecord {
            seqname: parts[0].to_string(), source: parts[1].to_string(), feature: parts[2].to_string(),
            start: parts[3].parse().unwrap_or(0), end: parts[4].parse().unwrap_or(0),
            score, strand, frame, attributes: parts[8].to_string()
        });
    }
    Ok(out)
}

pub fn build_genes(records: &[GtfRecord]) -> Vec<Gene> {
    use std::collections::HashMap;
    let mut genes: HashMap<String, Gene> = HashMap::new();
    for rec in records {
        if rec.feature != "exon" { continue; }
        let mut gene_id = String::new();
        let mut gene_name = None;
        let mut transcript_id = None;
        for kv in rec.attributes.split(';') {
            let kv = kv.trim();
            if kv.is_empty() { continue; }
            if let Some((k,v)) = kv.split_once(' ') {
                let v = v.trim_matches('"');
                match k {
                    "gene_id" => gene_id = v.to_string(),
                    "gene_name" => gene_name = Some(v.to_string()),
                    "transcript_id" => transcript_id = Some(v.to_string()),
                    _ => {}
                }
            }
        }
        if gene_id.is_empty() { continue; }
        let e = Exon{ gene_id: gene_id.clone(), transcript_id, start: rec.start, end: rec.end, strand: rec.strand };
        let g = genes.entry(gene_id.clone()).or_insert(Gene{ gene_id: gene_id.clone(), gene_name: gene_name.clone(), seqname: rec.seqname.clone(), exons: vec![] });
        g.exons.push(e);
    }
    genes.into_values().collect()
}
