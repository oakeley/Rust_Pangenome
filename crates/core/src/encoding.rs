use smallvec::SmallVec;

#[inline]
pub fn encode_kmers(seq: &[u8], k: usize) -> SmallVec<[u32; 2048]> {
    use smallvec::smallvec;
    if seq.len() < k { return smallvec![]; }
    let mut out: SmallVec<[u32; 2048]> = SmallVec::with_capacity(seq.len()-k+1);
    let mut val: u32 = 0;
    let mask: u32 = (1u32 << (2*k as u32)) - 1;
    let to2 = |c: u8| -> u32 { match c { b'A'|b'a'=>0, b'C'|b'c'=>1, b'G'|b'g'=>2, b'T'|b't'|b'U'|b'u'=>3, _=>0 } };
    for &c in &seq[..k] { val = (val<<2) | to2(c); }
    out.push(val);
    for &c in &seq[k..] { val = ((val<<2)|to2(c)) & mask; out.push(val); }
    out
}
pub fn reconstruct_seq_from_kmers(kmers: &[u32], k: usize) -> Vec<u8> {
    if kmers.is_empty() {
        return Vec::new();
    }
    
    // Simple reconstruction: take first k-mer and extend by one nucleotide per subsequent k-mer
    let mut result = Vec::new();
    
    // Decode first k-mer completely
    let first_kmer = kmers[0];
    for i in (0..k).rev() {
        let nt = (first_kmer >> (2 * i)) & 3;
        result.push(match nt {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            3 => b'T',
            _ => b'N',
        });
    }
    
    // For subsequent k-mers, just add the last nucleotide
    for &kmer in &kmers[1..] {
        let last_nt = kmer & 3;
        result.push(match last_nt {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            3 => b'T',
            _ => b'N',
        });
    }
    
    result
}
