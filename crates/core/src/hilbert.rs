#[derive(Debug, Clone)]
pub struct HilbertMapper {
    pub size: u16,
    pub total_positions: u32,
}

impl HilbertMapper {
    pub fn new(size: u16) -> Self {
        Self { size, total_positions: (size as u32) * (size as u32) }
    }
    pub fn index_to_xy(&self, mut idx: u32) -> (u16,u16) {
        let mut x: u32 = 0; let mut y: u32 = 0; let mut s: u32 = 1;
        while s < self.size as u32 {
            let rx = 1 & (idx / 2);
            let ry = 1 & (idx ^ rx);
            let (nx,ny) = Self::rot(s, x, y, rx, ry);
            x = nx + s*rx; y = ny + s*ry; idx /= 4; s *= 2;
        }
        (x as u16, y as u16)
    }
    fn rot(n: u32, x: u32, y: u32, rx: u32, ry: u32) -> (u32,u32) {
        if ry == 0 {
            if rx == 1 { return (n-1-x, n-1-y); }
            return (y,x);
        }
        (x,y)
    }
    pub fn position_to_xy(&self, pos: u64, genome_size: u64) -> (u16,u16) {
        let idx = ((pos as f64 / genome_size as f64) * self.total_positions as f64)
            .floor().min((self.total_positions-1) as f64) as u32;
        self.index_to_xy(idx)
    }
}
