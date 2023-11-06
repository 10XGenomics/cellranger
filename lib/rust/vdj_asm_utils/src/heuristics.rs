// Control over the algorithm.

pub struct Heuristics {
    pub free: bool,      // do denovo assembly?
    pub prim_trim: bool, // trim after primers?
}

impl Heuristics {
    pub fn new() -> Heuristics {
        Heuristics {
            free: false,
            prim_trim: true,
        }
    }
}

impl Default for Heuristics {
    fn default() -> Self {
        Self::new()
    }
}
