//
// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
//

/// Direction of motion in a DeBruijn graph
#[derive(Copy, Clone, Debug)]
pub enum Dir {
    Left,
    Right,
}

impl Dir {
    /// Return a fresh Dir with the opposite direction
    pub fn flip(&self) -> Dir {
        match *self {
            Dir::Left => Dir::Right,
            Dir::Right => Dir::Left,
        }
    }

    /// Return a fresh Dir opposite direction if do_flip == True
    pub fn cond_flip(&self, do_flip: bool) -> Dir {
        if do_flip {
            self.flip()
        } else {
            *self
        }
    }

    /// Pick between two alternatives, depending on the direction
    pub fn pick<T>(&self, if_left: T, if_right: T) -> T {
        match self {
            &Dir::Left => if_left,
            &Dir::Right => if_right,
        }
    }
}
