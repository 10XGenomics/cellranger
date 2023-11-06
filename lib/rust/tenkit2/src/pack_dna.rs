// Turning this into a grab bad of DNA utility functions, should rename.

// Pack DNA bases.  This in effect turns Ns into As.  Only implemented for certain
// kmers sizes.

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// K = 10
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn unpack_bases_10(y: &[u8; 3], x: &mut [u8; 10]) {
    for i in 0..3 {
        let mut yy = *y;
        for j in 0..4 {
            if 4 * i + 4 - j > 10 {
                continue;
            }
            let m = yy[i] % 4;
            match m {
                0 => {
                    x[4 * i + 4 - j - 1] = b'A';
                }
                1 => {
                    x[4 * i + 4 - j - 1] = b'C';
                }
                2 => {
                    x[4 * i + 4 - j - 1] = b'G';
                }
                _ => {
                    x[4 * i + 4 - j - 1] = b'T';
                }
            }
            yy[i] >>= 2;
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// K = 16
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn pack_bases_16x(x: &[u8], y: &mut [u8; 4]) {
    for i in 0..4 {
        y[i] = 0;
        for j in 0..4 {
            y[i] <<= 2;
            match x[4 * i + j] {
                b'C' => {
                    y[i] += 1;
                }
                b'G' => {
                    y[i] += 2;
                }
                b'T' => {
                    y[i] += 3;
                }
                _ => {}
            }
        }
    }
}

pub fn unpack_bases_16(y: &[u8; 4], x: &mut [u8; 16]) {
    for i in 0..4 {
        let mut yy = *y;
        for j in 0..4 {
            let m = yy[i] % 4;
            match m {
                0 => {
                    x[4 * i + 4 - j - 1] = b'A';
                }
                1 => {
                    x[4 * i + 4 - j - 1] = b'C';
                }
                2 => {
                    x[4 * i + 4 - j - 1] = b'G';
                }
                _ => {
                    x[4 * i + 4 - j - 1] = b'T';
                }
            }
            yy[i] >>= 2;
        }
    }
}

// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓
// K = 90
// ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓

pub fn pack_bases_80(x: &[u8; 80], y: &mut [u8; 20]) {
    for i in 0..20 {
        y[i] = 0;
        for j in 0..4 {
            y[i] <<= 2;
            match x[4 * i + j] {
                b'C' => {
                    y[i] += 1;
                }
                b'G' => {
                    y[i] += 2;
                }
                b'T' => {
                    y[i] += 3;
                }
                _ => {}
            }
        }
    }
}

pub fn unpack_bases_80(y: &[u8; 20], x: &mut [u8; 80]) {
    for i in 0..20 {
        let mut yy = *y;
        for j in 0..4 {
            let m = yy[i] % 4;
            match m {
                0 => {
                    x[4 * i + 4 - j - 1] = b'A';
                }
                1 => {
                    x[4 * i + 4 - j - 1] = b'C';
                }
                2 => {
                    x[4 * i + 4 - j - 1] = b'G';
                }
                _ => {
                    x[4 * i + 4 - j - 1] = b'T';
                }
            }
            yy[i] >>= 2;
        }
    }
}
pub fn reverse_complement(x: &mut [u8]) {
    x.reverse();
    for entry in x {
        *entry = match entry {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => *entry,
        }
    }
}
