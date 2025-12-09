#![deny(missing_docs)]
pub(super) fn pack_bases_80(x: &[u8; 80], y: &mut [u8; 20]) {
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

pub(super) fn unpack_bases_80(y: &[u8; 20], x: &mut [u8; 80]) {
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
