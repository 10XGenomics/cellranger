extern crate bincode;
extern crate byteorder;
extern crate itertools;
#[macro_use]
extern crate log;
extern crate regex;
extern crate rust_htslib;
#[macro_use]
extern crate serde_derive;
extern crate fxhash;
extern crate num;

pub mod lane;
pub mod chunk_bam;
pub mod collections;
pub mod geometry;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
