// Copyright 2014-2015 Johannes Köster, Peer Aramillo Irizar.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::borrow::Borrow;

use bit_set::BitSet;

#[derive(Debug, PartialEq, Eq)]
pub struct Alphabet {
    pub symbols: BitSet,
}

impl Alphabet {
    pub fn new<C, T>(symbols: T) -> Self
    where
        C: Borrow<u8>,
        T: IntoIterator<Item = C>,
    {
        let mut s = BitSet::new();
        s.extend(symbols.into_iter().map(|c| *c.borrow() as usize));

        Alphabet { symbols: s }
    }

    pub fn insert(&mut self, a: u8) {
        self.symbols.insert(a as usize);
    }

    pub fn is_word<C, T>(&self, text: T) -> bool
    where
        C: Borrow<u8>,
        T: IntoIterator<Item = C>,
    {
        text.into_iter()
            .all(|c| self.symbols.contains(*c.borrow() as usize))
    }

    pub fn max_symbol(&self) -> Option<u8> {
        self.symbols.iter().max().map(|a| a as u8)
    }

    pub fn len(&self) -> usize {
        self.symbols.len()
    }

    pub fn is_empty(&self) -> bool {
        self.symbols.is_empty()
    }

    pub fn intersection(&self, other: &Alphabet) -> Self {
        return Alphabet {
            symbols: self.symbols.intersection(&other.symbols).collect(),
        };
    }

    pub fn difference(&self, other: &Alphabet) -> Self {
        return Alphabet {
            symbols: self.symbols.difference(&other.symbols).collect(),
        };
    }

    pub fn union(&self, other: &Alphabet) -> Self {
        return Alphabet {
            symbols: self.symbols.union(&other.symbols).collect(),
        };
    }
}

pub fn english_ascii_lower_alphabet() -> Alphabet {
    Alphabet::new(&b"abcdefghijklmnopqrstuvwxyz"[..])
}

pub fn english_ascii_upper_alphabet() -> Alphabet {
    Alphabet::new(&b"ABCDEFGHIJKLMNOPQRSTUVWXYZ"[..])
}
