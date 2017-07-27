//
// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.
//

/*
#[cfg(test)]
mod tests {

    use croaring::Bitmap;

    #[test]
    fn play() {

        let mut rb1 = Bitmap::create();
        rb1.add(1);
        rb1.add(2);
        rb1.add(3);
        rb1.add(4);
        rb1.add(5);
        rb1.add(100);
        rb1.add(1000);
        rb1.run_optimize();

        let mut rb2 = Bitmap::create();
        rb2.add(3);
        rb2.add(4);
        rb2.add(1000);
        rb2.run_optimize();

        let mut rb3 = Bitmap::create();

        assert_eq!(rb1.cardinality(), 7);
        assert!(rb1.contains(3));

        rb1.and_inplace(&rb2);
        rb3.add(5);
        rb3.or_inplace(&rb1);

        let mut rb4 = Bitmap::fast_or(&[&rb1, &rb2, &rb3]);

        rb1.and_inplace(&rb2);
        println!("{:?}", rb1);

        rb3.add(5);
        rb3.or_inplace(&rb1);

        println!("{:?}", rb1);

        rb3.add(5);
        rb3.or_inplace(&rb1);

        println!("{:?}", rb3.as_slice());
        println!("{:?}", rb3);
        println!("{:?}", rb4);

        rb4 = Bitmap::fast_or(&[&rb1, &rb2, &rb3]);

        println!("{:?}", rb4);
    }
}
*/
