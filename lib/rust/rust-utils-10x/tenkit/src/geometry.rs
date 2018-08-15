

use num::Num;

pub struct Point2D<T: Num + Copy + Clone> {
    x: T,
    y: T,
}

impl<T: Num + Copy + Clone> Point2D<T> {
    pub fn new(x: T, y: T) -> Self {
        Point2D {x, y}
    }

    pub fn origin() -> Self {
        Point2D {
            x: T::zero(),
            y: T::zero()
        }
    }

    #[inline(always)]
    pub fn squared_distance(&self, other: &Self) -> T {
        let dx = self.x - other.x;
        let dy = self.y - other.y;
        dx * dx + dy * dy
    }

    #[inline(always)]
    pub fn quadrance(&self, other: &Self) -> T {
        self.squared_distance(other)
    }
}