#[derive(Debug, Clone, Copy)]
struct Point {
    x: f64,
    y: f64,
}

impl Point {
    fn new(x: f64, y: f64) -> Self {
        Point { x, y }
    }
}

#[derive(Debug)]
struct Boundary {
    min_point: Point,
    max_point: Point,
}

impl Boundary {
    fn new(min_point: Point, max_point: Point) -> Boundary {
        Boundary {
            min_point,
            max_point,
        }
    }

    fn contains(&self, point: &Point) -> bool {
        point.x >= self.min_point.x
            && point.x < self.max_point.x
            && point.y >= self.min_point.y
            && point.y < self.max_point.y
    }
    fn intersects(&self, boundary: &Boundary) -> bool {
        self.min_point.x < boundary.max_point.x
            && self.max_point.x >= boundary.min_point.x
            && self.min_point.y < boundary.max_point.y
            && self.max_point.y >= boundary.min_point.y
    }

    fn mid_point(&self) -> Point {
        Point::new(
            (self.min_point.x + self.max_point.x) / 2.0,
            (self.min_point.y + self.max_point.y) / 2.0,
        )
    }

    fn nw(&self) -> Self {
        Self {
            min_point: Point::new(self.min_point.x, self.min_point.y),
            max_point: self.mid_point(),
        }
    }
    fn ne(&self) -> Self {
        let mid_point = self.mid_point();

        Self {
            min_point: Point::new(mid_point.x, self.min_point.y),
            max_point: Point::new(self.max_point.x, mid_point.y),
        }
    }

    fn sw(&self) -> Self {
        let mid_point = self.mid_point();

        Self {
            min_point: Point::new(self.min_point.x, mid_point.y),
            max_point: Point::new(mid_point.x, self.max_point.y),
        }
    }

    fn se(&self) -> Self {
        Self {
            min_point: self.mid_point(),
            max_point: Point::new(self.max_point.x, self.max_point.y),
        }
    }

    fn divide(&self) -> (Self, Self, Self, Self) {
        (self.nw(), self.ne(), self.sw(), self.se())
    }
}

#[derive(Debug)]
struct Quadtree<T> {
    boundary: Boundary,
    capacity: usize,
    points: Vec<Point>,
    interests: Vec<T>,
    divided: bool,
    nw: Option<Box<Quadtree<T>>>,
    ne: Option<Box<Quadtree<T>>>,
    sw: Option<Box<Quadtree<T>>>,
    se: Option<Box<Quadtree<T>>>,
}

impl<T> Quadtree<T> {
    pub fn new(boundary: Boundary, capacity: usize) -> Quadtree<T> {
        Quadtree {
            boundary,
            capacity,
            points: Vec::with_capacity(capacity),
            interests: Vec::with_capacity(capacity),
            divided: false,
            nw: None,
            ne: None,
            sw: None,
            se: None,
        }
    }
    fn locate_subboundary_to_put(&mut self, point: &Point) -> &mut Option<Box<Quadtree<T>>> {
        if self.nw.as_ref().unwrap().boundary.contains(point) {
            &mut self.nw
        } else if self.ne.as_ref().unwrap().boundary.contains(point) {
            &mut self.ne
        } else if self.sw.as_ref().unwrap().boundary.contains(point) {
            &mut self.sw
        } else if self.se.as_ref().unwrap().boundary.contains(point) {
            &mut self.se
        } else {
            unreachable!()
        }
    }

    pub fn insert(&mut self, point: Point, object: T) {
        // If the given point is outside the boundary, return.

        if !self.boundary.contains(&point) {
            return;
        }

        if !self.divided && (self.points.len() < self.capacity) {
            self.points.push(point);
            self.interests.push(object);
        } else {
            if !self.divided {
                self.subdivide();
            }

            // if divided, it has sub-quadtree and you can try insert harmlessly!
            let sub_boundary = self.locate_subboundary_to_put(&point);
            sub_boundary.as_mut().unwrap().insert(point, object);
        }
    }

    fn subdivide(&mut self) {
        let (nw_boundary, ne_boundary, sw_boundary, se_boundary) = self.boundary.divide();

        self.nw = Some(Box::new(Quadtree::new(nw_boundary, self.capacity)));
        self.ne = Some(Box::new(Quadtree::new(ne_boundary, self.capacity)));
        self.sw = Some(Box::new(Quadtree::new(sw_boundary, self.capacity)));
        self.se = Some(Box::new(Quadtree::new(se_boundary, self.capacity)));

        self.divided = true;

        // Reinsert points
        for (point, object) in std::mem::take(&mut self.points)
            .into_iter()
            .zip(std::mem::take(&mut self.interests).into_iter())
        {
            self.insert(point, object);
        }
    }

    fn query<'a: 'b, 'b>(&'a self, range: &Boundary, found: &mut Vec<(Point, &'b T)>) {
        if !self.boundary.intersects(range) {
            return;
        }

        for (point, object) in self.points.iter().zip(self.interests.iter()) {
            if range.contains(point) {
                found.push((*point, object));
            }
        }

        if self.divided {
            if let Some(ref nw) = self.nw {
                nw.query(range, found);
            }
            if let Some(ref ne) = self.ne {
                ne.query(range, found);
            }
            if let Some(ref sw) = self.sw {
                sw.query(range, found);
            }
            if let Some(ref se) = self.se {
                se.query(range, found);
            }
        }
    }
}

#[cfg(test)]
mod test {
    use crate::{Boundary, Point, Quadtree};

    #[test]
    fn test_contains() {
        //GIVEN
        let boundary = Boundary::new(Point::new(0.0, 0.0), Point::new(100.0, 100.0));
        //WHEN
        let res = boundary.contains(&Point { x: 0.5, y: 2.0 });
        //THEN
        assert!(res);
    }
    #[test]
    fn test() {
        let mut quadtree: Quadtree<&str> = Quadtree::new(
            Boundary::new(Point::new(0.0, 0.0), Point::new(100.0, 100.0)),
            4,
        );

        // Inserting some points
        quadtree.insert(Point::new(30.0, 30.0), "A");
        quadtree.insert(Point::new(10.0, 50.0), "B");
        quadtree.insert(Point::new(70.0, 20.0), "C");
        quadtree.insert(Point::new(80.0, 80.0), "D");

        quadtree.insert(Point::new(80.0, 90.0), "E");

        // Querying points in a range
        let range = Boundary::new(Point::new(50.0, 50.0), Point::new(100.0, 100.0));
        let mut found_points = Vec::new();
        quadtree.query(&range, &mut found_points);

        assert!(quadtree.divided);
        assert_eq!(found_points.len(), 2);
    }
}
