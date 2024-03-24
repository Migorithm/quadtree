#[derive(Debug, Clone, Copy)]
struct Coordinate {
    x: f64,
    y: f64,
}

impl Coordinate {
    fn new(x: f64, y: f64) -> Self {
        Coordinate { x, y }
    }
}

#[derive(Debug)]
struct Boundary {
    top_left_coor: Coordinate,
    bottom_right_coor: Coordinate,
}

impl Boundary {
    fn new(top_left_coor: Coordinate, bottom_right_coor: Coordinate) -> Boundary {
        Boundary {
            top_left_coor,
            bottom_right_coor,
        }
    }

    fn contains(&self, point: &Coordinate) -> bool {
        point.x >= self.top_left_coor.x
            && point.x < self.bottom_right_coor.x
            && point.y >= self.top_left_coor.y
            && point.y < self.bottom_right_coor.y
    }
    fn intersects(&self, boundary: &Boundary) -> bool {
        self.top_left_coor.x < boundary.bottom_right_coor.x
            && self.bottom_right_coor.x >= boundary.top_left_coor.x
            && self.top_left_coor.y < boundary.bottom_right_coor.y
            && self.bottom_right_coor.y >= boundary.top_left_coor.y
    }

    fn disect(&self) -> (Self, Self, Self, Self) {
        let mid_coor = Coordinate::new(
            (self.top_left_coor.x + self.bottom_right_coor.x) / 2.0,
            (self.top_left_coor.y + self.bottom_right_coor.y) / 2.0,
        );

        let nw = Boundary::new(self.top_left_coor, mid_coor);
        let ne = Boundary::new(
            Coordinate::new(mid_coor.x, self.top_left_coor.y),
            Coordinate::new(self.bottom_right_coor.x, mid_coor.y),
        );
        let sw = Boundary::new(
            Coordinate::new(self.top_left_coor.x, mid_coor.y),
            Coordinate::new(mid_coor.x, self.bottom_right_coor.y),
        );

        let se = Boundary::new(mid_coor, self.bottom_right_coor.clone());

        (nw, ne, sw, se)
    }
}

#[derive(Debug)]
struct Quadtree<T> {
    boundary: Boundary,
    capacity: usize,
    coordinates: Vec<Coordinate>,
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
            coordinates: Vec::with_capacity(capacity),
            interests: Vec::with_capacity(capacity),
            divided: false,
            nw: None,
            ne: None,
            sw: None,
            se: None,
        }
    }
    fn find_subquadtree_for_newly_added_coordinate(
        &mut self,
        point: &Coordinate,
    ) -> &mut Option<Box<Quadtree<T>>> {
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

    pub fn insert(&mut self, coor: Coordinate, interest: T) {
        if !self.boundary.contains(&coor) {
            return;
        }

        if self.coordinates.len() >= self.capacity {
            self.subdivide();
        }

        if self.divided {
            let sub_quadtree = self.find_subquadtree_for_newly_added_coordinate(&coor);
            sub_quadtree.as_mut().unwrap().insert(coor, interest);
        } else {
            self.coordinates.push(coor);
            self.interests.push(interest)
        }
    }

    fn subdivide(&mut self) {
        let (nw_boundary, ne_boundary, sw_boundary, se_boundary) = self.boundary.disect();

        self.nw = Some(Box::new(Quadtree::new(nw_boundary, self.capacity)));
        self.ne = Some(Box::new(Quadtree::new(ne_boundary, self.capacity)));
        self.sw = Some(Box::new(Quadtree::new(sw_boundary, self.capacity)));
        self.se = Some(Box::new(Quadtree::new(se_boundary, self.capacity)));

        self.divided = true;

        // At this point, `self(quadtree)` is not leaf node anymore. Reinsert its interests and coordinates into subquadtree
        std::mem::take(&mut self.coordinates)
            .into_iter()
            .zip(std::mem::take(&mut self.interests))
            .for_each(|(coor, interest)| self.insert(coor, interest))
    }

    fn query<'a: 'b, 'b>(&'a self, range: &Boundary, found: &mut Vec<(Coordinate, &'b T)>) {
        if !self.boundary.intersects(range) {
            return;
        }

        for (point, object) in self.coordinates.iter().zip(self.interests.iter()) {
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
    use crate::{Boundary, Coordinate, Quadtree};

    #[test]
    fn test_contains() {
        //GIVEN
        let boundary = Boundary::new(Coordinate::new(0.0, 0.0), Coordinate::new(100.0, 100.0));
        //WHEN
        let res = boundary.contains(&Coordinate { x: 0.5, y: 2.0 });
        //THEN
        assert!(res);
    }

    #[test]
    fn test_insert() {
        //GIVEN
        let mut qt = Quadtree::<&str>::new(
            Boundary::new(Coordinate::new(0.0, 0.0), Coordinate::new(100.0, 100.0)),
            4,
        );

        //WHEN
        qt.insert(Coordinate::new(30.0, 20.0), "1");
        qt.insert(Coordinate::new(10.0, 50.0), "2");
        qt.insert(Coordinate::new(60.0, 60.0), "3");
        qt.insert(Coordinate::new(70.0, 90.0), "4");
        qt.insert(Coordinate::new(100.0, 100.0), "5");

        //THEN
        assert!(qt.divided);
        assert!(qt.interests.is_empty());
        assert!(qt.coordinates.is_empty());
        assert_eq!(qt.se.unwrap().interests.len(), 3);
    }

    #[test]
    fn test_query() {
        let mut quadtree: Quadtree<&str> = Quadtree::new(
            Boundary::new(Coordinate::new(0.0, 0.0), Coordinate::new(100.0, 100.0)),
            4,
        );

        // Inserting some points
        quadtree.insert(Coordinate::new(30.0, 30.0), "A");
        quadtree.insert(Coordinate::new(10.0, 50.0), "B");
        quadtree.insert(Coordinate::new(70.0, 20.0), "C");
        quadtree.insert(Coordinate::new(80.0, 80.0), "D");
        quadtree.insert(Coordinate::new(80.0, 90.0), "E");

        // Querying points in a range
        let range = Boundary::new(Coordinate::new(50.0, 50.0), Coordinate::new(100.0, 100.0));
        let mut found_points = Vec::new();
        quadtree.query(&range, &mut found_points);

        assert!(quadtree.divided);
        assert_eq!(found_points.len(), 2);
    }
}
