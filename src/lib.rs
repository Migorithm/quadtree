#[derive(Debug, Clone, Copy)]
pub struct Coordinate {
    x: f64,
    y: f64,
}

impl Coordinate {
    pub fn new(x: f64, y: f64) -> Self {
        Coordinate { x, y }
    }

    fn distance(&self, other: &Coordinate) -> f64 {
        ((self.x - other.x).powi(2) + (self.y - other.y).powi(2)).sqrt()
    }
}

#[derive(Debug)]
pub struct Boundary {
    top_left_coor: Coordinate,
    bottom_right_coor: Coordinate,
}

impl Boundary {
    pub fn new(top_left_coor: Coordinate, bottom_right_coor: Coordinate) -> Boundary {
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

    // used from subdivide function!
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

        let se = Boundary::new(mid_coor, self.bottom_right_coor);

        (nw, ne, sw, se)
    }
}

#[derive(Debug)]
pub struct Quadtree<T> {
    boundary: Boundary,
    capacity: usize,
    coordinates: Vec<Coordinate>,
    interests: Vec<T>,

    children: Vec<Quadtree<T>>,
}

impl<T> Quadtree<T> {
    pub fn new(boundary: Boundary, capacity: usize) -> Quadtree<T> {
        Quadtree {
            boundary,
            capacity,
            coordinates: Vec::with_capacity(capacity),
            interests: Vec::with_capacity(capacity),
            children: Vec::with_capacity(4),
        }
    }
    fn find_subquadtree_for_newly_added_coordinate(
        &mut self,
        point: &Coordinate,
    ) -> &mut Quadtree<T> {
        for child in self.children.iter_mut() {
            if child.boundary.contains(point) {
                return child;
            }
        }
        unreachable!()
    }

    pub fn insert(&mut self, coor: Coordinate, interest: T) {
        if !self.boundary.contains(&coor) {
            return;
        }

        if self.coordinates.len() >= self.capacity {
            self.subdivide();
        }

        if !self.children.is_empty() {
            let sub_quadtree = self.find_subquadtree_for_newly_added_coordinate(&coor);
            sub_quadtree.insert(coor, interest);
        } else {
            self.coordinates.push(coor);
            self.interests.push(interest)
        }
    }

    // the following devides self into four sub boundaries
    fn subdivide(&mut self) {
        let (nw_boundary, ne_boundary, sw_boundary, se_boundary) = self.boundary.disect();

        self.children.extend(vec![
            Quadtree::new(nw_boundary, self.capacity),
            Quadtree::new(ne_boundary, self.capacity),
            Quadtree::new(sw_boundary, self.capacity),
            Quadtree::new(se_boundary, self.capacity),
        ]);

        // At this point, `self(quadtree)` is not leaf node anymore. Reinsert its interests and coordinates into subquadtree
        std::mem::take(&mut self.coordinates)
            .into_iter()
            .zip(std::mem::take(&mut self.interests))
            .for_each(|(coor, interest)| self.insert(coor, interest))
    }

    // To find interests in specific boundary use the following
    pub fn query<'a: 'b, 'b>(&'a self, range: &Boundary, found: &mut Vec<(Coordinate, &'b T)>) {
        if !self.boundary.intersects(range) {
            return;
        }

        for (point, object) in self.coordinates.iter().zip(self.interests.iter()) {
            if range.contains(point) {
                found.push((*point, object));
            }
        }

        if !self.children.is_empty() {
            for child in self.children.iter() {
                child.query(range, found);
            }
        }
    }

    // recursive function
    fn search<'a: 'b, 'b>(&'a self, distances: &mut Vec<(&'b T, f64)>, query_point: &Coordinate) {
        if self.boundary.contains(query_point) && self.children.is_empty() {
            for (coordinate, interest) in self.coordinates.iter().zip(self.interests.iter()) {
                let distance = coordinate.distance(query_point);
                distances.push((interest, distance));
            }
        } else {
            for child in self.children.iter() {
                child.search(distances, query_point);
            }
        }
    }

    pub fn find_nearest_neighbors(&self, query_point: &Coordinate, k: usize) -> Vec<&T> {
        let mut distances = vec![];

        self.search(&mut distances, query_point);
        distances.sort_by(|(_, dis1), (_, dis2)| dis1.partial_cmp(dis2).unwrap());
        distances
            .into_iter()
            .take(k)
            .map(|(c, _)| c)
            .collect::<Vec<_>>()
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
        qt.insert(Coordinate::new(90.0, 90.0), "5");

        //THEN

        assert!(!qt.children.is_empty());
        assert!(qt.interests.is_empty());
        assert!(qt.coordinates.is_empty());
        assert_eq!(qt.children.get(3).unwrap().interests.len(), 3);
    }

    #[test]
    fn test_query() {
        //GIVEN
        let mut quadtree: Quadtree<&str> = Quadtree::new(
            Boundary::new(Coordinate::new(0.0, 0.0), Coordinate::new(100.0, 100.0)),
            4,
        );

        //WHEN
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

        //THEN
        assert!(!quadtree.children.is_empty());
        assert_eq!(found_points.len(), 2);
    }

    #[test]
    fn test_k_nearest() {
        //GIVEN
        let mut quadtree: Quadtree<&str> = Quadtree::new(
            Boundary::new(Coordinate::new(0.0, 0.0), Coordinate::new(5000.0, 5000.0)),
            4,
        );

        // Inserting some points
        quadtree.insert(Coordinate::new(30.0, 30.0), "A");
        quadtree.insert(Coordinate::new(10.0, 50.0), "B");
        quadtree.insert(Coordinate::new(70.0, 20.0), "C");
        quadtree.insert(Coordinate::new(80.0, 80.0), "D");
        quadtree.insert(Coordinate::new(80.0, 90.0), "E");
        quadtree.insert(Coordinate::new(60.0, 90.0), "F");
        quadtree.insert(Coordinate::new(2500.0, 2700.0), "G");
        quadtree.insert(Coordinate::new(1700.0, 1500.0), "H");
        quadtree.insert(Coordinate::new(4993.0, 4999.0), "I");
        quadtree.insert(Coordinate::new(4993.0, 4330.0), "J");
        quadtree.insert(Coordinate::new(4993.0, 4500.0), "K");

        //WHEN
        let query_point = Coordinate::new(4999.0, 4950.0);
        let interests = quadtree.find_nearest_neighbors(&query_point, 3);

        //THEN
        assert_eq!(interests.len(), 3);
        let mut expected = vec!["I", "J", "K"];
        for interest in interests {
            assert!(expected.contains(interest));
            let pos = expected.iter().position(|x| x == interest);
            expected.remove(pos.unwrap());
        }
        assert!(expected.is_empty());
    }

    #[test]
    fn test_insert_million() {
        use rand::Rng;

        //GIVEN
        // million record
        let mut rng = rand::thread_rng();
        let million_record = (0..100000000)
            .map(|_| (rng.gen_range(0.0..4999.9), rng.gen_range(0.0..4999.9)))
            .collect::<Vec<_>>();

        let mut quadtree: Quadtree<&str> = Quadtree::new(
            Boundary::new(Coordinate::new(0.0, 0.0), Coordinate::new(5000.0, 5000.0)),
            4,
        );
        let instance = std::time::Instant::now();

        // WHEN
        for (x, y) in million_record {
            quadtree.insert(Coordinate { x, y }, "A");
        }
        let second_instant = std::time::Instant::now();

        // THEN
        println!(
            "Time passed : {}",
            second_instant.duration_since(instance).as_millis()
        );
    }
}
