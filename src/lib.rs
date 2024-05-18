use std::{
    collections::BinaryHeap,
    fmt::Debug,
    iter::zip,
    ops::{Add, Deref, Div, Mul, Sub},
};

use serde::{Deserialize, Serialize};

#[derive(Debug, Default, Clone, Copy, Serialize, Deserialize)]
#[repr(transparent)]
pub struct NonNanFloat(f64);

impl PartialEq for NonNanFloat {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        if self.0.is_nan() {
            other.0.is_nan()
        } else {
            self.0 == other.0
        }
    }
}

impl PartialEq<f64> for NonNanFloat {
    #[inline]
    fn eq(&self, other: &f64) -> bool {
        if self.0.is_nan() {
            other.is_nan()
        } else {
            self.0 == *other
        }
    }
}

impl PartialOrd for NonNanFloat {
    #[inline]
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.0.partial_cmp(&other.0)
    }
}

impl Eq for NonNanFloat {}

impl Ord for NonNanFloat {
    #[inline]
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.partial_cmp(other).unwrap()
    }
}

impl Sub for NonNanFloat {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self(self.0 - other.0)
    }
}

impl Sub<f64> for NonNanFloat {
    type Output = Self;

    fn sub(self, other: f64) -> Self {
        Self(self.0 - other)
    }
}

impl Add for NonNanFloat {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        NonNanFloat(self.0 + other.0)
    }
}

impl Add<f64> for NonNanFloat {
    type Output = Self;

    fn add(self, other: f64) -> Self {
        NonNanFloat(self.0 + other)
    }
}

impl Mul for NonNanFloat {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        NonNanFloat(self.0 * other.0)
    }
}

impl Mul<f64> for NonNanFloat {
    type Output = Self;

    fn mul(self, other: f64) -> Self {
        NonNanFloat(self.0 * other)
    }
}

impl Div for NonNanFloat {
    type Output = Self;

    fn div(self, other: Self) -> Self {
        NonNanFloat(self.0 / other.0)
    }
}

impl Div<f64> for NonNanFloat {
    type Output = Self;

    fn div(self, other: f64) -> Self {
        NonNanFloat(self.0 / other)
    }
}

impl Deref for NonNanFloat {
    type Target = f64;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl NonNanFloat {
    fn new(value: f64) -> Self {
        assert!(!value.is_nan());
        NonNanFloat(value)
    }
}

impl From<f64> for NonNanFloat {
    fn from(value: f64) -> Self {
        NonNanFloat::new(value)
    }
}

#[derive(Debug, Default, Clone, Copy, Serialize, Deserialize)]
pub struct Coordinate {
    pub longitude: NonNanFloat,
    pub latitude: NonNanFloat,
}

pub trait TCoordianteFloat {
    fn into(self) -> NonNanFloat;
}
impl TCoordianteFloat for f64 {
    fn into(self) -> NonNanFloat {
        NonNanFloat::new(self)
    }
}
impl TCoordianteFloat for NonNanFloat {
    fn into(self) -> NonNanFloat {
        self
    }
}

impl Coordinate {
    pub fn new<F: TCoordianteFloat>(longitude: F, latitude: F) -> Self {
        Coordinate {
            longitude: longitude.into(),
            latitude: latitude.into(),
        }
    }

    fn distance(&self, other: &Coordinate) -> f64 {
        ((self.longitude - other.longitude).powi(2) + (self.latitude - other.latitude).powi(2))
            .sqrt()
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Boundary {
    pub top_left_coor: Coordinate,
    pub bottom_right_coor: Coordinate,
}

impl Boundary {
    pub fn new(top_left_coor: Coordinate, bottom_right_coor: Coordinate) -> Boundary {
        Boundary {
            top_left_coor,
            bottom_right_coor,
        }
    }

    fn contains(&self, point: &Coordinate) -> bool {
        point.longitude >= self.top_left_coor.longitude
            && point.longitude < self.bottom_right_coor.longitude
            && point.latitude >= self.top_left_coor.latitude
            && point.latitude < self.bottom_right_coor.latitude
    }
    fn intersects(&self, boundary: &Boundary) -> bool {
        self.top_left_coor.longitude < boundary.bottom_right_coor.longitude
            && self.bottom_right_coor.longitude >= boundary.top_left_coor.longitude
            && self.top_left_coor.latitude < boundary.bottom_right_coor.latitude
            && self.bottom_right_coor.latitude >= boundary.top_left_coor.latitude
    }

    // used from subdivide function!
    fn disect(&self) -> (Self, Self, Self, Self) {
        let mid_coor = Coordinate::new(
            (self.top_left_coor.longitude + self.bottom_right_coor.longitude) / 2.0,
            (self.top_left_coor.latitude + self.bottom_right_coor.latitude) / 2.0,
        );

        let nw = Boundary::new(self.top_left_coor, mid_coor);
        let ne = Boundary::new(
            Coordinate::new(mid_coor.longitude, self.top_left_coor.latitude),
            Coordinate::new(self.bottom_right_coor.longitude, mid_coor.latitude),
        );
        let sw = Boundary::new(
            Coordinate::new(self.top_left_coor.longitude, mid_coor.latitude),
            Coordinate::new(mid_coor.longitude, self.bottom_right_coor.latitude),
        );

        let se = Boundary::new(mid_coor, self.bottom_right_coor);

        (nw, ne, sw, se)
    }

    fn distance(&self, point: &Coordinate) -> NonNanFloat {
        if self.contains(point) {
            return NonNanFloat(0.0);
        }

        let dx = if point.longitude < self.top_left_coor.longitude {
            self.top_left_coor.longitude - point.longitude
        } else if point.longitude >= self.bottom_right_coor.longitude {
            point.longitude - self.bottom_right_coor.longitude
        } else {
            NonNanFloat(0.0)
        };

        let dy = if point.latitude < self.top_left_coor.latitude {
            self.top_left_coor.latitude - point.latitude
        } else if point.latitude >= self.bottom_right_coor.latitude {
            point.latitude - self.bottom_right_coor.latitude
        } else {
            NonNanFloat(0.0)
        };

        NonNanFloat((dx.powi(2) + dy.powi(2)).sqrt())
    }
}

#[derive(Debug, Serialize, Deserialize)]
pub struct Quadtree<T: Ord> {
    pub boundary: Boundary,
    pub capacity: usize,
    pub coordinates: Vec<Coordinate>,
    pub interests: Vec<T>,

    children: Vec<Quadtree<T>>,
}

impl<T: Ord + Debug> Quadtree<T> {
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
    fn search<'a>(
        &'a self,
        nearest_neighbers: &mut BinaryHeap<(NonNanFloat, &'a T)>,
        query_point: &Coordinate,
        k: usize,
    ) {
        if nearest_neighbers.len() == k {
            if let Some((max_distance, _)) = nearest_neighbers.peek() {
                if self.boundary.distance(query_point) >= *max_distance {
                    return;
                }
            }
        }

        for (coord, v) in zip(&self.coordinates, &self.interests) {
            let distance = coord.distance(query_point);
            if nearest_neighbers.len() < k {
                nearest_neighbers.push((NonNanFloat::new(distance), v));
            } else if let Some((max_distance, _)) = nearest_neighbers.peek() {
                if distance >= **max_distance {
                    continue;
                }
                if nearest_neighbers.len() == k {
                    nearest_neighbers.pop();
                }
                nearest_neighbers.push((NonNanFloat::new(distance), v));
            }
        }

        for child in self.children.iter() {
            child.search(nearest_neighbers, query_point, k);
        }
    }

    pub fn find_nearest_neighbors(&self, query_point: &Coordinate, k: usize) -> Vec<&T> {
        let mut nearest_neighbors = BinaryHeap::new();
        self.search(&mut nearest_neighbors, query_point, k);

        nearest_neighbors
            .into_sorted_vec()
            .iter()
            .map(|(_, v)| *v)
            .collect()
    }
}

#[cfg(test)]
mod test {
    use crate::{Boundary, Coordinate, NonNanFloat, Quadtree};

    #[test]
    fn non_nan_float() {
        //GIVEN
        let a = NonNanFloat(10.0);
        let b = NonNanFloat(20.0);
        let c = NonNanFloat(30.0);
        let d = NonNanFloat(40.0);
        let nan = NonNanFloat(f64::NAN);
        let another_nan = NonNanFloat(f64::NAN);
        let plain_float = 10.0;

        //WHEN
        let res = a + b;
        let res2 = c - d;
        let res3 = a * b;
        let res4 = c / d;
        let res5 = a > b;
        let res6 = a < b;
        let res7 = a == b;
        let res8 = nan == another_nan;
        let res9 = a == plain_float;
        let res10 = a + plain_float;
        let res11 = a - plain_float;
        let res12 = a * plain_float;
        let res13 = a / plain_float;

        //THEN
        assert_eq!(res, NonNanFloat(30.0));
        assert_eq!(res2, NonNanFloat(-10.0));
        assert_eq!(res3, NonNanFloat(200.0));
        assert_eq!(res4, NonNanFloat(0.75));
        assert!(!res5);
        assert!(res6);
        assert!(!res7);
        assert!(res8);
        assert!(res9);
        assert_eq!(res10, NonNanFloat(20.0));
        assert_eq!(res11, NonNanFloat(0.0));
        assert_eq!(res12, NonNanFloat(100.0));
        assert_eq!(res13, NonNanFloat(1.0));
    }

    #[test]
    fn test_boundary_contains() {
        //GIVEN
        let boundary = Boundary::new(Coordinate::new(0.0, 0.0), Coordinate::new(100.0, 100.0));
        //WHEN
        let res = boundary.contains(&Coordinate::new(0.5, 2.0));
        //THEN
        assert!(res);
    }

    #[test]
    fn test_boundary_distance() {
        //GIVEN
        let boundary = Boundary::new(Coordinate::new(0.0, 0.0), Coordinate::new(100.0, 100.0));

        // WHEN
        // The pair means (coordinate, expected distance)
        let cases = vec![
            // point is inside the boundary
            (Coordinate::new(0.5, 2.0), 0.0),
            (Coordinate::new(0.0, 0.0), 0.0),
            (Coordinate::new(100.0, 100.0), 0.0),
            (Coordinate::new(0.0, 100.0), 0.0),
            (Coordinate::new(100.0, 0.0), 0.0),
            (Coordinate::new(50.0, 50.0), 0.0),
            (Coordinate::new(50.0, 0.0), 0.0),
            (Coordinate::new(0.0, 50.0), 0.0),
            (Coordinate::new(100.0, 50.0), 0.0),
            (Coordinate::new(50.0, 100.0), 0.0),
            // point is outside the boundary
            (Coordinate::new(0.0, 200.0), 100.0),
            (Coordinate::new(200.0, 0.0), 100.0),
            (Coordinate::new(200.0, 200.0), f64::sqrt(20000.0)),
            (Coordinate::new(50.0, 200.0), 100.0),
            (Coordinate::new(200.0, 50.0), 100.0),
            (Coordinate::new(200.0, 100.0), 100.0),
            (Coordinate::new(100.0, 200.0), 100.0),
            (Coordinate::new(150.0, 150.0), f64::sqrt(5000.0)),
            (Coordinate::new(101.0, 103.0), f64::sqrt(10.0)),
            // real data
            (Coordinate::new(37.512428, 127.054513), 27.054513),
        ];

        //THEN
        for (coor, expected) in cases {
            let res = boundary.distance(&coor);
            assert_eq!(*res, expected);
        }
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
        struct TestCase<'a> {
            points: Vec<(f64, f64, &'a str)>,
            query_point: (f64, f64),
            k: usize,
            expected: Vec<&'a str>,
        }

        //GIVEN
        let mut cases = vec![
            TestCase {
                points: vec![
                    (30.0, 30.0, "A"),
                    (10.0, 50.0, "B"),
                    (70.0, 20.0, "C"),
                    (80.0, 80.0, "D"),
                    (80.0, 90.0, "E"),
                    (60.0, 90.0, "F"),
                    (2500.0, 2700.0, "G"),
                    (1700.0, 1500.0, "H"),
                    (4993.0, 4999.0, "I"),
                    (4993.0, 4330.0, "J"),
                    (4993.0, 4500.0, "K"),
                ],
                query_point: (4999.0, 4950.0),
                k: 8,
                expected: vec!["I", "J", "K", "G", "H", "E", "D", "F"],
            },
            TestCase {
                points: vec![
                    (30.0, 30.0, "A"),
                    (10.0, 50.0, "B"),
                    (70.0, 20.0, "C"),
                    (80.0, 80.0, "D"),
                    (80.0, 90.0, "E"),
                    (60.0, 90.0, "F"),
                    (2500.0, 2700.0, "G"),
                    (1700.0, 1500.0, "H"),
                    (4993.0, 4999.0, "I"),
                    (4993.0, 4330.0, "J"),
                    (4993.0, 4500.0, "K"),
                ],
                query_point: (4999.0, 4950.0),
                k: 5,
                expected: vec!["I", "J", "K", "G", "H"],
            },
            TestCase {
                points: vec![
                    (30.0, 30.0, "A"),
                    (10.0, 50.0, "B"),
                    (70.0, 20.0, "C"),
                    (80.0, 80.0, "D"),
                    (80.0, 90.0, "E"),
                    (60.0, 90.0, "F"),
                    (2500.0, 2700.0, "G"),
                    (1700.0, 1500.0, "H"),
                    (4993.0, 4999.0, "I"),
                    (4993.0, 4330.0, "J"),
                    (4993.0, 4500.0, "K"),
                ],
                query_point: (4999.0, 4950.0),
                k: 3,
                expected: vec!["I", "J", "K"],
            },
            TestCase {
                points: vec![(30.0, 30.0, "A"), (10.0, 50.0, "B"), (70.0, 20.0, "C")],
                query_point: (4999.0, 4950.0),
                k: 5,
                expected: vec!["A", "B", "C"],
            },
            TestCase {
                points: Vec::new(),
                query_point: (4999.0, 4950.0),
                k: 5,
                expected: Vec::new(),
            },
        ];

        //THEN
        for case in cases.iter_mut() {
            let mut quadtree: Quadtree<&str> = Quadtree::new(
                Boundary::new(Coordinate::new(0.0, 0.0), Coordinate::new(5000.0, 5000.0)),
                4,
            );

            for (longitude, latitude, interest) in &case.points {
                quadtree.insert(Coordinate::new(*longitude, *latitude), interest);
            }

            let mut interests = quadtree.find_nearest_neighbors(
                &Coordinate::new(case.query_point.0, case.query_point.1),
                case.k,
            );

            interests.sort();
            case.expected.sort();

            assert_eq!(interests.len(), case.expected.len());
            for (interest, expected) in interests.iter().zip(case.expected.iter()) {
                assert_eq!(*interest, expected);
            }
        }
    }

    #[test]
    fn test_insert_million() {
        use rand::Rng;

        //GIVEN
        // million record
        let mut rng = rand::thread_rng();
        let million_record = (0..1000000)
            .map(|_| (rng.gen_range(0.0..4999.9), rng.gen_range(0.0..4999.9)))
            .collect::<Vec<_>>();

        let mut quadtree: Quadtree<&str> = Quadtree::new(
            Boundary::new(Coordinate::new(0.0, 0.0), Coordinate::new(5000.0, 5000.0)),
            4,
        );
        let instance = std::time::Instant::now();

        // WHEN
        for (longitude, latitude) in million_record {
            quadtree.insert(Coordinate::new(longitude, latitude), "A");
        }
        let second_instant = std::time::Instant::now();

        // THEN
        println!(
            "Time passed : {}",
            second_instant.duration_since(instance).as_millis()
        );
    }

    #[test]
    fn test_near_5_from_million_records() {
        use rand::Rng;

        //GIVEN
        // million record
        let mut rng = rand::thread_rng();
        let million_record = (0..1000000)
            .map(|_| (rng.gen_range(-180.0..180.0), rng.gen_range(-90.0..90.0)))
            .collect::<Vec<_>>();

        let mut quadtree: Quadtree<String> = Quadtree::new(
            Boundary::new(Coordinate::new(-180.0, -90.0), Coordinate::new(180.0, 90.0)),
            20,
        );

        for (longitude, latitude) in million_record {
            quadtree.insert(
                Coordinate::new(longitude, latitude),
                format!("long : {longitude}, lat: {latitude}"),
            );
        }

        // WHEN
        let instance = std::time::Instant::now();

        let k_business = quadtree.find_nearest_neighbors(&Coordinate::new(127.0, 38.1), 5);

        // THEN
        let second_instant = std::time::Instant::now();
        println!(
            "Time passed : {}",
            second_instant.duration_since(instance).as_millis()
        );
        assert_eq!(k_business.len(), 5);
    }
}
