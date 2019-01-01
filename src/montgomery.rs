use num_bigint::BigUint;
use num_integer::Integer;
use num_traits::identities::{One, Zero};
use num_traits::pow::Pow;

use crate::global;
use crate::ffield::inverse;

#[derive(Debug, Clone, PartialEq)]
struct Curve {
    a: BigUint,
    b: BigUint,
}

impl Curve {
    pub fn new(a: BigUint, b: BigUint) -> Curve {
        Curve {
            a, b,
        }
    }

    pub fn contains(&self, p: &Point) -> bool {
        let left = &self.b * (p.y.pow(2u32));
        let right = p.x.pow(3u32) + &self.a * p.x.pow(2u32) + &p.x;

        left.mod_floor(&*global::p) == right.mod_floor(&*global::p)
    }

    pub fn recover(p: &Point, q: &ProjectivePoint, o: &ProjectivePoint) -> Point {
        let v1 = &p.x * &q.z;
        let v2 = &q.x + &v1;
        let v3 = &q.x + inverse(&v1);
        let v3 = &v3 * &v3;
        let v3 = &v3 * &o.x;
        let v1 = BigUint::from(2u32) * &p.curve.a * &q.z;
        let v2 = &v2 + &v1;
        let v4 = &p.x * &q.x;
        let v4 = &v4 + &q.z;
        let v2 = &v2 * &v4;
        let v1 = &v1 * &q.z;
        let v2 = &v2 + inverse(&v1);
        let v2 = &v2 * &o.z;
        let y  = &v2 + inverse(&v3);
        let v1 = BigUint::from(2u32) * &p.curve.b * &p.y;
        let v1 = &v1 * &q.z;
        let v1 = &v1 * &o.z;
        let x  = &v1 * &q.x;
        let z  = &v1 * &q.z;

        Point {
            curve: p.curve.clone(),
            x: x.mod_floor(&*global::p),
            y: y.mod_floor(&*global::p),
            z: z.mod_floor(&*global::p)
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
struct Point {
    curve: Curve,
    x: BigUint,
    y: BigUint,
    z: BigUint,
}

impl Point {
    pub fn new(curve: Curve, x: BigUint, y: BigUint) -> Point {
        Point { curve, x, y, z: 1u32.into() }
    }

    pub fn multiply(&self, k: &BigUint) -> Point {
        let (x0, x1) = Point::ladder(&k, self.projectivize());
        let q = Curve::recover(self, &x0, &x1);
        return q;
    }

    fn ladder(k: &BigUint, p: ProjectivePoint) -> (ProjectivePoint, ProjectivePoint) {
        let mut x0 = p.clone();
        let mut x1 = p.double();
        let l = k.bits();
        assert!(k >> (l-1) == BigUint::one());

        // BUG! Rust can't do backwards ranges...
        for i in (0..(l-1)).rev() {
            let b = (k >> i) & BigUint::one();
            if b == BigUint::zero() {
                let temp = x0.add(&x1, &p);
                x1 = temp;
                x0 = x0.double();
            } else {
                x0 = x0.add(&x1, &p);
                x1 = x1.double();
            }
        }

        return (x0, x1);
    }

    pub fn projectivize(&self) -> ProjectivePoint {
        ProjectivePoint {
            curve: self.curve.clone(),
            x: self.x.clone(),
            z: 1u32.into(),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
struct ProjectivePoint {
    curve: Curve,
    x: BigUint,
    z: BigUint,
}

impl ProjectivePoint {

    pub fn double(&self) -> ProjectivePoint {
        let v1 = &self.x + &self.z;
        let v1 = &v1 * &v1;
        let v2 = &self.x + inverse(&self.z);
        let v2 = &v2 * &v2;
        let x: BigUint = &v1 * &v2;

        let v1 = &v1 + inverse(&v2);
        // Here was a bug! I forgot that division -> Modulo Ring
        let v3 = (&self.curve.a + 2u32) * &*global::inv_4 * &v1;
        let v3 = &v3 + &v2;
        let z: BigUint = &v1 * &v3;

        assert!(x != BigUint::zero());
        assert!(z != BigUint::zero());
        return ProjectivePoint {
            curve: self.curve.clone(),
            x: x.mod_floor(&*global::p),
            z: z.mod_floor(&*global::p),
        }
    }

    pub fn add(&self, other: &ProjectivePoint, orig: &ProjectivePoint) -> ProjectivePoint {
        let v0 = &self.x + &self.z;
        let v1 = &other.x + inverse(&other.z);
        let v1 = &v1 * &v0;
        let v0 = &self.x + inverse(&self.z);
        let v2 = &other.x + &other.z;
        let v2 = &v2 * &v0;
        let v3 = &v1 + &v2;
        let v3 = &v3 * &v3;
        let v4 = &v1 + inverse(&v2);
        let v4 = &v4 * &v4;
        let x = &orig.z * &v3;
        let z = &orig.x * &v4;

        assert!(x != BigUint::zero());
        assert!(z != BigUint::zero());
        return ProjectivePoint {
            curve: self.curve.clone(),
            x: x.mod_floor(&*global::p),
            z: z.mod_floor(&*global::p)
        }
    }
}


#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn check_mul() {
        let x = BigUint::parse_bytes(b"1610136167417732114551750628834547689599377638968928848446790657003564444920832601989878342871943872732914644562225752648144282986871867284681507748894542", 10).unwrap();
        let y = BigUint::parse_bytes(b"4451648534257700053946276422085028026878964190648779111476216123161445920291343445824029807401583577610995159754431055928148978190346828059269568767034510", 10).unwrap();

        let otherX = BigUint::parse_bytes(b"5078812358276015976382593126191062178939149646412005198913702881619839138992439581471890018490803712280132371742602680177163514165386685881446765077133001", 10).unwrap();
        let otherY = BigUint::parse_bytes(b"2224187386808147719950316248531353129580936339598092469828138282231361456048530546430333214876724108025609984131177898306863977372254974363123638850862735", 10).unwrap();

        let curve = Curve::new(0u32.into(), 1u32.into());
        let point = Point::new(curve.clone(), x, y);
        assert!(curve.contains(&point));
        let other_point = Point::new(curve.clone(), otherX, otherY);
        assert!(curve.contains(&other_point));

        let multiplied = point.multiply(&BigUint::from(15u32));

        println!("X: {}\nY: {}\n\nX: {}\nY: {}", multiplied.x, multiplied.y, other_point.x, other_point.y);

        assert_eq!(other_point, multiplied);
    }
}
