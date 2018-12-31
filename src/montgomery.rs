use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::identities::{One, Zero};

use crate::global;

#[derive(Debug, Clone)]
struct Curve {
    a: BigInt,
}

impl Curve {
    pub fn new(a: BigInt) -> Curve {
        Curve {
            a
        }
    }

    pub fn recover(p: &Point, x0: &ProjectivePoint, x1: &ProjectivePoint) -> Point {
        Point {
            curve: p.curve.clone(),
            x: 0.into(),
            y: 0.into(),
        }
    }
}

#[derive(Debug, Clone)]
struct Point {
    curve: Curve,
    x: BigInt,
    y: BigInt,
}

impl Point {
    pub fn new(curve: Curve, x: BigInt, y: BigInt) -> Point {
        Point { curve, x, y, }
    }

    pub fn multiply(&self, k: BigInt) -> Point {
        let (x0, x1) = Point::ladder(k, self.projectivize());
        let q = Curve::recover(self, &x0, &x1);
        return q;
    }

    fn ladder(k: BigInt, p: ProjectivePoint) -> (ProjectivePoint, ProjectivePoint) {
        let mut x1 = p.double();
        let mut x0 = p.clone();
        let l = k.bits();

        for i in (l-2)..0 {
            let b = (&k >> i) & BigInt::one();
            if b == BigInt::zero() {
                x0 = x0.double();
                x1 = x1.add(&x0, &p);
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
            z: 0.into(),
        }
    }
}

#[derive(Debug, Clone)]
struct ProjectivePoint {
    curve: Curve,
    x: BigInt,
    z: BigInt,
}

impl ProjectivePoint {

    pub fn double(&self) -> ProjectivePoint {
        let v1 = &self.x + &self.z;
        let v1 = &v1 * &v1;
        let v2 = &self.x - &self.z;
        let v2 = &v2 * &v2;
        let x: BigInt = &v1 * &v2;

        let v1 = &v1 - &v2;
        let v3 = (&self.curve.a + 2)/4 * &v1;
        let v3 = &v3 + &v2;
        let z: BigInt = &v1 * &v3;

        return ProjectivePoint {
            curve: self.curve.clone(),
            x: x.mod_floor(&*global::p),
            z: z.mod_floor(&*global::p),
        }
    }

    pub fn add(&self, other: &ProjectivePoint, orig: &ProjectivePoint) -> ProjectivePoint {
        let v0 = &self.x + &self.z;
        let v1 = &other.x - &other.z;
        let v1 = &v0 * &v1;
        let v0 = &self.x - &self.z;
        let v2 = &other.x + &other.z;
        let v2 = &v2 * &v0;
        let v3 = &v1 + &v2;
        let v3 = &v3 * &v3;
        let v4 = &v1 - &v2;
        let v4 = &v4 * &v4;
        let x = &orig.z * &v3;
        let z = &orig.x * &v4;
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
        let p = Point::new(Curve::new(1.into()), 2.into());
    }
}
