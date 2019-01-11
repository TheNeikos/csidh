use crate::galois::{GaloisElement, LargeUint};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Curve {
    pub a: GaloisElement,
    b: GaloisElement,
}

impl Curve {
    pub fn new(a: LargeUint, b: LargeUint) -> Curve {
        Curve {
            a: GaloisElement::from_large_uint(a),
            b: GaloisElement::from_large_uint(b),
        }
    }

    pub fn right_side(&self, x: &GaloisElement) -> GaloisElement {
        let mut ret = *x;
        ret.square();
        let t = self.a * *x;
        ret.add_from(&t);
        ret.add_from(&crate::global::GAL_1);
        ret.mul_with(x);
        return ret;
    }

    pub fn contains(&self, p: &Point) -> bool {
        let left = (p.y * p.y) * self.b;
        let right = self.right_side(&p.x);

        left == right
    }

    fn recover(p: &Point, q: &ProjectivePoint, o: &ProjectivePoint) -> Point {
        let v1 = p.x * q.z;
        let v2 = q.x + v1;
        let v3 = q.x - v1;
        let v3 = v3 * v3;
        let v3 = v3 * o.x;
        let v1 = q.z * GaloisElement::from_u64(2) * p.curve.a;
        let v2 = v2 + v1;
        let v4 = p.x * q.x;
        let v4 = v4 + q.z;
        let v2 = v2 * v4;
        let v1 = v1 * q.z;
        let v2 = v2 - v1;
        let v2 = v2 * o.z;
        let y  = v2 - v3;
        let v1 = p.y * GaloisElement::from_u64(2) * p.curve.b;
        let v1 = v1 * q.z;
        let v1 = v1 * o.z;
        let x  = v1 * q.x;
        let z  = v1 * q.z;

        Point {
            curve: p.curve,
            x,
            y,
            z,
        }
    }

    pub fn isogeny(k: u64, ker: &ProjectivePoint, point: &ProjectivePoint, a: &ProjectivePoint)
        -> (ProjectivePoint, ProjectivePoint)
    {
        let mut t = [ker.z, ker.x, ker.x, ker.z];
        let mut tmp0;
        let mut tmp1;
        let mut q = ProjectivePoint::new({
            let x = ker.x * point.x;
            tmp0 = ker.z * point.z;
            x - tmp0
        }, {
            let z = ker.z * point.x;
            tmp0 = ker.x * point.z;
            z - tmp0
        });
        let mut m = [*ker; 3];
        m[1] = ker.double(&a.x);

        for i in 1..(k as usize / 2) {
            if i >= 2 {
                m[i % 3] = m[(i -1) % 3].add(ker, &m[(i - 2) % 3]);
            }

            tmp0 = m[i % 3].x * t[0];
            tmp1 = m[i % 3].z * t[1];
            t[0] = tmp0 + tmp1;

            t[1].mul_with(&m[i % 3].x);

            tmp0 = m[i % 3].z * t[2];
            tmp1 = m[i % 3].x * t[3];
            t[2] = tmp0 + tmp1;

            t[3].mul_with(&m[i % 3].z);


            tmp0 = point.x * m[i % 3].x;
            tmp1 = point.z * m[i % 3].z;
            tmp0.sub_from(&tmp1);
            q.x.mul_with(&tmp0);

            tmp0 = point.x * m[i % 3].z;
            tmp1 = point.z * m[i % 3].x;
            tmp0.sub_from(&tmp1);
            q.z.mul_with(&tmp0);
        }

        t[0].mul_with(&{t[1]});
        t[0].add_from(&{t[0]});

        t[1].square();

        t[2].mul_with(&{t[3]});
        t[2].add_from(&{t[2]});

        t[3].square();

        tmp0 = t[1] * t[2];
        tmp1 = t[0] * t[3];
        tmp0.sub_from(&tmp1);
        tmp0.mul_with(&a.z);
        tmp1 = tmp0 + tmp0;
        tmp0.add_from(&tmp1);

        tmp1 = t[1] * t[3];
        tmp1.mul_with(&a.x);

        let mut a_x = tmp1 - tmp0;

        t[3].square();
        let mut a_z = t[3];

        q.x.square();
        q.z.square();

        // We normalize here right away

        a_z.inverse();
        a_x.mul_with(&a_z);

        let new_a = ProjectivePoint::new(a_x, GaloisElement::from_u64(1));

        return (new_a, q);
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Point {
    curve: Curve,
    x: GaloisElement,
    y: GaloisElement,
    z: GaloisElement,
}

impl Point {
    pub fn new(curve: Curve, x: LargeUint, y: LargeUint) -> Point {
        Point {
            x: GaloisElement::from_large_uint(x),
            y: GaloisElement::from_large_uint(y),
            z: GaloisElement::from_u64(1),
            curve,
        }
    }

    pub fn multiply(&self, k: &LargeUint) -> Point {
        let (x0, x1) = self.projectivize().ladder(&self.curve.a, &k);
        let q = Curve::recover(self, &x0, &x1);
        return q.unproject();
    }

    fn projectivize(&self) -> ProjectivePoint {

        let zero = GaloisElement::from_u64(0);

        if self.x == zero || self.z == zero {
            ProjectivePoint {
                x: GaloisElement::from_u64(1),
                z: GaloisElement::from_u64(0),
            }
        } else {
            ProjectivePoint {
                x: self.x.clone(),
                z: GaloisElement::from_u64(1),
            }
        }
    }

    fn unproject(self) -> Point {
        assert!(self.z != GaloisElement::from_u64(0));

        let x = self.x / self.z;
        let y = self.y / self.z;
        let z = self.z / self.z;

        Point {
            x, y, z,
            curve: self.curve,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ProjectivePoint {
    pub x: GaloisElement,
    pub z: GaloisElement,
}

impl ProjectivePoint {
    pub fn new(x: GaloisElement, z: GaloisElement) -> ProjectivePoint {
        ProjectivePoint {
            x, z
        }
    }

    pub fn is_infinity(&self) -> bool {
        self.z == GaloisElement::from_u64(0)
    }

    pub fn ladder(&self, a: &GaloisElement, k: &LargeUint) -> (ProjectivePoint, ProjectivePoint) {
        let mut x0 = self.clone();
        let mut x1 = self.double(a);

        if k == &LargeUint::from_u64(1) {
            return (x0, x1);
        }

        let l = k.bits();

        // BUG! Rust can't do backwards ranges...
        for i in (0..=(l-2)).rev() {
            if !k.bit(i) {
                let temp = x0.add(&x1, self);
                x1 = temp;
                x0 = x0.double(a);
            } else {
                x0 = x0.add(&x1, self);
                x1 = x1.double(a);
            }
        }

        return (x0, x1);
    }


    pub fn add(&self, other: &ProjectivePoint, orig: &ProjectivePoint) -> ProjectivePoint {
        let v0 = self.x + self.z;
        let v1 = other.x - other.z;
        let v1 = v1 * v0;
        let v0 = self.x - self.z;
        let v2 = other.x + other.z;
        let v2 = v2 * v0;
        let v3 = v1 + v2;
        let v3 = v3 * v3;
        let v4 = v1 - v2;
        let v4 = v4 * v4;
        let x = orig.z * v3;
        let z = orig.x * v4;

        return ProjectivePoint {
            x,
            z,
        }
    }

    pub fn double(&self, a: &GaloisElement) -> ProjectivePoint {
        let v1 = self.x + self.z;
        let v1 = v1 * v1;
        let v2 = self.x - self.z;
        let v2 = v2 * v2;
        let x = v1 * v2;
        let v1 = v1 - v2;
        // Here was a bug! I forgot that division -> Modulo Ring
        let a_2 = *a + GaloisElement::from_u64(2);
        let v3 = v1 * (a_2 / GaloisElement::from_u64(4));
        let v3 = v3 + v2;
        let z = v1 * v3;

        return ProjectivePoint {
            x,
            z,
        }
    }
}


#[cfg(test)]
mod test {
    use super::*;
    use crate::global;

    #[test]
    fn check_mul() {
        let x = LargeUint::parse_bytes(b"2051044887188588280366899510711463515184102432059522841387541984999186019238289110841661333718393379209806643406155944602233875537370058705956384966209858");
        let y = LargeUint::parse_bytes(b"2999054700883294606115636709285947688603015463995111523694534197644452886751843273757676343103953201273958036952062931228773567734286840492294219977378136");

        let other_x = LargeUint::parse_bytes(b"1254817631949275079030490581963578364746575569014839158947538007979236709253796922466332191140273712204313677321924940880514829958528954596325165920058277");
        let other_y = LargeUint::parse_bytes(b"2381495309685763751265865484184529659090354786855457591442552214156841700513768692570497752099605704710183797526595611214891101033449784504091079214700929");

        let a = 0u32.into();
        let b = 1u32.into();
        let curve = Curve::new(a, b);
        let point = Point::new(curve, x, y);
        assert!(curve.contains(&point));
        let other_point = Point::new(curve, other_x, other_y);
        assert!(curve.contains(&other_point));

        let multiplied = point.multiply(&LargeUint::from_u64(9u64)).unproject();

        println!("X: {}\nY: {}\nZ: {}\n\nX: {}\nY: {}\nZ: {}",
                 multiplied.x.into_large_uint(), multiplied.y.into_large_uint(), multiplied.z.into_large_uint(),
                 other_point.x.into_large_uint(), other_point.y.into_large_uint(), other_point.z.into_large_uint());

        assert_eq!(other_point, multiplied);
    }
}
