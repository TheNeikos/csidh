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
        let inv_b = GaloisElement::from_u64(1) / self.b;
        let two_x = p.x * p.x;
        let x_c = p.x * inv_b;

        let left = (p.y * p.y) * self.b * inv_b;
        let right = two_x * x_c + two_x * inv_b * self.a + x_c;

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

    pub fn isogeny(&self, l: u64, ker: &ProjectivePoint, point: &ProjectivePoint) -> (GaloisElement, ProjectivePoint) {
        let d = (l - 1) / 2;
        let one = GaloisElement::from_u64(1);

        let mut rho = GaloisElement::from_u64(0);
        let mut rho_tilde = GaloisElement::from_u64(0);
        let mut pi = GaloisElement::from_u64(1);

        for i in 1..d {
            let x_k = ker.ladder(&LargeUint::from_u64(i)).0.x;
            rho = rho + x_k;
            rho_tilde = rho_tilde + one / x_k;
            pi = pi * x_k;
        }

        let six = GaloisElement::from_u64(6);
        let new_a = (six * rho - six * rho_tilde + self.a) * pi * pi;

        let mut x = GaloisElement::from_u64(1);
        let mut z = GaloisElement::from_u64(1);

        let x_min_z = point.x - point.z;
        let x_plu_z = point.x + point.z;

        for i in 1..d {
            let p_i = point.ladder(&LargeUint::from_u64(i)).0;

            let xi_min_zi = p_i.x - p_i.z;
            let xi_plu_zi = p_i.x + p_i.z;

            let pair_1 = x_min_z * xi_plu_zi;
            let pair_2 = x_plu_z * xi_min_zi;

            x = x * (pair_1 + pair_2);
            z = z * (pair_1 - pair_2);

        }

        x = x * x * point.x;
        z = z * z * point.z;

        let new_p = ProjectivePoint::new(*self, x, z);

        return (new_a, new_p);
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
        let (x0, x1) = self.projectivize().ladder(&k);
        let q = Curve::recover(self, &x0, &x1);
        return q.unproject();
    }

    fn projectivize(&self) -> ProjectivePoint {

        let zero = GaloisElement::from_u64(0);

        if self.x == zero || self.z == zero {
            ProjectivePoint {
                curve: self.curve,
                x: GaloisElement::from_u64(1),
                z: GaloisElement::from_u64(0),
            }
        } else {
            ProjectivePoint {
                curve: self.curve,
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

#[derive(Debug, Clone, PartialEq)]
pub struct ProjectivePoint {
    curve: Curve,
    x: GaloisElement,
    z: GaloisElement,
}

impl ProjectivePoint {
    pub fn new(curve: Curve, x: GaloisElement, z: GaloisElement) -> ProjectivePoint {
        ProjectivePoint {
            x, z, curve,
        }
    }

    pub fn set_curve(&mut self, curve: Curve) {
        self.curve = curve;
    }

    pub fn is_infinity(&self) -> bool {
        self.z == GaloisElement::from_u64(0)
    }

    pub fn ladder(&self, k: &LargeUint) -> (ProjectivePoint, ProjectivePoint) {
        let mut x0 = self.clone();
        let mut x1 = self.double();

        if k == &LargeUint::from_u64(1) {
            return (x0, x1);
        }

        let l = k.bits();

        // BUG! Rust can't do backwards ranges...
        for i in (0..=(l-2)).rev() {
            if k.bit(i) {
                let temp = x0.add(&x1, self);
                x1 = temp;
                x0 = x0.double();
            } else {
                x0 = x0.add(&x1, self);
                x1 = x1.double();
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
            curve: self.curve,
            x,
            z,
        }
    }

    pub fn double(&self) -> ProjectivePoint {
        let v1 = self.x + self.z;
        let v1 = v1 * v1;
        let v2 = self.x - self.z;
        let v2 = v2 * v2;
        let x = v1 * v2;
        let v1 = v1 - v2;
        // Here was a bug! I forgot that division -> Modulo Ring
        let a_2 = self.curve.a + GaloisElement::from_u64(2);
        let v3 = v1 * (a_2 / GaloisElement::from_u64(4));
        let v3 = v3 + v2;
        let z = v1 * v3;

        return ProjectivePoint {
            curve: self.curve,
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
                 multiplied.x, multiplied.y, multiplied.z,
                 other_point.x, other_point.y, other_point.z);

        assert_eq!(other_point, multiplied);
    }
}
