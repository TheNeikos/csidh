use num_bigint::BigUint;
use num_traits::identities::{One, Zero};

use crate::ffield::{Field, FieldElement};

use crate::montgomery::Curve as MontgomeryCurve;

#[derive(Debug, Clone, PartialEq)]
pub struct EdwardsCurve {
    pub field: Field,
    pub a: FieldElement,
    pub d: FieldElement,
}

impl EdwardsCurve {
    pub fn to_montgomery(&self) -> MontgomeryCurve {
        let a = &self.field.get(2u32) * (&self.a + &self.d);
        let c = &self.a - &self.d;

        MontgomeryCurve::new(self.field.clone(), a, c);
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct EdwardsPoint {
    curve: Curve,
    x: FieldElement,
    y: FieldElement,
    z: FieldElement,
}
