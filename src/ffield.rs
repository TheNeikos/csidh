use crate::global;
use num_bigint::BigUint;
use num_integer::Integer;

use std::fmt;
use std::ops::{Add, Sub, Mul, Div, BitXor};
use std::rc::Rc;

#[derive(Debug, PartialEq, Clone)]
pub struct FieldElement {
    order: Rc<BigUint>,
    val: BigUint,
}

impl FieldElement {
    pub fn is_square(&self) -> bool {
        let p_min_1_half = ((&*self.order) - BigUint::from(1u32)) / BigUint::from(2u32);
        let pow = self.val.modpow(&p_min_1_half, &self.order);
        return pow == 1u32.into();
    }

    fn mul_inverse(&self) -> FieldElement {
        let p_min_2 = &*self.order - BigUint::from(2u32);
        FieldElement {
            order: self.order.clone(),
            val: self.val.modpow(&p_min_2, &self.order),
        }
    }

    fn add_inverse(&self) -> FieldElement {
        FieldElement {
            order: self.order.clone(),
            val: &*self.order - &self.val,
        }
    }
}

impl fmt::Display for FieldElement {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.val)
    }
}

impl Add for &FieldElement {
    type Output = FieldElement;

    fn add(self, other: &FieldElement) -> FieldElement {
        let val = &self.val + &other.val;
        FieldElement {
            val: val.mod_floor(&self.order),
            order: self.order.clone(),
        }
    }
}

impl Add<FieldElement> for FieldElement {
    type Output = FieldElement;

    fn add(self, other: FieldElement) -> FieldElement {
        let val = &self.val + &other.val;
        FieldElement {
            val: val.mod_floor(&self.order),
            order: self.order.clone(),
        }
    }
}

impl Add<&FieldElement> for FieldElement {
    type Output = FieldElement;

    fn add(self, other: &FieldElement) -> FieldElement {
        let val = &self.val + &other.val;
        FieldElement {
            val: val.mod_floor(&self.order),
            order: self.order.clone(),
        }
    }
}

impl Mul for &FieldElement {
    type Output = FieldElement;

    fn mul(self, other: &FieldElement) -> FieldElement {
        let val = &self.val * &other.val;
        FieldElement {
            val: val.mod_floor(&self.order),
            order: self.order.clone(),
        }
    }
}

impl Mul for FieldElement {
    type Output = FieldElement;

    fn mul(self, other: FieldElement) -> FieldElement {
        let val = &self.val * &other.val;
        FieldElement {
            val: val.mod_floor(&self.order),
            order: self.order.clone(),
        }
    }
}

impl Mul<&FieldElement> for FieldElement {
    type Output = FieldElement;

    fn mul(self, other: &FieldElement) -> FieldElement {
        let val = &self.val * &other.val;
        FieldElement {
            val: val.mod_floor(&self.order),
            order: self.order.clone(),
        }
    }
}

impl Sub for &FieldElement {
    type Output = FieldElement;

    fn sub(self, other: &FieldElement) -> FieldElement {
        let val = self + &other.add_inverse();
        FieldElement {
            val: val.val.mod_floor(&self.order),
            order: self.order.clone(),
        }
    }
}

impl Div for &FieldElement {
    type Output = FieldElement;

    fn div(self, other: &FieldElement) -> FieldElement {
        let val = self * &other.mul_inverse();
        FieldElement {
            val: val.val.mod_floor(&self.order),
            order: self.order.clone(),
        }
    }
}

impl Mul<&BigUint> for &FieldElement {
     type Output = FieldElement;

    fn mul(self, other: &BigUint) -> FieldElement {
        let val = &self.val * other;
        FieldElement {
            val: val.mod_floor(&self.order),
            order: self.order.clone(),
        }
    }
}

impl Mul<&BigUint> for FieldElement {
     type Output = FieldElement;

    fn mul(self, other: &BigUint) -> FieldElement {
        let val = &self.val * other;
        FieldElement {
            val: val.mod_floor(&self.order),
            order: self.order.clone(),
        }
    }
}

impl Mul<BigUint> for FieldElement {
     type Output = FieldElement;

    fn mul(self, other: BigUint) -> FieldElement {
        let val = &self.val * other;
        FieldElement {
            val: val.mod_floor(&self.order),
            order: self.order.clone(),
        }
    }
}

impl Mul<BigUint> for &FieldElement {
     type Output = FieldElement;

    fn mul(self, other: BigUint) -> FieldElement {
        let val = &self.val * other;
        FieldElement {
            val: val.mod_floor(&self.order),
            order: self.order.clone(),
        }
    }
}

impl BitXor for &FieldElement {
     type Output = FieldElement;

    fn bitxor(self, other: &FieldElement) -> FieldElement {
        let val = &self.val ^ &other.val;
        FieldElement {
            val: val,
            order: self.order.clone(),
        }
    }
}

#[derive(Debug, PartialEq, Clone)]
pub struct Field {
    order: Rc<BigUint>,
}

impl Field {
    pub fn new(order: BigUint) -> Field {
        Field {
            order: Rc::new(order)
        }
    }

    pub fn get<S: Into<BigUint>>(&self, v: S) -> FieldElement {
        FieldElement {
            order: self.order.clone(),
            val: v.into().mod_floor(&self.order),
        }
    }

    pub fn order(&self) -> &BigUint {
        &self.order
    }
}

#[cfg(test)]
mod test {
    use num_bigint::BigUint;
    use num_integer::Integer;
    use crate::global;
    use super::*;

    #[test]
    fn check_fields() {
        let field = Field::new(global::p.clone());
        let five = field.get(5u32);
        let two = field.get(2u32);

        let seven = &five + &two;
        assert_eq!(seven, field.get(7u32));

        let three = &five - &two;
        assert_eq!(three, field.get(3u32));

        let ten = &five * &two;
        assert_eq!(ten, field.get(10u32));

        let new_five = &ten / &two;
        assert_eq!(five, new_five);
    }

    #[test]
    fn sanity_checks() {
        let field = Field::new(BigUint::from(43u32));
        let one = field.get(1u32);
        let five = field.get(36u32);
        assert_eq!(&five / &five, one);
    }
}

