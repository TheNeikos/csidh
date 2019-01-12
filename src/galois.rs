use std::ops::{Add, Sub, Mul, Div};
use rand::{CryptoRng, Rng};

use crate::global::*;

pub const LIMBS: usize = 8;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LargeUint {
    pub elements: [u64; LIMBS],
}

impl LargeUint {
    pub fn new() -> LargeUint {
        LargeUint {
            elements: [0; LIMBS],
        }
    }

    pub fn from_u64(u: u64) -> LargeUint {
        LargeUint {
            elements: [u, 0, 0, 0, 0, 0, 0, 0],
        }
    }

    pub fn as_bytes(&self) -> Vec<u8> {
        use byteorder::{ByteOrder, LittleEndian};
        let mut bytes = vec![0; 8 * LIMBS];

        LittleEndian::write_u64_into(&self.elements, &mut bytes[..]);
        return bytes;
    }

    pub fn parse_bytes(s: &[u8]) -> LargeUint {
        s.iter().fold(LargeUint::new(), |mut acc, x| {
            acc.mul_with_u64(10);
            acc.add_from(&LargeUint::from_u64((x - b'0') as u64));
            acc
        })
    }

    pub fn add_from(&mut self, other: &LargeUint) -> bool {
        let mut carry: bool = false;
        for i in 0..LIMBS {
            let (temp, c) = self.elements[i].overflowing_add(carry as u64);
            carry = c;
            let (res, c) = temp.overflowing_add(other.elements[i]);
            self.elements[i] = res;
            carry |= c;
        }

        return carry;
    }

    fn sub_from(&mut self, other: &LargeUint) -> bool {
        let mut carry: bool = false;
        for i in 0..LIMBS {
            let (temp, c) = self.elements[i].overflowing_sub(carry as u64);
            carry = c;
            let (res, c) = temp.overflowing_sub(other.elements[i]);
            self.elements[i] = res;
            carry |= c;
        }

        return carry;
    }

    pub fn mul_with_u64(&mut self, other: u64) {
        let mut c = 0u64;
        for i in 0..LIMBS {
            let t = self.elements[i] as u128 * other as u128 + c as u128;
            c = (t >> 64) as u64;
            self.elements[i] = t as u64;
        }
    }

    pub fn bits(&self) -> u64 {
        for i in (0..LIMBS).rev() {
            if self.elements[i] == 0 {
                continue;
            }
            let zeros = self.elements[i].leading_zeros();
            return (64 - zeros as u64) + (i * 64) as u64;
        }
        return 0;
    }

    pub fn bit(&self, i: u64) -> bool {
        (self.elements[i as usize / 64] >> i % 64) & 1 == 1
    }
}

impl std::convert::From<u32> for LargeUint {
    fn from(u: u32) -> LargeUint {
        LargeUint::from_u64(u as u64)
    }
}

impl std::fmt::Display for LargeUint {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "[")?;
        for i in 0..LIMBS {
            write!(f, "0x{:016x}", self.elements[i])?;
            if i != LIMBS - 1 {
                write!(f, " ")?;
            }
        }
        write!(f, "]")
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GaloisElement {
    pub elements: [u64; LIMBS],
}

impl GaloisElement {
    pub fn from_u64(u: u64) -> GaloisElement {
        let lu = LargeUint {
            elements: [u, 0, 0, 0, 0, 0, 0, 0]
        };
        GaloisElement::from_large_uint(lu)
    }

    pub fn from_large_uint(lu: LargeUint) -> GaloisElement {
        let mut t = GaloisElement {
            elements: lu.elements,
        };

        t.mul_with(&R_SQUARED_MOD_P);
        return t;
    }

    pub fn into_large_uint(&self) -> LargeUint {
        let mut s = *self;
        s.mul_with(&LUINT_1);
        LargeUint {
            elements: s.elements,
        }
    }

    fn into_large_uint_priv(&self) -> LargeUint {
        LargeUint {
            elements: self.elements,
        }
    }

    pub fn random_element<R: Rng + CryptoRng>(rng: &mut R) -> GaloisElement {
        loop {
            let mut elems = [0u64; LIMBS];
            rng.fill(&mut elems);
            let m = (1u64 << PBITS % 64) - 1;
            elems[LIMBS - 1] &= m;

            for i in (0..LIMBS).rev() {
                if elems[i] < P.elements[i] {
                    return GaloisElement { elements: elems };
                } else if elems[i] > P.elements[i] {
                    break;
                }
            }
        }
    }

    pub fn sub_from(&mut self, other: &GaloisElement) -> bool {
        let mut s = self.into_large_uint_priv();
        let o = other.into_large_uint_priv();
        let r = s.sub_from(&o);
        if r {
            s.add_from(&P_INT);
        }
        self.elements = s.elements;
        return r;
    }

    pub fn add_from(&mut self, other: &GaloisElement) -> bool {
        let mut s = self.into_large_uint_priv();
        let o = other.into_large_uint_priv();
        let r = s.add_from(&o);
        self.elements = s.elements;
        self.reduce_once();
        return r;
    }

    pub fn mul_with(&mut self, other: &GaloisElement) {
        let mut temp = [0u64; LIMBS + 1];
        for k in 0..LIMBS {
            let r = |i| -> usize { (k + i) % (LIMBS + 1) };

            let m: u64 = INV_MIN_P_MOD_R.wrapping_mul(self.elements[k].wrapping_mul(other.elements[0])
                                                 .wrapping_add(temp[r(0)]));
            let mut carry = false;
            let mut other_carry = false;
            for i in 0..LIMBS {
                let u: u128 = m as u128 * P.elements[i] as u128;

                let (res, c) = temp[r(i)].overflowing_add(other_carry as u64);
                other_carry = c;
                temp[r(i)] = res;

                let (res, c) = temp[r(i)].overflowing_add(u as u64);
                other_carry |= c;
                temp[r(i)] = res;

                let (res, c) = temp[r(i+1)].overflowing_add(carry as u64);
                carry = c;
                temp[r(i+1)] = res;

                let (res, c) = temp[r(i+1)].overflowing_add((u >> 64) as u64);
                carry |= c;
                temp[r(i+1)] = res;
            }
            temp[r(LIMBS)] += other_carry as u64;

            let mut carry = false;
            let mut other_carry = false;
            for i in 0..LIMBS {
                let u: u128 = self.elements[k] as u128 * other.elements[i] as u128;

                let (res, c) = temp[r(i)].overflowing_add(other_carry as u64);
                other_carry = c;
                temp[r(i)] = res;

                let (res, c) = temp[r(i)].overflowing_add(u as u64);
                other_carry |= c;
                temp[r(i)] = res;

                let (res, c) = temp[r(i+1)].overflowing_add(carry as u64);
                carry = c;
                temp[r(i+1)] = res;

                let (res, c) = temp[r(i+1)].overflowing_add((u >> 64) as u64);
                carry |= c;
                temp[r(i+1)] = res;
            }
            temp[r(LIMBS)] += other_carry as u64;
        }

        for i in 0..LIMBS {
            self.elements[i] = temp[(LIMBS + i) % (LIMBS + 1)];
        }

        self.reduce_once();
    }

    pub fn square(&mut self) -> GaloisElement {
        self.mul_with(&{*self});
        *self
    }

    fn pow(&mut self, exp: &LargeUint) {
        let mut prev: GaloisElement = *self;
        *self = GAL_1;
        for k in 0..LIMBS {
            let mut t = exp.elements[k];
            for _ in 0..64 {
                if (t & 1) != 0 {
                    self.mul_with(&prev);
                }
                prev.square();
                t >>= 1;
            }
        }
    }

    pub fn inverse(&mut self) {
        self.pow(&P_MINUS_2);
    }

    pub fn is_square(&self) -> bool {
        let mut t = *self;
        t.pow(&P_MINUS_1_HALVES);
        t == GAL_1
    }

    fn reduce_once(&mut self) {
        let mut temp = self.clone();
        if !temp.sub_from(&P) {
            *self = temp;
        }
    }
}

impl Add for GaloisElement {
    type Output = GaloisElement;

    fn add(mut self, other: GaloisElement) -> GaloisElement {
        self.add_from(&other);
        self
    }
}

impl std::fmt::Display for GaloisElement {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "[")?;
        for i in 0..LIMBS {
            write!(f, "0x{:016x}", self.elements[i])?;
            if i != LIMBS - 1 {
                write!(f, " ")?;
            }
        }
        write!(f, "] mod p")
    }
}

impl Sub for GaloisElement {
    type Output = GaloisElement;

    fn sub(mut self, other: GaloisElement) -> GaloisElement {
        self.sub_from(&other);
        self
    }
}

impl Mul for GaloisElement {
    type Output = GaloisElement;

    fn mul(mut self, other: GaloisElement) -> GaloisElement {
        self.mul_with(&other);
        self
    }
}

impl Div for GaloisElement {
    type Output = GaloisElement;

    fn div(mut self, mut other: GaloisElement) -> GaloisElement {
        &other.inverse();
        self.mul_with(&other);
        self
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn check_add() {
        let mut one = LargeUint {
            elements: [1, 0, 0, 0, 0, 0, 0, 0]
        };

        let two = LargeUint {
            elements: [2, 0, 0, 0, 0, 0, 0, 0]
        };

        let three = LargeUint {
            elements: [3, 0, 0, 0, 0, 0, 0, 0]
        };

        one.add_from(&two);
        assert_eq!(one, three);

        let one = LargeUint {
            elements: [1, 0, 0, 0, 0, 0, 0, 0]
        };

        let mut max_one = LargeUint {
            elements: [u64::max_value(), 0, 0, 0, 0, 0, 0, 0]
        };

        let max_two = LargeUint {
            elements: [0, 1, 0, 0, 0, 0, 0, 0]
        };

        max_one.add_from(&one);
        assert_eq!(max_one, max_two);
    }

    #[test]
    fn check_sub() {
        let one = LargeUint {
            elements: [1, 0, 0, 0, 0, 0, 0, 0]
        };

        let two = LargeUint {
            elements: [2, 0, 0, 0, 0, 0, 0, 0]
        };

        let mut three = LargeUint {
            elements: [3, 0, 0, 0, 0, 0, 0, 0]
        };

        three.sub_from(&two);
        assert_eq!(one, three);
    }

    #[test]
    fn check_mul() {
        let mut two = GaloisElement::from_u64(2);
        let four = GaloisElement::from_u64(4);

        two.mul_with(&two.clone());
        assert_eq!(two, four);
    }

    #[test]
    fn check_enc_dec() {
        let one = GaloisElement::from_u64(1);
        assert_ne!(one.elements[0], 1);

        let one = one.into_large_uint();
        assert_eq!(one.elements[0], 1);
    }

    #[test]
    fn check_pow() {
        let mut two = GaloisElement::from_u64(2);
        let four = GaloisElement::from_u64(4);
        two.square();
        assert_eq!(two, four);

        let mut two = GaloisElement::from_u64(2);
        let two_c = GaloisElement::from_u64(2);
        let one = LargeUint::from_u64(1);
        two.pow(&one);
        assert_eq!(two, two_c);

        let mut two = GaloisElement::from_u64(2);
        let three = LargeUint::from_u64(3);
        let eight = GaloisElement::from_u64(8);
        two.pow(&three);
        assert_eq!(two, eight);
    }

    #[test]
    fn check_inv() {
        let mut two = GaloisElement::from_u64(2);
        let o_two = GaloisElement::from_u64(2);
        let one = GaloisElement::from_u64(1);
        two.inverse();
        two.mul_with(&o_two);
        assert_eq!(two, one);
    }

    #[test]
    fn check_square() {
        let four = GaloisElement::from_u64(4);
        let one = GaloisElement::from_u64(2);
        assert!(four.is_square());
        assert!(!one.is_square());
    }

    #[test]
    fn check_add_impl() {
        let one = GaloisElement::from_u64(1);
        let two = GaloisElement::from_u64(2);
        assert_eq!(one + one, two);
    }

    #[test]
    fn check_sub_impl() {
        let one = GaloisElement::from_u64(1);
        let two = GaloisElement::from_u64(2);
        assert_eq!(two - one, one);
    }

    #[test]
    fn check_mul_impl() {
        let two = GaloisElement::from_u64(2);
        let four = GaloisElement::from_u64(4);
        assert_eq!(two * two, four);
    }

    #[test]
    fn check_div_impl() {
        let one = GaloisElement::from_u64(1);
        let two = GaloisElement::from_u64(2);
        assert_eq!(two / one, two);

        let six = GaloisElement::from_u64(6);
        let two = GaloisElement::from_u64(2);
        let three = GaloisElement::from_u64(3);
        assert_eq!(six / two, three);
    }

    #[test]
    fn check_complex_impl() {
        let one = GaloisElement::from_u64(1);
        let two = GaloisElement::from_u64(2);
        assert_eq!(((one - two) * two + one + one + two)/ two, one);
    }

    #[test]
    fn check_parse() {
        let one = LargeUint::from_u64(1);
        let one_parsed = LargeUint::parse_bytes(b"1");
        assert_eq!(one, one_parsed);

        let foo = LargeUint::from_u64(12314123);
        let foo_parsed = LargeUint::parse_bytes(b"12314123");
        assert_eq!(foo, foo_parsed);
    }

    #[test]
    fn check_bits() {
        let one = LargeUint::from_u64(1);
        assert_eq!(one.bits(), 1);

        let one = LargeUint {
            elements: [0, 1, 1, 1, 1, 1, 0, 2]
        };

        assert_eq!(one.bits(), 7 * 64 + 2);
    }
}
