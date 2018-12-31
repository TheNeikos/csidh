use crate::global;
use num_bigint::BigInt;

pub fn is_square(num: &BigInt) -> bool {
    // Eulers Criterion
    let pow = num.modpow(&global::p_min_1_half, &global::p);
    return pow == 1u32.into();
}

pub fn inverse(num: &BigInt) -> BigInt {
    // Uses the fact that the a^(p-2) = 1 mod p
    num.modpow(&global::p_min_2, &global::p)
}

#[cfg(test)]
mod test {
    use num::{BigInt, Integer};
    use crate::global;
    use super::{is_square, inverse};

    #[test]
    fn check_squares() {
        for i in 1..=16 {
            assert!(is_square(&(i * i).into()));
        }
    }

    #[test]
    fn check_inv() {
        for i in 1..=16 {
            let i = i.into();
            let a = inverse(&i);
            assert!((&a * &i).mod_floor(&*global::p) == BigInt::from(1));
        }
    }
}

