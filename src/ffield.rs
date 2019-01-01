use crate::global;
use num_bigint::BigUint;

pub fn is_square(num: &BigUint) -> bool {
    // Eulers Criterion
    let pow = num.modpow(&global::p_min_1_half, &global::p);
    return pow == 1u32.into();
}

pub fn inverse(num: &BigUint) -> BigUint {
    // Uses the fact that the a^(p-2) = 1 mod p
    num.modpow(&global::p_min_2, &global::p)
}

#[cfg(test)]
mod test {
    use num_bigint::BigUint;
    use num_integer::Integer;
    use crate::global;
    use super::{is_square, inverse};

    #[test]
    fn check_squares() {
        for i in 1u32..=16 {
            assert!(is_square(&(i * i).into()));
        }
    }

    #[test]
    fn check_inv() {
        for i in 1u32..=16 {
            let i = i.into();
            let a = inverse(&i);
            assert!((&a * &i).mod_floor(&*global::p) == BigUint::from(1u32));
        }
    }
}

