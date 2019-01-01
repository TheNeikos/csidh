use crate::global;
use rand::prelude::*;
use num_bigint::{BigUint, RandBigInt};

pub fn action(_a: BigUint, private: &[BigUint]) -> BigUint {
    let any_zero = || { private.iter().any(|x| x != &BigUint::from(0u32)) };
    let mut rng = thread_rng();

    while any_zero() {
        let _x = rng.gen_biguint_range(&BigUint::from(1u32), &global::p);
    }

    0u32.into()
}
