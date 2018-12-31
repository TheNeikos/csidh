use crate::global;
use rand::prelude::*;
use num_bigint::{BigInt, RandBigInt};

pub fn action(a: BigInt, private: &[BigInt]) -> BigInt {
    let any_zero = || { private.iter().any(|x| x != &BigInt::from(0)) };
    let mut rng = thread_rng();

    while any_zero() {
        let x = rng.gen_bigint_range(&BigInt::from(1), &global::p);
    }

    0.into()
}
