use lazy_static::lazy_static;
use num_bigint::BigInt;

lazy_static! {
    pub static ref p: BigInt = {
        BigInt::parse_bytes(
            b"65b48e8f740f89bffc8ab0d15e3e4c4ab42d083aedc8\
            8c425afbfcc69322c9cda7aac6c567f35507516730cc1\
            f0b4f25c2721bf457aca8351b81b90533c6c87b",
            16
        ).unwrap()
    };

    pub static ref p_min_1_half: BigInt = {
        ((&*p) - BigInt::from(1u32)) / BigInt::from(2u32)
    };

    pub static ref p_min_2: BigInt = {
        (&*p) - BigInt::from(2u32)
    };

    pub static ref l: Vec<BigInt> = {
        let mut vec = vec![];
        let first_primes = primal::Primes::all().skip(1).take(73).map(Into::into);
        vec.extend(first_primes);
        vec.push(587u32.into());
        vec
    };
}
