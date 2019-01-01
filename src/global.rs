use lazy_static::lazy_static;
use num_bigint::BigUint;
use num_integer::Integer;

lazy_static! {
    pub static ref p: BigUint = {
        BigUint::parse_bytes( b"65b48e8f740f89bffc8ab0d15e3e4c4ab42d083aedc88c425afbfcc69322c9cda7aac6c567f35507516730cc1f0b4f25c2721bf457aca8351b81b90533c6c87b", 16).unwrap()
    };

    pub static ref p_min_1_half: BigUint = {
        ((&*p) - BigUint::from(1u32)) / BigUint::from(2u32)
    };

    pub static ref p_min_2: BigUint = {
        (&*p) - BigUint::from(2u32)
    };

    pub static ref inv_4: BigUint = {
        BigUint::from(4u32).modpow(&p_min_2, &p)
    };

    pub static ref l: Vec<BigUint> = {
        let mut vec = vec![];
        let first_primes = primal::Primes::all().skip(1).take(73).map(Into::into);
        vec.extend(first_primes);
        vec.push(587u32.into());
        vec
    };
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn check_globals() {
        let other_p = BigUint::from(4u32) * l.iter().product::<BigUint>() - BigUint::from(1u32);
        assert!(*p == other_p);

        assert_eq!((&*inv_4 * BigUint::from(4u32)).mod_floor(&p), BigUint::from(1u32));
    }
}
