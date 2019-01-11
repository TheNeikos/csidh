use rand::prelude::*;

use crate::global;
use crate::montgomery::{Curve, ProjectivePoint};
use crate::galois::{LargeUint, GaloisElement};

pub fn action(curve: &Curve, private: &[i8]) -> LargeUint {
    let mut rng = thread_rng();

    let mut k = [LargeUint::from_u64(4) ;2];

    let mut e = [[0u8; global::NUM_PRIMES]; 2];

    for i in 0..global::NUM_PRIMES {
        let t = private[i];

        if t > 0 {
            e[0][i] = t as u8;
            e[1][i] = 0;
            k[1].mul_with_u64(global::PRIMES[i]);
        } else if t < 0 {
            e[1][i] = -t as u8;
            e[0][i] = 0;
            k[0].mul_with_u64(global::PRIMES[i]);
        } else {
            e[0][i] = 0;
            e[1][i] = 0;
            k[0].mul_with_u64(global::PRIMES[i]);
            k[1].mul_with_u64(global::PRIMES[i]);
        }
    }

    let mut p_curve = ProjectivePoint::new(curve.a, GaloisElement::from_u64(1));
    let mut done = [false; 2];

    loop {
        assert_eq!(p_curve.z, GaloisElement::from_u64(1));

        let x = GaloisElement::random_element(&mut rng);
        let sign = (!Curve::right_side(&p_curve.x, &x).is_square()) as usize;

        if done[sign] {
            continue;
        }

        let p = ProjectivePoint::new(x, GaloisElement::from_u64(1));

        let mut p = p.ladder2(&p_curve, &k[sign]);

        done[sign] = true;

        for i in (0..global::NUM_PRIMES).rev() {
            if e[sign][i] != 0 {
                let mut cof = LargeUint::from_u64(1);
                for j in 0..i {
                    if e[sign][j] != 0 {
                        cof.mul_with_u64(global::PRIMES[j]);
                    }
                }

                let kernel = p.ladder2(&p_curve, &cof);
                if !kernel.is_infinity() {
                    Curve::isogeny(&mut p_curve, &mut p, &kernel, global::PRIMES[i]);
                    e[sign][i] -= 1;
                    if e[sign][i] == 0 {
                        k[sign].mul_with_u64(global::PRIMES[i]);
                    }
                }
            }

            done[sign] &= e[sign][i] == 0;
        }

        p_curve.z.inverse();
        p_curve.x.mul_with(&p_curve.z);
        p_curve.z = GaloisElement::from_u64(1);

        if done[1] && done[0] {
            break;
        }
    }

    return p_curve.x.into_large_uint();
}
