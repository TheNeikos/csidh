use rand::prelude::*;

use crate::global;
use crate::montgomery::{Curve, ProjectivePoint};
use crate::galois::{LargeUint, GaloisElement};

use flamer::flame;

#[flame]
pub fn action(curve: &Curve, private: &[i32]) -> LargeUint {
    let mut rng = thread_rng();

    let mut k = [LargeUint::from_u64(4) ;2];

    let mut e = [[0; global::NUM_PRIMES]; 2];

    for i in 0..global::NUM_PRIMES {
        let t = ((private[i/2] << i % 2 * 4) >> 4) as i8;

        if t > 0 {
            e[0][i] = t;
            e[1][i] = 0;
            k[1].mul_with_u64(global::PRIMES[i]);
        } else if t < 0 {
            e[0][i] = -t;
            e[1][i] = 0;
            k[0].mul_with_u64(global::PRIMES[i]);
        } else {
            e[0][i] = 0;
            e[1][i] = 0;
            k[0].mul_with_u64(global::PRIMES[i]);
            k[1].mul_with_u64(global::PRIMES[i]);
        }
    }

    let mut curve = *curve;
    let mut done = [false; 2];

    loop {
        let x = GaloisElement::random_element(&mut rng);
        let sign = curve.right_side(&x).is_square() as usize;

        if done[sign] {
            continue;
        }

        let p = ProjectivePoint::new(curve, x, GaloisElement::from_u64(1));

        let (mut p, _) = p.ladder(&k[sign]);

        done[sign] = true;

        for i in (0..global::NUM_PRIMES).rev() {
            println!("Checking {}", i);
            if e[sign][i] != 0 {
                let mut cof = LargeUint::from_u64(1);
                for j in 0..i {
                    if e[sign][j] != 0 {
                        cof.mul_with_u64(global::PRIMES[j]);
                    }
                }

                let (kernel, _) = p.ladder(&cof);
                if !kernel.is_infinity() {
                    let (new_a, new_p) = curve.isogeny(global::PRIMES[i], &kernel, &p);
                    curve.a = new_a;
                    p = new_p;
                    e[sign][i] -= 1;
                    if e[sign][i] == 0 {
                        k[sign].mul_with_u64(global::PRIMES[i]);
                    }
                }
            }

            done[sign] &= e[sign][i] == 0;
        }

        if done[1] && done[0] {
            break;
        }
    }

    return curve.a.into_large_uint();;

    // while !all_zero(&private) {
    //     println!("{:?}", private);
    //     let x = GaloisElement::random_element(&mut rng);
    //     let val = x * x * x + curve.a * x * x + x;
    //     let square_sign = if val.is_square() { 1i32 } else { -1i32 };

    //     let all_with_sign: Vec<usize> = private.iter().enumerate()
    //         .filter(|&(_, x)| x.signum() == square_sign)
    //         .map(|(i,_)| i )
    //         .collect();

    //     let mut k = all_with_sign.iter()
    //         .map(|&i| global::PRIMES[i])
    //         .fold(LargeUint::from_u64(1), |mut acc, x| { acc.mul_with_u64(x); acc });

    //     if k == LargeUint::from_u64(0u64) {
    //         continue;
    //     }

    //     let p = ProjectivePoint::new(curve, x, GaloisElement::from_u64(1));

    //     let (mut p, _) = p.ladder(&(&*global::p_plus_1 / &k));

    //     for &i in all_with_sign.iter().rev() {
    //         println!("Checking: {}", i);

    //         let li = &global::PRIMES[i];
    //         let k_div_li = &k / li;
    //         let point_k = p.ladder(&k_div_li).0;

    //         if point_k.is_infinity() {
    //             continue;
    //         }

    //         k = k_div_li;
    //         private[i] -= square_sign;
    //         let iso = curve.isogeny(*li, &point_k);
    //         curve.a = iso.get_a();
    //         p = iso.evaluate(&p);
    //         p.set_curve(curve);
    //     }
    // }

    // curve.a.into_large_uint()
}
