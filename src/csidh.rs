use rand::prelude::*;

use crate::global;
use crate::montgomery::{Curve, ProjectivePoint};
use crate::galois::{LargeUint, GaloisElement};

use flamer::flame;

#[flame]
pub fn action(curve: &Curve, private: &[i32]) -> LargeUint {
    let mut curve = curve.clone();
    let mut private = Vec::from(private);
    let all_zero = |private: &[i32]| { private.iter().all(|&x| x == 0i32) };
    let mut rng = thread_rng();

    while !all_zero(&private) {
        println!("{:?}", private);
        let x = GaloisElement::random_element(&mut rng);
        let val = x * x * x + curve.a * x * x + x;
        let square_sign = if val.is_square() { 1i32 } else { -1i32 };

        let all_with_sign: Vec<usize> = private.iter().enumerate()
            .filter(|&(_, x)| x.signum() == square_sign)
            .map(|(i,_)| i )
            .collect();

        let mut k = all_with_sign.iter()
            .map(|&i| global::LS[i])
            .fold(LargeUint::from_u64(1), |mut acc, x| { acc.mul_with_u64(x); acc });

        if k == LargeUint::from_u64(0u64) {
            continue;
        }

        let p = ProjectivePoint::new(curve.clone(), x, GaloisElement::from_u64(1));

        let (mut p, _) = p.ladder(&(&*global::p_plus_1 / &k));

        for &i in all_with_sign.iter().rev() {
            println!("Checking: {}", i);

            let li = &global::LS[i];
            let k_div_li = &k / li;
            let point_k = p.ladder(&k_div_li).0;

            if point_k.is_infinity() {
                continue;
            }

            k = k_div_li;
            private[i] -= square_sign;
            let iso = curve.isogeny(*li, &point_k);
            curve.a = iso.get_a();
            p = iso.evaluate(&p);
            p.set_curve(curve.clone());
        }
    }

    curve.a.into_large_uint()
}
