use rand::distributions::{Distribution, Uniform};
use csidh::montgomery::Curve;
use csidh::ffield::Field;
use csidh::global;
use csidh::csidh::action;

fn main() {
    let between = Uniform::from(-5..=5);
    let mut rng = rand::thread_rng();
    let mut a_secret = vec![0i32; 74];

    for i in 0..74 {
        a_secret[i] = between.sample(&mut rng);
    }

    println!("Alice's Secret Key: {:?}", a_secret);

    let field = Field::new(global::p.clone());
    let curve = Curve::new(field, 0u32.into(), 1u32.into());

    let alice_public = action(&curve, &a_secret[..]);

    println!("Alice's public key is: {}", alice_public);

    let mut b_secret = vec![0i32; 74];

    for i in 0..74 {
        b_secret[i] = between.sample(&mut rng);
    }

    println!("Bob's Secret Key: {:?}", b_secret);

    let field = Field::new(global::p.clone());
    let curve = Curve::new(field, 0u32.into(), 1u32.into());

    let bob_public = action(&curve, &b_secret[..]);

    println!("Bob's public key is: {}", bob_public);

    let field = Field::new(global::p.clone());
    let alices_curve = Curve::new(field.clone(), alice_public, 1u32.into());
    let bobs_curve = Curve::new(field, bob_public, 1u32.into());

    let bob_shared = action(&alices_curve, &b_secret[..]);
    let alice_shared = action(&bobs_curve, &a_secret[..]);

    println!("Bob's Shared key: {}", bob_shared);
    println!("Alice's Shared key: {}", alice_shared);
}
