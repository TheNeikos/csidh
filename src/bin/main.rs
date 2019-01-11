use rand::distributions::{Distribution, Uniform};
use csidh::montgomery::Curve;
use csidh::global;
use csidh::csidh::action;

fn main() {
    let between = Uniform::from(-5..=5);
    let mut rng = rand::thread_rng();
    let mut a_secret = vec![0i8; global::NUM_PRIMES];

    for i in 0..global::NUM_PRIMES {
        a_secret[i] = between.sample(&mut rng);
    }

    println!("Alice's Secret Key: {:?}", a_secret);

    let curve = Curve::new(0u32.into(), 1u32.into());

    let alice_public = action(&curve, &a_secret[..]);

    println!("Alice's public key is: {}", alice_public);

    let mut b_secret = vec![0i8; global::NUM_PRIMES];

    for i in 0..global::NUM_PRIMES {
        b_secret[i] = between.sample(&mut rng);
    }

    println!("Bob's Secret Key: {:?}", b_secret);

    let bob_public = action(&curve, &b_secret[..]);

    println!("Bob's public key is: {}", bob_public);

    let alice_curve = Curve::new(alice_public, 1u32.into());
    let bob_curve = Curve::new(bob_public, 1u32.into());

    let alice_shared = action(&bob_curve, &a_secret[..]);
    let bob_shared = action(&alice_curve, &b_secret[..]);

    println!("Alice's Shared key: {}", alice_shared);
    println!("Bob's Shared key:   {}", bob_shared);

    if bob_shared == alice_shared {
        println!("========================");
        println!("========================");
        println!("=====SUCCESSSSS!!!!=====");
        println!("========================");
        println!("========================");
    } else {
        println!("Try again");
    }
}
