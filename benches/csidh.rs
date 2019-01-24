#![feature(test)]
extern crate test;
use test::Bencher;
use csidh;

#[bench]
fn speed(b: &mut Bencher) {
    let mut rng = rand::thread_rng();
    let private = csidh::CsidhPrivateKey::generate_new(&mut rng);
    b.iter(|| {
        let public = private.get_public_key();
    });
}
