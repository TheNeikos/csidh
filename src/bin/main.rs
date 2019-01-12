use csidh::CsidhPrivateKey;

fn main() {
    let mut rng = rand::thread_rng();

    let a_private = CsidhPrivateKey::generate_new(&mut rng);
    let b_private = CsidhPrivateKey::generate_new(&mut rng);

    let a_public = a_private.get_public_key();
    let b_public = b_private.get_public_key();

    let a_shared = a_private.get_shared_secret(&b_public);
    let b_shared = b_private.get_shared_secret(&a_public);

    assert_eq!(a_shared, b_shared);
    println!("You have a common secret!: \n{:?}", a_shared);
}
