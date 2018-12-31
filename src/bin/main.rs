use csidh::global;

fn main() {
    let p = &*global::p;
    println!("{}", p - 2u32 * p);
}
