#[macro_use]
extern crate cfg_if;
extern crate web_sys;
extern crate wasm_bindgen;
extern crate csidh;
extern crate rand;

use wasm_bindgen::prelude::*;
use csidh::{CsidhPrivateKey};

cfg_if! {
    // When the `console_error_panic_hook` feature is enabled, we can call the
    // `set_panic_hook` function to get better error messages if we ever panic.
    if #[cfg(feature = "console_error_panic_hook")] {
        extern crate console_error_panic_hook;
        use console_error_panic_hook::set_once as set_panic_hook;
    } else {
        #[inline]
        fn set_panic_hook() {}
    }
}

cfg_if! {
    // When the `wee_alloc` feature is enabled, use `wee_alloc` as the global
    // allocator.
    if #[cfg(feature = "wee_alloc")] {
        extern crate wee_alloc;
        #[global_allocator]
        static ALLOC: wee_alloc::WeeAlloc = wee_alloc::WeeAlloc::INIT;
    }
}

fn now() -> f64 {
    web_sys::window()
        .expect("should have a Window")
        .performance()
        .expect("should have a Performance")
        .now()
}

fn write_str(text: &str) -> Result<(), JsValue> {
    let window = web_sys::window().expect("should have a Window");
    let document = window.document().expect("should have a Document");
    let p: web_sys::Element = document.create_element("p")?.into();
    p.set_inner_html(&text);

    let body = document.body().expect("should have a body");
    let body: &web_sys::Node = body.as_ref();
    body.append_child(&p)?;
    Ok(())
}

// Called by our JS entry point to run the example.
#[wasm_bindgen]
pub fn run() -> Result<(), JsValue> {
    set_panic_hook();

    let mut rng = rand::thread_rng();

    let start = now();
    let a_private = CsidhPrivateKey::generate_new(&mut rng);
    let b_private = CsidhPrivateKey::generate_new(&mut rng);
    let end = now();

    write_str(&format!("Generated private keys in {}ms", end - start));

    let start = now();
    let a_public = a_private.get_public_key();
    let end = now();

    write_str(&format!("Generated Alice's Public Key in {}ms", end - start));

    let start = now();
    let b_public = b_private.get_public_key();
    let end = now();

    write_str(&format!("Generated Bob's Public Key in {}ms", end - start));

    let start = now();
    let a_shared = a_private.get_shared_secret(&b_public);
    let end = now();

    write_str(&format!("Generated Alice's Shared Secret in {}ms:<br><pre>{:?}</pre>", end - start, a_shared));

    let start = now();
    let b_shared = b_private.get_shared_secret(&a_public);
    let end = now();

    write_str(&format!("Generated Bob's Shared Secret in {}ms:<br><pre>{:?}</pre>", end - start, b_shared));

    assert_eq!(a_shared, b_shared);
    Ok(())
}
