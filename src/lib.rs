//! # A Rust implementation of the CSIDH Algorithm
//!
//! *DO NOT USE THIS CODE IN ANY WAY THAT REQUIRES IT TO BE CRYPTOGRAPHICALLY SECURE*
//!
//! ## Usage
//!
//! The CSIDH (pronounced sea-side) algorithm is commutative, this means that you can use it to
//! compute a shared secret over an insecure channel with a Diffie-Helmann key exchange. The
//! following snippet shows the basic gist of it.
//!
//! ```rust,no_run
//! # use csidh::CsidhPrivateKey;
//! let mut rng = rand::thread_rng();
//!
//! let a_private = CsidhPrivateKey::generate_new(&mut rng);
//! let b_private = CsidhPrivateKey::generate_new(&mut rng);
//!
//! let a_public = a_private.get_public_key();
//! let b_public = b_private.get_public_key();
//!
//! let a_shared = a_private.get_shared_secret(&b_public);
//! let b_shared = b_private.get_shared_secret(&a_public);
//!
//! assert_eq!(a_shared, b_shared);
//! ```
//!

mod global;
mod galois;
mod csidh;
mod montgomery;

pub use crate::csidh::{CsidhPrivateKey, CsidhPublicKey};
