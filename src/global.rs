use crate::galois::{LargeUint, GaloisElement};

pub static PBITS: u64 = 511;

// -p^-1 mod 2^64
pub static INV_MIN_P_MOD_R: u64 = 0x66c1301f632e294d;

// 1
pub static LUINT_1: GaloisElement = GaloisElement {
    elements: [1, 0, 0, 0, 0, 0, 0 ,0]
};

pub static P_INT: LargeUint = LargeUint {
    elements: [
        0x1b81b90533c6c87b, 0xc2721bf457aca835, 0x516730cc1f0b4f25, 0xa7aac6c567f35507,
        0x5afbfcc69322c9cd, 0xb42d083aedc88c42, 0xfc8ab0d15e3e4c4a, 0x65b48e8f740f89bf
    ]
};

pub static P_MINUS_2: LargeUint = LargeUint {
    elements: [
        0x1b81b90533c6c879, 0xc2721bf457aca835, 0x516730cc1f0b4f25, 0xa7aac6c567f35507,
        0x5afbfcc69322c9cd, 0xb42d083aedc88c42, 0xfc8ab0d15e3e4c4a, 0x65b48e8f740f89bf,
    ]
};

pub static P_MINUS_1_HALVES: LargeUint = LargeUint {
    elements: [
        0x8dc0dc8299e3643d, 0xe1390dfa2bd6541a, 0xa8b398660f85a792, 0xd3d56362b3f9aa83,
        0x2d7dfe63499164e6, 0x5a16841d76e44621, 0xfe455868af1f2625, 0x32da4747ba07c4df,
    ]
};

pub static P: GaloisElement = GaloisElement {
    elements: [
        0x1b81b90533c6c87b, 0xc2721bf457aca835, 0x516730cc1f0b4f25, 0xa7aac6c567f35507,
        0x5afbfcc69322c9cd, 0xb42d083aedc88c42, 0xfc8ab0d15e3e4c4a, 0x65b48e8f740f89bf
    ]
};

// (2^512)^2 mod p
pub static R_SQUARED_MOD_P: GaloisElement = GaloisElement {
    elements: [
        0x36905b572ffc1724, 0x67086f4525f1f27d, 0x4faf3fbfd22370ca, 0x192ea214bcc584b1,
        0x5dae03ee2f5de3d0, 0x1e9248731776b371, 0xad5f166e20e4f52d, 0x4ed759aea6f3917e,
    ]
};

// 2^512 mod p
pub static GAL_1: GaloisElement = GaloisElement {
    elements: [
        0xc8fc8df598726f0a, 0x7b1bc81750a6af95, 0x5d319e67c1e961b4, 0xb0aa7275301955f1,
        0x4a080672d9ba6c64, 0x97a5ef8a246ee77b, 0x06ea9e5d4383676a, 0x3496e2e117e0ec80,
    ]
};

pub static LS: [u64; 74] = [
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
    101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193,
    197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307,
    311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 587
];

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn check_bits() {
        assert_eq!(PBITS, P.into_large_uint().bits());
    }
}

