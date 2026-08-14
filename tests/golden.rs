use rng_compat_r::RVersion;

mod common;

const CASES: &[(&str, RVersion, &str)] = &[
    (
        "3.5.0",
        RVersion::R3_5,
        include_str!("fixtures/r-3.5.0.txt"),
    ),
    (
        "3.6.0",
        RVersion::R3_6,
        include_str!("fixtures/r-3.6.0.txt"),
    ),
    (
        "4.0.0",
        RVersion::R4_0,
        include_str!("fixtures/r-4.0.0.txt"),
    ),
    (
        "4.5.0",
        RVersion::R4_5,
        include_str!("fixtures/r-4.5.0.txt"),
    ),
    (
        "4.6.0",
        RVersion::R4_6,
        include_str!("fixtures/r-4.6.0.txt"),
    ),
];

#[test]
fn complete_r_sequences_and_states_match() {
    for &(label, version, fixture_text) in CASES {
        common::assert_fixture(label, version, fixture_text);
    }
}
