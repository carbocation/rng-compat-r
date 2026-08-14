use rng_compat_r::pnorm;

mod common;

const FIXTURE: &str = include_str!("fixtures/pnorm-r-4.6.0.txt");
const MAX_PNORM_ULPS: u64 = 6;

#[test]
fn pnorm_central_intermediate_and_extreme_tails_match_r_4_6() {
    let fixture = common::Fixture::parse(FIXTURE);
    assert_eq!(fixture.field("target_version"), "4.6.0");
    let sqrt_32 = 32.0_f64.sqrt();
    let standard = [
        -1e100,
        -50.0,
        -40.0,
        -38.5,
        -38.4674,
        -38.0,
        -10.0,
        -8.2924,
        -8.0,
        -sqrt_32,
        -1.0,
        -0.674_489_75,
        -0.1,
        0.0,
        0.1,
        0.674_489_75,
        1.0,
        sqrt_32,
        8.0,
        8.2924,
        10.0,
        38.0,
        38.4674,
        38.5,
        40.0,
        50.0,
        1e100,
    ];
    assert_case(
        "standard_lower_plain_bits",
        standard.map(|x| pnorm(x, 0.0, 1.0, true, false)),
        &fixture,
    );
    assert_case(
        "standard_upper_plain_bits",
        standard.map(|x| pnorm(x, 0.0, 1.0, false, false)),
        &fixture,
    );
    assert_case(
        "standard_lower_log_bits",
        standard.map(|x| pnorm(x, 0.0, 1.0, true, true)),
        &fixture,
    );
    assert_case(
        "standard_upper_log_bits",
        standard.map(|x| pnorm(x, 0.0, 1.0, false, true)),
        &fixture,
    );

    let shifted = [-100.0, -7.25, -1.0, 2.5, 6.0, 9.75, 100.0];
    assert_case(
        "shifted_lower_plain_bits",
        shifted.map(|x| pnorm(x, 2.5, 3.25, true, false)),
        &fixture,
    );
    assert_case(
        "shifted_upper_plain_bits",
        shifted.map(|x| pnorm(x, 2.5, 3.25, false, false)),
        &fixture,
    );
    assert_case(
        "shifted_lower_log_bits",
        shifted.map(|x| pnorm(x, 2.5, 3.25, true, true)),
        &fixture,
    );
    assert_case(
        "shifted_upper_log_bits",
        shifted.map(|x| pnorm(x, 2.5, 3.25, false, true)),
        &fixture,
    );
}

#[test]
fn pnorm_matches_r_edge_case_semantics() {
    assert!(pnorm(f64::NAN, 0.0, 1.0, true, false).is_nan());
    assert!(pnorm(0.0, f64::NAN, 1.0, true, false).is_nan());
    assert!(pnorm(0.0, 0.0, f64::NAN, true, false).is_nan());
    assert!(pnorm(0.0, 0.0, -1.0, true, false).is_nan());
    assert!(pnorm(f64::INFINITY, f64::INFINITY, 1.0, true, false).is_nan());

    assert_eq!(pnorm(-1.0, 0.0, 0.0, true, false), 0.0);
    assert_eq!(pnorm(0.0, 0.0, 0.0, true, false), 1.0);
    assert_eq!(pnorm(-1.0, 0.0, 0.0, false, false), 1.0);
    assert_eq!(pnorm(0.0, 0.0, 0.0, false, true), f64::NEG_INFINITY);
    assert_eq!(pnorm(7.0, 0.0, f64::INFINITY, true, false), 0.5);
    assert_eq!(
        pnorm(f64::NEG_INFINITY, 0.0, 1.0, true, true),
        f64::NEG_INFINITY
    );
    assert_eq!(pnorm(f64::INFINITY, 0.0, 1.0, false, false), 0.0);
}

fn assert_case<const N: usize>(name: &str, actual: [f64; N], fixture: &common::Fixture<'_>) {
    let expected = fixture.hex_u64(name);
    assert_eq!(actual.len(), expected.len(), "{name}");
    for (index, (&actual, &expected)) in actual.iter().zip(&expected).enumerate() {
        let actual = actual.to_bits();
        assert!(
            actual == expected || actual.abs_diff(expected) <= MAX_PNORM_ULPS,
            "{name}[{index}] differs by more than {MAX_PNORM_ULPS} ULPs: {actual:016x} != {expected:016x}"
        );
    }
}
