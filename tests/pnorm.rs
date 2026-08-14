use rng_compat_r::pnorm;

mod common;

const FIXTURE: &str = include_str!("fixtures/pnorm-r-4.6.0.txt");

#[test]
fn pnorm_central_intermediate_and_extreme_tails_match_r_4_6() {
    common::assert_pnorm_fixture(FIXTURE);
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
