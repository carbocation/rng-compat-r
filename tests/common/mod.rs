#![allow(dead_code)]

use std::collections::HashMap;

use rng_compat_r::{pnorm, RRng, RUniformKind, RVersion};

const MAX_NORMAL_ULPS: u64 = 5;
const MAX_PNORM_ULPS: u64 = 6;

pub fn assert_fixture(label: &str, version: RVersion, fixture_text: &str) {
    let fixture = Fixture::parse(fixture_text);
    assert_eq!(fixture.field("target_version"), label);

    let mut rng = RRng::from_seed(42).with_version(version);

    let actual_uniforms: Vec<_> = (0..100).map(|_| rng.runif().to_bits()).collect();
    assert_eq!(actual_uniforms, fixture.hex_u64("runif_bits"), "R {label}");
    assert_eq!(
        rng.random_seed(),
        fixture.i32s("state_after_runif"),
        "R {label}"
    );

    let actual_normals: Vec<_> = (0..100).map(|_| rng.rnorm(0.0, 1.0).to_bits()).collect();
    let expected_normals = fixture.hex_u64("rnorm_bits");
    for (index, (&actual, &expected)) in actual_normals.iter().zip(&expected_normals).enumerate() {
        let ulps = actual.abs_diff(expected);
        assert!(
            ulps <= MAX_NORMAL_ULPS,
            "R {label} rnorm[{index}] differs by {ulps} ULPs: {actual:016x} != {expected:016x}"
        );
    }
    assert_eq!(
        rng.random_seed(),
        fixture.i32s("state_after_rnorm"),
        "R {label}"
    );

    assert_eq!(
        rng.permutation(20).unwrap(),
        fixture.usizes("permutation"),
        "R {label}"
    );
    assert_eq!(
        rng.random_seed(),
        fixture.i32s("state_after_permutation"),
        "R {label}"
    );
}

pub fn assert_lecuyer_fixture(fixture_text: &str) {
    let fixture = Fixture::parse(fixture_text);
    assert_eq!(fixture.field("target_version"), "4.6.0");

    let mut rng = RRng::from_seed_with_kind(42, RUniformKind::LecuyerCmrg);
    assert_eq!(rng.random_seed(), fixture.i32s("initial_state"));
    assert_eq!(
        rng.next_rng_stream().unwrap().random_seed(),
        fixture.i32s("next_stream")
    );
    assert_eq!(
        rng.next_rng_substream().unwrap().random_seed(),
        fixture.i32s("next_substream")
    );

    let uniforms: Vec<_> = (0..100).map(|_| rng.runif().to_bits()).collect();
    assert_eq!(uniforms, fixture.hex_u64("runif_bits"));
    assert_eq!(rng.random_seed(), fixture.i32s("state_after_runif"));

    let normals: Vec<_> = (0..100).map(|_| rng.rnorm(0.0, 1.0).to_bits()).collect();
    for (index, (&actual, expected)) in normals
        .iter()
        .zip(fixture.hex_u64("rnorm_bits"))
        .enumerate()
    {
        assert!(
            actual.abs_diff(expected) <= MAX_NORMAL_ULPS,
            "rnorm[{index}] differs by more than {MAX_NORMAL_ULPS} ULPs: {actual:016x} != {expected:016x}"
        );
    }
    assert_eq!(rng.random_seed(), fixture.i32s("state_after_rnorm"));

    assert_eq!(rng.permutation(20).unwrap(), fixture.usizes("permutation"));
    assert_eq!(rng.random_seed(), fixture.i32s("state_after_permutation"));
}

pub fn assert_pnorm_fixture(fixture_text: &str) {
    let fixture = Fixture::parse(fixture_text);
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
    assert_pnorm_case(
        "standard_lower_plain_bits",
        standard.map(|x| pnorm(x, 0.0, 1.0, true, false)),
        &fixture,
    );
    assert_pnorm_case(
        "standard_upper_plain_bits",
        standard.map(|x| pnorm(x, 0.0, 1.0, false, false)),
        &fixture,
    );
    assert_pnorm_case(
        "standard_lower_log_bits",
        standard.map(|x| pnorm(x, 0.0, 1.0, true, true)),
        &fixture,
    );
    assert_pnorm_case(
        "standard_upper_log_bits",
        standard.map(|x| pnorm(x, 0.0, 1.0, false, true)),
        &fixture,
    );

    let shifted = [-100.0, -7.25, -1.0, 2.5, 6.0, 9.75, 100.0];
    assert_pnorm_case(
        "shifted_lower_plain_bits",
        shifted.map(|x| pnorm(x, 2.5, 3.25, true, false)),
        &fixture,
    );
    assert_pnorm_case(
        "shifted_upper_plain_bits",
        shifted.map(|x| pnorm(x, 2.5, 3.25, false, false)),
        &fixture,
    );
    assert_pnorm_case(
        "shifted_lower_log_bits",
        shifted.map(|x| pnorm(x, 2.5, 3.25, true, true)),
        &fixture,
    );
    assert_pnorm_case(
        "shifted_upper_log_bits",
        shifted.map(|x| pnorm(x, 2.5, 3.25, false, true)),
        &fixture,
    );
}

fn assert_pnorm_case<const N: usize>(name: &str, actual: [f64; N], fixture: &Fixture<'_>) {
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

pub struct Fixture<'a> {
    fields: HashMap<&'a str, &'a str>,
}

impl<'a> Fixture<'a> {
    pub fn parse(text: &'a str) -> Self {
        let fields = text
            .lines()
            .filter(|line| !line.is_empty() && !line.starts_with('#'))
            .map(|line| line.split_once('=').expect("fixture line must contain '='"))
            .collect();
        Self { fields }
    }

    pub fn field(&self, name: &str) -> &'a str {
        self.fields[name]
    }

    pub fn hex_u64(&self, name: &str) -> Vec<u64> {
        self.field(name)
            .split(',')
            .map(|value| u64::from_str_radix(value, 16).unwrap())
            .collect()
    }

    pub fn i32s(&self, name: &str) -> Vec<i32> {
        self.field(name)
            .split(',')
            .map(|value| value.parse().unwrap())
            .collect()
    }

    pub fn usizes(&self, name: &str) -> Vec<usize> {
        self.field(name)
            .split(',')
            .map(|value| value.parse().unwrap())
            .collect()
    }
}
