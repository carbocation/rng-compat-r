#![allow(dead_code)]

use std::collections::HashMap;

use rng_compat_r::{RRng, RVersion};

const MAX_NORMAL_ULPS: u64 = 5;

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
