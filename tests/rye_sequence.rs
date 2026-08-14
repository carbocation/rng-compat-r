use rng_compat_r::{RRng, RVersion};

mod common;

const FIXTURE: &str = include_str!("fixtures/rye-r-4.6.0.txt");

#[test]
fn rye_optimizer_rng_call_order_matches_r() {
    let fixture = common::Fixture::parse(FIXTURE);
    assert_eq!(
        fixture.field("source_commit"),
        "2541fe03a163da92ce9279d87bf75a089b5ebf60"
    );

    let mut rng = RRng::from_seed(42).with_version(RVersion::R4_6);
    let mut alpha_indices = Vec::with_capacity(50);
    let mut alpha_normals = Vec::with_capacity(50);
    let mut weight_indices = Vec::with_capacity(50);
    let mut weight_normals = Vec::with_capacity(50);
    let mut acceptance_uniforms = Vec::with_capacity(50);

    for _ in 0..50 {
        alpha_indices.push(rng.sample_index(7) + 1);
        alpha_normals.push(rng.rnorm(0.0, 1.0).to_bits());
        weight_indices.push(rng.sample_index(11) + 1);
        weight_normals.push(rng.rnorm(0.0, 1.0).to_bits());
        // Rye calls pnorm() immediately before this; it consumes no RNG state.
        acceptance_uniforms.push(rng.runif().to_bits());
    }

    assert_eq!(alpha_indices, fixture.usizes("alpha_indices"));
    assert_normals_close(&alpha_normals, &fixture.hex_u64("alpha_normal_bits"));
    assert_eq!(weight_indices, fixture.usizes("weight_indices"));
    assert_normals_close(&weight_normals, &fixture.hex_u64("weight_normal_bits"));
    assert_eq!(
        acceptance_uniforms,
        fixture.hex_u64("acceptance_uniform_bits")
    );
    assert_eq!(rng.random_seed(), fixture.i32s("final_state"));
}

fn assert_normals_close(actual: &[u64], expected: &[u64]) {
    assert_eq!(actual.len(), expected.len());
    for (index, (&actual, &expected)) in actual.iter().zip(expected).enumerate() {
        assert!(
            actual.abs_diff(expected) <= 5,
            "normal {index} differs by more than 5 ULPs: {actual:016x} != {expected:016x}"
        );
    }
}
