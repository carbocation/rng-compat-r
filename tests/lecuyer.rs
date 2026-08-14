use rng_compat_r::{RRng, RRngError, RUniformKind};

mod common;

const FIXTURE: &str = include_str!("fixtures/lecuyer-r-4.6.0.txt");

#[test]
fn lecuyer_sequence_state_and_stream_jumps_match_r_4_6() {
    let fixture = common::Fixture::parse(FIXTURE);
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
            actual.abs_diff(expected) <= 5,
            "rnorm[{index}] differs by more than 5 ULPs: {actual:016x} != {expected:016x}"
        );
    }
    assert_eq!(rng.random_seed(), fixture.i32s("state_after_rnorm"));

    assert_eq!(rng.permutation(20).unwrap(), fixture.usizes("permutation"));
    assert_eq!(rng.random_seed(), fixture.i32s("state_after_permutation"));
}

#[test]
fn lecuyer_state_round_trip_and_allocation_free_export() {
    let mut rng = RRng::from_seed_with_kind(-765_432, RUniformKind::LecuyerCmrg);
    for _ in 0..257 {
        let _ = rng.rnorm(3.0, 0.25);
        let _ = rng.sample_index(31);
    }

    assert_eq!(rng.random_seed_len(), 7);
    let mut output = [i32::MAX; 9];
    assert_eq!(rng.write_random_seed(&mut output).unwrap(), 7);
    assert_eq!(&output[..7], rng.random_seed());
    assert_eq!(&output[7..], &[i32::MAX; 2]);

    let mut restored = RRng::from_random_seed(&output[..7]).unwrap();
    assert_eq!(restored.uniform_kind(), RUniformKind::LecuyerCmrg);
    for _ in 0..1_000 {
        assert_eq!(restored.runif().to_bits(), rng.runif().to_bits());
    }

    assert_eq!(
        rng.write_random_seed(&mut [0; 6]).unwrap_err(),
        RRngError::OutputTooSmall {
            required: 7,
            actual: 6
        }
    );
    assert_eq!(
        RRng::from_seed(42).next_rng_stream().unwrap_err(),
        RRngError::StreamOperationRequiresLecuyer
    );
}
