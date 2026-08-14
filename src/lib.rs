// rng-compat-r: reproduce R's default random-number behavior in Rust.
// Copyright (C) 2026 carbocation/rng-compat-r contributors
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.

//! R-compatible random-number generation.
//!
//! This crate implements the generator stacks used by R 3.5 through R 4.6:
//! Mersenne Twister and L'Ecuyer-CMRG uniforms, inversion normals, and both
//! rejection and pre-R-3.6 rounding samplers.
//!
//! The serialized state is compatible with R's `.Random.seed` integer vector.

use std::error::Error;
use std::fmt;

mod lecuyer;
mod math;

pub use math::{pnorm, pnorm_with_mode};

const MT_LEN: usize = 624;
const MT_M: usize = 397;
const MATRIX_A: u32 = 0x9908_b0df;
const UPPER_MASK: u32 = 0x8000_0000;
const LOWER_MASK: u32 = 0x7fff_ffff;
const INV_2_POW_32: f64 = 2.328_306_436_538_696_3e-10;
const INV_U32_MAX: f64 = 2.328_306_437_080_797e-10;
const MAX_R_SAMPLE_N: u64 = 4_500_000_000_000_000;

/// The uniform generator represented by an [`RRng`].
#[derive(Clone, Copy, Debug, Default, Eq, Hash, PartialEq)]
#[non_exhaustive]
pub enum RUniformKind {
    /// R's default MT19937 generator.
    #[default]
    MersenneTwister,
    /// Pierre L'Ecuyer's MRG32k3a generator used by R for parallel streams.
    LecuyerCmrg,
}

impl RUniformKind {
    const fn mode_code(self) -> i32 {
        match self {
            Self::MersenneTwister => 3,
            Self::LecuyerCmrg => 7,
        }
    }
}

/// Mathematical-library selection for normal and probability calculations.
#[derive(Clone, Copy, Debug, Default, Eq, Hash, PartialEq)]
#[non_exhaustive]
pub enum MathMode {
    /// Use the target platform's math library, like a locally compiled R.
    #[default]
    Platform,
    /// Use the crate's pinned pure-Rust math implementation.
    ///
    /// For finite inputs under IEEE-754 round-to-nearest semantics, this mode
    /// has bit-stable results across x86-64 and AArch64 and does not depend on
    /// FMA availability. NaN payload bits are not part of the contract.
    Deterministic,
}

/// R compatibility mode.
///
/// The uniform and normal algorithms are the same for every listed version.
/// R 3.6.0 changed uniform integer sampling from rounding to rejection.
#[derive(Clone, Copy, Debug, Default, Eq, Hash, PartialEq)]
#[non_exhaustive]
pub enum RVersion {
    /// R 3.5.x behavior, including the legacy rounding sampler.
    R3_5,
    /// R 3.6.x behavior, with rejection sampling.
    R3_6,
    /// R 4.0.x behavior, with rejection sampling.
    R4_0,
    /// R 4.5.x behavior, with rejection sampling.
    R4_5,
    /// R 4.6.x behavior, with rejection sampling.
    #[default]
    R4_6,
}

impl RVersion {
    const fn uses_rejection_sampling(self) -> bool {
        !matches!(self, Self::R3_5)
    }
}

/// An error while validating state or sampling arguments.
#[derive(Clone, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum RRngError {
    /// R reserves the minimum signed integer as `NA_integer_`.
    InvalidSeed(i32),
    /// `.Random.seed` did not contain the number of integers required by its generator.
    InvalidSeedLength { expected: usize, actual: usize },
    /// The encoded uniform generator is not supported.
    UnsupportedUniformKind(i32),
    /// The encoded normal generator is not inversion.
    UnsupportedNormalKind(i32),
    /// The encoded sample generator is not rounding or rejection.
    UnsupportedSampleKind(i32),
    /// The encoded binomial generator is unsupported.
    UnsupportedBinomialKind(i32),
    /// The Mersenne Twister cursor is outside R's valid range.
    InvalidPosition(i32),
    /// A required generator-state component contains no set bits.
    AllZeroState,
    /// A L'Ecuyer-CMRG state word was outside its component modulus.
    InvalidLecuyerState,
    /// A stream operation was requested for a non-L'Ecuyer generator.
    StreamOperationRequiresLecuyer,
    /// An output slice was too short for the serialized state.
    OutputTooSmall { required: usize, actual: usize },
    /// A population size was zero or exceeded R's supported range.
    InvalidPopulationSize(usize),
    /// Sampling without replacement requested more values than exist.
    SampleLargerThanPopulation { size: usize, population: usize },
}

impl fmt::Display for RRngError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            Self::InvalidSeed(seed) => write!(f, "{seed} is not a valid R integer seed"),
            Self::InvalidSeedLength { expected, actual } => {
                write!(
                    f,
                    "invalid .Random.seed length: expected {expected}, got {actual}"
                )
            }
            Self::UnsupportedUniformKind(kind) => {
                write!(f, "unsupported R uniform RNG kind {kind}")
            }
            Self::UnsupportedNormalKind(kind) => {
                write!(f, "unsupported R normal RNG kind {kind}")
            }
            Self::UnsupportedSampleKind(kind) => {
                write!(f, "unsupported R sample RNG kind {kind}")
            }
            Self::UnsupportedBinomialKind(kind) => {
                write!(f, "unsupported R binomial RNG kind {kind}")
            }
            Self::InvalidPosition(position) => {
                write!(f, "invalid Mersenne Twister position {position}")
            }
            Self::AllZeroState => f.write_str("a required generator-state component is all zero"),
            Self::InvalidLecuyerState => f.write_str("invalid L'Ecuyer-CMRG state"),
            Self::StreamOperationRequiresLecuyer => {
                f.write_str("RNG stream jumps require L'Ecuyer-CMRG")
            }
            Self::OutputTooSmall { required, actual } => write!(
                f,
                "state output is too short: need {required} integers, got {actual}"
            ),
            Self::InvalidPopulationSize(size) => {
                write!(f, "population size {size} is outside R's supported range")
            }
            Self::SampleLargerThanPopulation { size, population } => write!(
                f,
                "cannot sample {size} values without replacement from {population} values"
            ),
        }
    }
}

impl Error for RRngError {}

#[derive(Clone, Debug)]
#[allow(clippy::large_enum_variant)] // Keep MT state inline so construction and cloning do not allocate.
enum UniformState {
    MersenneTwister {
        words: [u32; MT_LEN],
        position: usize,
    },
    LecuyerCmrg(lecuyer::LecuyerState),
}

/// An R-compatible random-number generator and its serializable state.
///
/// `RRng` is [`Send`] and [`Sync`]. Random draws require `&mut self`, so native
/// parallel code should move independent generators to worker threads. A
/// shared mutable generator still requires synchronization.
#[derive(Clone, Debug)]
pub struct RRng {
    uniform: UniformState,
    version: RVersion,
    math_mode: MathMode,
}

impl RRng {
    /// Reproduce `set.seed(seed)` using R's historical seed scrambling.
    ///
    /// # Panics
    ///
    /// Panics for `i32::MIN`, which R reserves as `NA_integer_`. Use
    /// [`Self::try_from_seed`] to handle that case without panicking.
    #[must_use]
    pub fn from_seed(seed: i32) -> Self {
        Self::try_from_seed(seed).expect("i32::MIN is R's NA_integer_, not a valid seed")
    }

    /// Checked form of [`Self::from_seed`].
    ///
    /// # Errors
    ///
    /// Returns [`RRngError::InvalidSeed`] for `i32::MIN`.
    pub fn try_from_seed(seed: i32) -> Result<Self, RRngError> {
        Self::try_from_seed_with_kind(seed, RUniformKind::MersenneTwister)
    }

    /// Reproduce `RNGkind(kind); set.seed(seed)` for a supported generator.
    ///
    /// # Panics
    ///
    /// Panics for `i32::MIN`. Use [`Self::try_from_seed_with_kind`] for
    /// checked input.
    #[must_use]
    pub fn from_seed_with_kind(seed: i32, kind: RUniformKind) -> Self {
        Self::try_from_seed_with_kind(seed, kind)
            .expect("i32::MIN is R's NA_integer_, not a valid seed")
    }

    /// Checked form of [`Self::from_seed_with_kind`].
    pub fn try_from_seed_with_kind(seed: i32, kind: RUniformKind) -> Result<Self, RRngError> {
        if seed == i32::MIN {
            return Err(RRngError::InvalidSeed(seed));
        }
        let mut scrambled = seed as u32;
        for _ in 0..50 {
            scrambled = r_lcg(scrambled);
        }

        let uniform = match kind {
            RUniformKind::MersenneTwister => {
                // R initializes 625 words. The first is the cursor and is then
                // overwritten with 624 by FixupSeeds(); retain that step.
                scrambled = r_lcg(scrambled);
                let mut words = [0_u32; MT_LEN];
                for word in &mut words {
                    scrambled = r_lcg(scrambled);
                    *word = scrambled;
                }
                UniformState::MersenneTwister {
                    words,
                    position: MT_LEN,
                }
            }
            RUniformKind::LecuyerCmrg => {
                UniformState::LecuyerCmrg(lecuyer::LecuyerState::from_scrambled_seed(scrambled))
            }
        };

        Ok(Self {
            uniform,
            version: RVersion::default(),
            math_mode: MathMode::default(),
        })
    }

    /// Restore a supported Inversion-normal `.Random.seed` vector.
    ///
    /// Mersenne Twister and L'Ecuyer-CMRG uniform states are accepted.
    pub fn from_random_seed(seed: &[i32]) -> Result<Self, RRngError> {
        if seed.is_empty() {
            return Err(RRngError::InvalidSeedLength {
                expected: 1,
                actual: 0,
            });
        }

        let modes = seed[0];
        if modes < 0 {
            return Err(RRngError::UnsupportedUniformKind(modes));
        }

        let uniform_kind = modes % 100;
        let normal_kind = modes % 10_000 / 100;
        let sample_kind = modes % 100_000 / 10_000;
        let binomial_kind = modes / 100_000;

        if normal_kind != 4 {
            return Err(RRngError::UnsupportedNormalKind(normal_kind));
        }
        if !matches!(sample_kind, 0 | 1) {
            return Err(RRngError::UnsupportedSampleKind(sample_kind));
        }
        if binomial_kind != 0 {
            return Err(RRngError::UnsupportedBinomialKind(binomial_kind));
        }

        let uniform = match uniform_kind {
            3 => {
                const EXPECTED: usize = MT_LEN + 2;
                if seed.len() != EXPECTED {
                    return Err(RRngError::InvalidSeedLength {
                        expected: EXPECTED,
                        actual: seed.len(),
                    });
                }
                let position = seed[1];
                if !(1..=MT_LEN as i32).contains(&position) {
                    return Err(RRngError::InvalidPosition(position));
                }
                let mut words = [0_u32; MT_LEN];
                for (word, &serialized) in words.iter_mut().zip(&seed[2..]) {
                    *word = serialized as u32;
                }
                if words.iter().all(|&word| word == 0) {
                    return Err(RRngError::AllZeroState);
                }
                UniformState::MersenneTwister {
                    words,
                    position: position as usize,
                }
            }
            7 => {
                const EXPECTED: usize = lecuyer::STATE_LEN + 1;
                if seed.len() != EXPECTED {
                    return Err(RRngError::InvalidSeedLength {
                        expected: EXPECTED,
                        actual: seed.len(),
                    });
                }
                UniformState::LecuyerCmrg(lecuyer::LecuyerState::from_serialized(&seed[1..])?)
            }
            _ => return Err(RRngError::UnsupportedUniformKind(uniform_kind)),
        };

        Ok(Self {
            uniform,
            version: if sample_kind == 0 {
                RVersion::R3_5
            } else {
                RVersion::R4_6
            },
            math_mode: MathMode::default(),
        })
    }

    /// Select the algorithms associated with an R release.
    ///
    /// Changing the version does not re-seed the generator, matching the fact
    /// that these modes share the same uniform-generator initialization.
    #[must_use]
    pub const fn with_version(mut self, version: RVersion) -> Self {
        self.version = version;
        self
    }

    /// Change the compatibility mode without changing the uniform state.
    pub const fn set_version(&mut self, version: RVersion) {
        self.version = version;
    }

    /// Return the selected compatibility mode.
    #[must_use]
    pub const fn version(&self) -> RVersion {
        self.version
    }

    /// Select platform-native or cross-platform deterministic math.
    ///
    /// This setting is not encoded by R's `.Random.seed`; restoring a state
    /// therefore selects [`MathMode::Platform`] until explicitly changed.
    #[must_use]
    pub const fn with_math_mode(mut self, mode: MathMode) -> Self {
        self.math_mode = mode;
        self
    }

    /// Change mathematical evaluation without changing RNG state.
    pub const fn set_math_mode(&mut self, mode: MathMode) {
        self.math_mode = mode;
    }

    /// Return the selected mathematical evaluation mode.
    #[must_use]
    pub const fn math_mode(&self) -> MathMode {
        self.math_mode
    }

    /// Return the selected uniform generator.
    #[must_use]
    pub const fn uniform_kind(&self) -> RUniformKind {
        match self.uniform {
            UniformState::MersenneTwister { .. } => RUniformKind::MersenneTwister,
            UniformState::LecuyerCmrg(_) => RUniformKind::LecuyerCmrg,
        }
    }

    /// Return the exact number of integers needed by [`Self::write_random_seed`].
    #[must_use]
    pub const fn random_seed_len(&self) -> usize {
        match self.uniform {
            UniformState::MersenneTwister { .. } => MT_LEN + 2,
            UniformState::LecuyerCmrg(_) => lecuyer::STATE_LEN + 1,
        }
    }

    /// Serialize the state in R's `.Random.seed` representation.
    #[must_use]
    pub fn random_seed(&self) -> Vec<i32> {
        let mut result = vec![0; self.random_seed_len()];
        self.write_random_seed(&mut result)
            .expect("freshly allocated state has the exact required length");
        result
    }

    /// Write the R `.Random.seed` representation without allocating.
    ///
    /// Returns the number of integers written. Extra output capacity is left
    /// unchanged.
    pub fn write_random_seed(&self, output: &mut [i32]) -> Result<usize, RRngError> {
        let required = self.random_seed_len();
        if output.len() < required {
            return Err(RRngError::OutputTooSmall {
                required,
                actual: output.len(),
            });
        }
        let sample_mode = i32::from(self.version.uses_rejection_sampling());
        output[0] = 10_000 * sample_mode + 400 + self.uniform_kind().mode_code();
        match &self.uniform {
            UniformState::MersenneTwister { words, position } => {
                output[1] = *position as i32;
                for (target, &word) in output[2..required].iter_mut().zip(words) {
                    *target = word as i32;
                }
            }
            UniformState::LecuyerCmrg(state) => {
                state.write_serialized(&mut output[1..required]);
            }
        }
        Ok(required)
    }

    /// Return a copy advanced by R's `nextRNGStream` jump (`2^127` draws).
    pub fn next_rng_stream(&self) -> Result<Self, RRngError> {
        let UniformState::LecuyerCmrg(state) = &self.uniform else {
            return Err(RRngError::StreamOperationRequiresLecuyer);
        };
        let mut next = self.clone();
        next.uniform = UniformState::LecuyerCmrg(state.next_stream());
        Ok(next)
    }

    /// Return a copy advanced by R's `nextRNGSubStream` jump (`2^76` draws).
    pub fn next_rng_substream(&self) -> Result<Self, RRngError> {
        let UniformState::LecuyerCmrg(state) = &self.uniform else {
            return Err(RRngError::StreamOperationRequiresLecuyer);
        };
        let mut next = self.clone();
        next.uniform = UniformState::LecuyerCmrg(state.next_substream());
        Ok(next)
    }

    /// Generate one uniform value exactly as R's default `runif(1)` does.
    #[must_use]
    pub fn runif(&mut self) -> f64 {
        let value = match &mut self.uniform {
            UniformState::MersenneTwister { words, position } => {
                f64::from(mt_next_u32(words, position)) * INV_2_POW_32
            }
            UniformState::LecuyerCmrg(state) => state.next_uniform(),
        };
        if value <= 0.0 {
            0.5 * INV_U32_MAX
        } else if 1.0 - value <= 0.0 {
            1.0 - 0.5 * INV_U32_MAX
        } else {
            value
        }
    }

    /// Generate one uniform value over `[min, max)`, following `runif()`.
    #[must_use]
    pub fn runif_range(&mut self, min: f64, max: f64) -> f64 {
        if !min.is_finite() || !max.is_finite() || max < min {
            return f64::NAN;
        }
        if min == max {
            return min;
        }
        min + (max - min) * self.runif()
    }

    /// Generate one inversion-based normal variate as R's `rnorm()` does.
    #[must_use]
    pub fn rnorm(&mut self, mean: f64, sd: f64) -> f64 {
        if mean.is_nan() || sd.is_nan() || sd < 0.0 || sd == f64::INFINITY {
            return f64::NAN;
        }
        if sd == 0.0 || mean.is_infinite() {
            return mean;
        }

        const BIG: f64 = 134_217_728.0; // 2^27
        let upper = (BIG * self.runif()) as u32;
        let p = (f64::from(upper) + self.runif()) / BIG;
        mean + sd * math::qnorm_standard(p, self.math_mode)
    }

    /// Sample an idiomatic zero-based index from `0..population`.
    ///
    /// This consumes exactly the same uniforms as `sample.int(population, 1)`.
    /// Panics when `population` is zero or larger than R supports; use
    /// [`Self::try_sample_index`] for checked input.
    #[must_use]
    pub fn sample_index(&mut self, population: usize) -> usize {
        self.try_sample_index(population)
            .expect("sample_index population must be in 1..=4.5e15")
    }

    /// Checked form of [`Self::sample_index`].
    pub fn try_sample_index(&mut self, population: usize) -> Result<usize, RRngError> {
        Self::validate_nonempty_population(population)?;
        Ok(self.unif_index(population))
    }

    /// Reproduce uniform `sample.int(n, size, replace)`.
    ///
    /// Values are one-based, just as in R. The no-replacement algorithm uses
    /// R's partial Fisher-Yates procedure and therefore preserves RNG state.
    pub fn sample_int(
        &mut self,
        population: usize,
        size: usize,
        replace: bool,
    ) -> Result<Vec<usize>, RRngError> {
        if population as u64 > MAX_R_SAMPLE_N || (population == 0 && size > 0) {
            return Err(RRngError::InvalidPopulationSize(population));
        }
        if !replace && size > population {
            return Err(RRngError::SampleLargerThanPopulation { size, population });
        }

        if replace || size < 2 {
            return Ok((0..size).map(|_| self.unif_index(population) + 1).collect());
        }

        let mut values: Vec<usize> = (0..population).collect();
        let mut remaining = population;
        let mut result = Vec::with_capacity(size);
        for _ in 0..size {
            let index = self.unif_index(remaining);
            result.push(values[index] + 1);
            remaining -= 1;
            values[index] = values[remaining];
        }
        Ok(result)
    }

    /// Sample one-based values without replacement, matching R.
    pub fn sample_without_replacement(
        &mut self,
        population: usize,
        size: usize,
    ) -> Result<Vec<usize>, RRngError> {
        self.sample_int(population, size, false)
    }

    /// Return a one-based random permutation, matching `sample.int(n)`.
    pub fn permutation(&mut self, population: usize) -> Result<Vec<usize>, RRngError> {
        self.sample_int(population, population, false)
    }

    fn validate_nonempty_population(population: usize) -> Result<(), RRngError> {
        if population == 0 || population as u64 > MAX_R_SAMPLE_N {
            Err(RRngError::InvalidPopulationSize(population))
        } else {
            Ok(())
        }
    }

    fn unif_index(&mut self, population: usize) -> usize {
        if !self.version.uses_rejection_sampling() {
            return self.rounding_index(population);
        }

        let bits = usize::BITS - (population - 1).leading_zeros();
        loop {
            let candidate = self.random_bits(bits) as usize;
            if candidate < population {
                return candidate;
            }
        }
    }

    fn rounding_index(&mut self, population: usize) -> usize {
        let uniform = if population > i32::MAX as usize {
            const TWO_POW_25: f64 = 33_554_432.0;
            ((TWO_POW_25 * self.runif()).floor() + self.runif()) / TWO_POW_25
        } else {
            self.runif()
        };
        (population as f64 * uniform).floor() as usize
    }

    fn random_bits(&mut self, bits: u32) -> u64 {
        let mut value = 0_u64;
        let mut generated = 0;
        while generated <= bits {
            let chunk = (self.runif() * 65_536.0).floor() as u64;
            value = 65_536 * value + chunk;
            generated += 16;
        }
        value & ((1_u64 << bits) - 1)
    }
}

const fn r_lcg(seed: u32) -> u32 {
    seed.wrapping_mul(69_069).wrapping_add(1)
}

fn mt_next_u32(words: &mut [u32; MT_LEN], position: &mut usize) -> u32 {
    if *position >= MT_LEN {
        mt_twist(words);
        *position = 0;
    }

    let mut value = words[*position];
    *position += 1;
    value ^= value >> 11;
    value ^= (value << 7) & 0x9d2c_5680;
    value ^= (value << 15) & 0xefc6_0000;
    value ^= value >> 18;
    value
}

fn mt_twist(words: &mut [u32; MT_LEN]) {
    for index in 0..(MT_LEN - MT_M) {
        let combined = (words[index] & UPPER_MASK) | (words[index + 1] & LOWER_MASK);
        words[index] =
            words[index + MT_M] ^ (combined >> 1) ^ if combined & 1 == 0 { 0 } else { MATRIX_A };
    }
    for index in (MT_LEN - MT_M)..(MT_LEN - 1) {
        let combined = (words[index] & UPPER_MASK) | (words[index + 1] & LOWER_MASK);
        words[index] = words[index + MT_M - MT_LEN]
            ^ (combined >> 1)
            ^ if combined & 1 == 0 { 0 } else { MATRIX_A };
    }
    let combined = (words[MT_LEN - 1] & UPPER_MASK) | (words[0] & LOWER_MASK);
    words[MT_LEN - 1] =
        words[MT_M - 1] ^ (combined >> 1) ^ if combined & 1 == 0 { 0 } else { MATRIX_A };
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn set_seed_matches_r_initial_state_prefix() {
        let state = RRng::from_seed(42).random_seed();
        assert_eq!(state.len(), 626);
        assert_eq!(&state[..4], &[10_403, 624, 507_561_766, 1_260_545_903]);
    }

    #[test]
    fn save_restore_round_trip() {
        let mut original = RRng::from_seed(-123_456);
        for _ in 0..917 {
            let _ = original.runif();
        }
        let state = original.random_seed();
        let mut restored = RRng::from_random_seed(&state).unwrap();
        assert_eq!(restored.random_seed(), state);
        for _ in 0..1_000 {
            assert_eq!(restored.runif().to_bits(), original.runif().to_bits());
        }
    }

    #[test]
    fn save_restore_property_across_seeds_versions_and_operations() {
        let seeds = [i32::MIN + 1, -1_234_567, -1, 0, 1, 42, i32::MAX];
        let versions = [RVersion::R3_5, RVersion::R3_6, RVersion::R4_6];

        for seed in seeds {
            for version in versions {
                let mut original = RRng::from_seed(seed).with_version(version);
                for iteration in 0..37 {
                    let _ = original.permutation(2 + iteration % 23).unwrap();
                    let _ = original.rnorm(-3.5, 2.25);
                    let _ = original.runif();
                }

                let state = original.random_seed();
                let mut restored = RRng::from_random_seed(&state).unwrap();
                assert_eq!(restored.random_seed(), state);
                for _ in 0..128 {
                    assert_eq!(restored.runif().to_bits(), original.runif().to_bits());
                }
            }
        }
    }

    #[test]
    fn version_changes_only_sampler_mode() {
        let modern = RRng::from_seed(42);
        let legacy = modern.clone().with_version(RVersion::R3_5);
        assert_eq!(modern.random_seed()[0], 10_403);
        assert_eq!(legacy.random_seed()[0], 403);
        assert_eq!(&modern.random_seed()[1..], &legacy.random_seed()[1..]);
    }

    #[test]
    fn allocation_free_mt_export_preserves_extra_capacity() {
        let rng = RRng::from_seed(42);
        let mut output = [i32::MAX; MT_LEN + 4];
        assert_eq!(rng.write_random_seed(&mut output).unwrap(), MT_LEN + 2);
        assert_eq!(&output[..MT_LEN + 2], rng.random_seed());
        assert_eq!(&output[MT_LEN + 2..], &[i32::MAX; 2]);
        assert_eq!(
            rng.write_random_seed(&mut output[..MT_LEN + 1])
                .unwrap_err(),
            RRngError::OutputTooSmall {
                required: MT_LEN + 2,
                actual: MT_LEN + 1,
            }
        );
    }

    #[test]
    fn rng_is_send_and_sync() {
        const fn assert_send_sync<T: Send + Sync>() {}
        assert_send_sync::<RRng>();
    }

    #[test]
    fn r_edge_cases_do_not_consume_state() {
        assert_eq!(
            RRng::try_from_seed(i32::MIN).unwrap_err(),
            RRngError::InvalidSeed(i32::MIN)
        );

        let mut rng = RRng::from_seed(42);
        let initial = rng.random_seed();
        assert_eq!(rng.runif_range(2.0, 2.0), 2.0);
        assert_eq!(rng.rnorm(7.0, 0.0), 7.0);
        assert_eq!(rng.sample_int(0, 0, false).unwrap(), Vec::<usize>::new());
        assert_eq!(rng.random_seed(), initial);

        assert!(rng.runif_range(f64::NEG_INFINITY, 1.0).is_nan());
        assert!(rng.rnorm(0.0, f64::INFINITY).is_nan());
        assert_eq!(rng.random_seed(), initial);
    }

    #[test]
    fn imported_state_modes_are_validated() {
        let state = RRng::from_seed(42).random_seed();

        let mut invalid = state.clone();
        invalid[0] = 10_400;
        assert_eq!(
            RRng::from_random_seed(&invalid).unwrap_err(),
            RRngError::UnsupportedUniformKind(0)
        );

        let mut invalid = state.clone();
        invalid[1] = 0;
        assert_eq!(
            RRng::from_random_seed(&invalid).unwrap_err(),
            RRngError::InvalidPosition(0)
        );

        let mut invalid = state;
        invalid[2..].fill(0);
        assert_eq!(
            RRng::from_random_seed(&invalid).unwrap_err(),
            RRngError::AllZeroState
        );
    }
}
