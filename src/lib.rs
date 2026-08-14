// rng-compat-r: reproduce R's default random-number behavior in Rust.
// Copyright (C) 2026 carbocation/rng-compat-r contributors
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.

//! R-compatible random-number generation.
//!
//! This crate implements the default generator stack used by modern R:
//! Mersenne Twister uniforms, inversion normals, and the rejection sampler
//! introduced in R 3.6.0. It can also select the pre-3.6 rounding sampler.
//!
//! The serialized state is compatible with R's `.Random.seed` integer vector.

use std::error::Error;
use std::fmt;

const MT_LEN: usize = 624;
const MT_M: usize = 397;
const MATRIX_A: u32 = 0x9908_b0df;
const UPPER_MASK: u32 = 0x8000_0000;
const LOWER_MASK: u32 = 0x7fff_ffff;
const INV_2_POW_32: f64 = 2.328_306_436_538_696_3e-10;
const INV_U32_MAX: f64 = 2.328_306_437_080_797e-10;
const MAX_R_SAMPLE_N: u64 = 4_500_000_000_000_000;

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
    /// `.Random.seed` did not contain the 626 integers required by MT19937.
    InvalidSeedLength { expected: usize, actual: usize },
    /// The encoded uniform generator is not Mersenne Twister.
    UnsupportedUniformKind(i32),
    /// The encoded normal generator is not inversion.
    UnsupportedNormalKind(i32),
    /// The encoded sample generator is not rounding or rejection.
    UnsupportedSampleKind(i32),
    /// The encoded binomial generator is unsupported.
    UnsupportedBinomialKind(i32),
    /// The Mersenne Twister cursor is outside R's valid range.
    InvalidPosition(i32),
    /// The Mersenne Twister state contains no set bits.
    AllZeroState,
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
            Self::AllZeroState => f.write_str("Mersenne Twister state is all zero"),
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

/// R's default random-number generator and its serializable state.
#[derive(Clone, Debug)]
pub struct RRng {
    mt: [u32; MT_LEN],
    position: usize,
    version: RVersion,
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
        if seed == i32::MIN {
            return Err(RRngError::InvalidSeed(seed));
        }
        let mut scrambled = seed as u32;

        for _ in 0..50 {
            scrambled = r_lcg(scrambled);
        }

        // R initializes 625 words here. The first is the cursor and is then
        // overwritten with 624 by FixupSeeds(); retain that discarded step.
        scrambled = r_lcg(scrambled);

        let mut mt = [0_u32; MT_LEN];
        for word in &mut mt {
            scrambled = r_lcg(scrambled);
            *word = scrambled;
        }

        Ok(Self {
            mt,
            position: MT_LEN,
            version: RVersion::default(),
        })
    }

    /// Restore an R Mersenne-Twister/Inversion `.Random.seed` vector.
    ///
    /// The vector must contain the mode word, cursor, and all 624 MT words.
    pub fn from_random_seed(seed: &[i32]) -> Result<Self, RRngError> {
        const EXPECTED: usize = MT_LEN + 2;
        if seed.len() != EXPECTED {
            return Err(RRngError::InvalidSeedLength {
                expected: EXPECTED,
                actual: seed.len(),
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

        if uniform_kind != 3 {
            return Err(RRngError::UnsupportedUniformKind(uniform_kind));
        }
        if normal_kind != 4 {
            return Err(RRngError::UnsupportedNormalKind(normal_kind));
        }
        if !matches!(sample_kind, 0 | 1) {
            return Err(RRngError::UnsupportedSampleKind(sample_kind));
        }
        if binomial_kind != 0 {
            return Err(RRngError::UnsupportedBinomialKind(binomial_kind));
        }

        let position = seed[1];
        if !(1..=MT_LEN as i32).contains(&position) {
            return Err(RRngError::InvalidPosition(position));
        }

        let mut mt = [0_u32; MT_LEN];
        for (word, &serialized) in mt.iter_mut().zip(&seed[2..]) {
            *word = serialized as u32;
        }
        if mt.iter().all(|&word| word == 0) {
            return Err(RRngError::AllZeroState);
        }

        Ok(Self {
            mt,
            position: position as usize,
            version: if sample_kind == 0 {
                RVersion::R3_5
            } else {
                RVersion::R4_6
            },
        })
    }

    /// Select the algorithms associated with an R release.
    ///
    /// Changing the version does not re-seed the generator, matching the fact
    /// that these modes share the same Mersenne Twister initialization.
    #[must_use]
    pub const fn with_version(mut self, version: RVersion) -> Self {
        self.version = version;
        self
    }

    /// Change the compatibility mode without changing the MT state.
    pub const fn set_version(&mut self, version: RVersion) {
        self.version = version;
    }

    /// Return the selected compatibility mode.
    #[must_use]
    pub const fn version(&self) -> RVersion {
        self.version
    }

    /// Serialize the state in R's `.Random.seed` representation.
    #[must_use]
    pub fn random_seed(&self) -> Vec<i32> {
        let mode = if self.version.uses_rejection_sampling() {
            10_403
        } else {
            403
        };
        let mut result = Vec::with_capacity(MT_LEN + 2);
        result.push(mode);
        result.push(self.position as i32);
        result.extend(self.mt.iter().map(|&word| word as i32));
        result
    }

    /// Generate one uniform value exactly as R's default `runif(1)` does.
    #[must_use]
    pub fn runif(&mut self) -> f64 {
        let value = f64::from(self.next_u32()) * INV_2_POW_32;
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
        mean + sd * qnorm_standard(p)
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

    fn next_u32(&mut self) -> u32 {
        if self.position >= MT_LEN {
            self.twist();
        }

        let mut value = self.mt[self.position];
        self.position += 1;
        value ^= value >> 11;
        value ^= (value << 7) & 0x9d2c_5680;
        value ^= (value << 15) & 0xefc6_0000;
        value ^= value >> 18;
        value
    }

    fn twist(&mut self) {
        for index in 0..(MT_LEN - MT_M) {
            let combined = (self.mt[index] & UPPER_MASK) | (self.mt[index + 1] & LOWER_MASK);
            self.mt[index] = self.mt[index + MT_M]
                ^ (combined >> 1)
                ^ if combined & 1 == 0 { 0 } else { MATRIX_A };
        }
        for index in (MT_LEN - MT_M)..(MT_LEN - 1) {
            let combined = (self.mt[index] & UPPER_MASK) | (self.mt[index + 1] & LOWER_MASK);
            self.mt[index] = self.mt[index + MT_M - MT_LEN]
                ^ (combined >> 1)
                ^ if combined & 1 == 0 { 0 } else { MATRIX_A };
        }
        let combined = (self.mt[MT_LEN - 1] & UPPER_MASK) | (self.mt[0] & LOWER_MASK);
        self.mt[MT_LEN - 1] =
            self.mt[MT_M - 1] ^ (combined >> 1) ^ if combined & 1 == 0 { 0 } else { MATRIX_A };
        self.position = 0;
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

// AS 241 as used by R's qnorm.c. rnorm's combined uniform never reaches the
// asymptotic r > 27 branch, so the two rational regions are sufficient here.
#[allow(clippy::excessive_precision)]
fn qnorm_standard(p: f64) -> f64 {
    if p == 0.0 {
        return f64::NEG_INFINITY;
    }
    if p == 1.0 {
        return f64::INFINITY;
    }
    if !(0.0..=1.0).contains(&p) || p.is_nan() {
        return f64::NAN;
    }

    let q = p - 0.5;
    if q.abs() <= 0.425 {
        let r = 0.180_625 - q * q;
        let numerator = horner(
            r,
            &[
                2509.0809287301226727,
                33430.575583588128105,
                67265.770927008700853,
                45921.953931549871457,
                13731.693765509461125,
                1971.5909503065514427,
                133.14166789178437745,
                3.387132872796366608,
            ],
        );
        let denominator = horner(
            r,
            &[
                5226.495278852854561,
                28729.085735721942674,
                39307.89580009271061,
                21213.794301586595867,
                5394.1960214247511077,
                687.1870074920579083,
                42.313330701600911252,
                1.0,
            ],
        );
        return q * numerator / denominator;
    }

    let tail = if q > 0.0 { (0.5 - p) + 0.5 } else { p };
    let mut r = (-tail.ln()).sqrt();
    let mut value = if r <= 5.0 {
        r -= 1.6;
        let numerator = horner(
            r,
            &[
                7.7454501427834140764e-4,
                0.0227238449892691845833,
                0.24178072517745061177,
                1.27045825245236838258,
                3.64784832476320460504,
                5.7694972214606914055,
                4.6303378461565452959,
                1.42343711074968357734,
            ],
        );
        let denominator = horner(
            r,
            &[
                1.05075007164441684324e-9,
                5.475938084995344946e-4,
                0.0151986665636164571966,
                0.14810397642748007459,
                0.68976733498510000455,
                1.6763848301838038494,
                2.05319162663775882187,
                1.0,
            ],
        );
        numerator / denominator
    } else {
        r -= 5.0;
        let numerator = horner(
            r,
            &[
                2.01033439929228813265e-7,
                2.71155556874348757815e-5,
                0.0012426609473880784386,
                0.026532189526576123093,
                0.29656057182850489123,
                1.7848265399172913358,
                5.4637849111641143699,
                6.6579046435011037772,
            ],
        );
        let denominator = horner(
            r,
            &[
                2.04426310338993978564e-15,
                1.4215117583164458887e-7,
                1.8463183175100546818e-5,
                7.868691311456132591e-4,
                0.0148753612908506148525,
                0.13692988092273580531,
                0.59983220655588793769,
                1.0,
            ],
        );
        numerator / denominator
    };

    if q < 0.0 {
        value = -value;
    }
    value
}

fn horner(x: f64, coefficients: &[f64]) -> f64 {
    let mut result = coefficients[0];
    for coefficient in &coefficients[1..] {
        result = result * x + coefficient;
    }
    result
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
