// rng-compat-r: reproduce R's random-number behavior in Rust.
// Copyright (C) 2026 carbocation/rng-compat-r contributors
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.

// Pierre L'Ecuyer's MRG32k3a generator and R's stream jumps.

use crate::{r_lcg, RRngError};

pub(crate) const STATE_LEN: usize = 6;

const M1: u64 = 4_294_967_087;
const M2: u64 = 4_294_944_443;
const NORM: f64 = 2.328_306_549_295_728e-10;

const A1P76: [[u64; 3]; 3] = [
    [82_758_667, 1_871_391_091, 4_127_413_238],
    [3_672_831_523, 69_195_019, 1_871_391_091],
    [3_672_091_415, 3_528_743_235, 69_195_019],
];
const A2P76: [[u64; 3]; 3] = [
    [1_511_326_704, 3_759_209_742, 1_610_795_712],
    [4_292_754_251, 1_511_326_704, 3_889_917_532],
    [3_859_662_829, 4_292_754_251, 3_708_466_080],
];
const A1P127: [[u64; 3]; 3] = [
    [2_427_906_178, 3_580_155_704, 949_770_784],
    [226_153_695, 1_230_515_664, 3_580_155_704],
    [1_988_835_001, 986_791_581, 1_230_515_664],
];
const A2P127: [[u64; 3]; 3] = [
    [1_464_411_153, 277_697_599, 1_610_723_613],
    [32_183_930, 1_464_411_153, 1_022_607_788],
    [2_824_425_944, 32_183_930, 2_093_834_863],
];

#[derive(Clone, Debug)]
pub(crate) struct LecuyerState {
    words: [u32; STATE_LEN],
}

impl LecuyerState {
    pub(crate) fn from_scrambled_seed(mut seed: u32) -> Self {
        let mut words = [0_u32; STATE_LEN];
        for word in &mut words {
            loop {
                seed = r_lcg(seed);
                if u64::from(seed) < M2 {
                    *word = seed;
                    break;
                }
            }
        }
        Self { words }
    }

    pub(crate) fn from_serialized(serialized: &[i32]) -> Result<Self, RRngError> {
        debug_assert_eq!(serialized.len(), STATE_LEN);
        let mut words = [0_u32; STATE_LEN];
        for (word, &value) in words.iter_mut().zip(serialized) {
            *word = value as u32;
        }

        if words[..3].iter().all(|&word| word == 0) || words[3..].iter().all(|&word| word == 0) {
            return Err(RRngError::AllZeroState);
        }
        if words[..3].iter().any(|&word| u64::from(word) >= M1)
            || words[3..].iter().any(|&word| u64::from(word) >= M2)
        {
            return Err(RRngError::InvalidLecuyerState);
        }

        Ok(Self { words })
    }

    pub(crate) fn write_serialized(&self, output: &mut [i32]) {
        debug_assert_eq!(output.len(), STATE_LEN);
        for (target, &word) in output.iter_mut().zip(&self.words) {
            *target = word as i32;
        }
    }

    pub(crate) fn next_uniform(&mut self) -> f64 {
        let mut p1 =
            1_403_580_i64 * i64::from(self.words[1]) - 810_728_i64 * i64::from(self.words[0]);
        let mut k = p1 / M1 as i64;
        p1 -= k * M1 as i64;
        if p1 < 0 {
            p1 += M1 as i64;
        }
        self.words[0] = self.words[1];
        self.words[1] = self.words[2];
        self.words[2] = p1 as u32;

        let mut p2 =
            527_612_i64 * i64::from(self.words[5]) - 1_370_589_i64 * i64::from(self.words[3]);
        k = p2 / M2 as i64;
        p2 -= k * M2 as i64;
        if p2 < 0 {
            p2 += M2 as i64;
        }
        self.words[3] = self.words[4];
        self.words[4] = self.words[5];
        self.words[5] = p2 as u32;

        let difference = if p1 > p2 {
            p1 - p2
        } else {
            p1 - p2 + M1 as i64
        };
        difference as f64 * NORM
    }

    pub(crate) fn next_stream(&self) -> Self {
        self.jump(&A1P127, &A2P127)
    }

    pub(crate) fn next_substream(&self) -> Self {
        self.jump(&A1P76, &A2P76)
    }

    fn jump(&self, first: &[[u64; 3]; 3], second: &[[u64; 3]; 3]) -> Self {
        let mut words = [0_u32; STATE_LEN];
        multiply_mod(first, &self.words[..3], M1, &mut words[..3]);
        multiply_mod(second, &self.words[3..], M2, &mut words[3..]);
        Self { words }
    }
}

fn multiply_mod(matrix: &[[u64; 3]; 3], input: &[u32], modulus: u64, output: &mut [u32]) {
    for (row, target) in matrix.iter().zip(output) {
        let mut value = 0_u128;
        for (&coefficient, &word) in row.iter().zip(input) {
            value = (value + u128::from(coefficient) * u128::from(word)) % u128::from(modulus);
        }
        *target = value as u32;
    }
}
