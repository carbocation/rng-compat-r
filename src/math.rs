// rng-compat-r: reproduce R's random-number behavior in Rust.
// Copyright (C) 2026 carbocation/rng-compat-r contributors
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.

// Mathematical functions ported from R's nmath library.

use crate::MathMode;

#[derive(Clone, Copy)]
pub(crate) struct MathBackend(MathMode);

impl MathBackend {
    pub(crate) const fn new(mode: MathMode) -> Self {
        Self(mode)
    }

    fn log(self, value: f64) -> f64 {
        match self.0 {
            MathMode::Platform => value.ln(),
            MathMode::Deterministic => libm::log(value),
        }
    }

    fn exp(self, value: f64) -> f64 {
        match self.0 {
            MathMode::Platform => value.exp(),
            MathMode::Deterministic => libm::exp(value),
        }
    }

    fn log1p(self, value: f64) -> f64 {
        match self.0 {
            MathMode::Platform => value.ln_1p(),
            MathMode::Deterministic => libm::log1p(value),
        }
    }

    fn sqrt(self, value: f64) -> f64 {
        match self.0 {
            MathMode::Platform => value.sqrt(),
            MathMode::Deterministic => libm::sqrt(value),
        }
    }
}

/// Evaluate R's `pnorm5(x, mean, sd, lower_tail, log_p)` algorithm.
///
/// This uses the platform mathematical library, as a locally compiled R does.
/// Use [`pnorm_with_mode`] with [`MathMode::Deterministic`] when identical
/// finite results across supported architectures matter more than matching the
/// last few bits of a particular R binary.
#[must_use]
pub fn pnorm(x: f64, mean: f64, sd: f64, lower_tail: bool, log_p: bool) -> f64 {
    pnorm_with_mode(x, mean, sd, lower_tail, log_p, MathMode::Platform)
}

/// Evaluate R's `pnorm5` algorithm with an explicit math implementation.
#[must_use]
pub fn pnorm_with_mode(
    x: f64,
    mean: f64,
    sd: f64,
    lower_tail: bool,
    log_p: bool,
    mode: MathMode,
) -> f64 {
    if x.is_nan() || mean.is_nan() || sd.is_nan() {
        return x + mean + sd;
    }
    if !x.is_finite() && mean == x {
        return f64::NAN;
    }
    if sd <= 0.0 {
        if sd < 0.0 {
            return f64::NAN;
        }
        return boundary(x < mean, lower_tail, log_p);
    }

    let standardized = (x - mean) / sd;
    if !standardized.is_finite() {
        return boundary(x < mean, lower_tail, log_p);
    }

    pnorm_standard(standardized, lower_tail, log_p, MathBackend::new(mode))
}

fn boundary(left_of_mean: bool, lower_tail: bool, log_p: bool) -> f64 {
    let probability_is_zero = left_of_mean == lower_tail;
    match (probability_is_zero, log_p) {
        (true, true) => f64::NEG_INFINITY,
        (true, false) => 0.0,
        (false, true) => 0.0,
        (false, false) => 1.0,
    }
}

#[allow(clippy::excessive_precision)]
fn pnorm_standard(x: f64, lower_tail: bool, log_p: bool, math: MathBackend) -> f64 {
    const A: [f64; 5] = [
        2.2352520354606839287,
        161.02823106855587881,
        1067.6894854603709582,
        18154.981253343561249,
        0.065682337918207449113,
    ];
    const B: [f64; 4] = [
        47.20258190468824187,
        976.09855173777669322,
        10260.932208618978205,
        45507.789335026729956,
    ];
    const C: [f64; 9] = [
        0.39894151208813466764,
        8.8831497943883759412,
        93.506656132177855979,
        597.27027639480026226,
        2494.5375852903726711,
        6848.1904505362823326,
        11602.651437647350124,
        9842.7148383839780218,
        1.0765576773720192317e-8,
    ];
    const D: [f64; 8] = [
        22.266688044328115691,
        235.38790178262499861,
        1519.377599407554805,
        6485.558298266760755,
        18615.571640885098091,
        34900.952721145977266,
        38912.003286093271411,
        19685.429676859990727,
    ];
    const P: [f64; 6] = [
        0.21589853405795699,
        0.1274011611602473639,
        0.022235277870649807,
        0.001421619193227893466,
        2.9112874951168792e-5,
        0.02307344176494017303,
    ];
    const Q: [f64; 5] = [
        1.28426009614491121,
        0.468238212480865118,
        0.0659881378689285515,
        0.00378239633202758244,
        7.29751555083966205e-5,
    ];
    const SQRT_32: f64 = 5.6568542494923801952;
    const INV_SQRT_2_PI: f64 = 0.39894228040143267794;

    let y = x.abs();
    let (lower, upper) = if y <= 0.674_489_75 {
        let (xnum, xden) = if y > f64::EPSILON * 0.5 {
            let xsq = x * x;
            let mut xnum = A[4] * xsq;
            let mut xden = xsq;
            for index in 0..3 {
                xnum = (xnum + A[index]) * xsq;
                xden = (xden + B[index]) * xsq;
            }
            (xnum, xden)
        } else {
            (0.0, 0.0)
        };
        let temp = x * (xnum + A[3]) / (xden + B[3]);
        let lower = 0.5 + temp;
        let upper = 0.5 - temp;
        if log_p {
            (math.log(lower), math.log(upper))
        } else {
            (lower, upper)
        }
    } else if y <= SQRT_32 {
        let mut xnum = C[8] * y;
        let mut xden = y;
        for index in 0..7 {
            xnum = (xnum + C[index]) * y;
            xden = (xden + D[index]) * y;
        }
        let temp = (xnum + C[7]) / (xden + D[7]);
        pnorm_del(x, y, temp, log_p, math)
    } else if (log_p && y < 1.0e170)
        || (lower_tail && -38.4674 < x && x < 8.2924)
        || (!lower_tail && -8.2924 < x && x < 38.4674)
    {
        let xsq = 1.0 / (x * x);
        let mut xnum = P[5] * xsq;
        let mut xden = xsq;
        for index in 0..4 {
            xnum = (xnum + P[index]) * xsq;
            xden = (xden + Q[index]) * xsq;
        }
        let temp = (INV_SQRT_2_PI - xsq * (xnum + P[4]) / (xden + Q[4])) / y;
        pnorm_del(x, x, temp, log_p, math)
    } else if x > 0.0 {
        if log_p {
            (0.0, f64::NEG_INFINITY)
        } else {
            (1.0, 0.0)
        }
    } else if log_p {
        (f64::NEG_INFINITY, 0.0)
    } else {
        (0.0, 1.0)
    };

    if lower_tail {
        lower
    } else {
        upper
    }
}

fn pnorm_del(
    sign_source: f64,
    split_source: f64,
    temp: f64,
    log_p: bool,
    math: MathBackend,
) -> (f64, f64) {
    let xsq = (split_source * 16.0).trunc() * 0.0625;
    let del = (split_source - xsq) * (split_source + xsq);
    let (small, complement) = if log_p {
        let exponent = -xsq * (xsq * 0.5) - del * 0.5;
        let small = exponent + math.log(temp);
        let complement = math.log1p(-math.exp(-xsq * (xsq * 0.5)) * math.exp(-del * 0.5) * temp);
        (small, complement)
    } else {
        let small = math.exp(-xsq * (xsq * 0.5)) * math.exp(-del * 0.5) * temp;
        (small, 1.0 - small)
    };

    if sign_source > 0.0 {
        (complement, small)
    } else {
        (small, complement)
    }
}

// AS 241 as used by R's qnorm.c. rnorm's combined uniform never reaches the
// asymptotic r > 27 branch, so the two rational regions are sufficient here.
#[allow(clippy::excessive_precision)]
pub(crate) fn qnorm_standard(p: f64, mode: MathMode) -> f64 {
    if p == 0.0 {
        return f64::NEG_INFINITY;
    }
    if p == 1.0 {
        return f64::INFINITY;
    }
    if !(0.0..=1.0).contains(&p) || p.is_nan() {
        return f64::NAN;
    }

    let math = MathBackend::new(mode);
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
    let mut r = math.sqrt(-math.log(tail));
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
