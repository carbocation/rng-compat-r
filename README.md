# rng-compat-r

`rng-compat-r` is a pure Rust implementation of R's random-number behavior.
Its purpose is exact cross-language reproducibility: matching values, sampling
decisions, stream jumps, and `.Random.seed` state without embedding R.

```rust
use rng_compat_r::{pnorm, MathMode, RRng, RVersion};

let mut rng = RRng::from_seed(42)
    .with_version(RVersion::R4_6)
    .with_math_mode(MathMode::Deterministic);

let uniform = rng.runif();
let normal = rng.rnorm(0.0, 1.0);
let zero_based_index = rng.sample_index(20);
let probability = pnorm(normal, 0.0, 1.0, true, false);

let state: Vec<i32> = rng.random_seed();
let restored = RRng::from_random_seed(&state)?;
# Ok::<(), rng_compat_r::RRngError>(())
```

## Compatibility

The crate covers:

- R's historical `set.seed()` scrambling and MT19937 initialization;
- default Mersenne Twister `runif()` output;
- default inversion-based `rnorm()` output;
- R-compatible `pnorm(x, mean, sd, lower_tail, log_p)`, including logged tails
  and R's non-finite and zero-scale behavior;
- modern rejection sampling from R 3.6.0 onward;
- explicit compatibility modes for R 3.5, 3.6, 4.0, 4.5, and 4.6;
- the pre-3.6 rounding sampler for migration fixtures;
- `sample.int()` with and without replacement, including full permutations;
- Mersenne Twister and L'Ecuyer-CMRG `.Random.seed` serialization and
  restoration;
- R-compatible `nextRNGStream` and `nextRNGSubStream` jumps;
- allocation-free `.Random.seed` export through `write_random_seed()`.

`sample_index()` is the idiomatic zero-based operation. `sample_int()`,
`sample_without_replacement()`, and `permutation()` return R-compatible one-based
values. A full permutation deliberately consumes one sampling operation per
element, matching `sample(seq(n))` rather than `sample.int(n, 1)`.

Mersenne-Twister/Inversion and L'Ecuyer-CMRG/Inversion combinations are accepted
on state import. Weighted sampling and other R generator kinds are not
implemented.

## L'Ecuyer-CMRG streams

```rust
use rng_compat_r::{RRng, RUniformKind};

let stream = RRng::from_seed_with_kind(42, RUniformKind::LecuyerCmrg);
let next_stream = stream.next_rng_stream()?;
let next_substream = stream.next_rng_substream()?;

let mut state = [0_i32; 7];
let written = next_stream.write_random_seed(&mut state)?;
assert_eq!(written, 7);
# Ok::<(), rng_compat_r::RRngError>(())
```

The jump methods follow R's `parallel::nextRNGStream()` (`2^127` draws) and
`parallel::nextRNGSubStream()` (`2^76` draws). They return a new generator and
leave the source stream unchanged, matching the state-to-state behavior of the
R functions.

## Floating-point modes

Normal results preserve R's AS 241 coefficients and evaluation order, but exact
floating-point results can differ between separately compiled R binaries.
`MathMode::Platform`, the default, uses the target's math library and is the
closest match to a locally compiled R. The Rust results are bit-identical to the
installed R 4.5 build; the official R 4.6.1 macOS build differs by at most five
ULPs in the 100-normal fixtures. The R 4.6.1 `pnorm` fixtures allow six ULPs for
the same compiler and math-library reason.

`MathMode::Deterministic` uses the exactly pinned pure-Rust `libm` 0.2.16
implementation for `log`, `sqrt`, `exp`, and `log1p`. For finite results under
IEEE-754 round-to-nearest semantics, its bit patterns are part of the
`rng-compat-r` 0.1 compatibility contract across x86-64 (including AVX2,
AVX-512, and FMA-enabled targets), AArch64, and Rust compiler upgrades. Rust
floating-point expressions do not implicitly contract into FMA instructions;
the pinned golden tests detect any compiler change that would alter a result.
NaN payload bits and processes that modify the floating-point rounding or
flush-to-zero environment are outside this contract.

The selected math mode is deliberately not encoded in `.Random.seed`, because
R has no corresponding field. Restoring a state selects `MathMode::Platform`;
call `with_math_mode(MathMode::Deterministic)` again when restoring a
deterministic computation. Both modes consume exactly two uniforms per normal,
so serialized state, sampled indices, and later RNG output remain exact even
when a normal differs by a few ULPs.

Golden fixtures cover R 3.5, 3.6, 4.0, 4.5, and 4.6. The historical fixtures
record their generator runtime and use R's `RNGversion()` compatibility mode;
the 4.6 fixture was generated directly with official R 4.6.1. CI additionally
runs the generator under R 3.5.3, 3.6.3, 4.0.5, and 4.6.1.

R 4.6.1 fixtures separately cover `pnorm`, L'Ecuyer-CMRG state and stream
jumps, and deterministic-math bit patterns. Uniforms, sampled indices,
permutations, stream states, and every serialized state are compared exactly.

## Provenance and license

The algorithms are ports of the GPL-licensed implementations in R's
[`RNG.c`](https://github.com/wch/r-source/blob/trunk/src/main/RNG.c),
[`random.c`](https://github.com/wch/r-source/blob/trunk/src/main/random.c),
[`snorm.c`](https://github.com/wch/r-source/blob/trunk/src/nmath/snorm.c), and
[`qnorm.c`](https://github.com/wch/r-source/blob/trunk/src/nmath/qnorm.c),
[`pnorm.c`](https://github.com/wch/r-source/blob/trunk/src/nmath/pnorm.c), and
[`rngstream.c`](https://github.com/wch/r-source/blob/trunk/src/library/parallel/src/rngstream.c).

This project is licensed under GPL-2.0-or-later. See `LICENSE`.
