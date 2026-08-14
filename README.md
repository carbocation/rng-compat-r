# rng-compat-r

`rng-compat-r` is a pure Rust implementation of R's default random-number
behavior. Its purpose is exact cross-language reproducibility: matching values,
sampling decisions, and `.Random.seed` state without embedding R.

```rust
use rng_compat_r::{RRng, RVersion};

let mut rng = RRng::from_seed(42).with_version(RVersion::R4_5);

let uniform = rng.runif();
let normal = rng.rnorm(0.0, 1.0);
let zero_based_index = rng.sample_index(20);

let state: Vec<i32> = rng.random_seed();
let restored = RRng::from_random_seed(&state)?;
# Ok::<(), rng_compat_r::RRngError>(())
```

## Compatibility

The initial implementation covers:

- R's historical `set.seed()` scrambling and MT19937 initialization;
- default Mersenne Twister `runif()` output;
- default inversion-based `rnorm()` output;
- modern rejection sampling from R 3.6.0 onward;
- explicit compatibility modes for R 3.5, 3.6, 4.0, 4.5, and 4.6;
- the pre-3.6 rounding sampler for migration fixtures;
- `sample.int()` with and without replacement, including full permutations;
- `.Random.seed` serialization and restoration.

`sample_index()` is the idiomatic zero-based operation. `sample_int()`,
`sample_without_replacement()`, and `permutation()` return R-compatible one-based
values. A full permutation deliberately consumes one sampling operation per
element, matching `sample(seq(n))` rather than `sample.int(n, 1)`.

Only the default Mersenne-Twister/Inversion combinations are accepted on state
import. Weighted sampling and other R generator kinds are not yet implemented.

Normal results preserve R's AS 241 coefficients and evaluation order, but exact
floating-point contraction is compiler-dependent. The Rust results were
bit-identical to the installed R 4.5 build; the official R 4.6.1 macOS build,
compiled with a newer Apple toolchain, differed by at most five ULPs in the
100-value fixture. Normal generation still consumes exactly two uniforms, so
`.Random.seed`, sampled indices, and later RNG output remain exact. The test
suite enforces this small normal tolerance while requiring bit identity for
uniforms and exact equality for every serialized state and discrete result.

Golden fixtures cover R 3.5, 3.6, 4.0, 4.5, and 4.6. The historical fixtures
record their generator runtime and use R's `RNGversion()` compatibility mode;
the 4.6 fixture was generated directly with official R 4.6.1. CI additionally
runs the generator under R 3.5.3, 3.6.3, 4.0.5, and 4.6.1.

The Rye regression fixture follows the optimizer call order at Rye commit
`539b818c9b6010e65b63d829a6bf775c1d10f962`: full alpha permutation, normal,
full weight permutation, normal, then acceptance uniform. It intentionally does
not replace either full permutation with a single-index sample.

## Provenance and license

The algorithms are ports of the GPL-licensed implementations in R's
[`RNG.c`](https://github.com/wch/r-source/blob/trunk/src/main/RNG.c),
[`random.c`](https://github.com/wch/r-source/blob/trunk/src/main/random.c),
[`snorm.c`](https://github.com/wch/r-source/blob/trunk/src/nmath/snorm.c), and
[`qnorm.c`](https://github.com/wch/r-source/blob/trunk/src/nmath/qnorm.c).

This project is licensed under GPL-2.0-or-later. See `LICENSE`.
