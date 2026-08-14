use std::env;
use std::path::PathBuf;
use std::process::Command;

use rng_compat_r::RVersion;

mod common;

const CASES: &[(&str, RVersion)] = &[
    ("3.5.0", RVersion::R3_5),
    ("3.6.0", RVersion::R3_6),
    ("4.0.0", RVersion::R4_0),
    ("4.5.0", RVersion::R4_5),
    ("4.6.0", RVersion::R4_6),
];

#[test]
fn compare_directly_with_installed_r() {
    let requested = env::var("RNG_COMPAT_R_VERSION").ok();
    let required = env::var_os("RNG_COMPAT_REQUIRE_R").is_some();
    let generator = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tools")
        .join("generate-fixture.R");

    for &(label, version) in CASES {
        if requested.as_deref().is_some_and(|target| target != label) {
            continue;
        }

        let output = match Command::new("Rscript")
            .arg("--vanilla")
            .arg(&generator)
            .arg(label)
            .output()
        {
            Ok(output) => output,
            Err(error) if !required => {
                eprintln!("skipping live R checks: {error}");
                return;
            }
            Err(error) => panic!("Rscript is required for this test: {error}"),
        };

        if !output.status.success() {
            if requested.is_none() {
                // An older local R cannot emulate a future RNGversion. It has
                // already checked every supported target encountered so far.
                continue;
            }
            panic!(
                "R fixture generator failed for {label}: {}",
                String::from_utf8_lossy(&output.stderr)
            );
        }

        let fixture = String::from_utf8(output.stdout).expect("R output must be UTF-8");
        common::assert_fixture(label, version, &fixture);
    }
}
