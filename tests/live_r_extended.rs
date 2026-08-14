use std::env;
use std::path::{Path, PathBuf};
use std::process::Command;

mod common;

#[test]
fn compare_pnorm_and_lecuyer_directly_with_r_4_6() {
    let tools = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tools");
    let Some(pnorm_fixture) = run_generator(&tools.join("generate-pnorm-fixture.R")) else {
        return;
    };
    common::assert_pnorm_fixture(&pnorm_fixture);

    let Some(lecuyer_fixture) = run_generator(&tools.join("generate-lecuyer-fixture.R")) else {
        return;
    };
    common::assert_lecuyer_fixture(&lecuyer_fixture);
}

fn run_generator(generator: &Path) -> Option<String> {
    let required = env::var_os("RNG_COMPAT_REQUIRE_R").is_some();
    let output = match Command::new("Rscript")
        .arg("--vanilla")
        .arg(generator)
        .output()
    {
        Ok(output) => output,
        Err(error) if !required => {
            eprintln!("skipping live R 4.6 checks: {error}");
            return None;
        }
        Err(error) => panic!("Rscript is required for this test: {error}"),
    };

    if !output.status.success() {
        if !required {
            eprintln!(
                "skipping live R 4.6 check from {}: {}",
                generator.display(),
                String::from_utf8_lossy(&output.stderr)
            );
            return None;
        }
        panic!(
            "R fixture generator {} failed: {}",
            generator.display(),
            String::from_utf8_lossy(&output.stderr)
        );
    }

    Some(String::from_utf8(output.stdout).expect("R output must be UTF-8"))
}
