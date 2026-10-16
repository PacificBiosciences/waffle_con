//! Fast gate for performance work. It runs the same cases as `benches/consensus_bench.rs`.
//!
//! * Correctness: each case runs once and its `Debug` output is hashed.
//!   A speed-up must not change the `total digest` line.
//! * Speed: each case repeats for about `BUDGET_SECS` (at least `MIN_REPS` times).
//!   The minimum wall time is reported, because it is the least noisy estimate.
//!
//! Usage: `cargo run --release --example perf_gate [-- <name filter>]`

use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
use std::time::{Duration, Instant};

use waffle_con::cdwfa_config::CdwfaConfigBuilder;
use waffle_con::consensus::ConsensusDWFA;
use waffle_con::dual_consensus::DualConsensusDWFA;
use waffle_con::example_gen::{generate_dual_test, generate_test};

const BUDGET_SECS: f64 = 0.3;
const MIN_REPS: usize = 3;

struct Case {
    label: String,
    dual: bool,
    min_count: u64,
    dataset: Vec<Vec<u8>>,
}

fn build_cases() -> Vec<Case> {
    let alphabet_size = 4;
    let mut cases = vec![];
    for dual in [false, true] {
        for sl in [50, 200, 1000] {
            let num_differences = (sl / 50).max(1);
            for ns in [10, 30] {
                for er in [0.0, 0.01] {
                    let (label, dataset) = if dual {
                        let (_a1, _a2, d) = generate_dual_test(alphabet_size, sl, ns, er, num_differences);
                        (format!("dual_{alphabet_size}x{sl}x{ns}_{er}_diff{num_differences}"), d)
                    } else {
                        let (_c, d) = generate_test(alphabet_size, sl, ns, er);
                        (format!("consensus_{alphabet_size}x{sl}x{ns}_{er}"), d)
                    };
                    cases.push(Case { label, dual, min_count: (ns as u64) / 4, dataset });
                }
            }
        }
    }
    cases
}

/// Runs one case and returns the `Debug` text of its result.
fn run_case(case: &Case) -> String {
    let config = CdwfaConfigBuilder::default().min_count(case.min_count).build().unwrap();
    if case.dual {
        let mut dwfa = DualConsensusDWFA::with_config(config).unwrap();
        for s in case.dataset.iter() {
            dwfa.add_sequence(s).unwrap();
        }
        format!("{:?}", dwfa.consensus().unwrap())
    } else {
        let mut dwfa = ConsensusDWFA::with_config(config).unwrap();
        for s in case.dataset.iter() {
            dwfa.add_sequence(s).unwrap();
        }
        format!("{:?}", dwfa.consensus().unwrap())
    }
}

fn digest(text: &str) -> u64 {
    let mut hasher = DefaultHasher::new();
    text.hash(&mut hasher);
    hasher.finish()
}

fn main() {
    let filter = std::env::args().nth(1).unwrap_or_default();
    let cases: Vec<Case> = build_cases().into_iter().filter(|c| c.label.contains(&filter)).collect();

    let mut total_hasher = DefaultHasher::new();
    let mut total_min = Duration::ZERO;
    for case in cases.iter() {
        let output = run_case(case);
        let case_digest = digest(&output);
        case_digest.hash(&mut total_hasher);

        let mut best = Duration::MAX;
        let mut reps = 0;
        let start = Instant::now();
        while reps < MIN_REPS || start.elapsed().as_secs_f64() < BUDGET_SECS {
            let t = Instant::now();
            std::hint::black_box(run_case(case));
            best = best.min(t.elapsed());
            reps += 1;
        }
        total_min += best;
        println!("{:<34} {:>12.3} ms  reps={:<5} digest={:016x}", case.label, best.as_secs_f64() * 1e3, reps, case_digest);
    }
    println!("sum of minimums: {:.3} ms", total_min.as_secs_f64() * 1e3);
    println!("total digest: {:016x}", total_hasher.finish());
}
