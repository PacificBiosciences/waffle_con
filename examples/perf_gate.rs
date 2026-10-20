//! Fast gate for performance work. It runs the same cases as `benches/consensus_bench.rs`.
//!
//! * Correctness: each case runs once and its result is hashed.
//!   A speed-up must not change the `total digest` line.
//! * Speed: each case repeats for about `BUDGET_SECS` (at least `MIN_REPS` times).
//!   The minimum wall time is reported, because it is the least noisy estimate.
//!   The timed runs only compute the result. They do not format it.
//!
//! Usage: `cargo run --release --example perf_gate [-- <name filter>]`

use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
use std::time::{Duration, Instant};

use waffle_con::cdwfa_config::CdwfaConfigBuilder;
use waffle_con::consensus::{Consensus, ConsensusDWFA};
use waffle_con::dual_consensus::{DualConsensus, DualConsensusDWFA};
use waffle_con::dynamic_wfa::{DWFALite, DWFALiteConfig, DWFALiteConfigBuilder};
use waffle_con::example_gen::{generate_dual_test, generate_finalize_pair, generate_test};

const BUDGET_SECS: f64 = 0.3;
const MIN_REPS: usize = 3;

enum Input {
    Consensus { min_count: u64, dataset: Vec<Vec<u8>> },
    Dual { min_count: u64, dataset: Vec<Vec<u8>> },
    Finalize { config: DWFALiteConfig, baseline: Vec<u8>, other: Vec<u8> },
}

enum Output {
    Consensus(Vec<Consensus>),
    Dual(Vec<DualConsensus>),
    Finalize(DWFALite),
}

struct Case {
    label: String,
    input: Input,
}

fn build_cases() -> Vec<Case> {
    let alphabet_size = 4;
    let mut cases = vec![];
    for dual in [false, true] {
        for sl in [50, 200, 1000] {
            let num_differences = (sl / 50).max(1);
            for ns in [10, 30] {
                for er in [0.0, 0.01] {
                    let min_count = (ns as u64) / 4;
                    cases.push(if dual {
                        let (_a1, _a2, dataset) = generate_dual_test(alphabet_size, sl, ns, er, num_differences);
                        Case {
                            label: format!("dual_{alphabet_size}x{sl}x{ns}_{er}_diff{num_differences}"),
                            input: Input::Dual { min_count, dataset },
                        }
                    } else {
                        let (_c, dataset) = generate_test(alphabet_size, sl, ns, er);
                        Case {
                            label: format!("consensus_{alphabet_size}x{sl}x{ns}_{er}"),
                            input: Input::Consensus { min_count, dataset },
                        }
                    });
                }
            }
        }
    }
    let (sl, er) = (4000, 0.01);
    for gap in [1000, 2000, 3000] {
        let config = DWFALiteConfigBuilder::default().max_edit_distance(Some(sl + gap)).build().unwrap();
        let (baseline, other) = generate_finalize_pair(alphabet_size, sl, er, gap);
        cases.push(Case {
            label: format!("finalize_{alphabet_size}x{sl}_{er}_gap{gap}"),
            input: Input::Finalize { config, baseline, other },
        });
    }
    cases
}

/// Runs one case and returns its result.
/// # Arguments
/// * `input` - the sequences and settings for one case
fn run_case(input: &Input) -> Output {
    match input {
        Input::Consensus { min_count, dataset } => {
            let config = CdwfaConfigBuilder::default().min_count(*min_count).build().unwrap();
            let mut dwfa = ConsensusDWFA::with_config(config).unwrap();
            for s in dataset.iter() {
                dwfa.add_sequence(s).unwrap();
            }
            Output::Consensus(dwfa.consensus().unwrap())
        },
        Input::Dual { min_count, dataset } => {
            let config = CdwfaConfigBuilder::default().min_count(*min_count).build().unwrap();
            let mut dwfa = DualConsensusDWFA::with_config(config).unwrap();
            for s in dataset.iter() {
                dwfa.add_sequence(s).unwrap();
            }
            Output::Dual(dwfa.consensus().unwrap())
        },
        Input::Finalize { config, baseline, other } => {
            let mut dwfa = DWFALite::new(config.clone());
            dwfa.update(baseline, other).unwrap();
            dwfa.finalize(baseline, other).unwrap();
            Output::Finalize(dwfa)
        },
    }
}

/// Hashes the text of a result. A `DWFALite` is hashed through its public getters only,
/// so that a change to its private fields does not change the digest.
/// # Arguments
/// * `output` - the result of one case
fn digest(output: &Output) -> u64 {
    let text = match output {
        Output::Consensus(result) => format!("{result:?}"),
        Output::Dual(result) => format!("{result:?}"),
        Output::Finalize(dwfa) => format!(
            "{:?}",
            (dwfa.edit_distance(), dwfa.wavefront(), dwfa.state(), dwfa.maximum_baseline_distance(), dwfa.maximum_other_distance())
        ),
    };
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
        let case_digest = digest(&run_case(&case.input));
        case_digest.hash(&mut total_hasher);

        let mut best = Duration::MAX;
        let mut reps = 0;
        let start = Instant::now();
        while reps < MIN_REPS || start.elapsed().as_secs_f64() < BUDGET_SECS {
            let t = Instant::now();
            std::hint::black_box(run_case(&case.input));
            best = best.min(t.elapsed());
            reps += 1;
        }
        total_min += best;
        println!("{:<34} {:>12.3} ms  reps={:<5} digest={:016x}", case.label, best.as_secs_f64() * 1e3, reps, case_digest);
    }
    println!("sum of minimums: {:.3} ms", total_min.as_secs_f64() * 1e3);
    println!("total digest: {:016x}", total_hasher.finish());
}
