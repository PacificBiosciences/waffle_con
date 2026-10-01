
use criterion::{black_box, criterion_group, criterion_main, Criterion};

use waffle_con::cdwfa_config::CdwfaConfigBuilder;
use waffle_con::consensus::ConsensusDWFA;
use waffle_con::dual_consensus::DualConsensusDWFA;
use waffle_con::example_gen::{generate_dual_test, generate_test};

pub fn bench_consensus(c: &mut Criterion) {
    let alphabet_size = 4;
    let seq_lens = [50, 200, 1000];
    let num_samples = [10, 30];
    let error_rates = [0.0, 0.01];

    let mut benchmark_group = c.benchmark_group("consensus-group");
    benchmark_group.sample_size(10);

    for &sl in seq_lens.iter() {
        for &ns in num_samples.iter() {
            // require 25% of reads to go forth
            let config = CdwfaConfigBuilder::default()
                .min_count((ns as u64) / 4)
                .build().unwrap();
            for &er in error_rates.iter() {
                let (_consensus, dataset) = generate_test(alphabet_size, sl, ns, er);
                /*
                // uncomment to print out the strings, mostly for initial testing
                println!("{consensus:?}");
                for d in dataset.iter() {
                    println!("{d:?}");
                }
                */
                let test_label = format!("consensus_{alphabet_size}x{sl}x{ns}_{er}");
                benchmark_group.bench_function(&test_label, |b| b.iter(|| {
                    black_box({
                        let mut consensus_dwfa = ConsensusDWFA::with_config(config.clone()).unwrap();
                        for s in dataset.iter() {
                            consensus_dwfa.add_sequence(s).unwrap();
                        }
                        let resolved_consensus = consensus_dwfa.consensus().unwrap();
                        
                        // this was an initial sanity check we did as a basic test
                        // assert_eq!(resolved_consensus[0].sequence(), &consensus);
                        
                        resolved_consensus
                    });
                }));
            }
        }
    }

    benchmark_group.finish();
}

pub fn bench_dual_consensus(c: &mut Criterion) {
    let alphabet_size = 4;
    let seq_lens = [50, 200, 1000];
    let num_samples = [10, 30];
    let error_rates = [0.0, 0.01];

    let mut benchmark_group = c.benchmark_group("dual-consensus-group");
    benchmark_group.sample_size(10);

    for &sl in seq_lens.iter() {
        // About 2% of positions differ between the two true alleles, and at least one.
        let num_differences = (sl / 50).max(1);
        for &ns in num_samples.iter() {
            // Require 25% of one allele's reads. Each allele contributes `ns` reads.
            let config = CdwfaConfigBuilder::default()
                .min_count((ns as u64) / 4)
                .build().unwrap();
            for &er in error_rates.iter() {
                let (_allele1, _allele2, dataset) = generate_dual_test(alphabet_size, sl, ns, er, num_differences);
                let test_label = format!("dual_{alphabet_size}x{sl}x{ns}_{er}_diff{num_differences}");
                benchmark_group.bench_function(&test_label, |b| b.iter(|| {
                    black_box({
                        let mut consensus_dwfa = DualConsensusDWFA::with_config(config.clone()).unwrap();
                        for s in dataset.iter() {
                            consensus_dwfa.add_sequence(s).unwrap();
                        }
                        consensus_dwfa.consensus().unwrap()
                    });
                }));
            }
        }
    }

    benchmark_group.finish();
}

criterion_group!(benches, bench_consensus, bench_dual_consensus);
criterion_main!(benches);