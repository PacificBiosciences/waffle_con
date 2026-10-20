mod common;

use waffle_con::cdwfa_config::{CdwfaConfig, ConsensusCost};
use waffle_con::consensus::{Consensus, ConsensusAssignment, ConsensusDWFA};
use waffle_con::dual_consensus::{DualConsensus, DualConsensusDWFA, SequenceAssignment};

fn scores(reads: &[&[u8]], haplotype: &[u8], cost: ConsensusCost) -> Vec<usize> {
    reads
        .iter()
        .map(|read| {
            let distance = common::edit_distance(read, haplotype, None, false);
            match cost {
                ConsensusCost::L1Distance => distance,
                ConsensusCost::L2Distance => distance * distance,
            }
        })
        .collect()
}

fn expected_consensus(
    sequence: &[u8],
    cost: ConsensusCost,
    scores: Vec<usize>,
    assignments: Option<Vec<ConsensusAssignment>>,
) -> Consensus {
    Consensus::new(sequence.to_vec(), cost, scores, assignments).unwrap()
}

fn single_consensus(reads: &[&[u8]], config: CdwfaConfig) -> Vec<Consensus> {
    let mut engine = ConsensusDWFA::with_config(config).unwrap();
    for &read in reads {
        engine.add_sequence(read).unwrap();
    }
    engine.consensus().unwrap()
}

fn dual_consensus(reads: &[&[u8]], config: CdwfaConfig) -> Vec<DualConsensus> {
    let mut engine = DualConsensusDWFA::with_config(config).unwrap();
    for &read in reads {
        engine.add_sequence(read).unwrap();
    }
    engine.consensus().unwrap()
}

fn assert_dual_results(actual: &[DualConsensus], expected: &[DualConsensus]) {
    // DualConsensus::PartialEq does not compare scores1 or scores2.
    assert_eq!(actual, expected);
    for (actual, expected) in actual.iter().zip(expected) {
        assert_eq!(actual.is_dual(), expected.is_dual());
        assert_eq!(actual.scores1(), expected.scores1());
        assert_eq!(actual.scores2(), expected.scores2());
    }
}

#[test]
fn single_consensus_reports_finalized_global_l1_and_l2_scores() {
    let reads: [&[u8]; 6] = [b"AC", b"AC", b"AC", b"AC", b"AC", b"CA"];
    for cost in [ConsensusCost::L1Distance, ConsensusCost::L2Distance] {
        // AC -> CA takes two edits, although its live prefix distance is one.
        // Five exact reads make AC strictly better than a one-edit compromise
        // under both objectives (2 < 6 for L1; 4 < 6 for L2).
        let expected_scores = scores(&reads, b"AC", cost);
        let outlier_score = if cost == ConsensusCost::L1Distance {
            2
        } else {
            4
        };
        assert_eq!(expected_scores, [0, 0, 0, 0, 0, outlier_score]);
        let actual = single_consensus(
            &reads,
            CdwfaConfig {
                consensus_cost: cost,
                ..Default::default()
            },
        );
        assert_eq!(
            actual,
            [expected_consensus(b"AC", cost, expected_scores, None)]
        );
    }
}

#[test]
fn dual_search_without_a_supported_split_reports_finalized_global_scores() {
    let reads: [&[u8]; 6] = [b"AC", b"AC", b"AC", b"AC", b"AC", b"CA"];
    for cost in [ConsensusCost::L1Distance, ConsensusCost::L2Distance] {
        for weighted_by_ed in [false, true] {
            let expected_scores = scores(&reads, b"AC", cost);
            let outlier_score = if cost == ConsensusCost::L1Distance {
                2
            } else {
                4
            };
            assert_eq!(expected_scores, [0, 0, 0, 0, 0, outlier_score]);
            let expected = DualConsensus::new(
                expected_consensus(b"AC", cost, expected_scores.clone(), None),
                None,
                vec![SequenceAssignment::Consensus1; reads.len()],
                expected_scores.into_iter().map(Some).collect(),
                vec![None; reads.len()],
            )
            .unwrap();
            let actual = dual_consensus(
                &reads,
                CdwfaConfig {
                    consensus_cost: cost,
                    weighted_by_ed,
                    ..Default::default()
                },
            );
            assert_dual_results(&actual, &[expected]);
        }
    }
}

#[test]
fn single_search_keeps_sibling_alleles_independent_and_honors_first_seen_return_limit() {
    // The first nominated allele is deliberately not lexicographically first.
    let reads: [&[u8]; 2] = [b"ATGT", b"ACGT"];
    for cost in [ConsensusCost::L1Distance, ConsensusCost::L2Distance] {
        let config = CdwfaConfig {
            consensus_cost: cost,
            min_count: 1,
            max_return_size: 1,
            ..Default::default()
        };
        let first_scores = scores(&reads, b"ATGT", cost);
        let second_scores = scores(&reads, b"ACGT", cost);
        assert_eq!(first_scores, [0, 1]);
        assert_eq!(second_scores, [1, 0]);
        assert_eq!(
            single_consensus(&reads, config.clone()),
            [expected_consensus(
                b"ATGT",
                cost,
                first_scores.clone(),
                None
            )],
        );

        // With room for both answers, the API returns sequence-sorted results.
        // Extending/finalizing the first branch must not overwrite its sibling.
        assert_eq!(
            single_consensus(
                &reads,
                CdwfaConfig {
                    max_return_size: 2,
                    ..config
                }
            ),
            [
                expected_consensus(b"ACGT", cost, second_scores, None),
                expected_consensus(b"ATGT", cost, first_scores, None),
            ],
        );
    }
}

#[test]
fn dual_common_prefix_split_preserves_alleles_and_equal_score_read() {
    let reads: [&[u8]; 7] = [
        b"ATGT", b"ATGT", b"ATGT", b"ACGT", b"ACGT", b"ACGT", b"AGGT",
    ];
    for cost in [ConsensusCost::L1Distance, ConsensusCost::L2Distance] {
        for weighted_by_ed in [false, true] {
            let scores1 = scores(&reads, b"ACGT", cost);
            let scores2 = scores(&reads, b"ATGT", cost);
            assert_eq!(scores1, [1, 1, 1, 0, 0, 0, 1]);
            assert_eq!(scores2, [0, 0, 0, 1, 1, 1, 1]);
            // C and T have three supporting reads each; G has only one.
            // AGGT is equally close to both alleles and belongs to neither
            // allele's assigned-read score list, but appears in both arrays.
            let expected = DualConsensus::new(
                expected_consensus(b"ACGT", cost, vec![0; 3], None),
                Some(expected_consensus(b"ATGT", cost, vec![0; 3], None)),
                vec![
                    SequenceAssignment::Consensus2,
                    SequenceAssignment::Consensus2,
                    SequenceAssignment::Consensus2,
                    SequenceAssignment::Consensus1,
                    SequenceAssignment::Consensus1,
                    SequenceAssignment::Consensus1,
                    SequenceAssignment::EqualScore,
                ],
                scores1.into_iter().map(Some).collect(),
                scores2.into_iter().map(Some).collect(),
            )
            .unwrap();
            let actual = dual_consensus(
                &reads,
                CdwfaConfig {
                    consensus_cost: cost,
                    weighted_by_ed,
                    ..Default::default()
                },
            );
            assert_dual_results(&actual, &[expected]);
        }
    }
}

#[test]
fn dual_minority_support_at_or_below_count_and_ceiling_fraction_thresholds() {
    let reads: [&[u8]; 7] = [b"A", b"A", b"A", b"A", b"A", b"C", b"C"];
    // At 1/4 of seven reads, ceil(7/4) = 2: the minority qualifies.
    // At 3/8, ceil(21/8) = 3: the same two reads do not qualify.
    // These exact binary fractions do not depend on fractional vote rounding.
    for (min_count, min_af, splits) in [
        (1, 0.25, true),
        (1, 0.375, false),
        (2, 0.0, true),
        (3, 0.0, false),
    ] {
        for cost in [ConsensusCost::L1Distance, ConsensusCost::L2Distance] {
            for weighted_by_ed in [false, true] {
                let scores1 = scores(&reads, b"A", cost);
                let scores2 = scores(&reads, b"C", cost);
                assert_eq!(scores1, [0, 0, 0, 0, 0, 1, 1]);
                assert_eq!(scores2, [1, 1, 1, 1, 1, 0, 0]);
                let expected = if splits {
                    DualConsensus::new(
                        expected_consensus(b"A", cost, vec![0; 5], None),
                        Some(expected_consensus(b"C", cost, vec![0; 2], None)),
                        vec![
                            SequenceAssignment::Consensus1,
                            SequenceAssignment::Consensus1,
                            SequenceAssignment::Consensus1,
                            SequenceAssignment::Consensus1,
                            SequenceAssignment::Consensus1,
                            SequenceAssignment::Consensus2,
                            SequenceAssignment::Consensus2,
                        ],
                        scores1.into_iter().map(Some).collect(),
                        scores2.into_iter().map(Some).collect(),
                    )
                    .unwrap()
                } else {
                    DualConsensus::new(
                        expected_consensus(b"A", cost, scores1.clone(), None),
                        None,
                        vec![SequenceAssignment::Consensus1; reads.len()],
                        scores1.into_iter().map(Some).collect(),
                        vec![None; reads.len()],
                    )
                    .unwrap()
                };
                let actual = dual_consensus(
                    &reads,
                    CdwfaConfig {
                        consensus_cost: cost,
                        min_count,
                        min_af,
                        weighted_by_ed,
                        ..Default::default()
                    },
                );
                assert_dual_results(&actual, &[expected]);
            }
        }
    }
}

#[test]
fn single_cap_distinguishes_a_tracked_error_from_a_halted_outlier() {
    let reads: [&[u8]; 5] = [
        b"TTTTTTTT",
        b"ACGTACGT",
        b"ACGTACGT",
        b"ACGTACGT",
        b"ACCTACGT",
    ];
    // The first read is six edits away; the last is only one edit away.
    assert_eq!(
        scores(&reads, b"ACGTACGT", ConsensusCost::L1Distance),
        [6, 0, 0, 0, 1]
    );
    for cost in [ConsensusCost::L1Distance, ConsensusCost::L2Distance] {
        let capped_score = if cost == ConsensusCost::L1Distance {
            2
        } else {
            4
        };
        let expected = expected_consensus(
            b"ACGTACGT",
            cost,
            vec![capped_score, 0, 0, 0, 1],
            Some(vec![
                ConsensusAssignment::EditDistanceLimit,
                ConsensusAssignment::Included,
                ConsensusAssignment::Included,
                ConsensusAssignment::Included,
                ConsensusAssignment::Included,
            ]),
        );
        let actual = single_consensus(
            &reads,
            CdwfaConfig {
                consensus_cost: cost,
                max_edit_distance: Some(2),
                ..Default::default()
            },
        );
        assert_eq!(actual, [expected]);
    }
}

#[test]
fn unweighted_dual_cap_retains_complete_pruned_capped_and_tied_score_arrays() {
    let reads: [&[u8]; 8] = [
        b"ATGT", b"ATGT", b"ATGT", b"GGGG", b"ACGT", b"ACGT", b"ACGT", b"AGGT",
    ];
    assert_eq!(
        scores(&reads, b"ACGT", ConsensusCost::L1Distance),
        [1, 1, 1, 3, 0, 0, 0, 1]
    );
    assert_eq!(
        scores(&reads, b"ATGT", ConsensusCost::L1Distance),
        [0, 0, 0, 3, 1, 1, 1, 1]
    );
    for cost in [ConsensusCost::L1Distance, ConsensusCost::L2Distance] {
        let capped_score = if cost == ConsensusCost::L1Distance {
            2
        } else {
            4
        };
        // Exact allele reads lose tracking on the opposite allele at the SNP.
        // GGGG remains tracked on both until it hits the cap; AGGT stays tied.
        let expected = DualConsensus::new(
            expected_consensus(b"ACGT", cost, vec![0; 3], None),
            Some(expected_consensus(b"ATGT", cost, vec![0; 3], None)),
            vec![
                SequenceAssignment::Consensus2,
                SequenceAssignment::Consensus2,
                SequenceAssignment::Consensus2,
                SequenceAssignment::EditDistanceLimit,
                SequenceAssignment::Consensus1,
                SequenceAssignment::Consensus1,
                SequenceAssignment::Consensus1,
                SequenceAssignment::EqualScore,
            ],
            vec![
                None,
                None,
                None,
                Some(capped_score),
                Some(0),
                Some(0),
                Some(0),
                Some(1),
            ],
            vec![
                Some(0),
                Some(0),
                Some(0),
                Some(capped_score),
                None,
                None,
                None,
                Some(1),
            ],
        )
        .unwrap();
        let actual = dual_consensus(
            &reads,
            CdwfaConfig {
                consensus_cost: cost,
                weighted_by_ed: false,
                max_edit_distance: Some(2),
                dual_max_ed_delta: 0,
                ..Default::default()
            },
        );
        assert_dual_results(&actual, &[expected]);
    }
}
