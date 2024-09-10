/*!
This module provides access to the PriorityConsensusDWFA, which iteratively uses a DualConsensusDWFA on multiple provided sequences.
In the classic paired example, it will first generate DualConsensus on the first sequences in the pair and then further split on the second sequence (if possible).
Additionally, the reported consensus would be a pair of consensus sequences.

# Example usage
```rust
use waffle_con::consensus::Consensus;
use waffle_con::priority_consensus::PriorityConsensusDWFA;
use waffle_con::cdwfa_config::ConsensusCost;

let sequences_chains = [
    vec![b"TCCGT".to_vec(), b"TCCGT".to_vec()], // con 1
    vec![b"TCCGT".to_vec(), b"TCCGT".to_vec()], // con 1
    vec![b"TCCGT".to_vec(), b"TCCGT".to_vec()], // con 1
    vec![b"TCCGT".to_vec(), b"ACGGT".to_vec()], // con 2
    vec![b"TCCGT".to_vec(), b"ACGGT".to_vec()], // con 2
    vec![b"TCCGT".to_vec(), b"ACGGT".to_vec()], // con 2
    vec![b"ACGT".to_vec(), b"ACCCGGTT".to_vec()], // con 3
    vec![b"ACGT".to_vec(), b"ACCCGGTT".to_vec()], // con 3
    vec![b"ACGT".to_vec(), b"ACCCGGTT".to_vec()], // con 3
];

// add all the sequences
let mut cdwfa: PriorityConsensusDWFA = Default::default();
for s_chain in sequences_chains.iter() {
    let ref_chain: Vec<&[u8]> = s_chain.iter().map(|c| c.as_slice()).collect();
    cdwfa.add_sequence_chain(ref_chain).unwrap();
}

// run consensus and check results
let consensuses = cdwfa.consensus().unwrap();
assert_eq!(consensuses.consensuses(), &[
    // these are in alphabetically ordered by the chains; in the example these chains are pairs (i.e. length 2)
    vec![
        Consensus::new(b"ACGT".to_vec(), ConsensusCost::L1Distance, vec![0; 3]),
        Consensus::new(b"ACCCGGTT".to_vec(), ConsensusCost::L1Distance, vec![0; 3])
    ],
    vec![
        Consensus::new(b"TCCGT".to_vec(), ConsensusCost::L1Distance, vec![0; 6]), // this is shared between consensus 1 and 2, so it has costs for both
        Consensus::new(b"ACGGT".to_vec(), ConsensusCost::L1Distance, vec![0; 3])
    ],
    vec![
        Consensus::new(b"TCCGT".to_vec(), ConsensusCost::L1Distance, vec![0; 6]), // this is shared between consensus 1 and 2, so it has costs for both
        Consensus::new(b"TCCGT".to_vec(), ConsensusCost::L1Distance, vec![0; 3])
    ]
]);
assert_eq!(consensuses.sequence_indices(), &[
    2, 2, 2, 1, 1, 1, 0, 0, 0
]);
```
*/

use itertools::Itertools;
use log::debug;
use rustc_hash::FxHashSet as HashSet;
use simple_error::bail;

use crate::cdwfa_config::{CdwfaConfig, ConsensusCost};
use crate::consensus::Consensus;
use crate::dual_consensus::DualConsensusDWFA;

/// Contains a final multi-consensus result
#[derive(Debug, PartialEq)]
pub struct PriorityConsensus {
    /// Results for the consensuses
    consensuses: Vec<Vec<Consensus>>,
    /// For each input sequence set, this is the index of the assigned consensus
    sequence_indices: Vec<usize>
}

impl PriorityConsensus {
    pub fn new(consensuses: Vec<Vec<Consensus>>, sequence_indices: Vec<usize>) -> PriorityConsensus {
        PriorityConsensus {
            consensuses,
            sequence_indices
        }
    }

    // Getters
    pub fn consensuses(&self) -> &[Vec<Consensus>] {
        &self.consensuses
    }

    pub fn sequence_indices(&self) -> &[usize] {
        &self.sequence_indices
    }
}

/// Core utility that will generate a consensus sequence using the DWFA approach.
/// This approach is a wrapper for MultiConsensusDWFA that allows for a consensus to get split based on multiple sequences in a priority order.
/// For example, this may be used to first split sequences by homo-polymer compressed sequences and then the full-length sequences.
/// All consensuses are returned as a "chain" of consensuses that a user can then make decisions with.
#[derive(Debug, Default)]
pub struct PriorityConsensusDWFA<'a>  {
    /// Contains all the sequences that have been added to this consensus so far.
    sequences: Vec<Vec<&'a [u8]>>,
    /// Approximate offsets into the consensus that sequence starts, this will get adjust at run-time. If None, then assume start.
    offsets: Vec<Vec<Option<usize>>>,
    /// Optional start seeds for running consensus
    seed_groups: Vec<Option<usize>>,
    /// The config for this consensus run
    config: CdwfaConfig,
    /// The alphabet we will use for consensus building
    alphabet: HashSet<u8>
}

impl<'a> PriorityConsensusDWFA<'a> {
    /// Creates a new instance of ConsensusDWFA and performs sanity checks.
    /// # Arguments
    /// * `config` - the type of consensus score we want to use
    /// # Errors
    /// * None so far
    pub fn with_config(config: CdwfaConfig) -> Result<PriorityConsensusDWFA<'a>, Box<dyn std::error::Error>> {
        Ok(PriorityConsensusDWFA {
            sequences: vec![],
            offsets: vec![],
            seed_groups: vec![],
            config,
            alphabet: Default::default()
        })
    }

    /// Wrapper function for adding an unseeded sequence chain to the list.
    /// # Arguments
    /// * `sequences` - the new sequences to add
    /// # Errors
    /// * Same as `add_seeded_sequence(...)`
    pub fn add_sequence_chain(&mut self, sequences: Vec<&'a [u8]>) -> Result<(), Box<dyn std::error::Error>> {
        let offsets = vec![None; sequences.len()];
        self.add_seeded_sequence_chain(sequences, offsets, None)
    }

    /// Adds a new sequence chain to the list with an optional seed group
    /// # Arguments
    /// * `sequences` - the new sequences to add
    /// * `offset` - the approximate offset of the sequence, use None if at the start
    /// * `seed_group` - an optional starting seed grouping for this sequence
    /// # Errors
    /// * None so far
    pub fn add_seeded_sequence_chain(&mut self, sequences: Vec<&'a [u8]>, offsets: Vec<Option<usize>>, seed_group: Option<usize>) -> Result<(), Box<dyn std::error::Error>> {
        // check for the errors that can break assumptions later
        if sequences.is_empty() {
            bail!("Must provide a non-empty sequences Vec");
        }
        
        if !self.sequences.is_empty() && self.sequences[0].len() != sequences.len() {
            bail!("Expected sequences Vec of length {}, but got one of length {}", self.sequences[0].len(), sequences.len());
        }

        // add ever character to our alphabet
        for sequence in sequences.iter() {
            self.alphabet.extend(sequence.iter().cloned());
        }
        
        // make sure we didn't add the wildcard by accident
        // note: it's easier to just remove it once here than do something fancy with the above iteration
        if let Some(wc) = self.config.wildcard {
            self.alphabet.remove(&wc);
        }

        // now add the sequences
        self.sequences.push(sequences);
        self.offsets.push(offsets);
        self.seed_groups.push(seed_group);
        Ok(())
    }

    /// The core function that gets called after adding all the sequence sets we care about
    /// # Errors
    /// * if the DWFA throws errors
    pub fn consensus(&self) -> Result<PriorityConsensus, Box<dyn std::error::Error>> {
        // initialize the first set with all the sequences
        let max_split_level = self.sequences[0].len();
        let mut to_split: Vec<Vec<bool>> = vec![];
        let mut split_levels: Vec<usize> = vec![];
        let mut consensus_chains: Vec<Vec<Consensus>> = vec![];

        // figure out the list of initial seeds, typicall all are "None"
        let initial_group_keys: HashSet<Option<usize>> = self.seed_groups.iter().cloned().collect();
        for igk in initial_group_keys.iter() {
            // build a vec of which ones have this particular seed and add it to the list
            let initial_group_set: Vec<bool> = self.seed_groups.iter()
                .map(|sg| sg == igk)
                .collect();

            // pushing this grouping, and indicate we are looking at the first set of sequences (0)
            to_split.push(initial_group_set);
            split_levels.push(0);
            consensus_chains.push(vec![]);
        }

        let mut consensuses = vec![];
        let mut assignments = vec![];
        while let Some(include_set) = to_split.pop() {
            // also pop the current split level
            let current_split_level = split_levels.pop().unwrap();
            let mut current_consensus_chain = consensus_chains.pop().unwrap();

            // build a dual consensus DWFA for these sequences
            let mut dc_dwfa = DualConsensusDWFA::with_config(self.config.clone())?;
            debug!("Calling Dual at level {current_split_level} with: {include_set:?}");

            for (&include, (seq, offset)) in include_set.iter().zip(self.sequences.iter().zip(self.offsets.iter())) {
                if include {
                    dc_dwfa.add_sequence_offset(seq[current_split_level], offset[current_split_level])?;
                }
            }

            // now solve it
            let dc_result = dc_dwfa.consensus()?;
            if dc_result.len() > 1 {
                debug!("Multiple dual consensuses detected, arbitrarily selecting first option.");
            }
            
            let chosen_result = &dc_result[0];
            debug!("Parsing result with {} consensuses...", if chosen_result.is_dual() { 2 } else { 1 });

            /*
            // this will print the consensus sequences
            println!(">c1\n{}", std::str::from_utf8(chosen_result.consensus1().sequence()).unwrap());
            if chosen_result.is_dual() {
                println!(">c2\n{}", std::str::from_utf8(chosen_result.consensus2().unwrap().sequence()).unwrap());
            }

            let index1: usize = 257; // this one is going into UNKNOWN eventually, why are they splitting?
            let index2 = 301; // this one is resolving to D7
            // /*
            if include_set[index1] && include_set[index2] {
                debug!("BOTH INSIDE {index1} & {index2}");
            }
            // */

            // this was a debugging where we would score the full WFA ed of each sequence
            let mut ic_index = 0;
            for (i, seq_chain) in self.sequences.iter().enumerate() {
                if include_set[i] {
                    let require_both_end = false; // set to true for end-to-end ED
                    let s1 = crate::sequence_alignment::wfa_ed_config(chosen_result.consensus1().sequence(), seq_chain[current_split_level], require_both_end, self.config.wildcard);
                    let s2 = if chosen_result.is_dual() {
                        crate::sequence_alignment::wfa_ed_config(chosen_result.consensus2().unwrap().sequence(), seq_chain[current_split_level], require_both_end, self.config.wildcard)
                    } else {
                        usize::MAX
                    };
                    debug!("\tseq_{i} => {s1} | {s2} => is_ci = {}; {:?} | {:?}", chosen_result.is_consensus1()[ic_index], chosen_result.scores1()[ic_index], chosen_result.scores2()[ic_index]);
                    ic_index += 1;
                    /*
                    if i == index1 || i == index2 {
                        println!(">seq_{i}\n{}", std::str::from_utf8(seq_chain[current_split_level]).unwrap());
                    }
                    // */
                }
            }
            // */

            if chosen_result.is_dual() {
                // consensus sequences actually don't matter at this point, we only care about how they were split
                let is_c1 = chosen_result.is_consensus1();
                let mut is_c1_index: usize = 0;

                // these are the new clusters
                let mut assign1 = vec![false; self.sequences.len()];
                let mut assign2 = vec![false; self.sequences.len()];

                for (i, &included) in include_set.iter().enumerate() {
                    if included {
                        // this one was part of the include_set
                        if is_c1[is_c1_index] {
                            // this one should go to the first one
                            assign1[i] = true;
                        } else {
                            assign2[i] = true;
                        }
                        is_c1_index += 1;
                    }
                }

                // make sure we have written the correct number of things
                assert_eq!(is_c1.len(), is_c1_index);

                // now add both new sets for re-splitting
                to_split.push(assign1);
                split_levels.push(current_split_level); // level has not increased since we found a split
                consensus_chains.push(current_consensus_chain.clone()); // consensus chain thus far is shared

                to_split.push(assign2);
                split_levels.push(current_split_level); // level has not increased since we found a split
                consensus_chains.push(current_consensus_chain); // consensus chain thus far is shared
            } else {
                let new_split_level = current_split_level + 1;
                current_consensus_chain.push(chosen_result.consensus1().clone());

                if new_split_level == max_split_level {
                    // we found no split at the final level, time to return the consensus
                    // consensuses.push(chosen_result.consensus1().clone());
                    consensuses.push(current_consensus_chain);
                    assignments.push(include_set);
                } else {
                    // we want to try to split at the next split level
                    to_split.push(include_set);
                    split_levels.push(new_split_level);

                    // we found a consensus for this level, save it to the chain and add it back in to the looping
                    consensus_chains.push(current_consensus_chain);
                }
            }
        }

        let (consensuses, sequence_indices) = if consensuses.len() > 1 {
            // probably need to sort and correct ordering, etc.
            let mut indices = vec![usize::MAX; self.sequences.len()];
            let sorted_cons = consensuses.into_iter().zip(assignments)
                .sorted_by(|a, b| {
                    let chain_a = a.0.iter().map(|c| c.sequence()).collect::<Vec<&[u8]>>();
                    let chain_b = b.0.iter().map(|c| c.sequence()).collect::<Vec<&[u8]>>();
                    chain_a.cmp(&chain_b)
                })
                .enumerate()
                .map(|(con_index, (consensus, assignments))| {
                    for (i, &assigned) in assignments.iter().enumerate() {
                        if assigned {
                            // make sure we have not assigned this one
                            assert_eq!(indices[i], usize::MAX);
                            indices[i] = con_index;
                        }
                    }
                    consensus
                })
                .collect();
            
            (sorted_cons, indices)
        } else {
            // single value, pass through
            (consensuses, vec![0; self.sequences.len()])
        };

        Ok(PriorityConsensus {
            consensuses,
            sequence_indices
        })
    }

    // getters
    pub fn sequences(&self) -> &[Vec<&'a [u8]>] {
        &self.sequences
    }

    pub fn consensus_cost(&self) -> ConsensusCost {
        self.config.consensus_cost
    }

    pub fn alphabet(&self) -> &HashSet<u8> {
        &self.alphabet
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    use std::path::PathBuf;

    use crate::cdwfa_config::CdwfaConfigBuilder;

    #[derive(Debug, serde::Deserialize)]
    struct MultiRecord {
        consensus: usize,
        edits: usize,
        sequence: String
    }

    /// Wrapper test function that loads a multi test from a csv file.
    /// Expected columns are "consensus" (integer >= 1), "edits" (u64), and "sequence" (String).
    /// The "sequence" column is a ";" separated series of sequences.
    /// The "edits" column is largely ignored except for marking the true consensus. The exact ED testing is in other tooling.
    /// Returns a tuple of (sequence_chains, PriorityConsensus).
    /// # Arguments
    /// * `filename` - the file path to load
    /// * `include_consensus` - if True, it will load the first consensus read into the sequences
    /// * `cost_mode` - the cost mode getting tested
    fn load_priority_csv_test(filename: &std::path::Path, include_consensus: bool, cost_mode: ConsensusCost) -> (Vec<Vec<Vec<u8>>>, PriorityConsensus) {
        // contains the consensus sequences
        let mut consensuses: Vec<Vec<Vec<u8>>> = vec![];
        // contains the edit distances for mappings the sequences to the consensus
        // let mut edit_distances: Vec<Vec<usize>> = vec![];
        
        // contains the sequences to get into the consensus algorithm
        let mut sequence_chains: Vec<Vec<Vec<u8>>> = vec![];
        // contains the index of the consensus they _should_ map to
        let mut sequence_indices: Vec<usize> = vec![];
        
        let mut csv_reader = csv::ReaderBuilder::new()
            .has_headers(true)
            .from_path(filename)
            .unwrap();
        for row in csv_reader.deserialize() {
            let record: MultiRecord = row.unwrap();
            assert!(record.consensus >= 1, "consensus column must be >= 1");
            let con_index = record.consensus - 1;
            let edits = match cost_mode {
                ConsensusCost::L1Distance => record.edits,
                ConsensusCost::L2Distance => record.edits.pow(2)
            };
            let sequence_chain: Vec<Vec<u8>> = record.sequence
                .split(';')
                .map(|s| s.as_bytes().to_vec())
                .collect();

            if con_index >= consensuses.len() {
                // an index is beyond our current length, extend
                consensuses.resize(con_index+1, vec![]);
                // edit_distances.resize(con_index+1, vec![]);
            }

            // check if this is the consensus
            if edits == 0 && consensuses[con_index].is_empty() {
                // this is the consensus
                consensuses[con_index] = sequence_chain.clone();
                // everything after here is for normal reads, so early exit if we don't include the read
                if !include_consensus {
                    continue;
                }
            }

            // edit_distances[con_index].push(edits);
            
            sequence_chains.push(sequence_chain);
            sequence_indices.push(con_index);
        }

        // make sure all consensuses are non-empty
        assert!(consensuses.iter().all(|v| !v.is_empty()));
        // make sure all consensuses have at least 1 read, how else would we ever find it?
        assert!(sequence_chains.iter().all(|edv| !edv.is_empty()));

        // figure out the consensus order which is just the sorted order
        for c in consensuses.iter() {
            println!("{c:?}");
        }

        for i in 0..consensuses.len() {
            for j in (i+1)..consensuses.len() {
                println!("{}.cmp({}) = {:?}", i, j, consensuses[i].cmp(&consensuses[j]));
            }
        }

        // future: I feel like this is over-complicated
        // first figure out the correct index order
        let argsort: Vec<usize> = (0..consensuses.len())
            .sorted_by_key(|&i| &consensuses[i])
            .collect();
        // now convert that into a lookup table we can use for ordering them
        let arg_lookup: Vec<usize> = (0..consensuses.len())
            .sorted_by_key(|&i| argsort[i])
            .collect();

        // reorder everything based on the argsort lookup
        let consensuses: Vec<Vec<Vec<u8>>> = consensuses.into_iter().enumerate()
            .sorted_by_key(|(i, _v)| arg_lookup[*i])
            .map(|(_i, v)| v)
            .collect();
        /*
        let edit_distances: Vec<Vec<usize>> = edit_distances.into_iter().enumerate()
            .sorted_by_key(|(i, _v)| arg_lookup[*i])
            .map(|(_i, v)| v)
            .collect();
        */
        let sequence_indices: Vec<usize> = sequence_indices.into_iter()
            .map(|i| arg_lookup[i])
            .collect();

        // finally make the multi-consensus result
        let con_vec = consensuses.into_iter()
            .map(|con_chain| {
                con_chain.into_iter().map(|con| {
                    Consensus::new(con, cost_mode, vec![])
                })
                .collect()
            })
            .collect();

        let consensus = PriorityConsensus {
            consensuses: con_vec,
            sequence_indices
        };

        (sequence_chains, consensus)
    }

    /// Entry point for most file-based tests.
    /// # Arguments
    /// * `filename` - the test file to load, will be a csv
    /// * `include_consensus` - if True, this will include the consensus sequence line as a read
    /// * `opt_config` - optional alternate config, otherwise a default is created
    fn run_test_file(filename: &str, include_consensus: bool, opt_config: Option<CdwfaConfig>) {
        let config = opt_config.unwrap_or(
            CdwfaConfigBuilder::default()
                .wildcard(Some(b'*'))
                .build().unwrap()
        );

        let (sequences, expected_consensus) = load_priority_csv_test(&PathBuf::from(filename), include_consensus, config.consensus_cost);

        // set queue size large since we do not want to worry about filtering for this test
        let mut consensus_dwfa = PriorityConsensusDWFA::with_config(config).unwrap();
        for sequence_chain in sequences.iter() {
            let vec_ref: Vec<&[u8]> = sequence_chain.iter()
                .map(|s| s.as_slice())
                .collect();
            consensus_dwfa.add_sequence_chain(vec_ref).unwrap();
        }

        // make sure our alphabet has four symbols, ACGT
        assert_eq!(consensus_dwfa.alphabet().len(), 4);

        // now check that the consensus is the same as our sequence
        let consensus = consensus_dwfa.consensus().unwrap();
        // assert_eq!(consensus, expected_consensus); // we don't load edit distances, so we can't directly compare this
        assert_eq!(consensus.sequence_indices(), expected_consensus.sequence_indices()); // check consensus assignment
        assert_eq!(consensus.consensuses().len(), expected_consensus.consensuses().len()); // check consensus length
        for (cc, ecc) in consensus.consensuses().iter().zip(expected_consensus.consensuses().iter()) {
            assert_eq!(cc.len(), ecc.len()); // check chain length
            for (c, ec) in cc.iter().zip(ecc.iter()) {
                assert_eq!(c.sequence(), ec.sequence()); // check equality of the sequences in the chain
            }
        }
    }

    // these tests are just copied from the others to verify nothing breaks with single chains
    #[test]
    fn test_csv_dual_001() {
        run_test_file("./tests/dual_001.csv", true, None);
    }

    #[test]
    fn test_multi_exact_001() {
        run_test_file("./tests/multi_exact_001.csv", true, None);
    }

    #[test]
    fn test_multi_exact_002() {
        run_test_file("./tests/multi_exact_002.csv", true, None);
    }

    #[test]
    fn test_multi_err_001() {
        run_test_file("./tests/multi_err_001.csv", false, None);
    }

    #[test]
    fn test_multi_err_002() {
        run_test_file("./tests/multi_err_002.csv", false, None);
    }


    #[test]
    fn test_multi_samesplit_001() {
        // this one tests four sequences that all have a unique symbol [ACGT] at the same position, so we must split into all four at once
        run_test_file("./tests/multi_samesplit_001.csv", true, None);
    }

    #[test]
    fn test_multi_postcon_001() {
        // this one tests a scenario where the consensus will get split off correctly, but we need to re-run to get the best overall consensus sequence
        run_test_file("./tests/multi_postcon_001.csv", true, Some(
            CdwfaConfigBuilder::default()
                .wildcard(Some(b'*'))
                .min_count(2) // has to be lower here due to the errors we inject
                .build().unwrap()
        ));
    }

    // now begins the multi-tests
    #[test]
    fn test_single_sequence() { // simple double chain with a single entry
        let sequence = b"ACGTACGTACGT";
        let mut consensus_dwfa = PriorityConsensusDWFA::default();
        consensus_dwfa.add_sequence_chain(vec![sequence, sequence]).unwrap();

        // make sure our alphabet has four symbols, ACGT
        assert_eq!(consensus_dwfa.alphabet().len(), 4);

        // now check that the consensus is the same as our sequence
        let consensus = consensus_dwfa.consensus().unwrap();
        assert_eq!(consensus, PriorityConsensus {
            consensuses: vec![
                vec![Consensus::new(
                    sequence.to_vec(),
                    ConsensusCost::L1Distance,
                    vec![0]
                ); 2]
            ],
            sequence_indices: vec![0]
        });
    }

    #[test]
    fn doc_test_example() {
        let sequences_chains = [
            vec![b"TCCGT".to_vec(), b"TCCGT".to_vec()], // con 1
            vec![b"TCCGT".to_vec(), b"TCCGT".to_vec()], // con 1
            vec![b"TCCGT".to_vec(), b"TCCGT".to_vec()], // con 1
            vec![b"TCCGT".to_vec(), b"ACGGT".to_vec()], // con 2
            vec![b"TCCGT".to_vec(), b"ACGGT".to_vec()], // con 2
            vec![b"TCCGT".to_vec(), b"ACGGT".to_vec()], // con 2
            vec![b"ACGT".to_vec(), b"ACCCGGTT".to_vec()], // con 3
            vec![b"ACGT".to_vec(), b"ACCCGGTT".to_vec()], // con 3
            vec![b"ACGT".to_vec(), b"ACCCGGTT".to_vec()], // con 3
        ];

        // add all the sequences
        let mut cdwfa: PriorityConsensusDWFA = Default::default();
        for s_chain in sequences_chains.iter() {
            let ref_chain: Vec<&[u8]> = s_chain.iter().map(|c| c.as_slice()).collect();
            cdwfa.add_sequence_chain(ref_chain).unwrap();
        }

        // run consensus and check results
        let consensuses = cdwfa.consensus().unwrap();
        assert_eq!(consensuses.consensuses(), &[
            // these are in alphabetically order by the consensus content
            vec![
                Consensus::new(b"ACGT".to_vec(), ConsensusCost::L1Distance, vec![0; 3]),
                Consensus::new(b"ACCCGGTT".to_vec(), ConsensusCost::L1Distance, vec![0; 3])
            ],
            vec![
                Consensus::new(b"TCCGT".to_vec(), ConsensusCost::L1Distance, vec![0; 6]),
                Consensus::new(b"ACGGT".to_vec(), ConsensusCost::L1Distance, vec![0; 3])
            ],
            vec![
                Consensus::new(b"TCCGT".to_vec(), ConsensusCost::L1Distance, vec![0; 6]),
                Consensus::new(b"TCCGT".to_vec(), ConsensusCost::L1Distance, vec![0; 3])
            ]
        ]);
        assert_eq!(consensuses.sequence_indices(), &[
            2, 2, 2, 1, 1, 1, 0, 0, 0
        ]);
    }

    #[test]
    fn test_priority_001() {
        run_test_file("./tests/priority_001.csv", true, None);
    }

    #[test]
    fn test_priority_002() {
        run_test_file("./tests/priority_002.csv", true, None);
    }
    
    #[test]
    fn test_priority_003() {
        run_test_file("./tests/priority_003.csv", true, None);
    }
}