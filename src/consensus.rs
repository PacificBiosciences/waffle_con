
/*!
This module provides access to the ConsensusDWFA, which generates the single best consensus for a set of sequences.

# Example usage
```rust
use waffle_con::consensus::ConsensusDWFA;

let sequences = [
    b"ACGT".to_vec(),
    b"ACCGT".to_vec(), // this should be the consensus
    b"ACCCGT".to_vec()
];

// add all the sequences
let mut cdwfa: ConsensusDWFA = Default::default();
for s in sequences.iter() {
    cdwfa.add_sequence(s).unwrap();
}

// run consensus and check the results
let consensuses = cdwfa.consensus().unwrap();
assert_eq!(consensuses.len(), 1);
assert_eq!(consensuses[0].sequence(), sequences[1]);
assert_eq!(consensuses[0].scores(), &[1, 0, 1]);
```
*/

use log::{debug, trace};
use priority_queue::PriorityQueue;
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};
use simple_error::bail;
use std::cmp::Reverse;

use crate::cdwfa_config::{CdwfaConfig, ConsensusCost};
use crate::dynamic_wfa::DWFALite;
use crate::pqueue_tracker::PQueueTracker;

type NodePriority = (Reverse<usize>, usize);

/// Contains a final consensus result
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Consensus {
    /// The generated consensus
    sequence: Vec<u8>,
    /// The consensus scoring model
    consensus_cost: ConsensusCost,
    /// Vector of the scores from the consensus to each sequence
    scores: Vec<usize>
}

impl Consensus {
    /// Constructor
    pub fn new(sequence: Vec<u8>, consensus_cost: ConsensusCost, scores: Vec<usize>) -> Consensus {
        Consensus {
            sequence,
            consensus_cost,
            scores
        }
    }

    // Getters
    pub fn sequence(&self) -> &[u8] {
        &self.sequence
    }

    pub fn consensus_cost(&self) -> ConsensusCost {
        self.consensus_cost
    }

    pub fn scores(&self) -> &[usize] {
        &self.scores
    }
}

/// Core utility that will generate a consensus sequence using the DWFA approach.
/// For now, it will assume that all sequences are full length representations
#[derive(Debug, Default)]
pub struct ConsensusDWFA<'a>  {
    /// Contains all the sequences that have been added to this consensus so far.
    sequences: Vec<&'a [u8]>,
    /// Last offset into the consensus that the sequence is allowed to start, this may get adjusted at run-time. If None, then assume start.
    offsets: Vec<Option<usize>>,
    /// The config for this consensus run
    config: CdwfaConfig,
    /// The alphabet we will use for consensus building
    alphabet: HashSet<u8>
}

impl<'a> ConsensusDWFA<'a> {
    /// Creates a new instance of ConsensusDWFA and performs sanity checks.
    /// # Arguments
    /// * `config` - the type of consensus score we want to use
    /// # Errors
    /// * None so far
    pub fn with_config(config: CdwfaConfig) -> Result<ConsensusDWFA<'a>, Box<dyn std::error::Error>> {
        Ok(ConsensusDWFA {
            sequences: vec![],
            offsets: vec![],
            config,
            alphabet: Default::default()
        })
    }

    /// Adds a new sequence to the list without an offset.
    /// # Arguments
    /// * `sequence` - the new sequence to add
    /// # Errors
    /// * None so far
    pub fn add_sequence(&mut self, sequence: &'a [u8]) -> Result<(), Box<dyn std::error::Error>> {
        self.add_sequence_offset(sequence, None)
    }

    /// Adds a new sequence to the list 
    /// # Arguments
    /// * `sequence` - the new sequence to add
    /// * `last_offset` - an optional offset into the consensus that is the last allowed offset for the sequence to start
    /// # Errors
    /// * None so far
    pub fn add_sequence_offset(&mut self, sequence: &'a [u8], last_offset: Option<usize>) -> Result<(), Box<dyn std::error::Error>> {
        // add ever character to our alphabet
        self.alphabet.extend(sequence.iter().cloned());
        
        // make sure we didn't add the wildcard by accident
        // note: it's easier to just remove it once here than do something fancy with the above iteration
        if let Some(wc) = self.config.wildcard {
            self.alphabet.remove(&wc);
        }

        // now add the sequence
        self.sequences.push(sequence);
        self.offsets.push(last_offset);
        Ok(())
    }

    /// The core function that gets called after adding all the sequences we care about
    /// # Errors
    /// * if the DWFA throws errors
    pub fn consensus(&self) -> Result<Vec<Consensus>, Box<dyn std::error::Error>> {
        // initialize everything
        let mut maximum_error = usize::MAX;
        let mut nodes_explored: usize = 0;
        let mut nodes_ignored: usize = 0;
        let mut peak_queue_size: usize = 0;
        let mut farthest_consensus: usize = 0;
        let mut last_constraint: u64 = 0;

        let offset_window = self.config.offset_window; // the window around the provided offset that we check
        let offset_compare_length = self.config.offset_compare_length; // the number of bases we check to figure out the real start

        let offsets: Vec<Option<usize>> = if self.config.auto_shift_offsets {
            // first figure out if we need to address the offsets for this re-run of consensus
            let mut min_offset = usize::MAX;
            let mut start_sequence_found = false;
            for offset in self.offsets.iter() {
                match offset {
                    Some(o) => min_offset = min_offset.min(*o),
                    None => start_sequence_found = true
                };
            }

            if !start_sequence_found {
                debug!("No start sequence detected, shifting all offsets by {min_offset}");
                self.offsets.iter()
                    .map(|o| {
                        let value = o.unwrap();
                        if value == min_offset {
                            // equal to minimum, flag as no offset
                            None
                        } else {
                            // past minimum, so down-shift it
                            Some(value - min_offset)
                        }
                    })
                    .collect::<Vec<Option<usize>>>()
            } else {
                self.offsets.clone()
            }
        } else {
            self.offsets.clone()
        };

        debug!("Offsets: {:?}", offsets);

        // build up the list of sizes where we need to activate one or more sequences
        let mut activate_points: HashMap<usize, Vec<usize>> = Default::default();
        let mut max_activate = 0;
        let mut initially_active: usize = 0;
        for (seq_index, offset) in offsets.iter().enumerate() {
            if let Some(last_offset) = offset {
                // we have an approx, calculate when it can get activated and add it to our list
                let activate_length = last_offset + offset_compare_length;
                let entry = activate_points.entry(activate_length).or_default();
                entry.push(seq_index);

                if activate_length > max_activate {
                    max_activate = activate_length;
                }
            } else {
                initially_active += 1;
            }
        }

        if initially_active == 0 {
            bail!("Must have at least one initial offset of None to see the consensus.");
        }
        
        let max_queue_size = self.config.max_queue_size;
        let max_capacity_per_size = self.config.max_capacity_per_size;
        let initial_size = self.sequences.iter().map(|s| s.len()).max().unwrap();
        let mut pqueue_tracker = PQueueTracker::with_capacity(initial_size, max_capacity_per_size);

        let initial_node = ConsensusNode::new_root_node(&offsets, self.config.wildcard, self.config.allow_early_termination)?;
        let initial_priority = initial_node.priority(self.consensus_cost());
        
        // start the priority queue, which defaults to bigger is better so we need a Reverse since want smaller costs
        let mut pqueue: PriorityQueue<ConsensusNode, NodePriority> = PriorityQueue::new();
        pqueue_tracker.insert(initial_node.consensus().len());
        pqueue.push(initial_node, initial_priority);

        let mut ret = vec![];

        // the way this will work is that we will eventually find one or more answers and anything worse will get drained off until no possibilities remain
        while !pqueue.is_empty() {
            // this just tracks how large our actual queue gets
            peak_queue_size = peak_queue_size.max(pqueue.len());

            // first, check if we need to restrict our pqueue threshold
            while (pqueue_tracker.len() > max_queue_size || last_constraint >= self.config.max_nodes_wo_constraint) && pqueue_tracker.threshold() < farthest_consensus {
                pqueue_tracker.increment_threshold();
                last_constraint = 0;
            }

            // get the top node and check if it's bad
            let (top_node, top_priority) = pqueue.pop().unwrap();
            let top_len = top_node.consensus().len();
            pqueue_tracker.remove(top_len);

            trace!("Pop: {:?} => {:?}", top_priority, top_node.consensus);
        
            if top_priority.0.0 > maximum_error || top_priority.1 < pqueue_tracker.threshold() || pqueue_tracker.at_capacity(top_len) {
                trace!("\tIgnored");
                nodes_ignored += 1;
                continue;
            }

            // mark this one as explored and updated farthest if necessary
            farthest_consensus = farthest_consensus.max(top_len);
            nodes_explored += 1;
            last_constraint += 1;
            pqueue_tracker.process(top_len)?;

            // now check if this node has reached the end
            if top_node.reached_end(&self.sequences, self.config.allow_early_termination) {
                // this node *might* be done, we need to finalize it to be sure
                // we also need to do it in a clone since it might not be finalized and we need to keep extending
                let mut finalized_node = top_node.clone();
                finalized_node.finalize(&self.sequences)?;

                // get the new finalized score, it may have changed
                let finalized_score = finalized_node.total_cost(self.config.consensus_cost);

                // check first if it's stricly BETTER than anything so far
                if finalized_score < maximum_error {
                    // this score is better than anything we have seen so far, clear out previous results if we have any
                    maximum_error = finalized_score;
                    ret.clear();
                }

                // now check if it's as good as anything so far; the above check will not nullify this one
                if finalized_score <= maximum_error && ret.len() < self.config.max_return_size {
                    // this consensus is at least as good as the best, so add it as a result
                    ret.push(Consensus::new(
                        finalized_node.consensus().to_vec(),
                        self.config.consensus_cost,
                        finalized_node.costs(self.config.consensus_cost)
                    ));
                }
            }

            // this fetches the list of options according to the WFA so far
            // NOTE: this CAN include the wildcard, but only if the wildcard is the only character
            let extension_candidates = top_node.get_extension_candidates(&self.sequences, self.config.wildcard);
            let max_observed = extension_candidates.values().cloned().max_by(|a, b| a.total_cmp(b))
                // if no observation, then just use min count
                .unwrap_or(self.config.min_count as f64);
            
            // the active threshold is the minimum of 1) the configured minimum count OR 2) the highest count we observed
            let active_threshold = (self.config.min_count as f64).min(max_observed);
            trace!("\tcandidates = {:?}", extension_candidates);

            // pull out the full list of viable candidates
            let passing_candidates: Vec<u8> = extension_candidates.into_iter()
                .filter_map(|(symbol, count)| {
                    if count >= active_threshold {
                        Some(symbol)
                    } else {
                        None
                    }
                }).collect();
            
            let mut new_nodes = vec![];
            if passing_candidates.is_empty() {
                if top_len < max_activate {
                    bail!("Encountered coverage gap: consensus is length {} with no candidates, but sequences activate at {}", top_len, max_activate);
                } else {
                    // no extensions remain, should just happen at the end
                }
            } else if passing_candidates.len() == 1 {
                // we have only one extension, we can extend in place without cloning
                let mut new_node = top_node;
                new_node.push(&self.sequences, passing_candidates[0])?;
                new_nodes.push(new_node);
            } else {
                // we have 2+ viable, we have to clone them
                for symbol in passing_candidates.into_iter() {
                    let mut new_node = top_node.clone();
                    new_node.push(&self.sequences, symbol)?;
                    new_nodes.push(new_node);
                }
            }

            // for each need node, do any activations and then add it to the queue
            for mut new_node in new_nodes.into_iter() {
                // check if we need to activate any strings
                let opt_activate_list = activate_points.get(&new_node.consensus().len());
                if let Some(activate_list) = opt_activate_list {
                    assert!(!activate_list.is_empty());
                    for &seq_index in activate_list.iter() {
                        new_node.activate_sequence(self.sequences[seq_index], seq_index, offset_window, offset_compare_length, self.config.wildcard, self.config.allow_early_termination)?;
                    }
                }

                // get the new cost and put it in the queue
                let new_priority = new_node.priority(self.consensus_cost());
                trace!("\tPush {:?} => {:?}", new_priority, new_node.consensus);
                pqueue_tracker.insert(new_node.consensus().len());
                pqueue.push(new_node, new_priority);
            }
        }

        assert_eq!(pqueue_tracker.len(), 0);

        // sort these by the sequence so we always have a fixed order
        ret.sort_by(|c1, c2| c1.sequence().cmp(c2.sequence()));

        debug!("nodes_explored: {nodes_explored}");
        debug!("nodes_ignored: {nodes_ignored}");
        debug!("peak_queue_size: {peak_queue_size}");
        Ok(ret)
    }

    // getters
    pub fn sequences(&self) -> &[&'a [u8]] {
        &self.sequences
    }

    pub fn consensus_cost(&self) -> ConsensusCost {
        self.config.consensus_cost
    }

    pub fn alphabet(&self) -> &HashSet<u8> {
        &self.alphabet
    }
}

/// Wrapper for a node containing a partial consensus as well as the DWFA tracking for that node
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
struct ConsensusNode {
    /// The consensus sequence so far
    consensus: Vec<u8>,
    /// The DWFAs that are tracked for each sequence; these are only None if they have not been activated yet due to an offset
    dwfas: Vec<Option<DWFALite>>
}

impl ConsensusNode {
    /// Constructor for a new consensus search root node
    /// # Arguments
    /// * `offsets` - a set of offsets into the sequences where the approximate starts are
    /// * `wildcard` - an optional wildcard symbol that will match anything
    /// * `allow_early_termination` - if true, then it will allow the consensus to go beyond the provided baseline sequences without penalty
    /// # Errors
    /// * if DWFA construction fails
    fn new_root_node(offsets: &[Option<usize>], wildcard: Option<u8>, allow_early_termination: bool) -> Result<ConsensusNode, Box<dyn std::error::Error>> {
        let dwfas: Vec<Option<DWFALite>> = offsets.iter()
            .map(|offset| match offset {
                // we have an offset, so do not create a map
                Some(_o) => None,
                // we don't have an offset, so this is active from the start
                None => Some(DWFALite::new(wildcard, allow_early_termination))
            })
            .collect();

        let active_count = dwfas.iter().filter(|d| d.is_some()).count();
        if active_count == 0 {
            bail!("Root consensus node has no initially active sequences.");
        }

        Ok(ConsensusNode {
            consensus: vec![],
            dwfas
        })
    }

    /// This will activate the DWFAs for a particular sequence.
    /// # Arguments
    /// * `sequence` - the sequence that is getting activated
    /// * `seq_index` - the sequence index, which will map to the DWFAs
    /// * `offset_window` - the window we are searching for a best match
    /// * `offset_compare_length` - the amount of bases we are comparing
    /// * `wildcard` - optional wildcard for scoring, passed into DWFA
    /// * `allow_early_termination` - enables sequences to end partway through the consensus, passed into DWFA
    fn activate_sequence(&mut self, sequence: &[u8], seq_index: usize, offset_window: usize, offset_compare_length: usize, wildcard: Option<u8>, allow_early_termination: bool) -> Result<(), Box<dyn std::error::Error>> {
        // make sure everything is currently inactive
        assert!(self.dwfas[seq_index].is_none());

        // build up the search space
        let con_len = self.consensus.len();
        
        // figure out how far back we search
        let start_delta = offset_window + offset_compare_length; // we search from current back to the offset_window
        let start_position = con_len.saturating_sub(start_delta);
        let end_position = con_len.saturating_sub(offset_compare_length);

        // figure out which offset has the best score; assume the middle of the offset window is the best
        let mut best_offset = con_len.saturating_sub(offset_compare_length + offset_window / 2);
        let mut min_ed = crate::sequence_alignment::wfa_ed_config(&self.consensus[best_offset..], &sequence[0..offset_compare_length], false, wildcard);
        
        // now check all the rest around this position
        for p in start_position..end_position {
            let ed = crate::sequence_alignment::wfa_ed_config(&self.consensus[p..], &sequence[0..offset_compare_length], false, wildcard);
            if ed < min_ed {
                min_ed = ed;
                best_offset = p;
            }
        }

        // now set up the DWFA with the best offset
        let mut new_dwfa = DWFALite::new(wildcard, allow_early_termination);
        new_dwfa.set_offset(best_offset);
        new_dwfa.update(sequence, &self.consensus)?;
        self.dwfas[seq_index] = Some(new_dwfa);

        Ok(())
    }

    /// Adds a new symbol to the DWFAs
    /// # Arguments
    /// * `symbol` - the symbol to extend each DWFA with
    /// # Errors
    /// * if any of the extensions fail
    fn push(&mut self, baseline_sequences: &[&[u8]], symbol: u8) -> Result<(), Box<dyn std::error::Error>> {
        self.consensus.push(symbol);
        for (&baseline, opt_dwfa) in baseline_sequences.iter().zip(self.dwfas.iter_mut()) {
            if let Some(dwfa) = opt_dwfa.as_mut() {
                dwfa.update(baseline, &self.consensus)?;
            }
        }
        Ok(())
    }

    /// Finalizes all of the contained DWFAs
    /// # Errors
    /// * if any of the individual finalizes fail
    fn finalize(&mut self, baseline_sequences: &[&[u8]]) -> Result<(), Box<dyn std::error::Error>> {
        for (&baseline, opt_dwfa) in baseline_sequences.iter().zip(self.dwfas.iter_mut()) {
            if let Some(dwfa) = opt_dwfa.as_mut() {
                dwfa.finalize(baseline, &self.consensus)?;
            } else {
                bail!("Finalize called on DWFA that was never initialized.");
            }
        }
        Ok(())
    }

    /// Returns a vec of all of the scores for this node
    fn costs(&self, consensus_cost: ConsensusCost) -> Vec<usize> {
        self.dwfas.iter()
            .map(|opt_d| {
                if let Some(d) = opt_d.as_ref() {
                    match consensus_cost {
                        ConsensusCost::L1Distance => d.edit_distance(),
                        ConsensusCost::L2Distance => d.edit_distance().pow(2)
                    }
                } else {
                    // it's not initialized, here we will just report 0
                    0
                }
            })
            .collect()
    }

    /// Returns the total score for the node
    fn total_cost(&self, consensus_cost: ConsensusCost) -> usize {
        self.costs(consensus_cost).iter().sum()
    }

    /// Returns the node priority.
    /// Currently, this is based on 1) lowest cost and 2) consensus length.
    /// # Arguments
    /// * `consensus_cost` - cost model to evaluate the cost
    fn priority(&self, consensus_cost: ConsensusCost) -> NodePriority {
        (
            Reverse(self.total_cost(consensus_cost)),
            self.consensus.len()
        )
    }

    /// Returns true if any of the dwfas are at the end of their respective baseline sequences
    /// # Arguments
    /// * `baseline_sequences` - the sequences that are fixed that we want the consensus of
    /// * `require_all` - if True, requires all of the sequences to be at their end (this is useful for allow_early_termination)
    fn reached_end(&self, baseline_sequences: &[&[u8]], require_all: bool) -> bool {
        let mut end_iter = baseline_sequences.iter().zip(self.dwfas.iter())
            .map(|(&baseline, opt_dwfa)| {
                if let Some(dwfa) = opt_dwfa.as_ref() {
                    dwfa.reached_baseline_end(baseline)
                } else {
                    false
                }
            });
        
        // same iterator, different check
        if require_all {
            end_iter.all(|b| b)
        } else {
            end_iter.any(|b| b)
        }
    }

    /// Returns a hashmap of the extension candidates with their votes.
    /// The votes are fractional if a sequence has multiple equally likely options.
    /// Note that if 
    /// # Arguments
    /// * `baseline_sequences` - the sequences that are fixed that we want the consensus of
    /// * `wildcard` - an optional wildcard character, will be removed from return set unless it is the only value in it
    fn get_extension_candidates(&self, baseline_sequences: &[&[u8]], wildcard: Option<u8>) -> HashMap<u8, f64> {
        let mut candidates: HashMap<u8, f64> = Default::default();
        for (baseline_seq, opt_dwfa) in baseline_sequences.iter().zip(self.dwfas.iter()) {
            if let Some(dwfa) = opt_dwfa.as_ref() {
                // get the candidates and the total observation weight
                let cand = dwfa.get_extension_candidates(baseline_seq, &self.consensus);
                let vote_split = cand.values().sum::<usize>() as f64;
                
                // iterate over each candidate and scale it by the occurrences count / total weight
                for (&c, &occ) in cand.iter() {
                    let entry = candidates.entry(c).or_insert(0.0);
                    *entry += occ as f64 / vote_split;
                }
            }
        }

        if let Some(wc) = wildcard {
            if candidates.len() > 1 {
                // we have 2+ candidates, we can safely remove the wildcard and have something left over
                candidates.remove(&wc);
            }
        }

        candidates
    }

    // getters
    fn consensus(&self) -> &[u8] {
        &self.consensus
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::cdwfa_config::CdwfaConfigBuilder;

    #[test]
    fn test_single_sequence() {
        let sequence = b"ACGTACGTACGT";
        let mut consensus_dwfa = ConsensusDWFA::default();
        consensus_dwfa.add_sequence(sequence).unwrap();

        // make sure our alphabet has four symbols, ACGT
        assert_eq!(consensus_dwfa.alphabet().len(), 4);

        // now check that the consensus is the same as our sequence
        let consensus = consensus_dwfa.consensus().unwrap();
        assert_eq!(consensus, vec![Consensus {
            sequence: sequence.to_vec(),
            consensus_cost: ConsensusCost::L1Distance,
            scores: vec![0]
        }]);
    }

    #[test]
    fn test_dual_sequence() {
        let sequence = b"ACGTACGTACGT";
        let sequence2 = b"ACGTACCTACGT";
        let mut consensus_dwfa = ConsensusDWFA::default();
        consensus_dwfa.add_sequence(sequence).unwrap();
        consensus_dwfa.add_sequence(sequence2).unwrap();

        // make sure our alphabet has four symbols, ACGT
        assert_eq!(consensus_dwfa.alphabet().len(), 4);

        // now check that the consensus is a tie between the two we put it, each 1 ED total away
        let consensus = consensus_dwfa.consensus().unwrap();

        // sequence2 is alphabetically before sequence 1, so it will come first
        assert_eq!(consensus, vec![
            Consensus {
                sequence: sequence2.to_vec(),
                consensus_cost: ConsensusCost::L1Distance,
                scores: vec![1, 0]
            },
            Consensus {
                sequence: sequence.to_vec(),
                consensus_cost: ConsensusCost::L1Distance,
                scores: vec![0, 1]
            },
        ]);
    }

    #[test]
    fn test_trio_sequence() {
        let sequence = b"ACGTACGTACGT";
        let sequence2 = b"ACGTACCTACGT";
        let mut consensus_dwfa = ConsensusDWFA::default();

        // add the first sequence twice, it will be the consensus
        consensus_dwfa.add_sequence(sequence).unwrap();
        consensus_dwfa.add_sequence(sequence).unwrap();
        consensus_dwfa.add_sequence(sequence2).unwrap();

        // make sure our alphabet has four symbols, ACGT
        assert_eq!(consensus_dwfa.alphabet().len(), 4);

        // now check that the consensus found sequence
        let consensus = consensus_dwfa.consensus().unwrap();

        // sequence2 is alphabetically before sequence 1, so it will come first
        assert_eq!(consensus, vec![
            Consensus {
                sequence: sequence.to_vec(),
                consensus_cost: ConsensusCost::L1Distance,
                scores: vec![0, 0, 1]
            }
        ]);
    }

    #[test]
    fn test_complicated() {
        let expected_consensus = b"ACGTACGTACGT";
        let sequences = [
            // b"ACGTACG-TACGT"
            // b"AC-TACGGTACGT"
            b"ACTACGGTACGT",
            // b"ACGTAAG-TCCGT"
            b"ACGTAAGTCCGT",
            // b"AAGTACG-TACGT"
            b"AAGTACGTACGT"
        ];
        
        // build a consensus over all inputs
        let mut consensus_dwfa = ConsensusDWFA::default();
        for &sequence in sequences.iter() {
            consensus_dwfa.add_sequence(sequence).unwrap();
        }

        // make sure our alphabet has four symbols, ACGT
        assert_eq!(consensus_dwfa.alphabet().len(), 4);

        // now check that the consensus found the right thing and only that
        let consensus = consensus_dwfa.consensus().unwrap();
        assert_eq!(consensus.len(), 1);
        assert_eq!(expected_consensus, consensus[0].sequence());
    }

    #[test]
    fn test_wildcards() {
        let expected_consensus = b"ACGTACGTACGT";
        let sequences = [
            // C insertion and tail missing
            b"ACGTACCGT****".to_vec(),
            // C>T and head+tail missing
            b"**GTATGTAC**".to_vec(),
            // just head missing
            b"****ACGTACGT".to_vec()
        ];

        // build a consensus over all inputs
        let mut consensus_dwfa = ConsensusDWFA::with_config(
            CdwfaConfigBuilder::default()
                .wildcard(Some(b'*'))
                .build().unwrap()
        ).unwrap();
        for sequence in sequences.iter() {
            consensus_dwfa.add_sequence(sequence).unwrap();
        }

        // make sure our alphabet has four symbols, ACGT
        assert_eq!(consensus_dwfa.alphabet().len(), 4);

        // now check that the consensus found the right thing and only that
        let consensus = consensus_dwfa.consensus().unwrap();
        assert_eq!(consensus.len(), 1);
        assert_eq!(expected_consensus, consensus[0].sequence());
        assert_eq!(&[1, 1, 0], consensus[0].scores());
    }

    #[test]
    fn test_all_wildcards() {
        let _expected_consensus = b"ACGTACGTACGT";
        let actual_consensus = b"*CGTACG*ACG*";
        let sequences = [
            // 1 ins
            b"*CGTAACG*ACG*".to_vec(),
            // nothing other than *
            b"*CGTACG*ACG*".to_vec(),
            // C>T
            b"*CGTACG*ATG*".to_vec()
        ];

        // build a consensus over all inputs
        let mut consensus_dwfa = ConsensusDWFA::with_config(
            CdwfaConfigBuilder::default()
                .wildcard(Some(b'*'))
                .build().unwrap()
        ).unwrap();
        for sequence in sequences.iter() {
            consensus_dwfa.add_sequence(sequence).unwrap();
        }

        // make sure our alphabet has four symbols, ACGT
        assert_eq!(consensus_dwfa.alphabet().len(), 4);

        // now check that the consensus found the right thing and only that
        let consensus = consensus_dwfa.consensus().unwrap();
        assert_eq!(consensus.len(), 1);
        assert_eq!(actual_consensus, consensus[0].sequence());
        assert_eq!(&[1, 0, 1], consensus[0].scores());
    }

    #[test]
    fn test_allow_early_termination_costs() {
        let expected_consensus = b"ACGT";

        // first verify baseline is not generating the full thing due to error
        let mut consensus_dwfa = ConsensusDWFA::with_config(
            CdwfaConfigBuilder::default()
                .wildcard(Some(b'*'))
                .build().unwrap()
        ).unwrap();
        for i in 1..(expected_consensus.len()+1) {
            consensus_dwfa.add_sequence(&expected_consensus[0..i]).unwrap();
        }
        // this first approach generated multiple possible ones that are in the middle
        let consensus = consensus_dwfa.consensus().unwrap();
        assert_eq!(consensus, [
            Consensus { sequence: vec![65, 67], consensus_cost: ConsensusCost::L1Distance, scores: vec![1, 0, 1, 2] }, 
            Consensus { sequence: vec![65, 67, 71], consensus_cost: ConsensusCost::L1Distance, scores: vec![2, 1, 0, 1] }
        ]);

        // second, verify that allowing early termination fixes it to the original result
        let mut consensus_dwfa = ConsensusDWFA::with_config(
            CdwfaConfigBuilder::default()
                .wildcard(Some(b'*'))
                .allow_early_termination(true)
                .build().unwrap()
        ).unwrap();
        for i in 1..(expected_consensus.len()+1) {
            consensus_dwfa.add_sequence(&expected_consensus[0..i]).unwrap();
        }
        // this first approach generated multiple possible ones that are in the middle
        let consensus = consensus_dwfa.consensus().unwrap();
        assert_eq!(consensus, [
            Consensus { sequence: expected_consensus.to_vec(), consensus_cost: ConsensusCost::L1Distance, scores: vec![0; 4] },
        ]);
    }

    #[test]
    fn test_offset_windows() {
        // this test is _only_ for verifying that our offsets are getting correctly adjust to match
        let expected_consensus = b"ACGTACGTACGTACGT";
        let sequences = [
            // no offset
            b"ACGTACGTACGTACGT".to_vec(),
            // offset 4
            b"ACGTACGTACGT".to_vec(),
            // offset 6
            b"GTACGTACGT".to_vec(),
        ];

        // make one offset exactly match and another be right on our boundary (1-shift here)
        let offsets = [None, Some(4), Some(7)];

        // build a consensus over all inputs
        let mut consensus_dwfa = ConsensusDWFA::with_config(
            CdwfaConfigBuilder::default()
                .offset_window(1) // both are exactly one off
                .offset_compare_length(4) // semi-arbitrary, just needs to be short here given our seq-len
                .build().unwrap()
        ).unwrap();
        for (sequence, offset) in sequences.iter().zip(offsets.iter()) {
            consensus_dwfa.add_sequence_offset(sequence, offset.clone()).unwrap();
        }

        // make sure our alphabet has four symbols, ACGT
        assert_eq!(consensus_dwfa.alphabet().len(), 4);

        // now check that the consensus found the right thing and only that
        let consensus = consensus_dwfa.consensus().unwrap();
        assert_eq!(consensus.len(), 1);
        assert_eq!(expected_consensus, consensus[0].sequence());
        assert_eq!(&[0, 0, 0], consensus[0].scores());
    }

    #[test]
    fn test_offset_gap_err() {
        // this test is _only_ for verifying that our offsets are getting correctly adjust to match
        let sequences = [
            // no offset
            b"ACGTACGTACGTACGT".to_vec(),
            // massive offset
            b"ACGTACGTACGTACGT".to_vec(),
        ];

        // one starts and one is way out there
        let offsets = [None, Some(1000)];

        // build a consensus over all inputs
        let mut consensus_dwfa = ConsensusDWFA::with_config(
            CdwfaConfigBuilder::default()
                .offset_window(1) // both are exactly one off
                .offset_compare_length(4) // semi-arbitrary, just needs to be short here given our seq-len
                .build().unwrap()
        ).unwrap();
        for (sequence, offset) in sequences.iter().zip(offsets.iter()) {
            consensus_dwfa.add_sequence_offset(sequence, offset.clone()).unwrap();
        }

        // make sure our alphabet has four symbols, ACGT
        assert_eq!(consensus_dwfa.alphabet().len(), 4);

        // now check that the consensus found the right thing and only that
        let consensus_err = consensus_dwfa.consensus();
        assert!(consensus_err.is_err());
        let error = consensus_err.err().unwrap();
        assert_eq!(error.to_string(), "Finalize called on DWFA that was never initialized.");
    }
}