
/*!
This module provides access to the DualConsensusDWFA, which generates the two best consensuses for a set of sequences.
Note that it can also generate a single consensus if no splits are identified.

# Example usage
```rust
use waffle_con::consensus::Consensus;
use waffle_con::dual_consensus::DualConsensusDWFA;
use waffle_con::cdwfa_config::ConsensusCost;

let sequences = [
    b"TCCGT".to_vec(),
    b"ACCGT".to_vec(), // consensus 1
    b"ACCGT".to_vec(), // consensus 1
    b"ACCAT".to_vec(),
    b"CCGTAAT".to_vec(),
    b"CGTAAAT".to_vec(),
    b"CGTAAT".to_vec(), // consensus 2
    b"CGTAAT".to_vec(), // consensus 2
];

// add all the sequences
let mut cdwfa: DualConsensusDWFA = Default::default();
for s in sequences.iter() {
    cdwfa.add_sequence(s).unwrap();
}

// run consensus and check the results
let consensuses = cdwfa.consensus().unwrap();
assert_eq!(consensuses.len(), 1);
assert_eq!(consensuses[0].consensus1(), &Consensus::new(sequences[1].clone(), ConsensusCost::L1Distance, vec![1, 0, 0, 1]));
assert_eq!(consensuses[0].consensus2().unwrap(), &Consensus::new(sequences[6].clone(), ConsensusCost::L1Distance, vec![1, 1, 0, 0]));
assert_eq!(consensuses[0].is_consensus1(), &[true, true, true, true, false, false, false, false]);
```
*/

use log::{debug, trace, warn};
use priority_queue::PriorityQueue;
use rustc_hash::{FxHashMap as HashMap, FxHashSet as HashSet};
use simple_error::{bail, SimpleError};
use std::cmp::Reverse;

use crate::cdwfa_config::{CdwfaConfig, ConsensusCost};
use crate::consensus::Consensus;
use crate::dynamic_wfa::DWFALite;
use crate::pqueue_tracker::PQueueTracker;

type NodePriority = (Reverse<usize>, usize);

/// Contains a final multi-consensus result
#[derive(Debug)]
pub struct DualConsensus {
    /// The first consensus
    consensus1: Consensus,
    /// The second consensus - this one can be optional when only one consensus is identified
    consensus2: Option<Consensus>,
    /// same length as input sequences; for each entry, if True, then the corresponding input sequence matches consensus #1; otherwise, it matches #2
    is_consensus1: Vec<bool>,
    /// The scores when compared to consensus1, will be None if it stopped tracking
    scores1: Vec<Option<usize>>,
    /// The scores when compared to consensus2, will be None if it stopped tracking
    scores2: Vec<Option<usize>>
}

impl PartialEq for DualConsensus {
    fn eq(&self, other: &Self) -> bool {
        self.consensus1 == other.consensus1 && 
            self.consensus2 == other.consensus2 && 
            self.is_consensus1 == other.is_consensus1 

            // for now, we do not care about these in the impl
            //&& self.scores1 == other.scores1 && self.scores2 == other.scores2
    }
}

impl DualConsensus {
    /// Generic constructor for outside usage. This does not force consensus order like `from_node(...)`.
    pub fn new(
        consensus1: Consensus, consensus2: Option<Consensus>, 
        is_consensus1: Vec<bool>, scores1: Vec<Option<usize>>, scores2: Vec<Option<usize>>
    ) -> Result<Self, SimpleError> {

        if is_consensus1.len() != scores1.len() || is_consensus1.len() != scores2.len() {
            bail!("is_consensus1, scores1, and scores2 must all be the same length");
        }

        Ok(Self {
            consensus1,
            consensus2,
            is_consensus1,
            scores1,
            scores2,
        })
    }

    /// Constructor from a DualConsensusNode.
    /// Of note, this will re-order dual alleles alpha-numerically, allowing for predictable outputs.
    /// # Arguments
    /// * `finalized_node` - the node that we are converting to a consensus
    /// * `consensus_cost` - the cost model that gets propated
    fn from_node(finalized_node: &DualConsensusNode, consensus_cost: ConsensusCost) -> DualConsensus {
        // figure out the best consensus and score for each read; note that this assumes no ties
        let (best_consensus_index, best_consensus_score) = finalized_node.costs(consensus_cost);

        // check if we need to swap the order
        let swap_order = finalized_node.is_dual && (finalized_node.consensus2 < finalized_node.consensus1);

        // now reformat the above information such that we can build out the multi-consensus return values
        let mut is_consensus1: Vec<bool> = Default::default();
        let mut consensus_scores: Vec<Vec<usize>> = vec![vec![]; 2];
        for (best_con_index, best_con_score) in best_consensus_index.into_iter()
            .zip(best_consensus_score.into_iter()) {
            // these MUST be equal length
            // store that this sequence index matches the particular consensus index
            assert!(best_con_index <= 1);
            // toggle the consensus assignment if we swap_order
            is_consensus1.push((best_con_index == 0) ^ swap_order);
            // also store this score for this consensus index
            consensus_scores[best_con_index].push(best_con_score);
        }

        // now we can store the consensus sequences as well as the corresponding indices in the final output
        let c1 = Consensus::new(finalized_node.consensus1.clone(), consensus_cost, consensus_scores[0].clone());
        let c2 = Consensus::new(finalized_node.consensus2.clone(), consensus_cost, consensus_scores[1].clone());

        // reformat the actual consensus assignments based on swappage, and build result
        let (consensus1, consensus2) = if swap_order {
            assert!(finalized_node.is_dual);
            (c2, Some(c1))
        } else {
            (c1, if finalized_node.is_dual { Some(c2) } else { None })
        };

        // now save the scores also
        let (s1, s2) = finalized_node.full_cost(consensus_cost);
        let (scores1, scores2) = if swap_order {
            (s2, s1)
        } else {
            (s1, s2)
        };

        DualConsensus {
            consensus1,
            consensus2,
            is_consensus1,
            scores1,
            scores2
        }
    }

    // Returns true if this is a dual consensus result
    pub fn is_dual(&self) -> bool {
        self.consensus2.is_some()
    }

    // Getters
    pub fn consensus1(&self) -> &Consensus {
        &self.consensus1
    }

    pub fn consensus2(&self) -> Option<&Consensus> {
        self.consensus2.as_ref()
    }

    pub fn is_consensus1(&self) -> &[bool] {
        &self.is_consensus1
    }

    pub fn scores1(&self) -> &[Option<usize>] {
        &self.scores1
    }

    pub fn scores2(&self) -> &[Option<usize>] {
        &self.scores2
    }
}

/// Core utility that will generate a consensus sequence using the DWFA approach.
/// For now, it will assume that all sequences are full length representations
#[derive(Debug, Default)]
pub struct DualConsensusDWFA<'a>  {
    /// Contains all the sequences that have been added to this consensus so far.
    sequences: Vec<&'a [u8]>,
    /// Approximate offsets into the consensus that sequence starts, this will get adjust at run-time. If None, then assume start.
    offsets: Vec<Option<usize>>,
    /// The config for this consensus run
    config: CdwfaConfig,
    /// The alphabet we will use for consensus building
    alphabet: HashSet<u8>
}

impl<'a> DualConsensusDWFA<'a> {
    /// Creates a new instance of ConsensusDWFA with the specified config.
    /// # Arguments
    /// * `config` - the type of consensus score we want to use
    /// # Errors
    /// * None so far
    pub fn with_config(config: CdwfaConfig) -> Result<DualConsensusDWFA<'a>, Box<dyn std::error::Error>> {
        Ok(DualConsensusDWFA {
            config,
            ..Default::default()
        })
    }

    /// Adds a new sequence to the list with no offset.
    /// # Arguments
    /// * `sequence` - the new sequence to add
    /// # Errors
    /// * None so far
    pub fn add_sequence(&mut self, sequence: &'a [u8]) -> Result<(), Box<dyn std::error::Error>> {
        self.add_sequence_offset(sequence, None)
    }

    /// Adds a new sequence to the list with no offset.
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

    // The core function that gets called after adding all the sequences we care about
    /// # Errors
    /// * if the DWFA throws errors
    pub fn consensus(&self) -> Result<Vec<DualConsensus>, Box<dyn std::error::Error>> {
        // initialize everything
        let mut maximum_error = usize::MAX;
        let mut nodes_explored: usize = 0;
        let mut nodes_ignored: usize = 0;
        let mut peak_queue_size: usize = 0;
        let mut farthest_consensus_single: usize = 0;
        let mut farthest_consensus_dual: usize = 0;
        let mut single_last_constraint: u64 = 0;
        let mut dual_last_constraint: u64 = 0;

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
        let mut initially_active: usize = 0;
        for (seq_index, offset) in offsets.iter().enumerate() {
            if let Some(last_offset) = offset {
                // we have an approx, calculate when it can get activated and add it to our list
                let activate_length = last_offset + offset_compare_length;
                let entry = activate_points.entry(activate_length).or_default();
                entry.push(seq_index);
            } else {
                initially_active += 1;
            }
        }

        if initially_active == 0 {
            bail!("Must have at least one initial offset of None to see the consensus.");
        }
        
        // we now track singletons and dual nodes separately so we can guarantee that we return _something_ even if final dual nodes are imbalanced
        let max_queue_size = self.config.max_queue_size;
        let max_capacity_per_size = self.config.max_capacity_per_size;
        let initial_size = self.sequences.iter().map(|s| s.len()).max().unwrap();
        let mut single_tracker = PQueueTracker::with_capacity(initial_size, max_capacity_per_size);
        let mut dual_tracker = PQueueTracker::with_capacity(initial_size, max_capacity_per_size);

        let initial_node = DualConsensusNode::new_root_node(&offsets, self.config.wildcard, self.config.allow_early_termination)?;
        let initial_priority = initial_node.priority(self.consensus_cost());

        // start the priority queue, which defaults to bigger is better so we need a Reverse since want smaller costs
        let mut pqueue: PriorityQueue<DualConsensusNode, NodePriority> = PriorityQueue::new();
        single_tracker.insert(initial_node.max_consensus_length());
        pqueue.push(initial_node, initial_priority);

        let mut ret: Vec<DualConsensus> = vec![];
        let mut last_indiv = 0;

        // we need to do min_count based on the number of *active* sequences when we have really high coverage
        // dynamic min count which is the max of either the raw minimum OR the number of sequences * min_af
        let full_min_count = self.config.min_count.max(
            (self.config.min_af * self.sequences.len() as f64).ceil() as u64
        );
        let mut total_active_count = vec![
            initially_active
        ];
        let mut active_min_count = vec![
            self.config.min_count.max(
                (self.config.min_af * initially_active as f64).ceil() as u64
            )
        ];

        // the way this will work is that we will eventually find one or more answers and anything worse will get drained off until no possibilities remain
        let mut iteration = 0;
        while !pqueue.is_empty() {
            // this just tracks how large our actual queue gets
            peak_queue_size = peak_queue_size.max(pqueue.len());

            // first, check if we need to restrict our pqueue threshold
            while (single_tracker.len() > max_queue_size || single_last_constraint >= self.config.max_nodes_wo_constraint) && single_tracker.threshold() < farthest_consensus_single {
                single_tracker.increment_threshold();
                single_last_constraint = 0;
            }
            while (dual_tracker.len() > max_queue_size || dual_last_constraint >= self.config.max_nodes_wo_constraint) && dual_tracker.threshold() < farthest_consensus_dual {
                dual_tracker.increment_threshold();
                dual_last_constraint = 0;
            }

            // get the top node and check if it's bad
            let (top_node, top_cost) = pqueue.pop().unwrap();
            let top_len = top_node.max_consensus_length();

            let (threshold_cutoff, at_capacity) = if top_node.is_dual {
                trace!(
                    "{} -> \"{}\" + \"{}\", {:?}", 
                    top_cost.0.0,
                    top_node.consensus1.len(), 
                    top_node.consensus2.len(),
                    top_node.costs(self.config.consensus_cost)
                );
                dual_tracker.remove(top_len);
                (dual_tracker.threshold(), dual_tracker.at_capacity(top_len))
            } else {
                trace!("{} -> \"{}\"", top_cost.0.0, top_node.consensus1.len());
                single_tracker.remove(top_len);
                (single_tracker.threshold(), single_tracker.at_capacity(top_len))
            };
        
            // check if this node is worse than the best full solution OR
            if top_cost.0.0 > maximum_error || 
                // if it's shorter than the pqueue threshold OR
                top_len < threshold_cutoff ||
                // if we have reached process capacity OR
                at_capacity ||
                // if it's a dual node that is no longer containing enough on each group
                top_node.is_dual_imbalanced(active_min_count[top_len] as usize) {
                nodes_ignored += 1;
                trace!("\tignored {} || {} || {}", top_cost.0.0 > maximum_error, top_len < threshold_cutoff, top_node.is_dual_imbalanced(active_min_count[top_node.max_consensus_length()] as usize));
                continue;
            }

            // mark this one as explored and updated farthest if necessary
            if top_node.is_dual {
                farthest_consensus_dual = farthest_consensus_dual.max(top_len);
                dual_last_constraint += 1;
                dual_tracker.process(top_len)?;
            } else {
                farthest_consensus_single = farthest_consensus_single.max(top_len);
                single_last_constraint += 1;
                single_tracker.process(top_len)?;
            }
            nodes_explored += 1;

            if !top_node.is_dual {
                last_indiv = last_indiv.max(top_len);
            }

            if iteration % 1000 == 0 {
                trace!("i: {}, t: {}, s: {}, d: {}", iteration, pqueue.len(), single_tracker.len(), dual_tracker.len());
                trace!(
                    "Handling: cost={}, dual={}, h1_len={}, h2_len={}; s_t: {}, d_t: {}",
                    top_node.total_cost(self.config.consensus_cost),
                    top_node.is_dual,
                    top_node.consensus1.len(),
                    top_node.consensus2.len(),
                    single_tracker.threshold(),
                    dual_tracker.threshold()
                );
            }

            let l = 50.min(top_node.consensus1.len());
            if l > 0 {
                trace!("\t{}..{}", 
                    std::str::from_utf8(&top_node.consensus1[..l]).unwrap(), 
                    std::str::from_utf8(&top_node.consensus1[(top_node.consensus1.len()-l)..]).unwrap());
            }
            if top_node.is_dual {
                let l = 50.min(top_node.consensus2.len());
                if l > 0 {
                    trace!("\t{}..{}", 
                        std::str::from_utf8(&top_node.consensus2[..l]).unwrap(),
                        std::str::from_utf8(&top_node.consensus2[(top_node.consensus2.len()-l)..]).unwrap());
                }
            }
            iteration += 1;
            /*
            if iteration == 40000 {
                panic!("max");
            }
            */

            // now check if this node has reached the end
            if top_node.reached_all_end(&self.sequences, self.config.allow_early_termination) {
                // this node *might* be done, we need to finalize it to be sure
                // we also need to do it in a clone since it might not be finalized and we need to keep extending
                let mut finalized_node = top_node.clone();
                finalized_node.finalize(&self.sequences)?;

                // check if this is a dual node, but without enough support on each allele
                let imbalanced = if finalized_node.is_dual {
                    let (best_consensus_index, _best_consensus_score) = finalized_node.costs(self.consensus_cost());
                    let counts1 = best_consensus_index.iter()
                        .filter(|&&v| v == 0)
                        .count();
                    let counts2 = best_consensus_index.len() - counts1;

                    // if either does not have enough support, we will ignore it
                    (counts1 as u64) < full_min_count || (counts2 as u64) < full_min_count
                } else {
                    // non-dual cannot be imbalanced
                    false
                };

                if !imbalanced {
                    // get the new finalized score, it may have changed
                    let finalized_score = finalized_node.total_cost(self.config.consensus_cost);
                    
                    // check first if it's strictly BETTER than anything so far
                    if finalized_score < maximum_error {
                        // this score is better than anything we have seen so far, clear out previous results if we have any
                        maximum_error = finalized_score;
                        trace!("\tMaximum error set to {maximum_error}");
                        ret.clear();
                    }

                    // now check if it's as good as anything so far; the above check will not nullify this one
                    if finalized_score <= maximum_error && ret.len() < self.config.max_return_size {
                        // this consensus is at least as good as the best, so add it as a result
                        let dual_con_result = DualConsensus::from_node(
                            &finalized_node,
                            self.config.consensus_cost
                        );
                        trace!("\tadding to ret");//: {dual_con_result:?}");
                        trace!("\tcon1: {}", std::str::from_utf8(dual_con_result.consensus1().sequence())?);
                        trace!("\tcon2: {}", std::str::from_utf8(dual_con_result.consensus2().unwrap().sequence())?);
                        ret.push(dual_con_result);
                    }
                } else {
                    // this node is not balanced correctly for us to return it as a solution
                    // note that it may still be a candidate for extension though
                    trace!("Finalized node is imbalanced, ignoring.");
                }
            }

            // everything past this point is about extending the top_node further

            // check if we need to update the active
            if active_min_count.len() == top_len+1 {
                // we need to copy and update
                let current_active = total_active_count[top_len];
                let new_additions = match activate_points.get(&top_len) {
                    Some(v) => v.len(),
                    None => 0
                };
                let new_total = current_active + new_additions;
                total_active_count.push(new_total);
                // debug!("total_active_count[{}] = {}", total_active_count.len() - 1, new_total);
                
                let new_min_af = self.config.min_count.max(
                    (self.config.min_af * new_total as f64).ceil() as u64
                );
                active_min_count.push(new_min_af);
                // debug!("active_min_count[{}] = {}", total_active_count.len() - 1, new_min_af);
            }

            // this fetches the list of options according to the WFA so far
            // NOTE: this CAN include the wildcard, but only if the wildcard is the only character
            let weighted_by_ed = self.config.weighted_by_ed;
            let extension_candidates1 = top_node.get_extension_candidates(&self.sequences, self.config.wildcard, true, weighted_by_ed);
            // let min_count1 = active_min_count[top_len];
            let min_count1 = self.config.min_count.max(
                (self.config.min_af * extension_candidates1.values().sum::<f64>()).ceil() as u64
            );
            let max_observed1 = extension_candidates1.values().cloned().max_by(|a, b| a.total_cmp(b))
                // if no observation, then just use min count
                .unwrap_or(min_count1 as f64);
            
            // the active threshold is the minimum of 1) the configured minimum count OR 2) the highest count we observed
            let active_threshold1 = (min_count1 as f64).min(max_observed1);

            if top_node.is_dual {
                // get the second candidate set also
                let extension_candidates2 = top_node.get_extension_candidates(&self.sequences, self.config.wildcard, false, weighted_by_ed);
                // let min_count2 = active_min_count[top_len];
                let min_count2 = self.config.min_count.max(
                    (self.config.min_af * extension_candidates2.values().sum::<f64>()).ceil() as u64
                );
                let max_observed2 = extension_candidates2.values().cloned().max_by(|a, b| a.total_cmp(b))
                    // if no observation, then just use min count
                    .unwrap_or(min_count2 as f64);
            
                // the active threshold is the minimum of 1) the configured minimum count OR 2) the highest count we observed
                let active_threshold2 = (min_count2 as f64).min(max_observed2);

                // one thing we have to handle here is the situation where alleles are different length
                // this means one might want more extensions and the other is all done
                // so lets check if either allele is ready to be done
                let is_con1_finalized = top_node.reached_consensus_end(&self.sequences, true, self.config.allow_early_termination);
                let is_con2_finalized = top_node.reached_consensus_end(&self.sequences, false, self.config.allow_early_termination);

                trace!("\tec1: {extension_candidates1:?}, {active_threshold1}");
                trace!("\tec2: {extension_candidates2:?}, {active_threshold2}");

                // now create adjust extension lists that can have None
                let opt_ec1: Vec<Option<u8>> = {
                    let mut v = vec![];
                    // if consensus 1 has at least one final candidate OR the candidate list is empty OR it's already locked
                    // THEN we need the no-extend option None
                    if is_con1_finalized || extension_candidates1.is_empty() || top_node.is_con1_locked() {
                        v.push(None);
                    }

                    // do not add extensions if we are sequence locked
                    if !top_node.is_con1_locked() {
                        v.extend(extension_candidates1.iter()
                            .filter_map(|(&symbol, &count)| {
                                // if this character occurs fewer than the threshold times, remove it from the candidate set
                                if count < active_threshold1 {
                                    None
                                } else {
                                    // double Some is because of filter_map
                                    Some(Some(symbol))
                                }
                            })
                        );
                    }
                    v
                };

                let opt_ec2: Vec<Option<u8>> = {
                    let mut v = vec![];
                    // if consensus 2 has at least one final candidate OR the candidate list is empty OR it's already locked
                    // then we need the no-extend option None
                    if is_con2_finalized || extension_candidates2.is_empty() || top_node.is_con2_locked() {
                        v.push(None);
                    }

                    // do not add extensions if we are sequence locked
                    if !top_node.is_con2_locked() {
                        v.extend(extension_candidates2.iter().filter_map(|(&symbol, &count)| {
                            // if this character occurs fewer than the threshold times, remove it from the candidate set
                            if count < active_threshold2 {
                                None
                            } else {
                                // double Some is because of filter_map
                                Some(Some(symbol))
                            }
                        })
                    );
                    }
                    v
                };

                // make sure we always have something in each
                assert!(!opt_ec1.is_empty() && !opt_ec2.is_empty());

                for opt_can1 in opt_ec1.iter() {
                    for opt_can2 in opt_ec2.iter() {
                        if opt_can1.is_none() && opt_can2.is_none() {
                            // we have to extend at least one of them or else it's the exact same node
                            continue;
                        }

                        trace!("\tExtending with: {opt_can1:?} + {opt_can2:?}");

                        // pair wise add each extension; above check enforces that at least one of these `if` statements are executed
                        let mut new_node = top_node.clone();
                        if let Some(c1) = opt_can1 {
                            // we have a base to add
                            new_node.push(&self.sequences, *c1, true)?;
                        } else {
                            // we don't, lock the sequence in place so we do not get duplicate nodes
                            new_node.lock_sequence(true);
                        }
                        if let Some(c2) = opt_can2 {
                            // we have a base to add
                            new_node.push(&self.sequences, *c2, false)?;
                        } else {
                            // we don't, lock the sequence in place so we do not get duplicate nodes
                            new_node.lock_sequence(false);
                        }

                        // check if we need to activate any strings
                        let opt_activate_list = activate_points.get(&new_node.max_consensus_length());
                        if let Some(activate_list) = opt_activate_list {
                            assert!(!activate_list.is_empty());
                            for &seq_index in activate_list.iter() {
                                new_node.activate_sequence(self.sequences[seq_index], seq_index, offset_window, offset_compare_length, self.config.wildcard, self.config.allow_early_termination)?;
                            }
                        }

                        // check if we can prune anything after extensions
                        new_node.prune_dwfa(self.config.dual_max_ed_delta)?;

                        // get the new cost and put it in the queue
                        let new_priority = new_node.priority(self.consensus_cost());
                        assert!(new_node.is_dual);
                        dual_tracker.insert(new_node.max_consensus_length()); // top_node is already dual
                        assert!(pqueue.push(new_node.clone(), new_priority).is_none());
                    }
                }
            } else {
                // this is currently a non-dual node
                // first, handle the option where it stays as non-dual
                trace!("\tec1: {extension_candidates1:?}, {active_threshold1}");
                for (&symbol, &count) in extension_candidates1.iter() {
                    if count < active_threshold1 {
                        // this symbol was not observed enough times, skip it to reduce overhead
                        continue;
                    }

                    // clone and add the symbol
                    let mut new_node = top_node.clone();
                    new_node.push(&self.sequences, symbol, true)?;

                    // check if we need to activate any strings
                    let opt_activate_list = activate_points.get(&new_node.max_consensus_length());
                    if let Some(activate_list) = opt_activate_list {
                        assert!(!activate_list.is_empty());
                        for &seq_index in activate_list.iter() {
                            new_node.activate_sequence(self.sequences[seq_index], seq_index, offset_window, offset_compare_length, self.config.wildcard, self.config.allow_early_termination)?;
                        }
                    }

                    // get the new cost and put it in the queue
                    let new_priority = new_node.priority(self.consensus_cost());
                    assert!(!new_node.is_dual);
                    single_tracker.insert(new_node.max_consensus_length());
                    assert!(pqueue.push(new_node.clone(), new_priority).is_none());
                }

                // now handle dual-node generation
                let mut num_passing = 0;
                let sorted_candidates = {
                    let mut sc: Vec<_> = extension_candidates1.iter()
                        .filter_map(|(&c, &count)| {
                            // when it comes to dual splitting nodes, we care about hitting the minimum count
                            if self.config.wildcard.is_some() && c == self.config.wildcard.unwrap() {
                                // this is the wildcard, we do not count it for the purpose of dual splitting
                                None
                            } else {
                                if count >= min_count1 as f64 {
                                    num_passing += 1;
                                }
                                Some((Reverse(count), c))
                            }
                        })
                        .collect();
                    sc.sort_by(|a, b| a.partial_cmp(b).unwrap());
                    sc
                };

                if num_passing > 1 {
                    // c1 in outer loop
                    for (i, &(_order, c1)) in sorted_candidates.iter().enumerate() {
                        // c2 in inner loop, must start after c1
                        for &(_order2, c2) in sorted_candidates[(i+1)..].iter() {
                            // sanity checks to make sure we never try to make a dual node where one extension is the wildcard
                            assert!(self.config.wildcard.is_none() || c1 != self.config.wildcard.unwrap());
                            assert!(self.config.wildcard.is_none() || c2 != self.config.wildcard.unwrap());
                            
                            // convert the clone into a dual node with the candidate; c1 should be the "major" candidate due to the above sort
                            let mut new_node = top_node.clone();
                            new_node.activate_dual(&self.sequences, c1, c2)?;

                            // check if we need to activate any strings
                            let opt_activate_list = activate_points.get(&new_node.max_consensus_length());
                            if let Some(activate_list) = opt_activate_list {
                                assert!(!activate_list.is_empty());
                                for &seq_index in activate_list.iter() {
                                    new_node.activate_sequence(self.sequences[seq_index], seq_index, offset_window, offset_compare_length, self.config.wildcard, self.config.allow_early_termination)?;
                                }
                            }

                            // now prune
                            new_node.prune_dwfa(self.config.dual_max_ed_delta)?;

                            // get the new cost and put it in the queue
                            let new_priority = new_node.priority(self.consensus_cost());
                            assert!(new_node.is_dual);
                            dual_tracker.insert(new_node.max_consensus_length());
                            assert!(pqueue.push(new_node.clone(), new_priority).is_none());
                        }
                    }
                }
            }

            // TODO: remove this eventually, it's a summation check
            assert_eq!(pqueue.len(), single_tracker.unfiltered_len() + dual_tracker.unfiltered_len());
        }

        assert_eq!(single_tracker.len(), 0);
        assert_eq!(dual_tracker.len(), 0);

        if ret.len() > 1 {
            // sort these by the sequence so we always have a fixed order
            ret.sort_by(|c1, c2| {
                let empty = vec![];
                
                // first consensus key
                let s1 = c1.consensus1.sequence();
                let s2 = match c1.consensus2.as_ref() {
                    Some(s) => s.sequence(),
                    None => &empty
                };
                let k1 = (s1, s2);

                // second consensus key
                let s1 = c2.consensus1.sequence();
                let s2 = match c2.consensus2.as_ref() {
                    Some(s) => s.sequence(),
                    None => &empty
                };
                let k2 = (s1, s2);
                k1.cmp(&k2)
            });
        }

        // TODO: we can hit this if we remove all final solution because they are imbalanced
        // if we find that happens in practice, then we need to somehow distinguish dual and singleton nodes to prevent the scenario where
        //     we remove all non-dual nodes, but then all the dual nodes fail the final check
        // assert!(!ret.is_empty());
        if ret.is_empty() {
            warn!("No consensus found that reached end, is there a gap between input sequences?");

            // TODO: how do we want to handle this long-term? this returns an empty string consensus
            let no_offsets = vec![None; self.sequences.len()]; // we need these to get costs of 0
            let root_node = DualConsensusNode::new_root_node(&no_offsets, self.config.wildcard, self.config.allow_early_termination)?;
            ret.push(DualConsensus::from_node(&root_node, self.consensus_cost()));
        }

        debug!("nodes_explored: {nodes_explored}");
        debug!("nodes_ignored: {nodes_ignored}");
        debug!("peak_queue_size: {peak_queue_size}");
        debug!("last_indiv: {last_indiv}");

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
struct DualConsensusNode {
    /// if True, then this node is tracking two consensuses
    is_dual: bool,
    /// if True, we are not allowed to modify consensus1 anymore
    is_con1_locked: bool,
    /// if True, we are not allowed to modify consensus2 anymore
    is_con2_locked: bool,
    /// The primary consensus sequence
    consensus1: Vec<u8>,
    /// The secondary consensus sequence
    consensus2: Vec<u8>,
    /// The set of DWFAs for consensus1; these are options because we stop tracking once the scores diverge
    dwfas1: Vec<Option<DWFALite>>,
    /// The set of DWFAs for consensus2; these are options because we stop tracking once the scores diverge
    dwfas2: Vec<Option<DWFALite>>,
}

impl DualConsensusNode {
    /// Constructor for a new consensus search root node.
    /// Note that initial it is not a dual node, it will become that when divergence is detected.
    /// # Arguments
    /// * `offsets` - a set of offsets into the sequences where the approximate starts are
    /// * `wildcard` - an optional wildcard symbol that will match anything
    /// * `allow_early_termination` - if true, then it will allow the consensus to go beyond the provided baseline sequences without penalty
    /// # Errors
    /// * if DWFA construction fails
    fn new_root_node(offsets: &[Option<usize>], wildcard: Option<u8>, allow_early_termination: bool) -> Result<DualConsensusNode, Box<dyn std::error::Error>> {
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
            bail!("Root DualConsensus node has no initially active sequences.");
        }

        // a root node will have one empty consensus, and one corresponding vec of empty DWFAs
        Ok(DualConsensusNode {
            is_dual: false,
            is_con1_locked: false,
            is_con2_locked: false,
            consensus1: vec![],
            consensus2: vec![],
            dwfas1: dwfas,
            dwfas2: vec![None; offsets.len()]
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
        // figure out whether we need to just do con1 or both
        let activators = if self.is_dual {
            vec![(&mut self.dwfas1, &self.consensus1), (&mut self.dwfas2, &self.consensus2)]
        } else {
            vec![(&mut self.dwfas1, &self.consensus1)]
        };

        for (dwfas, consensus) in activators.into_iter() {
            // make sure everything is currently inactive
            assert!(dwfas[seq_index].is_none());
    
            let con_len = consensus.len();
            let start_delta = offset_window + offset_compare_length;
            let start_position = con_len.saturating_sub(start_delta);
            let end_position = con_len.saturating_sub(offset_compare_length);

            // figure out which offset has the best score; assume the middle of the offset window is the best
            let mut best_offset = con_len.saturating_sub(offset_compare_length + offset_window / 2);
            let mut min_ed = crate::sequence_alignment::wfa_ed_config(&consensus[best_offset..], &sequence[0..offset_compare_length], false, wildcard);
            
            for p in start_position..end_position {
                let ed = crate::sequence_alignment::wfa_ed_config(&consensus[p..], &sequence[0..offset_compare_length], false, wildcard);
                if ed < min_ed {
                    min_ed = ed;
                    best_offset = p;
                }
            }

            // now set up the DWFA with the best offset
            let mut new_dwfa = DWFALite::new(wildcard, allow_early_termination);
            new_dwfa.set_offset(best_offset);
            new_dwfa.update(sequence, consensus)?;
            dwfas[seq_index] = Some(new_dwfa);
        }

        Ok(())
    }

    /// Adds a new symbol to the DWFAs
    /// # Arguments
    /// * `baseline_sequences` - the sequences that are getting compared against in the DWFAs
    /// * `symbol` - the symbol to extend each DWFA with
    /// * `is_consensus1` - if True, this will check consensus1 DWFAs, otherwise it checks those for consensus2
    /// # Errors
    /// * if any of the extensions fail
    fn push(&mut self, baseline_sequences: &[&[u8]], symbol: u8, is_consensus1: bool) -> Result<(), Box<dyn std::error::Error>> {
        if is_consensus1 && self.is_con1_locked {
            bail!("Consensus 1 is locked, cannot modify");
        } else if  !is_consensus1 && self.is_con2_locked {
            bail!("Consensus 2 is locked, cannot modify");
        }

        // get the relevant DWFAs
        let dwfa_iter = if is_consensus1 {
            self.dwfas1.iter_mut()
        } else {
            self.dwfas2.iter_mut()
        };

        // get the relevant consensus
        let consensus_seq = if is_consensus1 {
            &mut self.consensus1
        } else {
            &mut self.consensus2
        };

        // first extend the consensus
        consensus_seq.push(symbol);

        // now adjust the corresponding DWFA
        for (&baseline, opt_dwfa) in baseline_sequences.iter().zip(dwfa_iter) {
            if let Some(dwfa) = opt_dwfa {
                dwfa.update(baseline, consensus_seq)?;
            }
        }
        Ok(())
    }

    /// Converts this node into a dual node and extends with the new different symbols.
    /// # Arguments
    /// * `baseline_sequences` - the sequences that are getting compared against in the DWFAs
    /// * `symbol1` - the symbol to extend the first consensus with
    /// * `symbol2` - the symbol to extend the second consensus with
    /// # Errors
    /// * if this is already a dual node
    /// * if symbol1 == symbol2, they must be different or it's not a dual node
    /// * if any of the extensions fail
    fn activate_dual(&mut self, baseline_sequences: &[&[u8]], symbol1: u8, symbol2: u8) -> Result<(), Box<dyn std::error::Error>> {
        if self.is_dual {
            bail!("Cannot activate dual on a dual node");
        }
        self.is_dual = true;

        if symbol1 == symbol2 {
            bail!("Cannot activate dual mode with the same extension symbols");
        }

        // clone everything in 1 into 2
        self.consensus2 = self.consensus1.clone();
        self.dwfas2 = self.dwfas1.clone();

        // now extend both consensus with the new symbols
        self.push(baseline_sequences, symbol1, true)?;
        self.push(baseline_sequences, symbol2, false)?;

        Ok(())
    }

    /// Returns true if this is a dual node AND one of the two tracked has fewer counts than the minimum
    /// # Arguments
    /// * `min_count` - the minimum count required to have a dual allele result
    fn is_dual_imbalanced(&self, min_count: usize) -> bool {
        if self.is_dual {
            let counts1 = self.dwfas1.iter()
                .filter_map(|d| d.as_ref())
                .count();
            let counts2 = self.dwfas2.iter()
                .filter_map(|d| d.as_ref())
                .count();

            // return true if either is less than the minimum
            counts1 < min_count || counts2 < min_count

            /*
            let (c1, c2) = self.full_cost(ConsensusCost::L1Distance);
            let mut counts1 = 0.0;
            let mut counts2 = 0.0;
            let equality_score = 0.5;
            for (opt_v1, opt_v2) in c1.into_iter().zip(c2.into_iter()) {
                if let (Some(v1), Some(v2)) = (opt_v1, opt_v2) {
                    // both are present, compare
                    match v1.cmp(&v2) {
                        std::cmp::Ordering::Less => counts1 += 1.0,
                        std::cmp::Ordering::Equal => {
                            counts1 += equality_score;
                            counts2 += equality_score;
                        },
                        std::cmp::Ordering::Greater => counts2 += 1.0,
                    };
                } else if opt_v1.is_some() {
                    counts1 += 1.0;
                } else if opt_v2.is_some() {
                    counts2 += 1.0;
                } else {
                    // both are None, I think leave this empty for now
                }
            }
            
            // return true if either is less than the minimum
            counts1 < (min_count as f64) || counts2 < (min_count as f64)
            */
        } else {
            false
        }
    }

    /// This will go through the pairs of DWFA and remove any that have diverge from equality.
    /// This is a no-op if the node is not dual yet.
    /// # Arguments
    /// * `ed_delta` - the maximum difference between the edit distances before the higher one gets dropped from tracking
    fn prune_dwfa(&mut self, ed_delta: usize) -> Result<(), Box<dyn std::error::Error>> {
        if self.is_dual {
            for (opt_d1, opt_d2) in self.dwfas1.iter_mut().zip(self.dwfas2.iter_mut()) {
                if let (Some(d1), Some(d2)) = (opt_d1.as_ref(), opt_d2.as_ref()) {
                    if d1.edit_distance() + ed_delta < d2.edit_distance() {
                        // d1 is much less than d2; drop d2
                        *opt_d2 = None;
                    } else if d2.edit_distance() + ed_delta < d1.edit_distance() {
                        // d2 is much less than d1; drop d1
                        *opt_d1 = None;
                    }
                }
            }
        }
        Ok(())
    }

    /// Adds an internal lock that prevents any further additions to the sequence
    pub fn lock_sequence(&mut self, is_consensus1: bool) {
        if is_consensus1 {
            self.is_con1_locked = true;
        } else {
            self.is_con2_locked = true;
        }
    }

    /// Finalizes all of the contained DWFAs
    /// # Errors
    /// * if any of the individual finalizes fail
    fn finalize(&mut self, baseline_sequences: &[&[u8]]) -> Result<(), Box<dyn std::error::Error>> {
        // we will always loop over the first
        for (&baseline, (opt_dwfa1, opt_dwfa2)) in baseline_sequences.iter()
            .zip(self.dwfas1.iter_mut().zip(self.dwfas2.iter_mut())) {
            let mut seq_finalized = false;
            if let Some(dwfa) = opt_dwfa1 {
                dwfa.finalize(baseline, &self.consensus1)?;
                seq_finalized = true;
            }

            // if this is a dual node, loop over those also
            if self.is_dual {
                if let Some(dwfa) = opt_dwfa2 {
                    dwfa.finalize(baseline, &self.consensus2)?;
                    seq_finalized = true;
                }
            }

            if !seq_finalized {
                bail!("Finalize called on DWFA that was never initialized.");
            }
        }

        // if we are finalizing, we are definitely locking
        self.is_con1_locked = true;
        self.is_con2_locked = true;
        
        Ok(())
    }

    /// Returns a tuple of scores for this node of form (indices, scores).
    /// It will pick the minimum cost for any sequence where multiple possible consensuses are still being tracked.
    /// The minimum cost consensus index will go into `indices` and the actual minimum score into `scores`.
    fn costs(&self, consensus_cost: ConsensusCost) -> (Vec<usize>, Vec<usize>) {
        // we need to collect the minimum cost across multiple DWFAs for a given sequence
        let num_sequences = self.dwfas1.len();
        let mut best_consensus_index = vec![usize::MAX; num_sequences];
        let mut best_consensus_score = vec![usize::MAX; num_sequences];

        // iterate over the two DWFA collections
        let mut dwfa_iters = [self.dwfas1.iter(), self.dwfas2.iter()];
        for (con_index, dwfa_iter) in dwfa_iters.iter_mut().enumerate() {
            for (seq_index, opt_d) in dwfa_iter.enumerate() {
                if let Some(d) = opt_d {
                    let score = match consensus_cost {
                        ConsensusCost::L1Distance => d.edit_distance(),
                        ConsensusCost::L2Distance => d.edit_distance().pow(2)
                    };

                    // if this score is BETTER, then replace the previous
                    if score < best_consensus_score[seq_index] {
                        best_consensus_score[seq_index] = score;
                        best_consensus_index[seq_index] = con_index;
                    }
                } else {
                    // we can ignore these, they do not have a DWFA for one reason or another
                    // the assertion below verifies that at least one per sequence is present though
                }
            }
        }

        // make sure all of them match 0 or 1
        // assert!(best_consensus_index.iter().all(|&i| i == 0 || i == 1));
        
        // the above no longer holds when a sequence has not been activated yet, replace all of the usize::MAX with 0
        for (bci, bcs) in best_consensus_index.iter().zip(best_consensus_score.iter_mut()) {
            if *bci == usize::MAX && *bcs == usize::MAX {
                *bcs = 0;
            }
        }

        (best_consensus_index, best_consensus_score)
    }

    /// Returns the total score for the node
    fn total_cost(&self, consensus_cost: ConsensusCost) -> usize {
        let (_best_indices, best_costs) = self.costs(consensus_cost);
        best_costs.iter().sum()
    }

    /// Returns the full set of tracked costs for the two consensuses.
    fn full_cost(&self, consensus_cost: ConsensusCost) -> (Vec<Option<usize>>, Vec<Option<usize>>) {
        let s1 = self.dwfas1.iter()
            .map(|opt_dwfa| {
                opt_dwfa.as_ref().map(|d| match consensus_cost {
                    ConsensusCost::L1Distance => d.edit_distance(),
                    ConsensusCost::L2Distance => d.edit_distance().pow(2)
                })
            }).collect();
        let s2 = self.dwfas2.iter()
            .map(|opt_dwfa| {
                opt_dwfa.as_ref().map(|d| match consensus_cost {
                    ConsensusCost::L1Distance => d.edit_distance(),
                    ConsensusCost::L2Distance => d.edit_distance().pow(2)
                })
            }).collect();
        
        (s1, s2)
    }

    /// Returns the node priority.
    /// Currently, this is based on 1) lowest cost and 2) consensus length.
    /// # Arguments
    /// * `consensus_cost` - cost model to evaluate the cost
    fn priority(&self, consensus_cost: ConsensusCost) -> NodePriority {
        (
            Reverse(self.total_cost(consensus_cost)),
            self.max_consensus_length()
        )
    }

    /// Returns true if each consensus has at least one DWFA at the end of their respective baseline sequences.
    /// # Arguments
    /// * `baseline_sequences` - the sequences that are fixed that we want the consensus of
    /// * `require_all` - if True, requires *all* of the sequences to be at their end in at least one of the consensuses (this is useful for allow_early_termination)
    fn reached_all_end(&self, baseline_sequences: &[&[u8]], require_all: bool) -> bool {
        // create an iterator that determines if we're at the end
        let mut iter_map = baseline_sequences.iter()
            .zip(self.dwfas1.iter().zip(self.dwfas2.iter()))
            .map(|(&baseline, (opt_dwfa1, opt_dwfa2))| {
                // assert!(opt_dwfa1.is_some() || opt_dwfa2.is_some());
                // these can be None if either A) this one has not started or B) it has started, but dropped off due to high ED
                // in either case, it would default to NOT at end (i.e., false)
                let p1 = opt_dwfa1.as_ref().map(|d| d.reached_baseline_end(baseline));
                let p2 = opt_dwfa2.as_ref().map(|d| d.reached_baseline_end(baseline));

                // at least one of them needs to be at the end to pass; untracked does not count
                p1.unwrap_or(false) || p2.unwrap_or(false)
            });
        
        // handle iterator appropriately
        if require_all {
            iter_map.all(|b| b)
        } else {
            iter_map.any(|b| b)
        }
    }

    /// For the given consensus, returns true if at least one *active* DWFA at the end of the baseline sequences.
    /// If this is not a dual node, then it will always return false for consensus 2.
    /// # Arguments
    /// * `baseline_sequences` - the sequences that are fixed that we want the consensus of
    /// * `is_consensus1` - if True, this will check consensus1 DWFAs, otherwise it checks those for consensus2
    /// * `require_all` - if True, this requires all active DWFAs to be at the end
    fn reached_consensus_end(&self, baseline_sequences: &[&[u8]], is_consensus1: bool, require_all: bool) -> bool {
        // TODO: do we need to refactor this to look at both consensuses at once?
        let dwfa_iter = if is_consensus1 {
            self.dwfas1.iter()
        } else {
            if !self.is_dual {
                // this is not a dual node, but the user is asking for consensus 2
                return false;
            }

            self.dwfas2.iter()
        };
        
        // if we require all, then we do not want to penalize; so it should be True
        // if not, then we don't want to count it; so it should be False
        let value_for_inactive = require_all;

        // iterate over each DWFA and check if it's at the end
        let mut iter_map = baseline_sequences.iter().zip(dwfa_iter)
            .map(|(&baseline, opt_dwfa)| {
                opt_dwfa.as_ref()
                    .map(|d| d.reached_baseline_end(baseline))
                    .unwrap_or(value_for_inactive)
            });
        
        if require_all {
            iter_map.all(|b| b)
        } else {
            iter_map.any(|b| b)
        }
    }

    /// Returns a hashmap of the extension candidates with their votes.
    /// The votes are fractional if a sequence has multiple equally likely options.
    /// # Arguments
    /// * `baseline_sequences` - the sequences that are fixed that we want the consensus of
    /// * `wildcard` - an optional wildcard character, will be removed from return set unless it is the only value in it
    /// * `is_consensus1` - if True, this will check consensus1 DWFAs, otherwise it checks those for consensus2
    /// * `weighted_by_ed` - if True, then the weights are scaled based on the dual weights comparison; e.g. if ed1 = 2 and ed2 = 4, then weights for consensus 1 are ~1/3 and for consensus 2 are ~2/3 of the total
    fn get_extension_candidates(&self, baseline_sequences: &[&[u8]], wildcard: Option<u8>, is_consensus1: bool, weighted_by_ed: bool) -> HashMap<u8, f64> {
        // get the relevant DWFAs
        let dwfa_iter = if is_consensus1 {
            self.dwfas1.iter()
        } else {
            self.dwfas2.iter()
        };

        // get the relevant consensus
        let consensus_seq = if is_consensus1 {
            &self.consensus1
        } else {
            &self.consensus2
        };

        let weights = if weighted_by_ed {
            self.get_ed_weights(is_consensus1, weighted_by_ed)
        } else {
            // these weights will always factor in current scoring
            vec![1.0; self.dwfas1.len()]
        };

        // now pull out all the candidates and count how many times they were voted on
        let mut candidates: HashMap<u8, f64> = Default::default();
        for ((&baseline_seq, opt_dwfa), &weight) in baseline_sequences.iter().zip(dwfa_iter).zip(weights.iter()) {
            if weight > 0.0 {
                if let Some(dwfa) = opt_dwfa {
                    // get the candidates and the total observation weight
                    let cand = dwfa.get_extension_candidates(baseline_seq, consensus_seq);
                    let vote_split = cand.values().sum::<usize>() as f64;
                    
                    // iterate over each candidate and scale it by the occurrences count / total weight
                    for (&c, &occ) in cand.iter() {
                        let entry = candidates.entry(c).or_insert(0.0);
                        *entry += weight * occ as f64 / vote_split;
                    }
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

    /// This returns a series of weights that are scaled based on the edit distances to con1 and con2.
    /// For example, if ed1 = 2 and ed2 = 5, then a consensus1 scaled ed would be ~ 5 / (2+5).
    /// Conversely, the weight for consensus2 would be ~ 2 / (2+5).
    /// If `weight_by_ed` is False, this will instead just return 0.0, 0.5, or 1.0; giving full weight when one consensus is better.
    /// # Arguments
    /// * `is_consensus` - If True, then return weights for consensus 1, otherwise consensus 2.
    /// * `weight_by_ed` - If True, then weights are scaled by the relative edit distances. If False, then weights are either 0.0 or 1.0 typically, with 0.5 when equal.
    fn get_ed_weights(&self, is_consensus1: bool, weight_by_ed: bool) -> Vec<f64> {
        if self.is_dual {
            let min_ed = 0.5; // prevent divide by zero while still enabling a scaling
            // TODO: should we treat these as 1.0, 0.5, or 0.0?
            let equality_score = 0.5; // using 0.0 means don't let something that isn't for sure mapped here have a vote
            self.dwfas1.iter().zip(self.dwfas2.iter())
                .map(|(d1, d2)| {
                    let c1 = d1.as_ref().map(|d| (d.edit_distance() as f64).max(min_ed));
                    let c2 = d2.as_ref().map(|d| (d.edit_distance() as f64).max(min_ed));
                    
                    if let (Some(v1), Some(v2)) = (c1, c2) {
                        if weight_by_ed {
                            // numerator is the opposite entry; i.e for con1, if ed1 = 1 and ed2 = 9, we should get 90% of weight
                            let numer = if is_consensus1 { v2 } else { v1 };
                            numer / (v1 + v2)
                        } else if v1 == v2 {
                            // equal weights
                            equality_score
                        } else if (is_consensus1 && v1 < v2) || (!is_consensus1 && v2 < v1) {
                            // the consensus we are looking at is less than the other, so give it full weight
                            1.0
                        } else {
                            // the consensus we are looking at is more than the other, so it gets no weight
                            0.0
                        }
                    } else if (c1.is_some() && is_consensus1) || (c2.is_some() && !is_consensus1) {
                        // the one we are checking is present, but the other is not
                        1.0
                    } else {
                        // everything else gets no weight
                        0.0
                    }
                }).collect()
        } else {
            // this one isn't dual currently, just return all full weight
            vec![1.0; self.dwfas1.len()]
        }
    }

    // getters
    fn max_consensus_length(&self) -> usize {
        self.consensus1.len().max(self.consensus2.len())
    }

    fn is_con1_locked(&self) -> bool {
        self.is_con1_locked
    }

    fn is_con2_locked(&self) -> bool {
        self.is_con2_locked
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::path::PathBuf;

    use crate::cdwfa_config::CdwfaConfigBuilder;

    // first some more targeted tests
    #[test]
    fn test_get_ed_weights() {
        let vec_sequences = vec![
            b"ACGT".to_vec(),
            b"CGTA".to_vec()
        ];
        let sequences: Vec<&[u8]> = vec_sequences.iter().map(|v| v.as_slice()).collect();
        let offsets = vec![None; sequences.len()];

        let mut node = DualConsensusNode::new_root_node(&offsets, None, true).unwrap();
        node.activate_dual(sequences.as_slice(), b'A', b'C').unwrap();

        let weights1 = node.get_ed_weights(true, true);
        assert_eq!(weights1, vec![1.0 / 1.5, 0.5 / 1.5]);
        let weights2 = node.get_ed_weights(false, true);
        assert_eq!(weights2, vec![0.5 / 1.5, 1.0 / 1.5]);

        let weights1 = node.get_ed_weights(true, false);
        assert_eq!(weights1, vec![1.0, 0.0]);
        let weights2 = node.get_ed_weights(false, false);
        assert_eq!(weights2, vec![0.0, 1.0]);
    }

    // below here are mostly end-to-end DualConsensusDWFA tests

    #[derive(Debug, serde::Deserialize)]
    struct DualRecord {
        consensus: usize,
        edits: usize,
        sequence: String
    }

    /// Wrapper test function that loads a dual test from a csv file.
    /// Expected columns are "consensus" (1 or 2), "edits" (u64), and "sequence" (String).
    /// Returns a tuple of (sequences, DualConsensus).
    /// # Arguments
    /// * `filename` - the file path to load
    /// * `include_consensus` - if True, it will load the first consensus read into the sequences
    /// * `cost_mode` - the cost mode getting tested
    fn load_dual_csv_test(filename: &std::path::Path, include_consensus: bool, cost_mode: ConsensusCost) -> (Vec<Vec<u8>>, DualConsensus) {
        let mut sequences = vec![];
        let mut is_consensus1 = vec![];
        let mut ed1 = vec![];
        let mut ed2 = vec![];

        let mut con1: Option<Vec<u8>> = None;
        let mut con2: Option<Vec<u8>> = None;

        let mut csv_reader = csv::ReaderBuilder::new()
            .has_headers(true)
            .from_path(filename)
            .unwrap();
        for row in csv_reader.deserialize() {
            let record: DualRecord = row.unwrap();

            let is_con1 = record.consensus == 1;
            let edits = match cost_mode {
                ConsensusCost::L1Distance => record.edits,
                ConsensusCost::L2Distance => record.edits.pow(2)
            };
            let sequence = record.sequence.as_bytes().to_vec();

            if is_con1 {
                if con1.is_none() && edits == 0 {
                    // this is the first 0 ED for consensus 1
                    con1 = Some(sequence.clone());
                    if !include_consensus {
                        continue;
                    }
                }
                ed1.push(edits);
            } else {
                if con2.is_none() && edits == 0 {
                    // this is the first 0 ED for consensus 2
                    con2 = Some(sequence.clone());
                    if !include_consensus {
                        continue;
                    }
                }
                ed2.push(edits);
            }

            is_consensus1.push(is_con1);
            sequences.push(sequence);
        }

        // make sure that either we do not have consensus 2 OR consensus 1 comes before consensus 2
        assert!(con2.is_none() || con1.as_ref().unwrap() < con2.as_ref().unwrap());

        let consensus1 = Consensus::new(con1.unwrap(), cost_mode, ed1);
        let consensus2 = con2.map(|c2| Consensus::new(c2, cost_mode, ed2));
        let consensus = DualConsensus {
            consensus1,
            consensus2,
            is_consensus1,
            scores1: vec![None; sequences.len()],
            scores2: vec![None; sequences.len()]
        };

        (sequences, consensus)
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

        let (sequences, expected_consensus) = load_dual_csv_test(&PathBuf::from(filename), include_consensus, config.consensus_cost);

        // set queue size large since we do not want to worry about filtering for this test
        let mut consensus_dwfa = DualConsensusDWFA::with_config(config).unwrap();
        for sequence in sequences.iter() {
            consensus_dwfa.add_sequence(sequence).unwrap();
        }

        // make sure our alphabet has four symbols, ACGT
        assert_eq!(consensus_dwfa.alphabet().len(), 4);

        // now check that the consensus is the same as our sequence
        let consensus = consensus_dwfa.consensus().unwrap();
        assert_eq!(consensus, vec![expected_consensus]);
    }

    #[test]
    fn test_single_sequence() {
        let sequence = b"ACGTACGTACGT";
        let mut consensus_dwfa = DualConsensusDWFA::default();
        consensus_dwfa.add_sequence(sequence).unwrap();

        // make sure our alphabet has four symbols, ACGT
        assert_eq!(consensus_dwfa.alphabet().len(), 4);

        // now check that the consensus is the same as our sequence
        let consensus = consensus_dwfa.consensus().unwrap();
        assert_eq!(consensus, vec![DualConsensus {
            consensus1: Consensus::new(
                sequence.to_vec(),
                ConsensusCost::L1Distance,
                vec![0]
            ),
            consensus2: None,
            is_consensus1: vec![true],
            // these are not checked
            scores1: vec![],
            scores2: vec![]
        }]);
    }

    #[test]
    fn test_trio_sequence() {
        let sequence = b"ACGTACGTACGT";
        let sequence2 = b"ACGTACCTACGT";
        let mut consensus_dwfa = DualConsensusDWFA::default();

        // add the first sequence twice, it will be the consensus
        consensus_dwfa.add_sequence(sequence).unwrap();
        consensus_dwfa.add_sequence(sequence).unwrap();
        consensus_dwfa.add_sequence(sequence2).unwrap();

        // make sure our alphabet has four symbols, ACGT
        assert_eq!(consensus_dwfa.alphabet().len(), 4);

        // now check that the consensus found sequence
        let consensus = consensus_dwfa.consensus().unwrap();

        // sequence2 is alphabetically before sequence 1, so it will come first
        assert_eq!(consensus, vec![DualConsensus {
            consensus1: Consensus::new(
                sequence.to_vec(),
                ConsensusCost::L1Distance,
                vec![0, 0, 1]
            ),
            consensus2: None,
            is_consensus1: vec![true, true, true],
            // these are not checked
            scores1: vec![],
            scores2: vec![]
        }]);
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
        let mut consensus_dwfa = DualConsensusDWFA::default();
        for &sequence in sequences.iter() {
            consensus_dwfa.add_sequence(sequence).unwrap();
        }

        // make sure our alphabet has four symbols, ACGT
        assert_eq!(consensus_dwfa.alphabet().len(), 4);

        // now check that the consensus found the right thing and only that
        let consensus = consensus_dwfa.consensus().unwrap();
        // assert_eq!(consensus.len(), 1);
        assert_eq!(consensus, vec![DualConsensus {
            consensus1: Consensus::new(expected_consensus.to_vec(), ConsensusCost::L1Distance, vec![2, 2, 1]),
            consensus2: None,
            is_consensus1: vec![true; 3],
            // these are not checked
            scores1: vec![],
            scores2: vec![]
        }]);
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
        let mut consensus_dwfa = DualConsensusDWFA::with_config(
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
        // assert_eq!(consensus.len(), 1);
        assert_eq!(consensus, vec![DualConsensus {
            consensus1: Consensus::new(expected_consensus.to_vec(), ConsensusCost::L1Distance, vec![1, 1, 0]),
            consensus2: None,
            is_consensus1: vec![true; 3],
            // these are not checked
            scores1: vec![],
            scores2: vec![]
        }]);
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
        let mut consensus_dwfa = DualConsensusDWFA::with_config(
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
        // assert_eq!(consensus.len(), 1);
        assert_eq!(consensus, vec![DualConsensus {
            consensus1: Consensus::new(actual_consensus.to_vec(), ConsensusCost::L1Distance, vec![1, 0, 1]),
            consensus2: None,
            is_consensus1: vec![true; 3],
            // these are not checked
            scores1: vec![],
            scores2: vec![]
        }]);
    }

    #[test]
    fn test_dual_sequence() {
        let sequence     = b"ACGT";
        let alt_sequence = b"AGGT"; // single snp test
        let mut consensus_dwfa = DualConsensusDWFA::with_config(
            CdwfaConfigBuilder::default()
                .min_count(1)
                .build().unwrap()
        ).unwrap();
        consensus_dwfa.add_sequence(sequence).unwrap();
        consensus_dwfa.add_sequence(alt_sequence).unwrap();

        // make sure our alphabet has four symbols, ACGT
        assert_eq!(consensus_dwfa.alphabet().len(), 4);

        // now check that the consensus is the same as our sequence
        let consensus = consensus_dwfa.consensus().unwrap();
        assert_eq!(consensus, vec![DualConsensus {
            consensus1: Consensus::new(sequence.to_vec(), ConsensusCost::L1Distance, vec![0]),
            consensus2: Some(Consensus::new(alt_sequence.to_vec(), ConsensusCost::L1Distance, vec![0])),
            is_consensus1: vec![true, false],
            // these are not checked
            scores1: vec![],
            scores2: vec![]
        }]);
    }

    #[test]
    fn test_dual_unequal_001() {
        let sequence     = b"ACGT";
        let alt_sequence = b"AGGTA"; // single snp test
        let mut consensus_dwfa = DualConsensusDWFA::with_config(
            CdwfaConfigBuilder::default()
                .min_count(1)
                .build().unwrap()
        ).unwrap();
        consensus_dwfa.add_sequence(sequence).unwrap();
        consensus_dwfa.add_sequence(alt_sequence).unwrap();

        // make sure our alphabet has four symbols, ACGT
        assert_eq!(consensus_dwfa.alphabet().len(), 4);

        // now check that the consensus is the same as our sequence
        let consensus = consensus_dwfa.consensus().unwrap();
        assert_eq!(consensus, vec![DualConsensus {
            consensus1: Consensus::new(sequence.to_vec(), ConsensusCost::L1Distance, vec![0]),
            consensus2: Some(Consensus::new(alt_sequence.to_vec(), ConsensusCost::L1Distance, vec![0])),
            is_consensus1: vec![true, false],
            // these are not checked
            scores1: vec![],
            scores2: vec![]
        }]);
    }

    #[test]
    fn test_dual_unequal_002() {
        let sequence     = b"ACGTA";
        let alt_sequence = b"AGGT"; // single snp test
        let mut consensus_dwfa = DualConsensusDWFA::with_config(
            CdwfaConfigBuilder::default()
                .min_count(1)
                .build().unwrap()
        ).unwrap();
        consensus_dwfa.add_sequence(sequence).unwrap();
        consensus_dwfa.add_sequence(alt_sequence).unwrap();

        // make sure our alphabet has four symbols, ACGT
        assert_eq!(consensus_dwfa.alphabet().len(), 4);

        // now check that the consensus is the same as our sequence
        let consensus = consensus_dwfa.consensus().unwrap();
        assert_eq!(consensus, vec![DualConsensus {
            consensus1: Consensus::new(sequence.to_vec(), ConsensusCost::L1Distance, vec![0]),
            consensus2: Some(Consensus::new(alt_sequence.to_vec(), ConsensusCost::L1Distance, vec![0])),
            is_consensus1: vec![true, false],
            // these are not checked
            scores1: vec![],
            scores2: vec![]
        }]);
    }

    #[test]
    fn test_dual_noise_before_variation() {
        // this test includes a noisy C insertion; it will cause a dual node split that get's corrected later when the true variation is reached
        let con1 = b"ACGTACGTACGT";
        let con2 = b"ACGTACGTCCCT"; // two SNPs towards the end that will be consistent
        let sequences = [
            b"ACGTACGTACGT".to_vec(),
            b"ACCGTACGTACGT".to_vec(), // noisy C insert
            b"ACGTACGTACGT".to_vec(),
            b"ACGTACGTCCCT".to_vec(),
            b"ACGTACGTCCCT".to_vec(),
            b"ACCGTACGTCCCT".to_vec()  // noisy C insert
        ];

        // set queue size large since we do not want to worry about filtering for this test
        let mut consensus_dwfa = DualConsensusDWFA::with_config(
            CdwfaConfigBuilder::default()
                .min_count(1)
                .max_queue_size(1000)
                .build().unwrap()
        ).unwrap();
        for sequence in sequences.iter() {
            consensus_dwfa.add_sequence(sequence).unwrap();
        }

        // make sure our alphabet has four symbols, ACGT
        assert_eq!(consensus_dwfa.alphabet().len(), 4);

        // now check that the consensus is the same as our sequence
        let consensus = consensus_dwfa.consensus().unwrap();
        assert_eq!(consensus, vec![DualConsensus {
            consensus1: Consensus::new(con1.to_vec(), ConsensusCost::L1Distance, vec![0, 1, 0]),
            consensus2: Some(Consensus::new(con2.to_vec(), ConsensusCost::L1Distance, vec![0, 0, 1])),
            is_consensus1: vec![true, true, true, false, false, false],
            // these are not checked
            scores1: vec![],
            scores2: vec![]
        }]);
    }

    #[test]
    fn test_multi_extension() {
        // test the scenario where 3+ extensions are identified
        let con1 = b"ACGTACGTACGT";
        let con2 = b"ACGTACGTCCCT"; // two SNPs towards the end that will be consistent
        let sequences = [
            // matches con 1
            b"ACGTACGTACGT".to_vec(),
            b"ACGTACGTACGT".to_vec(), 
            b"ACGTACGTGCGT".to_vec(),// A read as G, should create multi-extension for exploration
            // matches con2
            b"ACGTACGTCCCT".to_vec(), 
            b"ACGTACGTCCCT".to_vec(),
            b"ACGTACGTGCCT".to_vec()  // first C read as G, should create multi-extension for exploration
        ];

        // set queue size large since we do not want to worry about filtering for this test
        let mut consensus_dwfa = DualConsensusDWFA::with_config(
            CdwfaConfigBuilder::default()
                .min_count(1)
                .max_queue_size(1000)
                .build().unwrap()
        ).unwrap();
        for sequence in sequences.iter() {
            consensus_dwfa.add_sequence(sequence).unwrap();
        }

        // make sure our alphabet has four symbols, ACGT
        assert_eq!(consensus_dwfa.alphabet().len(), 4);

        // now check that the consensus is the same as our sequence
        let consensus = consensus_dwfa.consensus().unwrap();
        assert_eq!(consensus, vec![DualConsensus {
            consensus1: Consensus::new(con1.to_vec(), ConsensusCost::L1Distance, vec![0, 0, 1]),
            consensus2: Some(Consensus::new(con2.to_vec(), ConsensusCost::L1Distance, vec![0, 0, 1])),
            is_consensus1: vec![true, true, true, false, false, false],
            // these are not checked
            scores1: vec![],
            scores2: vec![]
        }]);
    }

    #[test]
    fn test_equal_options() {
        // test the scenario multiple possible dual consensuses are possible
        let con1 = b"ACGTACGTACGT"; // 00
        let con2 = b"ACGTCCGTCCGT"; // 11
        let con3 = b"ACGTACGTCCGT"; // 01
        let con4 = b"ACGTCCGTACGT"; // 10
        let sequences = [
            con1.to_vec(),
            con2.to_vec(),
            con3.to_vec(),
            con4.to_vec()
        ];

        // set queue size large since we do not want to worry about filtering for this test
        let mut consensus_dwfa = DualConsensusDWFA::with_config(
            CdwfaConfigBuilder::default()
                .min_count(1)
                .max_queue_size(1000)
                .build().unwrap()
        ).unwrap();
        for sequence in sequences.iter() {
            consensus_dwfa.add_sequence(sequence).unwrap();
        }

        // make sure our alphabet has four symbols, ACGT
        assert_eq!(consensus_dwfa.alphabet().len(), 4);

        // this result actually has 6 equally possible dual consensuses, each of which is ed = 2 in total
        let consensus = consensus_dwfa.consensus().unwrap();
        assert_eq!(consensus.len(), 6);
        for dc in consensus.iter() {
            // make sure it's dual
            assert!(dc.consensus2.is_some());
            // make sure the edits are 2
            assert_eq!(2, dc.consensus1.scores().iter().sum::<usize>()+dc.consensus2.as_ref().unwrap().scores().iter().sum::<usize>());
            
            // future: if we need to test the exact result, we're going to have to add it by hand later
        }
    }

    #[test]
    fn test_tail_extension() {
        // test the scenario we have a simple tail extension
        // for now, this will NOT create a dual solution, but will only return one possible solution
        // future: we may want it to somehow figure out that this is dual, but this seems like a minor edge case for now
        let con1 = b"ACGT";
        let con2 = b"ACGTT"; // one consensus is 1 bp longer in the tail
        let sequences = [
            con1.to_vec(),
            con2.to_vec()
        ];

        // set queue size large since we do not want to worry about filtering for this test
        let mut consensus_dwfa = DualConsensusDWFA::with_config(
            CdwfaConfigBuilder::default()
                .min_count(1)
                .max_queue_size(1000)
                .build().unwrap()
        ).unwrap();
        for sequence in sequences.iter() {
            consensus_dwfa.add_sequence(sequence).unwrap();
        }

        // make sure our alphabet has four symbols, ACGT
        assert_eq!(consensus_dwfa.alphabet().len(), 4);

        // now check that the consensus is the same as our sequence
        let consensus = consensus_dwfa.consensus().unwrap();
        assert_eq!(consensus, vec![DualConsensus {
            consensus1: Consensus::new(con1.to_vec(), ConsensusCost::L1Distance, vec![0, 1]),
            consensus2: None,
            is_consensus1: vec![true, true],
            // these are not checked
            scores1: vec![],
            scores2: vec![]
        }, 
        DualConsensus {
            consensus1: Consensus::new(con2.to_vec(), ConsensusCost::L1Distance, vec![1, 0]),
            consensus2: None,
            is_consensus1: vec![true, true],
            // these are not checked
            scores1: vec![],
            scores2: vec![]
        }]);
    }

    #[test]
    fn test_csv_dual_001() {
        run_test_file("./tests/dual_001.csv", true, None);
    }

    #[test]
    fn test_dual_max_ed_delta() {
        // this is the same test as above BUT we have restricted the dual_max_ed_delta to too low, resulting in one read getting mis-assigned
        // this is not normal, but we are testing that the functionality is correct
        let cost_mode = ConsensusCost::L1Distance;
        let (sequences, expected_consensus) = load_dual_csv_test(&PathBuf::from("./tests/dual_001.csv"), true, cost_mode);
        
        // third read gets swapped and has a worse ED than before
        let expected_consensus = vec![
            DualConsensus { 
                consensus1: Consensus::new(
                    expected_consensus.consensus1.sequence().to_vec(), 
                    ConsensusCost::L1Distance, 
                    vec![0, 4, 4, 2] // delete the third entry here
                ),
                consensus2: Some(Consensus::new(
                    expected_consensus.consensus2.as_ref().unwrap().sequence().to_vec(), 
                    ConsensusCost::L1Distance,
                    vec![3, 0, 0, 0, 0, 0] // shift it here, with a worse ED
                )), 
                is_consensus1: vec![true, true, false, true, true, false, false, false, false, false], // mark third from true -> false
                // these are not checked
                scores1: vec![],
                scores2: vec![]
            }
        ];

        // set queue size large since we do not want to worry about filtering for this test
        let mut consensus_dwfa = DualConsensusDWFA::with_config(
            CdwfaConfigBuilder::default()
                .wildcard(Some(b'*'))
                .dual_max_ed_delta(0) // set this intentionally low to bundle a sequence incorrectly in our result
                .build().unwrap()
        ).unwrap();
        for sequence in sequences.iter() {
            consensus_dwfa.add_sequence(sequence).unwrap();
        }

        // make sure our alphabet has four symbols, ACGT
        assert_eq!(consensus_dwfa.alphabet().len(), 4);

        // now check that the consensus is the same as our sequence
        let consensus = consensus_dwfa.consensus().unwrap();
        assert_eq!(consensus, expected_consensus);
    }

    #[test]
    fn test_csv_length_gap_001() {
        run_test_file("./tests/length_gap_001.csv", false, Some(
            CdwfaConfigBuilder::default()
                .wildcard(Some(b'*'))
                .min_count(2)
                .dual_max_ed_delta(5)
                .max_queue_size(1000)
                .consensus_cost(ConsensusCost::L2Distance)
                .build().unwrap()
        ));
    }

    #[test]
    fn test_csv_early_termination_001() {
        run_test_file("./tests/dual_early_termination_001.csv", true, Some(
            CdwfaConfigBuilder::default()
                .wildcard(Some(b'*'))
                .allow_early_termination(true)
                .build().unwrap()
        ));
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
        let mut consensus_dwfa = DualConsensusDWFA::with_config(
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
        assert!(!consensus[0].is_dual());
        assert_eq!(expected_consensus, consensus[0].consensus1().sequence());
        assert_eq!(&[0, 0, 0], consensus[0].consensus1().scores());
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
        let mut consensus_dwfa = DualConsensusDWFA::with_config(
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