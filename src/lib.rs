/*!
# waffle_con
This library provides access to consensus algorithms based on building a consensus with dynamic WFA exploration.

Key benefits:
* Ability to generate different consensus types: single consensus, dual consensus (e.g., diplotype), or multi-consensus (e.g., unknown number of consensuses)
* Generally efficient for low noise inputs
* Does not rely on a sequence backbone

Performance notes:
* The underlying algorithm scales with WFA, so high error rates will increase compute time and memory consumption
* Certain branching patterns are more expensive than others, we plan to improve this in the future

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

/// Configuration for ConsensusDWFA
pub mod cdwfa_config;
/// Main functionality for the consensus component
pub mod consensus;
/// Functionality for a 2-way consensus
pub mod dual_consensus;
/// Main functionality for a dynamic WFA, which will just focus on a single string
pub mod dynamic_wfa;
/// Utility for generating examples
pub mod example_gen;
/// Main functionality for a multiple-consensus component
pub mod multi_consensus;
/// Utility for tracking a pqueue without force deleting items
pub mod pqueue_tracker;
/// Main functionality for a priority multiple-consensus component, allowing differentiation on multiple sequences in priority order
pub mod priority_consensus;
/// Basic pair-wise alignment utilities
pub mod sequence_alignment;
