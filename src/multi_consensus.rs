
/*!
Initially, this module provided access to a MultiConsensus approach. 
This approach has been deprecated for now (see PriorityConsensus), but the output type `MultiConsensus` is still useful in many contexts.
*/

use crate::consensus::Consensus;

/// Contains a final multi-consensus result
#[derive(Debug, PartialEq)]
pub struct MultiConsensus {
    /// Results for the consensuses
    consensuses: Vec<Consensus>,
    /// For each input sequence, this is the index of the assigned consensus
    sequence_indices: Vec<usize>
}

impl MultiConsensus {
    /// General Constructor for MultiConsensus.
    /// Of note, this will re-order alleles alpha-numerically, allowing for predictable outputs.
    /// This re-ordering will alter the sequence_indices to match the new order.
    /// # Arguments
    /// * `consensuses` - the actual consensus sequences with scores
    /// * `sequence_indices` - a mapping point out which consensus the sequence ended up in
    pub fn new(mut consensuses: Vec<Consensus>, sequence_indices: Vec<usize>) -> MultiConsensus {
        // sort the consensuses
        let mut ordered_indices = (0..consensuses.len()).collect::<Vec<usize>>();
        ordered_indices.sort_by(|i, j| {
            let ci = consensuses[*i].sequence();
            let cj = consensuses[*j].sequence();
            ci.cmp(cj)
        });

        // now create the reverse lookup table
        let mut reverse_lookup = vec![usize::MAX; consensuses.len()];
        for (new_index, &old_index) in ordered_indices.iter().enumerate() {
            reverse_lookup[old_index] = new_index;
        }

        // we can now just re-sort the consensuses in place, this is generally better than cloning and such
        // TODO: is there some way to re-use the ordered indices above? not sure there is that is more efficient
        consensuses.sort_by(|c1, c2| {
            c1.sequence().cmp(c2.sequence())
        });

        // re-map the sequence indices
        let sequence_indices: Vec<usize> = sequence_indices.iter()
            .map(|&bci| reverse_lookup[bci])
            .collect();

        MultiConsensus {
            consensuses,
            sequence_indices
        }
    }

    // Getters
    pub fn consensuses(&self) -> &[Consensus] {
        &self.consensuses
    }

    pub fn sequence_indices(&self) -> &[usize] {
        &self.sequence_indices
    }
}

#[cfg(test)]
mod tests {
    use crate::cdwfa_config::ConsensusCost;

    use super::*;

    #[test]
    fn test_multiconsensus_sort() {
        let consensuses = vec![
            Consensus::new(b"ACGT".to_vec(), ConsensusCost::L1Distance, vec![0]),
            Consensus::new(b"TGCA".to_vec(), ConsensusCost::L1Distance, vec![0]),
            Consensus::new(b"AAAA".to_vec(), ConsensusCost::L1Distance, vec![0]),
        ];
        let sequence_indices = vec![
            2, 0, 1
        ];
        let multicon = MultiConsensus::new(consensuses, sequence_indices);

        // these should now be sorted by sequence with the sequence_indices adjusted to match the new order
        assert_eq!(multicon, MultiConsensus {
            consensuses: vec![
                Consensus::new(b"AAAA".to_vec(), ConsensusCost::L1Distance, vec![0]),
                Consensus::new(b"ACGT".to_vec(), ConsensusCost::L1Distance, vec![0]),
                Consensus::new(b"TGCA".to_vec(), ConsensusCost::L1Distance, vec![0]),
            ], 
            sequence_indices: vec![0, 1, 2]
        });
    }
}