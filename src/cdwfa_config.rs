
/*!
Contains configuration information for the consensus DWFA algorithm.
Typical usage is to the use the builder to construct the config, e.g.
```
use waffle_con::cdwfa_config::{CdwfaConfig, CdwfaConfigBuilder, ConsensusCost};
let config: CdwfaConfig = CdwfaConfigBuilder::default()
    .consensus_cost(ConsensusCost::L2Distance)
    .wildcard(Some(b'N'))
    .build()
    .unwrap();
```
*/

use simple_error::bail;

use crate::dynamic_wfa::DWFALiteConfig;

/// Enumeration of difference scoring types for a consensus.
/// Initially just using L1 distance, which is the sum of edit distance across all inputs.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub enum ConsensusCost {
    /// Minimizes the total edit distance across all sequences
    #[default]
    L1Distance,
    /// Minimizes the square of the edit distance across all sequences
    L2Distance
}

/**
Contains configuration information for the consensus DWFA algorithm.
Typical usage is to the use the builder to construct the config, e.g.
```
use waffle_con::cdwfa_config::{CdwfaConfig, CdwfaConfigBuilder, ConsensusCost};
let config: CdwfaConfig = CdwfaConfigBuilder::default()
    .consensus_cost(ConsensusCost::L2Distance)
    .wildcard(Some(b'N'))
    .build()
    .unwrap();
```
*/
#[derive(derive_builder::Builder, Clone, Debug)]
#[builder(default)]
pub struct CdwfaConfig {
    /// The consensus scoring cost
    pub consensus_cost: ConsensusCost,
    /// Maximum queue size, which controls how many active branches we allow during exploration
    pub max_queue_size: usize,
    /// Maximum capacity, controls how many nodes of each length we process
    pub max_capacity_per_size: usize,
    /// Maximum return size, which controls how many equal returns we track
    pub max_return_size: usize,
    /// Maximum explored nodes without a constraint, this prevents hyper-branching in truly ambiguous regions
    pub max_nodes_wo_constraint: u64,
    /// Minimum number of occurrences of a candidate extension to get used (largest is always used regardless)
    pub min_count: u64,
    /// Minimum fraction of sequences of a candidate extension to get used
    pub min_af: f64,
    /// For dual/multi-consensus, this will weight the nominated extension by the current edit distance, "accelerating" convergence
    pub weighted_by_ed: bool,
    /// Enables an optional wildcard character that will match anything
    pub wildcard: Option<u8>,
    /// Dual mode DWFA edit distance pruning threshold; if the DWFAs for the two options diverge more than this threshold, the worse one stops getting tracked
    pub dual_max_ed_delta: usize,
    // if true, then this will not penalize input sequences that are shorter than the final consensus
    pub allow_early_termination: bool,
    /// if true, this will automatically shift offsets downwards if nothing starts at "0"
    pub auto_shift_offsets: bool,
    /// The number of bases before the last_offset to search for an optimal start point
    pub offset_window: usize,
    /// The number of bases to use in the comparison for calculating best optimal start point
    pub offset_compare_length: usize,
    /// Optional hard maximum edit distance for one read. `None` disables the hard cap.
    pub max_edit_distance: Option<usize>,
    /// Optional edit-distance cap as a fraction of the baseline read length. `None` disables the fractional cap.
    /// When both caps are set, the tighter one is used.
    pub max_edit_distance_fraction: Option<f64>,
}

impl Default for CdwfaConfig {
    fn default() -> Self {
        Self { 
            // L1 v. L2 is an open question
            consensus_cost: ConsensusCost::L1Distance,
            // 20 is relatively small, but this seems to work out in practice for our low-error sequences
            max_queue_size: 20,
            // set to the same by default
            max_capacity_per_size: 20,
            // Realistically, anything more than 10 is not particularly useful
            max_return_size: 10,
            // lower values help constrain hyper-branching scenarios, we probably should not go below 10 though
            max_nodes_wo_constraint: 1000,
            // 3 seems reasonable
            min_count: 3,
            // by default, we will just let the raw count work
            min_af: 0.0,
            // by default, it's not clear we want to do this
            weighted_by_ed: false,
            // by default, we likely do not want a wildcard symbol
            wildcard: None,
            // if the options have diverged by 20, seems like one is a clear candidate
            dual_max_ed_delta: 20,
            // to keep with our existing tests, we will not allow this by default
            allow_early_termination: false,
            // someone might not want this, but most will
            auto_shift_offsets: true,
            // these were just what we started with
            offset_window: 50,
            offset_compare_length: 50,
            // by default, do not cap how far a read may diverge
            max_edit_distance: None,
            max_edit_distance_fraction: None,
        }
    }
}

impl CdwfaConfig {
    /// Edit-distance cap for a baseline of `baseline_len` bases.
    /// Uses the tighter of [`Self::max_edit_distance`] and `floor(max_edit_distance_fraction * baseline_len)`.
    /// Returns `None` when neither limit is set.
    /// # Errors
    /// * if `max_edit_distance_fraction` is negative or non-finite
    pub fn max_edit_distance_for(&self, baseline_len: usize) -> Result<Option<usize>, Box<dyn std::error::Error>> {
        // first, resolve the dynamic max ED
        let dynamic_max = if let Some(fraction) = self.max_edit_distance_fraction {
            if !fraction.is_finite() || fraction < 0.0 {
                bail!("max_edit_distance_fraction must be finite and non-negative");
            }
            Some((fraction * baseline_len as f64).floor() as usize)
        } else {
            None
        };

        // then, resolve the final max ED
        Ok(match (self.max_edit_distance, dynamic_max) {
            // if both are set, use the tighter of the two
            (Some(hard_max), Some(dynamic_max)) => Some(hard_max.min(dynamic_max)),
            // otherwise, use any that are set
            (hard_max, dynamic_max) => hard_max.or(dynamic_max),
        })
    }

    /// Builds the fixed DWFA options for one baseline read of `baseline_len` bases.
    /// # Errors
    /// * if `max_edit_distance_fraction` is negative or non-finite
    pub fn dwfa_lite_config_for(&self, baseline_len: usize) -> Result<DWFALiteConfig, Box<dyn std::error::Error>> {
        Ok(DWFALiteConfig {
            wildcard: self.wildcard,
            allow_early_termination: self.allow_early_termination,
            max_edit_distance: self.max_edit_distance_for(baseline_len)?,
        })
    }
}