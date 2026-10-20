use rustc_hash::FxHashMap as HashMap;
use simple_error::bail;

use super::dwfa_config::DWFALiteConfig;
use super::votes::{ExtensionVotes, VoteScratch};
use super::wavefront::{
    grow_wavefront, grow_wide_wavefront, grow_without_matches, match_length, WavefrontGrower, FRONTIER_BITS,
};

/// Lifecycle of a [`DWFALite`].
/// A DWFA starts [`Active`](DWFALiteState::Active). 
/// Finalizing it or passing the edit-distance cap moves it out of that state.
#[derive(Clone, Copy, Debug, Default, Eq, Hash, PartialEq)]
pub enum DWFALiteState {
    /// The DWFA can still be extended.
    #[default]
    Active,
    /// `finalize` has completed and the DWFA can no longer be extended.
    Finalized,
    /// Another edit would pass `max_edit_distance`. The distance stays at the cap and this DWFA stops voting.
    ExceededEditDistanceLimit,
}

/// The core dynamic WFA structure for a lite implementation.
/// It is structured such that all the sequences being built are maintained **outside** of this struct (hence the "lite").
/// Essentially, it is keep the wavefront information, but all sequences are updated and tracked elsewhere.
/// As a result, all updates must pass the sequences so the information can get tracked.
/// If the outside sequences are changed in any way other than appending, then this may become desynchronized.
/// Conceptually, if this is a 2D grid, the baseline sequence will go from top to bottom (y-axis) and the other sequence will go from left to right (x-axis).
/// This means each character we add to `other_seq` will add a new _column_ to the grid.
#[derive(Debug, Eq, Hash, PartialEq)]
pub struct DWFALite {
    /// Fixed matching and edit-distance options for this DWFA.
    config: DWFALiteConfig,
    /// Whether this DWFA is still active, finalized, or frozen at its edit-distance cap.
    state: DWFALiteState,
    /// The minimum edit distance from `other_seq` to `baseline_seq` allowing for a large port of the tail of the `baseline_seq` to get ignored.
    edit_distance: usize,
    /// This is the core wavefront for the current `edit_distance`.
    /// It should always be length 2*`edit_distance`-1.
    /// For our mental model, the median diagonal means we have used the same number of bases in our baseline and other.
    /// Anything below the median diagonal means we have used more baseline than other.
    /// Anything above the median diagonal means we have used more other than baseline.
    /// The value stored will always be the number of bases consumed in `other_seq`.
    /// After a proper character addition, _at least one_ of these should always equal the length of `other_seq` - `offset`.
    wavefront: Vec<usize>,
    // this is an offset into `other_seq` that the baseline starts
    offset: usize,
    /// Cached maximum of `wavefront`. Values only grow within one edit distance.
    max_wavefront: usize,
    /// Cached farthest position reached in the baseline, see `maximum_baseline_distance()`.
    max_baseline: usize,
    /// Length of `other_seq` at the last extension.
    /// A diagonal that stopped before this position stopped on a mismatch or the baseline end.
    /// Because `other_seq` only grows, only diagonals at this frontier can extend again.
    extended_len: usize,
    /// Bit `i` is set when `wavefront[i] + offset == extended_len`, i.e., diagonal `i` is at the frontier.
    /// `None` when unknown (after `set_offset()`) or when the wavefront is wider than `FRONTIER_BITS`.
    frontier: Option<u128>,
}

/// Spare wavefront capacity that a clone gets. Clones are usually extended at once, and each edit adds two diagonals.
/// With spare capacity, the next few edits grow the wavefront in place instead of in a new allocation.
const CLONE_SPARE_DIAGONALS: usize = 8;

impl Clone for DWFALite {
    fn clone(&self) -> Self {
        let mut wavefront = Vec::with_capacity(self.wavefront.len() + CLONE_SPARE_DIAGONALS);
        wavefront.extend_from_slice(&self.wavefront);
        DWFALite {
            config: self.config.clone(),
            state: self.state,
            edit_distance: self.edit_distance,
            wavefront,
            offset: self.offset,
            max_wavefront: self.max_wavefront,
            max_baseline: self.max_baseline,
            extended_len: self.extended_len,
            frontier: self.frontier,
        }
    }

    /// Reuses the wavefront allocation of `self`.
    /// # Arguments
    /// * `source` - the DWFA to copy
    fn clone_from(&mut self, source: &Self) {
        self.wavefront.clear();
        self.wavefront.reserve(source.wavefront.len() + CLONE_SPARE_DIAGONALS);
        self.wavefront.extend_from_slice(&source.wavefront);
        self.config.clone_from(&source.config);
        self.state = source.state;
        self.edit_distance = source.edit_distance;
        self.offset = source.offset;
        self.max_wavefront = source.max_wavefront;
        self.max_baseline = source.max_baseline;
        self.extended_len = source.extended_len;
        self.frontier = source.frontier;
    }
}

impl Default for DWFALite {
    fn default() -> Self {
        DWFALite {
            config: DWFALiteConfig::default(),
            state: DWFALiteState::Active,
            edit_distance: 0,
            wavefront: vec![0],
            offset: 0,
            max_wavefront: 0,
            max_baseline: 0,
            extended_len: 0,
            // the single diagonal is at `other_seq` position 0, which is the frontier
            frontier: Some(1),
        }
    }
}
    

impl DWFALite {
    /// Creates a new dynamic WFA instance from fixed options.
    /// # Arguments
    /// * `config` - wildcard, early termination, and the absolute edit-distance cap
    pub fn new(config: DWFALiteConfig) -> DWFALite {
        DWFALite {
            config,
            ..Default::default()
        }
    }

    /// Sets an offset into `other_seq` where the baseline starts.
    /// This will effectively ignore the first `offset` characters of `other_seq`, as if they were not present.
    /// # Arguments
    /// * `offset` - the number of characters to ignore
    pub fn set_offset(&mut self, offset: usize) {
        self.offset = offset;
        // the frontier is relative to the offset, so the next extension must visit every diagonal
        self.extended_len = 0;
        self.frontier = None;
    }

    /// This is the main function to extend the WFA with a new symbol.
    /// This function will handle all of the extending and potential edit distance increases as well.
    /// # Arguments
    /// * `baseline_seq` - the baseline sequence, theoretically fixed
    /// * `other_seq` - the other sequence, typically getting updates
    /// # Errors
    /// * If this function is called after `finalize()` has been called.
    #[inline]
    pub fn update(&mut self, baseline_seq: &[u8], other_seq: &[u8]) -> Result<usize, Box<dyn std::error::Error>> {
        // try the fast path, fallback to update_general if it fails
        match self.advance_single_frontier(baseline_seq, other_seq) {
            Some(_) => Ok(self.edit_distance),
            None => self.update_general(baseline_seq, other_seq),
        }
    }

    /// `update()` and then `vote()`. When the fast path of `update()` applies, the vote is the next baseline symbol
    /// on the frontier diagonal, so it needs no second look at the wavefront.
    /// # Arguments
    /// * `baseline_seq` - the baseline sequence, theoretically fixed
    /// * `other_seq` - the other sequence, typically getting updates
    /// * `weight` - the total weight of this DWFA's vote
    /// * `votes` - the votes of all DWFAs so far
    /// * `scratch` - buffers for a DWFA with more than one frontier diagonal
    /// # Errors
    /// * If this function is called after `finalize()` has been called.
    #[inline(always)]
    pub(crate) fn update_and_vote(
        &mut self, baseline_seq: &[u8], other_seq: &[u8], weight: f64, votes: &mut ExtensionVotes, scratch: &mut VoteScratch
    ) -> Result<usize, Box<dyn std::error::Error>> {
        match self.advance_single_frontier(baseline_seq, other_seq) {
            Some(next_offset) => {
                // single wavefront case, just add weight once
                if let Some(&symbol) = baseline_seq.get(next_offset) {
                    votes.add_single(symbol, weight);
                }
                Ok(self.edit_distance)
            },
            None => {
                // general case, we may have multiple frontier wavefronts
                let edit_distance = self.update_general(baseline_seq, other_seq)?;
                self.vote(baseline_seq, other_seq, weight, votes, scratch);
                Ok(edit_distance)
            },
        }
    }

    /// The most common case of `update()`, without the general code. It applies when one symbol was appended,
    /// exactly one diagonal is at the frontier, and that diagonal matches the new symbol. Then the diagonal moves
    /// one step, and nothing else changes. It returns the new baseline offset of the diagonal.
    /// In all other cases, this returns `None` and changes nothing.
    /// For speed, it only handles a frontier diagonal in the low 64 bits of the bitmask.
    /// # Arguments
    /// * `baseline_seq` - the baseline sequence, theoretically fixed
    /// * `other_seq` - the other sequence, typically getting updates
    #[inline(always)]
    fn advance_single_frontier(&mut self, baseline_seq: &[u8], other_seq: &[u8]) -> Option<usize> {
        // the last update must have reached the end of `other_seq`, so that the frontier diagonal holds the maximum
        let bits = self.frontier?;
        let low_bits = bits as u64;

        // TODO: we could use u128 short circuiting here, but initial tests showed it was slightly slower
        //       could be worth revisiting in the future

        // short circuit if our assumptions are violated
        if self.state != DWFALiteState::Active || (bits >> 64) != 0 || !low_bits.is_power_of_two() ||
            other_seq.len() != self.extended_len + 1 || self.offset + self.max_wavefront != self.extended_len {
            // if any condition is violated, we have to take the general path
            return None;
        }
        // the frontier diagonal is at `extended_len`, so it holds the cached maximum; this saves a dependent load
        let i = low_bits.trailing_zeros() as usize;
        let d = self.max_wavefront;
        debug_assert_eq!(d, self.wavefront[i]);

        // check that the new symbol matches the frontier diagonal
        let baseline_offset = d + self.edit_distance - i;
        let symbol = *baseline_seq.get(baseline_offset)?;
        if symbol != other_seq[d + self.offset] && self.config.wildcard != Some(symbol) {
            // it doesn't, so we have to take the general path
            return None;
        }

        // update the wavefront and cached maxima
        self.wavefront[i] = d + 1;
        self.max_wavefront = d + 1;
        self.max_baseline = self.max_baseline.max(baseline_offset + 1);
        self.extended_len = other_seq.len();
        Some(baseline_offset + 1)
    }

    /// `update()` for all cases.
    /// # Arguments
    /// * `baseline_seq` - the baseline sequence, theoretically fixed
    /// * `other_seq` - the other sequence, typically getting updates
    /// # Errors
    /// * If this function is called after `finalize()` has been called.
    #[inline(never)]
    fn update_general(&mut self, baseline_seq: &[u8], other_seq: &[u8]) -> Result<usize, Box<dyn std::error::Error>> {
        match self.state {
            DWFALiteState::Finalized => {
                bail!("Cannot push more bases after finalizing a DWFA");
            },
            DWFALiteState::ExceededEditDistanceLimit => return Ok(self.edit_distance),
            DWFALiteState::Active => {}
        };

        // maximally extend the diagonals at the frontier; the others cannot move
        self.extend(baseline_seq, other_seq);
        
        // check how it looks
        let mut maximum_distance = self.maximum_other_distance();
        while maximum_distance < other_seq.len() && !(self.config.allow_early_termination && self.reached_baseline_end(baseline_seq)){
            if self.edit_distance_at_limit() {
                // another edit would pass the cap; freeze and stop contributing
                self.state = DWFALiteState::ExceededEditDistanceLimit;
                return Ok(self.edit_distance);
            }
            // increase the edit distance, re-extension happens automatically
            self.increase_edit_distance(baseline_seq, other_seq)?;

            // recalculate this
            maximum_distance = self.maximum_other_distance();
        }

        // final assertion just to make sure we don't break anything
        assert!(
            maximum_distance == other_seq.len() || 
            (self.config.allow_early_termination && self.maximum_baseline_distance() == baseline_seq.len())
        );
        Ok(self.edit_distance)
    }

    /// This will take the current wavefront and try to extend each one along it's current diagonal.
    /// It also refreshes the cached maxima.
    /// This should also be stable such that calling it after finalizing has no impact.
    /// # Arguments
    /// * `baseline_seq` - the baseline sequence, theoretically fixed
    /// * `other_seq` - the other sequence, typically getting updates
    fn extend(&mut self, baseline_seq: &[u8], other_seq: &[u8]) {
        // `i` is the index in the wavefront
        // `i // 2` is always the middle diagonal
        // anything less than that has used more baseline symbols
        // anything more than that has used more other symbols
        // example at edit distance 1: baseline offset would equal [d+1, d, d-1], which corresponds to
        // a skipped base in primary, equal usage, and a skipped base in other

        // if we truncate a wavefront (e.g., a bounded size); then `i` below will be too small for every wavefront from the start
        //     that gets truncated; this means we want the below formula to be `*d + self.edit_distance - (i + truncate_prefix)`
        //     where `truncate_prefix` is the number we have cut off the front (in total).
        //     IMPORTANT: this formula will need to get updated everywhere that `baseline_offset` is calculated

        let edit_distance = self.edit_distance;
        let offset = self.offset;
        let wildcard = self.config.wildcard;
        let other_len = other_seq.len();

        match self.frontier {
            // a shorter `other_seq` breaks the append-only contract, so that takes the full path below
            Some(mut bits) if other_len >= self.extended_len => {
                // only the frontier diagonals can move; values only grow, so the cached maxima stay valid lower bounds
                let mut new_bits: u128 = 0;
                while bits != 0 {
                    // get the index of the lowest set bit, and unset it
                    let i = bits.trailing_zeros() as usize;
                    bits &= bits - 1;

                    // get the diagonal and extend it by the match length
                    let d = &mut self.wavefront[i];
                    *d += match_length(baseline_seq, other_seq, *d + edit_distance - i, *d + offset, wildcard);

                    // if the diagonal is at the end, then it is still frontier
                    if *d + offset == other_len {
                        new_bits |= 1 << i;
                    }
                    self.max_wavefront = self.max_wavefront.max(*d);
                    self.max_baseline = self.max_baseline.max(*d + edit_distance - i);
                }
                self.frontier = Some(new_bits);
            },
            // general case when the frontier cannot short circuit
            _ => {
                let frontier = if other_len < self.extended_len { 0 } else { self.extended_len };
                let mut max_wavefront = 0;
                let mut max_baseline = 0;
                let mut new_bits: u128 = 0;

                // go through every diagonal and extend it
                for (i, d) in self.wavefront.iter_mut().enumerate() {
                    if *d + offset >= frontier {
                        *d += match_length(baseline_seq, other_seq, *d + edit_distance - i, *d + offset, wildcard);
                    }
                    if i < FRONTIER_BITS && *d + offset == other_len {
                        new_bits |= 1 << i;
                    }
                    max_wavefront = max_wavefront.max(*d);
                    max_baseline = max_baseline.max(*d + edit_distance - i);
                }
                self.max_wavefront = max_wavefront;
                self.max_baseline = max_baseline;
                self.frontier = (self.wavefront.len() <= FRONTIER_BITS).then_some(new_bits);
            },
        };
        self.extended_len = other_len;
    }

    /// This will increase the edit distance for this DWFA and create a new larger wavefront.
    /// It also extends each new diagonal and refreshes the cached maxima, as `extend()` does.
    /// # Arguments
    /// * `baseline_seq` - the baseline sequence, theoretically fixed
    /// * `other_seq` - the other sequence, typically getting updates
    /// # Errors
    /// * If the DWFA is already finalized
    fn increase_edit_distance(&mut self, baseline_seq: &[u8], other_seq: &[u8]) -> Result<(), Box<dyn std::error::Error>> {
        match self.state {
            DWFALiteState::Active => {},
            DWFALiteState::Finalized => {
                bail!("Cannot increase edit distance after finalizing a DWFA");
            },
            DWFALiteState::ExceededEditDistanceLimit => {
                bail!("Cannot increase edit distance after exceeding the edit distance limit");
            },
        };

        // first, increase the distance we're at
        self.edit_distance += 1;

        // grow and re-extend; the bitmask only pays off for narrow wavefronts, so wide ones skip that work
        let (edit_distance, offset, wildcard) = (self.edit_distance, self.offset, self.config.wildcard);
        let tracked = self.wavefront.len() + 2 <= FRONTIER_BITS;
        let grown = if tracked {
            grow_wavefront(&mut self.wavefront, WavefrontGrower::new(baseline_seq, other_seq, edit_distance, offset, wildcard))
        } else {
            grow_wide_wavefront(&mut self.wavefront, baseline_seq, other_seq, edit_distance, offset, wildcard)
        };
        self.max_wavefront = grown.max_wavefront;
        self.max_baseline = grown.max_baseline;
        self.frontier = tracked.then_some(grown.frontier);
        self.extended_len = other_seq.len();
        Ok(())
    }

    /// Applies all remaining edits of `finalize()` at once. Call it only when no diagonal can extend again, i.e.,
    /// each diagonal is at or past the end of `other_seq`. Then each edit raises both cached maxima by exactly one,
    /// so the number of edits up to the baseline end or the cap is known, and `grow_without_matches()` has a closed
    /// form for the wavefront. The result is the same as from that many calls of `increase_edit_distance()`.
    /// # Arguments
    /// * `baseline_len` - the length of the baseline sequence
    fn finish_without_matches(&mut self, baseline_len: usize) {
        let mut steps = baseline_len.saturating_sub(self.max_baseline);
        if let Some(limit) = self.config.max_edit_distance {
            steps = steps.min(limit.saturating_sub(self.edit_distance));
        }
        if steps == 0 {
            return;
        }
        grow_without_matches(&mut self.wavefront, self.edit_distance, steps);
        self.edit_distance += steps;
        self.max_wavefront += steps;
        self.max_baseline += steps;
        // this only runs for wavefronts that are too wide for the frontier bitmask
        self.frontier = None;
    }

    /// This function signals that base insertion into `other_seq` is completed.
    /// This will trigger the algorithm to make sure we have reached the end of both the `baseline_seq`, potentially increasing edit distance further.
    /// Note: other_seq should already be at the end UNLESS we allow_early_termination
    /// # Arguments
    /// * `baseline_seq` - the baseline sequence, theoretically fixed
    /// * `other_seq` - the other sequence, typically getting updates
    /// # Errors
    /// * If the edit distance cannot be increased further and it needs to be.
    pub fn finalize(&mut self, baseline_seq: &[u8], other_seq: &[u8]) -> Result<(), Box<dyn std::error::Error>> {
        match self.state {
            DWFALiteState::Finalized => {
                bail!("Cannot finalize a DWFA twice.");
            },
            DWFALiteState::ExceededEditDistanceLimit => return Ok(()),
            DWFALiteState::Active => {},
        };
        while self.maximum_baseline_distance() < baseline_seq.len() {
            if self.edit_distance_at_limit() {
                self.state = DWFALiteState::ExceededEditDistanceLimit;
                return Ok(());
            }
            // Once every diagonal is at or past the end of `other_seq`, no diagonal can extend again, and the remaining
            // edits have a closed form. Values only grow, so this stays true. The scan usually stops at once, on a low
            // diagonal that has not reached the end. In the worst case, it adds one pass per edit to the growth pass.
            if self.wavefront.len() > FRONTIER_BITS && self.wavefront.iter().all(|&d| d + self.offset >= other_seq.len()) {
                self.finish_without_matches(baseline_seq.len());
            } else {
                // while we have not reached the end of the primary sequence
                self.increase_edit_distance(baseline_seq, other_seq)?;
            }
        }
        self.state = DWFALiteState::Finalized;
        Ok(())
    }

    /// True when the next edit would pass the cap resolved at construction.
    fn edit_distance_at_limit(&self) -> bool {
        self.config.max_edit_distance.is_some_and(|limit| self.edit_distance >= limit)
    }

    /// Helper function that will determine the farthest distance reached into the `baseline_seq` so far.
    pub fn maximum_baseline_distance(&self) -> usize {
        // this is the maximum of `d + edit_distance - i` over the wavefront, cached by `extend()`
        self.max_baseline
    }

    /// Helper function that will determine the farthest distance reached into the `other_seq` so far.
    /// After a public function call, this should be the `other_seq` length unless the edit-distance limit was reached.
    pub fn maximum_other_distance(&self) -> usize {
        // other distance is directly tracked in our wavefront, cached by `extend()`
        self.offset + self.max_wavefront
    }

    /// Returns true if the farther wavefront in the baseline is at the end
    /// # Arguments
    /// * `baseline_seq` - the baseline sequence, theoretically fixed
    pub fn reached_baseline_end(&self, baseline_seq: &[u8]) -> bool {
        self.maximum_baseline_distance() == baseline_seq.len()
    }

    /// This will return the set of candidate extensions that do not require increasing the edit distance.
    /// It also includes how many times that character was counted in the event of multiple possible extension points.
    /// Only active alignments nominate candidates.
    /// # Arguments
    /// * `baseline_seq` - the baseline sequence, theoretically fixed
    /// * `other_seq` - the other sequence, typically getting updates
    #[allow(dead_code)]
    #[deprecated(note = "use `extension_symbols` instead")]
    pub(super) fn get_extension_candidates(&self, baseline_seq: &[u8], other_seq: &[u8]) -> HashMap<u8, usize> {
        if self.state != DWFALiteState::Active {
            return Default::default();
        }
        let mut ret: HashMap<u8, usize> = Default::default();
        for (i, &d) in self.wavefront.iter().enumerate() {
            let other_offset = d + self.offset;
            if other_offset == other_seq.len() {
                let offset = d + self.edit_distance - i;
                if offset < baseline_seq.len() {
                    // ret.insert(baseline_seq[offset]);
                    let entry = ret.entry(baseline_seq[offset]).or_insert(0);
                    *entry += 1;
                }
            }
        }
        ret
    }

    /// Appends the symbols that extend this DWFA without increasing the edit distance to `symbols`.
    /// The order is the wavefront order, which is also the insertion order of `get_extension_candidates()`.
    /// A symbol is repeated once for each diagonal that nominates it.
    /// Inactive alignments leave the output buffer unchanged.
    /// # Arguments
    /// * `baseline_seq` - the baseline sequence, theoretically fixed
    /// * `other_seq` - the other sequence, typically getting updates
    /// * `symbols` - output buffer; existing contents are kept
    pub(crate) fn extension_symbols(&self, baseline_seq: &[u8], other_seq: &[u8], symbols: &mut Vec<u8>) {
        if self.state != DWFALiteState::Active {
            return;
        }
        match self.frontier {
            // the bitmask holds exactly the diagonals with `d + offset == other_seq.len()`, in wavefront order
            Some(mut bits) if other_seq.len() == self.extended_len => {
                while bits != 0 {
                    let i = bits.trailing_zeros() as usize;
                    bits &= bits - 1;
                    let offset = self.wavefront[i] + self.edit_distance - i;
                    if offset < baseline_seq.len() {
                        symbols.push(baseline_seq[offset]);
                    }
                }
            },
            _ => {
                for (i, &d) in self.wavefront.iter().enumerate() {
                    if d + self.offset == other_seq.len() {
                        let offset = d + self.edit_distance - i;
                        if offset < baseline_seq.len() {
                            symbols.push(baseline_seq[offset]);
                        }
                    }
                }
            },
        };
    }

    /// Adds the vote of this DWFA to `votes`, like `votes.add()` with the output of `extension_symbols()`.
    /// Most DWFAs have one diagonal at the frontier, and those vote without the `symbols` buffer.
    /// Inactive alignments leave existing votes unchanged.
    /// # Arguments
    /// * `baseline_seq` - the baseline sequence, theoretically fixed
    /// * `other_seq` - the other sequence, typically getting updates
    /// * `weight` - the total weight of this DWFA's vote
    /// * `votes` - the votes of all DWFAs so far
    /// * `scratch` - buffers for a DWFA with more than one frontier diagonal
    #[inline(always)]
    pub(crate) fn vote(&self, baseline_seq: &[u8], other_seq: &[u8], weight: f64, votes: &mut ExtensionVotes, scratch: &mut VoteScratch) {
        if self.state != DWFALiteState::Active {
            return;
        }
        // try to short-circuit the extension symbols computation
        if let Some(bits) = self.frontier {
            let low_bits = bits as u64;
            if other_seq.len() == self.extended_len && (bits >> 64) == 0 && low_bits & low_bits.wrapping_sub(1) == 0 {
                // zero or one diagonal at the frontier, in the low 64 bits
                if low_bits != 0 {
                    let i = low_bits.trailing_zeros() as usize;
                    let offset = self.wavefront[i] + self.edit_distance - i;
                    if offset < baseline_seq.len() {
                        votes.add_single(baseline_seq[offset], weight);
                    }
                }
                return;
            }
        }

        // otherwise, compute the full set of extension symbols
        scratch.symbols.clear();
        self.extension_symbols(baseline_seq, other_seq, &mut scratch.symbols);
        votes.add(&scratch.symbols, weight, &mut scratch.counts);
    }

    // Getters below
    pub fn edit_distance(&self) -> usize {
        self.edit_distance
    }

    pub fn wavefront(&self) -> &[usize] {
        &self.wavefront
    }

    pub fn max_edit_distance(&self) -> Option<usize> {
        self.config.max_edit_distance
    }

    pub fn state(&self) -> DWFALiteState {
        self.state
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::dwfa_config::DWFALiteConfigBuilder;
    use super::super::wavefront::FRONTIER_BITS;

    /// `finalize()` with the closed form against growth one edit at a time, with caps and offsets
    #[test]
    fn test_bulk_finalize_matches_individual_growth() {
        for cap in [None, Some(63), Some(64), Some(65), Some(200), Some(800)] {
            for offset in [0, 7] {
                let baseline = vec![b'A'; 800];
                let other = vec![b'C'; offset + 30];
                let mut actual = DWFALite::new(dwfa_config(None, false, cap));
                actual.set_offset(offset);
                actual.update(&baseline, &other).unwrap();
                let mut expected = actual.clone();
                while expected.max_baseline < baseline.len() && !expected.edit_distance_at_limit() {
                    expected.increase_edit_distance(&baseline, &other).unwrap();
                }
                expected.state = if expected.max_baseline < baseline.len() {
                    DWFALiteState::ExceededEditDistanceLimit
                } else {
                    DWFALiteState::Finalized
                };
                actual.finalize(&baseline, &other).unwrap();
                assert_eq!(actual, expected, "cap {cap:?}, offset {offset}");
            }
        }
    }

    /// Edit-distance DP table of `baseline` (rows) against `other` (columns); a baseline `wildcard` matches anything.
    /// # Arguments
    /// * `baseline` - the baseline sequence, one row per symbol
    /// * `other` - the other sequence, one column per symbol
    /// * `wildcard` - an optional wildcard symbol; a baseline symbol equal to it matches anything
    fn reference_table(baseline: &[u8], other: &[u8], wildcard: Option<u8>) -> Vec<Vec<usize>> {
        let mut table = vec![vec![0; other.len() + 1]; baseline.len() + 1];
        for (k, row) in table.iter_mut().enumerate() {
            row[0] = k;
        }
        for (m, cell) in table[0].iter_mut().enumerate() {
            *cell = m;
        }
        for k in 1..=baseline.len() {
            for m in 1..=other.len() {
                let is_match = baseline[k-1] == other[m-1] || wildcard == Some(baseline[k-1]);
                let diagonal = table[k-1][m-1] + usize::from(!is_match);
                table[k][m] = diagonal.min(table[k-1][m] + 1).min(table[k][m-1] + 1);
            }
        }
        table
    }

    #[test]
    fn test_random_against_reference() {
        use rand::{Rng, SeedableRng};
        let mut rng = rand::rngs::StdRng::seed_from_u64(7);
        let mut max_seen = 0;
        for trial in 0..300 {
            let wildcard = if trial % 3 == 0 { Some(b'*') } else { None };
            let alphabet: &[u8] = if wildcard.is_some() { b"ACGT*" } else { b"ACGT" };
            let baseline: Vec<u8> = (0..rng.gen_range(0..200)).map(|_| alphabet[rng.gen_range(0..alphabet.len())]).collect();
            // related pairs keep the distance small; unrelated pairs push it past the frontier bitmask width
            let mut other: Vec<u8> = if trial % 2 == 0 {
                baseline.iter().filter(|_| rng.gen_range(0..10) != 0).cloned().collect()
            } else {
                (0..rng.gen_range(0..200)).map(|_| b"ACGT"[rng.gen_range(0..4)]).collect()
            };
            if other.is_empty() {
                other.push(b'A');
            }
            let offset = if trial % 5 == 0 { rng.gen_range(0..other.len()) } else { 0 };

            let table = reference_table(&baseline, &other[offset..], wildcard);
            let mut dwfa = DWFALite::new(dwfa_config(wildcard, false, None));
            if offset > 0 {
                dwfa.set_offset(offset);
            }
            // the first update after an offset sees the whole prefix at once, like `activate_sequence`
            let first_len = if offset > 0 { offset + 1 } else { 1 };
            for l in first_len..=other.len() {
                let expected = (0..=baseline.len()).map(|k| table[k][l - offset]).min().unwrap();
                assert_eq!(dwfa.update(&baseline, &other[..l]).unwrap(), expected, "trial {trial}, length {l}");
                // clones must behave the same as the original
                if l % 7 == 0 {
                    dwfa = dwfa.clone();
                }
            }
            // `finalize` stops when any diagonal reaches the baseline end, also before the end of `other`.
            // So it returns the best cost to the baseline end over every prefix of `other`, but never less than before.
            let update_ed = dwfa.edit_distance();
            dwfa.finalize(&baseline, &other).unwrap();
            let best_to_baseline_end = *table[baseline.len()].iter().min().unwrap();
            assert_eq!(dwfa.edit_distance(), update_ed.max(best_to_baseline_end), "trial {trial}, finalized");
            max_seen = max_seen.max(dwfa.edit_distance());
        }
        // make sure the wide-wavefront fallback was exercised
        assert!(2 * max_seen + 1 > FRONTIER_BITS);
    }

    /// The wavefront algorithm before the frontier and growth changes: grow by scattering into a zeroed wavefront,
    /// then extend every diagonal one symbol at a time.
    /// # Arguments
    /// * `baseline` - the baseline sequence, theoretically fixed
    /// * `other` - the other sequence
    /// * `wildcard` - an optional wildcard symbol; a baseline symbol equal to it matches anything
    fn reference_wavefront(baseline: &[u8], other: &[u8], wildcard: Option<u8>) -> (usize, Vec<usize>, usize) {
        let extend = |wavefront: &mut Vec<usize>, ed: usize| {
            for (i, d) in wavefront.iter_mut().enumerate() {
                while *d + ed - i < baseline.len() && *d < other.len() &&
                    (baseline[*d + ed - i] == other[*d] || wildcard == Some(baseline[*d + ed - i])) {
                    *d += 1;
                }
            }
        };
        let grow = |wavefront: &mut Vec<usize>| {
            let mut new_wavefront = vec![0; wavefront.len() + 2];
            for (i, &d) in wavefront.iter().enumerate() {
                new_wavefront[i] = new_wavefront[i].max(d);
                new_wavefront[i + 1] = new_wavefront[i + 1].max(d + 1);
                new_wavefront[i + 2] = new_wavefront[i + 2].max(d + 1);
            }
            *wavefront = new_wavefront;
        };
        let max_baseline = |wavefront: &[usize], ed: usize| wavefront.iter().enumerate().map(|(i, &d)| d + ed - i).max().unwrap();

        let (mut ed, mut wavefront) = (0, vec![0]);
        extend(&mut wavefront, ed);
        while *wavefront.iter().max().unwrap() < other.len() {
            ed += 1;
            grow(&mut wavefront);
            extend(&mut wavefront, ed);
        }
        let update_ed = ed;
        while max_baseline(&wavefront, ed) < baseline.len() {
            ed += 1;
            grow(&mut wavefront);
            extend(&mut wavefront, ed);
        }
        (update_ed, wavefront, ed)
    }

    /// Like the `finalize` benchmark: `other` is `baseline` without one stretch, with a few errors.
    /// Large gaps give wide wavefronts. The whole `other` is pushed at once, then `finalize` adds the gap.
    #[test]
    fn test_gapped_pairs_against_reference() {
        use rand::{Rng, SeedableRng};
        let mut rng = rand::rngs::StdRng::seed_from_u64(13);
        let mut max_seen = 0;
        for trial in 0..60 {
            let wildcard = if trial % 4 == 0 { Some(b'*') } else { None };
            let alphabet: &[u8] = if wildcard.is_some() { b"ACGT*" } else { b"ACGT" };
            let len = rng.gen_range(2..600);
            let baseline: Vec<u8> = (0..len).map(|_| alphabet[rng.gen_range(0..alphabet.len())]).collect();
            let gap = rng.gen_range(1..len);
            let start = rng.gen_range(0..=len - gap);
            let mut other: Vec<u8> = baseline[..start].iter().chain(baseline[start + gap..].iter()).cloned().collect();
            for symbol in other.iter_mut() {
                if rng.gen_range(0..50) == 0 {
                    *symbol = b"ACGT"[rng.gen_range(0..4)];
                }
            }
            if other.is_empty() {
                other.push(b'A');
            }
            // also swap the roles, so that `other` is the longer sequence
            let (baseline, other) = if trial % 2 == 0 { (baseline, other) } else { (other, baseline) };

            let (expected_update, expected_wavefront, expected_final) = reference_wavefront(&baseline, &other, wildcard);
            let table = reference_table(&baseline, &other, wildcard);
            assert_eq!(expected_update, (0..=baseline.len()).map(|k| table[k][other.len()]).min().unwrap(), "trial {trial}");

            let mut dwfa = DWFALite::new(dwfa_config(wildcard, false, None));
            assert_eq!(dwfa.update(&baseline, &other).unwrap(), expected_update, "trial {trial}, update");
            dwfa.finalize(&baseline, &other).unwrap();
            assert_eq!(dwfa.edit_distance(), expected_final, "trial {trial}, finalized");
            assert_eq!(dwfa.wavefront(), &expected_wavefront[..], "trial {trial}, wavefront");
            assert_eq!(dwfa.maximum_other_distance(), *expected_wavefront.iter().max().unwrap(), "trial {trial}");
            max_seen = max_seen.max(dwfa.edit_distance());
        }
        // make sure that both the narrow and the wide growth were exercised
        assert!(2 * max_seen + 1 > 4 * FRONTIER_BITS);
    }

    #[test]
    #[allow(deprecated)]
    fn test_extension_symbols_match_candidates() {
        let baseline = b"ACGTTACGGA";
        let other = b"ACTTAC";
        let mut dwfa = DWFALite::default();
        for l in 0..other.len() {
            dwfa.update(baseline, &other[..(l+1)]).unwrap();
        }
        let mut symbols = vec![];
        dwfa.extension_symbols(baseline, other, &mut symbols);
        let mut counts: HashMap<u8, usize> = Default::default();
        for &s in symbols.iter() {
            *counts.entry(s).or_insert(0) += 1;
        }
        assert_eq!(counts, dwfa.get_extension_candidates(baseline, other));
    }

    #[test]
    #[allow(deprecated)]
    fn test_finalized_alignments_do_not_nominate_candidates() {
        let baseline = b"AC";
        let other = b"C";
        let mut dwfa = DWFALite::default();
        dwfa.update(baseline, other).unwrap();
        assert_eq!(dwfa.get_extension_candidates(baseline, other), HashMap::from_iter([(b'A', 1), (b'C', 1)]));

        dwfa.finalize(baseline, other).unwrap();
        assert_eq!(dwfa.state(), DWFALiteState::Finalized);
        assert!(dwfa.get_extension_candidates(baseline, other).is_empty());
    }

    #[test]
    fn test_finalized_alignments_preserve_symbol_buffer() {
        let baseline = b"AC";
        let other = b"C";
        let mut dwfa = DWFALite::default();
        dwfa.update(baseline, other).unwrap();
        let mut symbols = vec![0, 255];
        dwfa.extension_symbols(baseline, other, &mut symbols);
        assert_eq!(symbols, vec![0, 255, b'C', b'A']);

        dwfa.finalize(baseline, other).unwrap();
        assert_eq!(dwfa.state(), DWFALiteState::Finalized);
        dwfa.extension_symbols(baseline, other, &mut symbols);
        assert_eq!(symbols, vec![0, 255, b'C', b'A']);
    }

    #[test]
    fn test_finalized_alignments_preserve_existing_votes() {
        let baseline = b"AC";
        let other = b"CA";
        let mut dwfa = DWFALite::new(dwfa_config(None, true, None));
        dwfa.update(baseline, other).unwrap();
        let mut votes = ExtensionVotes::default();
        votes.add_single(b'G', 0.25);
        let mut scratch = VoteScratch::default();
        dwfa.vote(baseline, other, 0.5, &mut votes, &mut scratch);
        let mut candidates = vec![];
        let expected = HashMap::from_iter([(b'G', 0.25), (b'C', 0.5)]);
        votes.candidates(None, &mut candidates);
        assert_eq!(HashMap::from_iter(candidates.iter().copied()), expected);

        dwfa.finalize(baseline, other).unwrap();
        assert_eq!(dwfa.state(), DWFALiteState::Finalized);
        dwfa.vote(baseline, other, 0.5, &mut votes, &mut scratch);
        votes.candidates(None, &mut candidates);
        assert_eq!(HashMap::from_iter(candidates.iter().copied()), expected);
    }

    fn dwfa_config(wildcard: Option<u8>, allow_early_termination: bool, max_edit_distance: Option<usize>) -> DWFALiteConfig {
        DWFALiteConfigBuilder::default()
            .wildcard(wildcard)
            .allow_early_termination(allow_early_termination)
            .max_edit_distance(max_edit_distance)
            .build()
            .unwrap()
    }
}
