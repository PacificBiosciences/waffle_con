
use rustc_hash::FxHashMap as HashMap;
use simple_error::bail;

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

/// Fixed options for a [`DWFALite`].
/// These values are chosen when the DWFA is constructed and do not change as it is extended.
///
/// ```
/// use waffle_con::dynamic_wfa::{DWFALite, DWFALiteConfigBuilder};
/// let config = DWFALiteConfigBuilder::default()
///     .wildcard(Some(b'N'))
///     .allow_early_termination(true)
///     .max_edit_distance(Some(500))
///     .build()
///     .unwrap();
/// let dwfa = DWFALite::new(config);
/// assert_eq!(dwfa.max_edit_distance(), Some(500));
/// ```
#[derive(derive_builder::Builder, Clone, Debug, Eq, Hash, PartialEq)]
#[builder(default)]
pub struct DWFALiteConfig {
    /// Optional wildcard symbol that matches anything.
    pub wildcard: Option<u8>,
    /// When true, reaching the end of the baseline does not penalize a longer other sequence.
    pub allow_early_termination: bool,
    /// Absolute edit-distance cap. `None` means unlimited.
    pub max_edit_distance: Option<usize>,
}

// explicitly set our default for now, even if it's derivable
#[allow(clippy::derivable_impls)]
impl Default for DWFALiteConfig {
    fn default() -> Self {
        DWFALiteConfig {
            wildcard: None,
            allow_early_termination: false,
            max_edit_distance: None,
        }
    }
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

/// The widest wavefront that the `frontier` bitmask can track.
const FRONTIER_BITS: usize = u128::BITS as usize;

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
        match self.advance_single_frontier(baseline_seq, other_seq) {
            Some(_) => Ok(self.edit_distance),
            None => self.update_general(baseline_seq, other_seq),
        }
    }

    /// The most common case of `update()`, without the general code. It applies when one symbol was appended,
    /// exactly one diagonal is at the frontier, and that diagonal matches the new symbol. Then the diagonal moves
    /// one step, and nothing else changes. It returns the new baseline offset of the diagonal.
    /// In all other cases, this returns `None` and changes nothing.
    /// For speed, it only handles a frontier diagonal in the low 64 bits of the bitmask.
    #[inline(always)]
    fn advance_single_frontier(&mut self, baseline_seq: &[u8], other_seq: &[u8]) -> Option<usize> {
        let bits = self.frontier?;
        let low_bits = bits as u64;
        // the last update must have reached the end of `other_seq`, so that the frontier diagonal holds the maximum
        if self.state != DWFALiteState::Active || (bits >> 64) != 0 || !low_bits.is_power_of_two() ||
            other_seq.len() != self.extended_len + 1 || self.offset + self.max_wavefront != self.extended_len {
            return None;
        }
        let i = low_bits.trailing_zeros() as usize;
        // the frontier diagonal is at `extended_len`, so it holds the cached maximum; this saves a dependent load
        let d = self.max_wavefront;
        debug_assert_eq!(d, self.wavefront[i]);
        let baseline_offset = d + self.edit_distance - i;
        let symbol = *baseline_seq.get(baseline_offset)?;
        if symbol != other_seq[d + self.offset] && self.config.wildcard != Some(symbol) {
            return None;
        }
        self.wavefront[i] = d + 1;
        self.max_wavefront = d + 1;
        self.max_baseline = self.max_baseline.max(baseline_offset + 1);
        self.extended_len = other_seq.len();
        Some(baseline_offset + 1)
    }

    /// `update()` for all cases.
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
                    let i = bits.trailing_zeros() as usize;
                    bits &= bits - 1;
                    let d = &mut self.wavefront[i];
                    *d += match_length(baseline_seq, other_seq, *d + edit_distance - i, *d + offset, wildcard);
                    if *d + offset == other_len {
                        new_bits |= 1 << i;
                    }
                    self.max_wavefront = self.max_wavefront.max(*d);
                    self.max_baseline = self.max_baseline.max(*d + edit_distance - i);
                }
                self.frontier = Some(new_bits);
            },
            _ => {
                let frontier = if other_len < self.extended_len { 0 } else { self.extended_len };
                let mut max_wavefront = 0;
                let mut max_baseline = 0;
                let mut new_bits: u128 = 0;
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
            // while we have not reached the end of the primary sequence
            self.increase_edit_distance(baseline_seq, other_seq)?;
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
    /// # Arguments
    /// * `baseline_seq` - the baseline sequence, theoretically fixed
    /// * `other_seq` - the other sequence, typically getting updates
    pub fn get_extension_candidates(&self, baseline_seq: &[u8], other_seq: &[u8]) -> HashMap<u8, usize> {
        if self.state == DWFALiteState::ExceededEditDistanceLimit {
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
    /// # Arguments
    /// * `baseline_seq` - the baseline sequence, theoretically fixed
    /// * `other_seq` - the other sequence, typically getting updates
    /// * `symbols` - output buffer; existing contents are kept
    pub(crate) fn extension_symbols(&self, baseline_seq: &[u8], other_seq: &[u8], symbols: &mut Vec<u8>) {
        if self.state == DWFALiteState::ExceededEditDistanceLimit {
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

/// Extension votes from several DWFAs.
/// It gives the same map as voting into a `HashMap<u8, f64>` with `entry()`, but it does one map insertion for each
/// symbol instead of one map lookup for each DWFA.
#[derive(Debug, Default)]
pub(crate) struct ExtensionVotes {
    /// Votes in order of the first vote for each symbol. The alphabet is small, so a linear search is fast.
    votes: Vec<(u8, f64)>,
}

impl ExtensionVotes {
    /// Adds the votes of one DWFA. Each symbol gets `weight * count / symbols.len()`.
    /// The symbols get their first votes in the same order as when voting from the map of
    /// `DWFALite::get_extension_candidates()`. Usually all symbols are equal, and then no map is built.
    /// # Arguments
    /// * `symbols` - the output of `DWFALite::extension_symbols()` for one DWFA
    /// * `weight` - the total weight of this DWFA's vote
    pub(crate) fn add(&mut self, symbols: &[u8], weight: f64) {
        let Some(&first) = symbols.first() else {
            return;
        };
        let vote_split = symbols.len() as f64;
        if symbols.iter().all(|&s| s == first) {
            self.add_vote(first, weight * symbols.len() as f64 / vote_split);
        } else {
            // Colliding keys make the iteration order of the final map depend on insertion order.
            // Build the same per-DWFA map as before, so the symbols get their first votes in the same order.
            let mut counts: HashMap<u8, usize> = Default::default();
            for &s in symbols.iter() {
                *counts.entry(s).or_insert(0) += 1;
            }
            for (&c, &occ) in counts.iter() {
                self.add_vote(c, weight * occ as f64 / vote_split);
            }
        }
    }

    fn add_vote(&mut self, symbol: u8, vote: f64) {
        match self.votes.iter_mut().find(|(s, _)| *s == symbol) {
            Some((_, total)) => *total += vote,
            None => self.votes.push((symbol, 0.0 + vote)),
        };
    }

    /// Converts the votes into a map. Keys are inserted in order of their first vote, as `entry()` would do.
    pub(crate) fn into_map(self) -> HashMap<u8, f64> {
        let mut candidates: HashMap<u8, f64> = Default::default();
        for (symbol, total) in self.votes.into_iter() {
            candidates.insert(symbol, total);
        }
        candidates
    }
}

/// The cached values after a wavefront grows.
struct GrownWavefront {
    max_wavefront: usize,
    max_baseline: usize,
    /// Bit `i` is set when diagonal `i` is at the end of `other_seq`. Only filled when the frontier is tracked.
    frontier: u128,
}

impl GrownWavefront {
    fn new() -> GrownWavefront {
        GrownWavefront { max_wavefront: 0, max_baseline: 0, frontier: 0 }
    }
}

// Grows a wavefront by two diagonals for an edit distance that is one larger, then extends each diagonal.
// Each new diagonal sits at a position that was never compared, so all of them are visited.
// Each new value is the maximum of up to three old values:
// * deletion (skipping) of a base in `baseline_seq` from diagonal `i`; this does not change the distance into `other_seq`
// * mismatch from diagonal `i - 1`; this progresses both sequences
// * insertion of a base into `baseline_seq` from diagonal `i - 2`; this progresses `other_seq`
// Diagonal 0 has only a deletion, diagonal 1 has no insertion, and the two new diagonals at the end have no deletion.

/// Grows a wavefront that fits the frontier bitmask, in place, and fills `GrownWavefront::frontier`.
/// It is filled from the end, so that `wavefront[i - 2..=i]` still holds the old values when diagonal `i` is computed.
/// The `Vec` grows with amortized doubling, and a clone has spare capacity, so most calls do not allocate.
#[inline(always)]
fn grow_wavefront(wavefront: &mut Vec<usize>, mut grower: WavefrontGrower<'_>) -> GrownWavefront {
    let old_len = wavefront.len();
    debug_assert!(old_len > 0 && old_len + 2 <= FRONTIER_BITS);
    wavefront.resize(old_len + 2, 0);
    let wavefront = &mut wavefront[..old_len + 2];
    wavefront[old_len + 1] = grower.visit(old_len + 1, wavefront[old_len - 1] + 1);
    let insertion = if old_len >= 2 { wavefront[old_len - 2] + 1 } else { 0 };
    wavefront[old_len] = grower.visit(old_len, (wavefront[old_len - 1] + 1).max(insertion));
    if old_len >= 2 {
        // the old values of `i` and `i - 1` carry over from the last step
        let mut deletion = wavefront[old_len - 1];
        let mut mismatch = wavefront[old_len - 2];
        for i in (2..old_len).rev() {
            let insertion = wavefront[i - 2];
            wavefront[i] = grower.visit(i, deletion.max(mismatch + 1).max(insertion + 1));
            deletion = mismatch;
            mismatch = insertion;
        }
        wavefront[1] = grower.visit(1, wavefront[1].max(wavefront[0] + 1));
    }
    wavefront[0] = grower.visit(0, wavefront[0]);
    grower.grown
}

/// Grows a wide wavefront (without the frontier bitmask) in two passes over a new allocation. The result is the
/// same as `grow_wavefront()`. Wide wavefronts come from large edit distances, e.g., in `finalize()`, where most
/// diagonals are past the end of a sequence and can never extend again. A separate pass lets the gather use vector
/// instructions, and the extension pass skips those diagonals with a well-predicted branch.
fn grow_wide_wavefront(
    wavefront: &mut Vec<usize>, baseline_seq: &[u8], other_seq: &[u8], edit_distance: usize, offset: usize, wildcard: Option<u8>
) -> GrownWavefront {
    let old = &wavefront[..];
    let old_len = old.len();
    debug_assert!(old_len >= 2);
    let mut new_wavefront: Vec<usize> = Vec::with_capacity(old_len + 2);

    // gather, see `grow_wavefront()` for the three moves
    new_wavefront.push(old[0]);
    new_wavefront.push(old[1].max(old[0] + 1));
    new_wavefront.extend(old.windows(3).map(|w| w[2].max(w[1] + 1).max(w[0] + 1)));
    new_wavefront.push((old[old_len - 1] + 1).max(old[old_len - 2] + 1));
    new_wavefront.push(old[old_len - 1] + 1);

    // extend and get the maxima; `match_length()` is 0 for a diagonal at or past the end of either sequence
    let (baseline_len, other_len) = (baseline_seq.len(), other_seq.len());
    let mut grown = GrownWavefront::new();
    for (i, d) in new_wavefront.iter_mut().enumerate() {
        if *d + offset < other_len && *d + edit_distance - i < baseline_len {
            *d += match_length(baseline_seq, other_seq, *d + edit_distance - i, *d + offset, wildcard);
        }
        grown.max_wavefront = grown.max_wavefront.max(*d);
        grown.max_baseline = grown.max_baseline.max(*d + edit_distance - i);
    }
    *wavefront = new_wavefront;
    grown
}

/// The inputs that each diagonal needs while a narrow wavefront grows. A method, unlike a closure, can be forced inline.
struct WavefrontGrower<'a> {
    baseline_seq: &'a [u8],
    other_seq: &'a [u8],
    edit_distance: usize,
    offset: usize,
    wildcard: Option<u8>,
    grown: GrownWavefront,
}

impl<'a> WavefrontGrower<'a> {
    #[inline(always)]
    fn new(baseline_seq: &'a [u8], other_seq: &'a [u8], edit_distance: usize, offset: usize, wildcard: Option<u8>) -> Self {
        WavefrontGrower {
            baseline_seq,
            other_seq,
            edit_distance,
            offset,
            wildcard,
            grown: GrownWavefront::new(),
        }
    }

    /// Extends diagonal `i` from `d`, updates the cached values, and returns the extended `d`.
    #[inline(always)]
    fn visit(&mut self, i: usize, d: usize) -> usize {
        let d = d + match_length(self.baseline_seq, self.other_seq, d + self.edit_distance - i, d + self.offset, self.wildcard);
        if d + self.offset == self.other_seq.len() {
            self.grown.frontier |= 1 << i;
        }
        self.grown.max_wavefront = self.grown.max_wavefront.max(d);
        self.grown.max_baseline = self.grown.max_baseline.max(d + self.edit_distance - i);
        d
    }
}

/// Returns how many symbols match from `baseline_seq[baseline_offset..]` and `other_seq[other_offset..]`.
/// A baseline symbol equal to `wildcard` matches anything.
#[inline]
fn match_length(baseline_seq: &[u8], other_seq: &[u8], baseline_offset: usize, other_offset: usize, wildcard: Option<u8>) -> usize {
    let (Some(baseline), Some(other)) = (baseline_seq.get(baseline_offset..), other_seq.get(other_offset..)) else {
        return 0;
    };
    match wildcard {
        Some(wc) => baseline.iter().zip(other.iter())
            .take_while(|&(&b, &o)| b == o || b == wc)
            .count(),
        None => common_prefix_length(baseline, other),
    }
}

/// Returns the length of the common prefix of `a` and `b`, comparing 8 bytes at a time.
#[inline]
fn common_prefix_length(a: &[u8], b: &[u8]) -> usize {
    let n = a.len().min(b.len());
    let mut i = 0;
    while i + 8 <= n {
        let x = u64::from_le_bytes(a[i..i + 8].try_into().unwrap());
        let y = u64::from_le_bytes(b[i..i + 8].try_into().unwrap());
        let diff = x ^ y;
        if diff != 0 {
            // little-endian load, so the lowest set bit is the first mismatched byte
            return i + (diff.trailing_zeros() / 8) as usize;
        }
        i += 8;
    }
    while i < n && a[i] == b[i] {
        i += 1;
    }
    i
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cdwfa_config::CdwfaConfig;

    #[test]
    fn test_common_prefix_length() {
        assert_eq!(common_prefix_length(b"", b"ACGT"), 0);
        assert_eq!(common_prefix_length(b"ACGTACGTACGTA", b"ACGTACGTACGTA"), 13);
        for mismatch in 0..20 {
            let a = vec![b'A'; 20];
            let mut b = a.clone();
            b[mismatch] = b'C';
            assert_eq!(common_prefix_length(&a, &b), mismatch);
            assert_eq!(common_prefix_length(&a[..mismatch], &b), mismatch);
        }
    }

    /// Edit-distance DP table of `baseline` (rows) against `other` (columns); a baseline `wildcard` matches anything.
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
    fn test_extension_votes_match_entry_voting() {
        // `C` and `G` collide in a small FxHashMap, so key order depends on insertion order; compare order and values
        use rand::{Rng, SeedableRng};
        let mut rng = rand::rngs::StdRng::seed_from_u64(11);
        let alphabet = b"ACGT*N";
        for _ in 0..2000 {
            let mut expected: HashMap<u8, f64> = Default::default();
            let mut votes = ExtensionVotes::default();
            for _ in 0..rng.gen_range(0..8) {
                let symbols: Vec<u8> = (0..rng.gen_range(0..4)).map(|_| alphabet[rng.gen_range(0..alphabet.len())]).collect();
                let weight = [1.0, 0.5, 0.25, 0.9][rng.gen_range(0..4)];
                // the voting code before `ExtensionVotes`
                let mut cand: HashMap<u8, usize> = Default::default();
                for &s in symbols.iter() {
                    *cand.entry(s).or_insert(0) += 1;
                }
                let vote_split = cand.values().sum::<usize>() as f64;
                for (&c, &occ) in cand.iter() {
                    *expected.entry(c).or_insert(0.0) += weight * occ as f64 / vote_split;
                }
                votes.add(&symbols, weight);
            }
            let actual = votes.into_map();
            let to_bits = |m: &HashMap<u8, f64>| m.iter().map(|(&k, &v)| (k, v.to_bits())).collect::<Vec<_>>();
            assert_eq!(to_bits(&actual), to_bits(&expected));
        }
    }

    #[test]
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

    fn dwfa_config(wildcard: Option<u8>, allow_early_termination: bool, max_edit_distance: Option<usize>) -> DWFALiteConfig {
        DWFALiteConfigBuilder::default()
            .wildcard(wildcard)
            .allow_early_termination(allow_early_termination)
            .max_edit_distance(max_edit_distance)
            .build()
            .unwrap()
    }

    /// Resolves a consensus cap and stores that absolute maximum on a new DWFA.
    fn dwfa_with_consensus_cap(
        baseline_len: usize,
        hard_max: Option<usize>,
        fraction: Option<f64>,
        floor: Option<usize>,
    ) -> Result<DWFALite, Box<dyn std::error::Error>> {
        let max_edit_distance = CdwfaConfig {
            max_edit_distance: hard_max,
            max_edit_distance_fraction: fraction,
            min_edit_distance: floor,
            ..Default::default()
        }.max_edit_distance_for(baseline_len)?;
        Ok(DWFALite::new(dwfa_config(None, false, max_edit_distance)))
    }

    #[test]
    fn test_new() {
        let dwfa = DWFALite::default();
        assert_eq!(dwfa.edit_distance(), 0);
        assert_eq!(dwfa.wavefront(), &[0]);
    }

    #[test]
    fn test_exact_match() {
        let sequence = b"ACGTACGTACGT";
        let mut other_seq = vec![];
        let mut dwfa = DWFALite::default();
        for &c in sequence.iter() {
            other_seq.push(c);
            assert_eq!(dwfa.update(sequence, &other_seq).unwrap(), 0);
        }
    }

    #[test]
    fn test_simple_mismatch() {
        let sequence =     b"ACGTACGTACGT";
        let alt_sequence = b"ACGTACCTACGT";
        let mut dwfa = DWFALite::default();
        for l in 0..alt_sequence.len() {
            dwfa.update(sequence, &alt_sequence[..(l+1)]).unwrap();
        }
        assert_eq!(dwfa.edit_distance(), 1);
    }

    #[test]
    fn test_simple_insertion() {
        let sequence =     b"ACGTACGTACGT";
        let alt_sequence = b"ACGTACIGTACGT";
        let mut dwfa = DWFALite::default();
        for l in 0..alt_sequence.len() {
            dwfa.update(sequence, &alt_sequence[..(l+1)]).unwrap();
        }
        assert_eq!(dwfa.edit_distance(), 1);
    }

    #[test]
    fn test_simple_deletion() {
        let sequence =     b"ACGTACGTACGT";
        let alt_sequence = b"ACGTACTACGT";
        let mut dwfa = DWFALite::default();
        for l in 0..alt_sequence.len() {
            dwfa.update(sequence, &alt_sequence[..(l+1)]).unwrap();
        }
        assert_eq!(dwfa.edit_distance(), 1);
    }

    #[test]
    fn test_complex_001() {
        let sequence =     b"ACGTACGTACGT";
        let alt_sequence = b"ACTACGCACGGGT";
        let mut dwfa = DWFALite::default();
        for l in 0..alt_sequence.len() {
            dwfa.update(sequence, &alt_sequence[..(l+1)]).unwrap();
        }
        assert_eq!(dwfa.edit_distance(), 4);
    }

    #[test]
    fn test_complex_002() {
        //modified_seq has 2 separate deletions, 1 2bp insertion, and 1 mismatch
        let sequence     = b"AACGGATCAAGCTTACCAGTATTTACGT";
        let alt_sequence = b"AACGGACAAAAGCTTACCTGTATTACGT";
        
        let mut dwfa = DWFALite::default();
        /*
        for l in 0..alt_sequence.len() {
            dwfa.update(sequence, &alt_sequence[..(l+1)]).unwrap();
        }
        */
        // test a shortcut
        dwfa.update(sequence, alt_sequence).unwrap();
        assert_eq!(dwfa.edit_distance(), 5);
    }

    #[test]
    fn test_big_insertion() {
        // one big insertion in the middle
        let sequence     = b"AACGGATTTTACGT";
        let alt_sequence = b"AACGGATAAAAGCTTACCTGTTTTACGT";
        
        let mut dwfa = DWFALite::default();
        for l in 0..alt_sequence.len() {
            dwfa.update(sequence, &alt_sequence[..(l+1)]).unwrap();
        }
        assert_eq!(dwfa.edit_distance(), alt_sequence.len() - sequence.len());
    }

    #[test]
    fn test_big_deletion() {
        // one big deletion in the middle
        let sequence     = b"ATTTTTTTTTTAAAAAAAAAA";
        let alt_sequence = b"AAAAAAAAAAA";
        
        let mut dwfa = DWFALite::default();
        for l in 0..alt_sequence.len() {
            dwfa.update(sequence, &alt_sequence[..(l+1)]).unwrap();
        }
        assert_eq!(dwfa.edit_distance(), sequence.len() - alt_sequence.len());
    }

    #[test]
    fn test_required_finalize() {
        // the ALT is a lot smaller than the baseline, so we need to run a finalize to get the full ED
        let sequence     = b"ATTTTTTTTTTA";
        let alt_sequence = b"AA";
        
        let mut dwfa = DWFALite::default();
        for l in 0..alt_sequence.len() {
            dwfa.update(sequence, &alt_sequence[..(l+1)]).unwrap();
        }

        // here it has only compared AT to AA, so ED=1
        assert_eq!(dwfa.edit_distance(), 1);

        // after finalizing, it will do end-to-end comparison
        dwfa.finalize(sequence, alt_sequence).unwrap();
        assert_eq!(dwfa.edit_distance(), sequence.len() - alt_sequence.len());
    }

    #[test]
    fn test_cloning() {
        let sequence     = b"AAAAAAA";
        let alt_sequence = b"AAACAAA";

        let mut dwfa = DWFALite::default();
        let mut dwfa2 = dwfa.clone();
        for l in 0..alt_sequence.len() {
            dwfa.update(sequence, &sequence[..(l+1)]).unwrap();
            dwfa2.update(sequence, &alt_sequence[..(l+1)]).unwrap();

            if sequence[l] == alt_sequence[l] {
                // same sequence still
                assert_eq!(dwfa, dwfa2);
            } else {
                // should have different sequences now
                assert_ne!(dwfa, dwfa2);

                // re-clone the first one
                dwfa2 = dwfa.clone();
            }
        }

        // in the end, both should exactly match due to cloning
        assert_eq!(dwfa.edit_distance(), 0);
        assert_eq!(dwfa2.edit_distance(), 0);
    }

    #[test]
    fn test_wildcards_001() {
        //modified_seq has several wildcards, but otherwise an exact match
        let consensus= b"AACGGATCAAGCTTACCAGTATTTACGT";
        let baseline = b"*ACGGATCAA**TTACCA*TATTTACG*";
        
        let mut dwfa = DWFALite::new(dwfa_config(Some(b'*'), false, None));
        dwfa.update(baseline, consensus).unwrap();
        assert_eq!(dwfa.edit_distance(), 0);
    }

    #[test]
    fn test_wildcards_002() {
        //modified_seq has several wildcards, and we added 1 del, 1 SNV, and 1 ins
        let consensus= b"AACGGATCAAGCTTACCAGTATTTACGT";
        let baseline = b"*ACGATCAA**TATACCA*TATCTACG*";
        
        let mut dwfa = DWFALite::new(dwfa_config(Some(b'*'), false, None));
        dwfa.update(baseline, consensus).unwrap();
        assert_eq!(dwfa.edit_distance(), 3);
    }

    #[test]
    fn test_early_termination_001() {
        let consensus = b"ACGTACGT";
        let baseline = b"ACGT";
        let mut dwfa = DWFALite::new(dwfa_config(None, true, None));
        dwfa.update(baseline, consensus).unwrap();
        assert_eq!(dwfa.edit_distance(), 0);
    }

    #[test]
    fn test_big_early_termination() {
        let c1 = "AGCCCATTCTGGCCCCTTCCCCACATGCCAGGACAATGTAGTCCTTGTCACCAATCTGGGCAGTCAGAGTTGGGTCAGTGGGGGACACGGGATTATGGGCAAGGGTAACTGACATCTGCTCAGCCTCAACGTACCCGTCTCAAATGCGGCCAGGCGGTGGGGTAAGCAGGAATGAGGCAGGGGTGGGGTTGCCCTGAGGAGGATGATCCCAACGAGGGCGTGAGCAGGGGACCCGAGTTGGAACTACCACATTGCTTTATTGTACATTAGAGCCTCTGGCTAGGGAGCAGGCTGGGGACTAGGTACCCCATTCTAGCGGGGCACAGCACAAAGCTCATAGGGGGATGGGGTCACCAGAAAGCTGACGACACGAGAGTGGCTGGGCCGGGGCTGTCCGGCGGCCACGGAGAAGCTGAAGTGCTGCAGCAGGGAGGTGAAGAAGAGGAAGAGCTCCATGCGGGCCAGGGGCTCCCCGAGGCATGCACGGCGGCCTGTGGGGAGGGGAGGGGCGTCAGTGAGCCTGGCTCCTGGGTGATACCCCTGCAAGACTCCACGGAAGGGGACAGGGAGCCGGGCTCCCCACAGGCACCTGCTGAGAAAGGCAGGAAGGCCTCCGGCTTCACAAAGTGGCCCTGGGCATCCAGGAAGTGTTCGGGGTGGAAGCGGAAGGGCTTCTCCCAGACGGCCTCATCCTTCAGCACCGATGACAGGTTGGTGATGAGTGTCGTTCCCTGGGCAGGAGATGCAGGGTGAGAGTGGGGACTGGACTCTAGGATGCTGGGACCCCTGCCACCAAACACACGGGGGACACACACTGCCTGGCACACAGCTGGACTCTGTCAACTAGTCCTGCGCCCGAGAAGCTCCACAGTACCCTCTCCGACCCCACAGCAGGGCGCAGTCACACCTCTCAGAGGCACCCACACTGCCCCCTCTCCCTGCAGGCGCTGGGTCCTCCAACATTCTGGCAGGTCCTGGTTTGTCTCCCCACTAGACGGGGGCTCTGGATGGACAGGCCAGCCCTGCCTATACTCTGGACCCCCCACCCAAGTGGGGACAGTCAGTGTGGTGGCATTGAGGACTAGGTGGCCAGGGTTCCTAGAGTGGGCCCACCTGGCAGTAGCCATGCTGGGGCTATCACCAGGGGCTGGTGCTGAGCTGGGGTGAGGAGGGCGCCAGGCCTACCTTAGGGATGCGGAAGCCCTGTACTTCGATGTCACGGGATGTCATATGGGTCACACCCAGGGGGACGATGTCCCCAAAGCGCTGCACCTCATGAATCACGGCAGTGGTGTAGGGCATGTGAGCCTGGTCACCCATCTCTGGTCGCCGCACCTGCCCTATCACGTCGTCGATCTCCTGTTGGACACGGACTGGACAGACATGCGTCCCCACAATGGGTCAGCACCCAGGGGACACTCTCCTTCCTCCTGTGTTGGAGGAAGTTAGGCTTACAGGAGCCTGGCCACGCCTGTGCTGGAAGCCCCGGGTGTCCCAGCTAAGCCCAGGGGCCCCCAGCTGTACCCTTCCTCCCTCAGTCCCTGCCTTGGGCCCCAGCTGGGCTCACGCTGCACATCCAGGTGTAGGATCATGAGCAGGAGGCCCCAGGCCAGCGTGGTCAAGGTGGTCACCATCCCGGCAAGGAACAGGTTACCCACCACTATGCGCAGGTTCTCATCATTGAAGCTGCTCTCAGGGCTCCCCTTGGCCTGAGCAGGGCCGAGAGGATACTCAGGGGATAGAACGGGGTAGCCCCCAAATGACCTCCAATTCTGCACCTGTCAGCCCAGATGCGGCTCGCCGGGTGATGCACTGGTCCAACCTTTTGCCCAGCCTCCCCTCATTCCTCCTGGGACGTTCAACCCACCACCCTTGCCCCCCACCGTGGCAGCCACTCTCACCTTCTCCTTCTTTGCCAGGAAGGCCTCAGTCAGGTCTCGGGGTGGCTGGGCTGGGTCCCAGGTCATCCTGTGCTCAGTTAGCAGCTCATCCAGCTGGGTCAGGAAAGCCTTTTGGAAGCGTAGGACCTTGCCAGCCAGCGCTGGGATGTGCGGGAGGACGGGGACAGCATTCAGCACCTACACCAGACAGAACCGGGTCTCAATCCTTCCTGTGCTCTGCGTTCATCTGGACCAGTCTCAGGCCCCAGCCATCTCCAGGAAGACCCAGGGCCTGCCTGTCCTTACCACTGACCTCACCAAGTCCCTCCCCAAGTGCCAGCCTCCACCCTCTCTCTCCTTGCCCAGAGGAGAAACCTAAAATCGAAATCTCCAACGTGGACGGGGGTACAGAGTCCTTGGCCTCTCCTGGTGCCCCCTGACCCGGGCACACCTCTCCCACGACCATGTCTGAGATGTCCCCTCCTCCTCCAGGCCCTTCTTACAGTGGGGTCTCCTGGAATGTCCTTTCCCAAACCCATCTACGCAAATCCTGCCCTTCGGAGGCCCCAGTCCAGCCCCGGCACCTCTCAGGAGCTCGCCCTGCAAAGACCCTTGCTCCGCACCTCGCGCAGGAAGCCCGACTCCTCCTTCGATCCCTCCCTGAGCTAGGTCCAGCAGCCTGAGGAAGCGAGGGTCGTCGTACTCGAAGCGGCGCCCGCAGGTGAGGGAGGCGATCACGTTGCTCACGGCTTTGTCCAAGAGACCGTTGGGGCGAAAGGGGCGTCCTGGGGGTGGGAGATGCGGGTAAGGGGTCGCCTTCTCCGTCCCCCGCCTTCCCAGTTCCCGCTGTGTGCCCTTCTGCCCATCACCCACCGGCTTGGTCGGCGAAGGCGGCACAAAGGCAGGCGGCCTCCTCGGTCACCCACTGCTCCAGCGACTTCTTGCCCAGGCCCAAGTTGCGCAAGGTGGACACGGAGAAGCGCCTCTGCTCGCGCCACGCGGGCCCATAGCGCGACAGGATCACCCCTGTGGGCGGGACGGACACGTGGGCGTTGCCATGAAGGCCTTGGCCCCACCCTCCGCCACCCACTCCAACCCTGGCGCTCCACAAGGTCTCCCGCAGTCCCTAGCCCGGTCCAGCTGGGCACAGGGCCCACTCTTTGCTCACCCACATTGCTCCCCTGCCTGGGGCGGGGTTTGGCCCCACCTCGTCTCTGCCCACCCTGACCACCTTTCCACTCAAGGAAGATCCCGCCCGTCCCGCCCACACTGAGCCCGCAGCATAGGCGCGGTCCCCGCCACCGCCACTTCGACGCATCAGCCTCGCCCACCGGGCTTCTGGCGGGTCTGGGCAGTAGCCCCGCCCCCTCCCAGCCCACAGACTCGCACCTCCCCCGTGCAGGTGGTTTCCTGGCCCACTGTCCTCAGCCCACTCGCTGGCCTTTATCTCTGTTTCACGTCCAGGACCCCACGCCCTGTCGGCGCTGCTTGGGCTACGGTCACTGTCCACCCGGGGCCCACGGAAACGCGGTCTCTGTCCCCCACCGCCGCTTGCCTTGGGAACGCGGCCCGAAGCCCAGGACCTGGTAGATGGGCGCAGGCGGGCGGTCGGCCGTGTCCTCGCCGCGGGTCACCATCGCCTCGCGCACGGCCGCCAGCCCATTGAGCACGACCACCGGCGTCCAGGCCAGCTGCAGGCTGAACACGTCCCCGAAGCGGCGCCGCAACTGCAGAGGGAGGGTCAGGGCCTCTTGTCAAGCCAGGATCCCCCCAGACTACAGGTCCTAGTCCTATTTGAACCTTGGACGACCCCCGGGGCTACCAGGAGTGAGCAGGTGGAAGGAGGAGACCCAGCCTCCTGATCCTGGGGCGGGGGTGGGGGTCACACCTTCTGTGATGGAGGAACTCAGTTTGGATGCGTCACCCAGGTATGACCTTGCAAGAGTCACCAAAATTGCCGAGAGGCCCCAGTTAGCATCCCATTCCCAGATGATGGTCCATGCCGGTGAGCAGTGAGGCCCGAGGACCCACAGTGCAAAAGGTTTGAACCGGGTCACTGCACCCCCTTCATCCTCGATTTCGTGATTTAAACGGCACTCAGGACTAACTCATCTTCCATTCCCAAGGCCTTTCCTTCTGGTGTCAGCAGAAGGGACTTTGTACTCCATAACATATGTTGCCCAATGGGCTTGCATGCCCACTGCCAAGTCCAGCTCCACCTCCAGGCCCTTGCCCTACTCTTCCTTGGCCTTTGGAAAATCCAGTCCTTCATGCCATGTATAAATGTCCTTCCCCAGGACGTCCCCCAAACCTGCTTCCCCTTCTCAGCCTGGCTTCTGATCCAGCCTGTGGTTTAACCCACCACCCATGTTTGCTGGTGGTGGGGCATCCTCAGGACCTCTGCCGCCCTCCAGGACCTCCTCCCTCACCTGGTCGAAGCAGTATGGTGTGTTCTGGAAGTCCACATGCAGCAAGGTTGCCCAGCCCGGGCAGTGGCAGGGGACCTGGCGGGTAGCGTGCAGCCCAGCGTTGGTGCCGGTGCATCAGGTCCACCAGGAGCAGGAAGATGGCCACTATCATGGCCAGGGGCACCAGTGCTTCTAGCCCCATGGCTGCCTCACTACCAACTGGGCTCCTCTGGACACACCTGGCACCCCCACCCCACCAGGCACAGAGGACCAGGCAGGACACTCTCAGCACACCGAGCGCGTGACCCTTCCCTTATAAAGGGAGCTGATGATGGCCTTCGCCCTCTGCTGTGAGTGAACCTGCTGTGTTGACTGTGCTGCCAGTGGCAGAGTCAGGCCAGGGTGGGTATGGGCTGCTCCAGAGGTCCTTGCCGCTGCTTCCTGCTCCAGGCCCTTACCCAGGGTAGGGTGGTAGAAAGGCCTGGTCGGAGAAGTCACCCCCTCTCCCCACTCCAAGCTCCCCAAGCCCACACAGGCTTCTGGGATAACCAGGGTCTCAGTGGACCCGGCCATCCACCTCCCAGCTAGGCTCATACACCGTAATGTAGTCACAACCCCTCCTCCAGAACATGGCCTTGCCCTTTCCCTACCCCCACCTGCCCACTCCAGAGTGACCTTCAGCACCCTTATCTGTCACTGGCACTTACCTGGGGCCTTAGAGCTCCTGATGATGAGTGGCATCATGGGCCTGGTCCCTTCACTTCACCTTGCACTCTTGACATGCACAGACGCTATGCACACACCTGATGGTGCACAGATCTCTTGTCCACTCCCAGACACTTGTCCACTTGTTCACACTTGCAGGGACACGATTACACATGCAGAAAATCACCCACACAAAGACAATATTCACACATACACAGACTCACACTGACACTCAGGGCACACATTCTCTCTCACACACACCAGTCACACACACATACAGACCCGGCACCAAGTACCCCACTTCCCAGCCATGCCCAAGGTTTCCTGGATGGGACCTCTCCTGTCCAGAGGCTGCTCCCAGTGAGCCTCAAAGCTGTCACGTGGATCCCAGCTCAGCCCACATTCTGGGCTCTGGCCGGGCCATGGCTTCCTGTTTGCAACAGGGCTGTTCCCAGAGCTCCCAGTTGGTAGCCTGAAGGCCCTTGCCCCAGCCTGTGACAGCATCCTCCAGGGCTGCCTGAGGGTCGTCATTCTCCACTGCTTCCTGGCCTCCATGTTTCTGATTAGAAATCTGGTGGAAACATTATGGAGGATCCTTTATTTAGGATATGTTGCTTTTTTATTTTTATTTTTTCTTTAGACAGGGTCTCACTCTGTTGCCCGGGCCGGAGTGCAGTGGCAGGATCACGGCTCACTGCAATCTCAACATCAAGTGGACCTCCTGCCTCCCAAGTAGCTGGGACTACAGGCACCACCGAGCCCAAATAATTTTTTTTTTGAGACGGAGTTTTGCTCTGTCGCCCAGGTGGGAGTGCAATGATGCGATCTCGGCTCACTGCAACCTCCACCTCCAGGGTTCAAGCGATTCTCCTGCCTCAGCCTCCCAAGTAGCTGGGATTACAGGTGCCCACCACCATGCCTGGCTGATTTTTTGTA";
        let seq_23 = "AGCCCATTCTGGCCCCTTCCCCACATGCCAGGACAATGTAGTCCTTGTCACCAATCTGGGCAGTCAGAGTTGGGTCAGTGGGGGACATGGGATTATGGGCAAGGGTAACTGACATCTGCTCAGCCTCAACGTACCCGTCTCAAATGCGGCCAGGCGGTGGGGTAAGCAGGAATGAGGCAGGGGTGGGGTTGCCCTGAGGAGGATGATCCCAACGAGGGCGTGAGCAGGGGACCCGAGTTGGAACTACCACATTGCTTTATTGTACATTAGAGCCTCTGGCTAGGGAGCAGGCTGGGGACTAGGTACCCCATTCTAGCGGGGCACAGCACAAAGCTCGTAGGGGGATGGGGTCACCAGAAAGCTGACGACACGAGAGTGGCTGGGCCGGGGCTGTCCGGCGGCCACGGAGAAGCTGAAGTGCTGCAGCAGGGAGGTGAAGAAGAGGAAGAGCTCCATGCGGGCCAGGGGCTCCCCGAGGCATGCACGGCGGCCTGTGGGGAGGGGAGGGGCGTCAGTGAGCCTGGCTCCTGGGTGATACCCCTGCAAGACTCCACGGAAGGGGACAGGGAGCCGGGCTCCCCACAGGCACCTGCTGAGAAAGGCAGGAAGGCCTCCGGCTTCACAAAGTGGCCCTGGGCATCCAGGAAGTGT";

        // iterate, making sure everything is fine even as we go well beyond seq_23
        let mut dwfa = DWFALite::new(dwfa_config(None, true, None));
        for i in 0..c1.len() {
            dwfa.update(seq_23.as_bytes(), c1[0..(i+1)].as_bytes()).unwrap();
            assert!(dwfa.edit_distance() <= 2);
        }
        assert_eq!(dwfa.edit_distance(), 2);

        // make sure finalize does not break it
        dwfa.finalize(seq_23.as_bytes(), c1.as_bytes()).unwrap();
        assert_eq!(dwfa.edit_distance(), 2);
        
    }

    #[test]
    fn test_offsets() {
        // this sequence
        let consensus = b"ACGTACGT";
        let baseline =    b"GTACGT";
        let mut dwfa = DWFALite::new(dwfa_config(None, true, None));
        dwfa.set_offset(2);
        dwfa.update(baseline, consensus).unwrap();
        assert_eq!(dwfa.edit_distance(), 0);
    }

    #[test]
    fn test_exact_match_under_edit_distance_cap() {
        let sequence = b"ACGTACGT";
        let mut dwfa = dwfa_with_consensus_cap(sequence.len(), Some(2), None, None).unwrap();
        let other = &sequence[..4];
        assert_eq!(dwfa.update(sequence, other).unwrap(), 0);
        assert_eq!(dwfa.get_extension_candidates(sequence, other).get(&b'A'), Some(&1));
        assert_eq!(dwfa.state(), DWFALiteState::Active);
    }

    #[test]
    fn test_hard_edit_distance_cap() {
        let baseline = b"ACGTACGTACGT";
        let other = b"TTTTTTTTTTTT";
        let mut dwfa = dwfa_with_consensus_cap(baseline.len(), Some(3), None, None).unwrap();
        assert_eq!(dwfa.update(baseline, other).unwrap(), 3);
        assert_eq!(dwfa.state(), DWFALiteState::ExceededEditDistanceLimit);
        assert!(dwfa.get_extension_candidates(baseline, other).is_empty());

        let mut longer = other.to_vec();
        longer.push(b'A');
        assert_eq!(dwfa.update(baseline, &longer).unwrap(), 3);
        dwfa.finalize(baseline, &longer).unwrap();
        assert_eq!(dwfa.edit_distance(), 3);
        assert!(dwfa.get_extension_candidates(baseline, &longer).is_empty());
    }

    #[test]
    fn test_fractional_edit_distance_cap() {
        let baseline = vec![b'A'; 100];
        let mut other = vec![b'T'; 6];
        other.extend(std::iter::repeat(b'A').take(94));
        let mut dwfa = dwfa_with_consensus_cap(baseline.len(), None, Some(0.05), None).unwrap();
        assert_eq!(dwfa.max_edit_distance(), Some(5));
        assert_eq!(dwfa.update(&baseline, &other).unwrap(), 5);
        assert_eq!(dwfa.state(), DWFALiteState::ExceededEditDistanceLimit);
        assert!(dwfa.get_extension_candidates(&baseline, &other).is_empty());
    }

    #[test]
    fn test_edit_distance_cap_uses_tighter_limit() {
        let tighter_hard = dwfa_with_consensus_cap(100, Some(3), Some(0.05), None).unwrap();
        assert_eq!(tighter_hard.max_edit_distance(), Some(3));
        let tighter_fraction = dwfa_with_consensus_cap(100, Some(10), Some(0.05), None).unwrap();
        assert_eq!(tighter_fraction.max_edit_distance(), Some(5));
        let unlimited = dwfa_with_consensus_cap(100, None, None, None).unwrap();
        assert_eq!(unlimited.max_edit_distance(), None);
        assert!(dwfa_with_consensus_cap(100, None, Some(-0.1), None).is_err());
        assert!(dwfa_with_consensus_cap(100, None, Some(f64::NAN), None).is_err());
    }

    #[test]
    fn test_edit_distance_floor_raises_derived_cap() {
        // floor(0.05 * 19) == 0, so the floor is the resolved cap
        let raised = dwfa_with_consensus_cap(19, None, Some(0.05), Some(2)).unwrap();
        assert_eq!(raised.max_edit_distance(), Some(2));
        // a derived cap already above the floor stays at the derived value
        let unchanged = dwfa_with_consensus_cap(100, None, Some(0.05), Some(2)).unwrap();
        assert_eq!(unchanged.max_edit_distance(), Some(5));
        // a floor with no maximum does not create a cap
        let unlimited = dwfa_with_consensus_cap(19, None, None, Some(2)).unwrap();
        assert_eq!(unlimited.max_edit_distance(), None);
    }
}
