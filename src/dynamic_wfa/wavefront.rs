use std::collections::VecDeque;

/// Cached extrema after one growth.
/// `max_wavefront` is the farthest position reached after `offset`; add `offset` for the index in `other_seq`.
/// `max_baseline` is the farthest consumed position in `baseline_seq`.
/// `frontier_start..frontier_end` is the smallest half-open range containing equality lanes.
/// `frontier_start == frontier_end == 0` means none. Lanes past the end are not included.
pub(super) struct GrownWavefront {
    /// Farthest position reached after `offset`. The index in `other_seq` is `offset + max_wavefront`.
    pub(super) max_wavefront: usize,
    /// Farthest consumed position in `baseline_seq`.
    pub(super) max_baseline: usize,
    /// First diagonal at the end of `other_seq`, or `0` when there are none.
    pub(super) frontier_start: usize,
    /// One past the last diagonal at the end of `other_seq`, or `0` when there are none.
    pub(super) frontier_end: usize,
}

impl GrownWavefront {
    /// A growth that has not seen an equality lane yet.
    pub(super) fn new() -> GrownWavefront {
        GrownWavefront { max_wavefront: 0, max_baseline: 0, frontier_start: 0, frontier_end: 0 }
    }

    /// Includes diagonal `i` in the equality-lane bounds.
    fn record_equality_lane(&mut self, i: usize) {
        if self.frontier_end == 0 {
            self.frontier_start = i;
            self.frontier_end = i + 1;
        } else {
            self.frontier_start = self.frontier_start.min(i);
            self.frontier_end = self.frontier_end.max(i + 1);
        }
    }
}

// Grows a wavefront by two diagonals for an edit distance that is one larger, then extends each diagonal.
// Each new diagonal sits at a position that was never compared, so all of them are visited.
// Each new value is the maximum of up to three old values:
// * deletion (skipping) of a base in `baseline_seq` from diagonal `i`; this does not change the distance into `other_seq`
// * mismatch from diagonal `i - 1`; this progresses both sequences
// * insertion of a base into `baseline_seq` from diagonal `i - 2`; this progresses `other_seq`
// Diagonal 0 has only a deletion, diagonal 1 has no insertion, and the two new diagonals at the end have no deletion.

/// Grows a wavefront in place and records equality-lane bounds.
/// It is filled from the end, so that `wavefront[i - 2..=i]` still holds the old values when diagonal `i` is computed.
/// The `Vec` grows with amortized doubling, and a clone has spare capacity, so most calls do not allocate.
/// # Arguments
/// * `wavefront` - the wavefront to grow in place
/// * `grower` - the sequences and edit distance for the new diagonals
#[inline(always)]
pub(super) fn grow_wavefront(wavefront: &mut Vec<usize>, mut grower: WavefrontGrower<'_>) -> GrownWavefront {
    let old_len = wavefront.len();
    debug_assert!(old_len > 0);
    
    // Resize the wavefront to the new length.
    wavefront.resize(old_len + 2, 0);
    let wavefront = &mut wavefront[..old_len + 2];

    // visit the two new diagonals at the end.
    wavefront[old_len + 1] = grower.visit(old_len + 1, wavefront[old_len - 1] + 1);
    let insertion = if old_len >= 2 { wavefront[old_len - 2] + 1 } else { 0 };
    wavefront[old_len] = grower.visit(old_len, (wavefront[old_len - 1] + 1).max(insertion));

    // visit the old diagonals in reverse order.
    if old_len >= 2 {
        // the old values of `i` and `i - 1` carry over from the last step
        let mut deletion = wavefront[old_len - 1];
        let mut mismatch = wavefront[old_len - 2];

        // Fill right-to-left so reads of `i`, `i - 1`, and `i - 2` still see the old wavefront.
        for i in (2..old_len).rev() {
            let insertion = wavefront[i - 2];
            wavefront[i] = grower.visit(i, deletion.max(mismatch + 1).max(insertion + 1));
            deletion = mismatch;
            mismatch = insertion;
        }

        // visit next to first diagonal.
        wavefront[1] = grower.visit(1, wavefront[1].max(wavefront[0] + 1));
    }

    // visit first diagonal.
    wavefront[0] = grower.visit(0, wavefront[0]);
    grower.grown
}

/// Conceptually, this function plays out the wavefront edit distance extension.
/// However, it assumes that the wavefronts are all at the end of the `other_seq`, so it can only add edits to the baseline.
/// It uses some sliding window shortcuts to avoid computing all of the wavefront comparisons.
/// See the test `test_growth_without_matches()` for a more intuitive example of what this computes.
/// # Arguments
/// * `wavefront` - the wavefront to replace
/// * `edit_distance` - the edit distance of `wavefront` before these edits
/// * `steps` - how many edits to apply
pub(super) fn grow_without_matches(wavefront: &mut Vec<usize>, edit_distance: usize, steps: usize) {
    // Move the old allocation out, then reuse `wavefront` for the result.
    let old = std::mem::take(wavefront);
    // Indices of old diagonals that can still raise the value written at output diagonal `j`.
    // Diagonal `i` qualifies while `j` is from `i` through `i + steps`. It contributes `old[i] + (j - i)`:
    // one more for each diagonal it has moved. The front index makes the largest contribution.
    // Each old diagonal is added at the back once, when `j` lands on it.
    // It is removed from the back right away when a newer diagonal has an equal or larger baseline offset,
    // because that newer diagonal is at least as good here and at every later `j`.
    // It is removed from the front when `j` moves past `i + steps`: its value has stopped rising.
    let mut rising: VecDeque<usize> = VecDeque::new();
    // Indices of old diagonals whose added amount has capped.
    // Diagonal `i` qualifies while `j` is from `i + steps` through `i + 2 * steps`.
    // Every one of those positions gets the same contribution, `old[i] + steps`. The front index has the largest `old` value.
    // A diagonal is added at the back when `j` reaches `i + steps`, the first position where it has stopped rising.
    // It is removed from the back right away when a newer diagonal has an equal or larger `old` value,
    // because that newer flat contribution is at least as large for the rest of the window.
    // It is removed from the front when `j` moves past `i + 2 * steps`: it can no longer reach `j`.
    let mut constant: VecDeque<usize> = VecDeque::new();
    wavefront.reserve(old.len() + 2 * steps);
    for j in 0..old.len() + 2 * steps {
        // Add old diagonal `j` at the back. First remove any back entries it makes obsolete.
        if j < old.len() {
            let score = old[j] + edit_distance - j;
            while rising.back().is_some_and(|&i| old[i] + edit_distance - i <= score) {
                rising.pop_back();
            }
            rising.push_back(j);
        }
        // Add the diagonal that just stopped rising. First remove any back entries it makes obsolete.
        if j >= steps && j - steps < old.len() {
            let i = j - steps;
            while constant.back().is_some_and(|&k| old[k] <= old[i]) {
                constant.pop_back();
            }
            constant.push_back(i);
        }
        // Remove front entries that `j` has moved past: their value is no longer rising.
        while rising.front().is_some_and(|&i| i + steps < j) {
            rising.pop_front();
        }
        // Remove front entries that `j` has moved past: they can no longer reach this diagonal.
        while constant.front().is_some_and(|&i| i + 2 * steps < j) {
            constant.pop_front();
        }
        // The front of `rising` contributes `old[i] + (j - i)`. The front of `constant` contributes `old[i] + steps`.
        let rising_value = rising.front().map_or(0, |&i| old[i] + (j - i));
        let constant_value = constant.front().map_or(0, |&i| old[i] + steps);
        wavefront.push(rising_value.max(constant_value));
    }
}

/// The inputs that each diagonal needs while a narrow wavefront grows. A method, unlike a closure, can be forced inline.
pub(super) struct WavefrontGrower<'a> {
    baseline_seq: &'a [u8],
    other_seq: &'a [u8],
    edit_distance: usize,
    offset: usize,
    wildcard: Option<u8>,
    grown: GrownWavefront,
}

impl<'a> WavefrontGrower<'a> {
    /// Creates a grower for one edit.
    /// # Arguments
    /// * `baseline_seq` - the baseline sequence, theoretically fixed
    /// * `other_seq` - the other sequence, typically getting updates
    /// * `edit_distance` - the edit distance of the grown wavefront
    /// * `offset` - the offset into `other_seq` where the baseline starts
    /// * `wildcard` - an optional wildcard symbol
    #[inline(always)]
    pub(super) fn new(baseline_seq: &'a [u8], other_seq: &'a [u8], edit_distance: usize, offset: usize, wildcard: Option<u8>) -> Self {
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
    /// # Arguments
    /// * `i` - the index of the diagonal in the grown wavefront
    /// * `d` - the distance into `other_seq` before extending
    #[inline(always)]
    fn visit(&mut self, i: usize, d: usize) -> usize {
        let d = d + match_length(self.baseline_seq, self.other_seq, d + self.edit_distance - i, d + self.offset, self.wildcard);
        if d + self.offset == self.other_seq.len() {
            self.grown.record_equality_lane(i);
        }
        self.grown.max_wavefront = self.grown.max_wavefront.max(d);
        self.grown.max_baseline = self.grown.max_baseline.max(d + self.edit_distance - i);
        d
    }
}

/// Returns how many symbols match from `baseline_seq[baseline_offset..]` and `other_seq[other_offset..]`.
/// A baseline symbol equal to `wildcard` matches anything.
/// Returns zero when either offset is at or past its sequence end.
/// # Arguments
/// * `baseline_seq` - the baseline sequence, theoretically fixed
/// * `other_seq` - the other sequence, typically getting updates
/// * `baseline_offset` - the first baseline index to compare
/// * `other_offset` - the first other index to compare
/// * `wildcard` - an optional wildcard symbol; a baseline symbol equal to it matches anything
#[inline]
pub(super) fn match_length(baseline_seq: &[u8], other_seq: &[u8], baseline_offset: usize, other_offset: usize, wildcard: Option<u8>) -> usize {
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
/// # Arguments
/// * `a` - the first sequence
/// * `b` - the second sequence
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

    /// The closed form of `grow_without_matches()` against growth one edit at a time
    #[test]
    fn test_growth_without_matches() {
        use rand::{Rng, SeedableRng};
        let mut rng = rand::rngs::StdRng::seed_from_u64(51);
        for ed in [0, 1, 2, 63, 64, 129] {
            for steps in [1, 2, 3, 63, 64, 151] {
                let mut actual: Vec<usize> = (0..2 * ed + 1)
                    .map(|i: usize| i.saturating_sub(ed) + rng.gen_range(0..50)).collect();
                let mut expected = actual.clone();
                for _ in 0..steps {
                    let mut next = vec![0; expected.len() + 2];
                    for (i, &d) in expected.iter().enumerate() {
                        next[i] = next[i].max(d);
                        next[i + 1] = next[i + 1].max(d + 1);
                        next[i + 2] = next[i + 2].max(d + 1);
                    }
                    expected = next;
                }
                grow_without_matches(&mut actual, ed, steps);
                assert_eq!(actual, expected, "distance {ed}, steps {steps}");
            }
        }
    }

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
}