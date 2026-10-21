
/// Extension votes from several DWFAs.
/// Symbols stay in the order of their first vote. This vector is not sorted by weight.
#[derive(Clone, Debug, Default)]
pub(crate) struct ExtensionVotes {
    /// Votes in order of the first vote for each symbol. The alphabet is small, so a linear search is fast.
    votes: Vec<(u8, f64)>,
}

impl ExtensionVotes {
    /// Adds the votes of one DWFA. Each symbol gets `weight * count / symbols.len()`.
    /// A symbol's first vote is recorded in the order it appears in `symbols`.
    /// Usually all symbols are equal, and then the whole weight goes to that symbol.
    /// That path stays bit-identical to `add_single`.
    /// # Arguments
    /// * `weight` - the total weight of this DWFA's vote
    /// * `scratch` - reusable buffers; the VoteScratch.symbols are expected to be populated by the DWFA,
    ///   and the VoteScratch.counts are cleared and modified by this function
    pub(crate) fn add(&mut self, weight: f64, scratch: &mut VoteScratch) {
        let Some(&first) = scratch.symbols.first() else {
            // This DWFA nominates nothing, so its weight is dropped.
            return;
        };

        // All symbols are the same, so the whole weight goes to that symbol.
        // This stays bit-identical to `add_single`.
        if scratch.symbols.iter().all(|&s| s == first) {
            self.add_vote(first, weight);
            return;
        }

        // otherwise, we have to split the weight, keeping the first-seen order
        let VoteScratch { symbols, counts } = scratch;
        let vote_split = symbols.len() as f64;
        counts.clear();
        for &symbol in symbols.iter() {
            match counts.iter_mut().find(|(seen, _)| *seen == symbol) {
                Some((_, count)) => *count += 1,
                None => counts.push((symbol, 1)),
            }
        }
        for &(symbol, count) in counts.iter() {
            self.add_vote(symbol, weight * count as f64 / vote_split);
        }
    }

    /// Adds the vote of a DWFA that nominates one symbol, the same as `add` with a one-symbol `VoteSymbols` buffer.
    /// # Arguments
    /// * `symbol` - the symbol this DWFA nominates
    /// * `weight` - the total weight of this DWFA's vote
    #[inline]
    pub(crate) fn add_single(&mut self, symbol: u8, weight: f64) {
        self.add_vote(symbol, weight);
    }

    /// Adds one amount to the total for a symbol.
    /// # Arguments
    /// * `symbol` - the symbol receiving the vote
    /// * `vote` - the amount to add
    #[inline]
    fn add_vote(&mut self, symbol: u8, vote: f64) {
        match self.votes.iter_mut().find(|(s, _)| *s == symbol) {
            // An existing symbol keeps its place and gains `vote`. A new symbol is appended.
            Some((_, total)) => *total += vote,
            None => self.votes.push((symbol, vote)),
        };
    }

    /// Removes all votes and keeps the allocation.
    pub(crate) fn clear(&mut self) {
        self.votes.clear();
    }

    /// Writes the candidates into `out`, in first-seen order.
    /// When there are two or more candidates, `wildcard` is removed. One candidate is kept as it is,
    /// including when that candidate is the wildcard.
    /// # Arguments
    /// * `wildcard` - an optional wildcard symbol
    /// * `out` - output buffer; it is cleared first
    pub(crate) fn candidates(&self, wildcard: Option<u8>, out: &mut Vec<(u8, f64)>) {
        out.clear();
        // One symbol, or none: keep it, even when it is the wildcard.
        if self.votes.len() <= 1 {
            out.extend_from_slice(&self.votes);
            return;
        }
        // Symbols in `votes` are distinct, so dropping the wildcard leaves at least one.
        out.extend(self.votes.iter().copied().filter(|&(symbol, _)| wildcard != Some(symbol)));
    }
}

/// Distinct symbols expected in one mixed vote. DNA, a wildcard, and the ambiguity codes fit.
/// This is only a starting capacity.
const MIXED_SYMBOL_CAPACITY: usize = 16;

/// Symbols nominated by one DWFA. `VoteScratch::clear` empties this buffer and reuses it.
/// The one-symbol path does not use it.
pub(crate) type VoteSymbols = Vec<u8>;

/// First-seen symbol counts for one mixed vote. `ExtensionVotes::add` clears this buffer and reuses it.
pub(crate) type VoteCounts = Vec<(u8, usize)>;

/// Buffers for a DWFA vote with more than one frontier diagonal.
#[derive(Debug)]
pub(crate) struct VoteScratch {
    pub(super) symbols: VoteSymbols,
    pub(super) counts: VoteCounts,
}

impl VoteScratch {
    /// Clears both buffers and keeps their allocations.
    pub(crate) fn clear(&mut self) {
        self.symbols.clear();
        self.counts.clear();
    }
}

impl Default for VoteScratch {
    fn default() -> Self {
        Self {
            symbols: Vec::new(),
            counts: Vec::with_capacity(MIXED_SYMBOL_CAPACITY),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Counts `symbols` in first-seen order and adds `weight * count / symbols.len()` to `expected`.
    fn add_expected(expected: &mut Vec<(u8, f64)>, symbols: &[u8], weight: f64) {
        let Some(&first) = symbols.first() else {
            return;
        };
        if symbols.iter().all(|&symbol| symbol == first) {
            match expected.iter_mut().find(|(symbol, _)| *symbol == first) {
                Some((_, total)) => *total += weight,
                None => expected.push((first, weight)),
            }
            return;
        }
        let vote_split = symbols.len() as f64;
        let mut counts: VoteCounts = Vec::new();
        for &symbol in symbols.iter() {
            match counts.iter_mut().find(|(seen, _)| *seen == symbol) {
                Some((_, count)) => *count += 1,
                None => counts.push((symbol, 1)),
            }
        }
        for (symbol, count) in counts {
            let vote = weight * count as f64 / vote_split;
            match expected.iter_mut().find(|(seen, _)| *seen == symbol) {
                Some((_, total)) => *total += vote,
                None => expected.push((symbol, vote)),
            }
        }
    }

    /// Tests that voting is First-Seen Order (FSO).
    #[test]
    fn test_extension_votes_are_first_seen() {
        use rand::{Rng, SeedableRng};
        let mut rng = rand::rngs::StdRng::seed_from_u64(11);
        let mut scratch = VoteScratch::default();
        for trial in 0..4000 {
            let alphabet: &[u8] = if trial % 4 == 0 { b"ACGT*NRYKMSW" } else { b"ACGT*N" };
            let mut expected: Vec<(u8, f64)> = Vec::new();
            let mut votes = ExtensionVotes::default();
            for _ in 0..rng.gen_range(0..10) {
                scratch.clear();
                scratch.symbols.extend((0..rng.gen_range(0..10)).map(|_| alphabet[rng.gen_range(0..alphabet.len())]));
                let weight = [1.0, 0.5, 0.25, 0.9][rng.gen_range(0..4)];
                add_expected(&mut expected, &scratch.symbols, weight);
                votes.add(weight, &mut scratch);
            }
            let to_bits = |pairs: &[(u8, f64)]| pairs.iter().map(|&(k, v)| (k, v.to_bits())).collect::<Vec<_>>();
            let mut actual = vec![];
            votes.candidates(None, &mut actual);
            assert_eq!(to_bits(&actual), to_bits(&expected));

            // a wildcard is removed only when another candidate remains
            let wildcard = if trial % 2 == 0 { Some(b'*') } else { None };
            if let Some(wc) = wildcard {
                if expected.len() > 1 {
                    expected.retain(|&(symbol, _)| symbol != wc);
                }
            }
            votes.candidates(wildcard, &mut actual);
            assert_eq!(to_bits(&actual), to_bits(&expected));
        }
    }
}
