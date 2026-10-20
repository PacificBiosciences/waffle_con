use rustc_hash::FxHashMap as HashMap;

use super::dwfa_lite::DWFALite;

/// Extension votes from several DWFAs.
/// It gives the same map as voting into a `HashMap<u8, f64>` with `entry()`, but it does one map insertion for each
/// symbol instead of one map lookup for each DWFA.
#[derive(Clone, Debug, Default)]
pub(crate) struct ExtensionVotes {
    /// Votes in order of the first vote for each symbol. The alphabet is small, so a linear search is fast.
    votes: Vec<(u8, f64)>,
}

impl ExtensionVotes {
    /// Adds the votes of one DWFA. Each symbol gets `weight * count / symbols.len()`.
    /// The symbols get their first votes in the same order as when voting from the map of
    /// `DWFALite::get_extension_candidates()`. Usually all symbols are equal, and then no map is needed.
    /// # Arguments
    /// * `symbols` - the output of `DWFALite::extension_symbols()` for one DWFA
    /// * `weight` - the total weight of this DWFA's vote
    /// * `orders` - cached map orders
    pub(crate) fn add(&mut self, symbols: &[u8], weight: f64, orders: &mut MapOrders) {
        let Some(&first) = symbols.first() else {
            return;
        };
        let vote_split = symbols.len() as f64;
        if symbols.iter().all(|&s| s == first) {
            self.add_vote(first, weight * symbols.len() as f64 / vote_split);
            return;
        }
        // Colliding keys make the iteration order of the final map depend on insertion order.
        // The symbols get their first votes in the iteration order of the per-DWFA count map from before.
        let mut keys = [0u8; MAX_CACHED_KEYS];
        let mut counts = [0usize; MAX_CACHED_KEYS];
        let mut num_keys = 0;
        for &s in symbols.iter() {
            match keys[..num_keys].iter().position(|&k| k == s) {
                Some(index) => counts[index] += 1,
                None if num_keys < MAX_CACHED_KEYS => {
                    keys[num_keys] = s;
                    counts[num_keys] = 1;
                    num_keys += 1;
                },
                None => {
                    // too many symbols to cache; build the map
                    let mut count_map: HashMap<u8, usize> = Default::default();
                    for &s in symbols.iter() {
                        *count_map.entry(s).or_insert(0) += 1;
                    }
                    for (&c, &occ) in count_map.iter() {
                        self.add_vote(c, weight * occ as f64 / vote_split);
                    }
                    return;
                },
            };
        }
        let order = orders.counted_order(&keys[..num_keys]);
        for position in 0..num_keys {
            let index = MapOrders::index_at(order, position);
            self.add_vote(keys[index], weight * counts[index] as f64 / vote_split);
        }
    }

    /// Adds the vote of a DWFA that nominates one symbol, the same as `add(&[symbol], weight)`.
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
            Some((_, total)) => *total += vote,
            None => self.votes.push((symbol, 0.0 + vote)),
        };
    }

    /// Removes all votes and keeps the allocation.
    pub(crate) fn clear(&mut self) {
        self.votes.clear();
    }

    /// Converts the votes into a map. Keys are inserted in order of their first vote, as `entry()` would do.
    pub(crate) fn to_map(&self) -> HashMap<u8, f64> {
        let mut candidates: HashMap<u8, f64> = Default::default();
        for &(symbol, total) in self.votes.iter() {
            candidates.insert(symbol, total);
        }
        candidates
    }

    /// Writes the candidates into `out`, in the iteration order of `to_map()`.
    /// When there are two or more candidates, `wildcard` is removed. The callers iterated this map before,
    /// so the order keeps their tie-breaks. One candidate needs no map.
    /// # Arguments
    /// * `wildcard` - an optional wildcard symbol
    /// * `out` - output buffer; it is cleared first
    /// * `orders` - cached map orders
    pub(crate) fn candidates(&self, wildcard: Option<u8>, out: &mut Vec<(u8, f64)>, orders: &mut MapOrders) {
        out.clear();
        if self.votes.len() <= 1 {
            out.extend_from_slice(&self.votes);
            return;
        }
        if self.votes.len() > MAX_CACHED_KEYS {
            let mut map = self.to_map();
            if let Some(wc) = wildcard {
                map.remove(&wc);
            }
            out.extend(map.iter().map(|(&symbol, &total)| (symbol, total)));
            return;
        }
        // removing a key does not change the order of the other keys
        let mut keys = [0u8; MAX_CACHED_KEYS];
        for (key, &(symbol, _total)) in keys.iter_mut().zip(self.votes.iter()) {
            *key = symbol;
        }
        let order = orders.inserted_order(&keys[..self.votes.len()]);
        for position in 0..self.votes.len() {
            let (symbol, total) = self.votes[MapOrders::index_at(order, position)];
            if wildcard != Some(symbol) {
                out.push((symbol, total));
            }
        }
    }
}

/// The most keys of a map order that `MapOrders` caches. The keys and their count must fit into a `u64`.
const MAX_CACHED_KEYS: usize = 7;

/// The iteration orders of small, new `FxHashMap`s with `u8` keys, cached by the insertion order of the keys.
/// The iteration order of a new map only depends on its keys and their insertion order, so a cached order is exact.
/// A lookup here is much faster than to build a map for each search step.
#[derive(Debug, Default)]
pub(crate) struct MapOrders {
    /// Orders of a `HashMap<u8, f64>` filled with `insert()`, see `ExtensionVotes::to_map()`
    inserted: HashMap<u64, u64>,
    /// Orders of a `HashMap<u8, usize>` filled with `entry()`, see `ExtensionVotes::add()`
    counted: HashMap<u64, u64>,
}

impl MapOrders {
    /// Packs up to `MAX_CACHED_KEYS` keys and their count into one `u64`
    /// # Arguments
    /// * `keys` - the distinct keys, in insertion order
    fn pack(keys: &[u8]) -> u64 {
        debug_assert!(keys.len() <= MAX_CACHED_KEYS);
        keys.iter().enumerate().fold(keys.len() as u64, |packed, (index, &key)| packed | (u64::from(key) << (8 * (index + 1))))
    }

    /// Packs an iteration order: byte `position` is the index into the keys of the key at that position
    /// # Arguments
    /// * `indices` - the index of each key, in iteration order
    fn pack_order(indices: impl Iterator<Item = usize>) -> u64 {
        indices.enumerate().fold(0, |order, (position, index)| order | ((index as u64) << (8 * position)))
    }

    /// The index into the keys of the key at `position` in a packed order
    /// # Arguments
    /// * `order` - a packed iteration order
    /// * `position` - the position in that order
    #[inline(always)]
    fn index_at(order: u64, position: usize) -> usize {
        ((order >> (8 * position)) & 0xff) as usize
    }

    /// The packed iteration order of a new `HashMap<u8, f64>` after `insert()` of the distinct `keys` in order
    /// # Arguments
    /// * `keys` - the distinct keys, in insertion order
    fn inserted_order(&mut self, keys: &[u8]) -> u64 {
        *self.inserted.entry(Self::pack(keys)).or_insert_with(|| {
            let mut map: HashMap<u8, f64> = Default::default();
            for (index, &key) in keys.iter().enumerate() {
                map.insert(key, index as f64);
            }
            Self::pack_order(map.values().map(|&index| index as usize))
        })
    }

    /// The packed iteration order of a new `HashMap<u8, usize>` after `entry()` of the distinct `keys` in order
    /// # Arguments
    /// * `keys` - the distinct keys, in insertion order
    fn counted_order(&mut self, keys: &[u8]) -> u64 {
        *self.counted.entry(Self::pack(keys)).or_insert_with(|| {
            let mut map: HashMap<u8, usize> = Default::default();
            for (index, &key) in keys.iter().enumerate() {
                map.entry(key).or_insert(index);
            }
            Self::pack_order(map.values().copied())
        })
    }
}

/// Buffers for a DWFA vote with more than one frontier diagonal
#[derive(Debug, Default)]
pub(crate) struct VoteScratch {
    pub(super) symbols: Vec<u8>,
    pub(super) orders: MapOrders,
}

impl VoteScratch {
    /// The cached map orders, e.g., for `ExtensionVotes::candidates()`
    pub(crate) fn orders(&mut self) -> &mut MapOrders {
        &mut self.orders
    }
}

/// Reusable buffers that collect the extension votes of many DWFAs, so that a search step does not allocate.
#[derive(Debug, Default)]
pub(crate) struct VoteBuffers {
    votes: ExtensionVotes,
    scratch: VoteScratch,
}

impl VoteBuffers {
    /// Removes the votes of the last search step.
    pub(crate) fn clear(&mut self) {
        self.votes.clear();
    }

    /// Adds the vote of one DWFA, see `DWFALite::vote()`.
    /// # Arguments
    /// * `dwfa` - the DWFA casting the vote
    /// * `baseline_seq` - the baseline sequence, theoretically fixed
    /// * `other_seq` - the other sequence, typically getting updates
    /// * `weight` - the total weight of this DWFA's vote
    #[inline]
    pub(crate) fn add(&mut self, dwfa: &DWFALite, baseline_seq: &[u8], other_seq: &[u8], weight: f64) {
        dwfa.vote(baseline_seq, other_seq, weight, &mut self.votes, &mut self.scratch);
    }

    /// Writes the candidates into `out`, see `ExtensionVotes::candidates()`.
    /// # Arguments
    /// * `wildcard` - an optional wildcard symbol
    /// * `out` - output buffer; it is cleared first
    pub(crate) fn candidates(&mut self, wildcard: Option<u8>, out: &mut Vec<(u8, f64)>) {
        self.votes.candidates(wildcard, out, &mut self.scratch.orders);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rustc_hash::FxHashMap as HashMap;

    #[test]
    fn test_extension_votes_match_entry_voting() {
        // `C` and `G` collide in a small FxHashMap, so key order depends on insertion order; compare order and values
        use rand::{Rng, SeedableRng};
        let mut rng = rand::rngs::StdRng::seed_from_u64(11);
        // one cache for all trials, so that most orders come from the cache; the large alphabet passes `MAX_CACHED_KEYS`
        let mut orders = MapOrders::default();
        for trial in 0..4000 {
            let alphabet: &[u8] = if trial % 4 == 0 { b"ACGT*NRYKMSW" } else { b"ACGT*N" };
            let mut expected: HashMap<u8, f64> = Default::default();
            let mut votes = ExtensionVotes::default();
            for _ in 0..rng.gen_range(0..10) {
                let symbols: Vec<u8> = (0..rng.gen_range(0..10)).map(|_| alphabet[rng.gen_range(0..alphabet.len())]).collect();
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
                votes.add(&symbols, weight, &mut orders);
            }
            let actual = votes.to_map();
            let to_bits = |m: &HashMap<u8, f64>| m.iter().map(|(&k, &v)| (k, v.to_bits())).collect::<Vec<_>>();
            assert_eq!(to_bits(&actual), to_bits(&expected));

            // `candidates()` gives the iteration order of the map, without the wildcard when there are 2+ keys
            let wildcard = if trial % 2 == 0 { Some(b'*') } else { None };
            if let Some(wc) = wildcard {
                if expected.len() > 1 {
                    expected.remove(&wc);
                }
            }
            let mut candidates = vec![];
            votes.candidates(wildcard, &mut candidates, &mut orders);
            let candidate_bits = candidates.iter().map(|&(k, v)| (k, v.to_bits())).collect::<Vec<_>>();
            assert_eq!(candidate_bits, to_bits(&expected));
        }
    }
}
