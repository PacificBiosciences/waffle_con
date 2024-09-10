
use std::cmp::max;

/// Returns the full edit distance between two u8 Vecs by using a version of WFA.
/// # Arguments
/// * `v1` - the first Vec
/// * `v2` - the second Vec
/// # Examples
/// ```rust
/// use waffle_con::sequence_alignment::wfa_ed;
/// let v1: Vec<u8> = vec![0, 1, 2, 4, 5];
/// let v2: Vec<u8> = vec![0, 1, 3, 4, 5];
/// let v3: Vec<u8> = vec![1, 2, 3, 5];
/// assert_eq!(wfa_ed(&v1, &v1), 0);
/// assert_eq!(wfa_ed(&v1, &v2), 1);
/// assert_eq!(wfa_ed(&v1, &v3), 2);
/// ```
pub fn wfa_ed(v1: &[u8], v2: &[u8]) -> usize {
    wfa_ed_config(v1, v2, true, Some(b'*'))
}

/// Returns the edit distance between two u8 Vecs by using a version of WFA.
/// If `require_both_end` is true, it requires the full end-to-end edit distance. If false, then it only requires that `v2` is at the end.
/// # Arguments
/// * `v1` - the first Vec
/// * `v2` - the second Vec
/// * `require_both_end` - if true, it requires the full end-to-end edit distance; if false, then it only requires that `v2` is at the end
/// # Examples
/// ```rust
/// use waffle_con::sequence_alignment::wfa_ed_config;
/// let v1: Vec<u8> = vec![0, 1, 2, 4, 5];
/// let v2: Vec<u8> = vec![0, 1, 2, 4];
/// assert_eq!(wfa_ed_config(&v1, &v2, false, Some(b'*')), 0);
/// assert_eq!(wfa_ed_config(&v1, &v2, true, Some(b'*')), 1);
/// ```
pub fn wfa_ed_config(v1: &[u8], v2: &[u8], require_both_end: bool, wildcard: Option<u8>) -> usize {
    //we need the lengths to know where we are in the vecs
    let l1 = v1.len();
    let l2 = v2.len();

    //stores the next indices that should be compared
    let mut curr_wf: Vec<(usize, usize)> = vec![(0, 0)];
    let mut next_wf: Vec<(usize, usize)> = vec![(0, 0); 3];
    let mut edits = 0;

    //main idea is to iterate until we're at the end of BOTH vecs, this is guaranteed because i and j monotonically increase
    loop {
        //during each iteration, we go over all wavefronts; at iteration e, there are 2*e+1 current wavefronts that will generate 2*(e+1)+1 wavefronts
        //"e" in this context corresponds to the edit distance "edits"
        for (wf_index, &wf) in curr_wf.iter().enumerate() {
            let mut i = wf.0;
            let mut j = wf.1;

            // as long as the symbols match, keep moving along the diagonal
            while i < l1 && j < l2 && (v1[i] == v2[j] || wildcard.map_or(false, |w| v1[i] == w) || wildcard.map_or(false, |w| w == v2[j])) {
                i += 1;
                j += 1;
            }
            
            if (i == l1 || !require_both_end) && j == l2 {
                //we found the end, return the number of edits required to get here
                return edits;
            }
            else if i == l1 {
                //push the wavefront, but i cannot increase
                next_wf[wf_index] = max(next_wf[wf_index], (i, j));
                next_wf[wf_index+1] = max(next_wf[wf_index+1], (i, j+1));
                next_wf[wf_index+2] = max(next_wf[wf_index+2], (i, j+1));
            } else if j == l2 {
                //push the wavefront, but j cannot increase
                next_wf[wf_index] = max(next_wf[wf_index], (i+1, j));
                next_wf[wf_index+1] = max(next_wf[wf_index+1], (i+1, j));
                next_wf[wf_index+2] = max(next_wf[wf_index+2], (i, j));
            } else {
                //v1 and v2 do not match at i, j; add mismatch, insert, and del to the next wavefront
                next_wf[wf_index] = max(next_wf[wf_index], (i+1, j)); //v2 has a deletion relative to v1
                next_wf[wf_index+1] = max(next_wf[wf_index+1], (i+1, j+1)); //v2 has a mismatch relative to v1
                next_wf[wf_index+2] = max(next_wf[wf_index+2], (i, j+1)); //v2 has an insertion relative to v1
            }
        }

        //we finished this wave, increment the edit count and generate the buffer for the next wavefront
        edits += 1;
        curr_wf = next_wf;
        next_wf = vec![(0, 0); 3+2*edits];
    }
}
