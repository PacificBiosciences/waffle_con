mod common;

use common::edit_distance;
use waffle_con::dynamic_wfa::{DWFALite, DWFALiteConfig, DWFALiteState};

fn config(wildcard: Option<u8>, early: bool, cap: Option<usize>) -> DWFALiteConfig {
    DWFALiteConfig {
        wildcard,
        allow_early_termination: early,
        max_edit_distance: cap,
    }
}

fn assert_terminal(actual: &DWFALite, distance: usize, context: &str) {
    let cap = actual.max_edit_distance();
    let expected_state = if cap.is_some_and(|limit| distance > limit) {
        DWFALiteState::ExceededEditDistanceLimit
    } else {
        DWFALiteState::Finalized
    };
    assert_eq!(
        (actual.edit_distance(), actual.state()),
        (
            cap.map_or(distance, |limit| distance.min(limit)),
            expected_state
        ),
        "{context}: uncapped DP={distance}, cap={cap:?}"
    );
}

#[test]
fn exhaustive_short_endpoints_and_caps_match_dp() {
    let strings: Vec<Vec<u8>> = (0..=3)
        .flat_map(|length| {
            (0..3usize.pow(length)).map(move |mut encoded| {
                (0..length)
                    .map(|_| {
                        let byte = [0, 1, 255][encoded % 3];
                        encoded /= 3;
                        byte
                    })
                    .collect()
            })
        })
        .collect();

    for baseline in &strings {
        for suffix in &strings {
            for wildcard in [None, Some(255)] {
                for early in [false, true] {
                    // One independent DP calculation per endpoint, not per update or cap.
                    let distance = edit_distance(baseline, suffix, wildcard, early);
                    for offset in [0, 2] {
                        let mut other = [17, 18][..offset].to_vec();
                        other.extend_from_slice(suffix);
                        for cap in [None, Some(0), Some(1), Some(3)] {
                            for incremental in [false, true]
                                .into_iter()
                                .take(if suffix.is_empty() { 1 } else { 2 })
                            {
                                let mut actual = DWFALite::new(config(wildcard, early, cap));
                                actual.set_offset(offset);
                                if incremental {
                                    for end in offset..=other.len() {
                                        actual.update(baseline, &other[..end]).unwrap();
                                    }
                                } else {
                                    actual.update(baseline, &other).unwrap();
                                }
                                actual.finalize(baseline, &other).unwrap();
                                let expected_state = if cap.is_some_and(|limit| distance > limit) {
                                    DWFALiteState::ExceededEditDistanceLimit
                                } else {
                                    DWFALiteState::Finalized
                                };
                                assert_eq!(
                                    (actual.edit_distance(), actual.state(), actual.max_edit_distance()),
                                    (cap.map_or(distance, |limit| distance.min(limit)), expected_state, cap),
                                    "baseline={baseline:?}, suffix={suffix:?}, wildcard={wildcard:?}, early={early}, offset={offset}, cap={cap:?}, incremental={incremental}"
                                );
                            }
                        }
                    }
                }
            }
        }
    }
}

/// A scalar scatter recurrence. Every extension scans every lane byte by byte;
/// maxima are recomputed, and finalization grows one edit at a time.
#[derive(Clone)]
struct ScalarWavefront {
    config: DWFALiteConfig,
    state: DWFALiteState,
    score: usize,
    lanes: Vec<usize>,
    offset: usize,
}

impl ScalarWavefront {
    fn new(config: DWFALiteConfig, offset: usize) -> Self {
        Self {
            config,
            state: DWFALiteState::Active,
            score: 0,
            lanes: vec![0],
            offset,
        }
    }

    fn maxima(&self) -> (usize, usize) {
        let baseline = self
            .lanes
            .iter()
            .enumerate()
            .map(|(lane, &position)| position + self.score - lane)
            .max()
            .unwrap();
        let other = self.offset + self.lanes.iter().copied().max().unwrap();
        (baseline, other)
    }

    fn extend(&mut self, baseline: &[u8], other: &[u8]) {
        for (lane, position) in self.lanes.iter_mut().enumerate() {
            let mut row = *position + self.score - lane;
            let mut column = *position + self.offset;
            while row < baseline.len() && column < other.len() {
                if baseline[row] != other[column] && self.config.wildcard != Some(baseline[row]) {
                    break;
                }
                *position += 1;
                row += 1;
                column += 1;
            }
        }
    }

    fn grow(&mut self, baseline: &[u8], other: &[u8]) -> bool {
        if self
            .config
            .max_edit_distance
            .is_some_and(|limit| self.score >= limit)
        {
            self.state = DWFALiteState::ExceededEditDistanceLimit;
            return false;
        }
        let mut next = vec![0; self.lanes.len() + 2];
        // Scatter deletion, substitution and insertion independently from each old lane.
        for (lane, &position) in self.lanes.iter().enumerate() {
            next[lane] = next[lane].max(position);
            next[lane + 1] = next[lane + 1].max(position + 1);
            next[lane + 2] = next[lane + 2].max(position + 1);
        }
        self.lanes = next;
        self.score += 1;
        self.extend(baseline, other);
        true
    }

    fn update(&mut self, baseline: &[u8], other: &[u8]) {
        if self.state != DWFALiteState::Active {
            return;
        }
        self.extend(baseline, other);
        while self.maxima().1 < other.len()
            && !(self.config.allow_early_termination && self.maxima().0 >= baseline.len())
        {
            if !self.grow(baseline, other) {
                break;
            }
        }
    }

    fn finalize(&mut self, baseline: &[u8], other: &[u8]) {
        if self.state != DWFALiteState::Active {
            return;
        }
        loop {
            // Ordinary mode needs a single path to the joint terminal cell, not
            // two unrelated lanes furnishing the two coordinate maxima.
            let diagonal = (other.len() - self.offset) as isize - baseline.len() as isize;
            let terminal_lane = self.score as isize + diagonal;
            let complete = if self.config.allow_early_termination {
                self.maxima().0 >= baseline.len()
            } else {
                terminal_lane >= 0
                    && (terminal_lane as usize) < self.lanes.len()
                    && self.lanes[terminal_lane as usize] >= other.len() - self.offset
            };
            if complete {
                self.state = DWFALiteState::Finalized;
                return;
            }
            if !self.grow(baseline, other) {
                return;
            }
        }
    }
}

fn assert_state(
    actual: &DWFALite,
    expected: &ScalarWavefront,
    baseline: &[u8],
    other: &[u8],
    context: &str,
) {
    let (baseline_max, other_max) = expected.maxima();
    assert_eq!(
        (
            actual.edit_distance(),
            actual.state(),
            actual.wavefront(),
            actual.maximum_baseline_distance(),
            actual.maximum_other_distance(),
            actual.reached_baseline_end(baseline),
            actual.max_edit_distance(),
        ),
        (
            expected.score,
            expected.state,
            expected.lanes.as_slice(),
            baseline_max,
            other_max,
            baseline_max == baseline.len(),
            expected.config.max_edit_distance,
        ),
        "{context}: config={:?}, offset={}, baseline_len={}, other_len={}",
        expected.config,
        expected.offset,
        baseline.len(),
        other.len()
    );
}

fn update_pair(
    actual: &mut DWFALite,
    expected: &mut ScalarWavefront,
    baseline: &[u8],
    other: &[u8],
    context: &str,
) {
    expected.update(baseline, other);
    assert_eq!(
        actual.update(baseline, other).unwrap(),
        expected.score,
        "{context}: update return value"
    );
    assert_state(actual, expected, baseline, other, context);
}

fn check_clone_branches(
    source: &DWFALite,
    model: &ScalarWavefront,
    baseline: &[u8],
    other: &[u8],
    wider: &DWFALite,
    context: &str,
) -> Vec<(DWFALite, ScalarWavefront, Vec<u8>)> {
    // These are live public wavefront lengths, not assertions about spare capacity.
    let mut smaller = DWFALite::new(config(Some(42), true, Some(0)));
    assert!(smaller.wavefront().len() < source.wavefront().len());
    assert!(wider.wavefront().len() > source.wavefront().len());
    let mut larger = wider.clone();
    smaller.clone_from(source);
    larger.clone_from(source);
    let copies = [source.clone(), smaller, larger];
    let mut finished: Vec<(DWFALite, ScalarWavefront, Vec<u8>)> = Vec::with_capacity(copies.len());
    for (mut branch, byte) in copies.into_iter().zip([b'T', 2, 255]) {
        let mut expected = model.clone();
        assert_state(&branch, &expected, baseline, other, context);
        let mut branch_other = other.to_vec();
        branch_other.push(byte);
        update_pair(&mut branch, &mut expected, baseline, &branch_other, context);
        branch.finalize(baseline, &branch_other).unwrap();
        expected.finalize(baseline, &branch_other);
        assert_state(&branch, &expected, baseline, &branch_other, context);
        let distance = edit_distance(
            baseline,
            &branch_other[model.offset..],
            model.config.wildcard,
            model.config.allow_early_termination,
        );
        assert_terminal(&branch, distance, context);
        assert_state(source, model, baseline, other, context);
        // Later branches must not change earlier independently finalized copies.
        for (previous, previous_model, previous_other) in &finished {
            assert_state(previous, previous_model, baseline, previous_other, context);
        }
        finished.push((branch, expected, branch_other));
    }
    finished
}

fn compare_scenario(
    baseline: &[u8],
    suffix: &[u8],
    caps: &[Option<usize>],
    wider: &DWFALite,
    context: &str,
) {
    for wildcard in [None, Some(255)] {
        let mut baseline = baseline.to_vec();
        if wildcard.is_some() {
            for index in (0..baseline.len()).step_by(17) {
                baseline[index] = 255;
            }
        }
        for early in [false, true] {
            let distance = edit_distance(&baseline, suffix, wildcard, early);
            for offset in [0, 2] {
                let mut other = [17, 18][..offset].to_vec();
                other.extend_from_slice(suffix);
                for &cap in caps {
                    for chunk in [suffix.len().max(1), 1, 7] {
                        let options = config(wildcard, early, cap);
                        let mut actual = DWFALite::new(options.clone());
                        actual.set_offset(offset);
                        let mut expected = ScalarWavefront::new(options, offset);
                        assert_state(&actual, &expected, &baseline, &other[..offset], context);
                        // Include the valid offset-only prefix, then append whole, one,
                        // or seven bytes at a time without changing either sequence.
                        let mut end = offset;
                        loop {
                            update_pair(
                                &mut actual,
                                &mut expected,
                                &baseline,
                                &other[..end],
                                context,
                            );
                            if end == other.len() {
                                break;
                            }
                            end = (end + chunk).min(other.len());
                        }
                        // Branch once per chunking matrix. A zero-edit source has
                        // no genuinely smaller public wavefront to clone into.
                        let branches = if chunk == 7 && actual.wavefront().len() > 1 {
                            check_clone_branches(
                                &actual, &expected, &baseline, &other, wider, context,
                            )
                        } else {
                            Vec::new()
                        };
                        actual.finalize(&baseline, &other).unwrap();
                        expected.finalize(&baseline, &other);
                        assert_state(&actual, &expected, &baseline, &other, context);
                        assert_terminal(&actual, distance, context);
                        // Finalizing the source must not mutate any of its finished branches.
                        for (branch, branch_model, branch_other) in branches {
                            assert_state(&branch, &branch_model, &baseline, &branch_other, context);
                        }
                    }
                }
            }
        }
    }
}

fn wider_destination() -> DWFALite {
    let mut destination = DWFALite::new(config(Some(42), false, Some(200)));
    destination.update(&[0; 180], &[1; 180]).unwrap();
    destination
}

#[test]
fn boundary_wavefronts_and_clone_branches_match_scalar_reference() {
    let wider = wider_destination();
    for edits in [15, 16, 31, 32, 63, 64, 65, 127, 128] {
        let baseline = [vec![b'A'; edits], vec![b'G'; 80], vec![b'T']].concat();
        let other = [vec![b'C'; edits], vec![b'G'; 80]].concat();
        compare_scenario(
            &baseline,
            &other,
            &[None, Some(edits - 1), Some(edits), Some(edits + 1)],
            &wider,
            "boundary match run",
        );
    }
}

#[test]
fn mixed_bytes_indels_holes_and_long_finalization_match_scalar_reference() {
    let wider = wider_destination();
    let pattern: Vec<u8> = (0..160).map(|index| (index * 37) as u8).collect();
    // Nonuniform bytes make shifted lanes disagree, unlike the repeated-byte boundaries.
    let inserted = [
        pattern[..73].to_vec(),
        vec![254; 20],
        pattern[73..].to_vec(),
    ]
    .concat();
    compare_scenario(
        &pattern,
        &inserted,
        &[None, Some(19), Some(20), Some(21)],
        &wider,
        "interior insertion",
    );
    let deleted = [pattern[..73].to_vec(), pattern[93..].to_vec()].concat();
    compare_scenario(
        &pattern,
        &deleted,
        &[None, Some(19), Some(20), Some(21)],
        &wider,
        "interior deletion",
    );

    // A short other sequence leaves many no-match edits during finalization,
    // continuing well beyond the narrow/wide transition.
    compare_scenario(
        &pattern,
        &pattern[..9],
        &[None, Some(150), Some(151), Some(152)],
        &wider,
        "long deletion tail",
    );

    // Two terminal lanes separated by a hole exercise a noncontiguous frontier;
    // appending byte 2 in a clone then changes which lanes can extend.
    let baseline = [vec![0; 65], vec![1, 2, 1]].concat();
    let other = [vec![255; 65], vec![2, 1]].concat();
    let mut fixture = ScalarWavefront::new(config(None, false, None), 0);
    fixture.update(&baseline, &other);
    let terminal_lanes: Vec<usize> = fixture
        .lanes
        .iter()
        .enumerate()
        .filter_map(|(lane, &position)| (position == other.len()).then_some(lane))
        .collect();
    assert_eq!((fixture.score, terminal_lanes), (66, vec![65, 67]));
    compare_scenario(
        &baseline,
        &other,
        &[None, Some(65), Some(66), Some(67)],
        &wider,
        "separated terminal lanes",
    );
}
