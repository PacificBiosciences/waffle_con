#[cfg(test)]
use crate::cdwfa_config::CdwfaConfig;
use super::{DWFALite, DWFALiteConfig, DWFALiteConfigBuilder, DWFALiteState};

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
#[allow(deprecated)]
fn test_exact_match_under_edit_distance_cap() {
    let sequence = b"ACGTACGT";
    let mut dwfa = dwfa_with_consensus_cap(sequence.len(), Some(2), None, None).unwrap();
    let other = &sequence[..4];
    assert_eq!(dwfa.update(sequence, other).unwrap(), 0);
    assert_eq!(dwfa.get_extension_candidates(sequence, other).get(&b'A'), Some(&1));
    assert_eq!(dwfa.state(), DWFALiteState::Active);
}

#[test]
#[allow(deprecated)]
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
#[allow(deprecated)]
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
