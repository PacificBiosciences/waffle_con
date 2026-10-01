
use rand::distributions::Uniform;
use rand::{Rng, SeedableRng};

/// Creates a test set we can verify is working
/// # Arguments
/// * `alphabet_size` - the length of the alphabet, e.g. for DNA it's 4
/// * `seq_len` - the length of the sequences
/// * `num_samples` - the number of samples to generate from the consensus
/// * `error_rate` - overall error rate, assumes mismatch, insertion, and deletion are equally likely sub-components of this error rate
pub fn generate_test(alphabet_size: u8, seq_len: usize, num_samples: usize, error_rate: f64) -> (Vec<u8>, Vec<Vec<u8>>) {
    assert!(alphabet_size > 1);
    assert!((0.0..=1.0).contains(&error_rate));

    let mut rng = rand::rngs::StdRng::seed_from_u64(0);
    let base_distribution = Uniform::new(0, alphabet_size);

    let consensus: Vec<u8> = (0..seq_len)
        .map(|_i| rng.sample(base_distribution))
        .collect();
    let samples = sample_reads(&mut rng, &consensus, alphabet_size, num_samples, error_rate);

    (consensus, samples)
}

/// Creates a dual-allele test set.
/// Allele 1 is random. Allele 2 is that sequence with `num_differences` substitutions at distinct positions.
/// `num_samples` reads are then generated from each allele.
/// Reads from allele 1 come first, followed by reads from allele 2.
/// # Arguments
/// * `alphabet_size` - the length of the alphabet, e.g. for DNA it's 4
/// * `seq_len` - the length of both alleles
/// * `num_samples` - the number of reads to generate from each allele
/// * `error_rate` - per-base mismatch, insertion, and deletion rate applied independently to each read
/// * `num_differences` - how many positions differ between the two true alleles
pub fn generate_dual_test(
    alphabet_size: u8,
    seq_len: usize,
    num_samples: usize,
    error_rate: f64,
    num_differences: usize,
) -> (Vec<u8>, Vec<u8>, Vec<Vec<u8>>) {
    assert!(alphabet_size > 1);
    assert!(seq_len > 0);
    assert!((0.0..=1.0).contains(&error_rate));
    assert!(num_differences > 0 && num_differences <= seq_len);

    let mut rng = rand::rngs::StdRng::seed_from_u64(0);
    let base_distribution = Uniform::new(0, alphabet_size);

    let allele1: Vec<u8> = (0..seq_len)
        .map(|_i| rng.sample(base_distribution))
        .collect();

    // Pick distinct positions and change each one to a different symbol.
    let mut positions: Vec<usize> = (0..seq_len).collect();
    for i in 0..num_differences {
        let j = rng.gen_range(i..seq_len);
        positions.swap(i, j);
    }

    let mut allele2 = allele1.clone();
    for &position in positions.iter().take(num_differences) {
        let c = allele2[position];
        // Offset is at least 1, so this site always changes.
        let sub_offset = rng.gen_range(1..alphabet_size);
        allele2[position] = (c + sub_offset) % alphabet_size;
    }

    let mut reads = sample_reads(&mut rng, &allele1, alphabet_size, num_samples, error_rate);
    reads.extend(sample_reads(&mut rng, &allele2, alphabet_size, num_samples, error_rate));

    (allele1, allele2, reads)
}

/// Samples noisy reads from one consensus using mismatch, deletion, and insertion at equal rates.
fn sample_reads<R: Rng>(
    rng: &mut R,
    consensus: &[u8],
    alphabet_size: u8,
    num_samples: usize,
    error_rate: f64,
) -> Vec<Vec<u8>> {
    let base_distribution = Uniform::new(0, alphabet_size);
    let basem1_distribution = Uniform::new(0, alphabet_size - 1);
    let error_distribution = Uniform::new(0.0, 1.0);
    let error_type_distribution = Uniform::new(0, 3);

    (0..num_samples)
        .map(|_i| {
            let mut seq = vec![];
            let mut con_index = 0;
            while con_index < consensus.len() {
                let c = consensus[con_index];
                let is_error = rng.sample(error_distribution) < error_rate;
                if is_error {
                    let error_type = rng.sample(error_type_distribution);
                    match error_type {
                        0 => {
                            // substition
                            let sub_offset = rng.sample(basem1_distribution);
                            let alt_c = (c + sub_offset) % alphabet_size;
                            seq.push(alt_c);
                            con_index += 1;
                        },
                        1 => {
                            // deletion
                            con_index += 1;
                        },
                        2 => {
                            //insertion
                            let s = rng.sample(base_distribution);
                            seq.push(s);
                        },
                        _ => panic!("no impl")
                    }
                } else {
                    seq.push(c);
                    con_index += 1;
                }
            }
            seq
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dual_alleles_differ_at_requested_sites() {
        let (allele1, allele2, reads) = generate_dual_test(4, 50, 5, 0.0, 3);

        assert_eq!(allele1.len(), 50);
        assert_eq!(allele2.len(), 50);
        assert_eq!(
            allele1.iter().zip(allele2.iter()).filter(|(a, b)| a != b).count(),
            3
        );
        assert_eq!(reads.len(), 10);
        assert!(reads[..5].iter().all(|read| read == &allele1));
        assert!(reads[5..].iter().all(|read| read == &allele2));
        assert!(reads.iter().flatten().all(|&base| base < 4));
    }

    #[test]
    fn dual_generation_is_deterministic() {
        let first = generate_dual_test(4, 40, 4, 0.02, 2);
        let second = generate_dual_test(4, 40, 4, 0.02, 2);
        assert_eq!(first, second);
    }
}