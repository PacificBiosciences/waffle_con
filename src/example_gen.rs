
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
    let basem1_distribution = Uniform::new(0, alphabet_size-1);
    let error_distribution = Uniform::new(0.0, 1.0);
    let error_type_distribution = Uniform::new(0, 3);

    let consensus: Vec<u8> = (0..seq_len)
        .map(|_i| rng.sample(base_distribution))
        .collect();

    let samples: Vec<Vec<u8>> = (0..num_samples)
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
                            let alt_c = (c+sub_offset) % alphabet_size;
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
        .collect();

    (consensus, samples)
}