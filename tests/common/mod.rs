/// Wagner-Fischer edit distance, with wildcard matching only on the baseline.
/// Early termination consumes the baseline but may leave an other-sequence suffix.
pub fn edit_distance(baseline: &[u8], other: &[u8], wildcard: Option<u8>, early: bool) -> usize {
    let mut row: Vec<usize> = (0..=other.len()).collect();
    for (index, &base) in baseline.iter().enumerate() {
        let mut diagonal = row[0];
        row[0] = index + 1;
        for (column, &symbol) in other.iter().enumerate() {
            let above = row[column + 1];
            let substitution = usize::from(base != symbol && wildcard != Some(base));
            row[column + 1] = (diagonal + substitution)
                .min(above + 1)
                .min(row[column] + 1);
            diagonal = above;
        }
    }
    if early {
        *row.iter().min().unwrap()
    } else {
        row[other.len()]
    }
}
