
use log::trace;
use simple_error::bail;

/// Struct for tracking the length of results in our queue that are greater than some threshold.
/// This allows us to track how long a queue actually is when we know we will ignore events smaller than the threshold.
/// Currently, it does not contain the queue itself, which may be worth doing in the long term.
/// It also tracks the number of "processed" items separately, allowing a user to set a capacity at each size and report when that capacity is reached.
#[derive(Debug)]
pub struct PQueueTracker {
    /// The count of each haplotype size in the queue
    length_counts: Vec<usize>,
    /// The total number of haplotypes in the queue with length >= threshold
    total_count: usize,
    /// The minimum threshold for a haplotype to count
    threshold: usize,
    /// The count of each haplotype size that was processed
    processed_counts: Vec<usize>,
    /// The maximum number of elements allowed at each index
    capacity_per_size: usize
}

impl PQueueTracker {
    /// Creates a new tracker with a given initial size.
    /// # Arguments
    /// * `initial_size` - the initial size of the tracked elements
    /// * `capacity_per_size` - the maximum capacity at each element
    pub fn with_capacity(initial_size: usize, capacity_per_size: usize) -> PQueueTracker {
        PQueueTracker {
            length_counts: vec![0; initial_size],
            total_count: 0,
            threshold: 0,
            processed_counts: vec![0; initial_size],
            capacity_per_size
        }
    }

    /// Inserts a length to our tracker
    /// # Arguments
    /// * `value` - the length to track
    pub fn insert(&mut self, value: usize) {
        if value >= self.length_counts.len() {
            self.length_counts.resize(value+1, 0);
        }
        self.length_counts[value] += 1;
        if value >= self.threshold {
            self.total_count += 1;
        }
    }

    /// Removes a length from the tracker
    /// # Arguments
    /// * `value` - the length to remove from tracking
    pub fn remove(&mut self, value: usize) {
        assert!(self.length_counts[value] > 0);
        self.length_counts[value] -= 1;
        if value >= self.threshold {
            assert!(self.total_count > 0);
            self.total_count -= 1;
        }
    }

    /// Wrapper for increasing by one
    pub fn increment_threshold(&mut self) {
        self.increase_threshold(self.threshold+1);
    }

    /// Increased the threshold of what is included in our total count
    /// # Arguments
    /// * `new_threshold` - the new minimum threshold to track, must be >= current threshold
    pub fn increase_threshold(&mut self, new_threshold: usize) {
        assert!(new_threshold >= self.threshold);
        trace!("increase_threshold => {}, size = {}", self.threshold, self.total_count);
        for t in self.threshold..new_threshold {
            self.total_count -= self.length_counts[t];
        }
        self.threshold = new_threshold;
        trace!("increase_threshold => {}, size = {}", self.threshold, self.total_count);
    }

    /// Attempts to mark a new value as processed.
    /// # Arguments
    /// * `value` - the length to track
    pub fn process(&mut self, value: usize) -> Result<(), Box<dyn std::error::Error>> {
        if value >= self.processed_counts.len() {
            self.processed_counts.resize(value+1, 0);
        }

        if self.processed_counts[value] >= self.capacity_per_size {
            bail!("Capacity is full");
        }

        // increment the value
        self.processed_counts[value] += 1;
        Ok(())
    }

    /// Returns the number of processed elements at a particular level
    /// # Arguments
    /// * `value` - the length to track
    pub fn processed(&self, value: usize) -> usize {
        if value >= self.processed_counts.len() {
            0
        } else {
            self.processed_counts[value]
        }
    }

    /// Returns true if we are at capacity for a particular occupancy
    /// # Arguments
    /// * `value` - the length to track
    pub fn at_capacity(&self, value: usize) -> bool {
        self.processed(value) >= self.capacity_per_size
    }

    /// Returns the total number of haplotypes in the queue with length >= the internal threshold.
    pub fn len(&self) -> usize {
        self.total_count
    }

    /// Returns the total number of tracked elements, regardless of the threshold
    pub fn unfiltered_len(&self) -> usize {
        self.length_counts.iter().sum()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn threshold(&self) -> usize {
        self.threshold
    }

    /// Returns the number of tracked elements at a particular level
    /// # Arguments
    /// * `value` - the length to track
    pub fn occupancy(&self, value: usize) -> usize {
        if value >= self.length_counts.len() {
            0
        } else {
            self.length_counts[value]
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_capacity() {
        let mut capacity_tracker = PQueueTracker::with_capacity(0, 2);

        // currently empty
        assert!(!capacity_tracker.at_capacity(1));
        assert_eq!(capacity_tracker.processed(1), 0);
        
        // add 1 of 2
        capacity_tracker.process(1).unwrap();
        assert!(!capacity_tracker.at_capacity(1));
        assert_eq!(capacity_tracker.processed(1), 1);
        
        // add 2 if 2
        capacity_tracker.process(1).unwrap();
        assert!(capacity_tracker.at_capacity(1));
        assert_eq!(capacity_tracker.processed(1), 2);

        // it's full now, this should fail without changing occupancy
        assert!(capacity_tracker.process(1).is_err());
        assert_eq!(capacity_tracker.processed(1), 2); // should not have been modified when we failed to insert
    }
}