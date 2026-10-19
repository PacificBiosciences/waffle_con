
use std::collections::BinaryHeap;

use crate::cdwfa_config::ConsensusCost;

/// Search node that stores the rank [`PriorityQueue`] uses for ordering.
pub trait CachePriority {
    /// Caches the rank for `consensus_cost`, using `tie_breaker` when the other rank fields match.
    /// A lower tie-breaker was queued earlier and wins ties.
    /// # Arguments
    /// * `consensus_cost` - cost model used to score the node
    /// * `tie_breaker` - insertion order; a lower value is earlier and wins ties
    fn cache_priority(&mut self, consensus_cost: ConsensusCost, tie_breaker: u64);
}

/// Max-heap that caches a node's search priority before every push.
/// The tie-breaker is the number of nodes pushed so far, so earlier nodes win ties.
pub struct PriorityQueue<T> {
    /// The actual queue is stored as a binary heap.
    /// The type `T` must implement [`CachePriority`] and [`Ord`].
    heap: BinaryHeap<T>,
    /// Tie-breaker that the next [`push`](Self::push) will assign.
    next_id: u64,
    /// Cost model used to score the nodes.
    consensus_cost: ConsensusCost,
}

impl<T> PriorityQueue<T>
where
    T: CachePriority + Ord,
{
    /// Creates an empty queue that scores nodes with `consensus_cost`.
    pub fn new(consensus_cost: ConsensusCost) -> Self {
        Self {
            heap: BinaryHeap::new(),
            next_id: 0,
            consensus_cost,
        }
    }

    /// Tie-breaker that the next [`push`](Self::push) will assign.
    pub fn next_id(&self) -> u64 {
        self.next_id
    }

    /// Caches `node`'s priority from the current cost model, then pushes it.
    pub fn push(&mut self, mut node: T) {
        node.cache_priority(self.consensus_cost, self.next_id);
        self.next_id += 1;
        self.heap.push(node);
    }

    /// Removes the greatest node, or `None` when the queue is empty.
    pub fn pop(&mut self) -> Option<T> {
        self.heap.pop()
    }

    /// Number of nodes currently in the queue.
    pub fn len(&self) -> usize {
        self.heap.len()
    }

    /// Returns true when the queue contains no nodes.
    pub fn is_empty(&self) -> bool {
        self.heap.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use std::cmp::{Ordering, Reverse};

    use super::*;

    /// Item that implements [`CachePriority`] and [`Ord`], just for testing.
    #[derive(Debug, Eq, PartialEq)]
    struct Item {
        cost: usize,
        rank: Option<(Reverse<usize>, Reverse<u64>)>,
    }

    impl CachePriority for Item {
        fn cache_priority(&mut self, _consensus_cost: ConsensusCost, tie_breaker: u64) {
            self.rank = Some((Reverse(self.cost), Reverse(tie_breaker)));
        }
    }

    impl PartialOrd for Item {
        fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
            Some(self.cmp(other))
        }
    }

    impl Ord for Item {
        fn cmp(&self, other: &Self) -> Ordering {
            self.rank.cmp(&other.rank)
        }
    }

    /// Tests that [`push`](PriorityQueue::push) caches the priority before inserting the node.
    #[test]
    fn push_caches_priority_before_insert() {
        let mut queue = PriorityQueue::new(ConsensusCost::L1Distance);
        assert_eq!(queue.next_id(), 0);

        queue.push(Item { cost: 1, rank: None });
        queue.push(Item { cost: 1, rank: None });
        queue.push(Item { cost: 0, rank: None });

        assert_eq!(queue.len(), 3);
        assert_eq!(queue.next_id(), 3);

        let first = queue.pop().unwrap();
        let second = queue.pop().unwrap();
        let third = queue.pop().unwrap();
        assert!(queue.is_empty());

        assert_eq!(first.rank, Some((Reverse(0), Reverse(2))));
        assert_eq!(second.rank, Some((Reverse(1), Reverse(0))));
        assert_eq!(third.rank, Some((Reverse(1), Reverse(1))));
    }
}
