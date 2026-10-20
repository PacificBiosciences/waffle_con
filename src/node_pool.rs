
/// The most nodes that a pool keeps. A pool frees the nodes that it gets above this count.
const MAX_POOLED_NODES: usize = 64;

/// Recycles search nodes. A clone overwrites a recycled node with `Clone::clone_from`, so it can reuse the
/// buffers of that node instead of allocating new ones. This saves allocations only when `clone_from` of the
/// node type reuses its buffers.
pub(crate) struct NodePool<T> {
    nodes: Vec<T>,
}

impl<T: Clone> NodePool<T> {
    /// Creates an empty pool.
    pub(crate) fn new() -> Self {
        Self { nodes: Vec::new() }
    }

    /// Returns a clone of `source`. It is built in a recycled node when the pool has one.
    /// # Arguments
    /// * `source` - the node to clone
    pub(crate) fn clone_node(&mut self, source: &T) -> T {
        match self.nodes.pop() {
            Some(mut node) => {
                node.clone_from(source);
                node
            },
            None => source.clone(),
        }
    }

    /// Keeps `node` for a later clone.
    /// # Arguments
    /// * `node` - the node to keep
    pub(crate) fn recycle(&mut self, node: T) {
        if self.nodes.len() < MAX_POOLED_NODES {
            self.nodes.push(node);
        }
    }
}
