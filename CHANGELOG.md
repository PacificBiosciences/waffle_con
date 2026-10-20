# Unreleased
## Features
* Added optional edit-distance caps to the config so a read stops contributing once it diverges too far: a hard maximum, a fraction of the read length, and a floor on the resolved cap. By default, none of these are set, matching previous behavior. Set them to improve performance when highly divergent reads are part of the input set.
* Single consensus results now report which reads were halted at that cap
* Dual consensus results now report which allele each read belongs to, including reads that hit the cap on both alleles

## Bugs
* Correct global DWFA finalization to require the joint terminal endpoint, including mode-correct bulk finalization and endpoints reached exactly at an edit-distance cap. Early termination continues to allow an unused suffix of the other sequence.

## Performance
* Replaced the search hash queue with a binary heap ordered by cached priority
* `DWFALite` extends only the diagonals at the frontier after a push, and caches its maximum distances
* `DWFALite::update` has a fast path for one matching frontier diagonal
* `DWFALite` grows a wavefront in place, and a clone gets spare capacity for the next edits
* `DWFALite::finalize` applies the remaining edits at once when no diagonal can match again
* Consensus nodes collect their extension votes without hash maps or allocations
* The last child of a consensus node takes the parent node instead of a clone, and the dual search reuses the buffers of dropped nodes
* The consensus and dual consensus benchmarks are 9.9x faster than v0.4.4 (geometric mean). Equal-vote ties follow the order each symbol is first seen, so a capped search can differ from v0.4.4. Uncapped results are still sorted by sequence.

## Breaking changes
* `Consensus::new` now takes per-read assignments and returns a `Result`

# v0.4.4
## Bugs
* Fixed a panic caused by input sequences shorter than the offset compare length

# v0.4.3
## Bugs
* Fixed a panic caused by an unchecked `unwrap()` in a `trace!` statement

# v0.4.2
Initial release
