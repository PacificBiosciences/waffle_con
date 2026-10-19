# Unreleased
## Features
* Added optional edit-distance caps to the config so a read stops contributing once it diverges too far: a hard maximum, a fraction of the read length, and a floor on the resolved cap. By default, none of these are set, matching previous behavior. Set them to improve performance when highly divergent reads are part of the input set.
* Single consensus results now report which reads were halted at that cap
* Dual consensus results now report which allele each read belongs to, including reads that hit the cap on both alleles

## Performance
* Replaced the search hash queue with a binary heap ordered by cached priority

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
