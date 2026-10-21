# Unreleased
## Features
* Added optional edit-distance caps to the config so a read stops contributing once it diverges too far: a hard maximum, a fraction of the read length, and a floor on the resolved cap. By default, none of these are set, matching previous behavior. Set them to improve performance when highly divergent reads are part of the input set.
* Single consensus results report which reads were halted at that cap. `Consensus::assignments()` is `None` when every read was included.
* Dual consensus results report which allele each read belongs to. A read can tie, or hit the cap on both alleles.
* Numerous performance enhancements were incorporated. All consensus modes should run faster.

## Bugs
* Fixed global DWFA finalization so it requires the diagonal that reaches both sequence ends. It previously stopped when any diagonal reached the end of the baseline. Early termination still allows an unused suffix of the other sequence.

## Breaking changes
* `Consensus::new` now takes per-read assignments and returns a `Result`
* `DualConsensus::new` takes `Vec<SequenceAssignment>` instead of `Vec<bool>`, and `is_consensus1()` is replaced by `assignments()`
* `DWFALite::new` takes a `DWFALiteConfig` instead of a wildcard and an early-termination flag
* `DWFALite::update` returns `(edit_distance, Option<u8>)` instead of the edit distance alone. The optional symbol is the next baseline symbol when one frontier diagonal advanced.
* `DWFALite::get_extension_candidates` is no longer public
* The order of tie-breaking has changed, which may lead to slightly different results when there is not a dominant symbol

# v0.4.4
## Bugs
* Fixed a panic caused by input sequences shorter than the offset compare length

# v0.4.3
## Bugs
* Fixed a panic caused by an unchecked `unwrap()` in a `trace!` statement

# v0.4.2
Initial release
