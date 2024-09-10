[![Build status](https://github.com/PacificBiosciences/waffle_con/actions/workflows/test-ci.yml/badge.svg)](https://github.com/PacificBiosciences/waffle_con/actions)


# waffle_con
This crate contains our implementation of the Dynamic WFA (DWFA) consensus algorithms, or `waffle_con`.
The algorithms contained within were built to support the consensus steps of [pb-StarPhase](https://github.com/PacificBiosciences/pb-StarPhase).
Any issues related to StarPhase should be reported on the StarPhase GitHub page.

This crate contains functionality for:

* ConsensusDWFA - One input string per sequence, one output sequence
* DualConsensusDWFA - One input string per sequence, 1 or 2 output sequences
* PriorityConsensuDWFA - Multiple input strings per sequence (priority), 1+ output sequences

## Full documentation
`waffle_con` provides extensive in-line documentation.
A user-friendly HTML version can be generated via `cargo doc`.

## Methods
At a high level, this project provided many consensus algorithms that slowly build upon each other from a baseline single-consensus method.
The single consensus method is based on the idea of cost-based exploration of an assembly (or consensus) space.
"Cost" in this context is basically the edit distance between the assembled sequence and a set of inputs.
We use a dynamic WFA algorithm to both nominate and score the assembled sequences.
These sequences are then explored using a Dijkstra-like approach (least-cost-first).
Thus the core loop is:

* Pop a candidate from the priority queue - check if this candidate is finished and compare to current best results
* If not finished, each dynamic WFA nominates one or more extensions to the candidate
* Each candidate extension is added to the dynamic WFA for each sequence
* The total combined edit distance is used to score that candidate
* The candidate is placed into the min-cost priority queue

The algorithm can be extended to the "dual" option by allowing it to split into two candidate sequences when sufficient evidence is present (e.g., sufficient minor allele frequency and/or number of sequences).
Then the best score is used for each input sequence when calculating the edit distance cost.
This can be further split into a multi-consensus by repeatedly running dual-consensus in a binary-tree like system (e.g., repeatedly split the sequences into groups until they do not want to split further).
Finally, priority consensus is just multi-consensus on a chain of sequence inputs instead of a single sequence input.

## Limitations
`waffle_con` has been designed for the specific purpose of PGx consensus in StarPhase using long, accurate HiFi reads.
The underlying algorithms rely on a basic edit-distance wavefront algorithm, which scales with the total edit distance between two sequences.
Thus, high error or high divergence sequences are more likely to lead to longer run times.
Additionally, high error may cause the traversal algorithm to "get lost" in the search space due to weak consensus, which may ultimatley lead to lower quality consensus sequences.
Additionally, best results are with full-length input sequences.
Partial sequences require estimating start/stop positions, which injects more opportunities for error.

## Support information
The `waffle_con` crate is a pre-release software intended for research use **only** and not for use in diagnostic procedures. 
While efforts were made to ensure that `waffle_con` lives up to the quality that PacBio strives for, we make no warranty regarding this software.

As `waffle_con` is **not** covered by any service level agreement or the like, please do not contact a PacBio Field Applications Scientists or PacBio Customer Service for assistance with any `waffle_con` release. 
Please report all issues through GitHub instead. 
We make no warranty that any such issue will be addressed, to any extent or within any time frame.

### DISCLAIMER
THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.
