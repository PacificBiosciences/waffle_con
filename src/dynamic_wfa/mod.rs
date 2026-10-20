//! A dynamic wavefront alignment of one baseline against one other sequence.
//! The sequences live outside [`DWFALite`]; this module keeps the wavefront.

mod dwfa_config;
mod dwfa_lite;
mod votes;
mod wavefront;

#[cfg(test)]
mod tests;

pub use dwfa_config::{DWFALiteConfig, DWFALiteConfigBuilder};
pub use dwfa_lite::{DWFALite, DWFALiteState};
pub(crate) use votes::{ExtensionVotes, VoteBuffers, VoteScratch};
