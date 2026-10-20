/// Fixed options for a [`DWFALite`].
/// These values are chosen when the DWFA is constructed and do not change as it is extended.
///
/// ```
/// use waffle_con::dynamic_wfa::{DWFALite, DWFALiteConfigBuilder};
/// let config = DWFALiteConfigBuilder::default()
///     .wildcard(Some(b'N'))
///     .allow_early_termination(true)
///     .max_edit_distance(Some(500))
///     .build()
///     .unwrap();
/// let dwfa = DWFALite::new(config);
/// assert_eq!(dwfa.max_edit_distance(), Some(500));
/// ```
#[derive(derive_builder::Builder, Clone, Debug, Eq, Hash, PartialEq)]
#[builder(default)]
pub struct DWFALiteConfig {
    /// Optional wildcard symbol that matches anything.
    pub wildcard: Option<u8>,
    /// When true, reaching the end of the baseline does not penalize a longer other sequence.
    pub allow_early_termination: bool,
    /// Absolute edit-distance cap. `None` means unlimited.
    pub max_edit_distance: Option<usize>,
}

// explicitly set our default for now, even if it's derivable
#[allow(clippy::derivable_impls)]
impl Default for DWFALiteConfig {
    fn default() -> Self {
        DWFALiteConfig {
            wildcard: None,
            allow_early_termination: false,
            max_edit_distance: None,
        }
    }
}
