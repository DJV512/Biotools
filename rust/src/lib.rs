mod constants;
mod pairwise;

pub use pairwise::pairwise_protein;

#[cfg(feature = "python")]
use pyo3::prelude::*;

/// Globally align two protein sequences using affine gap penalties.
///
/// Inputs must be normalized sequences from validated Protein objects.
/// This internal function does not perform input validation.
///
/// Returns the two aligned sequences, with '-' characters marking gaps.
#[cfg(feature = "python")]
#[pyfunction]
fn _align_protein(seq1: &str, seq2: &str) -> (String, String) {
    pairwise_protein(seq1, seq2)
}

#[cfg(feature = "python")]
#[pymodule]
fn biotools_rust(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(_align_protein, m)?)?;
    Ok(())
}
