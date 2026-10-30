# Changelog

All notable changes to Biotools will be noted here.

## [0.2.0] - 2026-10-30

### Added
- Accelerated protein MSA using a Rust backend, substantially improving performance over the Python implementation (Python functions still available)
- Protein MSA can be tuned with custom affine --gap-open and --gap-extend penalties
- Introduced _non_positive_int to validate gap penalty arguments
- Argument --backend added to MSA, with choices "rust" and "python" (Rust defaults) 

### Changed
- Bundled JetBrains Mono Regular Font for cross-platform image generation, replacing the macOS-specific Menlo font dependency
- Updated the README

### Fixed
- Type _positive_int was updated so that floating point arguments return a helpful error message
- Changed Rust alignment scoring from f32 to f64 to improve numerical consistency with Python
- Added the proper repo URL to the README install instructions

### Testing
- Added Hypothesis-based randomized testing to validate consistency between the Rust and Python MSA implementation
- Started a testing module for the image renderer


## [0.1.0] - 2026-10-20

### Added

#### DNA sequence tools
- Transcribe to RNA
- Translate in three reading frames
- Generate reverse complement
- Calculate GC content
- Pairwise sequence alignment

#### Protein sequence tools
- Calculate net charge
- Calculate molecular weight
- Calculate extinction coefficient
- Search for and highlight N-linked glycosylation sites
- Compare two sequences and look for conserved motifs of varying lengths
- Pairwise sequence alignment (default implementation in Rust, Python functions available)
- Multiple sequence alignment (MSA) from FASTA files or FASTA-Formatted strings

#### General sequence tools
- Search for sequences and report their positions
- Returns the length of a sequence
- Reports the percent identity between all sequences from a FASTA file

#### Rendering
- Alignments output to terminal (when no --save argument provided)
- Generate publication-quality PNG and TIFF alignment images
    - Save Pairwise DNA and protein alignments
    - Save Multiple Sequence Alignment with proteins
    - Alignment row length can be tuned with --wrap
    - Argument --real-ruler displays amino acid residue number based on the first sequence, skipping gaps
- Argument --diff displays alignment differences relative to a reference sequence
- MSA optionally adds a consensus sequence with --consensus

#### Testing
- Comprehensive unit testing of all Python algorithm modules
- End to end testing of the MSA pipeline, verifying preservation of original sequences, uniform alignment lengths, and consistency between the Rust and Python pairwise_protein implementations