# chromatic 0.0.03 (03-10-2025)
## New features
- **Reference group Z-scoring**: `RunChromatic()` and `NormalizeStateMatrix()` now support baseline Z-scoring 
  by specifying `z_group_by` and `z_group_name`. This allows erosion scores to be interpreted relative to a chosen reference group 
  (e.g. control samples).
- **Flexible peak annotation**: Added `skip_annotation` option in `RunChromatic()` to allow users to provide 
  pre-annotated peaks with ChromHMM state labels, bypassing the annotation step.
- **Minimum overlap fraction**: `AnnotatePeaks()` and `RunChromatic()` now support filtering assignments 
  with `min_overlap_frac`, requiring a peak to overlap a ChromHMM state by at least a fraction of its length 
  (default = 0.25).
- **Direct peak input**: `RunChromatic()` accepts an optional `peaks_gr` argument, allowing users to provide 
  a custom peak set instead of extracting peaks from the Seurat object.

## Improvements
- Expanded output: `RunChromatic()` now returns multiple normalized versions of the chromatin state matrix 
  (`state_counts`, `state_frac`, `state_CLR`, `state_z`) for downstream analysis and benchmarking.
- Improved input validation in `AnnotatePeaks()` and `RunChromatic()` (checks for valid metadata columns, overlap thresholds, etc.).
- More informative progress messages for annotation, filtering, and scoring steps.

## Bug fixes
- Fixed handling of peaks with no overlaps in `AnnotatePeaks()` (previously dropped silently, now controlled by `keep_unannotated`).
- Corrected `ErosionScore()` to ensure consistent use of baseline scaling when reference group provided.


# chromatic 0.0.02 (27-09-2025)
## Added
- The chromatic tutorial in a human scATAC-seq dataset.

## Changes
- None
