

#' CalculateEpigenomeScores
#'
#' A wrapper function to compute chromatin-state–based scores from scATAC-seq data 
#' stored in a Seurat object. This function:
#' (1) annotates peaks with ChromHMM states,
#' (2) filters peaks and sets up a chromatin state matrix per cell,
#' (3) calculates an erosion score (active vs. repressive balance), and
#' (4) calculates an entropy-based plasticity score.
#'
#' @param seurat_obj Seurat object containing scATAC-seq data. Must have peaks as 
#'   features (rows) and cells as columns in the `assay` specified.
#' @param chromHMM_states A \code{GenomicRanges} object containing ChromHMM state annotations. 
#'   Should have a column with state labels specified by \code{state_col}.
#' @param stoplist (Optional) A \code{GenomicRanges} object of regions to exclude (e.g. blacklisted regions).
#'   If provided, peaks overlapping these regions will be removed.
#' @param remove_nonstandard_chromosomes Logical; if TRUE (default), non-standard chromosomes 
#'   (e.g. scaffolds) will be removed from both peaks and chromHMM annotations.
#' @param filter_features Logical; if TRUE (default), peaks are filtered to exclude features 
#'   with low coverage using \code{ExcludeUncommonPeaks}.
#' @param min_cells Integer; minimum number of cells required for a peak to be kept 
#'   (passed to \code{ExcludeUncommonPeaks}).
#' @param min_counts Integer; minimum total counts required for a peak to be kept 
#'   (passed to \code{ExcludeUncommonPeaks}).
#' @param state_signs Named vector indicating the “sign” (active/repressive) of each chromatin state.
#'   If NULL (default), will be generated automatically from \code{chromHMM_states} using 
#'   \code{ChromatinStateSigns()} and the patterns provided.
#' @param covariates (Optional) Character vector of covariate column names from 
#'   \code{seurat_obj@meta.data} to regress out from the scores (e.g. TSS.enrichment, nCount_ATAC).
#' @param state_col Character; name of the metadata column in \code{chromHMM_states} containing the state label.
#' @param active_patterns Character vector of patterns used to identify active states 
#'   (passed to \code{ChromatinStateSigns()}). Default includes "TssA", "TssFlnk", "Tx", "EnhA", "EnhG", "EnhWk".
#' @param repressive_patterns Character vector of patterns used to identify repressive states 
#'   (passed to \code{ChromatinStateSigns()}). Default includes "ReprPC", "Quies", "Het".
#' @param pseudocount Numeric; pseudocount to add before calculating fractions 
#'   (used by scoring functions). Default = 0.5.
#' @param assay Character; name of the Seurat assay containing scATAC data. Default = 'ATAC'.
#'
#' @details
#' This function orchestrates several steps:
#' - Annotates each ATAC peak with its overlapping ChromHMM state.
#' - Filters peaks to remove blacklisted/nonstandard regions and optionally low-coverage peaks.
#' - Constructs a cell-by-state matrix (fraction of accessibility per state per cell).
#' - Computes an erosion score (active vs. repressive chromatin balance).
#' - Computes an entropy-based plasticity score (degree of state heterogeneity per cell).
#'
#' The helper functions \code{AnnotatePeaks}, \code{ExcludeUncommonPeaks}, 
#' \code{ChromatinStateSigns}, \code{CalculateStateMatrix}, \code{ErosionScore}, and 
#' \code{EntropyScore} must be defined and return objects in the expected format.
#'
#' @return A list containing:
#' \describe{
#'   \item{peaks_gr}{\code{GenomicRanges} of the filtered and annotated peaks.}
#'   \item{state_matrix}{Matrix of chromatin-state fractions per cell (rows = cells, cols = states).}
#'   \item{erosion}{Data frame of erosion scores per cell (optionally covariate-regressed).}
#'   \item{entropy}{Data frame of entropy-based plasticity scores per cell (optionally covariate-regressed).}
#' }
#'
#' @examples
#' \dontrun{
#' scores <- CalculateEpigenomeScores(
#'     seurat_obj = atac_seurat,
#'     chromHMM_states = chromHMM_mouse,
#'     stoplist = blacklist_gr,
#'     covariates = c("TSS.enrichment","nCount_ATAC")
#' )
#' }
#'
#' @importFrom dplyr %>%
#' @importFrom IRanges subsetByOverlaps
#' @importFrom GenomeInfoDb keepStandardChromosomes, seqnames
#' @export
CalculateEpigenomeScores <- function(
    seurat_obj,
    chromHMM_states, 
    stoplist = NULL,
    remove_nonstandard_chromosomes = TRUE,
    filter_features = TRUE,
    min_cells = 100,
    min_counts = 100,
    state_signs = NULL,
    covariates = NULL,
    state_col = "name",
    active_patterns = c("TssA", "TssFlnk", "Tx", "EnhA", "EnhG", "EnhWk"),
    repressive_patterns = c("ReprPC", "Quies", "Het"),
    pseudocount = 0.5,
    assay = 'ATAC'
){

    # 0. Checks 

    # check that covariates exist in meta
    missing_covars <- setdiff(covariates, colnames(seurat_obj@meta.data))
    if(length(missing_covars) > 0){
        stop(paste("The following covariates are missing from meta:", 
                    paste(missing_covars, collapse=", ")))
    }

    #---------------------------------------------------------------
    # 1. Pre-processing
    #---------------------------------------------------------------
    
    DefaultAssay(seurat_obj) <- assay
    peaks_gr <- Signac::granges(seurat_obj)

    # filter peaks by stoplist
    if(!is.null(stoplist)){
        
        print("Filtering by stoplist regions")

        # TODO: Check that stoplist is valid format 
        # remove stoplist regions from peakset
        peaks_gr <- IRanges::subsetByOverlaps(peaks_gr, stoplist, invert = TRUE)
    }

    # only keep seqnames in common 
    common_seqnames <- as.character(GenomicRanges::intersect(GenomeInfoDb::seqnames(peaks_gr), GenomeInfoDb::seqnames(chromHMM_states)))
    print(common_seqnames)
    # peaks_gr <- subset(peaks_gr, seqnames %in% common_seqnames)
    # chromHMM_states <- subset(chromHMM_states, seqnames %in% common_seqnames)
    peaks_gr <- peaks_gr[GenomeInfoDb::seqnames(peaks_gr) %in% common_seqnames]
    chromHMM_states <- chromHMM_states[GenomeInfoDb::seqnames(chromHMM_states) %in% common_seqnames]

    # remove non-standard chromosomes:
    if(remove_nonstandard_chromosomes){
        print("Filtering nonstandard chromosomes")
        peaks_gr <- GenomeInfoDb::keepStandardChromosomes(peaks_gr)
        chromHMM_states <- GenomeInfoDb::keepStandardChromosomes(chromHMM_states)
    }

    print("Annotating peaks by overlapping with chromatin states")
    
    # TODO: impose a min overlap?
    peaks_gr <- AnnotatePeaks(
        peaks_gr, 
        chromHMM_states,
        state_col = state_col
    )

    # get the names of the peaks 
    # TODO: this should be changed based on how the Signac obj is setup
    peaks_names <- paste0(GenomeInfoDb::seqnames(peaks_gr), "-", BiocGenerics::start(peaks_gr), "-", BiocGenerics::end(peaks_gr))

    # filter the peaks matrix by peaks that we are keeping
    peaks_mat <- GetAssayData(seurat_obj, layer='counts', assay = assay)
    peaks_mat <- peaks_mat[peaks_names,]

    print('peaksmat')
    print(dim(peaks_mat))

    # get the peaks matrix
    if(filter_features){
        print("Excluding uncommon peaks")
        output <- ExcludeUncommonPeaks(
            peaks_gr,
            peaks_mat,
            min_counts = min_counts,
            min_cells = min_cells
        )
        peaks_mat <- output$peaks_mat 
        peaks_gr <- output$peaks_gr

    }

    #---------------------------------------------------------------
    # 2. Set up the matrix
    #---------------------------------------------------------------

    print("State matrix")

    state_matrix <- CalculateStateMatrix(
        peaks_mat = peaks_mat, 
        peaks_gr = peaks_gr
    )

    #---------------------------------------------------------------
    # 3. Calculate Erosion score 
    #---------------------------------------------------------------

    # get the state_signs vector
    if(is.null(state_signs)){
        state_signs <- ChromatinStateSigns(
            chromHMM_states = chromHMM_states,
            state_col = state_col,
            active_patterns = active_patterns,
            repressive_patterns = repressive_patterns
        )
    } 

    print("Calculating erosion score ")

    erosion_df <- ErosionScore(
        input = state_matrix, 
        meta = seurat_obj@meta.data,
        state_signs = state_signs,
        pseudocount = pseudocount,
        covariates = covariates
    )

    #---------------------------------------------------------------
    # 5. Calculate Entropy (Plasticity) score 
    #---------------------------------------------------------------

    print("Calculating entropy score ")

    entropy_df <- EntropyScore(
        input = state_matrix, 
        meta = seurat_obj@meta.data,
        pseudocount = pseudocount,
        covariates = covariates
    )

    output <- list(
        "peaks_gr" = peaks_gr,
        "state_matrix" = state_matrix,
        "erosion" = erosion_df,
        "entropy" = entropy_df
    )

    return(output)

}

#' Annotate Peaks with Chromatin States
#'
#' Assigns each peak in a \code{GRanges} object to the chromatin state with 
#' which it has the largest overlap from a ChromHMM \code{GRanges} annotation.
#' This function is typically used to annotate peaks before computing 
#' state-level metrics (e.g., erosion or entropy scores).
#'
#' @param peaks_gr A \code{GRanges} object of peak regions to annotate.
#' @param chromHMM_states A \code{GRanges} object of ChromHMM state annotations. 
#' Must contain a metadata column (default \code{"name"}) specifying state names.
#' @param state_col Character string specifying the metadata column in 
#' \code{chromHMM_states} containing state labels. Default = \code{"name"}.
#' @param keep_unannotated Logical; if \code{TRUE}, peaks without overlaps 
#' to any ChromHMM state are returned with \code{NA} in their \code{annotation} column.
#' If \code{FALSE} (default), unannotated peaks are dropped.
#' @param min_overlap Numeric; minimum overlap width (in bp) required for a peak 
#' to be assigned to a ChromHMM state. Peaks whose maximum overlap is below this 
#' threshold will be considered unannotated. Default = \code{1} (any overlap).
#'
#' @details
#' For each peak:
#' \enumerate{
#'   \item The function finds all overlaps between peaks and ChromHMM states.
#'   \item It computes the width of each overlap.
#'   \item It assigns to each peak the state with the largest overlap 
#'         (if \code{min_overlap} is met).
#'   \item It stores the overlap width in a new metadata column 
#'         (\code{overlap_width}).
#'   \item By default, peaks without any overlap are dropped 
#'         (set \code{keep_unannotated=TRUE} to retain them).
#' }
#'
#' @return A \code{GRanges} object of peaks with added metadata columns:
#' \itemize{
#'   \item \code{annotation} — the ChromHMM state label.
#'   \item \code{overlap_width} — width of overlap between the peak and 
#'         its assigned ChromHMM state.
#' }
#'
#' @examples
#' \dontrun{
#' annotated_peaks <- AnnotatePeaks(
#'   peaks_gr = peaks_gr,
#'   chromHMM_states = chromHMM_states,
#'   state_col = "name",
#'   keep_unannotated = TRUE,
#'   min_overlap = 50
#' )
#' }
#'
#' @seealso \code{\link{ChromatinStateSigns}}, \code{\link{ErosionScore}}
#' @export
AnnotatePeaks <- function(
    peaks_gr, 
    chromHMM_states,
    state_col = "name",
    keep_unannotated = FALSE,
    min_overlap = 1
){
    # ---- sanity checks ----
    if (!inherits(peaks_gr, "GRanges")) {
        stop("`peaks_gr` must be a GRanges object.")
    }
    if (!inherits(chromHMM_states, "GRanges")) {
        stop("`chromHMM_states` must be a GRanges object.")
    }
    if (!(state_col %in% colnames(mcols(chromHMM_states)))) {
        stop("Column '", state_col, "' not found in chromHMM_states metadata.")
    }
    if (!is.numeric(min_overlap) || length(min_overlap) != 1 || min_overlap < 0) {
        stop("`min_overlap` must be a single non-negative numeric value.")
    }

    # find overlaps between peaks and chromHMM states
    ov <- findOverlaps(peaks_gr, chromHMM_states)
    if (length(ov) == 0) {
        warning("No overlaps found between peaks and chromHMM states. Returning input peaks_gr unchanged.")
        peaks_gr$annotation <- NA_character_
        peaks_gr$overlap_width <- NA_integer_
        return(peaks_gr)
    }

    # compute overlap widths
    ov_df <- as.data.frame(ov)
    ov_df$width <- width(
        GenomicRanges::pintersect(peaks_gr[ov_df$queryHits],
                   chromHMM_states[ov_df$subjectHits])
    )

    # for each peak, choose the overlap with maximum width
    ov_best <- ov_df %>%
        dplyr::group_by(.data$queryHits) %>%
        dplyr::slice_max(.data$width, n = 1, with_ties = FALSE) %>%
        dplyr::ungroup()

    # apply minimum overlap filter: peaks below threshold treated as unannotated
    ov_best <- ov_best[ov_best$width >= min_overlap, ]

    # map peaks to state labels and overlap widths
    peak_to_state <- rep(NA_character_, length(peaks_gr))
    peak_to_width <- rep(NA_integer_, length(peaks_gr))

    peak_to_state[ov_best$queryHits] <- 
        mcols(chromHMM_states)[ov_best$subjectHits, state_col]
    peak_to_width[ov_best$queryHits] <- ov_best$width

    # annotate peaks
    peaks_gr$annotation <- peak_to_state
    peaks_gr$overlap_width <- peak_to_width

    # filter out peaks without annotation if requested
    if (!keep_unannotated) {
        peaks_gr <- peaks_gr[!is.na(peaks_gr$annotation)]
    }

    # optional: ensure peaks still overlap chromHMM_states (final cleanup)
    if (!keep_unannotated) {
        peaks_gr <- subsetByOverlaps(peaks_gr, chromHMM_states)
    }

    return(peaks_gr)
}

#' Exclude Lowly Detected or Low-Count Peaks
#'
#' Filters out peaks that are rarely observed across cells or have 
#' low total counts in a peak-by-cell matrix. This is a simple feature 
#' selection step prior to downstream analyses (e.g., entropy or erosion scores).
#'
#' @param peaks_gr A \code{GRanges} object containing the genomic coordinates 
#' of peaks. Rows should correspond to peaks in \code{peaks_mat}.
#' @param peaks_mat A numeric matrix of peak counts (rows = peaks, columns = cells).
#' @param min_cells Integer; minimum number of cells in which a peak must be detected 
#' (nonzero counts) to be retained. Default = \code{100}.
#' @param min_counts Integer; minimum total counts across all cells for a peak 
#' to be retained. Default = \code{100}.
#' @param verbose Logical; if \code{TRUE}, print a message summarizing how many peaks 
#' were retained/filtered. Default = \code{TRUE}.
#'
#' @details
#' The function applies two filters:
#' \enumerate{
#'   \item A peak must be present (nonzero) in at least 
#'         \code{min_cells} cells.
#'   \item A peak must have at least \code{min_counts} total counts 
#'         across all cells.
#' }
#' Peaks passing both filters are retained in the output. This step helps reduce 
#' noise and memory usage in large peak-by-cell matrices.
#'
#' @return A named \code{list} with:
#' \itemize{
#'   \item \code{peaks_gr} — the filtered \code{GRanges} object of peaks.
#'   \item \code{peaks_mat} — the filtered peak-by-cell count matrix.
#' }
#'
#' @examples
#' \dontrun{
#' filtered <- ExcludeUncommonPeaks(
#'   peaks_gr = peaks_gr,
#'   peaks_mat = peaks_mat,
#'   min_cells = 200,
#'   min_counts = 500,
#'   verbose = TRUE
#' )
#' filtered$peaks_gr
#' filtered$peaks_mat
#' }
#'
#' @seealso \code{\link{AnnotatePeaks}}
#' @export
ExcludeUncommonPeaks <- function(
    peaks_gr,
    peaks_mat,
    min_cells = 100,
    min_counts = 100,
    verbose = TRUE
){
    # ---- sanity checks ----
    if (!inherits(peaks_gr, "GRanges")) {
        stop("`peaks_gr` must be a GRanges object.")
    }
    if (!inherits(peaks_mat, 'Matrix')) {
        stop("`peaks_mat` must be a matrix (rows = peaks, cols = cells).")
    }
    if (nrow(peaks_mat) != length(peaks_gr)) {
        stop("Number of rows in `peaks_mat` must match length of `peaks_gr`.")
    }
    if (!is.numeric(min_cells) || length(min_cells) != 1 || min_cells < 0) {
        stop("`min_cells` must be a single non-negative numeric value.")
    }
    if (!is.numeric(min_counts) || length(min_counts) != 1 || min_counts < 0) {
        stop("`min_counts` must be a single non-negative numeric value.")
    }
    if (!is.logical(verbose) || length(verbose) != 1) {
        stop("`verbose` must be a single logical (TRUE/FALSE).")
    }

    # ---- filtering ----
    total_peaks <- nrow(peaks_mat)

    # total counts per peak
    peak_counts <- rowSums(peaks_mat)

    # number of cells with nonzero counts per peak
    cells_per_peak <- rowSums(peaks_mat > 0)

    # indices of peaks meeting both thresholds
    keep_peaks <- intersect(
        which(cells_per_peak >= min_cells),
        which(peak_counts >= min_counts)
    )

    # filter GRanges and matrix
    peaks_gr <- peaks_gr[keep_peaks]
    peaks_mat <- peaks_mat[keep_peaks, , drop = FALSE]

    if (verbose) {
        message(sprintf(
            "Filtered peaks: retained %d of %d peaks (%.1f%%).",
            length(keep_peaks),
            total_peaks,
            100 * length(keep_peaks) / total_peaks
        ))
    }

    # return list
    return(list(
        'peaks_gr' = peaks_gr,
        'peaks_mat' = peaks_mat
    ))
}


#' Assign Active/Repressive Signs to Chromatin States
#'
#' Generates a named vector of sign values (+1 = repressive, -1 = active, 0 = unclassified) 
#' for chromatin states based on user-specified patterns. This is typically used to 
#' weight states in functions such as \code{\link{ErosionScore}}.
#'
#' @param chromHMM_states A \code{GRanges} or similar object containing chromatin state 
#' annotations with a metadata column specifying state names.
#' @param state_col Character string specifying the metadata column name in 
#' \code{chromHMM_states@elementMetadata} containing the state names. Default = \code{"name"}.
#' @param active_patterns Character vector of regex patterns used to identify active states. 
#' Default = \code{c("TssA","TssFlnk","Tx","EnhA","EnhG","EnhWk")}.
#' @param repressive_patterns Character vector of regex patterns used to identify 
#' repressive states. Default = \code{c("ReprPC","Quies","Het")}.
#' @param error_if_unclassified Logical indicating whether to stop with an error if any 
#' states remain unclassified (sign = 0). Default = \code{FALSE} (issues a message instead).
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Extracts unique chromatin state names from \code{chromHMM_states}.
#'   \item Initializes all state signs to 0 (unclassified).
#'   \item Assigns \code{-1} for states matching any \code{active_patterns}.
#'   \item Assigns \code{+1} for states matching any \code{repressive_patterns}.
#' }
#'
#' The returned vector can be supplied to \code{\link{ErosionScore}} 
#' as the \code{state_signs} argument.
#'
#' @return A named numeric vector of signs with one entry per unique chromatin state:
#' \itemize{
#'   \item \code{-1}: Active state
#'   \item \code{+1}: Repressive state
#'   \item \code{0}: Unclassified (did not match any pattern)
#' }
#'
#' @examples
#' \dontrun{
#' state_signs <- ChromatinStateSigns(
#'   chromHMM_states = chromHMM_gr,
#'   state_col = "name",
#'   error_if_unclassified = TRUE
#' )
#' }
#'
#' @seealso \code{\link{ErosionScore}}
#' @export
ChromatinStateSigns <- function(
    chromHMM_states,
    state_col = "name",
    active_patterns = c("TssA", "TssFlnk", "Tx", "EnhA", "EnhG", "EnhWk"),
    repressive_patterns = c("ReprPC", "Quies", "Het"),
    error_if_unclassified = FALSE
){

    # extract all states
    all_states <- unique(as.character(chromHMM_states@elementMetadata[, state_col]))

    # initialize the sign vector
    state_signs <- rep(0, length(all_states))
    names(state_signs) <- all_states

    # assign active states (-1)
    for (pat in active_patterns) {
        state_signs[grepl(pat, all_states, ignore.case = TRUE)] <- -1
    }

    # assign repressive states (+1)
    for (pat in repressive_patterns) {
        state_signs[grepl(pat, all_states, ignore.case = TRUE)] <- +1
    }

    # sanity checks
    if (all(state_signs == 0)) {
        stop("No states matched any active or repressive patterns.")
    }

    unclassified_states <- names(state_signs)[state_signs == 0]
    if (length(unclassified_states) > 0) {
        if (error_if_unclassified) {
            stop("The following states were not classified: ",
                 paste(unclassified_states, collapse = ", "))
        } else {
            message("The following states were not classified (sign=0): ",
                    paste(unclassified_states, collapse = ", "))
        }
    }

    return(state_signs)
}



#' Assign Active/Repressive Signs to Chromatin States
#'
#' Generates a named vector of sign values (+1 = repressive, -1 = active, 0 = unclassified) 
#' for chromatin states based on user-specified patterns. This is typically used to 
#' weight states in functions such as \code{\link{ErosionScore}}.
#'
#' @param chromHMM_states A \code{GRanges} or similar object containing chromatin state 
#' annotations with a metadata column specifying state names.
#' @param state_col Character string specifying the metadata column name in 
#' \code{chromHMM_states@elementMetadata} containing the state names. Default = \code{"name"}.
#' @param active_patterns Character vector of regex patterns used to identify active states. 
#' Default = \code{c("TssA","TssFlnk","Tx","EnhA","EnhG","EnhWk")}.
#' @param repressive_patterns Character vector of regex patterns used to identify 
#' repressive states. Default = \code{c("ReprPC","Quies","Het")}.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Extracts unique chromatin state names from \code{chromHMM_states}.
#'   \item Initializes all state signs to 0 (unclassified).
#'   \item Assigns \code{-1} for states matching any \code{active_patterns}.
#'   \item Assigns \code{+1} for states matching any \code{repressive_patterns}.
#' }
#'
#' The returned vector can be supplied to \code{\link{ErosionScore}} 
#' as the \code{state_signs} argument.
#'
#' @return A named numeric vector of signs with one entry per unique chromatin state:
#' \itemize{
#'   \item \code{-1}: Active state
#'   \item \code{+1}: Repressive state
#'   \item \code{0}: Unclassified (did not match any pattern)
#' }
#'
#' @examples
#' \dontrun{
#' state_signs <- ChromatinStateSigns(
#'   chromHMM_states = chromHMM_gr,
#'   state_col = "name"
#' )
#' }
#'
#' @seealso \code{\link{ErosionScore}}
#' @export
ChromatinStateSigns <- function(
    chromHMM_states,
    state_col = "name",
    active_patterns = c("TssA", "TssFlnk", "Tx", "EnhA", "EnhG", "EnhWk"),
    repressive_patterns = c("ReprPC", "Quies", "Het")
){

    # extract all states
    all_states <- unique(as.character(chromHMM_states@elementMetadata[, state_col]))

    # initialize the sign vector
    state_signs <- rep(0, length(all_states))
    names(state_signs) <- all_states

    # assign active states (-1)
    for (pat in active_patterns) {
        state_signs[grepl(pat, all_states, ignore.case = TRUE)] <- -1
    }

    # assign repressive states (+1)
    for (pat in repressive_patterns) {
        state_signs[grepl(pat, all_states, ignore.case = TRUE)] <- +1
    }

    # sanity checks
    if (all(state_signs == 0)) {
        warning("No states matched any active or repressive patterns.")
    }

    unclassified_states <- names(state_signs)[state_signs == 0]
    if (length(unclassified_states) > 0) {
        message("The following states were not classified (sign=0): ",
                paste(unclassified_states, collapse = ", "))
    }

    return(state_signs)
}

#' Calculate per-cell chromatin state counts from annotated peaks
#'
#' This helper function aggregates a peak-by-cell matrix into a 
#' chromatin state–by-cell matrix based on the state annotations 
#' assigned to each peak. Each cell’s counts across peaks belonging 
#' to the same chromatin state are summed to produce a per-state 
#' count matrix, which is then transposed to return a cell-by-state 
#' matrix suitable for downstream analyses.
#'
#' @param peaks_mat A sparse or dense matrix of peak accessibility counts 
#'   (rows = peaks, columns = cells). Typically obtained from 
#'   \code{\link[Seurat]{GetAssayData}} on the ATAC assay.
#' @param peaks_gr A \code{\link[GenomicRanges]{GRanges}} object of the same 
#'   length as the number of rows in \code{peaks_mat}, containing at least 
#'   one column named \code{annotation} giving the chromatin state label 
#'   for each peak.
#'
#' @return A cell-by-state matrix of summed fragment counts per chromatin state 
#'   (rows = cells, columns = chromatin states).
#'
#' @details This function assumes that each row of \code{peaks_mat} corresponds 
#'   to the same peak as the corresponding entry in \code{peaks_gr}, and that 
#'   each peak has a valid chromatin state annotation in \code{peaks_gr$annotation}. 
#'   Peaks with the same state annotation are grouped and their counts summed 
#'   across cells. The resulting state-by-cell matrix is transposed to return 
#'   cell-by-state counts.
#'
#' @examples
#' # peaks_mat: peaks x cells matrix from ATAC assay
#' # peaks_gr: GRanges of peaks with 'annotation' column
#' state_matrix <- CalculateStateMatrix(peaks_mat, peaks_gr)
#'
#' @export
CalculateStateMatrix <- function(
    peaks_mat, 
    peaks_gr 
){ 

    if (nrow(peaks_mat) != length(peaks_gr)) {
        stop("Number of rows in peaks_mat (", nrow(peaks_mat),
             ") does not match length of peaks_gr (", length(peaks_gr), ").")
    }

    if (!"annotation" %in% colnames(mcols(peaks_gr))) {
        stop("peaks_gr must contain a column named 'annotation' with chromatin state labels.")
    }

    # set rownames to the state annotations 
    rownames(peaks_mat) <- as.character(peaks_gr$annotation) 
    
    # Sum across segments within each state: 
    groupings <- rownames(peaks_mat)
    state_matrix <- rowsum(peaks_mat, groupings)
    state_matrix <- t(state_matrix) 
    return(state_matrix) 
}


#' Calculate Erosion Scores from Chromatin State Matrices
#'
#' Computes cell-level epigenomic erosion scores from a chromatin state 
#' count matrix. The function performs fraction normalization, 
#' centered log-ratio (CLR) transformation, z-scoring, and applies a 
#' state-specific sign vector (+1 for repressive states, -1 for active states). 
#' Optionally, it can regress out covariates (e.g., TSS enrichment, 
#' nCount_ATAC) from the resulting erosion scores.
#'
#' @param input A numeric matrix of chromatin state counts with cells as rows 
#' and states as columns. Typically produced by \code{CalculateStateMatrix()}.
#' @param meta A data frame (such as \code{seurat_obj@meta.data}) containing 
#' cell-level metadata. Required if covariate regression is to be performed.
#' @param state_signs A named numeric vector indicating the sign (+1 or -1) 
#' for each chromatin state. Names must match the column names of \code{input}.
#' @param pseudocount A numeric value to add to counts to avoid division by 
#' zero. Default is 0.5.
#' @param covariates A character vector of column names in \code{meta} 
#' specifying which covariates to regress out of the erosion score. Default is \code{NULL}.
#'
#' @details
#' The erosion score is calculated as:
#' \enumerate{
#'   \item Normalize counts per cell to fractions.
#'   \item Apply log transform and center across states.
#'   \item Z-score each state across all cells.
#'   \item Multiply z-scored states by the sign vector (+1 for repressive, -1 for active).
#'   \item Sum across states to obtain one erosion score per cell.
#' }
#'
#' If \code{covariates} is supplied, a linear model is fit using the 
#' specified covariates as predictors of the erosion score. The residuals 
#' from this model are returned as a covariate-corrected erosion score.
#'
#' @return A data frame with:
#' \itemize{
#'   \item Z-scored values for each state per cell (columns = states).
#'   \item \code{erosion_score}: the raw erosion score per cell.
#'   \item \code{erosion_score_corrected} (if covariates given): residual 
#'   erosion score with covariates regressed out.
#' }
#'
#' @examples
#' \dontrun{
#' erosion_df <- ErosionScore(
#'   input = state_matrix,
#'   meta = seurat_obj@meta.data,
#'   state_signs = state_signs,
#'   covariates = c("TSS.enrichment", "nCount_ATAC")
#' )
#' }
#'
#' @seealso \code{\link{EntropyScore}}, \code{\link{CalculateStateMatrix}}
#' @export
ErosionScore <- function(
    input, 
    meta,
    state_signs,
    pseudocount = 0.5,
    covariates = NULL
){

    # calculate fraction
    state_frac <- (input + pseudocount) / rowSums(input + pseudocount)

    # centered log-ratio normalization (CLR)
    log_frac <- log(state_frac)
    log_frac_centered <- log_frac - rowMeans(log_frac)

    # Z-score each state (column) across all cells
    # (center to 0, scale to unit variance per state)
    state_z <- scale(log_frac_centered)  # gives a matrix same dim as log_frac_centered

    # Multiply each state by its sign (+1 for repressive, -1 for active)
    signed_z <- sweep(state_z, 2, state_signs[colnames(state_z)], "*")

    # Sum across states to get one erosion/entropy score per cell
    erosion_score <- rowSums(signed_z)

    out_df <- as.data.frame(signed_z)
    out_df$erosion_score <- erosion_score

    if(!is.null(covariates)){

        # build data frame for regression
        reg_df <- meta[, covariates, drop=FALSE]
        reg_df$erosion_score <- erosion_score

        # build formula: erosion_score ~ covar1 + covar2 + ...
        form <- as.formula(
            paste("erosion_score ~", paste(covariates, collapse = " + "))
        )

        # fit model
        model <- lm(form, data = reg_df)

        # store residuals as covariate-corrected erosion score
        out_df$erosion_score_corrected <- resid(model)
    }
    
    return(out_df)

}

#' Calculate Cell-Level Chromatin State Entropy (Plasticity) Scores
#'
#' Computes per-cell entropy scores from a chromatin state count matrix.  
#' This function normalizes counts to fractions and then calculates 
#' Shannon entropy for each cell to quantify the diversity of chromatin 
#' states (a proxy for epigenomic plasticity). Optionally, it can regress 
#' out covariates (e.g. sequencing depth, TSS enrichment) to return 
#' a covariate-corrected entropy score.
#'
#' @param input A numeric matrix of chromatin state counts 
#' (rows = cells, columns = states). Typically generated by \code{CalculateStateMatrix()}.
#' @param meta A data frame of cell-level metadata (e.g. \code{seurat_obj@meta.data}). 
#' Required if \code{covariates} is provided.
#' @param pseudocount Numeric pseudocount to add to all counts before normalization 
#' to avoid zeros. Default = 0.5.
#' @param covariates Optional character vector of column names in \code{meta} 
#' specifying covariates to regress out of the entropy score. Default = \code{NULL}.
#'
#' @details
#' For each cell:
#' \enumerate{
#'   \item Convert counts to per-cell fractions with pseudocount.
#'   \item Compute Shannon entropy of the fraction vector 
#'         (\code{entropy::entropy(p, unit="log2")}).
#'   \item Normalize the entropy score by dividing by \code{log2(n_states)} 
#'         (maximum possible entropy) to scale it between 0 and 1.
#' }
#'
#' If \code{covariates} is supplied, a linear model is fit with the specified 
#' covariates as predictors of the raw entropy score. The residuals from this model 
#' are returned as the covariate-corrected entropy score.
#'
#' @return A data frame with:
#' \itemize{
#'   \item \code{entropy}: raw Shannon entropy score per cell (in bits).
#'   \item \code{entropy_score_norm}: entropy normalized to 0–1 range by 
#'   dividing by \code{log2(n_states)}.
#'   \item \code{entropy_score_corrected} (if \code{covariates} given): 
#'   residual entropy score with covariates regressed out.
#' }
#'
#' @examples
#' \dontrun{
#' entropy_df <- EntropyScore(
#'   input = state_matrix,
#'   meta = seurat_obj@meta.data,
#'   covariates = c("TSS.enrichment", "nCount_ATAC")
#' )
#' }
#'
#' @seealso \code{\link{ErosionScore}}, \code{\link{CalculateStateMatrix}}
#' @export
EntropyScore <- function(
    input, 
    meta,
    pseudocount = 0.5,
    covariates = NULL
){

    # calculate fraction
    state_frac <- (input + pseudocount) / rowSums(input + pseudocount)

    # calculate entropy scores per cell
    entropy_score <- apply(
        input, 
        1, 
        function(p) {
            entropy::entropy(p, unit = "log2")  # Shannon entropy in bits
        }
    )
    
    n_states <- ncol(input)
    entropy_score_norm <- entropy_score / log2(n_states)

    out_df <- data.frame(
        "entropy" = entropy_score,
        "entropy_score_norm" = entropy_score_norm
    )

  if(!is.null(covariates)){

        # build data frame for regression
        reg_df <- meta[, covariates, drop=FALSE]
        reg_df$entropy_score <- entropy_score  # regress the raw entropy score

        # build formula: entropy_score ~ covar1 + covar2 + ...
        form <- as.formula(
            paste("entropy_score ~", paste(covariates, collapse = " + "))
        )

        # fit model
        model <- lm(form, data = reg_df)

        # store residuals as covariate-corrected entropy score
        out_df$entropy_score_corrected <- resid(model)
    }

    return(out_df)

}
