#' RunChromatic
#'
#' A wrapper function to compute chromatin-state–based scores from scATAC-seq data 
#' stored in a Seurat object. This function orchestrates the full pipeline:
#' (1) annotates peaks with ChromHMM states (or accepts pre-annotated peaks),
#' (2) filters peaks (stoplist, nonstandard chromosomes, and low-coverage filtering),
#' (3) constructs and normalizes a chromatin state count/fraction matrix per cell,
#' (4) computes an erosion score (active vs. repressive balance),
#' and (5) computes an entropy-based plasticity score.
#'
#' @param seurat_obj Seurat object containing scATAC-seq data. Must have peaks as 
#'   features (rows) and cells as columns in the `assay` specified.
#' @param chromHMM_states A \code{GRanges} object containing ChromHMM state annotations. 
#'   Should have a column with state labels specified by \code{state_col}.
#' @param peaks_gr (Optional) A \code{GRanges} object of peaks. If \code{NULL} (default), 
#'   peaks are extracted from the input Seurat object. If provided, will be used directly.
#' @param stoplist (Optional) A \code{GRanges} object of regions to exclude (e.g. blacklisted regions).
#'   If provided, peaks overlapping these regions will be removed.
#' @param remove_nonstandard_chromosomes Logical; if TRUE (default), non-standard chromosomes 
#'   (e.g. scaffolds) will be removed from both peaks and chromHMM annotations.
#' @param min_overlap Numeric; minimum overlap width (bp) required for a peak to be assigned 
#'   to a ChromHMM state. Default = 50.
#' @param min_overlap_frac Numeric or \code{NULL}; minimum fraction of a peak's length 
#'   that must overlap a ChromHMM state for assignment. Default = 0.25.
#' @param filter_features Logical; if TRUE (default), peaks are filtered to exclude features 
#'   with low coverage using \code{ExcludeUncommonPeaks}.
#' @param skip_annotation Logical; if TRUE, skips ChromHMM annotation step and assumes that 
#'   the provided \code{peaks_gr} (or peaks from \code{seurat_obj}) already contain 
#'   an \code{annotation} column. Default = FALSE.
#' @param min_cells Integer; minimum number of cells required for a peak to be kept 
#'   (passed to \code{ExcludeUncommonPeaks}). Default = 100.
#' @param min_counts Integer; minimum total counts required for a peak to be kept 
#'   (passed to \code{ExcludeUncommonPeaks}). Default = 100.
#' @param state_signs Named vector indicating the “sign” (active/repressive) of each chromatin state.
#'   If NULL (default), will be generated automatically from \code{chromHMM_states} using 
#'   \code{ChromatinStateSigns()} and the provided patterns.
#' @param covariates (Optional) Character vector of covariate column names from 
#'   \code{seurat_obj@meta.data} to regress out from the scores 
#'   (e.g. TSS.enrichment, nCount_ATAC).
#' @param state_col Character; name of the metadata column in \code{chromHMM_states} 
#'   containing the state label. Default = \code{"name"}.
#' @param active_patterns Character vector of regex patterns used to identify active states 
#'   (passed to \code{ChromatinStateSigns()}). 
#'   Default includes "TssA", "TssFlnk", "Tx", "EnhA", "EnhG", "EnhWk".
#' @param repressive_patterns Character vector of regex patterns used to identify 
#'   repressive states (passed to \code{ChromatinStateSigns()}). 
#'   Default includes "ReprPC", "Quies", "Het".
#' @param pseudocount Numeric; pseudocount to add before calculating fractions 
#'   (used by scoring functions). Default = 0.5.
#' @param z_group_by (Optional) Column name in \code{seurat_obj@meta.data} specifying a group 
#'   to use as reference for baseline Z-scoring. If \code{NULL}, Z-scores are computed 
#'   across all cells.
#' @param z_group_name (Optional) Name of the group (within \code{z_group_by}) to use as 
#'   reference for baseline Z-scoring.
#' @param assay Character; name of the Seurat assay containing scATAC data. Default = 'ATAC'.
#'
#' @details
#' The workflow consists of:
#' \enumerate{
#'   \item Peak annotation with ChromHMM states (unless \code{skip_annotation=TRUE}).
#'   \item Filtering peaks using stoplists, standard chromosomes, and minimum overlap rules.
#'   \item Feature-level filtering for low-coverage peaks.
#'   \item Construction of a cell-by-state matrix and normalization via fractions, CLR, 
#'         and Z-scores (with optional reference group scaling).
#'   \item Calculation of the erosion score (balance between active and repressive states).
#'   \item Calculation of the entropy score (plasticity of chromatin states per cell).
#' }
#'
#' @return A list containing:
#' \describe{
#'   \item{peaks_gr}{\code{GenomicRanges} of the filtered and annotated peaks.}
#'   \item{state_counts}{Matrix of raw chromatin-state counts per cell.}
#'   \item{state_frac}{Matrix of state fractions per cell.}
#'   \item{state_CLR}{Matrix of CLR-normalized state values.}
#'   \item{state_z}{Matrix of Z-scored CLR values (optionally baseline-scaled).}
#'   \item{scores}{Data frame containing erosion and entropy scores per cell 
#'   (optionally covariate-regressed).}
#' }
#'
#' @examples
#' \dontrun{
#' output <- RunChromatic(
#'     seurat_obj = seurat_obj,
#'     chromHMM_states = chromHMM_states,
#' )
#' }
#'
#' @seealso \code{\link{AnnotatePeaks}}, \code{\link{NormalizeStateMatrix}}, 
#'   \code{\link{ErosionScore}}, \code{\link{EntropyScore}}
#' @export
RunChromatic <- function(
    seurat_obj,
    chromHMM_states,
    peaks_gr = NULL, 
    stoplist = NULL,
    remove_nonstandard_chromosomes = TRUE,
    min_overlap = 50,
    min_overlap_frac = 0.25,
    filter_features = TRUE,
    skip_annotation = FALSE,
    min_cells = 100,
    min_counts = 100,
    state_signs = NULL,
    covariates = NULL,
    state_col = "name",
    active_patterns = c("TssA", "TssFlnk", "Tx", "EnhA", "EnhG", "EnhWk"),
    repressive_patterns = c("ReprPC", "Quies", "Het"),
    pseudocount = 0.5,
    z_group_by = NULL,
    z_group_name = NULL,
    assay = 'ATAC'
){

    # 0. Checks 

    # check that covariates exist in meta
    missing_covars <- setdiff(covariates, colnames(seurat_obj@meta.data))
    if(length(missing_covars) > 0){
        stop(paste("The following covariates are missing from seurat_obj@meta.data slot:", 
                    paste(missing_covars, collapse=", ")))
    }

    #---------------------------------------------------------------
    # 1. Pre-processing
    #---------------------------------------------------------------
    
    DefaultAssay(seurat_obj) <- assay

    if(is.null(peaks_gr)){
        peaks_gr <- Signac::granges(seurat_obj)
    }

    if(skip_annotation){
        if(!("annotation" %in% names(mcols(peaks_gr)))){
            stop("Chromatin state annotations missing from granges(seurat_obj). Please run again with skip_annotation=FALSE.")
        }
    } else{

        # filter peaks by stoplist
        if(!is.null(stoplist)){
            print("Filtering by stoplist regions")

            # TODO: Check that stoplist is valid format 
            # remove stoplist regions from peakset
            peaks_gr <- IRanges::subsetByOverlaps(peaks_gr, stoplist, invert = TRUE)
        }

        # only keep seqnames in common 
        common_seqnames <- as.character(GenomicRanges::intersect(GenomeInfoDb::seqnames(peaks_gr), GenomeInfoDb::seqnames(chromHMM_states)))
        peaks_gr <- peaks_gr[as.character(GenomeInfoDb::seqnames(peaks_gr)) %in% common_seqnames]
        chromHMM_states <- chromHMM_states[as.character(GenomeInfoDb::seqnames(chromHMM_states)) %in% common_seqnames]

        # remove non-standard chromosomes:
        if(remove_nonstandard_chromosomes){
            print("Filtering nonstandard chromosomes")
            peaks_gr <- GenomeInfoDb::keepStandardChromosomes(peaks_gr)
            chromHMM_states <- GenomeInfoDb::keepStandardChromosomes(chromHMM_states)
        }

        print("Annotating peaks by overlapping with chromatin states")
        peaks_gr <- AnnotatePeaks(
            peaks_gr, 
            chromHMM_states,
            state_col = state_col,
            min_overlap = min_overlap,
            min_overlap_frac = min_overlap_frac,
            verbose = TRUE
        )

    }

    # TODO: this should be changed based on how the Signac obj is setup
    # Could detect this automatically
    peaks_names <- paste0(GenomeInfoDb::seqnames(peaks_gr), "-", BiocGenerics::start(peaks_gr), "-", BiocGenerics::end(peaks_gr))

    # filter the peaks matrix by peaks that we are keeping
    peaks_mat <- GetAssayData(seurat_obj, layer='counts', assay = assay)
    peaks_mat <- peaks_mat[peaks_names,]

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

    print("Calculating chromatin state matrix")
    state_mat <- CalculateStateMatrix(
        peaks_mat = peaks_mat, 
        peaks_gr = peaks_gr
    )

    print("Normalizing chromatin state matrix")
    mat_list <- NormalizeStateMatrix(
        state_mat = state_mat,
        meta = seurat_obj@meta.data,
        pseudocount = pseudocount,
        group_by = z_group_by,
        group_name = z_group_name
    )

    #---------------------------------------------------------------
    # 3. Calculate Erosion score 
    #---------------------------------------------------------------

    state_names <- unique(as.character(peaks_gr@elementMetadata[, 'annotation']))

    # get the state_signs vector
    if(is.null(state_signs)){
        state_signs <- ChromatinStateSigns(
            state_names = state_names,
            active_patterns = active_patterns,
            repressive_patterns = repressive_patterns
        )
    } 

    print("Calculating erosion score ")

    erosion_df <- ErosionScore(
        mat = mat_list$zscore, 
        meta = seurat_obj@meta.data,
        state_signs = state_signs,
        covariates = covariates
    )

    #---------------------------------------------------------------
    # 5. Calculate Entropy (Plasticity) score 
    #---------------------------------------------------------------

    print("Calculating entropy score ")

    entropy_df <- EntropyScore(
        mat = mat_list$frac, 
        meta = seurat_obj@meta.data,
        covariates = covariates
    )

    #---------------------------------------------------------------
    # 6. Set up the outputs 
    #---------------------------------------------------------------

    score_df <- cbind(erosion_df, entropy_df)

    output <- list(
        "peaks_gr" = peaks_gr,
        "state_counts" = state_mat,
        "state_frac" = mat_list$frac,
        "state_CLR" = mat_list$CLR,
        "state_z" = mat_list$zscore,
        "scores" = score_df
    )

    return(output)
}

#' Annotate Peaks with Chromatin States
#'
#' Assigns each peak in a \code{GRanges} object to the chromatin state with 
#' which it has the largest overlap from a ChromHMM \code{GRanges} annotation.
#' Supports filtering by both absolute and fractional overlap requirements.  
#'
#' @param peaks_gr A \code{GRanges} object of peak regions to annotate.
#' @param chromHMM_states A \code{GRanges} object of ChromHMM state annotations. 
#' Must contain a metadata column (default \code{"name"}) specifying state names.
#' @param state_col Character string specifying the metadata column in 
#' \code{chromHMM_states} containing state labels. Default = \code{"name"}.
#' @param keep_unannotated Logical; if \code{TRUE}, peaks without qualifying overlaps 
#' to any ChromHMM state are retained with \code{NA} in their \code{annotation} column.
#' If \code{FALSE} (default), unannotated peaks are dropped.
#' @param min_overlap Numeric; minimum overlap width (in bp) required for a peak 
#' to be assigned to a ChromHMM state. Default = \code{1} (any overlap).
#' @param min_overlap_frac Optional numeric; minimum fraction of a peak's length that 
#' must overlap a ChromHMM state for assignment. Must be between 0 and 1. 
#' If \code{NULL} (default), no fraction filter is applied.
#' @param verbose Logical; if \code{TRUE} (default), prints a summary message about 
#' peaks removed due to overlap thresholds.
#'
#' @details
#' For each peak:
#' \enumerate{
#'   \item The function finds all overlaps between peaks and ChromHMM states.
#'   \item It computes the absolute overlap width and fractional overlap relative 
#'         to peak length.
#'   \item It assigns each peak to the state with the largest overlap, subject to 
#'         \code{min_overlap} and \code{min_overlap_frac} thresholds.
#'   \item Overlap statistics are stored in new metadata columns.
#'   \item Peaks that fail the filters are optionally retained (\code{keep_unannotated}) 
#'         or dropped (default).
#' }
#'
#' @return A \code{GRanges} object of peaks with added metadata columns:
#' \itemize{
#'   \item \code{annotation} — the ChromHMM state label.
#'   \item \code{overlap_width} — width of overlap between the peak and 
#'         its assigned ChromHMM state.
#'   \item \code{overlap_frac} — fraction of the peak length overlapping the state.
#' }
#'
#' @seealso \code{\link{ChromatinStateSigns}}, \code{\link{ErosionScore}}
#' @export
AnnotatePeaks <- function(
    peaks_gr, 
    chromHMM_states,
    state_col = "name",
    keep_unannotated = FALSE,
    min_overlap = 1,
    min_overlap_frac = NULL,
    verbose = TRUE
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
    if (!is.null(min_overlap_frac) && (!is.numeric(min_overlap_frac) || min_overlap_frac <= 0 || min_overlap_frac > 1)) {
        stop("`min_overlap_frac` must be a numeric between 0 and 1.")
    }

    # find overlaps between peaks and chromHMM states
    ov <- findOverlaps(peaks_gr, chromHMM_states)
    if (length(ov) == 0) {
        warning("No overlaps found between peaks and chromHMM states. Returning input peaks_gr unchanged.")
        peaks_gr$annotation <- NA_character_
        peaks_gr$overlap_width <- NA_integer_
        return(peaks_gr)
    }

    total_peaks <- length(peaks_gr)

    # compute overlap widths
    ov_df <- as.data.frame(ov)
    ov_df$width <- width(
        GenomicRanges::pintersect(peaks_gr[ov_df$queryHits],
                   chromHMM_states[ov_df$subjectHits])
    )

    # compute peak lengths (for fraction filtering)
    peak_lengths <- width(peaks_gr)
    ov_df$peak_length <- peak_lengths[ov_df$queryHits]
    ov_df$frac <- ov_df$width / ov_df$peak_length

    # for each peak, choose the overlap with maximum width
    ov_best <- ov_df %>%
        dplyr::group_by(.data$queryHits) %>%
        dplyr::slice_max(.data$width, n = 1, with_ties = FALSE) %>%
        dplyr::ungroup()

    # apply minimum overlap filter
    before_filter <- nrow(ov_best)
    
    # apply filters
    ov_best <- ov_best[ov_best$width >= min_overlap, ]
    if (!is.null(min_overlap_frac)) {
        ov_best <- ov_best[ov_best$frac >= min_overlap_frac, ]
    }

    after_filter <- nrow(ov_best)

    removed <- before_filter - after_filter
    perc_removed <- round((removed / total_peaks) * 100, 2)

    if (verbose) {
        msg <- paste0(
            "AnnotatePeaks: Removed ", removed, " peaks (", perc_removed, "%) ",
            "that did not meet min_overlap = ", min_overlap
        )
        if (!is.null(min_overlap_frac)) {
            msg <- paste0(msg, " and/or min_overlap_frac = ", min_overlap_frac)
        }
        message(msg)
    }

    # map peaks to state labels and overlap widths
    peak_to_state <- rep(NA_character_, length(peaks_gr))
    peak_to_width <- rep(NA_integer_, length(peaks_gr))
    peak_to_frac<- rep(NA_integer_, length(peaks_gr))

    peak_to_state[ov_best$queryHits] <- 
        mcols(chromHMM_states)[ov_best$subjectHits, state_col]
    peak_to_width[ov_best$queryHits] <- ov_best$width
    peak_to_frac[ov_best$queryHits] <- ov_best$frac

    # annotate peaks
    peaks_gr$annotation <- peak_to_state
    peaks_gr$overlap_width <- peak_to_width
    peaks_gr$overlap_frac <- peak_to_frac

    # filter out peaks without annotation if requested
    if (!keep_unannotated) {
        peaks_gr <- peaks_gr[!is.na(peaks_gr$annotation)]
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
    peak_counts <- Matrix::rowSums(peaks_mat)

    # number of cells with nonzero counts per peak
    cells_per_peak <- Matrix::rowSums(peaks_mat > 0)

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
#' @param state_names Character string containing a unique list of chromatin state names.
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
    state_names,
    active_patterns = c("TssA", "TssFlnk", "Tx", "EnhA", "EnhG", "EnhWk"),
    repressive_patterns = c("ReprPC", "Quies", "Het"),
    error_if_unclassified = FALSE
){

    # initialize the sign vector
    state_signs <- rep(0, length(state_names))
    names(state_signs) <- state_names

    # assign active states (-1)
    for (pat in active_patterns) {
        state_signs[grepl(pat, state_names, ignore.case = TRUE)] <- -1
    }

    # assign repressive states (+1)
    for (pat in repressive_patterns) {
        state_signs[grepl(pat, state_names, ignore.case = TRUE)] <- +1
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


#' Normalize a chromatin state counts matrix
#'
#' This function takes a cell-by-chromatin-state counts matrix and produces
#' three normalized versions: fractions (per cell), centered log-ratio (CLR),
#' and z-scored CLR values. Optionally, normalization can be performed relative
#' to a reference cell population defined in the metadata.
#'
#' @param state_mat A numeric matrix of chromatin state counts with cells as rows
#'   and chromatin states as columns.
#' @param meta A data frame containing cell metadata. Must have rownames matching
#'   \code{rownames(state_mat)}.
#' @param pseudocount A numeric value added to counts prior to fraction
#'   calculation to avoid division by zero. Default is 0.5.
#' @param group_by (Optional) A column name in \code{meta} specifying the variable
#'   used to define a reference cell group.
#' @param group_name (Optional) The value within \code{group_by} that defines
#'   the reference group of cells.
#' @param baseline_mu (Optional) A numeric vector of baseline means per chromatin
#'   state. If provided, used for reference-based z-scoring.
#' @param baseline_sigma (Optional) A numeric vector of baseline standard
#'   deviations per chromatin state. If provided, used for reference-based
#'   z-scoring.
#'
#' @return A list with three elements:
#'   \describe{
#'     \item{\code{frac}}{Matrix of fractions of each chromatin state per cell.}
#'     \item{\code{CLR}}{Matrix of centered log-ratio (CLR) transformed fractions.}
#'     \item{\code{zscore}}{Matrix of z-scored CLR values, either relative to the
#'     full dataset or to a specified reference group.}
#'   }
#'
#' @examples
#' # Normalize a counts matrix without a reference group
#' out <- NormalizeStateMatrix(state_mat, meta)
#'
#' # Normalize relative to a reference group
#' out <- NormalizeStateMatrix(state_mat, meta, group_by = "cluster", group_name = "Excitatory")
#'
#' @export
NormalizeStateMatrix <- function(
  state_mat,
  meta,
  pseudocount = 0.5,
  group_by = NULL,
  group_name = NULL,
  baseline_mu = NULL,
  baseline_sigma = NULL 
){

    # TODO: check that group_by and group_name are valid (present in meta)
    if(!is.null(group_by) | !is.null(group_name)){
        if(!(group_by %in% colnames(meta))){
            stop(paste0("group_by '", group_by, "' not found in meta"))
        }
        if(!(group_name %in% meta[[group_by]])){
            stop(paste0("group_name '", group_name, "' not found in meta$", group_by))
        }
    }

    # calculate fraction
    state_frac <- (state_mat + pseudocount) / rowSums(state_mat + pseudocount)

    # centered log-ratio normalization (CLR)
    log_frac <- log(state_frac)
    log_frac_centered <- log_frac - rowMeans(log_frac)

    if(!is.null(group_by) & !is.null(group_name)){
        reference_cells <- meta %>% subset(get(group_by) == group_name) %>% rownames

        # 1. Subset the log-ratio matrix to ONLY the reference cells
        reference_matrix <- log_frac_centered[reference_cells, ]

        # 2. Calculate the mean (mu) and standard deviation (sigma) for each state (column)
        baseline_mu <- colMeans(reference_matrix)
        baseline_sigma <- apply(reference_matrix, 2, sd)
    }

    #---------------------------------------------------------------------
    # Z-scoring by the full dataset, or by a reference cell population?
    #---------------------------------------------------------------------
    if(!is.null(baseline_mu) & !is.null(baseline_sigma)){
        # Baseline Z-scoring (Z-score against a reference group)
        # Center the data using the reference mean
        mat_centered <- sweep(log_frac_centered, 2, baseline_mu, "-")
        
        # Scale the data using the reference standard deviation
        z_mat <- sweep(mat_centered, 2, baseline_sigma, "/")

    } else {
        # Standard Z-scoring (Z-score against the whole input population)
       z_mat <- scale(log_frac_centered)
    }

    # setup the output 
    output <- list(
      "frac" = state_frac,
      "CLR" = log_frac_centered,
      "zscore" = z_mat
    )

    return(output)
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
#' @param mat A numeric matrix of chromatin states as columns an cells as rows 
#' and states as columns. Typically produced by \code{CalculateStateMatrix()}.
#' @param meta A data frame (such as \code{seurat_obj@meta.data}) containing 
#' cell-level metadata. Required if covariate regression is to be performed.
#' @param state_signs A named numeric vector indicating the sign (+1 or -1) 
#' for each chromatin state. Names must match the column names of \code{input}.
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
#'   \item \code{erosion_score}: the raw erosion score per cell.
#'   \item \code{erosion_score_corrected} (if covariates given): residual 
#'   erosion score with covariates regressed out.
#' }.
#'
#' @seealso \code{\link{EntropyScore}}, \code{\link{CalculateStateMatrix}}
#' @export
ErosionScore <- function(
    mat, 
    meta,
    state_signs,
    covariates = NULL
){

    # Multiply each state by its sign (+1 for repressive, -1 for active)
    signed_mat <- sweep(mat, 2, state_signs[colnames(mat)], "*")
    
    # Sum across states to get one erosion score per cell
    erosion_score <- rowSums(signed_mat)

    # initialize output df
    out_df <- data.frame("erosion" = erosion_score) 

    if(!is.null(covariates)){
        # build data frame for regression and fit model
        reg_df <- meta[, covariates, drop=FALSE]
        reg_df$erosion <- erosion_score
        form <- as.formula(paste("erosion ~", paste(covariates, collapse = " + ")))
        model <- lm(form, data = reg_df)
        out_df$erosion_corrected <- resid(model)
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
#' @param mat A numeric matrix of chromatin state counts 
#' (rows = cells, columns = states). Typically generated by \code{CalculateStateMatrix()}.
#' @param meta A data frame of cell-level metadata (e.g. \code{seurat_obj@meta.data}). 
#' Required if \code{covariates} is provided.
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
#'   \item \code{entropy_norm}: entropy normalized to 0–1 range by 
#'   dividing by \code{log2(n_states)}.
#'   \item \code{entropy_corrected} (if \code{covariates} given): 
#'   residual entropy score with covariates regressed out.
#' }
#'
#' @seealso \code{\link{ErosionScore}}, \code{\link{CalculateStateMatrix}}
#' @export
EntropyScore <- function(
    mat, 
    meta,
    covariates = NULL
){

    # calculate entropy scores per cell
    entropy_score <- apply(
        mat, 
        1, 
        function(p) {
            entropy::entropy(p, unit = "log2")  # Shannon entropy in bits
        }
    )
    
    n_states <- ncol(mat)
    entropy_score_norm <- entropy_score / log2(n_states)

    out_df <- data.frame(
        "entropy" = entropy_score,
        "entropy_norm" = entropy_score_norm
    )

  if(!is.null(covariates)){

        # build data frame for regression
        reg_df <- meta[, covariates, drop=FALSE]
        reg_df$entropy_score <- entropy_score  # regress the raw entropy score

        # build formula: entropy_score ~ covar1 + covar2 + ...
        form <- as.formula(paste("entropy_score ~", paste(covariates, collapse = " + ")))
        model <- lm(form, data = reg_df)
        out_df$entropy_corrected <- resid(model)
    }

    return(out_df)

}


