#' Download Roadmap ChromHMM Segmentation as GRanges
#'
#' @param model ChromHMM model. One of 15, 18, or 25.
#' @param tissue Tissue code (e.g. "E075") or human-readable name (e.g. "Liver").
#' @param genome Genome build. One of "hg19" or "hg38".
#' @param mnemonics Logical. If TRUE (default), download the mnemonics version (state names).
#'                  If FALSE, download the raw numeric states version.
#' @param cache_dir Directory to cache downloads (default "~/.chromatic_cache").
#'
#' @return A GRanges object with ChromHMM state annotations.
#' @export
#'
#' @examples
#' chromHMM <- FetchChromatinStates(model = 18, tissue = "E075", genome = "hg38")
FetchChromatinStates <- function(model = 18,
                             tissue,
                             genome = c("hg19", "hg38"),
                             mnemonics = TRUE,
                             cache_dir = "~/.chromatic_cache") {

  genome <- match.arg(genome, c("hg19", "hg38"))

  # --- 1. Look up tissue code ---
  tissue_code <- tissue
  tissue_code <- toupper(tissue_code)
  
  # Example lookup (extend this table for full Roadmap):
  # TODO, fix this
  # tissue_map <- data.frame(
  #   code = c("E066", "E075"),
  #   name = c("Liver", "Brain"),
  #   stringsAsFactors = FALSE
  # )

  # is this needed?
  tissue_map <- chromatic::tissue_map

  if (!(toupper(tissue) %in% tissue_map$EID)) {
    # try match by mnemonic
    match_idx <- match(tolower(tissue), tolower(tissue_map$Mnemonic))
    if (!is.na(match_idx)) {
      tissue_code <- tissue_map$EID[match_idx]
    } else {
      stop("Unknown tissue name/code: ", tissue)
    }
  } else {
    tissue_code <- toupper(tissue)
  }

  # --- 2. Model URL parts ---
  if (model == 15) {
    model_dir <- "coreMarks/jointModel/final"
    model_str <- "15_coreMarks"
  } else if (model == 18) {
    model_dir <- "core_K27ac/jointModel/final"
    model_str <- "18_core_K27ac"
  } else if (model == 25) {
    model_dir <- "imputed12marks/jointModel/final"
    model_str <- "25_imputed12marks"
  } else {
    stop("Model must be 15, 18, or 25")
  }

  # --- 3. Genome suffix ---
  if (genome == "hg19") {
    genome_suffix <- ""
  } else if (genome == "hg38") {
    genome_suffix <- "_hg38lift"
  }

  # --- 4. Mnemonics vs raw ---
  mnemonic_str <- if (mnemonics) "_mnemonics" else ""

  # --- 5. Construct URL ---
  base_url <- "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels"
  file_name <- paste0(
    tissue_code, "_", model_str, genome_suffix, mnemonic_str, ".bed.gz"
  )
  url <- file.path(base_url, model_dir, file_name)

  # --- 6. Cache file ---
  cache_dir <- path.expand(cache_dir)
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
  cache_path <- file.path(cache_dir, file_name)

  if (!file.exists(cache_path)) {
    message("Downloading ChromHMM segmentation from: ", url)
    utils::download.file(url, cache_path, mode = "wb")
  } else {
    message("Using cached file: ", cache_path)
  }

  # --- 7. Import BED ---
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    stop("Package 'rtracklayer' is required. Please install it.")
  }
  gr <- rtracklayer::import(cache_path, format = "BED")

  message("Loaded ", length(gr), " regions for ", tissue_code, " (", model, "-state, ", genome, ")")

  return(gr)
}
