#' Fetch ChromHMM Chromatin States from Roadmap Epigenomics
#'
#' @param model Which ChromHMM model to use: 15, 18, or 25 states.
#' @param tissue Either a Roadmap Epigenome ID (e.g., "E066"), a tissue mnemonic (e.g., "LNG.IMR90"),
#' or a human-readable tissue name (e.g., "IMR90 fetal lung fibroblasts").
#' @param genome Genome assembly: "hg19" or "hg38".
#' @param mnemonics Logical, if TRUE returns the mnemonic states (default) instead of raw state numbers.
#' @param cache_dir Directory to cache downloaded files.
#' @param fuzzy_matching Logical, if TRUE allows partial/fuzzy matching with an interactive prompt.
#'
#' @return A `GRanges` object with ChromHMM states.
#' @export
FetchChromatinStates <- function(model = 18,
                                 tissue,
                                 genome = c("hg19", "hg38"),
                                 mnemonics = TRUE,
                                 cache_dir = "~/.chromatic_cache",
                                 fuzzy_matching = TRUE) {
  genome <- match.arg(genome, c("hg19", "hg38"))

  # --- 1. Match tissue code ---
  tissue_input <- tissue

  if (!fuzzy_matching) {
    ## Exact match only
    idx <- which(
      toupper(tissue_input) == toupper(tissue_map$EID) |
      tolower(tissue_input) == tolower(tissue_map$Mnemonic) |
      tolower(tissue_input) == tolower(tissue_map$Name)
    )
    if (length(idx) == 0) {
      stop("No exact match found for tissue: ", tissue_input)
    } else if (length(idx) > 1) {
      stop("Multiple exact matches found. Please be more specific.")
    }
  } else {
    ## Fuzzy / partial match with interactive selection
    eid_match      <- grepl(toupper(tissue_input), toupper(tissue_map$EID))
    mnemonic_match <- grepl(tolower(tissue_input), tolower(tissue_map$Mnemonic))
    name_match     <- grepl(tolower(tissue_input), tolower(tissue_map$Name))

    combined_match <- eid_match | mnemonic_match | name_match
    matches <- which(combined_match)

    if (length(matches) == 0) {
      stop("No match found for tissue: ", tissue_input)
    } else if (length(matches) > 1) {
      cur_map <- tissue_map[matches,]
      rownames(cur_map) <- 1:nrow(cur_map)
      message("Multiple tissues matched your input:")
      print(cur_map[, c("EID", "Mnemonic", "Name")])
      sel <- readline(prompt = "Enter the row number of the tissue you want: ")
      sel <- as.integer(sel)
      if (is.na(sel) || !(sel %in% seq_along(matches))) {
        stop("Invalid selection. Please run again.")
      }
      idx <- matches[sel]
    } else {
      idx <- matches
    }
  }

  tissue_code     <- tissue_map$EID[idx]
  tissue_name     <- tissue_map$Name[idx]
  tissue_mnemonic <- tissue_map$Mnemonic[idx]

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
  genome_suffix <- ifelse(genome == "hg38", "_hg38lift", "")

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

  # --- 8. Final message ---
  message("Loaded ", length(gr), " regions for ",
          tissue_code, " (", tissue_mnemonic, ", ", tissue_name, ") — ",
          model, "-state model, genome ", genome)

  return(gr)
}



# #' Fetch ChromHMM Chromatin States from Roadmap Epigenomics
# #'
# #' @param model Which ChromHMM model to use: 15, 18, or 25 states.
# #' @param tissue Either a Roadmap Epigenome ID (e.g., "E066"), a tissue mnemonic (e.g., "LNG.IMR90"),
# #' or a human-readable tissue name (e.g., "IMR90 fetal lung fibroblasts").
# #' @param genome Genome assembly: "hg19" or "hg38".
# #' @param mnemonics Logical, if TRUE returns the mnemonic states (default) instead of raw state numbers.
# #' @param cache_dir Directory to cache downloaded files.
# #'
# #' @return A `GRanges` object with ChromHMM states.
# #' @export
# FetchChromatinStates <- function(
#   model = 18,
#   tissue,
#   genome = c("hg19", "hg38"),
#   mnemonics = TRUE,
#   cache_dir = "~/.chromatic_cache"
# ) {
#   genome <- match.arg(genome, c("hg19", "hg38"))

#   # --- 1. Match tissue code ---
#   tissue_input <- tissue

#   # Columns to check:
#   eid_match     <- grepl(toupper(tissue_input), toupper(tissue_map$EID))
#   mnemonic_match<- grepl(tolower(tissue_input), tolower(tissue_map$Mnemonic))
#   name_match    <- grepl(tolower(tissue_input), tolower(tissue_map$Name))

#   combined_match <- eid_match | mnemonic_match | name_match
#   matches <- which(combined_match)

#   if (length(matches) == 0) {
#     stop("No match found for tissue: ", tissue_input)
#   } else if (length(matches) > 1) {
#     cur_map <- tissue_map[matches,]
#     rownames(cur_map) <- 1:nrow(cur_map)
#     message("Multiple tissues matched your input:")
#     print(cur_map[, c("EID", "Mnemonic", "Name")])
#     sel <- readline(prompt = "Enter the row number of the tissue you want: ")
#     sel <- as.integer(sel)
#     if (is.na(sel) || !(sel %in% seq_along(matches))) {
#       stop("Invalid selection. Please run again.")
#     }
#     idx <- matches[sel]
#   } else {
#     idx <- matches
#   }

#   tissue_code <- tissue_map$EID[idx]
#   tissue_name <- tissue_map$Name[idx]
#   tissue_mnemonic <- tissue_map$Mnemonic[idx]

#   # --- 2. Model URL parts ---
#   if (model == 15) {
#     model_dir <- "coreMarks/jointModel/final"
#     model_str <- "15_coreMarks"
#   } else if (model == 18) {
#     model_dir <- "core_K27ac/jointModel/final"
#     model_str <- "18_core_K27ac"
#   } else if (model == 25) {
#     model_dir <- "imputed12marks/jointModel/final"
#     model_str <- "25_imputed12marks"
#   } else {
#     stop("Model must be 15, 18, or 25")
#   }

#   # --- 3. Genome suffix ---
#   genome_suffix <- ifelse(genome == "hg38", "_hg38lift", "")

#   # --- 4. Mnemonics vs raw ---
#   mnemonic_str <- if (mnemonics) "_mnemonics" else ""

#   # --- 5. Construct URL ---
#   base_url <- "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels"
#   file_name <- paste0(
#     tissue_code, "_", model_str, genome_suffix, mnemonic_str, ".bed.gz"
#   )
#   url <- file.path(base_url, model_dir, file_name)

#   # --- 6. Cache file ---
#   cache_dir <- path.expand(cache_dir)
#   if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
#   cache_path <- file.path(cache_dir, file_name)

#   if (!file.exists(cache_path)) {
#     message("Downloading ChromHMM segmentation from: ", url)
#     utils::download.file(url, cache_path, mode = "wb")
#   } else {
#     message("Using cached file: ", cache_path)
#   }

#   # --- 7. Import BED ---
#   if (!requireNamespace("rtracklayer", quietly = TRUE)) {
#     stop("Package 'rtracklayer' is required. Please install it.")
#   }
#   gr <- rtracklayer::import(cache_path, format = "BED")

#   # --- 8. Final message ---
#   message("Loaded ", length(gr), " regions for ",
#           tissue_code, " (", tissue_mnemonic, ", ", tissue_name, ") — ",
#           model, "-state model, genome ", genome)

#   return(gr)
# }
