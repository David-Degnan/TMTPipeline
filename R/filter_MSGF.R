#' Filter msnid object to a specific FDR
#'
#' Filtering MSGF data. In this implementation the peptide level filter optimizes both ppm and
#' PepQValue thresholds to achieve maximum number of peptide identifications within
#' the given FDR constraints. https://github.com/vladpetyuk/PlexedPiper/blob/master/R/filter_msgf_data.R
#'
#' @param msnid (MSnID object) collated MSGF output. Required.
#' @param level (character) Level at which to perform FDR filter. Options are peptide and accession (protein). Default is peptide.
#' @param fdr.max (numeric) Maximum acceptable FDR rate. Default is 0.01.
#' @param n.iter.grid (numeric) Bumber of grid-distributed evaluation points. Default 500.
#' @param n.iter.nm (numeric) Number of iterations for Nelder-Mead optimization algorithm. Default is 100.
#' @return (MSnID object) Filtered MSGF output
#'
#' @importFrom MSnID MSnIDFilter MSnIDFilter optimize_filter mass_measurement_error apply_filter
#' @export
fdr_filter <- function(msnid,
                       level = "peptide",
                       fdr.max = 0.01,
                       n.iter.grid = 500,
                       n.iter.nm = 100) {

  ####################
  ## INPUT CHECKING ##
  ####################

  # Input should be msnid object
  if (!inherits(msnid, "MSnID")) {
    stop("msnid should be a MSnID object.")
  }

  # The msnid object should contain PepQValue
  if (!"PepQValue" %in% names(msnid)) {
    if ("MS-GF:PepQValue" %in% names(msnid)) {
      message("MS-GF:PepQValue column detected. Generating a PepQValue from it.")
      msnid$PepQValue <- msnid$`MS-GF:PepQValue`
    } else {
      stop("No PepQValue column detected.")
    }
  }

  # If using accession, trigger warning
  if (level == "accession") {
    message("Protein level filtering is not necessary, as that step is carried out in pmartR.")
  }

  ###################
  ## RUN ALGORITHM ##
  ###################

  # setup
  if (level == "peptide") {
    msnid$msmsScore <- -log10(msnid$PepQValue)
    msnid$absParentMassErrorPPM <- abs(mass_measurement_error(msnid))
    filtObj <- MSnIDFilter(msnid)
    filtObj$absParentMassErrorPPM <- list(comparison="<", threshold=10.0)
    filtObj$msmsScore <- list(comparison=">", threshold=2.0)
    method <- "Nelder-Mead"
  } else if (level == "accession") {
    filtObj <- MSnIDFilter(msnid)
    filtObj$peptides_per_1000aa <- list(comparison=">", threshold=1.0)
    method <- "SANN"
  }

  # step 1
  filtObj.grid <- optimize_filter(filtObj,
                                  msnid,
                                  fdr.max=fdr.max,
                                  method="Grid",
                                  level=level,
                                  n.iter=n.iter.grid)
  # step 2
  filtObj.nm <- optimize_filter(filtObj.grid,
                                msnid,
                                fdr.max=fdr.max,
                                method=method,
                                level=level,
                                n.iter=n.iter.nm)

  # Apply the filter
  filtApplied <- apply_filter(msnid, filtObj.nm)

  # Add an attributes
  attr(filtApplied, "TMTPipeline")$FDRFiltered <- TRUE

  return(filtApplied)

}


#' Remove decoys from an msnid object
#'
#' Filtering decoys from an msnid object AFTER the FDR filter has been applied.
#'
#' @param msnid (MSnID object) collated MSGF output. Required.
#' @return (MSnID object) Filtered MSGF output
#'
#' @export
decoy_filter <- function(msnid) {

  ####################
  ## INPUT CHECKING ##
  ####################

  # Input should be msnid object
  if (!inherits(msnid, "MSnID")) {
    stop("msnid should be a MSnID object.")
  }

  # Data should be already FDR filtered
  if (is.null(attr(msnid, "TMTPipeline")$FDRFiltered)) {
    stop("Run the fdr_filter first.")
  }

  ###################
  ## RETURN RESULT ##
  ###################

  # Apply and track filter
  filtApplied <- MSnID::apply_filter(msnid, "isDecoy == FALSE")
  attr(filtApplied, "TMTPipeline")$DecoyFiltered <- TRUE

  return(filtApplied)

}


#' Filter msnid by an ascore
#' 
#' Requires an extra ascore_file that can be pulled with get_Ascore. Any NA AScores
#' are maintained. 
#' 
#' @param msnid (MSnID object) collated MSGF output. Required.
#' @param ascore_file (data.frame) Add an AScore file to filter phosphoproteomics data by uploading
#'     this data.frame. Default is NULL.  
#' @param ascore_threshold (numeric) If phosphoproteomics, filter by an AScore threshold.
#'     Default is 17, see [here](https://www.nature.com/articles/nbt1240)
#' @return (MSnID object) Filtered MSGF output
#' 
#' @importFrom dplyr %>%
#' 
#' @export
ascore_filter <- function(msnid, 
                          ascore_file,
                          ascore_threshold = 17) {
  
  ####################
  ## INPUT CHECKING ##
  ####################
  
  # Input should be msnid object
  if (!inherits(msnid, "MSnID")) {
    stop("msnid should be a MSnID object.")
  }
  
  # ascore needs to be a data.frame
  if (!inherits(ascore_file, "data.frame")) {
    stop("ascore_file should be a data.frame")
  }
  
  # ascore needs the spectrumFile, Scan, OriginalSequence, and AScore columns 
  if (!all(c("spectrumFile", "Scan", "OriginalSequence", "AScore") %in% colnames(ascore_file))) {
    stop("ascore_file requires the spectrumFile, Scan, peptide, and AScore columns.")
  }
  
  # ascore_threshold should be a single numeric value
  if (length(ascore_threshold) != 1 || !is.numeric(ascore_threshold)) {
    stop("ascore_threshold should be a single numeric value.")
  }
  
  ################
  ## RUN FILTER ##
  ################
  
  # Pull relevant columns 
  ascore <- ascore_file %>%
    dplyr::select(spectrumFile, Scan, OriginalSequence, AScore) %>%
    dplyr::rename(peptide = OriginalSequence)
  
  # Combine and filter 
  filtered <- msnid@psms %>%
    dplyr::left_join(ascore, by = c("spectrumFile", "Scan", "peptide")) %>%
    dplyr::filter(AScore >= ascore_threshold)
  
  ###################
  ## RETURN RESULT ##
  ###################
  
  # Update msnid 
  msnid@psms <- filtered
  
  # Add TMTPipeline attributes
  attr(msnid, "TMTPipeline")$AScore_Filtered <- TRUE

  return(msnid)
  
}