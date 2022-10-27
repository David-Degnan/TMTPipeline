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

  return(apply_filter(msnid, filtObj.nm))

}


#' Filtering MSGF Data: https://github.com/vladpetyuk/PlexedPiper/blob/master/R/filter_msgf_data.R
#'
#' Filtering MSGF data. In this implementation the peptide level filter optimizes both ppm and
#' PepQValue thresholds to achieve maximum number of peptide identifications within
#' the given FDR constraints.
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
remove_decoys <- function(msnid) {

  return(NULL)


}
