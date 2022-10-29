#' Filter MASIC data by an interference score and remove all missing measurements
#'
#' @param masic An object of the masic data class. Required.
#' @param interference_score_threshold A numeric between 0-1 to view the maximum interference.
#'     The higher the number, the cleaner parent ion at MS1 level. Default is 0.5.
#' @param remove_all_missing A logical to indicate whether a dataset and scan number
#'     where there is no intensity/abundance at any ion intensity should be removed.
#'     Default is TRUE.
#' @importFrom dplyr %>%
#'
#' @return A filtered object of the masic data class
interference_filter <- function(masic,
                                interference_score_threshold = 0.5,
                                remove_all_missing = TRUE) {

  ##################
  ## CHECK INPUTS ##
  ##################

  if (!inherits(masic, "masic_data")) {
    stop("interference_filter requires a masic_data object.")
  }

  # Interference filter should be a value between 0 and 1
  if (!is.numeric(interference_score_threshold) || interference_score_threshold < 0 | interference_score_threshold > 1) {
    stop("intereference_score_threshold should be a numeric that range from 0 to 1.")
  }

  # Remove all missing should be a true or false
  if (remove_all_missing != TRUE & remove_all_missing != FALSE) {
    stop("remove_all_missing should be true or false.")
  }

  ##################
  ## APPLY FILTER ##
  ##################

  # Apply and track filter
  filtApplied <- masic[masic$InterferenceScore >= interference_score_threshold,]
  attr(filtApplied, "TMTPipeline")$InterferenceFiltered <- TRUE

  # Remove all missing if necessary
  if (remove_all_missing) {

    # Get intensity columns
    IntensColumns <- colnames(masic)[!grepl("_OriginalIntensity|_SignalToNoise|InterferenceScore|Dataset|ScanNumber", colnames(masic))]

    # Add a sum intensity
    class(filtApplied) <- c("data.table", "data.frame")

    filtApplied <- filtApplied %>%
      dplyr::mutate(Sum = rowSums(.[IntensColumns])) %>%
      dplyr::filter(Sum > 0) %>%
      dplyr::select(-Sum)

    class(filtApplied) <- c("data.table", "data.frame", "masic_data")

  }

  return(filtApplied)

}
