#' Filter MASIC data by an interference score and remove all missing measurements
#'
#' @param masic An object of the masic data class. Required.
#' @param interference_score_threshold A numeric between 0-1 to view the maximum interference.
#'     The higher the number, the cleaner parent ion at MS1 level. Default is 0.5.
#' @importFrom dplyr %>%
#'
#' @return A filtered object of the masic data class
#' @export
interference_filter <- function(masic,
                                interference_score_threshold = 0.5) {

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

  ##################
  ## APPLY FILTER ##
  ##################

  # Apply and track filter
  filtApplied <- masic[masic$InterferenceScore >= interference_score_threshold,]
  attr(filtApplied, "TMTPipeline")$InterferenceFiltered <- TRUE

  return(filtApplied)

}
