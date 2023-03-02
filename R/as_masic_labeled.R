#' Make a MASIC data object for labeled data
#'
#' @param masic (data.frame) The data.frame to be converted into a labeled masic data object. Required.
#' @return (masic_data_labeled) A labeled masic data object. 
#'
#' @importFrom dplyr %>%
#'
#' @export
as_masic_labeled <- function(masic) {

  #################
  ## CHECK INPUT ##
  #################

  # Ensure the input is a data frame
  if (!inherits(masic, "data.frame")) {
    stop("masic must be a data.frame")
  }

  # Ensure the ion channels exist and other important columns
  IonChannels <- c("Ion_126.128", "Ion_127.125", "Ion_127.131", "Ion_128.128",
                   "Ion_128.134", "Ion_129.131", "Ion_129.138", "Ion_130.135",
                   "Ion_130.141", "Ion_131.138", "Ion_131.144", "Ion_132.142",
                   "Ion_132.148", "Ion_133.145", "Ion_133.151", "Ion_134.148")
  ImportantCols <- c("Dataset", "ScanNumber", IonChannels, "InterferenceScore")

  # Make sure the data frame contains the Ion Channels
  ColCheck <-  ImportantCols %in% colnames(masic) %>% all()
  if (!ColCheck) {
    stop(
      paste("One of the required columns for a masic_data_labeled is missing. The required list is:",
            ImportantCols %>% paste0(collapse = ", "))
    )
  }

  ############################
  ## MAKE AND RETURN OBJECT ##
  ############################

  class(masic) <- c(class(masic), "masic_data_labeled")

  return(masic)

}
