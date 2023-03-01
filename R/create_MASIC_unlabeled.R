#' Create MASIC data when the proteomics dataset is unlabeled
#'
#' @param folder_path Path to the SICstats txt files. Required.
#' @param interference_score_threshold A numeric between 0-1 to view the maximum interference.
#'      The higher the number, the cleaner parent ion at MS1 level. Default is 0.5.
#'
#' @importFrom dplyr %>%
#'
#' @return An inteference score filtered masic_data_unlabeled object
#'
#' @export
create_MASIC_unlabeled <- function(folder_path,
                                   interference_score_threshold = 0.5) {

  ##################
  ## CHECK INPUTS ##
  ##################

  # Interference filter should be a value between 0 and 1
  if (!is.numeric(interference_score_threshold) || interference_score_threshold < 0 | interference_score_threshold > 1) {
    stop("intereference_score_threshold should be a numeric that range from 0 to 1.")
  }

  # Pull all SIC files
  SIC_Files <- list.files(folder_path, full.names = TRUE)

  # Ensure all SIC files have SICstats.txt in them
  SICStatsCheck <- lapply(SIC_Files, function(x) {grepl("SICstats.txt", x)}) %>% unlist() %>% all()
  if (!SICStatsCheck) {
    stop("All files in folder_path must be SICstats.txt files.")
  }

  # Ensure all names are unique
  if (length(unique(SIC_Files)) != length(SIC_Files)) {
    stop("All SICstats.txt files in folder_path must have unique names.")
  }

  ####################################
  ## PULL DATA FROM SIC STATS FILES ##
  ####################################

  # Pull and knit files together into one data frame
  masic <- do.call(rbind,
   lapply(SIC_Files, function(theFile) {
     message(theFile)
     data.table::fread(theFile) %>%
       dplyr::select(FragScanNumber, PeakArea, InterferenceScore) %>%
       dplyr::rename(ScanNumber = FragScanNumber, Intensity = PeakArea) %>%
       dplyr::mutate(Dataset = theFile %>% strsplit("/", fixed = T) %>% unlist() %>% tail(1) %>% gsub("_SICstats.txt", "", x = ., fixed = T)) %>%
       dplyr::relocate(Dataset) %>%
       dplyr::filter(InterferenceScore >= 0.5)
   })
  )

  # Add the correct class
  class(masic) <- c(class(masic), "masic_data_unlabeled")

  # Return the result
  return(masic)

}




