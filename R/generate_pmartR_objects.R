#' Generate an f_data (sample information object) from masic data and a metadata file
#'
#' @param masic A masic_data object that has been filtered by intereference_filter. Required.
#' @param metadata A data frame with the PlexNames, IonChannelNames, SampleName,
#'    and any other metadata to end up in the f_data file. Know that PlexNames and
#'    IonChannelNames are used to find the data in the masic table, and Sample Name
#'    is what the column names will be in the final table. Required.
#'
#' @return A pmartR f_data (sample information) object
#' @export
generate_f_data <- function(masic,
                            metadata) {

  ##################
  ## CHECK INPUTS ##
  ##################

  # masic data should be a masic data object
  if (!inherits(masic, "masic_data")) {
    stop("A masic_data object is required.")
  }

  # masic data is suggested to be interference filtered
  if (is.null(attr(masic, "TMTPipeline")$InterferenceFiltered)) {
    message("Applying an intereference filter is recommended before building pmartR objects")
  }

  # Assert that metadata has the required columns
  if (!all(c("PlexNames", "IonChannelNames", "SampleNames") %in% colnames(metadata))) {
    stop("PlexNames, IonChannelNames, and SampleNames must be in the metadata object.")
  }

  # Assert that all datasets have the proper plex names
  if (!all(unique(masic$Dataset) %in% unique(metadata$PlexNames))) {
    stop(paste0("PlexNames is missing ", unique(masic$Dataset)[unique(masic$Dataset) %in% unique(metadata$PlexNames)], sep = ", "))
  }

  # Get intensity column names
  IntensColumns <- colnames(masic)[!grepl("_OriginalIntensity|_SignalToNoise|InterferenceScore|Dataset|ScanNumber", colnames(masic))]

  # Assert that all ion names are included
  if (!all(IntensColumns %in% unique(metadata$IonChannelNames))) {
    stop(paste0("IonChannelNames is missing ", IntensColumns[IntensColumns %in% unique(metadata$IonChannelNames)], sep = ", "))
  }

  # Let the user know if blanks were detect
  if ("" %in% metadata$SampleNames) {
    message("Blanks detected in metadata. Those entries will be removed.")
    metadata <- metadata[metadata$SampleNames != "",]
  }

  # Assert that sample names are unique
  if (length(unique(metadata$SampleNames)) != length(metadata$SampleNames)) {
    stop("Duplicate names detected in SampleNames")
  }

  #############################
  ## BUILD THE F_DATA OBJECT ##
  #############################

  # The point is to check and clean the f_data file. This is all for now.
  f_data <- metadata

  # Add attributes
  attr(f_data, "valid_f_data") <- TRUE

  return(f_data)

}
