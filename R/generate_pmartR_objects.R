#' Generate an f_data (sample information object) from masic data and a metadata file
#'
#' @param masic A masic_data object that has been filtered by intereference_filter. Required.
#' @param metadata A data frame with the Plex Names, Ion Channel Names, Sample Name,
#'    and any other metadata to end up in the f_data file. Know that Plex Names and
#'    Ion Channel names are used to find the data in the masic table, and Sample Name
#'    is what the column names will be in the final table. Required.
#'
#' @return A pmartR f_data (sample information) object
