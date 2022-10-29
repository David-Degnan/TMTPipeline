#' Create an f_data (sample information object) from masic data and a metadata file
#'
#' @param masic A masic_data object that has been filtered by intereference_filter. Required.
#' @param metadata A data frame with the PlexNames, IonChannelNames, SampleName,
#'    and any other metadata to end up in the f_data file. Know that PlexNames and
#'    IonChannelNames are used to find the data in the masic table, and Sample Name
#'    is what the column names will be in the final table. Required.
#'
#' @return A pmartR f_data (sample information) object
#' @export
create_f_data <- function(masic,
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

#' Create an e_data (biomolecule expression or crosstab) and e_meta (biomolecule data) from masic data, msnid, and f_data
#'
#' @param masic A masic_data object that has been filtered by intereference_filter. Required.
#' @param msnid An msnid object that has been fdr and decoy filtered. Required.
#' @param f_data An f_data object from create_f_data. Required.
#'
#' @return A list with a pmartR e_data (biomolecule expression) object and a pmartR e_meta (biomolecule data) object
#' @export
create_e_objects <- function(masic,
                            msnid,
                            f_data) {

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

  # Input should be an msnid object
  if (!inherits(msnid, "MSnID")) {
    stop("msnid should be a MSnID object.")
  }

  # msnid should already be FDR filtered
  if (is.null(attr(msnid, "TMTPipeline")$FDRFiltered)) {
    stop("Run the fdr_filter on the msnid data first.")
  }

  # msnid should already be decoy filtered
  if (is.null(attr(msnid, "TMTPipeline")$DecoyFiltered)) {
    stop("Run the decoy_filter on the msnid data first.")
  }

  # f_data should have been created with the appropriate function
  if (is.null(attr(f_data, "valid_f_data"))) {
    stop("f_data should be created with create_f_data.")
  }

  #############################
  ## BUILD THE E_DATA OBJECT ##
  #############################

  # Pull the peptide identifications. Select only the required columns for pmartR.
  # Get the protein and contaminant names. If the protein matches the contaminant
  # name only, it is a contaminant. If it matches both the contaminant and protein
  # names, it is potentially a contaminant. Otherwise, it is not a contaminant.
  psm_data <- MSnID::psms(msnid) %>%
    dplyr::select(Dataset, Scan, peptide, Protein) %>%
    dplyr::rename(ScanNumber = Scan, Peptide = peptide) %>%
    dplyr::mutate(
      ContaminantName = ifelse(grepl("Contaminant", Protein), gsub("Contaminant_", "", Protein), NA),
      ProteinName = purrr::map(Protein, function(x) {
        ifelse(grepl("Contaminant_", x), NA, strsplit(x, "|", fixed = T) %>% unlist() %>% tail(1))
      }) %>% unlist()
    ) %>%
    dplyr::group_by(Dataset, ScanNumber, Peptide) %>%
    tidyr::nest() %>%
    dplyr::mutate(
      Protein = purrr::map(data, function(x) {x[,c("ContaminantName", "ProteinName")] %>% unlist() %>% unique() %>% .[!is.na(.)]}),
    ) %>%
    tidyr::unnest(Protein) %>%
    dplyr::mutate(
      Contaminant = purrr::map2(data, Protein, function(x, y) {
        ifelse(y %in% x$ContaminantName, ifelse(y %in% x$ProteinName, "Potential", "Yes"), "No")
      }) %>% unlist()
    ) %>%
    dplyr::select(-data) %>%
    dplyr::rename(PlexNames = Dataset)

  # Generate the emeta object
  e_meta <- psm_data[,c("Peptide", "Protein", "Contaminant")] %>% unique()

  class(masic) <- c("data.table", "data.frame")

  # Create the e_data object by combining the masic data, peptide spectrum matches,
  # and the f_data. It is definitely possible for their to be multiple peptide
  # identifications in a sample. For now, return the maximum value.
  e_data <- masic %>%
    dplyr::select(c(Dataset, ScanNumber, f_data$IonChannelNames)) %>%
    tidyr::pivot_longer(f_data$IonChannelNames) %>%
    dplyr::rename(PlexNames = Dataset, IonChannelNames = name) %>%
    dplyr::inner_join(psm_data[,c("PlexNames", "ScanNumber", "Peptide", "Protein")], by = c("PlexNames", "ScanNumber")) %>%
    dplyr::inner_join(f_data[,c("PlexNames", "IonChannelNames", "SampleNames")], by = c("PlexNames", "IonChannelNames")) %>%
    dplyr::group_by(SampleNames, Peptide) %>%
    dplyr::summarise(value = sum(value)) %>%
    tidyr::pivot_wider(values_from = value, names_from = SampleNames)

  e_data[is.na(e_data)] <- 0

  # Filter e_meta to identified peptides
  e_meta <- e_meta[e_meta$Peptide %in% e_data$Peptide,] %>% dplyr::ungroup()

  # Return objects
  return(list(data.frame(e_data),
              data.frame(e_meta)))

}

