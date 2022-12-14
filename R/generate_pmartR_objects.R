#' Create an f_data (sample information object) from masic data and a metadata file
#'
#' @param masic A masic_data object that has been filtered by intereference_filter. Required.
#' @param metadata A data frame with the IonChannelNames, PlexNames, and SampleNames,
#'    and any other metadata to end up in the f_data file. SampleNames
#'    will replace the IonChannelNames and PlexNames combinations.
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
  if (!all(c("IonChannelNames", "SampleNames", "PlexNames") %in% colnames(metadata))) {
    stop("IonChannelNames and SampleNames must be in the metadata object.")
  }

  # Get intensity column names
  IntensColumns <- colnames(masic)[!grepl("_OriginalIntensity|_SignalToNoise|InterferenceScore|Dataset|ScanNumber", colnames(masic))]

  # Assert that all ion names are included
  if (!all(IntensColumns %in% unique(metadata$IonChannelNames))) {
    stop(paste0("IonChannelNames is missing ", IntensColumns[IntensColumns %in% unique(metadata$IonChannelNames)], sep = ", "))
  }

  # Let the user know if blanks were detect
  if ("" %in% metadata$SampleNames | any(is.na(metadata$SampleNames))) {
    message("Blanks detected in metadata. Those entries will be removed.")
    metadata <- metadata[metadata$SampleNames != "" & !is.na(metadata$SampleNames),]
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
#' @param f_data An f_data object from create_f_data. Contains at a minimum the
#'     SampleNames, PlexNames, and IonChannelNames. Required.
#' @param plex_data Contains at a minimum the Dataset and PlexNames information. Both
#'     f_data and plex_data are mapped to their plex and ion channels for merging. Required.
#'
#' @return A list with a pmartR e_data (biomolecule expression) object and a pmartR e_meta (biomolecule data) object
#' @export
create_e_objects <- function(masic,
                             msnid,
                             f_data,
                             plex_data) {

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

  # plex_data can only contains Dataset and PlexNames
  if (!all(colnames(plex_data) %in% c("Dataset", "PlexNames"))) {
    stop("plex_data should contain Dataset and PlexNames")
  }

  # Assert that all datasets have the proper plex names
  if (!all(unique(masic$Dataset) %in% unique(plex_data$Dataset))) {
    stop(paste0("PlexNames in plex_data is missing ", unique(masic$Dataset)[unique(masic$Dataset) %in% unique(plex_data$PlexNames)], sep = ", "))
  }

  # Assert that plex names maps correctly
  if (!all(f_data$PlexNames %in% plex_data$PlexNames)) {
    stop("All PlexNames in f_data and plex_data should match.")
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
    dplyr::inner_join(plex_data, by = "Dataset") %>%
    dplyr::ungroup()

  # Generate the emeta object
  e_meta <- psm_data[,c("Peptide", "Protein", "Contaminant")] %>% unique()

  class(masic) <- c("data.table", "data.frame")

  # Create the e_data object by combining the masic data, peptide spectrum matches,
  # and the f_data. It is definitely possible for their to be multiple peptide
  # identifications in a sample. For now, return the maximum value.
  e_data <- masic %>%
    dplyr::select(c(Dataset, ScanNumber, f_data$IonChannelNames)) %>%
    dplyr::inner_join(plex_data, by = "Dataset") %>%
    tidyr::pivot_longer(f_data$IonChannelNames) %>%
    dplyr::rename(IonChannelNames = name) %>%
    dplyr::inner_join(f_data[,c("PlexNames", "IonChannelNames", "SampleNames")], by = c("PlexNames", "IonChannelNames")) %>%
    dplyr::inner_join(psm_data[,c("Peptide", "ScanNumber", "PlexNames", "Dataset")], by = c("ScanNumber", "PlexNames", "Dataset")) %>%
    dplyr::group_by(SampleNames, Peptide) %>%
    dplyr::summarise(value = sum(value, na.rm = T)) %>%
    tidyr::pivot_wider(values_from = value, names_from = SampleNames)

  e_data[is.na(e_data)] <- 0

  # Filter e_meta to identified peptides
  e_meta <- e_meta[e_meta$Peptide %in% e_data$Peptide,] %>% dplyr::ungroup()

  # Return objects
  return(list(data.frame(e_data),
              data.frame(e_meta)))

}

