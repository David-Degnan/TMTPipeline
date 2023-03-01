#' Main labeled pipeline script
#'
#' @param msnid (MSnID object) collated MSGF output. Required.
#' @param masic An object of the masic data class. Required.
#' @param metadata A data frame with the IonChannelNames, PlexNames, and SampleNames,
#'     and any other metadata to end up in the f_data file. SampleNames will
#'     replace the IonChannelNames and PlexNames combinations.
#' @param plex_data Contains at a minimum the Dataset and PlexNames information.
#'     Both f_data and plex_data are mapped to their plex and ion channels for merging. Required.
#' @param outpath The folder to write the e_data, e_meta, and f_data files to. Required.
#' @param interference_score_threshold A numeric between 0-1 to view the maximum interference.
#'     The higher the number, the cleaner parent ion at MS1 level. Default is 0.5.
#' @param ascore_threshold If phosphoproteomics, filter by an A score threshold.
#'     Default is NULL.
#'
#' @export
labeled_pipeline <- function(msnid,
                             masic,
                             metadata,
                             plex_data,
                             outpath = "~/Downloads/",
                             interference_score_threshold = 0.5,
                             a_score_threshold = NULL) {

  # 1. FDR filter the MS-GF data
  msgf_filtered <- fdr_filter(msnid)

  # 2. Decoy filter the MS-GF data
  msgf_nodecoy <- decoy_filter(msgf_filtered)

  # 3. Filter the masic data
  masic <- as_masic_labeled(masic)
  masic_filtered <- interference_filter_labeled(masic, interference_score_threshold)

  # 4. Generate the f_data, e_data, and e_meta
  f_data <- create_f_data(masic_filtered, metadata)
  e_objects <- create_e_objects_labeled(masic_filtered, msgf_nodecoy, f_data, plex_data)

  # Write results
  write.csv(e_objects[[1]], file.path(outpath, "e_data.csv"), quote = F, row.names = F)
  write.csv(e_objects[[2]], file.path(outpath, "e_meta.csv"), quote = F, row.names = F)
  write.csv(f_data, file.path(outpath, "f_data.csv"), quote = F, row.names = F)

}
