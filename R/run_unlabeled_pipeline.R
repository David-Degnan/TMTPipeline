#' Main unlabeled pipeline script
#'
#' @param msnid (MSnID object) collated MSGF output. Required.
#' @param folder_path Path to the SICstats txt files to generate the masic dataset. Required.
#' @param metadata A data frame with the FileNames and the SampleNames. Required.
#' @param outpath The folder to write the e_data, e_meta, and f_data files to. Required.
#' @param interference_score_threshold A numeric between 0-1 to view the maximum interference.
#'     The higher the number, the cleaner parent ion at MS1 level. Default is 0.5.
#'
#' @export
#' @export
unlabeled_pipeline <- function(msnid,
                               folder_path,
                               metadata,
                               outpath = "~/Downloads/",
                               interference_score_threshold = 0.5) {

  # 1. FDR filter the MS-GF data
  msgf_filtered <- fdr_filter(msnid)

  # 2. Decoy filter the MS-GF data
  msgf_nodecoy <- decoy_filter(msgf_filtered)

  # 3. Generate the masic data from the SICstats files
  masic <- create_MASIC_unlabeled(folder_path = folder_path, interference_score_threshold = interference_score_threshold)

  # 4. Generate the f_data, e_data, and e_meta
  f_data <- create_f_data(masic, metadata)
  e_objects <- create_e_objects_unlabeled(masic_filtered, msgf_nodecoy, f_data, plex_data)





}
