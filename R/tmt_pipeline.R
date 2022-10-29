#' Main pipeline script to run the TMT Pipeline
#'
#' @param msnid (MSnID object) collated MSGF output. Required.
#' @param masic An object of the masic data class. Required.
#' @param interference_score_threshold A numeric between 0-1 to view the maximum interference.
#'     The higher the number, the cleaner parent ion at MS1 level. Default is 0.5.
#' @param remove_all_missing A logical to indicate whether a dataset and scan number where
#'     there is no intensity/abundance at any ion intensity should be removed. Default is TRUE.
#' @param metadata A data frame with the PlexNames, IonChannelNames, SampleName,
#'    and any other metadata to end up in the f_data file. Know that PlexNames and
#'    IonChannelNames are used to find the data in the masic table, and Sample Name
#'    is what the column names will be in the final table. Required.
#' @param outpath The folder to write the e_data, e_meta, and f_data files to. Required,
#'
#' @export
tmt_pipeline <- function(msnid,
                         masic,
                         interference_score_threshold = 0.5,
                         remove_all_missing = TRUE,
                         metadata,
                         outpath = "~/Downloads/") {

  # 1. FDR filter the MS-GF data
  msgf_filtered <- fdr_filter(msnid)

  # 2. Decoy filter the MS-GF data
  msgf_nodecoy <- decoy_filter(msgf_filtered)

  # 3. Filter the masic data
  class(masic) <- c(class(masic), "masic_data")
  masic_filtered <- interference_filter(masic, interference_score_threshold, remove_all_missing)

  # 4. Generate the f_data, e_data, and e_emta
  f_data <- create_f_data(masic_filtered, metadata)
  e_objects <- create_e_objects(masic, msnid, f_data)

  # Write results
  write.csv(e_objects[[1]], file.path(outpath, "e_data.csv"), quote = F, row.names = F)
  write.csv(e_objects[[2]], file.path(outpath, "e_meta.csv"), quote = F, row.names = F)
  write.csv(f_Data, file.path(outpath, "f_data.csv"), quote = F, row.names = F)

}
