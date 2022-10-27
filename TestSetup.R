library(devtools)

load_all(".")

msnid <- readRDS("~/Desktop/OldProjects/TMT/4234_Tafesse/4234_Data_Package.RDS")$msnid

# Apply the filtering
msnid
msgf_filtered <- filter_fdr(msnid)
msgf_filtered

# Apply filtering to built object
msnid_built
built_filtered <- filter_msgf_data(msnid_built)
built_filtered

