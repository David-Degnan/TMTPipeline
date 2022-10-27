library(devtools)

load_all(".")

msnid <- readRDS("~/Desktop/OldProjects/TMT/4234_Tafesse/4234_Data_Package.RDS")$msnid

# Apply the fdr filtering
msnid
msgf_filtered <- fdr_filter(msnid)
msgf_filtered

# Apply the fdr filtering to built object
msnid_built
built_filtered <- fdr_filter(msnid_built)
built_filtered

# Apply the decoy filtering
msgf_nodecoy <- decoy_filter(msgf_filtered)
msgf_nodecoy
