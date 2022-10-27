library(devtools)

load_all(".")

msnid <- readRDS("~/Desktop/OldProjects/TMT/4234_Tafesse/4234_Data_Package.RDS")$msnid

####################
## MSGF Filtering ##
####################

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

# Apply the decoy filtering to built object
built_nodecoy <- decoy_filter(built_filtered)
built_nodecoy

#####################
## MASIC Filtering ##
#####################

masic <- readRDS("~/Desktop/OldProjects/TMT/4234_Tafesse/4234_Data_Package.RDS")$masic
class(masic) <- c(class(masic), "masic_data")

# Filter masic data
masic_filtered <- interference_filter(masic, 0.9, TRUE)

#############################
## Generate pmartR Objects ##
#############################


