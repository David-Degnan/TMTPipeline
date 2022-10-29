library(dplyr)
library(devtools)
load_all(".")

#######################
## Data Package 4234 ##
#######################

# Load files
msnid <- readRDS("~/Desktop/OldProjects/TMT/4234_Tafesse/4234_Data_Package.RDS")$msnid
masic <- readRDS("~/Desktop/OldProjects/TMT/4234_Tafesse/4234_Data_Package.RDS")$masic
class(masic) <- c(class(masic), "masic_data")
metadata = data.frame(
  "PlexNames" = rep(unique(masic$Dataset) %>% sort(), each = 16),
  "IonChannelNames" = rep(colnames(masic)[3:18], 9),
  "SampleNames" = c(
    "78679_EBC_B1_122_3", "78679_EBC_B1_081_6", "78679_EBC_B1_055_3", "78679_EBC_B1_055_5",
    "78679_EBC_B1_081_4","78679_EBC_B1_081_3", "78679_EBC_B1_122_5","78679_EBC_B1_055_1",
    "78679_EBC_B1_081_2","78679_EBC_B1_122_4", "78679_EBC_B1_081_5","78679_EBC_B1_055_4",
    "78679_EBC_B1_081_1","78679_EBC_B1_055_2", "78679_EBC_B1_055_6","78679_EBC_B1_Ref_1",
    "78679_EBC_B1_040_3","78679_EBC_B1_076_4", "78679_EBC_B1_040_5","78679_EBC_B1_040_4",
    "78679_EBC_B1_040_2","78679_EBC_B1_046_3", "78679_EBC_B1_046_5","78679_EBC_B1_046_1",
    "78679_EBC_B1_046_4","78679_EBC_B1_076_5", "78679_EBC_B1_040_6","78679_EBC_B1_046_6",
    "78679_EBC_B1_046_2","78679_EBC_B1_040_1", "78679_EBC_B1_076_6","78679_EBC_B1_Ref_2",
    "78679_EBC_B1_042_2","78679_EBC_B1_042_1", "78679_EBC_B1_066_3","78679_EBC_B1_042_4",
    "78679_EBC_B1_069_5","78679_EBC_B1_066_4", "78679_EBC_B1_069_4","78679_EBC_B1_066_2",
    "78679_EBC_B1_069_6","78679_EBC_B1_066_6", "78679_EBC_B1_066_1","78679_EBC_B1_042_6",
    "78679_EBC_B1_042_3","78679_EBC_B1_042_5", "78679_EBC_B1_066_5","78679_EBC_B1_Ref_3",
    "78679_EBC_B1_064_2","78679_EBC_B1_003_5", "78679_EBC_B1_014_6","78679_EBC_B1_014_5",
    "78679_EBC_B1_014_4","78679_EBC_B1_003_2", "78679_EBC_B1_064_1","78679_EBC_B1_014_1",
    "78679_EBC_B1_014_2","78679_EBC_B1_003_4", "78679_EBC_B1_003_6","78679_EBC_B1_003_3",
    "78679_EBC_B1_064_3","78679_EBC_B1_014_3", "78679_EBC_B1_003_1","78679_EBC_B1_Ref_4",
    "78679_EBC_B1_061_3","78679_EBC_B1_024_5", "78679_EBC_B1_061_2","78679_EBC_B1_061_6",
    "78679_EBC_B1_061_5","78679_EBC_B1_025_2", "78679_EBC_B1_061_4","78679_EBC_B1_025_5",
    "78679_EBC_B1_024_6","78679_EBC_B1_025_3", "78679_EBC_B1_024_4","78679_EBC_B1_061_1",
    "78679_EBC_B1_025_1","78679_EBC_B1_025_6", "78679_EBC_B1_025_4","78679_EBC_B1_Ref_5",
    "78679_EBC_B1_120_6","78679_EBC_B1_082_1", "78679_EBC_B1_005_3","78679_EBC_B1_005_2",
    "78679_EBC_B1_082_5","78679_EBC_B1_082_6", "78679_EBC_B1_120_3","78679_EBC_B1_120_5",
    "78679_EBC_B1_120_4","78679_EBC_B1_120_2", "78679_EBC_B1_082_2","78679_EBC_B1_082_4",
    "78679_EBC_B1_120_1","78679_EBC_B1_082_3", "78679_EBC_B1_005_1","78679_EBC_B1_Ref_6", "",
    "78679_EBC_B1_062_4","78679_EBC_B1_027_1", "78679_EBC_B1_062_1","78679_EBC_B1_027_3",
    "78679_EBC_B1_062_2","78679_EBC_B1_076_1", "78679_EBC_B1_027_6","78679_EBC_B1_062_3",
    "78679_EBC_B1_062_6","78679_EBC_B1_027_5", "78679_EBC_B1_122_2","78679_EBC_B1_062_5",
    "78679_EBC_B1_027_4","78679_EBC_B1_027_2", "78679_EBC_B1_Ref_7","","78679_EBC_B1_034_6",
    "78679_EBC_B1_047_3","78679_EBC_B1_047_6", "78679_EBC_B1_069_3","78679_EBC_B1_047_5",
    "78679_EBC_B1_047_4","78679_EBC_B1_047_1", "78679_EBC_B1_047_2","78679_EBC_B1_034_1",
    "78679_EBC_B1_034_2","78679_EBC_B1_064_4", "78679_EBC_B1_034_5","78679_EBC_B1_034_3", "",
    "78679_EBC_B1_Ref_8","78679_EBC_B1_088_5", "78679_EBC_B1_005_4","78679_EBC_B1_057_6",
    "78679_EBC_B1_088_6","78679_EBC_B1_057_4", "78679_EBC_B1_057_5","78679_EBC_B1_024_2",
    "78679_EBC_B1_088_3","78679_EBC_B1_088_4", "","78679_EBC_B1_057_2","78679_EBC_B1_057_1",
    "78679_EBC_B1_057_3","","78679_EBC_B1_088_1", "78679_EBC_B1_Ref_9"
  ),
  Type = rep(c(rep("Sample", 15), "Reference"), 9)
)
write.csv(metadata, "~/Desktop/OldProjects/TMT/Results/DataPackage4234/metadata.csv", quote = F, row.names = F)

# Run pipeline
tmt_pipeline(msnid, masic, metadata, "~/Desktop/OldProjects/TMT/Results/DataPackage4234/")
f_data <- read.csv("~/Desktop/OldProjects/TMT/Results/DataPackage4234/f_data.csv")
f_data$SampleNames <- paste0("X", f_data$SampleNames)
write.csv(f_data, "~/Desktop/OldProjects/TMT/Results/DataPackage4234/f_data.csv", quote = F, row.names = F)

e_data <- read.csv("~/Desktop/OldProjects/TMT/Results/DataPackage4234/e_data.csv")
e_meta <- read.csv("~/Desktop/OldProjects/TMT/Results/DataPackage4234/e_meta.csv")

pmartR::as.isobaricpepData(
  e_data, f_data, e_meta, "Peptide", "SampleNames", "Protein"
)

#######################
## DATA PACKAGE 4005 ##
#######################

# Load files
msnid <- readRDS("~/Desktop/OldProjects/TMT/4005_Tafess_LI/4005_Data_Package.RDS")$msnid
masic <- readRDS("~/Desktop/OldProjects/TMT/4005_Tafess_LI/4005_Data_Package.RDS")$masic
class(masic) <- c(class(masic), "masic_data")

metadata = data.frame(
  "PlexNames" = unique(masic$Dataset)[1:12],
  "IonChannelNames" = c(colnames(masic)[3:13], NA),
  "SampleNames" = c("Taffesse_LI_Mock_1",
                    "Taffesse_LI_SarsCov2_1",
                    "Taffesse_LI_Mock_2",
                    "Taffesse_LI_Mock_3",
                    "Taffesse_LI_SarsCov2_2",
                    "Taffesse_LI_SarsCov2_3",
                    "Taffesse_LI_Mock_4",
                    "Taffesse_LI_SarsCov2_4",
                    "Taffesse_LI_SarsCov2_5",
                    "Taffesse_LI_Mock_5",
                    "Reference",
                    "")
)
write.csv(metadata, "~/Desktop/OldProjects/TMT/Results/DataPackage4005/metadata.csv", quote = F, row.names = F)

# Run pipeline
tmt_pipeline(msnid, masic, metadata, "~/Desktop/OldProjects/TMT/Results/DataPackage4005/")

e_data <- read.csv("~/Desktop/OldProjects/TMT/Results/DataPackage4005/e_data.csv")
e_meta <- read.csv("~/Desktop/OldProjects/TMT/Results/DataPackage4005/e_meta.csv")

pmartR::as.isobaricpepData(
  e_data, f_data, e_meta, "Peptide", "SampleNames", "Protein"
)






#######################
## DATA PACKAGE 4235 ##
#######################

# Load files
msnid <- readRDS("~/Desktop/OldProjects/TMT/4235_Tafesse/4235_Data_Package.RDS")$msnid
masic <- readRDS("~/Desktop/OldProjects/TMT/4235_Tafesse/4235_Data_Package.RDS")$masic
class(masic) <- c(class(masic), "masic_data")



