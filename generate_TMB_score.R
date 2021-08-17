#!/usr/bin/Rscript
myArgs = commandArgs(trailingOnly=TRUE)

#load tidyverse package
library(tidyverse)

## usage:
## Rscript generate_TMB_score.R 
# --input_variants          - variant file
# --roi_size                - size of ROI (in bp) over which coding variants are called

# PGDx (2021)
# Laurel Keefer and Sam Angiuoli

# Reading user parameters
# input_variants
argVals <- getOpts(myArgs)

#check for variant input file
if(!argSpecified("--input_variants")){
  write("No variant file specified!", stderr())
  quit(status=1)
}
inVariants <- argVals["--input_variants"]
if(!file.exists(inVariants)){
  write(sprintf("Variants file does not exist!\n\t%s",inVariants), stderr())
  quit(status=1)
}


# check that ROI size is specified
if(!argSpecified("--roi_size")){
  write("No ROI size specified!  Please provide size of coding ROI in bp", stderr())
  quit(status=1)
}
ROI_size <- argVals["--roi_size"]
if(ROI_size%%1 !=0 ){
  write(sprintf("ROI size does not appear to be an integer\n\t%s",ROI_size), stderr())
  quit(status=1)
}else if(ROI_size == 0) {
  write(sprintf("ROI size is invalid!\n\t%s",ROI_size), stderr())
  quit(status=1)
}

vars <- read.delim(inVariants)

#filter variants
filtered_vars <- vars %>% filter(mut_reads >= 4) %>%
        filter(variant_allele_fraction > 0.05) %>%
        filter(coding_noncoding == "coding") %>%
        filter(germline_status != "germline") %>%
        filter(hotspot_status == "non_hotspot") %>%
        filter(quality_pass_fail == "pass")

#calculate TMB, transform to muts/Mbp
TMB_score <- nrow(filtered_vars) / (ROI_size / 1000000)

write.table(filtered_vars, "TMB_filtered_candidate_variants.txt", sep = "\t", row.names = F, quote = F)
write.table(data.frame(TMB_Score = TMB_score, n_Variants = nrow(filtered_vars)), file = "TMB_score.txt", row.names = F, sep = "\t", quote = F)
