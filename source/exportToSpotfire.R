# A R script to export identifications to Spotfire (or any other statistical software, actually)

library(data.table)
library(ggplot2)

projectPath <- "C:/Users/Nicolas Housset/Documents/RetentionTimeAnalysis"
load(file = paste0(projectPath,"/data/identified.RData"))

# The "most common peptide" notion is here project-based.
setkey(identified, l_projectid, modified_sequence)
countsPerProject <- unique(identified)
countsPerProject[, modified_sequence.f := factor(modified_sequence)]
nbProjPerPeptides <- summary(countsPerProject[, modified_sequence.f], maxsum = 500000)
nbProjPerPeptides[1]

# Create an alphabetical-based index
id_peptide <- 1:NROW(nbProjPerPeptides)
dt <- data.table(id_peptide)
dt[, modified_sequence := labels(nbProjPerPeptides)]
dt[, nbProj := -nbProjPerPeptides]
setkey(dt, nbProj)
# Here, the index will depend of the number of projects in which each peptide appear
dt[, rank_peptide := 1:NROW(nbProjPerPeptides)]
dt[, nbProj := -nbProj]
setkey(dt, modified_sequence)

setkey(identified, modified_sequence)
identified <- identified[dt]

write.csv(identified, file = paste0(projectPath,"/data/identified.csv"))

subs_id <- identified[rank_peptide == 1]

identified[, list(rank_peptide, nbProj)]
setkey(identified, rank_peptide)
identified