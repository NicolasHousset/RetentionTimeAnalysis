# A R script to export identifications to Spotfire (or any other statistical software, actually)

library(data.table)
library(ggplot2)

projectPath <- "C:/Users/Nicolas Housset/Documents/RetentionTimeAnalysis"
load(file = paste0(projectPath,"/data/identified.RData"))

# Adding the number of identifications of one peptide per LC-RUN
setkey(identified, l_lcrunid, modified_sequence)

numbIDperLCRun <- identified[, .N, by = c("l_lcrunid", "modified_sequence")]
setkey(numbIDperLCRun, l_lcrunid, modified_sequence)

identified <- identified[numbIDperLCRun]

# Calculation of the min and max of identified retention times per project
identified[, l_projectid := as.character(l_projectid)]
setkey(identified, l_projectid)
rowNumbers <- identified[,j=list(.I[which.max(rtsec)]), by=c("l_projectid")][,j=V1]
maxRTperProject <- identified[rowNumbers, list(l_projectid,rtsec)]
maxRTperProject[, rtsecMax := rtsec]
maxRTperProject[, rtsec := NULL]
setkey(maxRTperProject, l_projectid)

rowNumbers <- identified[,j=list(.I[which.min(rtsec)]), by=c("l_projectid")][,j=V1]
minRTperProject <- identified[rowNumbers, list(l_projectid,rtsec)]
minRTperProject[, rtsecMin := rtsec]
minRTperProject[, rtsec := NULL]
setkey(minRTperProject, l_projectid)

setkey(identified, l_projectid)
identified <- identified[minRTperProject][maxRTperProject]

identified[, rtsecNorm0 := 3000 * rtsec / rtsecMax]
identified[, rtsecNorm1 := 3000 * rtsec / rtsecMax - rtsecMin]
identified[, rtsecNorm2 := 3000 * (rtsec - rtsecMin) / rtsecMax]

# The "most common peptide" notion is here project-based.
setkey(identified, l_projectid, modified_sequence)
countsPerProject <- unique(identified)
countsPerProject[, modified_sequence.f := factor(modified_sequence)]
nbProjPerPeptides <- summary(countsPerProject[, modified_sequence.f], maxsum = 500000)

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

# Test code to check if everything went as expected
subs_id <- identified[rank_peptide == 1]

identified[, list(rank_peptide, nbProj)]
setkey(identified, rank_peptide)
identified