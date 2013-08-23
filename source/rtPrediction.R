# Again a reboot...
# We start to get used to it right ?
# Anyway...

library(data.table)
library(ggplot2)

projectPath <- "C:/Users/Nicolas Housset/Documents/RetentionTimeAnalysis"
load(file = paste0(projectPath,"/data/identified.RData"))

identified <- identified[l_instrumentid == 10]

# To filter out experiments with very high retention times 
identified[, l_projectid := as.character(l_projectid)]
rowNumbers <- identified[,j=list(.I[which.max(rtsec)]), by=c("l_projectid")][,j=V1]
maxRTperProject <- identified[rowNumbers, list(l_projectid,rtsec)]
maxRTperProject[, rtsecMax := rtsec]
maxRTperProject[, rtsec := NULL]
setkey(maxRTperProject, l_projectid)
setkey(identified, l_projectid)
identified <- identified[maxRTperProject]
identified <- identified[rtsecMax < 3001]

identified[, weird := (l_projectid > 2468 & l_projectid < 2476) |
           (l_projectid > 2461 & l_projectid < 2466)]

# Since we will calculate the median per protocol, let's put the weird stuff on a different protocol id
identified[weird == TRUE, l_protocolid := 99]

# Most common peptide notion depends on the protocol

setkey(identified, l_projectid, modified_sequence)
countsPerProject <- unique(identified)[, list(l_projectid,modified_sequence, l_protocolid)]
countsPerProject[, modified_sequence.f := factor(modified_sequence)]

nbProjPerPeptide <- summary(countsPerProject[l_protocolid == 1, modified_sequence.f], maxsum = 100000)
nbProjPerPeptide <- countsPerProject[, summary(modified_sequence.f, maxsum = 100), by = l_protocolid]

labels(nbProjPerPeptide[, V1])
rm(countsPerProject)

# The "most common peptide" notion is here project-based.
setkey(identified, l_projectid, modified_sequence)
countsPerProject <- unique(identified)[, list(l_projectid,modified_sequence)]
countsPerProject[, modified_sequence.f := factor(modified_sequence)]
nbProjPerPeptide <- summary(countsPerProject[, modified_sequence.f], maxsum = 1000000)
rm(countsPerProject)
# 677930 peptides (21/08/2013)

# Create an alphabetical-based index
id_peptide <- 1:NROW(nbProjPerPeptide)
dt <- data.table(id_peptide)
dt[, modified_sequence := labels(nbProjPerPeptide)]
dt[, nbProjPep := -nbProjPerPeptide]
setkey(dt, nbProjPep)
# Here, the index will depend of the number of projects in which each peptide appear
dt[, rank_peptide := 1:NROW(nbProjPerPeptide)]
dt[, nbProjPep := -nbProjPep]
setkey(dt, modified_sequence)

setkey(identified, modified_sequence)
identified <- identified[dt]

# We repeat this part on the protein level
setkey(identified, l_projectid, accession)
protsPerProject <- unique(identified)[, list(l_projectid, accession)]
protsPerProject[, accession.f := factor(accession)]
nbProjPerProtein <- summary(protsPerProject[, accession.f], maxsum = 500000)
rm(protsPerProject)
# 54402 proteins (21/08/2013)

# Create an alphabetical-based index
id_protein <- 1:NROW(nbProjPerProtein)
dt <- data.table(id_protein)
dt[, accession := labels(nbProjPerProtein)]
dt[, nbProjProt := -nbProjPerProtein]
setkey(dt, nbProjProt)
# Here, the index will depend of the number of projects in which each protein appear
dt[, rank_protein := 1:NROW(nbProjPerProtein)]
dt[, nbProjProt := -nbProjProt]
setkey(dt, accession)

setkey(identified, accession)
identified <- identified[dt]