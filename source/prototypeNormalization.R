# Something huge...

# How can we bring closer the retention times of different set of lc-run ? I'll try stuff on two sets of lc-run I found with Spotfire

library(data.table)
library(ggplot2)

projectPath <- "C:/Users/Nicolas Housset/Documents/RetentionTimeAnalysis"
load(file = paste0(projectPath,"/data/identified.RData"))

identified[, l_lcrunid := as.character(l_lcrunid)]
setkey(identified, l_lcrunid)

identified <- identified[c(as.character(87371:87376),
                     as.character(88556:88562))]

# The "most common peptide" notion is here lcrun-based.
setkey(identified, l_lcrunid, modified_sequence)
countsPerProject <- unique(identified)[, list(l_lcrunid,modified_sequence)]
countsPerProject[, modified_sequence.f := factor(modified_sequence)]
nbProjPerPeptide <- summary(countsPerProject[, modified_sequence.f], maxsum = 1000000)
rm(countsPerProject)
# 2774 peptides (22/08/2013)

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
setkey(identified, l_lcrunid, accession)
protsPerProject <- unique(identified)[, list(l_lcrunid, accession)]
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

setkey(identified, l_lcrunid, modified_sequence, rtsec)
# To remove rt that have been identified more than once (otherwise, index are altered)
identified <- unique(identified)
convenient_vector <- 1:4000

# Add an index : 1 for the first time a peptide is encountered in a LC-run, 2 the second time, etc...
# convenient_vector is automatically shrinked to the appropriate size : that is very convenient :)
identified[, index_rt1 := convenient_vector, by = c("l_lcrunid","modified_sequence")]
# Slightly different index : number of times the peptide is identified in the LC-run.
identified[, size_rt := .N, by = c("l_lcrunid", "modified_sequence")]

identified[,total_spectrum_intensity := -total_spectrum_intensity]
setkey(identified, l_lcrunid, modified_sequence, total_spectrum_intensity)
identified[, index_rt2 := convenient_vector, by = c("l_lcrunid","modified_sequence")]
identified[,total_spectrum_intensity := -total_spectrum_intensity]

identified[, grpLC := (as.numeric(l_lcrunid) > 88555)]

table(identified[, grpLC])

test <- identified[, nbId := .N * (index_rt2 == 1), by = c("grpLC", "modified_sequence")]
test <- identified[, nbId := .N, by = c("grpLC", "modified_sequence")]

table(test[, nbId, by = grpLC])

part1 <- identified[grpLC == FALSE]
part2 <- identified[grpLC == TRUE]

part1[, rt1 := quantile(rtsec, probs = 0.5), by = modified_sequence]
part2[, rt2 := quantile(rtsec, probs = 0.5), by = modified_sequence]

setkey(part1, modified_sequence)
setkey(part2, modified_sequence)
fusion1 <- unique(part1)[, j = list(sequence, rt1, modified_sequence)]
fusion2 <- unique(part2)[, j = list(sequence, rt2, modified_sequence)]

setkey(fusion1, sequence, modified_sequence)
setkey(fusion2, sequence, modified_sequence)

fusion <- merge(fusion1, fusion2)

plot_fusion <- ggplot(fusion, aes(x=rt1, y=rt2))
plot_fusion + geom_point()

fusion[, rt1_adjust := rt1 - 800]
fusion[, rt2_adjust := rt2 - 800]

fusion[,rtdiff := rt2 - rt1]

setkey(fusion, rtdiff)

plot_diff <- ggplot(fusion, aes(x=rt1,y=rtdiff))
plot_diff + geom_point()

fusion_2 <- fusion[rtdiff > 40]

plot_diff <- ggplot(fusion_2, aes(x=rt1,y=rtdiff))
plot_diff + geom_point()

fusion_2[,rt2trans := 0.951223 * rt2 - 56.760450]
plot_diff <- ggplot(fusion_2, aes(x=rt1,y=rt2trans))
plot_diff + geom_point()
