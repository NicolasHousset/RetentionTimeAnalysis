library(data.table)
library(ggplot2)

projectPath <- "C:/Users/Nicolas Housset/Documents/RetentionTimeAnalysis"
load(file = paste0(projectPath,"/data/identified.RData"))


# Calculation of the max of identified retention times per project
identified[, l_projectid := as.character(l_projectid)]

rowNumbers <- identified[,j=list(.I[which.max(rtsec)]), by=c("l_projectid")][,j=V1]
maxRTperProject <- identified[rowNumbers, list(l_projectid,rtsec)]
maxRTperProject[, rtsecMax := rtsec]
maxRTperProject[, rtsec := NULL]
setkey(maxRTperProject, l_projectid)

setkey(identified, l_projectid)
identified <- identified[maxRTperProject]
# To filter out experiments with very high retention times 
identified <- identified[rtsecMax < 3001]

# The "most common peptide" notion is here project-based.
setkey(identified, l_projectid, modified_sequence)
countsPerProject <- unique(identified)
countsPerProject[, modified_sequence.f := factor(modified_sequence)]
nbProjPerPeptide <- summary(countsPerProject[, modified_sequence.f], maxsum = 1000000)
rm(countsPerProject)
# 675331 peptides

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
protsPerProject <- unique(identified)
protsPerProject[, accession.f := factor(accession)]
nbProjPerProtein <- summary(protsPerProject[, accession.f], maxsum = 500000)
rm(protsPerProject)
# 35989 proteins for the instrument 10

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

identified_subs <- identified[nbProjPep > 19]
setkey(identified_subs, l_lcrunid, modified_sequence, rtsec)
# To remove rt that have been identified more than once (otherwise, index are altered)
identified_subs <- unique(identified_subs)
convenient_vector <- 1:4000

# Add an index : 1 for the first time a peptide is encountered in a LC-run, 2 the second time, etc...
# convenient_vector is automatically shrinked to the appropriate size : that is very convenient :)
identified_subs[, index_rt1 := convenient_vector, by = c("l_lcrunid","modified_sequence")]
# Slightly different index : number of times the peptide is identified in the LC-run.
identified_subs[, size_rt := .N, by = c("l_lcrunid", "modified_sequence")]

identified_subs[,total_spectrum_intensity := -total_spectrum_intensity]
setkey(identified_subs, l_lcrunid, modified_sequence, total_spectrum_intensity)
identified_subs[, index_rt2 := convenient_vector, by = c("l_lcrunid","modified_sequence")]
identified_subs[,total_spectrum_intensity := -total_spectrum_intensity]


# Statistics computed for all the rt measurements
setkey(identified_subs, l_instrumentid, modified_sequence)
identified_subs[, max_1 := quantile(rtsec, probs = 1.0), by = c("l_instrumentid", "modified_sequence")]
identified_subs[, min_1 := quantile(rtsec, probs = 0.0), by = c("l_instrumentid", "modified_sequence")]
identified_subs[, q975_1 := quantile(rtsec, probs = 0.975), by = c("l_instrumentid", "modified_sequence")]
identified_subs[, q025_1 := quantile(rtsec, probs = 0.025), by = c("l_instrumentid", "modified_sequence")]
identified_subs[, q75_1 := quantile(rtsec, probs = 0.75), by = c("l_instrumentid", "modified_sequence")]
identified_subs[, q25_1 := quantile(rtsec, probs = 0.25), by = c("l_instrumentid", "modified_sequence")]
identified_subs[, wid100_1 := max_1 - min_1]
identified_subs[, wid95_1 := q975_1 - q025_1]
identified_subs[, wid50_1 := q75_1 - q25_1]
identified_subs[, QCD_1 := (q75_1 - q25_1) / (q75_1 + q25_1)]


# Statistics computed on rt measurements where index_rt2 < 5 : 4 most intense spectrum of the peptide per LC-run

identified_subs[index_rt2 <5, max_2 := quantile(rtsec, probs = 1.0), by = c("l_instrumentid", "modified_sequence")]
identified_subs[index_rt2 <5, min_2 := quantile(rtsec, probs = 0.0), by = c("l_instrumentid", "modified_sequence")]
identified_subs[index_rt2 <5, q975_2 := quantile(rtsec, probs = 0.975), by = c("l_instrumentid", "modified_sequence")]
identified_subs[index_rt2 <5, q025_2 := quantile(rtsec, probs = 0.025), by = c("l_instrumentid", "modified_sequence")]
identified_subs[index_rt2 <5, q75_2 := quantile(rtsec, probs = 0.75), by = c("l_instrumentid", "modified_sequence")]
identified_subs[index_rt2 <5, q25_2 := quantile(rtsec, probs = 0.25), by = c("l_instrumentid", "modified_sequence")]
identified_subs[index_rt2 <5, wid100_2 := max_2 - min_2]
identified_subs[index_rt2 <5, wid95_2 := q975_2 - q025_2]
identified_subs[index_rt2 <5, wid50_2 := q75_2 - q25_2]
identified_subs[index_rt2 <5, QCD_2 := (q75_2 - q25_2) / (q75_2 + q25_2)]


setkey(identified_subs, l_lcrunid, modified_sequence, rtsec)

save(identified_subs, file = "C:/Users/Nicolas Housset/Documents/RetentionTimeAnalysis/data/identified_subs.RData", compression_level = 1)

write.csv(identified_subs, file = paste0(projectPath,"/data/identified_sample.csv"))







