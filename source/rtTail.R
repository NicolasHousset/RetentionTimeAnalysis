# Something new, I don't really know what, yet...
# Let me tell you a story, the tale of the tails : once upon a time...


library(data.table)
library(ggplot2)

projectPath <- "C:/Users/Nicolas Housset/Documents/RetentionTimeAnalysis"
load(file = paste0(projectPath,"/data/identified.RData"))

id_10 <- identified[l_instrumentid == 10]
rm(identified)
gc()

# Calculation of the max of identified retention times per project
id_10[, l_projectid := as.character(l_projectid)]

rowNumbers <- id_10[,j=list(.I[which.max(rtsec)]), by=c("l_projectid")][,j=V1]
maxRTperProject <- id_10[rowNumbers, list(l_projectid,rtsec)]
maxRTperProject[, rtsecMax := rtsec]
maxRTperProject[, rtsec := NULL]
setkey(maxRTperProject, l_projectid)

setkey(id_10, l_projectid)
id_10 <- id_10[maxRTperProject]
# To filter out experiments with very high retention times 
id_10 <- id_10[rtsecMax < 3001]

# The "most common peptide" notion is here project-based.
setkey(id_10, l_projectid, modified_sequence)
countsPerProject <- unique(id_10)
countsPerProject[, modified_sequence.f := factor(modified_sequence)]
nbProjPerPeptide <- summary(countsPerProject[, modified_sequence.f], maxsum = 500000)
# 270992 peptides for instrument 10

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

setkey(id_10, modified_sequence)
id_10 <- id_10[dt]

# We repeat this part on the protein level
setkey(id_10, l_projectid, accession)
protsPerProject <- unique(id_10)
protsPerProject[, accession.f := factor(accession)]
nbProjPerProtein <- summary(protsPerProject[, accession.f], maxsum = 500000)
# 25165 proteins for the instrument 10

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

setkey(id_10, accession)
id_10 <- id_10[dt]

id_10_subs <- id_10[nbProjPep > 9]
setkey(id_10_subs, l_lcrunid, modified_sequence, rtsec)
# To remove rt that have been identified more than once (otherwise, index are altered)
id_10_subs <- unique(id_10_subs)
convenient_vector <- 1:2000

# Add an index : 1 for the first time a peptide is encountered in a LC-run, 2 the second time, etc...
# convenient_vector is automatically shrinked to the appropriate size : that is very convenient :)
id_10_subs[, index_rt1 := convenient_vector, by = c("l_lcrunid","modified_sequence")]
# Slightly different index : number of times the peptide is identified in the LC-run.
id_10_subs[, size_rt := .N, by = c("l_lcrunid", "modified_sequence")]

id_10_subs[,total_spectrum_intensity := -total_spectrum_intensity]
setkey(id_10_subs, l_lcrunid, modified_sequence, total_spectrum_intensity)
id_10_subs[, index_rt2 := convenient_vector, by = c("l_lcrunid","modified_sequence")]
id_10_subs[,total_spectrum_intensity := -total_spectrum_intensity]


# Statistics computed for all the rt measurements
setkey(id_10_subs, modified_sequence)
id_10_subs[, max_1 := quantile(rtsec, probs = 1.0), by = modified_sequence]
id_10_subs[, min_1 := quantile(rtsec, probs = 0.0), by = modified_sequence]
id_10_subs[, q975_1 := quantile(rtsec, probs = 0.975), by = modified_sequence]
id_10_subs[, q025_1 := quantile(rtsec, probs = 0.025), by = modified_sequence]
id_10_subs[, q75_1 := quantile(rtsec, probs = 0.75), by = modified_sequence]
id_10_subs[, q25_1 := quantile(rtsec, probs = 0.25), by = modified_sequence]
id_10_subs[, wid100_1 := max_1 - min_1]
id_10_subs[, wid95_1 := q975_1 - q025_1]
id_10_subs[, wid50_1 := q75_1 - q25_1]
id_10_subs[, QCD_1 := (q75_1 - q25_1) / (q75_1 + q25_1)]


# Statistics computed on rt measurements where index_rt2 < 5 : 4 most intense spectrum of the peptide per LC-run

sub_sub <- id_10_subs[index_rt2 < 5]
setkey(sub_sub, modified_sequence)
id_10_subs[index_rt2 <5, max_2 := quantile(rtsec, probs = 1.0), by = modified_sequence]
id_10_subs[index_rt2 <5, min_2 := quantile(rtsec, probs = 0.0), by = modified_sequence]
id_10_subs[index_rt2 <5, q975_2 := quantile(rtsec, probs = 0.975), by = modified_sequence]
id_10_subs[index_rt2 <5, q025_2 := quantile(rtsec, probs = 0.025), by = modified_sequence]
id_10_subs[index_rt2 <5, q75_2 := quantile(rtsec, probs = 0.75), by = modified_sequence]
id_10_subs[index_rt2 <5, q25_2 := quantile(rtsec, probs = 0.25), by = modified_sequence]
id_10_subs[index_rt2 <5, wid100_2 := max_2 - min_2]
id_10_subs[index_rt2 <5, wid95_2 := q975_2 - q025_2]
id_10_subs[index_rt2 <5, wid50_2 := q75_2 - q25_2]
id_10_subs[index_rt2 <5, QCD_2 := (q75_2 - q25_2) / (q75_2 + q25_2)]


setkey(id_10_subs, l_lcrunid, modified_sequence, rtsec)

setkey(id_10_subs, l_lcrunid, modified_sequence)

write.csv(id_10_subs, file = paste0(projectPath,"/data/id_10_sample.csv"))





test_stats <- id_10_subs[, list(quantile(rtsec, probs = 1.0),
                    quantile(rtsec, probs = 0.0),
                    quantile(rtsec, probs = 0.975),
                    quantile(rtsec, probs = 0.025),
                    quantile(rtsec, probs = 0.75),
                    quantile(rtsec, probs = 0.25))
                    , by = c("modified_sequence")]

test_stats[, wid100 := V1 - V2]
test_stats[, wid95 := V3 - V4]
test_stats[, wid50 := V5 - V6]
test_stats[, QCD := (V5 - V6) / (V5 + V6)]

plot_stat <- ggplot(test_stats, aes(modified_sequence, QCD))
plot_stat + geom_point(alpha=1/2)

# Building of the aggregated dataset
setkey(id_10, modified_sequence)
pep_stat <- unique(id_10)[, list(modified_sequence,
                                 id_peptide,
                                 nbProjPep,
                                 rank_peptide)]

id_10[, N := as.numeric(N)]
setkey(id_10, N)

a <- as.data.frame(table(id_10[, N]))

ggplot(a, aes(Var1, Freq)) + geom_point()
