# Something new, I don't really know what, yet...
# Let me tell you a story, the tale of the tails : once upon a time...


library(data.table)
library(ggplot2)

projectPath <- "C:/Users/Nicolas Housset/Documents/RetentionTimeAnalysis"
load(file = paste0(projectPath,"/data/identified.RData"))

id_10 <- identified[l_instrumentid == 10]
rm(identified)
gc()

# Calculation of the min and max of identified retention times per project
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

id_10_subs <- id_10[nbProjPep > 39]
setkey(id_10_subs, l_lcrunid, modified_sequence, rtsec)
convenient_vector <- 1:2000

id_10_subs[, index_rt := convenient_vector, by = c("l_lcrunid","modified_sequence")]
id_10_subs
test_stats <- id_10_subs[, list(quantile(rtsec, probs = 1.0),
                    quantile(rtsec, probs = 0.0),
                    quantile(rtsec, probs = 0.975),
                    quantile(rtsec, probs = 0.025),
                    quantile(rtsec, probs = 0.75),
                    quantile(rtsec, probs = 0.25))
                    , by = c("modified_sequence")]

test_stats <- test_stats[, wid100 := V1 - V2]
test_stats <- test_stats[, wid95 := V3 - V4]
test_stats <- test_stats[, wid50 := V5 - V6]
test_stats <- test_stats[, IQR := (V3 - V4) / (V3 + V4)]

plot_stat <- ggplot(test_stats, aes(modified_sequence, IQR))
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
