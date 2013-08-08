projectPath <- "C:/Users/Nicolas Housset/Documents/RetentionTimeAnalysis"
library(data.table)
library(ggplot2)

load(file = paste0(projectPath,"/data/identified.RData"))

# We would like to count the number of times each peptide is identified
# However, we need to be careful that a peptide can be identified many times in a LC-run
setkey(identified, l_lcrunid, modified_sequence)
counts <- unique(identified)

setkey(identified, l_projectid, modified_sequence)
countsPerProject <- unique(identified)
countsPerProject[, modified_sequence.f := factor(modified_sequence)]
nbProjPerPeptides <- summary(countsPerProject[, modified_sequence.f], maxsum = 2000)

counts[, modified_sequence.f := factor(modified_sequence)]
liste_peptides <- summary(counts[, modified_sequence.f], maxsum = 2000)

setkey(identified,modified_sequence)
one_peptide <- identified[labels(nbProjPerPeptides)[1:30]]

table(one_peptide[, modified_sequence])

setkey(one_peptide, l_lcrunid, modified_sequence)
one_peptide_reduced <- one_peptide[one_peptide[,j=list(.I[which.max(total_spectrum_intensity)]), by=c("l_lcrunid","modified_sequence")][,j=V1]]

table(one_peptide_reduced[, modified_sequence])
setkey(countsPerProject, modified_sequence)
nbId_2 <- table(countsPerProject[labels(nbProjPerPeptides)[1:50]][,modified_sequence])


plot_all <- ggplot(one_peptide, aes(modified_sequence, rtsec, colour=factor(l_projectid)))
plot_all + geom_jitter(alpha = 1/2, size = 1.5) + ylim(0,3000) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_maxMS <- ggplot(one_peptide_reduced, aes(modified_sequence, rtsec, colour=factor(l_projectid)))
plot_maxMS + geom_jitter(alpha = 1/2, size = 1.5) + ylim(0,3000) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

one_peptide[, l_projectid := as.character(l_projectid)]
one_peptide_reduced[, l_projectid := as.character(l_projectid)]
one_peptide_reduced[, l_projectid.f := factor(l_projectid)]
liste_projects <- summary(one_peptide_reduced[, l_projectid.f])
setkey(one_peptide, l_projectid)
setkey(one_peptide_reduced, l_projectid)

plot_all <- ggplot(one_peptide[labels(liste_projects)[1:6]], aes(modified_sequence, rtsec, colour=factor(l_projectid)))
plot_all + geom_jitter(alpha = 1, size = 1.5) + ylim(0,3000) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_maxMS <- ggplot(one_peptide_reduced[labels(liste_projects)[1:6]], aes(modified_sequence, rtsec, colour=factor(l_projectid)))
plot_maxMS + geom_jitter(alpha = 1, size = 1.5) + ylim(0,3000) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# What if we remove the most common projects (high number of identifications)
# Hypothesis is that the peptides are identified maaaany times inside each LC-run
setkey(one_peptide, l_projectid)

one_peptide[, l_projectid:= as.character(l_projectid)]
setkey(one_peptide, l_projectid)
plot_all <- ggplot(one_peptide[labels(liste_projects)[7]], aes(modified_sequence, rtsec, colour=factor(l_projectid)))


nbId <- table(one_peptide[, l_projectid])

nbIdPerRun <- table(one_peptide_reduced[, l_projectid])

rapport <- nbId / nbIdPerRun

rapportdf <- as.data.frame(rapport)

labels(nbId)



plot <- ggplot(rapportdf,aes(Var1,Freq, colour=factor(l_instrumentid)))

plot+geom_point()


setkey(identified, l_lcrunid)
countsLC <- unique(identified)
nbLCPerProject <- table(countsLC[, l_projectid])
nbLCPerProject