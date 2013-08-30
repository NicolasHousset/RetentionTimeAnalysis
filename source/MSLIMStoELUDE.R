# R script to convert MS-LIMS modified sequence into ELUDE modified sequence
# Main difference is how modifications are handled : <> for MS-LIMS, [] for ELUDE
# ELUDE can handle arbitrary modifications but they are recognized with []

library(data.table)
library(ggplot2)

projectPath <- "C:/Users/Nicolas Housset/Documents/RetentionTimeAnalysis"


load(file = paste0(projectPath,"/data/identified_protocol.RData"))
setkey(identified_subs, modified_sequence)

rtPeptide <- identified_subs[, list(modified_sequence, l_protocolid, grpProj, q50_2)]
rtPeptide <- rtPeptide[l_protocolid=='11'][grpProj==3]

# For now, only keeping the non-modified peptides
rtPeptide[, matching := grepl("NH2-[ABCDEFGHIJKLMNOPQRSTUVWXYZ]*-COOH", modified_sequence)]
rtPeptide <- rtPeptide[matching==TRUE]
setkey(rtPeptide, modified_sequence)
rtPeptide <- unique(rtPeptide)
rtPeptide[, modified_sequence := sub("NH2-", "", modified_sequence)]
rtPeptide[, modified_sequence := sub("-COOH", "", modified_sequence)]
rtPeptide <- rtPeptide[, list(modified_sequence, q50_2)]

# rtPeptide[, q50_2 := q50_2 / 60]
# summary(rtPeptide, q50_2)
rtPeptideFlush <- sample(rtPeptide[, modified_sequence])
setkey(rtPeptide, modified_sequence)
rtPeptide <- rtPeptide[rtPeptideFlush]


rtPeptideTrain <- rtPeptide[1:1000]
rtPeptideTest <- rtPeptide[1:NROW(rtPeptide)]

write.table(rtPeptideTrain, file=paste0(projectPath, "/data/rtPeptideTrain.txt"), quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
write.table(rtPeptideTest, file=paste0(projectPath, "/data/rtPeptideTest.txt"), quote = FALSE, , sep="\t", row.names = FALSE, col.names = FALSE)

# elude-running <- "elude -v 5 -t "C:\Users\Nicolas Housset\Documents\RetentionTimeAnalysis\data\rtPeptideTrain.txt" -e "C:\Users\Nicolas Housset\Documents\RetentionTimeAnalysis\data\rtPeptideTest.txt" -s "C:\Users\Nicolas Housset\Documents\RetentionTimeAnalysis\data\modelTest" -r "C:\Users\Nicolas Housset\Documents\RetentionTimeAnalysis\data\retentionIndexTest" -y -p"

results <- data.table(read.table(file=paste0(projectPath, "/data/predictions.out"), header = TRUE, sep = "\t"))

ggplot(results, aes(Predicted_RT, Observed_RT)) + geom_point(alpha=(1/2))
ggplot(results, aes(Observed_RT, Predicted_RT - Observed_RT)) + geom_point(alpha=(1/2))
ggplot(rtPeptideTrain, aes("A",q50_2)) + geom_jitter(alpha=(1/2))

results[, modified_sequence := paste0("NH2-", Peptide, "-COOH")]

rtPeptide <- identified_subs[l_protocolid=='11'][grpProj==3]

