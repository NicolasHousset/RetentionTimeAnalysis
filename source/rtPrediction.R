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

# Very different retention times on those projects
identified[, weird := (l_projectid > 2468 & l_projectid < 2476) |
           (l_projectid > 2461 & l_projectid < 2466)]

# Since we will calculate the median per protocol, let's put the weird stuff on a different protocol id
identified[weird == TRUE, l_protocolid := 99]

# Most common peptide notion depends on the protocol

# First get all the different protocols
identified[, l_protocolid.f := factor(l_protocolid)]
list_protocol <- levels(unique(identified[,l_protocolid.f]))



######
identified[, l_lcrunid:= as.character(l_lcrunid)]
setkey(identified, l_lcrunid, modified_sequence)
countsPerProject <- unique(identified)[, list(l_lcrunid,modified_sequence, l_protocolid.f)]
countsPerProject[, modified_sequence.f := factor(modified_sequence)]

# Tip : factor of modified_sequence has been computed with all the protocols.
# When summarizing it per protocol, peptides not present will still be here but with a counting of 0
# The order of summary(modified_sequence.f) will be the same (alphabetical)

nbProjPerPeptide <- summary(countsPerProject[, modified_sequence.f], maxsum = 1000000)
id_peptide <- 1:NROW(nbProjPerPeptide)
dt <- data.table(id_peptide)
dt[, modified_sequence := labels(nbProjPerPeptide)]

for (i in 1:NROW(list_protocol)){
  protocolSummary <- summary(countsPerProject[list_protocol[i]][, modified_sequence.f], maxsum = 1000000)
  stringExpression <- paste0("dt[, nbProjPepProtocol",list_protocol[i]," := -protocolSummary]")
  print(stringExpression)
  eval(parse(text=stringExpression))
}

for(i in 1:NROW(list_protocol)){
  # Ordering peptides by descending frequency in one given protocol
  eval(parse(text=paste0("setkey(dt, nbProjPepProtocol",list_protocol[i],")")))
  # Saving the order in the current setting
  eval(parse(text=paste0("dt[, rank_peptideProtocol",list_protocol[i]," := 1:NROW(nbProjPerPeptide)]")))
  # If peptide is not counted in one given protocol, set the rank to 0
  eval(parse(text=paste0("dt[nbProjPepProtocol",list_protocol[i],"==0, rank_peptideProtocol",list_protocol[i]," := 0L]")))
  # Returning the number of projects to a positive value
  eval(parse(text=paste0("dt[, nbProjPepProtocol", list_protocol[i]," := -nbProjPepProtocol", list_protocol[i],"]")))
}  

setkey(dt, modified_sequence)
setkey(identified, modified_sequence)
identified <- identified[dt]

identified[, grpProj := 1]
identified[l_projectid > 2062, grpProj := 2]
identified[l_projectid > 2499, grpProj := 3]

# Create an alphabetical-based index
id_peptide <- 1:NROW(nbProjPerPeptide)
dt <- data.table(id_peptide)
dt[, modified_sequence := labels(nbProjPerPeptide)]
dt[, nbProjPep := -nbProjPerPeptide]
setkey(dt, nbProjPep)
# Here, the index will depend of the number of projects in which each peptide appear
dt[, rank_peptide := 1:NROW(nbProjPerPeptide)]
dt[, nbProjPep := -nbProjPep]
dt[, id_peptide := NULL]
setkey(dt, modified_sequence)

setkey(identified, modified_sequence)
identified <- identified[dt]
rm(countsPerProject)

# test <- identified[(l_protocolid!=5 & l_protocolid!=11) | (nbProjPepProtocol5>4 | nbProjPepProtocol11>4)]

setkey(identified, l_lcrunid, modified_sequence, rtsec)
# To remove rt that have been identified more than once (otherwise, index are altered)
identified_subs <- unique(identified)
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
setkey(identified_subs, l_instrumentid, modified_sequence, l_protocolid, grpProj)
identified_subs[, q975_1 := quantile(rtsec, probs = 0.975), by = c("l_instrumentid", "modified_sequence", "l_protocolid", "grpProj")]
identified_subs[, q025_1 := quantile(rtsec, probs = 0.025), by = c("l_instrumentid", "modified_sequence", "l_protocolid", "grpProj")]
identified_subs[, q75_1 := quantile(rtsec, probs = 0.75), by = c("l_instrumentid", "modified_sequence", "l_protocolid", "grpProj")]
identified_subs[, q50_1 := quantile(rtsec, probs = 0.50), by = c("l_instrumentid", "modified_sequence", "l_protocolid", "grpProj")]
identified_subs[, q25_1 := quantile(rtsec, probs = 0.25), by = c("l_instrumentid", "modified_sequence", "l_protocolid", "grpProj")]
identified_subs[, wid95_1 := q975_1 - q025_1]
identified_subs[, wid50_1 := q75_1 - q25_1]
identified_subs[, QCD_1 := (q75_1 - q25_1) / (q75_1 + q25_1)]

# Statistics computed on rt measurements where index_rt2 < 5 : 4 most intense spectrum of the peptide per LC-run

identified_subs[index_rt2 <3, q975_2 := quantile(rtsec, probs = 0.975), by = c("l_instrumentid", "modified_sequence", "l_protocolid", "grpProj")]
identified_subs[index_rt2 <3, q025_2 := quantile(rtsec, probs = 0.025), by = c("l_instrumentid", "modified_sequence", "l_protocolid", "grpProj")]
identified_subs[index_rt2 <3, q75_2 := quantile(rtsec, probs = 0.75), by = c("l_instrumentid", "modified_sequence", "l_protocolid", "grpProj")]
identified_subs[index_rt2 <3, q50_2 := quantile(rtsec, probs = 0.50), by = c("l_instrumentid", "modified_sequence", "l_protocolid", "grpProj")]
identified_subs[index_rt2 <3, q25_2 := quantile(rtsec, probs = 0.25), by = c("l_instrumentid", "modified_sequence", "l_protocolid", "grpProj")]
identified_subs[index_rt2 <3, wid95_2 := q975_2 - q025_2]
identified_subs[index_rt2 <3, wid50_2 := q75_2 - q25_2]
identified_subs[index_rt2 <3, QCD_2 := (q75_2 - q25_2) / (q75_2 + q25_2)]

save(identified_subs, file = paste0(projectPath,"/data/identified_protocol.RData"), compression_level=1)
write.csv(identified_subs, file = paste0(projectPath,"/data/identified_protocol.csv"))

