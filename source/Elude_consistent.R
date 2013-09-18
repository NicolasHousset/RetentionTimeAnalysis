# Trying to see what happens if we take more consistent data

library(data.table)
library(ggplot2)

projectPath <- "C:/Users/Nicolas Housset/Documents/RetentionTimeAnalysis"


load(file = paste0(projectPath,"/data/identified_protocol.RData"))
setkey(identified_subs, modified_sequence)
rtPeptide <- identified_subs[nbProjPepProtocol11 > 4][, list(l_projectid, l_lcrunid, modified_sequence, index_rt2, l_protocolid, grpProj, rtsec, q50_2)]
rtPeptide <- rtPeptide[l_protocolid=='11'][grpProj==2]
setkey(rtPeptide, modified_sequence)
rtPepUnique <- unique(rtPeptide)


# Yes, I want to get the modifications. And this is complicated.

mod_pep <- rtPepUnique[, modified_sequence]
modif <- gregexpr("[ACDEFGHIKLMNPQRSTVWY]{1}<[ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789]+[*]*>", mod_pep, ignore.case = TRUE)
modifN <- gregexpr("^[ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789]+[-]", mod_pep, ignore.case = TRUE)
modifC <- gregexpr("[-][ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789]+$", mod_pep, ignore.case = TRUE)

listmod <- vector("list", NROW(modif))
for (i in (1:NROW(modif))){
  listmod[[i]] <- vector("character",4)
  for (j in (1:NROW(modif[[i]]))){
    if(modif[[i]][[j]] > -1){
      listmod[[i]][[j]] <- substr(mod_pep[[i]], modif[[i]][[j]], modif[[i]][[j]]+attr(modif[[i]],"match.length")[[j]]-1)
    }
  }
}

listmodN <- vector("list", NROW(modifN))
for (i in (1:NROW(modifN))){
  listmodN[[i]] <- vector("character",1)
  for (j in (1:NROW(modifN[[i]]))){
    if(modifN[[i]][[j]] > -1){
      listmodN[[i]][[j]] <- substr(mod_pep[[i]], modifN[[i]][[j]], modifN[[i]][[j]]+attr(modifN[[i]],"match.length")[[j]]-1)
    }
  }
}

listmodC <- vector("list", NROW(modifC))
for (i in (1:NROW(modifC))){
  listmodC[[i]] <- vector("character",1)
  for (j in (1:NROW(modifC[[i]]))){
    if(modifC[[i]][[j]] > -1){
      listmodC[[i]][[j]] <- substr(mod_pep[[i]], modifC[[i]][[j]], modifC[[i]][[j]]+attr(modifC[[i]],"match.length")[[j]]-1)
    }
  }
}

rtPepUnique[, listMod := listmod]
rtPepUnique[, mod1 := as.character(lapply(listmod, "[[", 1))]
rtPepUnique[, mod2 := as.character(lapply(listmod, "[[", 2))]
rtPepUnique[, mod3 := as.character(lapply(listmod, "[[", 3))]
rtPepUnique[, mod4 := as.character(lapply(listmod, "[[", 4))]
rtPepUnique[, modN := as.character(lapply(listmodN, "[[", 1))]
rtPepUnique[, modC := as.character(lapply(listmodC, "[[", 1))]

setkey(rtPepUnique, modified_sequence)
setkey(rtPeptide, modified_sequence)
rtPeptide <- rtPepUnique[rtPeptide]

# Set some condition
rtPeptide <- rtPeptide[mod2==""]
rtPeptide <- rtPeptide[mod1=="" | mod1=="K<prop>" | mod1=="K<propC13>"]

rtPeptide[substr(modified_sequence, 1, 4) == "Ace-", test := paste0(substr(modified_sequence, 5,5),"[Ace]",substr(modified_sequence, 6, nchar(modified_sequence)))]
rtPeptide[substr(modified_sequence, 1, 4) == "NH2-", test := substr(modified_sequence, 5,nchar(modified_sequence))]
rtPeptide[substr(modified_sequence, 1, 5) == "prop-", test := paste0(substr(modified_sequence, 6,6),"[prop]",substr(modified_sequence, 7, nchar(modified_sequence)))]
rtPeptide[substr(modified_sequence, 1, 7) == "propC13-", test := paste0(substr(modified_sequence, 8,8),"[prop]",substr(modified_sequence, 9, nchar(modified_sequence)))]

rtPeptide[, test := sub("<","[",test)]
rtPeptide[, test := sub(">","]",test)]
rtPeptide[, test := sub("-COOH", "", test)]
rtPeptide <- rtPeptide[!is.na(test)]
# Sometimes there is a modification at the N-terminal and at the first amino acide : not sure if ELUDE can handle that.
# This is very uncommon but for now we'll just remove those cases.
rtPeptide[, error := grepl("\\]\\[", test)]
rtPeptide <- rtPeptide[error == FALSE]

# This line is useful if selecting customized project numbers, recomputes the median
# rtPeptide[index_rt2 <2, q50_3 := quantile(rtsec, probs = 0.50), by = c("modified_sequence")]

rtPeptide <- rtPeptide[index_rt2 <2] 
rtPeptide <- rtPeptide[, list(test, q50_2, rtsec)]

setkey(rtPeptide, test)
rtPeptideUnique <- unique(rtPeptide)

# Stratified sampling : we want to have  more peptides selected at early and late retention times, less in the middle
early <- as.character(7:10)
late <- as.character(18:24)
middle <- as.character(11:17)
earlyCoefficient <- 1
lateCoefficient <- 1
middleCoefficient <- 1
threshold <- 0.1

breaks <- (7:25)*100
rtPeptideUnique[, rtCent := cut(x=q50_2, breaks=((7:25)*100), labels=7:24)]
rtPeptideUnique[, rtCent := as.character(rtCent)]
# We generate a random number between 0 and 1 to decide if we include the peptide or not in the training sample
rtPeptideUnique[, trainRand := runif(NROW(rtPeptideUnique))]

setkey(rtPeptideUnique, rtCent)
rtPeptideUnique[, train := FALSE]
rtPeptideUnique[early, train := (trainRand < (threshold * earlyCoefficient))]
rtPeptideUnique[middle, train := (trainRand < (threshold * middleCoefficient))]
rtPeptideUnique[late, train := (trainRand < (threshold * lateCoefficient))]

ggplot(rtPeptideUnique[train==TRUE], aes(rtCent)) + geom_histogram()

rtPeptideTrain <- rtPeptideUnique[train == TRUE, list(test, q50_2)]
min_RT <- rtPeptideTrain[, min(q50_2)]
max_RT <- rtPeptideTrain[, max(q50_2)]
rtPeptideTrain[, q50_2 := (q50_2-min_RT)/(max_RT-min_RT)]
# modified_sequence is included because it's the key, but I suspect this kind of behavior to change, so check.
# We want to assess the performance of the algorithm on the actual retention time observed, not the median.
rtPeptideTest <- rtPeptide[rtPeptideUnique[train == FALSE, test], rtsec]
rtPeptideTest[, rtsec := (rtsec-min_RT)/(max_RT-min_RT)]
rtPeptideTestMedian <- unique(rtPeptide[rtPeptideUnique[train == FALSE, test], q50_2])
rtPeptideTestMedian[, q50_2 := (q50_2-min_RT)/(max_RT-min_RT)]
# To use with RTModel and RTPredict (OpenMS)
rtPeptideTestRTPred <- rtPeptideTestMedian[, test]

# To execute if we want to slightly randomize the retention times
setkey(rtPeptideTrain, modified_sequence)
rtPeptideTrain <- unique(rtPeptideTrain)
rtPeptideTrain[, simulated_RT := rnorm(NROW(rtPeptideTrain), mean = q50_2, sd = (0.05*q50_2))]
rtPeptideTrain <- rtPeptideTrain[, list(modified_sequence,simulated_RT)]
# ggplot(rtPeptideTrain, aes(q50_2, simulated_RT)) + geom_point(alpha=(1/2))


write.table(rtPeptideTrain, file=paste0(projectPath, "/data/rtPeptideTrain.txt"), quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
write.table(rtPeptideTest, file=paste0(projectPath, "/data/rtPeptideTest.txt"), quote = FALSE, , sep="\t", row.names = FALSE, col.names = FALSE)
write.table(rtPeptideTestMedian, file=paste0(projectPath, "/data/rtPeptideTestMedian.txt"), quote = FALSE, , sep="\t", row.names = FALSE, col.names = FALSE)
write.table(rtPeptideTestRTPred, file=paste0(projectPath, "/data/rtPeptideTestRTPred.txt"), quote = FALSE, , sep="\t", row.names = FALSE, col.names = FALSE)

# Calling ELUDE from R using the command shell, does it work ?

projectPath <- "C:/Users/Nicolas Housset/Documents/RetentionTimeAnalysis"
trainData <- "/data/rtPeptideTrain.txt"
testData <- "/data/rtPeptideTestMedian.txt"
saveModel <- "/data/modelTest.model"
saveIndex <- "/data/retentionIndexTest.index"
savePredict <- "/data/predictions.out"

trainData <- shQuote(paste0(projectPath, trainData))
testData <- shQuote(paste0(projectPath, testData))
saveModel <- shQuote(paste0(projectPath, saveModel))
saveIndex <- shQuote(paste0(projectPath, saveIndex))
savePredict <- shQuote(paste0(projectPath, savePredict))

verbFlag <- " -v "
trainFlag <- " -t "
testFlag <- " -e "
saveModelFlag <- " -s "
saveIndexFlag <- " -r "
savePredictFlag <- " -o "
testRTFlag <- " -g "
# noInSourceFlag <- " -y "
ignoreNewTestPTMFlag <- " -p "
verbLevel <- " 5"

eludePath <- "C:/Program Files (x86)/Elude"
strCommand <- paste0("cd ",shQuote(eludePath), " && elude ", verbFlag, verbLevel, trainFlag, trainData, testFlag,
                     testData, saveModelFlag, saveModel, saveIndexFlag, saveIndex,
                     savePredictFlag, savePredict, testRTFlag, ignoreNewTestPTMFlag)

# elude -v 4 -t "C:\Users\Nicolas Housset\Documents\RetentionTimeAnalysis\data\rtPeptideTrain.txt" -e "C:\Users\Nicolas Housset\Documents\RetentionTimeAnalysis\data\rtPeptideTest.txt" -s "C:\Users\Nicolas Housset\Documents\RetentionTimeAnalysis\data\modelTest" -r "C:\Users\Nicolas Housset\Documents\RetentionTimeAnalysis\data\retentionIndexTest" -o "C:\Users\Nicolas Housset\Documents\RetentionTimeAnalysis\data\predictions.out" -g -p

# Alright, the trick is to change the working directory to ELUDE's
shell(strCommand, translate = TRUE, wait = TRUE)

results <- data.table(read.table(file=paste0(projectPath, "/data/predictions.out"), header = TRUE, sep = "\t"))

ggplot(results, aes(Predicted_RT, Observed_RT)) + geom_point(alpha=(1/2))
ggplot(results, aes(Observed_RT, Predicted_RT - Observed_RT)) + geom_point(alpha=(1/2))
ggplot(results, aes(Predicted_RT, Predicted_RT - Observed_RT)) + geom_point(alpha=(1/4))

ggplot(rtPeptideTrain, aes("A",q50_2)) + geom_jitter(alpha=(1/2))

testData <- "/data/rtPeptideTestMedian.txt"
testData <- shQuote(paste0(projectPath, testData))

savePredict <- "/data/predictionsTestMedian.out"
savePredict <- shQuote(paste0(projectPath, savePredict))

loadModelFlag <- " -l "

# This time we apply the model on the testing data, but using the median
eludePath <- "C:/Program Files (x86)/Elude"
strCommand <- paste0("cd ",shQuote(eludePath), " && elude ", verbFlag, verbLevel, testFlag,
                     testData, loadModelFlag, saveModel,
                     savePredictFlag, savePredict, testRTFlag, ignoreNewTestPTMFlag)
shell(strCommand, translate = TRUE, wait = TRUE)

resultsTest <- data.table(read.table(file=paste0(projectPath, "/data/predictionsTestMedian.out"), header = TRUE, sep = "\t"))

ggplot(resultsTest, aes(Predicted_RT, Observed_RT)) + geom_point(alpha=(1/2))
resultsTest[, diff := Predicted_RT - Observed_RT]
setkey(resultsTest, Observed_RT)

for(i in 499:(NROW(resultsTest))){
  resultsTest[i, moving_average := resultsTest[(i-499):i, mean(diff)]]
  resultsTest[i, tmoving_average := resultsTest[(i-499):i, mean(diff, trim = 0.25)]]
  resultsTest[i, movq95 := resultsTest[(i-499):i, quantile(diff, probs = 0.975) - quantile(diff, probs = 0.025)]]
  resultsTest[i, movq90 := resultsTest[(i-499):i, quantile(diff, probs = 0.95) - quantile(diff, probs = 0.05)]]
  resultsTest[i, movq50 := resultsTest[(499):i, quantile(diff, probs = 0.75) - quantile(diff, probs = 0.25)]]
}
ggplot(resultsTest, aes(Observed_RT, moving_average)) + geom_point(alpha=(1/2)) + xlim(750,2500)+ylim(-400,100)
ggplot(resultsTest, aes(Observed_RT, tmoving_average)) + geom_point(alpha=(1/2)) + xlim(750,2500)+ylim(-400,100)
ggplot(resultsTest, aes(Observed_RT, movq95)) + geom_point(alpha=(1/2))
ggplot(resultsTest, aes(Observed_RT, movq90)) + geom_point(alpha=(1/2))
ggplot(resultsTest, aes(Observed_RT, movq50)) + geom_point(alpha=(1/2))

ggplot(resultsTest, aes(Observed_RT, movq95)) + geom_smooth(alpha=(1/2)) + xlim(750,2500)+ylim(-400,100)



resultsTrain[, bsv := 0]
resultsTrain[roundedDiff == -9 | roundedDiff == 9, bsv := 1]
setkey(resultsTrain, Observed_RT)
resultsTrain_1 <- resultsTrain[bsv==0]
resultsTrain_2 <- resultsTrain[bsv==1]
for(i in 20:(NROW(resultsTrain_1))){
  resultsTrain_1[i, moving_average := resultsTrain_1[(i-19):i, mean(diff)]]
  resultsTrain_1[i, movq95 := resultsTrain_1[(i-19):i, quantile(diff, probs = 0.975) - quantile(diff, probs = 0.025)]]
  resultsTrain_1[i, movq90 := resultsTrain_1[(i-19):i, quantile(diff, probs = 0.95) - quantile(diff, probs = 0.05)]]
  resultsTrain_1[i, movq50 := resultsTrain_1[(i-19):i, quantile(diff, probs = 0.75) - quantile(diff, probs = 0.25)]]
}
ggplot(resultsTrain_1, aes(Observed_RT, moving_average)) + geom_point(alpha=(1/2)) + xlim(750,2500)+ylim(-400,100)
ggplot(resultsTrain_1, aes(Observed_RT, movq95)) + geom_point(alpha=(1/2))
ggplot(resultsTrain_1, aes(Observed_RT, movq90)) + geom_point(alpha=(1/2))
ggplot(resultsTrain_1, aes(Observed_RT, movq50)) + geom_point(alpha=(1/2))

resultsTrain_1[, summary(diff)]

table(resultsTrain[, bsv])

ggplot(resultsTrain, aes(x = Observed_RT, colour = factor(bsv))) + geom_density(alpha=(1/2))

ggplot(resultsTrain, aes(x = Observed_RT, y = diff, colour = factor(bsv))) + geom_point(alpha=(1/2))


setkey(resultsTrain, diff)

setkey(results, Observed_RT)
results[, diff := Predicted_RT - Observed_RT]

# View by centile. Quick, but moving average is way better
results[, centile := 1:NROW(results)]
results[, centile := ceiling(100 * centile / NROW(results))]
results[, meanError := mean(diff), by = centile]
results[, q95_pred := quantile(diff, probs = 0.975) - quantile(diff, probs = 0.025), by = centile]
results[, q90_pred := quantile(diff, probs = 0.95) - quantile(diff, probs = 0.05), by = centile]
ggplot(results, aes(centile, meanError)) + geom_point(alpha=(1/2))
ggplot(results, aes(centile, q95_pred)) + geom_point(alpha=(1/2))
ggplot(results, aes(centile, q90_pred)) + geom_point(alpha=(1/2))

for(i in 1:(NROW(results)-999)){
  results[i, moving_average := results[i:(i+999), mean(diff)]]
  results[i, movq95 := results[i:(i+999), quantile(diff, probs = 0.975) - quantile(diff, probs = 0.025)]]
}

for(i in 500:(NROW(results))){
  results[i, movq975 := results[(i-499):i, quantile(diff, probs = 0.975)]]
  results[i, movq025 := results[(i-499):i, quantile(diff, probs = 0.025)]]
}
ggplot(results, aes(Observed_RT, moving_average)) + geom_point(alpha=(1/2)) + xlim(750,2500)+ylim(-400,100)
ggplot(results, aes(Observed_RT, movq95)) + geom_point(alpha=(1/2))
ggplot(results, aes(Observed_RT, movq90)) + geom_point(alpha=(1/2))
ggplot(results, aes(Observed_RT, movq50)) + geom_point(alpha=(1/2))


ggplot(results, aes(Observed_RT, movq975)) + geom_point(alpha=(1/2))
ggplot(results, aes(Observed_RT, movq025)) + geom_point(alpha=(1/2))
ggplot(results, aes(Observed_RT, list(movq975,movq025))) + geom_point(alpha=(1/2))

results[, movq95 := movq975 - movq025]

results[, part := 1]
graphDS <- results[, list(Observed_RT, Predicted_RT, movq95, part)][0]
listeVar <- list("movq95", "movq975", "movq025")
for(i in 1:3){
  eval(parse(text=paste0("graphPart <- results[, list(Observed_RT, Predicted_RT, ", listeVar[i], ")]")))
  eval(parse(text=paste0("graphPart[,", listeVar[[1]], " := ", listeVar[[i]], "]")))
  eval(parse(text=paste0("graphPart[, part := ", i, "]")))
  eval(parse(text=paste0("graphDS <- rbind(graphDS, graphPart[, list(Observed_RT, Predicted_RT, movq95, part)])"))) 
}
graphDS <- results[, list(evalsubstitute('Observed_RT'), Predicted_RT, movq975)]
graphDS[, movq025 := movq975]
graphDS[, movq975 := NULL]
graphDS[, part := 1]
graphDS <- rbind(graphDS, results[, list(Observed_RT, Predicted_RT, movq025, part)])

ggplot(graphDS, aes(Observed_RT, movq95, colour = factor(part, labels = c("95 % width", "0.975 quantile", "0.025 quantile")))) + geom_point(alpha=(1/2)) + xlim(750,2500)+ ylim(-625,625)
