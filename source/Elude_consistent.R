# Trying to see what happens if we take more consistent data

library(data.table)
library(ggplot2)

projectPath <- "C:/Users/Nicolas Housset/Documents/RetentionTimeAnalysis"


load(file = paste0(projectPath,"/data/identified_protocol.RData"))
setkey(identified_subs, modified_sequence)

rtPeptide <- identified_subs[nbProjPepProtocol11 > 4][, list(l_projectid, l_lcrunid, modified_sequence, index_rt2, l_protocolid, grpProj, rtsec, q50_2)]
rtPeptide <- rtPeptide[l_protocolid=='11'][grpProj==2]

# This line is useful if selecting customized project numbers, recomputes the median
# rtPeptide[index_rt2 <2, q50_3 := quantile(rtsec, probs = 0.50), by = c("modified_sequence")]

# For now, only keeping the non-modified peptides
rtPeptide[, matching := grepl("NH2-[ABCDEFGHIJKLMNOPQRSTUVWXYZ]*-COOH", modified_sequence)]
rtPeptide <- rtPeptide[matching==TRUE]
# setkey(rtPeptide, modified_sequence)
# rtPeptide <- unique(rtPeptide)
rtPeptide <- rtPeptide[index_rt2 <2] 
rtPeptide[, modified_sequence := sub("NH2-", "", modified_sequence)]
rtPeptide[, modified_sequence := sub("-COOH", "", modified_sequence)]
rtPeptide <- rtPeptide[, list(modified_sequence, q50_2, rtsec)]

# rtPeptide[, q50_2 := q50_2 / 60]
# summary(rtPeptide, q50_2)
setkey(rtPeptide, modified_sequence)
rtPeptideUnique <- unique(rtPeptide)

# Taking the peptides in a different order to randomize the training set
rtPeptideFlush <- sample(rtPeptideUnique[, modified_sequence])
setkey(rtPeptideUnique, modified_sequence)
rtPeptideUnique <- rtPeptideUnique[rtPeptideFlush]

# Stratified sampling : we want to have  more peptides selected at early and late retention times, less in the middle
early <- as.character(7:10)
late <- as.character(18:24)
middle <- as.character(11:17)
earlyCoefficient <- 2
lateCoefficient <- 4
middleCoefficient <- 0.4
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

rtPeptideTrain <- rtPeptideUnique[train == TRUE, list(modified_sequence, q50_2)]
# modified_sequence is included because it's the key, but I suspect this kind of behavior to change, so check.
# We want to assess the performance of the algorithm on the actual retention time observed, not the median.
rtPeptideTest <- rtPeptide[rtPeptideUnique[train == FALSE, modified_sequence], rtsec]
rtPeptideTestMedian <- unique(rtPeptide[rtPeptideUnique[train == FALSE, modified_sequence], q50_2])


setkey(rtPeptideTrain, modified_sequence)
rtPeptideTrain <- unique(rtPeptideTrain)
rtPeptideTrain[, simulated_RT := rnorm(NROW(rtPeptideTrain), mean = q50_2, sd = (0.05*q50_2))]
rtPeptideTrain <- rtPeptideTrain[, list(modified_sequence,simulated_RT)]
# ggplot(rtPeptideTrain, aes(q50_2, simulated_RT)) + geom_point(alpha=(1/2))

setkey(rtPeptideTrain, modified_sequence)
uniqueTrain <- unique(rtPeptideTrain)
write.table(rtPeptideTrain, file=paste0(projectPath, "/data/rtPeptideTrain.txt"), quote = FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
write.table(rtPeptideTest, file=paste0(projectPath, "/data/rtPeptideTest.txt"), quote = FALSE, , sep="\t", row.names = FALSE, col.names = FALSE)
write.table(rtPeptideTestMedian, file=paste0(projectPath, "/data/rtPeptideTestMedian.txt"), quote = FALSE, , sep="\t", row.names = FALSE, col.names = FALSE)

# Calling ELUDE from R using the command shell, does it work ?

projectPath <- "C:/Users/Nicolas Housset/Documents/RetentionTimeAnalysis"
trainData <- "/data/rtPeptideTrain.txt"
testData <- "/data/rtPeptideTest.txt"
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
verbLevel <- " 5"

eludePath <- "C:/Program Files (x86)/Elude"
strCommand <- paste0("cd ",shQuote(eludePath), " && elude ", verbFlag, verbLevel, trainFlag, trainData, testFlag,
                     testData, saveModelFlag, saveModel, saveIndexFlag, saveIndex,
                     savePredictFlag, savePredict, testRTFlag)

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
                     savePredictFlag, savePredict, testRTFlag, noInSourceFlag)
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

ggplot(results, aes(Observed_RT, moving_average)) + geom_point(alpha=(1/2)) + xlim(750,2500)+ylim(-400,100)
ggplot(results, aes(Observed_RT, movq95)) + geom_point(alpha=(1/2))
ggplot(results, aes(Observed_RT, movq90)) + geom_point(alpha=(1/2))
ggplot(results, aes(Observed_RT, movq50)) + geom_point(alpha=(1/2))

