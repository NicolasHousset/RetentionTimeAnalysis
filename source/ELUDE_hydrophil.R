projectPath <- "C:/Users/Nicolas Housset/Documents/RetentionTimeAnalysis"
testData <- "/data/DeNovoDigest.txt"
savePredict <- "/data/predictions.out"
loadmodel <- "/data/modelHydrophil.model"
testData <- shQuote(paste0(projectPath, testData))
savePredict <- shQuote(paste0(projectPath, savePredict))

verbFlag <- " -v "
testFlag <- " -e "
savePredictFlag <- " -o "
ignoreNewTestPTMFlag <- " -p "
verbLevel <- " 5"
loadModelFlag <- " -l "

# This time we apply the model on the testing data, but using the median
eludePath <- "C:/Program Files (x86)/Elude"
strCommand <- paste0("cd ",shQuote(eludePath), " && elude ", verbFlag, verbLevel, testFlag,
                     testData, loadModelFlag, saveModel,
                     savePredictFlag, savePredict, testRTFlag, ignoreNewTestPTMFlag)
shell(strCommand, translate = TRUE, wait = TRUE)