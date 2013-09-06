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

verbLevel <- " 5"

eludePath <- "C:/Program Files (x86)/Elude"
strCommand <- paste0("cd ",shQuote(eludePath), " && elude ", verbFlag, verbLevel, trainFlag, trainData, testFlag,
                     testData, saveModelFlag, saveModel, saveIndexFlag, saveIndex,
                     savePredictFlag, savePredict, testRTFlag)



# elude -v 4 -t "C:\Users\Nicolas Housset\Documents\RetentionTimeAnalysis\data\rtPeptideTrain.txt" -e "C:\Users\Nicolas Housset\Documents\RetentionTimeAnalysis\data\rtPeptideTest.txt" -s "C:\Users\Nicolas Housset\Documents\RetentionTimeAnalysis\data\modelTest" -r "C:\Users\Nicolas Housset\Documents\RetentionTimeAnalysis\data\retentionIndexTest" -o "C:\Users\Nicolas Housset\Documents\RetentionTimeAnalysis\data\predictions.out" -g -p

# Alright, the trick is to change the working directory to ELUDE's
shell(strCommand, translate = TRUE, wait = TRUE)
