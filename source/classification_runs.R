# Starting project : 2065
# Ending project : ???

# We restrict ourselves to Vanessa
# Classification in Low/High/Very high retention times
# The latter one is to put apart some weird experiments

projectPath <- "C:/Users/Nicolas Housset/Documents/RetentionTimeAnalysis"


load(file = paste0(projectPath,"/data/identified_protocol.RData"))

setkey(identified_subs, l_projectid)
test <- identified_subs[l_projectid >= "2065" & l_instrumentid == "10" & l_protocolid == "11"][, list(l_projectid, l_lcrunid, sequence, modified_sequence, index_rt2, l_protocolid, grpProj, rtsec, q50_2)]

setkey(test, l_projectid, l_lcrunid)
list_runs <- unique(test)[,l_projectid, l_lcrunid]

setkey(list_runs, l_projectid)


list_runs[as.character(2065:2079), group_run := 1L]
list_runs[as.character(2104), group_run := 2L]
list_runs[as.character(2127:2142), group_run := 3L]
list_runs[as.character(2150:2151), group_run := 4L]
list_runs[as.character(2154:2157), group_run := 5L]
list_runs[as.character(2169:2172), group_run := 5L] # repeat of 2154:2157
list_runs[as.character(2181:2189), group_run := 6L]
list_runs[as.character(2196:2197), group_run := 7L] # quite the unusual one, probably impossible to classify
list_runs[as.character(2210:2212), group_run := 8L]
list_runs[as.character(2229:2243), group_run := 9L]
list_runs[as.character(2247:2252), group_run := 10L] # Guess : only project 2250 (run 95548) is on the higher retention time
list_runs[as.character(2253), group_run := 11L] # No contaminant on this project ? Pattern to confirm : odd high, even low
list_runs[as.character(2263), group_run := 12L] # No contaminant on this project ? Pattern to confirm : odd high, even low
