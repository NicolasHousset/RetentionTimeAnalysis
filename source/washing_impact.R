# Suspicion that the two columns don't have the same gradient of B
# First column reaches 55 % of B at 30 minutes (1800 sec)
# Secund columns reaches 55 % of B at 33 minutes (2100 sec)

# We should see a shift in the retention times
# Typical example : project 2150, with 24 runs : from 93683 to 93706
# Odd runs have higher rt
# Even runs have lower rt
# A (quick) conclusion would be that even runs match to the first column, and odd runs to the secund column

# But I wonder if the rise in gradient even start at the same point ?



projectPath <- "C:/Users/Nicolas Housset/Documents/RetentionTimeAnalysis"


load(file = paste0(projectPath,"/data/identified_protocol.RData"))

setkey(identified_subs, l_projectid)
test <- identified_subs[as.character(2127:2142)][, list(l_projectid, l_lcrunid, sequence, modified_sequence, index_rt2, l_protocolid, grpProj, rtsec, q50_2)]
setkey(test,l_projectid)
ggplot(test[as.character(2127)], aes(rtsec)) + geom_histogram(binwidth=15) + xlim(1800,2600)

gradient_dt <- data.table(1:99)
for(i in 2127:2142){
  eval(parse(text=paste0("gradient_dt[, run_",i," := quantile(test[as.character(i)][,rtsec], probs = seq(0.01,0.99,0.01))]")))
}
list_var <- vector("list", 16)
for(i in 2127:2142){
  list_var[[i-2126]] <- paste0("run_",i)
}

gradient_dt2 <- data.table(melt(gradient_dt, id.vars = "V1", measure.vars = as.character(list_var), variable.name = "run_id", value.name = "RT"))
setkey(gradient_dt2, run_id)
gradient_dt2[as.character(list_var[c(1,3,5,7,9,11,13,15)]), group_run := "odd"]
gradient_dt2[as.character(list_var[c(2,4,6,8,10,12,14,16)]), group_run := "even"]

ggplot(gradient_dt2[as.character(list_var[c(1,3,5,7,9,11,13,15)])], aes(V1,RT, colour = run_id)) + geom_line() + xlim(30,70) + ylim(900,1600)
ggplot(gradient_dt2[as.character(list_var[c(2,4,6,8,10,12,14,16)])], aes(V1,RT, colour = run_id)) + geom_line() + xlim(30,70) + ylim(1200,1800)

setkey(gradient_dt2, group_run)
ggplot(gradient_dt2, aes(V1,RT, colour = group_run, shape = group_run)) + geom_point(size = 1.5) + xlim(1,99) + ylim(800,2200)

gradient_dt2["even", RT2 := RT - 90]
gradient_dt2["odd", RT2 := RT]
ggplot(gradient_dt2, aes(V1,RT2, colour = group_run, shape = group_run)) + geom_point(size = 1.5) + xlim(1,99) + ylim(800,2200)

setkey(gradient_dt2, V1, group_run)
a <- unique(gradient_dt2[, med_grp := quantile(RT, probs = 0.5), by = list(V1, group_run)])
setkey(a, group_run)
b <- a["odd"]
c <- a["even"]

d <- data.frame(b[, med_grp] - c[, med_grp])
test[,quant_0.05 := quantile(rtsec, probs = 0.05),by = l_lcrunid]


quantile(test[as.character(93683)][,rtsec], probs = seq(0.01,0.99,0.01))
test[,quantile(rtsec, probs = seq(0.01,0.99,0.01),by = l_lcrunid]

test[,quantile(rtsec, probs = 0.50),by = l_lcrunid]
