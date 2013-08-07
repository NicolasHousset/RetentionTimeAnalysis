load("C:/Users/Nicolas Housset/Documents/rt_project.RData")
library(data.table)

library(ggplot2)

library(RMySQL);
con <- dbConnect(MySQL(), user="nicolas", password="nicolas,13*", dbname="projects", host="muppet03.ugent.be");

# Selecting relevant information from the database
# Selection criteria are project numbers (to avoid huge database) and orbitrap instruments (id between 8 and 10)

rt_project1 <- data.table(dbGetQuery(con,
                                    "SELECT scanid, number, spectrumid, l_lcrunid, l_projectid, l_instrumentid, identified, score, identitythreshold, confidence, DB,
rtsec, total_spectrum_intensity, mass_to_charge, spectrum.charge, accession, start, end, sequence, modified_sequence FROM
(spectrum LEFT JOIN scan ON spectrum.spectrumid = scan.l_spectrumid 
LEFT JOIN identification ON spectrum.spectrumid = identification.l_spectrumid
LEFT JOIN spectrum_file ON spectrum.spectrumid = spectrum_file.l_spectrumid)
WHERE l_projectid BETWEEN 651 AND 750 AND l_instrumentid BETWEEN 8 AND 14;"));

save(rt_project1, file = "C:/Users/Nicolas Housset/Documents/rt_project1.RData", compression_level = 1)
rm(rt_project1)
gc()

rt_project2 <- data.table(dbGetQuery(con,
                                     "SELECT scanid, number, spectrumid, l_lcrunid, l_projectid, l_instrumentid, identified, score, identitythreshold, confidence, DB,
rtsec, total_spectrum_intensity, mass_to_charge, spectrum.charge, accession, start, end, sequence, modified_sequence FROM
(spectrum LEFT JOIN scan ON spectrum.spectrumid = scan.l_spectrumid 
LEFT JOIN identification ON spectrum.spectrumid = identification.l_spectrumid
LEFT JOIN spectrum_file ON spectrum.spectrumid = spectrum_file.l_spectrumid)
WHERE l_projectid BETWEEN 751 AND 850 AND l_instrumentid BETWEEN 8 AND 14;"));

save(rt_project2, file = "C:/Users/Nicolas Housset/Documents/rt_project2.RData", compression_level = 1)
rm(rt_project2)
gc()

rt_project3 <- data.table(dbGetQuery(con,
                                     "SELECT scanid, number, spectrumid, l_lcrunid, l_projectid, l_instrumentid, identified, score, identitythreshold, confidence, DB,
rtsec, total_spectrum_intensity, mass_to_charge, spectrum.charge, accession, start, end, sequence, modified_sequence FROM
(spectrum LEFT JOIN scan ON spectrum.spectrumid = scan.l_spectrumid 
LEFT JOIN identification ON spectrum.spectrumid = identification.l_spectrumid
LEFT JOIN spectrum_file ON spectrum.spectrumid = spectrum_file.l_spectrumid)
WHERE l_projectid BETWEEN 851 AND 950 AND l_instrumentid BETWEEN 8 AND 14;"));

save(rt_project3, file = "C:/Users/Nicolas Housset/Documents/rt_project3.RData", compression_level = 1)
rm(rt_project3)
gc()

rt_project4 <- data.table(dbGetQuery(con,
                                     "SELECT scanid, number, spectrumid, l_lcrunid, l_projectid, l_instrumentid, identified, score, identitythreshold, confidence, DB,
rtsec, total_spectrum_intensity, mass_to_charge, spectrum.charge, accession, start, end, sequence, modified_sequence FROM
(spectrum LEFT JOIN scan ON spectrum.spectrumid = scan.l_spectrumid 
LEFT JOIN identification ON spectrum.spectrumid = identification.l_spectrumid
LEFT JOIN spectrum_file ON spectrum.spectrumid = spectrum_file.l_spectrumid)
WHERE l_projectid BETWEEN 951 AND 1050 AND l_instrumentid BETWEEN 8 AND 14;"));

save(rt_project4, file = "C:/Users/Nicolas Housset/Documents/rt_project4.RData", compression_level = 1)
rm(rt_project4)
gc()

instruments <- data.table(dbGetQuery(con,
                                   "SELECT * FROM instrument"));
save(instruments, file="C:/Users/Nicolas Housset/Documents/instruments.RData")

protocols <- data.table(dbGetQuery(con,
                                   "SELECT * FROM protocol"));
save(protocols, file="C:/Users/Nicolas Housset/Documents/protocols.RData")
# Instruments 10, 13 and 14 are probably of better quality

projects <- data.table(dbGetQuery(con,
                                  "SELECT * from project"))
save(projects, file="C:/Users/Nicolas Housset/Documents/projects.RData")
write.csv(rt_project, file="C:/Users/Nicolas Housset/Documents/rt_project.csv")

# Working on identified peptides
load(file = "C:/Users/Nicolas Housset/Documents/rt_project3.RData")
identified <- rt_project3[(!is.na(modified_sequence))]

# Often the best solution is the easiest one. Here we go, with the unique function
setkey(identified, l_lcrunid, modified_sequence)
counts <- unique(identified)



counts[, modified_sequence.f := factor(modified_sequence)]

liste_peptides <- summary(counts[, modified_sequence.f], maxsum = 2000)


# Adding the protocol information to our data.
setkey(identified, l_projectid)
projects[, projectid := as.numeric(projectid)]
load(file="C:/Users/Nicolas Housset/Documents/projects.RData")
setkey(projects, projectid)

identified <- identified[projects, j=l_protocolid][identified]


setkey(identified,modified_sequence)

one_peptide <- identified[labels(liste_peptides)[1:20]]

load(file="C:/Users/Nicolas Housset/Documents/protocols.RData")

# Keeping the most intense spectrum intensity information
setkey(one_peptide,l_lcrunid,modified_sequence,total_spectrum_intensity)
one_peptide_reduced <- one_peptide[one_peptide[,j=list(.I[which.max(total_spectrum_intensity)]), by=c("l_lcrunid","modified_sequence")][,j=V1]]

one_peptide_reduced <- one_peptide_reduced[,nombre := .N, by = list(l_projectid, modified_sequence)][nombre >= 3]
# one_peptide_reduced <- one_peptide_reduced[l_instrumentid != 10]

plot_global <- ggplot(one_peptide_reduced, aes(modified_sequence, rtsec, colour=factor(l_projectid)))

plot_global + geom_jitter(alpha = 1) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + facet_grid(. ~ sequence)

plot_global + geom_jitter(alpha = 1/2, size = 1.5) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + facet_grid(. ~ l_instrumentid)

plot_global + geom_jitter(alpha = 1/2, size = 1.5) + ylim(0,3000) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Making the plots of retention times with fewer projects
one_peptide_reduced[, l_projectid.f := factor(l_projectid)]
liste_projects <- summary(one_peptide_reduced[, l_projectid.f])
setkey(one_peptide_reduced, l_projectid.f)

plot_global <- ggplot(one_peptide_reduced[labels(liste_projects)[1:10]], aes(modified_sequence, rtsec, colour=factor(l_projectid)))
plot_global + geom_jitter(alpha = 1/2, size = 1.5) + ylim(0,3000) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

one_peptide_reduced[, grouping := paste0(l_projectid,modified_sequence)]
# Binning the rt in a quick way
one_peptide_reduced[, rtint := floor(rtsec)]
one_peptide_reduced[, rttri := floor(rtsec/3)]
one_peptide_reduced[, rtdec := floor(rtsec/10)]
# Adding the shannon entropy to the data (needs the entropy package), with different bin sizes of retention times
one_peptide_reduced[, ShanEntropyInt := entropy(table(rtint),method="ML", unit="log2"), keyby = list(modified_sequence, l_projectid)]
one_peptide_reduced[, ShanEntropyTri := entropy(table(rttri),method="ML", unit="log2"), by = list(modified_sequence, l_projectid)]
one_peptide_reduced[, ShanEntropyDec := entropy(table(rtdec),method="ML", unit="log2"), by = list(modified_sequence, l_projectid)]

pep_stat <- one_peptide_reduced[,nombre := .N, by = list(l_projectid, modified_sequence)][nombre >= 10]

pepstat <- unique(setkey(one_peptide_reduced[, j=list(quant25 = quantile(rtsec, probs=0.25),
                                                      quant50 = quantile(rtsec, probs=0.5),
                                                      quant75 = quantile(rtsec, probs=0.75),
                                                      stdev = sd(rtsec),
                                                      MS2 = mean(total_spectrum_intensity),
                                                      l_instrumentid,
                                                      ShanEntropyInt,
                                                      ShanEntropyTri,
                                                      ShanEntropyDec), by = list(l_projectid, modified_sequence)],
                         l_projectid, modified_sequence))

pepstat[, QCD := (quant75-quant25)/(quant75+quant25)]

# Plotting log(QCD) and shannon entropy, coloured by projectid
plot_stat <- ggplot(pepstat, aes(log(QCD), ShanEntropyTri, colour=factor(l_projectid)))
plot_stat <- ggplot(pepstat, aes(log(MS2), log(QCD), colour=factor(l_projectid)))
plot_stat <- ggplot(pepstat, aes(log(MS2), ShanEntropyTri, colour=factor(l_projectid)))
plot_stat + geom_point(alpha=1/2)

# Limiting the plot on the most common projects
pepstat[, l_projectid.f := factor(l_projectid)]
liste_projects <- summary(pepstat[, l_projectid.f])
setkey(pepstat, l_projectid.f)
plot_stat <- ggplot(pepstat[(labels(liste_projects)[1:10])], aes(log(QCD), ShanEntropyTri, colour=factor(l_projectid)))
plot_stat <- ggplot(pepstat[(labels(liste_projects)[1:10])], aes(log(MS2), log(QCD), colour=factor(l_projectid)))
plot_stat <- ggplot(pepstat[(labels(liste_projects)[1:10])], aes(log(MS2), ShanEntropyTri, colour=factor(l_projectid)))

plot_stat + geom_point(alpha=1/2)


# Plotting relationships between log(QCD), log(MS2) and shannon entropy, coloured by instrument
plot_stat <- ggplot(pepstat, aes(log(QCD), ShanEntropyTri, colour=factor(l_instrumentid)))
plot_stat <- ggplot(pepstat, aes(log(MS2), log(QCD), colour=factor(l_instrumentid)))
plot_stat <- ggplot(pepstat, aes(log(MS2), ShanEntropyTri, colour=factor(l_instrumentid)))

plot_stat + geom_point(alpha=1/2)


# Limiting the plot on the most common projects
pepstat[, l_projectid.f := factor(l_projectid)]
liste_projects <- summary(pepstat[, l_projectid.f])
setkey(pepstat, l_projectid.f)
plot_stat <- ggplot(pepstat[(labels(liste_projects)[1:7])], aes(log(QCD), ShanEntropyDec, colour=factor(l_instrumentid)))
plot_stat <- ggplot(pepstat[(labels(liste_projects)[1:7])], aes(log(MS2), log(QCD), colour=factor(l_instrumentid)))
plot_stat <- ggplot(pepstat[(labels(liste_projects)[1:7])], aes(log(MS2), ShanEntropyInt, colour=factor(l_instrumentid)))

plot_stat + geom_point(alpha=1/2)


# Question : are the projects different within an instrument ?
pepstat8 <- pepstat[l_instrumentid==8]

# Plotting log(QCD) and shannon entropy, coloured by projectid
plot_stat <- ggplot(pepstat8, aes(log(QCD), ShanEntropyTri, colour=factor(l_projectid)))
plot_stat <- ggplot(pepstat8, aes(log(MS2), log(QCD), colour=factor(l_projectid)))
plot_stat <- ggplot(pepstat8, aes(log(MS2), ShanEntropyTri, colour=factor(l_projectid)))
plot_stat + geom_point(alpha=1/2)

# For the instrument 8, projects are still separated on those graphs


pepstat9 <- pepstat[l_instrumentid==9]

# Plotting log(QCD) and shannon entropy, coloured by projectid
plot_stat <- ggplot(pepstat9, aes(log(QCD), ShanEntropyTri, colour=factor(l_projectid)))
plot_stat <- ggplot(pepstat9, aes(log(MS2), log(QCD), colour=factor(l_projectid)))
plot_stat <- ggplot(pepstat9, aes(log(MS2), ShanEntropyTri, colour=factor(l_projectid)))
plot_stat + geom_point(alpha=1/2)

liste_projects9 <- summary(pepstat9[, l_projectid.f])
setkey(pepstat9, l_projectid.f)
plot_stat <- ggplot(pepstat9[(labels(liste_projects9)[1:5])], aes(log(QCD), ShanEntropyTri, colour=factor(l_projectid)))
plot_stat <- ggplot(pepstat9[(labels(liste_projects9)[1:5])], aes(log(MS2), log(QCD), colour=factor(l_projectid)))
plot_stat <- ggplot(pepstat9[(labels(liste_projects9)[1:5])], aes(log(MS2), ShanEntropyTri, colour=factor(l_projectid)))

# The separation is not so clear between projects within the instrument 9


pepstat10 <- pepstat[l_instrumentid==10]

# Plotting log(QCD) and shannon entropy, coloured by projectid
plot_stat <- ggplot(pepstat10, aes(log(QCD), ShanEntropyTri, colour=factor(l_projectid)))
plot_stat <- ggplot(pepstat10, aes(log(MS2), log(QCD), colour=factor(l_projectid)))
plot_stat <- ggplot(pepstat10, aes(log(MS2), ShanEntropyTri, colour=factor(l_projectid)))
plot_stat + geom_point(alpha=1/2)

liste_projects10 <- summary(pepstat10[, l_projectid.f])
setkey(pepstat10, l_projectid.f)
plot_stat <- ggplot(pepstat10[(labels(liste_projects10)[1:10])], aes(log(QCD), ShanEntropyTri, colour=factor(l_projectid)))
plot_stat <- ggplot(pepstat10[(labels(liste_projects10)[1:10])], aes(log(MS2), log(QCD), colour=factor(l_projectid)))
plot_stat <- ggplot(pepstat10[(labels(liste_projects10)[1:10])], aes(log(MS2), ShanEntropyTri, colour=factor(l_projectid)))

# There is good separation within the instrument 10


pepstat13 <- pepstat[l_instrumentid==13]

# Plotting log(QCD) and shannon entropy, coloured by projectid
plot_stat <- ggplot(pepstat13, aes(log(QCD), ShanEntropyTri, colour=factor(l_projectid)))
plot_stat <- ggplot(pepstat13, aes(log(MS2), log(QCD), colour=factor(l_projectid)))
plot_stat <- ggplot(pepstat13, aes(log(MS2), ShanEntropyTri, colour=factor(l_projectid)))
plot_stat + geom_point(alpha=1/2)

liste_projects13 <- summary(pepstat10[, l_projectid.f])
setkey(pepstat13, l_projectid.f)
plot_stat <- ggplot(pepstat10[(labels(liste_projects13)[1:10])], aes(log(QCD), ShanEntropyTri, colour=factor(l_projectid)))
plot_stat <- ggplot(pepstat10[(labels(liste_projects13)[1:10])], aes(log(MS2), log(QCD), colour=factor(l_projectid)))
plot_stat <- ggplot(pepstat10[(labels(liste_projects13)[1:10])], aes(log(MS2), ShanEntropyTri, colour=factor(l_projectid)))

# Not enough data for instrument 13


pepstat14 <- pepstat[l_instrumentid==14]

# Plotting log(QCD) and shannon entropy, coloured by projectid
plot_stat <- ggplot(pepstat14, aes(log(QCD), ShanEntropyTri, colour=factor(l_projectid)))
plot_stat <- ggplot(pepstat14, aes(log(MS2), log(QCD), colour=factor(l_projectid)))
plot_stat <- ggplot(pepstat14, aes(log(MS2), ShanEntropyTri, colour=factor(l_projectid)))
plot_stat + geom_point(alpha=1/2)

liste_projects14 <- summary(pepstat14[, l_projectid.f])
setkey(pepstat14, l_projectid.f)
plot_stat <- ggplot(pepstat10[(labels(liste_projects14)[1:10])], aes(log(QCD), ShanEntropyTri, colour=factor(l_projectid)))
plot_stat <- ggplot(pepstat10[(labels(liste_projects14)[1:10])], aes(log(MS2), log(QCD), colour=factor(l_projectid)))
plot_stat <- ggplot(pepstat10[(labels(liste_projects14)[1:10])], aes(log(MS2), ShanEntropyTri, colour=factor(l_projectid)))

# Good separation within instrument 14
