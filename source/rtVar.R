# Binning the rt in a quick way
one_peptide[, rtint := floor(rtsec)]
one_peptide[, rttri := floor(rtsec/3)]
one_peptide[, rtdec := floor(rtsec/10)]
# Adding the shannon entropy to the data (needs the entropy package), with different bin sizes of retention times
one_peptide[, ShanEntropyInt := entropy(table(rtint),method="ML", unit="log2"), keyby = list(modified_sequence, l_projectid)]
one_peptide[, ShanEntropyTri := entropy(table(rttri),method="ML", unit="log2"), by = list(modified_sequence, l_projectid)]
one_peptide[, ShanEntropyDec := entropy(table(rtdec),method="ML", unit="log2"), by = list(modified_sequence, l_projectid)]

pepstat <- unique(setkey(one_peptide[, j=list(quant25 = quantile(rtsec, probs=0.25),
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

one_peptide[, l_projectid.f := factor(l_projectid)]
liste_projects <- summary(one_peptide[, l_projectid.f], maxsum = 500)
setkey(pepstat, l_projectid)
plot_stat <- ggplot(pepstat[(labels(liste_projects)[100:500])], aes(log(QCD), ShanEntropyTri, colour=factor(l_instrumentid)))
plot_stat <- ggplot(pepstat[(labels(liste_projects)[100:500])], aes(log(MS2), log(QCD), colour=factor(l_projectid)))
plot_stat <- ggplot(pepstat[(labels(liste_projects)[50:500])], aes(log(MS2), ShanEntropyTri, colour=factor(l_projectid)))

plot_stat <- ggplot(pepstat, aes(modified_sequence, ShanEntropyDec, colour=factor(l_projectid)))
plot_stat + geom_jitter(alpha = 1/2, size = 1.5) + theme(axis.text.x = element_text(angle = 45, hjust = 1))


