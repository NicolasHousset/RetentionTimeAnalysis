# The purpose of this script is to extract relevant information concerning projects with id between 2000 and 3000.
# Reason behind this range of numbers is that projects are recent
# Projects are extracted by chunks of 100 to avoid overloading the RAM
# A db of identified peptides is built step by step

library(RMySQL);
library(data.table);

con <- dbConnect(MySQL(), group="MSDB", dbname="projects");

rt_project2m0 <- data.table(dbGetQuery(con,
                                     "SELECT scanid, number, spectrumid, l_lcrunid, l_projectid, l_instrumentid, identified, score, identitythreshold, confidence, DB,
rtsec, total_spectrum_intensity, mass_to_charge, spectrum.charge, accession, start, end, sequence, modified_sequence FROM
(spectrum LEFT JOIN scan ON spectrum.spectrumid = scan.l_spectrumid 
LEFT JOIN identification ON spectrum.spectrumid = identification.l_spectrumid
LEFT JOIN spectrum_file ON spectrum.spectrumid = spectrum_file.l_spectrumid)
WHERE l_projectid BETWEEN 2000 AND 2100 AND l_instrumentid BETWEEN 8 AND 14;"));

identified <- rt_project2m0[(!is.na(modified_sequence))]

save(rt_project2m0, file = "C:/Users/Nicolas Housset/Documents/rt_project2m0.RData", compression_level = 1)
rm(rt_project2m0)
gc()

rt_project2m1 <- data.table(dbGetQuery(con,
                                       "SELECT scanid, number, spectrumid, l_lcrunid, l_projectid, l_instrumentid, identified, score, identitythreshold, confidence, DB,
rtsec, total_spectrum_intensity, mass_to_charge, spectrum.charge, accession, start, end, sequence, modified_sequence FROM
(spectrum LEFT JOIN scan ON spectrum.spectrumid = scan.l_spectrumid 
LEFT JOIN identification ON spectrum.spectrumid = identification.l_spectrumid
LEFT JOIN spectrum_file ON spectrum.spectrumid = spectrum_file.l_spectrumid)
WHERE l_projectid BETWEEN 2101 AND 2200 AND l_instrumentid BETWEEN 8 AND 14;"));

identified <- rbind(identified,rt_project2m1[(!is.na(modified_sequence))])

save(rt_project2m1, file = "C:/Users/Nicolas Housset/Documents/rt_project2m1.RData", compression_level = 1)
rm(rt_project2m1)
gc()

rt_project2m2 <- data.table(dbGetQuery(con,
                                       "SELECT scanid, number, spectrumid, l_lcrunid, l_projectid, l_instrumentid, identified, score, identitythreshold, confidence, DB,
rtsec, total_spectrum_intensity, mass_to_charge, spectrum.charge, accession, start, end, sequence, modified_sequence FROM
(spectrum LEFT JOIN scan ON spectrum.spectrumid = scan.l_spectrumid 
LEFT JOIN identification ON spectrum.spectrumid = identification.l_spectrumid
LEFT JOIN spectrum_file ON spectrum.spectrumid = spectrum_file.l_spectrumid)
WHERE l_projectid BETWEEN 2201 AND 2300 AND l_instrumentid BETWEEN 8 AND 14;"));

identified <- rbind(identified,rt_project2m2[(!is.na(modified_sequence))])

save(rt_project2m2, file = "C:/Users/Nicolas Housset/Documents/rt_project2m2.RData", compression_level = 1)
rm(rt_project2m2)
gc()


rt_project2m3 <- data.table(dbGetQuery(con,
                                       "SELECT scanid, number, spectrumid, l_lcrunid, l_projectid, l_instrumentid, identified, score, identitythreshold, confidence, DB,
rtsec, total_spectrum_intensity, mass_to_charge, spectrum.charge, accession, start, end, sequence, modified_sequence FROM
(spectrum LEFT JOIN scan ON spectrum.spectrumid = scan.l_spectrumid 
LEFT JOIN identification ON spectrum.spectrumid = identification.l_spectrumid
LEFT JOIN spectrum_file ON spectrum.spectrumid = spectrum_file.l_spectrumid)
WHERE l_projectid BETWEEN 2301 AND 2400 AND l_instrumentid BETWEEN 8 AND 14;"));

identified <- rbind(identified,rt_project2m3[(!is.na(modified_sequence))])
save(rt_project2m3, file = "C:/Users/Nicolas Housset/Documents/rt_project2m3.RData", compression_level = 1)
rm(rt_project2m3)
gc()


rt_project2m4 <- data.table(dbGetQuery(con,
                                       "SELECT scanid, number, spectrumid, l_lcrunid, l_projectid, l_instrumentid, identified, score, identitythreshold, confidence, DB,
rtsec, total_spectrum_intensity, mass_to_charge, spectrum.charge, accession, start, end, sequence, modified_sequence FROM
(spectrum LEFT JOIN scan ON spectrum.spectrumid = scan.l_spectrumid 
LEFT JOIN identification ON spectrum.spectrumid = identification.l_spectrumid
LEFT JOIN spectrum_file ON spectrum.spectrumid = spectrum_file.l_spectrumid)
WHERE l_projectid BETWEEN 2401 AND 2500 AND l_instrumentid BETWEEN 8 AND 14;"));

identified <- rbind(identified,rt_project2m4[(!is.na(modified_sequence))])
save(rt_project2m4, file = "C:/Users/Nicolas Housset/Documents/rt_project2m4.RData", compression_level = 1)
rm(rt_project2m4)
gc()


rt_project2m5 <- data.table(dbGetQuery(con,
                                       "SELECT scanid, number, spectrumid, l_lcrunid, l_projectid, l_instrumentid, identified, score, identitythreshold, confidence, DB,
rtsec, total_spectrum_intensity, mass_to_charge, spectrum.charge, accession, start, end, sequence, modified_sequence FROM
(spectrum LEFT JOIN scan ON spectrum.spectrumid = scan.l_spectrumid 
LEFT JOIN identification ON spectrum.spectrumid = identification.l_spectrumid
LEFT JOIN spectrum_file ON spectrum.spectrumid = spectrum_file.l_spectrumid)
WHERE l_projectid BETWEEN 2501 AND 2600 AND l_instrumentid BETWEEN 8 AND 14;"));

identified <- rbind(identified,rt_project2m5[(!is.na(modified_sequence))])
save(rt_project2m5, file = "C:/Users/Nicolas Housset/Documents/rt_project2m5.RData", compression_level = 1)
rm(rt_project2m5)
gc()


rt_project2m6 <- data.table(dbGetQuery(con,
                                       "SELECT scanid, number, spectrumid, l_lcrunid, l_projectid, l_instrumentid, identified, score, identitythreshold, confidence, DB,
rtsec, total_spectrum_intensity, mass_to_charge, spectrum.charge, accession, start, end, sequence, modified_sequence FROM
(spectrum LEFT JOIN scan ON spectrum.spectrumid = scan.l_spectrumid 
LEFT JOIN identification ON spectrum.spectrumid = identification.l_spectrumid
LEFT JOIN spectrum_file ON spectrum.spectrumid = spectrum_file.l_spectrumid)
WHERE l_projectid BETWEEN 2601 AND 2700 AND l_instrumentid BETWEEN 8 AND 14;"));

identified <- rbind(identified,rt_project2m6[(!is.na(modified_sequence))])
save(rt_project2m6, file = "C:/Users/Nicolas Housset/Documents/rt_project2m6.RData", compression_level = 1)
rm(rt_project2m6)
gc()


rt_project2m7 <- data.table(dbGetQuery(con,
                                       "SELECT scanid, number, spectrumid, l_lcrunid, l_projectid, l_instrumentid, identified, score, identitythreshold, confidence, DB,
rtsec, total_spectrum_intensity, mass_to_charge, spectrum.charge, accession, start, end, sequence, modified_sequence FROM
(spectrum LEFT JOIN scan ON spectrum.spectrumid = scan.l_spectrumid 
LEFT JOIN identification ON spectrum.spectrumid = identification.l_spectrumid
LEFT JOIN spectrum_file ON spectrum.spectrumid = spectrum_file.l_spectrumid)
WHERE l_projectid BETWEEN 2701 AND 2800 AND l_instrumentid BETWEEN 8 AND 14;"));

identified <- rbind(identified,rt_project2m7[(!is.na(modified_sequence))])
save(rt_project2m7, file = "C:/Users/Nicolas Housset/Documents/rt_project2m7.RData", compression_level = 1)
rm(rt_project2m7)
gc()


rt_project2m8 <- data.table(dbGetQuery(con,
                                       "SELECT scanid, number, spectrumid, l_lcrunid, l_projectid, l_instrumentid, identified, score, identitythreshold, confidence, DB,
rtsec, total_spectrum_intensity, mass_to_charge, spectrum.charge, accession, start, end, sequence, modified_sequence FROM
(spectrum LEFT JOIN scan ON spectrum.spectrumid = scan.l_spectrumid 
LEFT JOIN identification ON spectrum.spectrumid = identification.l_spectrumid
LEFT JOIN spectrum_file ON spectrum.spectrumid = spectrum_file.l_spectrumid)
WHERE l_projectid BETWEEN 2801 AND 2900 AND l_instrumentid BETWEEN 8 AND 14;"));

identified <- rbind(identified,rt_project2m8[(!is.na(modified_sequence))])
save(rt_project2m8, file = "C:/Users/Nicolas Housset/Documents/rt_project2m8.RData", compression_level = 1)
rm(rt_project2m8)
gc()


rt_project2m9 <- data.table(dbGetQuery(con,
                                       "SELECT scanid, number, spectrumid, l_lcrunid, l_projectid, l_instrumentid, identified, score, identitythreshold, confidence, DB,
rtsec, total_spectrum_intensity, mass_to_charge, spectrum.charge, accession, start, end, sequence, modified_sequence FROM
(spectrum LEFT JOIN scan ON spectrum.spectrumid = scan.l_spectrumid 
LEFT JOIN identification ON spectrum.spectrumid = identification.l_spectrumid
LEFT JOIN spectrum_file ON spectrum.spectrumid = spectrum_file.l_spectrumid)
WHERE l_projectid BETWEEN 2901 AND 3000 AND l_instrumentid BETWEEN 8 AND 14;"));

identified <- rbind(identified,rt_project2m9[(!is.na(modified_sequence))])
save(rt_project2m9, file = "C:/Users/Nicolas Housset/Documents/rt_project2m9.RData", compression_level = 1)
rm(rt_project2m9)
gc()


save(identified, file = "C:/Users/Nicolas Housset/Documents/identified.RData", compression_level = 1)
