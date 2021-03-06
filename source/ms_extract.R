# The purpose of this script is to extract relevant information concerning projects with id between 1500 and 3000.
# Reason behind this range of numbers is that projects are recent
# Projects are extracted by chunks of 100 to avoid overloading the RAM
# A db of identified peptides is built step by step

library(RMySQL);
library(data.table);

con <- dbConnect(MySQL(), group="MSDB", dbname="projects");

# Since it is a function the data.table will be created in the environment of the function
sampleExtract <- function(saveName, projectStart, projectEnd, ...){
  projectPath <- "C:/Users/Nicolas Housset/Documents/RetentionTimeAnalysis";
  # Beginning of the SQL statement which contains the variable we want to extract
  varSQL <- "\"SELECT scanid, number, spectrumid, l_lcrunid, l_projectid, l_instrumentid, l_protocolid, l_userid, identified, score, identitythreshold, confidence, DB,
rtsec, total_spectrum_intensity, mass_to_charge, spectrum.charge, accession, start, end, sequence, modified_sequence, identification.description FROM
(spectrum LEFT JOIN scan ON spectrum.spectrumid = scan.l_spectrumid 
LEFT JOIN identification ON spectrum.spectrumid = identification.l_spectrumid
RIGHT JOIN project ON spectrum.l_projectid = project.projectid) 
WHERE l_projectid BETWEEN ";
  varSQL <- paste0(saveName,
                   " <- data.table(dbGetQuery(con,",
                   varSQL, projectStart, " AND ", projectEnd, " AND l_instrumentid BETWEEN 8 AND 14;\"))");
  
  # print(varSQL);
  eval(parse(text=varSQL));
  saveText <- paste0("save(",saveName, ", file = \"", projectPath, "/data/", saveName, ".RData\", compression_level = 1)");
  # print(saveText);
  eval(parse(text=saveText));
  return(NULL);
}

# From project 1501, extract and save in chunks of 100 projects
start <- 1500
for(i in 1:15){
  start <- start + 1
  end <- start + 99
  stringExpression <- paste0("sampleExtract(\"rt_project_part", i, "\", \"", start, "\", \"", end, "\")")
  start <- end
  print(stringExpression)
  eval(parse(text=stringExpression))
}

# Our database is divided into parts, we build step by step the identified peptides database
identified <- data.table(NULL)
projectPath <- "C:/Users/Nicolas Housset/Documents/RetentionTimeAnalysis"
for(i in 1:15){
  stringExpression <- paste0("load(file = \"", projectPath, "/data/rt_project_part", i, ".RData\")" )
  # print(stringExpression)
  eval(parse(text=stringExpression))
  stringExpression <- paste0("identified <- rbind(identified,rt_project_part", i, "[(!is.na(modified_sequence))])")
  # print(stringExpression)
  eval(parse(text=stringExpression))
  stringExpression <- paste0("rm(rt_project_part", i, ")")
  # print(stringExpression)
  eval(parse(text=stringExpression))
}
save(identified, file = "C:/Users/Nicolas Housset/Documents/RetentionTimeAnalysis/data/identified.RData", compression_level = 1)


#
#
#
#
# From there, do not run, previous code
rt_project2m0 <- data.table(dbGetQuery(con,
                                       "SELECT scanid, number, spectrumid, l_lcrunid, l_projectid, l_instrumentid, identified, score, identitythreshold, confidence, DB,
rtsec, total_spectrum_intensity, mass_to_charge, spectrum.charge, accession, start, end, sequence, modified_sequence FROM
(spectrum LEFT JOIN scan ON spectrum.spectrumid = scan.l_spectrumid 
                                       LEFT JOIN identification ON spectrum.spectrumid = identification.l_spectrumid
                                       LEFT JOIN spectrum_file ON spectrum.spectrumid = spectrum_file.l_spectrumid)
                                       WHERE l_projectid BETWEEN 2101 AND 2200 AND l_instrumentid BETWEEN 8 AND 14;"));

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
