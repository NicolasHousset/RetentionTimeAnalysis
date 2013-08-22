# Script to extract Enolase database

library(RMySQL);
library(data.table);

con <- dbConnect(MySQL(), group="ENOLASE", dbname="enolase_db");

projects <- data.table(dbGetQuery(con,"SELECT * FROM project"))

# OK, 31 projects, this is a little bit disappointing.