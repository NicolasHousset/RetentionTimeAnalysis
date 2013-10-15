# The problem : we want to synthetise artificial peptides based on the sequence NH2-...XXX...-COOH
# We want to cover the possible combinations of XXX (20*20*20 = 8000 combinations)
# To generate an artificial peptide, we need a "bit" as material. How many "bits" are necessary to cover 99 % of the peptides ?
# With 25000 bits -> around 7650 peptides are covered
# With 50000 bits -> around 7980 peptides are covered 



simulation <- data.table(1:800000)


simulation[, alea := runif(800000)]

simulation[, peptide := as.character(ceiling(alea*8000))]

simulation[, group := as.character(ceiling((1:800000)/50000))]
setkey(simulation, group, peptide)

uniquepep <- unique(simulation)

setkey(uniquepep, group)

coverage_list <- vector("list", 16)
for(i in 1:16){
  coverage_list[[i]] <- NROW(uniquepep[as.character(i)])
}
NROW(uniquepep[as.character(i)])

table(simulation[, group])
