library(Amelia)
Amp <- read.table(snakemake@input[["CpG"]], header=TRUE, sep="\t")

#Remove rows with more than 50% NA
Amp <- Amp[which(rowMeans(!is.na(Amp)) > 0.5), ]

idvars = c('X')
Amp2 <- amelia(Amp, m = 1, idvars = idvars)


imputed <- Amp2$imputations[[1]]
saveRDS(imputed, snakemake@output[["imputedCpG"]])
