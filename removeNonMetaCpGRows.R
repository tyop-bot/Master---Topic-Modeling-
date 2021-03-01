CpG <- readRDS(snakemake@input[["preCpG"]])
meta3 <- readRDS(snakemake@input[["postmeta"]])


#Removing from main DF columns that is not in metadatas rownames
drop <- c(rownames(meta3))
CpGDropped = CpG[,(names(CpG) %in% drop)]

saveRDS(CpGDropped, file=snakemake@output[["CpGReady4cis"]])
