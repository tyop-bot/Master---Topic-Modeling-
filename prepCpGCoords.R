CpGs <- readRDS(snakemake@input[["preCpG"]])
#Changed from read.table bcuz of BED being 0-based. This to avoid any problems since R is 1-based.
Coords <- read.table(snakemake@input[["coords"]])

library(tidyverse)
mergedCol <- unite(Coords, comb, V2, V3, sep = "-", remove = TRUE, na.rm = FALSE)
mergedChrDone <- unite(mergedCol, Coords, V1, comb, sep = ":", remove = TRUE, na.rm = FALSE)


names(CpGs)[names(CpGs) == "X"] <- "probes"
CpGs$probes <- as.character(CpGs$probes)

#Should not be needed cus IMPUTATION
#cpg <- CpG.sample.tab[complete.cases(CpG.sample.tab), ]

names(mergedChrDone)[names(mergedChrDone) == "V4"] <- "probes"
mergedChrDone$probes <-as.character(mergedChrDone$probes)

merged <- left_join(mergedChrDone, CpGs, by = c("probes"))

mergedSansprobes <- within(merged, rm("probes"))

merged2 <- mergedSansprobes[,-1]
rownames(merged2) <- mergedSansprobes[,1]

saveRDS(merged2, file=snakemake@output[["preCpG"]])

#will make df look like this:
#                   DO1274    DO1275
#chr1:15864-15865   0.8594239 0.8634571
#chr1:18826-18827   0.6293810 0.6655710
#chr1:29424-29425   0.2301290 0.2686210
#chr1:69590-69591   0.5033606 0.6797580
#chr1:135251-135252 0.8083419 0.6473794
#chr1:568535-568536 0.1705664 0.1788103
