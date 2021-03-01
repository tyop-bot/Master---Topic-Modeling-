library(Amelia)
library(tidyverse)


Amp <- read.table(snakemake@input[["CpG"]], header=TRUE, sep="\t")
Coords <- read.table(snakemake@input[["coords"]])

mergedCol <- unite(Coords, comb, V2, V3, sep = "-", remove = TRUE, na.rm = FALSE)
mergedChrDone <- unite(mergedCol, Coords, V1, comb, sep = ":", remove = TRUE, na.rm = FALSE)

names(Amp)[names(Amp) == "X"] <- "probes"
Amp$probes <- as.character(Amp$probes)


names(mergedChrDone)[names(mergedChrDone) == "V4"] <- "probes"
mergedChrDone$probes <-as.character(mergedChrDone$probes)

merged <- left_join(mergedChrDone, Amp, by = c("probes"))

mergedSansprobes <- within(merged, rm("probes"))

mergedSansprobes <- mergedSansprobes[which(rowMeans(!is.na(mergedSansprobes)) > 0.5), ]

idvars = c('Coords')

Amp2 <- amelia(mergedSansprobes, m = 1, idvars = idvars)

imputedMatrix <- Amp2$imputations[[1]]

imputedMatrix2 <- imputedMatrix[,-1]
rownames(imputedMatrix2) <- imputedMatrix[,1]

#############PrepMeta###############

samplemeta <- read.table(snakemake@input[["premeta"]], header=TRUE, sep="\t")


###OTHERS THAN BRCA::: _Subtype_mRNA', 'Subtype_Selected' & 'Subtype_Immune_Model_Based
metaframe <- samplemeta %>% select(icgc_donor_id, Subtype_Selected, Subtype_Immune_Model_Based)


###BRCA_FULL
#metaframe <- samplemeta %>% select(donor_id, PAM50, ER.Status)


#replace NA with string
metaframe.c <- metaframe %>% mutate_all(as.character)
metaframe.c <- metaframe.c %>% replace(is.na(.), "NA")

#If redundancy in data (E.g same donors multiple times)
metaSD <- metaframe.c %>% distinct(icgc_donor_id, .keep_all = TRUE)
#metaSD <- metaframe.c %>% distinct(donor_id, .keep_all = TRUE)

#set rownames
meta2 <- metaSD
rownames(meta2) <- metaSD[,1]


meta <- meta2 %>% select(-1)


###########RemoveNonMetaCpGRows############

drop <- c(rownames(meta))
CpG.sample.tab = imputedMatrix2[,(names(imputedMatrix2) %in% drop)]


############CISTOPIC########
suppressWarnings(library("cisTopic"))

#Name the PDF after what cancer type and TF:
PDF.name <- snakemake@output[[1]]

#Removing NAs by deletion: (CHECK IF SHOULD USE AMPUTATION!!!
#cpg <- CpG.sample.tab[complete.cases(CpG.sample.tab), ]
#Finding what to use from meta:
metacolnames <- colnames(meta)

cisTopicObject <- createcisTopicObject(is.acc=0.5, count.matrix=data.frame(CpG.sample.tab))
cisTopicObject <- addCellMetadata(cisTopicObject, cell.data = meta)
cisTopicObject <- runWarpLDAModels(cisTopicObject, topic=c(5:15,18,20,22,25), seed=123, nCores=20, addModels=FALSE)

pdf(PDF.name, height=7, width=10)
cisTopicObject <- selectModel(cisTopicObject, type='perplexity')
cisTopicObject <- runUmap(cisTopicObject, target ='cell')
cellTopicHeatmap(cisTopicObject, colorBy=c(metacolnames))
cellTopicHeatmap(cisTopicObject, method = "Probability", colorBy=c(metacolnames))

###
plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
     xaxt='n', yaxt='n', xlab='', ylab='')
text(1,4,"Colored by metadata", pos=4)

par(mfrow=c(1,2))
plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr=NULL, colorBy=c(metacolnames), cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE)
###
plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
     xaxt='n', yaxt='n', xlab='', ylab='')
text(1,4,"Topic by probability", pos=4)

par(mfrow=c(3,3))
plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr='Probability', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE)
###
plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
     xaxt='n', yaxt='n', xlab='', ylab='')
text(1,4,"Topic by Z-score", pos=4)

par(mfrow=c(3,3))
plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr='Z-score', colorBy=NULL, cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE)


cisTopicObject <- getRegionsScores(cisTopicObject, method='NormTop', scale=TRUE)
par(mfrow=c(3,5))
cisTopicObject <- binarizecisTopics(cisTopicObject, thrP=0.975, plot=TRUE)


cisTopicObject <- GREAT(cisTopicObject, genome='hg19', fold_enrichment=2, geneHits=1, sign=0.05, request_interval=10)
ontologyDotPlot(cisTopicObject, top=5, var.y='name', order.by='Binom_Adjp_BH')


# plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
#      xaxt='n', yaxt='n', xlab='', ylab='')
# topicNumber <- tail(strsplit(getwd() ,split="/")[[1]],1)
#text(1,4, getwd() , pos=4)

#should notify snakemake that pdf is output
dev.off()

saveRDS(cisTopicObject, file=snakemake@output[[2]])
#getBedFiles(cisTopicObject, path='CisTopicBed_BRCA_SNA')
