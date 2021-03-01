suppressWarnings(library("cisTopic"))
CpG.sample.tab <- readRDS(snakemake@input[["CpGReady4cis"]])
meta <- readRDS(snakemake@input[["meta"]])

####Drop CpGs
#drop <- c(rownames(meta2))
#CpGDropped = CpG[,(names(CpG) %in% drop)]

#Name the PDF after what cancer type and TF:
PDF.name <- snakemake@output[[1]]

#Removing NAs by deletion: (CHECK IF SHOULD USE AMPUTATION!!!
#cpg <- CpG.sample.tab[complete.cases(CpG.sample.tab), ]
#Finding what to use from meta:
metacolnames <- colnames(meta)

cisTopicObject <- createcisTopicObject(is.acc=0.5, count.matrix=data.frame(cpg))
cisTopicObject <- addCellMetadata(cisTopicObject, cell.data = meta)
cisTopicObject <- runWarpLDAModels(cisTopicObject, topic=c(5:15,18,20,22,25), seed=123, nCores=20, addModels=FALSE)

pdf(PDF.name)
cisTopicObject <- selectModel(cisTopicObject, type='derivative')
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





# plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
#      xaxt='n', yaxt='n', xlab='', ylab='')
# topicNumber <- tail(strsplit(getwd() ,split="/")[[1]],1)
#text(1,4, getwd() , pos=4)

#should notify snakemake that pdf is output
dev.off()

saveRDS(cisTopicObject, file=snakemake@output[[2]])
#getBedFiles(cisTopicObject, path='CisTopicBed_BRCA_SNA')

######
