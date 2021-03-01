library(rGREAT)
library(gridExtra)

args = commandArgs(trailingOnly=TRUE)
bedfile <- toString(args[1])
print(bedfile)

bedJob = submitGreatJob(bedfile, species="hg38")
enrichTable = getEnrichmentTables(bedJob)

df <- print(enrichTable[[1]])
# df <- enrichTable[[1]]

pdfname <- paste(bedfile,"_GREAT.pdf", sep="")
pdf(pdfname, height=4, width=40)




grid.table(df)
res = plotRegionGeneAssociationGraphs(bedJob)

plot(NA, xlim=c(0,10), ylim=c(0,10), bty='n',
     xaxt='n', yaxt='n', xlab='', ylab='')
text(4,9, bedfile, pos=4)

dev.off()
