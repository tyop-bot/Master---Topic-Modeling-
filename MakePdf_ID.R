#Create new PDF with name of folder to identify where the swarmplot is from
pdf("new.pdf")
plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n',
     xaxt='n', yaxt='n', xlab='', ylab='')
topicNumber <- tail(strsplit(getwd() ,split="/")[[1]],1)
text(1,4, topicNumber , pos=4)

#Merge new PDF with
