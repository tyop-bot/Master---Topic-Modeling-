library(tidyverse)

samplemeta <- read.table(snakemake@input[["premeta"]], header=TRUE, sep="\t")

#If lots of NA in col "Subtype_DNAmeth", then col not added to metaframe:

#BRCA-US: 'Subtype_mRNA' & 'Subtype_Immune_Model_Based'

metaframe <- samplemeta %>% select(icgc_donor_id, Subtype_mRNA, Subtype_Immune_Model_Based)

# if ((sum(is.na(samplemeta$Subtype_DNAmeth)))>(nrow(samplemeta)-10)) {
# metaframe <- samplemeta %>% select(icgc_donor_id, Subtype_Immune_Model_Based)
# } else {
# metaframe <- samplemeta %>% select(icgc_donor_id, Subtype_DNAmeth, Subtype_Immune_Model_Based)
# }

#replace NA with string
metaframe.c <- metaframe %>% mutate_all(as.character)
metaframe.c <- metaframe.c %>% replace(is.na(.), "NA")

#If redundancy in data (E.g same donors multiple times)
metaSD <- metaframe.c %>% distinct(icgc_donor_id, .keep_all = TRUE)

#set rownames
meta2 <- metaSD
rownames(meta2) <- metaSD[,1]
meta3 <- meta2 %>% select(-1)

saveRDS(meta3, file=snakemake@output[["postmeta"]])

##Not neeeded here, transport to after cpgs done!!!
#Removing from main DF columns that is not in metadatas rownames
#drop <- c(rownames(meta2))
#CpGDropped = CpG[,(names(CpG) %in% drop)]
