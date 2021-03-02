##Dict to go through matrices and appropriate transcription factors

projectTfDict = {
"BLCA-US": ["FOXA1", "GATA3"],
"BRCA-US": ["FOXA1", "ESR1", "GATA3", "RUNX3"],
"HNSC-US": ["SOX2"],
"KIRP-US": ["SPI1", "FLI1"],
"LGG-US": ["BHLHE40", "TCF12", "TFAP4", "ZEB1"],
"LUAD-US": ["FOXA2"],
"LUSC-US": ["SOX2", "TP63"],
"PRAD-US": ["TP63", "SNAI2"],
"SKCM-US": ["SPI1"],
"STAD-US": ["HNF4A", "GRHL2", "GATA3", "ELF3"],
"THCA-US": ["RUNX1", "SPI1", "PBX3"],
"UCEC-US": ["FOXA2"]
}


def findProjectandTF():
    for project in projectTfDict:
        #rename full probe list
        return(PROJECT=project)
        for tfactor in projectTfDict[project]:
            #return(project, tfactor)
            TFACTOR=tfactor
            return(PROJECT, TFACTOR)



print("Dette gikk bra")
print(TFACTOR)
print(findProjectandTF())
print(findProjectandTF())
#print(/storage/mathelierarea/processed/rozabl/Methylation_TF/results/20180802_ICGC_UniBind_all_corr/project/project_TRs/tfactor/project_tfactor.bed.hg19.wgEncodeHaibMethyl450CpgIslandDetails_emap.probes.bed)
#/storage/scratch/rozabl/reTFme/results/20190504_ICGC_UniBind_all_corr/{project}/{project}_TRs/{tfactor}/{project}_{tfactor}.bed.hg19.wgEncodeHaibMethyl450CpgIslandDetails_emap.probes.bed
