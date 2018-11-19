# Vasilisa A. Rudneva
# Nov 2018
# 
# A workflow to analyze DNA Methylation data

rm(list=ls())
dir=getwd()
setwd(dir)
path_to_idats="/Users/rudneva/Documents/InfantMB/data/IDATs"
masterTable=read.table("/Users/rudneva/Documents/InfantMB/data/SJYC07-Metadata.txt", header = T, sep="\t", stringsAsFactors = F)

libs=c("minfi", "limma", "minfiData", "stringr", "Rtsne", "weights")
lapply(libs, require, character.only = TRUE)

rgSet=CreateRGSet(path_to_idats)

material=masterTable[match(sampleNames(rgSet), masterTable$METH_450K),c("METH_450K", "Material")]; colnames(material)=c("id", "mat")

bVals=rgSetToBetasFiltering(rgSet, material)

Y=RuntSNE(bVals)

i.dbscan <- dbscan::dbscan(Y, eps = 10, minPts = 10)$cluster
names(i.dbscan) <- rownames(Y)
table(i.dbscan)

pdf("Infants_TSNE.All_Subgroups.pdf", height=10, width=8, onefile=TRUE, useDingbats=FALSE)
SUBGROUP = c("SHH"="red3", "Group_3"="darkgoldenrod2", "Group_4"="darkgreen")
plot(Y, cex=0.75, col=SUBGROUP[masterTable[match(rownames(Y), masterTable$METH_450K),]$Subgroup],las=2, 
          main=paste(dim(Y)[1], " Infant MB samples\n","11,664 most variable probes (SD > 0.25)", sep=""),
          pch=19)
legend("topright", legend = c("SHH", "Group_3", "Group_4"), col = SUBGROUP, pch=19)
dev.off()
