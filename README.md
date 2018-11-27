## Description
R code that reproduces figures from [Robinson GW &amp; Rudneva VA et al., Lancet Oncology 2018](https://www.thelancet.com/journals/lanonc/article/PIIS1470-2045(18)30204-3/fulltext)

## Workflow

View genetic alterations and download metadata from [PeCan](https://pecan.stjude.cloud/proteinpaint/study/MB-SJYC07):
`Data > Downloads > Matrix TSV`

Resulting file is [heatmap.tsv](https://github.com/VasilisaRudneva/InfantMB/blob/master/Data/heatmap.tsv).

### Preprocess metadata
```
require(data.table)
masterTable<-t(as.data.frame(fread("heatmap.tsv")))
masterTable=as.data.frame(masterTable, stringsAsFactors = F)
rownames(masterTable)=gsub("0;;", "", rownames(masterTable))
colnames(masterTable)=masterTable[1,]
masterTable=masterTable[-1,]
masterTable$PID=rownames(masterTable)
masterTable=masterTable[,-1]
```
### Set dependencies
```
rm(list=ls())
dir="/User/Documents/InfantMB/"
path_to_idats="/User/Documents/InfantMB/idats/"
setwd(dir)

libs=c("minfi", "limma", "minfiData", "stringr", "Rtsne", "weights")
lapply(libs, require, character.only = TRUE)
```
### Preprocess DNA methylation data
```
rgSet=CreateRGSet(path_to_idats)
material=masterTable[match(sampleNames(rgSet), masterTable$METH_450K),c("METH_450K", "Material")]
colnames(material)=c("id", "mat")
bVals=rgSetToBetasFiltering(rgSet, material)
```
### t-SNE structure
```
Y=RuntSNE(bVals)
color_scheme = c("SHH"="red3", "Group_3"="darkgoldenrod2", "Group_4"="darkgreen")
plot(Y, cex=0.75, col=color_scheme[masterTable[match(rownames(Y), masterTable$METH_450K),]$Subgroup], las=2, 
          main=paste(dim(Y)[1], " Infant MB samples\n", dim(bVals)[1]," most variable probes (SD > 0.25)", sep=""),
          pch=19)
legend("topright", legend = c("SHH", "Group_3", "Group_4"), col = color_scheme, pch=19)
```
### Oncoprint

### Defining clusters
```
i.dbscan <- dbscan::dbscan(Y, eps = 10, minPts = 10)$cluster
names(i.dbscan) <- rownames(Y)
table(i.dbscan)
```

