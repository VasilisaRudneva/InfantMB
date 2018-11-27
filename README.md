## Description
Code that reproduces figures from the [Robinson GW &amp; Rudneva VA et al., Lancet Oncology 2018](https://www.thelancet.com/journals/lanonc/article/PIIS1470-2045(18)30204-3/fulltext) paper

View genetic alterations and download metadata from [PeCan](https://pecan.stjude.cloud/proteinpaint/study/MB-SJYC07):
`Data > Downloads > Matrix TSV`

Preprocess metadata using R
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

