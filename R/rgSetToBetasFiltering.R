# Vasilisa A. Rudneva
# Nov 2018
# 
# This function takes an object of class MethylSet and a data.frame with samples material for batch effect correction as input
# > head(material)
#                 id    mat
#1 9297949037_R03C02 Frozen
#2 9297949037_R05C02 Frozen
#3 9305216163_R01C01 Frozen
#4 9305216153_R05C02 Frozen
#5 9297949040_R01C01 Frozen
#6 9297949147_R06C02 Frozen
#
# Preprocesses, removes batch effects, and filters the methylation data according to the array type
# Returs beta values

rgSetToBetasFiltering <- function(rgSet, material) {
  
  mSetRaw <- preprocessRaw(rgSet)
  
  Meth=log2(minfi::getMeth(mSetRaw))
  Unmeth=log2(minfi::getUnmeth(mSetRaw))
  
  data <- data.frame(id=colnames(Meth)); data$vector1 <- material[match(data$id, material$id),]
  dim(data); head(data)
  
  Meth=Meth[,as.vector(na.omit(data$vector1$id))]
  Unmeth=Unmeth[,as.vector(na.omit(data$vector1$id))]
  
  if (length(unique(data$vector1$mat))>1){
    
    Meth2=removeBatchEffect(Meth, batch = as.vector(data$vector1$mat))
    Unmeth2=removeBatchEffect(Unmeth, batch = as.vector(data$vector1$mat))
    
    Meth=2^Meth2
    Unmeth=2^Unmeth2
  }
  
  mSetRaw=MethylSet(Meth, Unmeth)
   
  # load appropriate library
  if (annotation(rgSet)[[1]] == "IlluminaHumanMethylationEPIC"){
    suppressMessages(library("IlluminaHumanMethylationEPICanno.ilm10b3.hg19"))
    suppressMessages(library("IlluminaHumanMethylationEPICmanifest")) 
    annThisArrayType = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b3.hg19)
  } else {
    suppressMessages(library("IlluminaHumanMethylation450kanno.ilmn12.hg19")) 
    suppressMessages(library("IlluminaHumanMethylation450kmanifest")) 
    annThisArrayType = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  }
  
  # Filtering (adopted from https://github.com/sirselim/illumina450k_filtering)
  # remove probes on the sex chromosomes
  keep <- !(featureNames(mSetRaw) %in% annThisArrayType$Name[annThisArrayType$chr %in% c("chrX","chrY")])
  table(keep)
  mSetRaw <- mSetRaw[keep,]; dim(mSetRaw)
  # exclude cross reactive probes
  xReactiveProbes <- read.csv(file="MethylationProbes_to_filter/450k/48639-non-specific-probes-Illumina450k.csv", stringsAsFactors=FALSE) 
  keep <- !(featureNames(mSetRaw) %in% xReactiveProbes$TargetID)
  table(keep)
  mSetRaw <- mSetRaw[keep,]
  dim(mSetRaw)
  # remove probes with any SNP
  snps=getSnpInfo(rgSet)
  remove=rownames(snps[!is.na(snps$SBE_rs) | !is.na(snps$CpG_rs),])
  keep <- !(featureNames(mSetRaw) %in% remove)
  table(keep)
  mSetRaw <- mSetRaw[keep,]
  dim(mSetRaw)
  # remove probes that are not uniquely mapped to the hg19 genome
  # BOWTIE2 multi-mapped
  multi.map <- read.csv('MethylationProbes_to_filter/450k/HumanMethylation450_15017482_v.1.1_hg19_bowtie_multimap.txt', head = F, as.is = T)
  multi.map.probes <- as.character(multi.map$V1)
  keep <- !(featureNames(mSetRaw) %in% multi.map.probes)
  table(keep)
  mSetRaw <- mSetRaw[keep,]
  dim(mSetRaw)
  
  # Additional probes for EPIC arrays
  if (annotation(rgSet)[[1]] == "IlluminaHumanMethylationEPIC"){
    suppressMessages(library("IlluminaHumanMethylationEPICanno.ilm10b3.hg19"))
    suppressMessages(library("IlluminaHumanMethylationEPICmanifest")) 
    # probes from Pidsley 2016 (EPIC)
    epic.cross1 <- read.csv('MethylationProbes_to_filter/EPIC/13059_2016_1066_MOESM1_ESM.csv', head = T)
    epic.variants1 <- read.csv('MethylationProbes_to_filter/EPIC/13059_2016_1066_MOESM4_ESM.csv', head = T)
    epic.variants2 <- read.csv('MethylationProbes_to_filter/EPIC/13059_2016_1066_MOESM5_ESM.csv', head = T)
    epic.variants3 <- read.csv('MethylationProbes_to_filter/EPIC/13059_2016_1066_MOESM6_ESM.csv', head = T)
    # additional filter probes
    epic.add.probes <- c(as.character(epic.cross1$X), as.character(epic.variants1$PROBE), as.character(epic.variants2$PROBE), 
                         as.character(epic.variants3$PROBE))
    # final list of unique probes
    epic.add.probes <- unique(epic.add.probes)
    
    keep <- !(featureNames(mSetRaw) %in% epic.add.probes)
    table(keep)
    mSetRaw <- mSetRaw[keep,]
    dim(mSetRaw)
  }
  
  bVals <- getBeta(mSetRaw, offset=100)
  names=colnames(bVals)
  
  return(bVals)
}
