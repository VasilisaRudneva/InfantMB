# Vasilisa A. Rudneva
# Nov 2018
# 
# This function downloads raw DNA Methylation data (idat files) from a directory and returns an object of class MethylSet 

CreateRGSet <- function(path_to_idats) {

  # Get a list of files from the directory
  filenames=list.files(path = path_to_idats, recursive = T, pattern ="_Red.idat")
  
  # Make targets
  targets=data.frame(matrix("", nrow = length(filenames), ncol = 3), stringsAsFactors = F)
  colnames(targets)=c("Basename", "folder", "ArrayTypes")
  for (i in 1:length(filenames)){
    this_line=filenames[i]
    this_res=unlist(strsplit(this_line, "/"))
    
    if (grepl("_Red.idat.gz",this_res[length(this_res)])){
      targets[i,"Basename"]=gsub("_Red.idat.gz", "", this_res[length(this_res)])
    } else {
      targets[i,"Basename"]=gsub("_Red.idat", "", this_res[length(this_res)])
    }

    targets[i,"Basename"]=gsub("_Red.idat", "", this_res[length(this_res)])
    targets[i,"folder"]=paste(as.vector(this_res[-length(this_res)]), collapse = "/")
    targets[i,"ArrayTypes"]="NA"
  }
  
  # Get array types
  capture.output(read.metharray.exp(base=paste0(path_to_idats, targets$folder), targets =targets, force=T), file="tmp.txt", type = c("output", "message"))
  tmp=read.table("tmp.txt", skip = 2, header = T, sep="\t", stringsAsFactors = F)
  targets$ArrayTypes=matrix(unlist(strsplit(tmp$array..........................size, " +")), ncol=3, byrow=TRUE)[,2]; rm(tmp)

  # Load different array types separately
  array_types=unique(targets$ArrayTypes)
  this_array_type=array_types[1]
  print(this_array_type)
  targets.tmp=targets[targets$ArrayTypes==this_array_type,]
  rgset1 = read.metharray.exp(base=paste0(path_to_idats, targets.tmp$folder), targets =targets.tmp, force=T)
  rgset1$ArrayTypes=rep(this_array_type, dim(rgset1)[2])

  this_array_type=array_types[2]
  targets.tmp=targets[targets$ArrayTypes==this_array_type,]
  rgset2 = read.metharray.exp(base=paste0(path_to_idats, targets.tmp$folder), targets =targets.tmp, force=T)
  rgset2$ArrayTypes=rep(this_array_type, dim(rgset2)[2])

  # Combine together different array types
  rgSet=combineArrays(rgset1, rgset2); rm(rgset1); rm(rgset2)
  
  return(rgSet)
 }
