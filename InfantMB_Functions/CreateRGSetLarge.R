# Vasilisa A. Rudneva
# Nov 2018
# 
# This function does the same as CreateRGSet but is more convenient for large datasets 

CreateRGSet <- function(path_to_idats) {

  # Get a list of files from the directory
  filenames=list.files(path = path_to_idats, recursive = T, pattern ="_Red.idat")
  
  # Make targets
  targets=data.frame(matrix("", nrow = length(filenames), ncol = 6), stringsAsFactors = F)
  colnames(targets)=c("Basename", "Sample_Name", "Sample_Plate", "filenames", "folder", "ArrayTypes")
  for (i in 1:length(filenames)){
    this_line=filenames[i]
    this_res=unlist(strsplit(this_line, "/"))
    
    if (grepl("_Red.idat.gz",this_res[length(this_res)])){
      targets[i,"Basename"]=gsub("_Red.idat.gz", "", this_res[length(this_res)])
      targets[i,"Sample_Name"]=gsub("_Red.idat.gz", "", this_res[length(this_res)])
      targets[i,"Sample_Plate"]=gsub("_Red.idat.gz", "", this_res[length(this_res)])
      targets[i,"filenames"]=gsub("_Red.idat.gz", "", this_res[length(this_res)])
    } else {
      targets[i,"Basename"]=gsub("_Red.idat", "", this_res[length(this_res)])
      targets[i,"Sample_Name"]=gsub("_Red.idat", "", this_res[length(this_res)])
      targets[i,"Sample_Plate"]=gsub("_Red.idat", "", this_res[length(this_res)])
      targets[i,"filenames"]=gsub("_Red.idat", "", this_res[length(this_res)])
    }
    targets[i,"folder"]=paste(as.vector(this_res[-length(this_res)]), collapse = "/")
    targets[i,"ArrayTypes"]="NA"
  }
  
  # Get array types
  capture.output(read.metharray.exp(base=paste0(path_to_idats, targets$folder), targets =targets, force=T), file="tmp.txt", type = c("output", "message"))
  tmp=read.table("tmp.txt", skip = 2, header = T, sep="\t", stringsAsFactors = F)
  targets$ArrayTypes=matrix(unlist(strsplit(tmp$array..........................size, " +")), ncol=3, byrow=TRUE)[,2]; rm(tmp)

  # Initiate with two first folders
  folder1=unique(targets$folder)[1]
  targets.tmp=targets[targets$folder==folder1,]
  rgSet1= read.metharray.exp(base=paste0(path_to_idats, targets.tmp$Sample_Plate, "/", folder1), 
                             targets =targets.tmp, force=T)

  folder2=unique(targets$folder)[2]
  targets.tmp=targets[targets$folder==folder2,]
  rgSet2= read.metharray.exp(base=paste0(path_to_idats, targets.tmp$Sample_Plate, "/", folder2), 
                             targets =targets.tmp, force=T)

  rgSetInit=combineArrays(rgSet1, rgSet2)
  rm(rgSet1)
  rm(rgSet2)

  len=length(unique(targets$folder))
  for (ind in 3:len){
    folder=unique(targets$folder)[ind]
    targets.tmp=targets[targets$folder==folder,]
    this_rgset = read.metharray.exp(base=paste0(path_to_idats, targets.tmp$Sample_Plate, "/", folder), 
                                    targets =targets.tmp, force=T)

    rgSetnew=combineArrays(rgSetInit, this_rgset); rm(this_rgset); rm(rgSetInit)
    rgSetInit=rgSetnew; rm(rgSetnew)
  }

  return(rgSetInit)
 }


