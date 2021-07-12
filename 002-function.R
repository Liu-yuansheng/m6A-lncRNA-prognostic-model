
##########################处理Json文件的函数##########################
#jsonFile是完整的文件路径
processingJsonFiles <- function(jsonFile){
  library(rjson)
  metadata_json_File <- fromJSON(file=jsonFile)
  json_File_Info <- data.frame(filesName = c(),TCGA_Barcode = c())
  for(i in 1:length(metadata_json_File)){
    TCGA_Barcode <- metadata_json_File[[i]][["associated_entities"]][[1]][["entity_submitter_id"]]
    file_name <- metadata_json_File[[i]][["file_name"]]
    json_File_Info <- rbind(json_File_Info,data.frame(filesName = file_name,TCGA_Barcode = TCGA_Barcode))
  }
  
  rownames(json_File_Info) <- json_File_Info[,1]
  json_File_Info <-json_File_Info[-1]
  return(json_File_Info)
}


############## 一个函数，通过gtf文件获取Ensemble_ID与基因名称的对应关系
get_map <- function(input) {
  if (is.character(input)) {
    if(!file.exists(input)) stop("Bad input file.")
    message("Treat input as file")
    input = data.table::fread(input, header = FALSE)
  } else {
    data.table::setDT(input)
  }
  input = input[input[[3]] == "gene", ]
  pattern_id = ".*gene_id \"([^;]+)\";.*"
  pattern_name = ".*gene_name \"([^;]+)\";.*"
  pattern_type = ".*gene_type \"([^;]+)\";.*"
  gene_id = sub(pattern_id, "\\1", input[[9]])
  gene_name = sub(pattern_name, "\\1", input[[9]])
  gene_type = sub(pattern_type, "\\1", input[[9]])
  EnsemblTOGenename <- data.frame(Ensembl = gene_id,
                                  gene_name = gene_name,
                                  gene_type = gene_type,
                                  stringsAsFactors = FALSE)
  return(EnsemblTOGenename)
}


############一个函数，一个Barcode向量，区分正常组织和肿瘤组织
getSampleInfo <- function(Barcode){
  TumorBarcode <- Barcode[grep("-01[A-D]-",Barcode)]
  NormalBarcode <- Barcode[-grep("-01[A-D]-",Barcode)]
  sampleNum <- data.frame(NormalBarcode = length(NormalBarcode),TumorBarcode = length(TumorBarcode))
  sampleList <- data.frame(ID = c(NormalBarcode,TumorBarcode),Type = c(rep("Normal",length(NormalBarcode)),
                                                                       rep("Tumor",length(TumorBarcode))))
  return(list(BarcodeSort = c(NormalBarcode,TumorBarcode),
              sampleNum = sampleNum ,sampleList = sampleList,
              TumorBarcode=TumorBarcode,NormalBarcode=NormalBarcode))
}
############ 上面函数返回一个重新排序的Barcode向量，正常组织在前



############获取生存时间，clindataM应为读入的tsv文件
getSurData <- function(clindataM){
  case_submitter_id <- c()
  SurvivalTime <- c()
  VitalStatus <-c()
  for(i in 1:nrow(clindataM)){
    if(clindataM[i,"vital_status"] == "Alive"){
      case_submitter_id = c(case_submitter_id,clindataM[i,"case_submitter_id"])
      SurvivalTime <- c(SurvivalTime,clindataM[i,"days_to_last_follow_up"])
      VitalStatus <- c(VitalStatus,0)
    }
    else if(clindataM[i,"vital_status"]== "Dead"){
      case_submitter_id = c(case_submitter_id,clindataM[i,"case_submitter_id"])
      SurvivalTime <- c(SurvivalTime,clindataM[i,"days_to_death"])
      VitalStatus <- c(VitalStatus,1)
    }else {print(paste("第",i,"行数据的生存状态为空值"))}
  }
  surData <- data.frame(case_submitter_id = case_submitter_id,SurvivalTime = SurvivalTime,VitalStatus = VitalStatus)
  return(surData)
}

