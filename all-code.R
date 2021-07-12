###------------001----------------
rm(list = ls())
options(stringsAsFactors = F)
# 在此之前用TCGAbiolinks获得转录组数据和clinical data
######----------------------1.定义相关文件路径-------------------------------
# 需要输入的文件，注意修改
#gene文件是我们要分析的基因，比如m6A,免疫调控相关的基因。gtf文件用于注释
inputPath <- ".\\02-inputFiles" #gene文件和gtf文件路径
RNASeqDataPath <- ".\\01-RawData\\FPKMandCounts\\" # TCGA 下载的转录组数据
jsonFilePath <- ".\\01-RawData\\metadata.cart.2021-07-08.json"# 转录组数据对应的json文件
gtfFilesPath <- ".\\02-inputFiles\\human.v33lift37.annotation.gtf"#基因注释文件
geneFilesPath <- ".\\02-inputFiles\\gene.txt" #要分析的基因，比如m6A RNA甲基化相关基因

# 临床数据文件地址：TCGA下载的临床数据
ClinicalDataPath <- ".\\01-RawData\\clinical.cart.2021-07-08\\clinical.tsv"

#############################上面代码地址注意修改

# setwd("./TCGA/Prognositc model analysis/")
## dir.create("00-RowData") 手动创建用存放数据
dir.create("03-AnalysisResult")
setwd(".\\03-AnalysisResult")
FolderName <- c("101-extractGenesExp","102-GenesDiff",
                "103-correlationAnalysis","104-tumorCluster",
                "105-PCA","106-ClinicalDataExtract","107-tumorClusterSurvivalAnalysis",
                "108-tumorClusterCorrelationAnalysis","109-singleFactorCOXAnlysis",
                "110-LassoRegression","111-RiskSurvivalCurve","112-ROC-Curve",
                "113-ClinicalCorrelationOfRisk","114-IndependentPrognostic",
                "115-ImmuneInfiltration"
)
mkdir <- function(FolderName){
  
  for (name in FolderName) {
    dir.create(name,showWarnings = F)
  }
}
mkdir(FolderName)
setwd("..\\")

#####################################################100-2
options(stringsAsFactors = F)
################################
# 分析结果输出路径
outputPath <- c(".\\03-AnalysisResult\\101-extractGenesExp\\",
               ".\\03-AnalysisResult\\102-GenesDiff\\",
               ".\\03-AnalysisResult\\103-CorrelationAnalysis\\",
               ".\\03-AnalysisResult\\104-tumorCluster\\",
               ".\\03-AnalysisResult\\105-PCA\\",
               ".\\03-AnalysisResult\\106-ClinicalDataExtract\\",
               ".\\03-AnalysisResult\\107-tumorClusterSurvivalAnalysis\\",
               ".\\03-AnalysisResult\\108-tumorClusterCorrelationAnalysis\\",
               ".\\03-AnalysisResult\\109-singleFactorCOXAnlysis\\",
               ".\\03-AnalysisResult\\110-LassoRegression\\",
               ".\\03-AnalysisResult\\111-RiskSurvivalCurve\\",
               ".\\03-AnalysisResult\\112-ROC-Curve\\",
               ".\\03-AnalysisResult\\113-ClinicalCorrelationOfRisk\\",
               ".\\03-AnalysisResult\\114-IndependentPrognostic\\")
####-------------------002-----------------

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

#####---------------------101----------------

options(stringsAsFactors = F)
############################将所有文件移动到同一个文件夹tempFiles下
if(file.exists("01-tempFiles")== FALSE){dir.create("01-tempFiles")}
txtfilepath <- dir(path = RNASeqDataPath,full.names = TRUE)
for(wd in txtfilepath){
  files <-dir(path = wd,pattern="gz$") #查看满足条件文件
  fromfilepath <- paste(wd,"\\",files,sep ="")
  tofilepath <- paste(".\\01-tempFiles\\",files,sep ="")
  file.copy(fromfilepath,tofilepath)
}

##processingJsonFiles函数在000-function代码中
jsonFileInfo <- processingJsonFiles(jsonFile = jsonFilePath)
write.table(jsonFileInfo,file = paste(outputPath[1],"JsonFileInfo.txt",sep = ""))


####ID转换,#gtf文件路径，注意更改
EnsemblTOGene <- get_map(gtfFilesPath) 
EnsemblTOGene$Ensembl <- substr(EnsemblTOGene[,"Ensembl"],1,15)#去掉版本号
colnames(EnsemblTOGene) <- c("Ensembl_ID","gene_id","gene_type")
write.table(EnsemblTOGene,file = paste(outputPath[1],"EnsemblTOGene.txt",sep = ""))

#可以将其保存为R对象，以后每次都不需要进行转换啦。
SampleFilesPath <- ".\\01-tempFiles"

###获取表达矩阵，样品信息已通过getSampleInfo函数排序
num = 0
setwd(SampleFilesPath)
if(file.exists("MANIFEST.txt")){file.remove("MANIFEST.txt")}

###下载的RNASeq数据包括了FPKM和counts的
###FPKM数据文件末尾都是htseq_fpkm.txt.gz，counts数据文件名末尾都是htseq_counts.gz
###先处理FPKM的数据
FPKM.FileNames<-dir(pattern="htseq_fpkm.tsv.gz$") #list.files()函数也行
geneExpMatrix.FPKM <- data.frame()
for(FPKM.File in FPKM.FileNames){
  #每一个循环读取一个文件
  print(paste("正在读入文件：",FPKM.File,sep = ""))
  OneSampleExp <- read.table(FPKM.File,header =FALSE)
  # 根据BarcodeFilesPath文件中文件名称与barcode对应关系，命名列名
  patientsID <- jsonFileInfo[FPKM.File,"TCGA_Barcode"]
  colnames(OneSampleExp) <- c("Ensembl_ID",patientsID)
  
  if (dim(geneExpMatrix.FPKM)[1]== 0){
    geneExpMatrix.FPKM <- OneSampleExp
  }else {geneExpMatrix.FPKM<- merge(geneExpMatrix.FPKM,OneSampleExp,by = "Ensembl_ID")}
  num = num +1
}
print(paste("读取文件完毕,共读入",num,"个文件",sep = ""))
###再处理counts的数据
num = 0
counts.FileNames<-dir(pattern="htseq.counts.gz$") #list.files()函数也行
geneExpMatrix.count <- data.frame()
for(count.File in counts.FileNames){
  #每一个循环读取一个文件
  print(paste("正在读入文件：",count.File,sep = ""))
  OneSampleExp <- read.table(count.File,header =FALSE)
  # 根据BarcodeFilesPath文件中文件名称与barcode对应关系，命名列名
  patientsID <- jsonFileInfo[count.File,"TCGA_Barcode"]
  colnames(OneSampleExp) <- c("Ensembl_ID",patientsID)
  
  if (dim(geneExpMatrix.count)[1]== 0){
    geneExpMatrix.count <- OneSampleExp
  }else {geneExpMatrix.count <- merge(geneExpMatrix.count,OneSampleExp,by = "Ensembl_ID")}
  num = num +1
}
print(paste("读取文件完毕,共读入",num,"个文件",sep = ""))
setwd("../") # 读完数据后返回最初工作目录


geneExpMatrix.FPKM <- geneExpMatrix.FPKM[grep("ENSG",geneExpMatrix.FPKM$Ensembl_ID),]
geneExpMatrix.count <- geneExpMatrix.count[grep("ENSG",geneExpMatrix.count$Ensembl_ID),]

library(stringr)
geneExpMatrix.FPKM$Ensembl_ID <- c(str_split(geneExpMatrix.FPKM$Ensembl_ID,'[.]',simplify = T)[,1])
geneExpMatrix.FPKM <- merge(EnsemblTOGene,geneExpMatrix.FPKM,by = "Ensembl_ID")

geneExpMatrix.count$Ensembl_ID <- c(str_split(geneExpMatrix.count$Ensembl_ID,'[.]',simplify = T)[,1])
geneExpMatrix.count <- merge(EnsemblTOGene,geneExpMatrix.count,by = "Ensembl_ID")

save(geneExpMatrix.FPKM,geneExpMatrix.count,file = paste0(outputPath[1],"COAD.RNAseq.Rdata"))


####重复基因取均值,提取lncRNA or mRNA(protein_coding)
library(limma)
proteinExpMatrix.FPKM  <- geneExpMatrix.FPKM[geneExpMatrix.FPKM$gene_type=="protein_coding",]
geneNames.FPKM <- proteinExpMatrix.FPKM$gene_id
proteinExpMatrix.FPKM  <-  proteinExpMatrix.FPKM[,-c(1:3)] %>% data.matrix()
rownames(proteinExpMatrix.FPKM ) <- as.vector(geneNames.FPKM)
proteinExpMatrix.FPKM  <- avereps(proteinExpMatrix.FPKM ) %>% as.data.frame()
dim(proteinExpMatrix.FPKM )

proteinExpMatrix.count  <- geneExpMatrix.count[geneExpMatrix.count$gene_type=="protein_coding",]
geneNames.count <- proteinExpMatrix.count$gene_id
proteinExpMatrix.count <-  proteinExpMatrix.count[,-c(1:3)] %>% data.matrix()
rownames(proteinExpMatrix.count) <- as.vector(geneNames.count)
proteinExpMatrix.count <- avereps(proteinExpMatrix.count) %>% as.data.frame()
dim(proteinExpMatrix.count)

###
# 从数据集中删除零方差的基因行，可以使用相同的apply 表达式，设置方差不等于零。
proteinExpMatrix.FPKM  <- proteinExpMatrix.FPKM [which(apply(proteinExpMatrix.FPKM , 1, var)!=0),]
dim(proteinExpMatrix.FPKM )
# 从数据集中删除零方差的基因行，可以使用相同的apply 表达式，设置方差不等于零。
proteinExpMatrix.count  <- proteinExpMatrix.count [which(apply(proteinExpMatrix.count , 1, var)!=0),]
dim(proteinExpMatrix.count )

##保存数据
##记录所有样本ID
Barcode <- intersect(colnames(proteinExpMatrix.FPKM),colnames(proteinExpMatrix.count))
sample.info <- getSampleInfo(Barcode)##Barcode
#保存文件
write.table(proteinExpMatrix.FPKM ,file = paste(outputPath[1],"proteinExpMatrix.FPKM.txt",sep = ""),
            sep = "\t",row.names = T)
#保存文件
write.table(proteinExpMatrix.count ,file = paste(outputPath[1],"proteinExpMatrix.count.txt",sep = ""),
            sep = "\t",row.names = T)
save(proteinExpMatrix.FPKM ,sample.info,file = paste0(outputPath[1],"COAD.proteinExpMatrix.FPKM.Rdata"))
save(proteinExpMatrix.count,sample.info,file = paste0(outputPath[1],"COAD.proteinExpMatrix.count.Rdata"))


####计算lncRNA的
####重复基因取均值,提取lncRNA or mRNA(protein_coding)
library(limma)
lncRNAExpMatrix.FPKM  <- geneExpMatrix.FPKM[geneExpMatrix.FPKM$gene_type=="lncRNA",]
geneNames.FPKM <- geneExpMatrix.FPKM$gene_id
lncRNAExpMatrix.FPKM  <-  lncRNAExpMatrix.FPKM[,-c(1:3)] %>% data.matrix()
rownames(lncRNAExpMatrix.FPKM ) <- as.vector(geneNames.FPKM)
lncRNAExpMatrix.FPKM  <- avereps(lncRNAExpMatrix.FPKM ) %>% as.data.frame()
dim(lncRNAExpMatrix.FPKM )

lncRNAExpMatrix.count  <- geneExpMatrix.count[geneExpMatrix.count$gene_type=="lncRNA",]
lncRNANames.count <- lncRNAExpMatrix.count$gene_id
lncRNAExpMatrix.count <-  lncRNAExpMatrix.count[,-c(1:3)] %>% data.matrix()
rownames(lncRNAExpMatrix.count) <- as.vector(lncRNANames.count)
lncRNAExpMatrix.count <- avereps(lncRNAExpMatrix.count) %>% as.data.frame()
dim(lncRNAExpMatrix.count)
###
# 从数据集中删除零方差的基因行，可以使用相同的apply 表达式，设置方差不等于零。
lncRNAExpMatrix.FPKM  <- lncRNAExpMatrix.FPKM [which(apply(lncRNAExpMatrix.FPKM , 1, var)!=0),]
dim(lncRNAExpMatrix.FPKM )
# 从数据集中删除零方差的基因行，可以使用相同的apply 表达式，设置方差不等于零。
lncRNAExpMatrix.count  <- lncRNAExpMatrix.count [which(apply(lncRNAExpMatrix.count , 1, var)!=0),]
dim(lncRNAExpMatrix.count )

##保存数据
##记录所有样本ID,这个非常重要，后边还会用到此功能，特别是今后在用101.1 脚本时候注意保存ID
Barcode <- intersect(colnames(lncRNAExpMatrix.FPKM),colnames(lncRNAExpMatrix.count))
sample.info <- getSampleInfo(Barcode)##Barcode
#保存文件
write.table(lncRNAExpMatrix.FPKM ,file = paste(outputPath[1],"lncRNAExpMatrix.FPKM.txt",sep = ""),
            sep = "\t",row.names = T)
#保存文件
write.table(lncRNAExpMatrix.count ,file = paste(outputPath[1],"lncRNAExpMatrix.count.txt",sep = ""),
            sep = "\t",row.names = T)
save(lncRNAExpMatrix.FPKM ,sample.info,file = paste0(outputPath[1],"COAD.lncRNAExpMatrix.FPKM.Rdata"))
save(lncRNAExpMatrix.count,sample.info,file = paste0(outputPath[1],"COAD.lncRNAExpMatrix.count.Rdata"))



######################################，此步骤其实m6A用得到，因为需要m6A gene与矩阵取交集，
# 而lncRNA直接在上边获取，不需要下边的步骤
options(stringsAsFactors = F)
library(limma)
#####################
#开始读入要分析的基因数据，比如免疫相关基因，能量代谢相关基因等。
gene <- readLines("02-inputFiles/Genes.txt")#最后一行为空行
# 取表达矩阵中的基因与我们要分析的基因的交集
yesgene <-intersect(gene,as.vector(row.names(proteinExpMatrix.count)))##proteinExpMatrix.FPKM也行
# 表达矩阵中没有的基因
nogene <- setdiff(gene,yesgene)

writeLines(c(paste("总共输入了",length(gene),"基因：",sep=""),
             paste("在表达矩阵中有的基因个数：",length(yesgene),sep=""),
             yesgene,paste("在表达矩阵中没有的基因个数：",length(nogene),sep=""),nogene),
           con=paste(outputPath[1],"filespec.txt",sep = ""),sep="\n")

####保存用于本次分析的表达数据,分两大部分查找m6A和LncRNA的各自矩阵
####m6A
#输入基因所有样本的表达矩阵
inputGeneAllSampleExpMatrix.count <- proteinExpMatrix.count[yesgene,]
write.table(inputGeneAllSampleExpMatrix.count,
            file=paste(outputPath[1],"inputGeneAllSampleExpMatrix.count.txt",sep = ""),
            sep="\t",row.names=T,quote=F)
#输入基因删除正常样本的表达矩阵
inputGeneOnlyTumorExpMatrix.count <- inputGeneAllSampleExpMatrix.count[,sample.info[["TumorBarcode"]]]
write.table(inputGeneOnlyTumorExpMatrix.count,
            file=paste0(outputPath[1],"inputGeneOnlyTumorExpMatrix.count.txt"),
            sep="\t",row.names=T,quote=F)
save(inputGeneAllSampleExpMatrix.count,file = paste0(outputPath[1],"COAD.m6AgeneExpMatrix.count.Rdata"))
save(inputGeneOnlyTumorExpMatrix.count,file = paste0(outputPath[1],"COAD.m6AgeneOnlyTumorExpMatrix.count.Rdata"))
####上面是count，下面是FPKM
inputGeneAllSampleExpMatrix.FPKM <- proteinExpMatrix.FPKM[yesgene,]
write.table(inputGeneAllSampleExpMatrix.FPKM,
            file=paste(outputPath[1],"inputGeneAllSampleExpMatrix.FPKM.txt",sep = ""),
            sep="\t",row.names=T,quote=F)
inputGeneOnlyTumorExpMatrix.FPKM <- inputGeneAllSampleExpMatrix.FPKM[,sample.info[["TumorBarcode"]]]
write.table(inputGeneOnlyTumorExpMatrix.FPKM,
            file=paste0(outputPath[1],"inputGeneOnlyTumorExpMatrix.FPKM.txt"),
            sep="\t",row.names=T,quote=F)
save(inputGeneOnlyTumorExpMatrix.FPKM,file = paste0(outputPath[1],"COAD.m6AgeneOnlyTumorExpMatrix.FPKM.Rdata"))
save(inputGeneOnlyTumorExpMatrix.FPKM,file = paste0(outputPath[1],"COAD.m6AgeneOnlyTumorExpMatrix.FPKM.Rdata"))

####LncRNA的
lncRNAAllSampleExpMatrix.count <- lncRNAExpMatrix.count
write.table(lncRNAAllSampleExpMatrix.count,
            file=paste(outputPath[1],"lncRNAAllSampleExpMatrix.count.txt",sep = ""),
            sep="\t",row.names=T,quote=F)
#输入基因删除正常样本的表达矩阵
lncRNAOnlyTumorExpMatrix.count <- lncRNAAllSampleExpMatrix.count[,sample.info[["TumorBarcode"]]]
write.table(lncRNAOnlyTumorExpMatrix.count,
            file=paste0(outputPath[1],"lncRNAOnlyTumorExpMatrix.count.txt"),
            sep="\t",row.names=T,quote=F)
save(lncRNAAllSampleExpMatrix.count,file = paste0(outputPath[1],"COAD.lncRNAALLExpMatrix.count.Rdata"))
save(lncRNAOnlyTumorExpMatrix.count,file = paste0(outputPath[1],"COAD.lncRNAOnlyTumorExpMatrix.count.Rdata"))
####上面是count，下面是FPKM
load(file = paste0(outputPath[1],"COAD.lncRNAExpMatrix.FPKM.Rdata"))
lncRNAAllSampleExpMatrix.FPKM <- lncRNAExpMatrix.FPKM
write.table(lncRNAAllSampleExpMatrix.FPKM,
            file=paste(outputPath[1],"lncRNAAllSampleExpMatrix.FPKM.txt",sep = ""),
            sep="\t",row.names=T,quote=F)
save(lncRNAAllSampleExpMatrix.FPKM,file = paste0(outputPath[1],"COAD.lncRNAAllSampleExpMatrix.FPKM.Rdata"))
lncRNAOnlyTumorExpMatrix.FPKM <- lncRNAAllSampleExpMatrix.FPKM[,sample.info[["TumorBarcode"]]]
write.table(lncRNAOnlyTumorExpMatrix.FPKM,
            file=paste0(outputPath[1],"lncRNAOnlyTumorExpMatrix.FPKM.txt"),
            sep="\t",row.names=T,quote=F)
save(lncRNAOnlyTumorExpMatrix.FPKM,file = paste0(outputPath[1],"COAD.lncRNAOnlyTumorExpMatrix.FPKM.Rdata"))

####上边最后得到lncRNA相关和m6A的CRC表达矩阵，刘元升，2021-07-08

# here, i used 101.1 R script :'m6A-lncRNA correlation analysis' finally, got the m6A relatied lncRNA of CRC matrix, for following analysis


# then, differential analysis
#####--------------------------102----------------------------------------
#差异分析

library("limma")
library( "edgeR" )
load(paste0(outputPath[1],"COAD.m6AlncRNAAllSampleExpMatrix.count.Rdata"))
dim(m6AlncRNAAllSampleExpMatrix.countexp)
##前面已经去除了0方差的基因，这里可以检查一下
proExpMat.count <- m6AlncRNAAllSampleExpMatrix.countexp[which(apply(m6AlncRNAAllSampleExpMatrix.countexp, 1, var)!=0),]
dim(proExpMat.count)
# 创建分组
Barcode <- colnames(proExpMat.count)
sample.info <- getSampleInfo(Barcode)##Barcode
sample.info$BarcodeSort <- grep("^TCGA",sample.info$BarcodeSort,value = T)  #去掉不知名的NA样本名
BarcodeSort <- sample.info$BarcodeSort ##正常样本ID在前，肿瘤样本ID在后
proExpMat.count <- proExpMat.count[,BarcodeSort]
conNum <- sample.info[["sampleNum"]][["NormalBarcode"]]#control组样品数目
treatNum <- sample.info[["sampleNum"]][["TumorBarcode"]]#tumor组样品数目
df <- as.matrix(proExpMat.count)
df <- apply(df,2,as.numeric)
rownames(df) <- rownames(proExpMat.count)
group <- c(rep(1,conNum),rep(2,treatNum))#1表示对照，2表示测试组(这里：肿瘤)
y <- DGEList(counts=df,group=group)
# 数据过滤
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
# 计算标准化因子
y <- calcNormFactors(y)
# 查看标准化因子
y$samples
# 计算离散度
y <- estimateDisp(y)
# 显著性检验
T2C.DiffExpGene <- exactTest(y)
# 获取排名靠前的基因，这里设置n=100000是为了输出所有基因
T2C.DiffExpGene <- topTags(T2C.DiffExpGene, n=100000)
# 转换为数据框类型
T2C.DiffExpGene <- as.data.frame(T2C.DiffExpGene)
# 将行名粘贴为数据框的第一列
T2C.DiffExpGene <- cbind(rownames(T2C.DiffExpGene),T2C.DiffExpGene)
# 指定列名
colnames(T2C.DiffExpGene) <- c("gene_id", "log2FoldChange", "log2CPM", "PValue", "FDR")

# 获取具有显著差异表达的结果
T2C.DiffExpGene$logP <- -log10(T2C.DiffExpGene$FDR)
sum(is.na(rownames(T2C.DiffExpGene)))
# 保存数据到本地
write.table(T2C.DiffExpGene, paste0(outputPath[2],"COAD_ALL_DEG.xls"), sep = "\t", col.names = TRUE,
            row.names = FALSE, quote = FALSE, na = "")

###绘制火山图
####可视化
library(TCGAbiolinks)
volFig <- pdf(paste0(outputPath[2],"COAD_ALL_DEG_vol.pdf"),width=10,height=6)
TCGAVisualize_volcano(T2C.DiffExpGene$log2FoldChange, T2C.DiffExpGene$FDR,
                      filename = volFig, xlab = "logFC",
                      names = rownames(T2C.DiffExpGene), show.names = "highlighted",
                      x.cut = 1, y.cut = 0.01, 
                      highlight = rownames(T2C.DiffExpGene)[which(abs(T2C.DiffExpGene$log2FoldChange) >=8)],
                      highlight.color = "orange")
dev.off()
###############

#####热图的绘制
#install.packages("pheatmap")
# load(paste0(outputPath[1],"COAD.proteinExpMatrix.FPKM.Rdata"))
# hdata <- proteinExpMatrix.FPKM[,BarcodeSort]##正常样本ID在前，肿瘤样本ID在后
hdata <- df
sum(is.na(rownames(hdata)))
yesgeneExp <- log2(hdata+1)
yesgeneDiff <- T2C.DiffExpGene
yesgeneDiff <- na.omit(yesgeneDiff)
rownames(yesgeneDiff)
write.table(yesgeneDiff,file = paste0(outputPath[2],"yesgeneDiff.txt"),sep = "\t")
##筛选差异基因条件根据结果数量决定
sig.gene <- yesgeneDiff$gene_id[abs(yesgeneDiff$log2FoldChange)>2 & yesgeneDiff$FDR <0.01]

volFig <- pdf(paste0(outputPath[2],"COAD_sigyesgeneDiff_vol.pdf"),width=10,height=6)
TCGAVisualize_volcano(T2C.DiffExpGene$log2FoldChange, T2C.DiffExpGene$FDR,
                      filename = volFig, xlab = "logFC",
                      names = rownames(T2C.DiffExpGene), show.names = "highlighted",
                      x.cut = 1, y.cut = 0.01, 
                      highlight = sig.gene,
                      highlight.color = "orange")
dev.off()


diffGeneExpHeatMapData <- hdata[sig.gene,]
diffGeneExpHeatMapData <- 4*diffGeneExpHeatMapData/(max(diffGeneExpHeatMapData)-min(diffGeneExpHeatMapData))-2

library(pheatmap)
Type <- c(rep("N",conNum),rep("T",treatNum))
names(Type)=colnames(diffGeneExpHeatMapData)
Type <- as.data.frame(Type)


pdf(paste(outputPath[2],"diffgeneHeatMap.pdf",sep = ""),height=10,width=20)
pheatmap(diffGeneExpHeatMapData, 
         annotation=Type, 
         color = colorRampPalette(c("#0509F2", "white", "#F00528"))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 10,
         fontsize_row=10,
         fontsize_col=3)
dev.off()
# 保存一下方便后边哟弄个到这个m6A-lncRNA矩阵
save(sample.info,yesgeneDiff,yesgeneExp,proExpMat.count,file = paste0(outputPath[1],"proExpMat.count.Rdata"))
###############------------------103--------------------

#################################103
#install.packages("corrplot")
# load(paste0(outputPath[1],"COAD.proteinExpMatrix.FPKM.Rdata"))

###删除正常样本
tur.yg.FPKM <- proExpMat.count[,sample.info[["TumorBarcode"]]]
tur.yg.FPKM <- apply(tur.yg.FPKM, 2, as.numeric)
tur.yg.FPKM <- log2(tur.yg.FPKM+1)
rownames(tur.yg.FPKM) <- rownames(proExpMat.count)
library(corrplot)
library(ggcorrplot)
cordata <- t(tur.yg.FPKM)
res1 <- cor_pmat(cordata)
corr <- cor(cordata,method = "spearman")
write.table(corr,file = paste0(outputPath[3],"corr.xls"),sep = "\t")
##绘图，基因很多，只做sig.gene
cd <- cordata[,sig.gene]
res2 <- cor_pmat(cd)
corr2 <- cor(cd,method = "spearman")
###根据基因数量调整图片大小
pdf(paste(outputPath[3],"corHeatFig.pdf",sep = ""),height=20,width=20)#保存图片的文件名称
ggcorrplot(corr2,method = "circle",hc.order = T,hc.method = "ward.D",
           outline.color = "white",ggtheme = theme_bw(),type = "up",
           p.mat = res2,insig = "blank",lab = T,lab_size = 2)
dev.off()


########################-------------104-------------------

#######################肿瘤分型######################
# load(paste0(outputPath[1],"COAD.proteinExpMatrix.FPKM.Rdata"))
###删除正常样本
tur.yg.FPKM <- proExpMat.count[,sample.info[["TumorBarcode"]]]   #此处要是只筛m6A或者lncRNA,则需要行筛选处加上yesgene


library(limma)

#聚类
maxK=9
library(ConsensusClusterPlus)
results <- ConsensusClusterPlus(as.matrix(tur.yg.FPKM),
                                maxK=maxK,
                                reps=50,
                                pItem=0.8,
                                pFeature=1,
                                title=outputPath[4],
                                clusterAlg="km",
                                distance="euclidean",
                                seed=123456,
                                plot="png")

#输出结果,根据结果可看出选择K=3最合适（CDF曲线下面积图在3处有拐点）
clusterNum <- 3  #分几类，根据判断标准判断
cluster <- results[[clusterNum]][["consensusClass"]]
write.table(cluster,file=paste(outputPath[4],"Clouster.txt",sep = ""),
            sep="\t",quote=F,col.names=F)

#########--------------------105-------

##############主成分分析#############
# load(paste0(outputPath[1],"COAD.proteinExpMatrix.FPKM.Rdata"))
###删除正常样本
tur.yg.FPKM <- proExpMat.count[,sample.info[["TumorBarcode"]]]

#pca analysis
### PCA分析
pcaData <- t(as.matrix(tur.yg.FPKM))
pcaData.class <- rownames(pcaData)
pcaData <- apply(pcaData, 2, as.numeric)
rownames(pcaData) <- colnames(tur.yg.FPKM)
pcaData.pca <- prcomp(pcaData, scale. = TRUE)                                  #PCA分析
write.table(predict(pcaData.pca),file=paste(outputPath[5],"PCAResult.txt",sep = "")
            ,quote=F,sep="\t")        #输出新表

#可视化
library(ggplot2)
clusterInfo <- as.data.frame(cluster)
group <- paste0("cluster",as.vector(clusterInfo[,1]))
pcaPredict <- predict(pcaData.pca)
PCA <- data.frame(PCA1 = pcaPredict[,1], PCA2 = pcaPredict[,2],group=group)
pdf(file=paste(outputPath[5],"CloustePCAFig.pdf",sep = ""),height=10,width=13)#保存输入出文件
ggplot(data = PCA, aes(PCA1, PCA2)) + geom_point(aes(color = group)) +
  theme_bw()+
  theme(plot.margin=unit(rep(1.5,4),'lines'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()
##################----------------106-----------
###################################106

options(stringsAsFactors = F)
ClinicalData <- as.data.frame(data.table::fread(ClinicalDataPath,))
ClinicalData <- ClinicalData[!duplicated(ClinicalData$case_submitter_id),]
names(ClinicalData)
pClinData <- ClinicalData[,c("case_submitter_id",
                             "gender",
                             "days_to_death",
                             "days_to_last_follow_up",
                             "vital_status",
                             "ajcc_pathologic_stage",
                             "ajcc_pathologic_t",
                             "ajcc_pathologic_n",
                             "ajcc_pathologic_m")]

surData <- getSurData(ClinicalData)
clindata <- merge(surData,pClinData,by = "case_submitter_id")
clindata$ajcc_pathologic_stage <- gsub("[A-D].*","",c(pClinData$ajcc_pathologic_stage))
clindata$ajcc_pathologic_t <- gsub("[a-c].*","",c(clindata$ajcc_pathologic_t))
clindata$ajcc_pathologic_n <- gsub("[a-c].*","",c(clindata$ajcc_pathologic_n))
clindata$ajcc_pathologic_m <- gsub("[a-c].*","",c(clindata$ajcc_pathologic_m))
clindata <- clindata[,c("case_submitter_id","SurvivalTime","VitalStatus",
                        "ajcc_pathologic_stage",
                        "ajcc_pathologic_t",
                        "ajcc_pathologic_n",
                        "ajcc_pathologic_m")]

colnames(clindata) <- c("Barcode","SurvivalTime","VitalStatus","Stage","T","N","M")
write.table(clindata,file=paste(outputPath[6],"clindata.txt",sep = ""),sep="\t",row.names=F,quote=F)
rio::export(clindata,file = paste(outputPath[6],"clindata.csv"))
save(clindata,file = paste(outputPath[6],"clindata.Rdata"))
#########
##ultiClinDataFN <- "ultiClinData.csv"
# 主要是查看一下2种方式下载的数据是否有区别，就生存时间会有区别

###去掉重复的数据

######上面代码产生的数据不用于分析
# load(paste0(outputPath[1],"COAD.proteinExpMatrix.FPKM.Rdata"))
# 此处还是用m6A-lncRNA的数据
###删除正常样本
tur.FPKM <- proExpMat.count[,sample.info[["TumorBarcode"]]]

##重新定义列名，使得ID与临床数据ID匹配
colnames(tur.FPKM) <- gsub("(.*?)-(.*?)-(.*?)-(.*?)-.*","\\1-\\2-\\3",colnames(tur.FPKM))
tur.FPKM <- rbind(Barcode=colnames(tur.FPKM),tur.FPKM)
tur.FPKM <- t(tur.FPKM)
tur.clin.FPKM <- merge(clindata,tur.FPKM,by ="Barcode")

###保存数据
save(tur.clin.FPKM,file = paste0(outputPath[6],"COAD.tur.clin.FPKM.Rdata"))

######---------------------107------------------------
##########107-tumorClusterSurvivalAnalysis
options(stringsAsFactors = F)

load(paste0(outputPath[6],"COAD.tur.clin.FPKM.Rdata"))
SurD <- tur.clin.FPKM

clusterNum = clusterNum #聚类分为几类 和004相关

classification <- as.data.frame(cluster)
classification$Barcode <- rownames(classification)
classification$Barcode <- gsub("(.*?)-(.*?)-(.*?)-(.*?)-.*","\\1-\\2-\\3",classification$Barcode)
classification <- classification[,c("Barcode","cluster")]


cluster.sur <- merge(SurD,classification,by = "Barcode")
###因为一个病人可能测序就几个数据，所以前3端ID可能有重复的，可通过第四个字段区别
cluster.sur <- cluster.sur[!duplicated(cluster.sur$Barcode),]
rownames(cluster.sur) <- cluster.sur$Barcode
cluster.sur$cluster <- paste0("Cluster ",cluster.sur$cluster)
cluster.sur$SurvivalTime <- as.numeric(cluster.sur$SurvivalTime)
cluster.sur <- na.omit(cluster.sur) 
library(dplyr)


cluster.sur <- select(cluster.sur, Barcode, SurvivalTime,VitalStatus, 
                      Stage,`T`,N,M,cluster,everything())
#保存数据
save(cluster.sur,file = paste0(outputPath[7],"COAD-TNM-Sur.Rdata"))

#install.packages("survival")
library(survival)
library(survminer)
diff <- survdiff(Surv(SurvivalTime, VitalStatus) ~ cluster,data = cluster.sur)

pValue <- 1-pchisq(diff$chisq,df=clusterNum-1)# df为自由度n-1
if(pValue<0.001){
  pValue <- signif(pValue,4)
  pValue <- format(pValue, scientific = TRUE)
}else{
  pValue <- round(pValue,3)
}

fit <- survfit(Surv(SurvivalTime, VitalStatus) ~ cluster, data = cluster.sur)
# summary(fit)# 查看5年生存率

pdf(file=paste(outputPath[7],"ClousterSurCurve.pdf",sep = ""),width = 10,height =10)
ggsurvplot(fit,palette = c("#FF008C","#0900FF","#00FF5E","#EAFF00"),
           xlab ="Day", ylim=c(0,1),
           risk.table = T,
           conf.int = TRUE,
           surv.median.line = "hv"
)
dev.off()

##############---------------------108------------------
##############################108
options(stringsAsFactors = F)

load(paste0(outputPath[7],"COAD-TNM-Sur.Rdata"))##来自107保存的数据

SurD <- cluster.sur
######删除缺失数据，事先看一下缺失数据是什么表示的，一般是"'--",但这里是"'--"
SurD <- SurD[which(SurD$Stage != "'--"),]
SurD <- SurD[which(SurD$T != "'--"),]
SurD <- SurD[which(SurD$N != "'--"),]
SurD <- SurD[which(SurD$N != "NX"),]
SurD <- SurD[which(SurD$M != "'--"),]

###如果MX少的话，可以删除。多的话不删
if(length(grep("MX",SurD$M)) < nrow(SurD)/10){
  SurD <- SurD[which(SurD$M != "MX"),]
}
indepData <- SurD[,1:7]
save(indepData,file = paste0(outputPath[108],"indepData.Rdata"))

### 临床相关性检验：分型与临床症状的相关性检验
library(dplyr)
#SurD <- select(SurD,cluster,everything())
colnames(SurD)[1:8]
##自己先查看一下列名，注意顺序
newLabels <- c("Barcode","SurvivalTime","VitalStatus")
for(i in 4:7){
  nameStat <- colnames(SurD)[i]
  tableStat <- table(SurD[,c(nameStat,"cluster")])
  pStat <- chisq.test(tableStat) # 注意检验方法
  pvalue <- pStat$p.value
  if(pvalue<0.001){
    newLabels <- c(newLabels,paste0(colnames(SurD)[i],"***"))
  }else if(pvalue<0.01){
    newLabels <- c(newLabels,paste0(colnames(SurD)[i],"**"))
  }else if(pvalue<0.05){
    newLabels <- c(newLabels,paste0(colnames(SurD)[i],"*"))
  }else{
    newLabels <- c(newLabels,colnames(SurD)[i])
  }
  print(paste(colnames(SurD)[i],pvalue,sep=" "))
}
colnames(SurD) <- c(newLabels,"cluster",colnames(SurD)[9:ncol(SurD)])

######################################################
#install.packages("pheatmap")
library(pheatmap)
SurD <- SurD[order(SurD$cluster),]
rownames(SurD) <- SurD$Barcode
st <- SurD[,c(3:8)]
sum(is.na(st$VitalStatus))
st$VitalStatus <- ifelse(st$VitalStatus ==0,"Alive","Dead")
ExpMatrix <- t(SurD[,sig.gene]) %>% as.data.frame()#
ExpMatrix <- data.matrix(ExpMatrix) %>% as.data.frame()
ExpMatrix <- log2(ExpMatrix+1)

pdf(paste(outputPath[8],"ClousterWithClinCor.pdf",sep = ""),height=8,width=10)
pheatmap(ExpMatrix, annotation = st, 
         color = colorRampPalette(c("#0509F2", "white", "#F00528"))(50),
         cluster_cols =F,
         fontsize=8,
         fontsize_row=8,
         scale="row",
         show_colnames=F,
         fontsize_col=3)
dev.off()
###############-----------------109------------------
##########109
options(stringsAsFactors = F)
##############################################################
### 读入数据并整理
load(paste0(outputPath[7],"COAD-TNM-Sur.Rdata"))##来自107保存的数据

coxSur.fpkm <- cluster.sur[,c("SurvivalTime","VitalStatus",sig.gene)]
coxSur.fpkm <- data.matrix(coxSur.fpkm) %>% as.data.frame()


library(survival)
library(forestplot)
options(forestplot_new_page = FALSE) # 在pdf文件中生成图形就在第一页，不声明会在第二页
clrs <- fpColors(box="green",line="darkblue", summary="royalblue")   #定义森林图颜色

coxResult <- data.frame()
for(i in colnames(coxSur.fpkm[,3:ncol(coxSur.fpkm)])){
  cox <- coxph(Surv(SurvivalTime,VitalStatus) ~ coxSur.fpkm[,i], data = coxSur.fpkm)
  coxSummary <- summary(cox)
  coxP <- coxSummary$coefficients[,"Pr(>|z|)"]
  coxResult <- rbind(coxResult,
                     cbind(Gene=i,
                           HR=coxSummary$conf.int[,"exp(coef)"],
                           HR.95L=coxSummary$conf.int[,"lower .95"],
                           HR.95H=coxSummary$conf.int[,"upper .95"],
                           pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}

write.table(coxResult,file =paste(outputPath[9],"coxResult.txt",sep = ""),sep = "\t")
write.table(coxResult,file=paste(outputPath[9],"coxResult.xls",sep = ""),sep="\t",row.names=F,quote=F)
#绘制森林图
rt<-read.table(paste(outputPath[9],"coxResult.xls",sep = ""),
               header=T,sep="\t",row.names=1,check.names=F)
data <- as.matrix(rt)
HR <- data[,1:3]
hr <- sprintf("%.3f",HR[,"HR"])
hrLow <- sprintf("%.3f",HR[,"HR.95L"])
hrHigh <- sprintf("%.3f",HR[,"HR.95H"])
pVal <- data[,"pvalue"]
pVal <- ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
tabletext <- 
  list(c("m6A-LncRNA", rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )         
pdf(file=paste(outputPath[9],"forest.pdf",sep = ""),
    width = 8,             
    height = 10,           
)
forestplot(tabletext, 
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(50, "mm"),
           xlog=T,
           lwd.ci=2,
           boxsize=0.4,
           xlab="Hazard ratio"
)
dev.off()
############-------------110------------
options(stringsAsFactors = F)
########################110-LassoRegression
### 读入数据并整理,也可以用COX回归的中合并的文件

###############################下面代码还没有修改
#install.packages("glmnet")
library("glmnet") # 
library("survival")
coxGene <-coxResult[which(coxResult$pvalue < 0.05),]$Gene
lassoData <- coxSur.fpkm[,c("SurvivalTime","VitalStatus",coxGene)]
lassoData <- lassoData[lassoData$SurvivalTime >0,]

x <- data.matrix(lassoData[,c(3:ncol(lassoData))])# 
y <- data.matrix(Surv(lassoData$SurvivalTime,lassoData$VitalStatus))

####################################################

#################################################
fit <- glmnet(x, y,family = "cox" ,maxit = 5000)# 

pdf(paste(outputPath[10],"lambda.pdf",sep = ""))
plot(fit, xvar = "lambda", label = TRUE)
dev.off()

cvfit <- cv.glmnet(x, y, family = "cox", maxit = 10000)
pdf(paste0(outputPath[10],"cvfit.pdf"))
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
dev.off()

coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene <- row.names(coef)[index]# lassoGene=row.names(coef)[index]

geneCoef <- cbind(Gene=lassoGene,Coef=actCoef)

write.table(geneCoef,file=paste0(outputPath[10],"geneCoef.txt"),sep="\t",quote=F,row.names=F)

riskScore <- predict(cvfit, newx = x, s = "lambda.min",type="response")
outCol <- c("SurvivalTime","VitalStatus",lassoGene)
risk <- as.vector(ifelse(riskScore>median(riskScore),"high","low"))
outTab <- cbind(lassoData[,outCol],riskScore=as.vector(riskScore),risk)
write.table(cbind(Barcode = rownames(outTab),outTab),
            file=paste0(outputPath[10],"lassoRisk.txt"),
            sep="\t",
            quote=F,
            row.names=F)

#########--------------111---------------
#################111
#install.packages("survival")

library(survival)
lr <- read.table(paste0(outputPath[10],"lassoRisk.txt"),header=T,sep="\t")
diff <- survdiff(Surv(SurvivalTime, VitalStatus) ~risk,data = lr)
pValue <- 1-pchisq(diff$chisq,df=1)
pValue <- round(pValue,3)
pValue <- signif(pValue,4)
pValue <- format(pValue, scientific = TRUE)

fit <- survfit(Surv(SurvivalTime, VitalStatus) ~ risk, data = lr)
summary(fit)    #查看五年生存率
pdf(file=paste0(outputPath[11],"RiskSurvival.pdf"),width=5.5,height=6)
ggsurvplot(fit,palette = c("#FF008C","#0900FF","#00FF5E","#EAFF00"),
           xlab ="Day", ylim=c(0,1),
           risk.table = T)
dev.off()
############-----------------------112----------------

###112
#install.packages("survival")
library(survivalROC)
lr <- read.table(paste0(outputPath[10],"lassoRisk.txt"),header=T,sep="\t",
                 check.names=F,row.names=1)    #读取lasso回归风险文件
lr$SurvivalTime <- lr$SurvivalTime/365
pdf(file=paste0(outputPath[12],"ROC.pdf"),width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc <- survivalROC(Stime=lr$SurvivalTime, status=lr$VitalStatus, marker = lr$riskScore, 
                   predict.time =5, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="False positive rate", ylab="True positive rate",
     main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
dev.off()
##############
#####---------------113——————————————————————————
###############113
#######数据整理
##需要一风险基因的表达矩阵，和一个临床数据
lassoRisk <- read.table(paste0(outputPath[10],"lassoRisk.txt"),header=T,sep="\t",
                        check.names=F)    #读取lasso回归风险文件
lassoRiskGroup <- lassoRisk[,c("Barcode","risk")]
#临床表型数据文件
load(paste0(outputPath[7],"COAD-TNM-Sur.Rdata"))
cs <- cluster.sur
lrcs <- merge(lassoRiskGroup,cs,by="Barcode")
names(lrcs)[1:9]
######显著性检验
lrcs <- select(lrcs,Barcode,SurvivalTime,VitalStatus,Stage,`T`,N,M,cluster,risk,everything())
lrcs <- arrange(lrcs,risk)
lrcs$VitalStatus <- ifelse(lrcs$VitalStatus==0,"Alive","Dead")

names(lrcs)[1:9]

newcoln = names(lrcs)[1:2]
for(i in names(lrcs)[3:8]){
  tableStat=table(lrcs[,c(i,"risk")])
  pStat=chisq.test(tableStat)# 注意检验方法
  pvalue=pStat$p.value
  if(pvalue<0.001){
    newcoln=c(newcoln,paste0(i,"***"))
  }else if(pvalue<0.01){
    newcoln=c(newcoln,paste0(i,"**"))
  }else if(pvalue<0.05){
    newcoln=c(newcoln,paste0(i,"*"))
  }else{
    newcoln=c(newcoln,i)
  }
  print(paste(i,pvalue,sep=" "))
}
newcoln <- c(newcoln,"risk")
colnames(lrcs)=c(newcoln,colnames(lrcs)[10:length(colnames(lrcs))])
rownames(lrcs)<-lrcs$Barcode
write.table(lrcs,file=paste0(outputPath[13],"RiskWtihClinSignificanceTest.txt"),
            sep="\t",row.names=F,quote=F)
##############################绘制热图######
Type <- lrcs[,3:9]

ExpMatrix2 <- t(lrcs[,sig.gene]) %>% as.data.frame()#
ExpMatrix2 <- data.matrix(ExpMatrix2) %>% as.data.frame()
ExpMatrix2 <- log2(ExpMatrix2+1)

library(pheatmap)
pdf(paste0(outputPath[13],"heatmap.pdf"),height=6,width=10)
pheatmap(ExpMatrix2, annotation=Type, 
         color = colorRampPalette(c("#0509F2", "white", "#F00528"))(50),
         cluster_cols =F,
         fontsize=7.5,
         fontsize_row=8,
         scale="row",
         show_colnames=F,
         fontsize_col=3)
dev.off()

###coxGene
ExpMatrix3 <- t(lrcs[,coxGene]) %>% as.data.frame()#
ExpMatrix3 <- data.matrix(ExpMatrix3) %>% as.data.frame()
ExpMatrix3 <- log2(ExpMatrix3 + 1)

library(pheatmap)
pdf(paste0(outputPath[13],"inputGene-heatmap.pdf"),height=6,width=10)
pheatmap(ExpMatrix3, annotation=Type, 
         color = colorRampPalette(c("#0509F2", "white", "#F00528"))(50),
         cluster_cols =F,
         fontsize=7.5,
         fontsize_row=8,
         scale="row",
         show_colnames=F,
         fontsize_col=3)
dev.off()
###---------------------114-----------------
###################114-IndependentPrognostic

#读取lasso回归风险文件,lassoRisk.txt文件中 
lassoRisk <- read.table(paste0(outputPath[10],"lassoRisk.txt"),header=T,sep="\t",
                        check.names=F)    #读取lasso回归风险文件  
load(paste0(outputPath[108],"indepData.Rdata"))
ilr <- merge(indepData,lassoRisk[,c("Barcode","riskScore")],by="Barcode")
names(ilr)
ilr$Stage <- as.numeric(as.factor(ilr$Stage))
ilr$T <- as.numeric(as.factor(ilr$T))
ilr$N <- as.numeric(as.factor(ilr$N))
ilr$M <- as.numeric(as.factor(ilr$M))
row.names(ilr) <- ilr$Barcode
#####写出数据用
write.csv(ilr,file = paste0(outputPath[14],"IndependentPrognosticData.csv"),row.names = F)
ilr <- ilr[,-1]

####################

options(stringsAsFactors = F)
#install.packages('survival')
#install.packages('forestplot')
library(survival)
library(forestplot)
options(forestplot_new_page = FALSE)
clrs <- fpColors(box="green",line="darkblue", summary="royalblue") #定义森林图颜色

names(ilr)
outTab=data.frame()
for(i in colnames(ilr[,3:ncol(ilr)])){
  cox <- coxph(Surv(SurvivalTime, VitalStatus) ~ ilr[,i], data = ilr)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(id=i,
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
write.table(outTab,file=paste0(outputPath[14],"uniCox.xls"),
            sep="\t",row.names=F,quote=F)

#绘制森林图
uniCox=read.table(paste0(outputPath[14],"uniCox.xls"),header=T,sep="\t",row.names=1,check.names=F)
data=as.matrix(uniCox)
HR=data[,1:3]
hr=sprintf("%.3f",HR[,"HR"])
hrLow=sprintf("%.3f",HR[,"HR.95L"])
hrHigh=sprintf("%.3f",HR[,"HR.95H"])
pVal=data[,"pvalue"]
pVal=ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
tabletext <- 
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) ) #定义图片文字
pdf(file=paste0(outputPath[14],"singleFactorForestFig.pdf"),
    width = 6,             #图片的宽度
    height = 4,            #图片的高度
)
forestplot(tabletext, 
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(50, "mm"),
           xlog=T,
           lwd.ci=2,
           boxsize=0.3,
           xlab="Hazard ratio"
)
dev.off()
#######
#install.packages('survival')
#install.packages('forestplot')

library(survival)
library(forestplot)
options(forestplot_new_page = FALSE)
clrs <- fpColors(box="red",line="darkblue", summary="royalblue")  #定义森林图颜色

multiCox <- coxph(Surv(SurvivalTime, VitalStatus) ~ ., data = ilr)# 多因素就变成点啦
multiCoxSum <- summary(multiCox)

outTab <- data.frame()
outTab <- cbind(
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])

write.table(outTab,file=paste0(outputPath[14],"multiCox.xls")
            ,sep="\t",row.names=T,quote=F)

#绘制森林图
data <- as.matrix(outTab)
HR <- data[,1:3]
hr <- sprintf("%.3f",HR[,"HR"])
hrLow <- sprintf("%.3f",HR[,"HR.95L"])
hrHigh <- sprintf("%.3f",HR[,"HR.95H"])
pVal <- data[,"pvalue"]
pVal <- ifelse(pVal<0.001, "<0.001", sprintf("%.3f", pVal))
tabletext <- 
  list(c(NA, rownames(HR)),
       append("pvalue", pVal),
       append("Hazard ratio",paste0(hr,"(",hrLow,"-",hrHigh,")")) )  #定义图片文字
pdf(file=paste0(outputPath[14],"multiCoxforest.pdf"),
    width = 6,             #图片的宽度
    height = 4,            #图片的高度
)
forestplot(tabletext, 
           rbind(rep(NA, 3), HR),
           col=clrs,
           graphwidth=unit(50, "mm"),
           xlog=T,
           lwd.ci=2,
           boxsize=0.3,
           xlab="Hazard ratio"
)
dev.off()


