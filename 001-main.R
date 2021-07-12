options(stringsAsFactors = F)
######----------------------1.定义相关文件路径-------------------------------
# 需要输入的文件，注意修改
dir.create("01-RawData")
dir.create("02-inputFiles")
#gene文件是我们要分析的基因，比如m6A,免疫调控相关的基因。gtf文件用于注释
inputPath <- ".\\02-inputFiles" #gene文件和gtf文件路径
RNASeqDataPath <- ".\\01-RawData\\FPKMandCounts\\" # TCGA 下载的转录组数据
jsonFilePath <- ".\\01-RawData\\metadata.cart.2020-09-09.json"# 转录组数据对应的json文件
gtfFilesPath <- ".\\02-inputFiles\\human.v33lift37.annotation.gtf"#基因注释文件
geneFilesPath <- ".\\02-inputFiles\\gene.txt" #要分析的基因，比如m6A RNA甲基化相关基因

# 临床数据文件地址：TCGA下载的临床数据
ClinicalDataPath <- ".\\01-RawData\\Clinical2020-09-09\\clinical.csv"

#############################上面代码地址注意修改

setwd(".\\")
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
outputPath<- c(".\\03-AnalysisResult\\101-extractGenesExp\\",
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
