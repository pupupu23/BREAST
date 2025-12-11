#设置镜像，
local({r <- getOption("repos")  
r["CRAN"] <- "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"   
options(repos=r)}) 



#2 安装bioconductor常用包
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")



p="argparse"

if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
  install.packages(p,  warn.conflicts = FALSE)
  suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
}

parser <- ArgumentParser(description='batch unvariate cox regression  gene expression')


parser$add_argument("-m", "--metadata", type="character",required=T,
                    help="input metadata file path with suvival time [required]",
                    metavar="metadata")
parser$add_argument( "-g", "--expset", type="character",required=T,
                     help="input gene expression set file [required]",
                     metavar="expset")
# parser$add_argument( "-g", "--geneinfo", type="character",required=F, default=NULL,
#                      help="input gene info to add gene Name to result [default NULL]",
#                      metavar="geneinfo")
# parser$add_argument( "-f", "--follow", type="character",required=F, default=30,
#                      help="follow-up ≤ the days were excluded from the survival analysis [default 30]",
#                      metavar="follow-up")

parser$add_argument( "-B", "--by", type="character",required=F,default="barcode",
                     help="input sample ID column name in metadata [default %(default)s]",
                     metavar="by")
parser$add_argument( "-t", "--time", type="character",required=F,default="TIME",
                     help="set suvival time column name in metadata [default TIME]",
                     metavar="time")
parser$add_argument( "-e", "--event", type="character",required=F, default="EVENT",
                     help="set event  column name in metadata [default EVENT]",
                     metavar="event")


parser$add_argument( "-l", "--pvalue", type="double",required=F, default=0.01,
                     help="pvalue cutoff to choose sig gene  [default 0.01]",
                     metavar="pvalue")
parser$add_argument( "-b", "--blocksize", type="integer",required=F, default=2,
                     help="Number of variables Parallel to test in each   [default 2]",
                     metavar="blocksize")
parser$add_argument( "--log2", action='store_true',
                     help="whether do log2 transfrom for expression data [optional, default: False]")
parser$add_argument( "-o", "--outdir", type="character", default=getwd(),
                     help="output file directory [default cwd]",
                     metavar="outdir")
parser$add_argument("-p", "--prefix", type="character", default="cox",
                    help="out file name prefix [default cox]",
                    metavar="prefix")
# parser$add_argument( "-H", "--height", type="double", default=5,
#                      help="the height of pic   inches  [default 5]",
#                      metavar="height")
# parser$add_argument("-W", "--width", type="double", default=5,
#                     help="the width of pic   inches [default 5]",
#                     metavar="width")
opt <- parser$parse_args()

#parser$print_help()
#set some reasonable defaults for the options that are needed,
#but were not specified.
if( !file.exists(opt$outdir) ){
  if( !dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE) ){
    stop(paste("dir.create failed: outdir=",opt$outdir,sep=""))
  }
}



#包的加载

package_list <- c("tidyverse","survival")
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p,  warn.conflicts = FALSE)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

package_list <- c("RegParallel")

for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(p)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}



# mycolorc=list(
#   brewer.pal(5, "Blues"),
#   brewer.pal(5, "BuGn"),
#   brewer.pal(5, "BuPu"),
#   brewer.pal(5, "Oranges"),
#   brewer.pal(5, "Reds"),
#   brewer.pal(5, "Greys"))
# mycolor=c(
#   brewer.pal(9, "Set1"),
#   brewer.pal(8, "Set2"),
#   brewer.pal(12, "Set3"),
#   brewer.pal(12, "Paired"),
#   brewer.pal(8, "Dark2"),
#   brewer.pal(8, "Accent"))
#-----------------------------------------------------------------
# reading data
#-----------------------------------------------------------------

metadata<-read.table(opt$metadata,header=T,check.names = F,stringsAsFactors = F,sep="\t",comment.char = "")
geneExpset<-read.table(opt$expset,header=T,check.names = F,row.names=1,stringsAsFactors = F,sep="\t",comment.char = "")


#合并表达数据,样本取交集

sampleID<-metadata[,opt$by]
geneExpset=geneExpset[,sampleID]

if(opt$log2){
  print("log2 transform:  \n")
  geneExpset=log2(geneExpset+1)
}

#合并临床数据与表达数据

index=ncol(metadata)+1

coxdata=cbind(metadata,t(geneExpset))


write.table(coxdata,file = paste0(opt$outdir,"/",opt$prefix,".metadata-exp.tsv"),sep="\t",quote = F,row.names = F)

#coxdata=lung


f=paste0("Surv(",opt$time,",",opt$event,")~ [*]")

coln=colnames(coxdata)
coln=gsub("-","xxxxx",coln)
colnames(coxdata)=coln

res5 <- RegParallel(
  data = coxdata,
  formula = f,
  FUN = function(formula, data)
    coxph(formula = formula,
          data = data,
          ties = 'breslow',
          singular.ok = TRUE),
          FUNtype = 'coxph',
          variables = colnames(coxdata)[index:ncol(coxdata)],
          blocksize = opt$blocksize,
          p.adjust = "none")
res5 <- res5[!is.na(res5$P),]
#gene_name=geneInfo[res5$Variable,]
#res5$gene.name<-gene_name$external_gene_name



res5=as.data.frame(res5)
res5$Variable=gsub("xxxxx","-",res5$Variable)
res5$Term=gsub("xxxxx","-",res5$Term)




#mycol<-c("Variable",	"Term",	"Beta",	"StandardError",	"Z",	"P",	"LRT",	"Wald",	"LogRank"	,"HR",	"HRlower",	"HRupper")
mycol<-c("Variable",	"Term",	"Beta"	,"HR",	"HRlower",	"HRupper",		"P")

res5<-res5[,mycol]


res5 <- res5[order(res5$P, decreasing = FALSE),]
final <- subset(res5, P < opt$pvalue)




write.table(res5,file = paste0(opt$outdir,"/",opt$prefix,".all.res.tsv"),sep="\t",quote = F,row.names = F)

write.table(final,file = paste0(opt$outdir,"/",opt$prefix,".sig.res.tsv"),sep="\t",quote = F,row.names = F)

