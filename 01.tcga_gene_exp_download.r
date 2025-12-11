
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


#######命令行参数设置################################################################################
parser <- ArgumentParser(description='TCGA_download')

parser$add_argument( "-p", "--project", type="character",required=T,
		help="input  project ID of TCGA,  for  example TCGA-KIRC[required]",
		metavar="project")
parser$add_argument( "-f", "--files.per.chunk", type="integer",required=F,default=10,
                     help="This will make the API method only download n (files.per.chunk) files at a time. This may reduce the download problems when the data size is too large [default 10]",
                     metavar="files.per.chunk")
parser$add_argument("-m", "--more.sample.info", action='store_true',
                    help="keep more sample info  [optional, default: False]")
parser$add_argument( "-o", "--outdir", type="character", default=getwd(),
		help="output file directory [default cwd]",
		metavar="outdir")

opt <- parser$parse_args()

if( !file.exists(opt$outdir) ){
	if( !dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE) ){
		stop(paste("dir.create failed: outdir=",opt$outdir,sep=""))
	}
}



###############################################################################################
#############################################################
package_list <- c("SummarizedExperiment","TCGAbiolinks")

for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(p)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}




####################################################




ensemblToName<-function(expr,gene.info){
  expr=as.data.frame(expr)
  #顺序保持一致
  gene.info<-as.data.frame(gene.info)
  row.names(gene.info)<-gene.info$gene_id
  gene.info<-gene.info[rownames(expr),]
  
  sample.id=colnames(expr)
  #添加gene name列  相同基因随机选取
  gene_name_rmdup=!duplicated(gene.info$gene_name)
  
  #用dplyr包，相同的 gene_name 求平均  
  
  #expr %>% group_by(gene_name) %>% summarise(across(everything(), list(sum)))->expr.res
  
  #expr.res<-as.data.frame(expr.res)
  
  expr.res=expr[gene_name_rmdup,]
  gene.info=gene.info[gene_name_rmdup,]
  rownames(expr.res)<-gene.info$gene_name
  
  expr.res
  
}

#######################################################

## 2.2 数据下载代码
#设置输出目录
##setwd(opt$outdir)
		
		
#######################################################################
#######################################################################
		
# 查询可以下载的数据,详情可参考：
clinical <- GDCquery_clinic(project = opt$project, type = "clinical")
dim(clinical)
###########################################
#根据临床信息，筛选一下样品
##########################################

# 临床数据输出表格
write.table(clinical, file = paste0(opt$outdir,"/",opt$project,'_clinical_index.tsv'), sep="\t",row.names =F,quote = F)

################################原始XML文件 最全的临床数据########################
query <- GDCquery(project = opt$project, 
data.category = "Clinical", 
file.type = "xml")
GDCdownload(query,directory = opt$project,files.per.chunk=opt$files.per.chunk, method='api')

##整理所有的数据并保存
clinical.info<-c("drug","follow_up","radiation","patient","stage_event","new_tumor_event","admin")
for(i in clinical.info){
clinical <- GDCprepare_clinic(query, clinical.info = i,directory = opt$project)
write.table(clinical, file = paste0(opt$outdir,"/",opt$project,'_clinical_xml_',i,'.tsv'), sep="\t", row.names =F,quote = F)
}
		

## 2.3 下载RNA-seq表达数据

##############################HTSeq - Counts######################
query <- GDCquery(project = opt$project,
	data.category = "Transcriptome Profiling",  
	data.type = "Gene Expression Quantification",
	workflow.type = "STAR - Counts",
	#sample.type = sample_type,
	legacy = FALSE)
files <- getResults(query,cols=c("cases"))
cat("Total files to download:", length(files),"\n")

# TP
dataSmTP <- TCGAquery_SampleTypes(barcode = files, typesample = "TP")
cat("Total TP samples to down:", length(dataSmTP),"\n")

# NT
dataSmNT <- TCGAquery_SampleTypes(barcode = files,typesample = "NT")
cat("Total NT samples to down:", length(dataSmNT),"\n")

GDCdownload(query = query,directory = opt$project,files.per.chunk=opt$files.per.chunk, method='api')
# 保存整理下载数据结果
gene.data <- GDCprepare(query = query, 
	save = TRUE, 
	directory =  opt$project,
	save.filename = paste0(opt$outdir,"/",opt$project,"_gene_expression_Counts.rda"))

sample.info.list=colData(gene.data)@listData
sample.info=NULL
if(opt$more.sample.info){
  rmcol<-c("treatments","primary_site","disease_type")
  sample.info<-as.data.frame(sample.info.list[!names(sample.info.list)%in%rmcol])
}else{
  sample.info=as.data.frame(sample.info.list[1:10])
}

write.table(sample.info, file = paste0(opt$outdir,"/",opt$project,'_sample_info.tsv'),sep="\t", row.names =F, quote = F)
		
gene.info=as.data.frame(rowRanges(gene.data))

write.table(gene.info, file = paste0(opt$outdir,"/",opt$project,'_gene_info.tsv'), sep="\t", row.names =F, quote = F)

# 表达量提取
# unstranded stranded_first stranded_second tpm_unstrand fpkm_unstrand fpkm_uq_unstrand
data_expr <- assay(gene.data,i = "unstranded")

data_expr<-ensemblToName(data_expr,gene.info)
write.table(data.frame(ID=rownames(data_expr),data_expr,check.names = F), file = paste0(opt$outdir,"/",opt$project,'_gene_expression_Counts.tsv'), sep="\t",row.names =F, quote = F)

# 表达量提取
data_expr_fpkm_uq <- assay(gene.data,i = "fpkm_uq_unstrand")
data_expr_fpkm_uq<-ensemblToName(data_expr_fpkm_uq,gene.info)
write.table(data.frame(ID=rownames(data_expr_fpkm_uq),data_expr_fpkm_uq,check.names = F), file = paste0(opt$outdir,"/",opt$project,'_gene_expression_FPKM-UQ.tsv'), sep="\t",row.names =F, quote = F)	
# 表达量提取
data_expr_fpkm <- assay(gene.data,i = "fpkm_unstrand")
data_expr_tpm <- assay(gene.data,i = "tpm_unstrand")
data_expr_fpkm<-ensemblToName(data_expr_fpkm,gene.info)
write.table(data.frame(ID=rownames(data_expr_fpkm),data_expr_fpkm,check.names = F), file = paste0(opt$outdir,"/",opt$project,'_gene_expression_FPKM.tsv'), sep="\t",row.names =F, quote = F)
data_expr_tpm<-ensemblToName(data_expr_tpm,gene.info)
write.table(data.frame(ID=rownames(data_expr_tpm),data_expr_tpm,check.names = F), file = paste0(opt$outdir,"/",opt$project,'_gene_expression_TPM.tsv'), sep="\t",row.names =F, quote = F)

	ncRNA_type<-c("3prime_overlapping_ncrna",
			"antisense",
			"lncRNA",
			"macro_lncRNA",
			"processed_transcript",
			"sense_intronic",
			"sense_overlapping"
	)
    
    protein_coding_gene_type<-c("IG_C_gene",
                                "IG_D_gene",
                                "IG_J_gene",
                                "IG_V_gene",
                                "polymorphic_pseudogene",
                                "protein_coding",
                                "TR_C_gene",
                                "TR_D_gene",
                                "TR_J_gene",
                                "TR_V_gene"
    )
    
	#	protein_coding_gene<-subset(gene.info,gene_type == "protein_coding")
protein_coding_gene<-gene.info[gene.info$gene_type %in% protein_coding_gene_type,]
lncRNA_gene<-gene.info[gene.info$gene_type %in% ncRNA_type,]

#counts
data_expr_protein_coding<-data_expr[row.names(data_expr)%in%protein_coding_gene$gene_name,]
data_expr_lncRNA<-data_expr[row.names(data_expr)%in%lncRNA_gene$gene_name,]


write.table(data.frame(ID=rownames(data_expr_protein_coding),data_expr_protein_coding,check.names = F), file = paste0(opt$outdir,"/",opt$project,'_gene_expression_Counts_protein_coding.tsv'), sep="\t",row.names =F, quote = F)
write.table(data.frame(ID=rownames(data_expr_lncRNA),data_expr_lncRNA,check.names = F), file = paste0(opt$outdir,"/",opt$project,'_gene_expression_Counts_lncRNA.tsv'), sep="\t", row.names =F, quote = F)


# FPKM
data_expr_fpkm_protein_coding<-data_expr_fpkm[row.names(data_expr_fpkm)%in%protein_coding_gene$gene_name,]
data_expr_fpkm_lncRNA<-data_expr_fpkm[row.names(data_expr_fpkm)%in%lncRNA_gene$gene_name,]


write.table(data.frame(ID=rownames(data_expr_fpkm_protein_coding),data_expr_fpkm_protein_coding,check.names = F), file = paste0(opt$outdir,"/",opt$project,'_gene_expression_FPKM_protein_coding.tsv'), sep="\t", row.names =F, quote = F)
write.table(data.frame(ID=rownames(data_expr_fpkm_lncRNA),data_expr_fpkm_lncRNA,check.names = F), file = paste0(opt$outdir,"/",opt$project,'_gene_expression_FPKM_lncRNA.tsv'), row.names =F, sep="\t", quote = F)
# FPKM-UQ
data_expr_fpkm_protein_coding_uq<-data_expr_fpkm_uq[row.names(data_expr_fpkm_uq)%in%protein_coding_gene$gene_name,]
data_expr_fpkm_lncRNA_uq<-data_expr_fpkm_uq[row.names(data_expr_fpkm_uq)%in%lncRNA_gene$gene_name,]


write.table(data.frame(ID=rownames(data_expr_fpkm_protein_coding_uq),data_expr_fpkm_protein_coding_uq,check.names = F), file = paste0(opt$outdir,"/",opt$project,'_gene_expression_FPKM-UQ_protein_coding.tsv'), sep="\t", row.names =F, quote = F)
write.table(data.frame(ID=rownames(data_expr_fpkm_lncRNA_uq),data_expr_fpkm_lncRNA_uq,check.names = F), file = paste0(opt$outdir,"/",opt$project,'_gene_expression_FPKM-UQ_lncRNA.tsv'), row.names =F, sep="\t", quote = F)


# TPM
data_expr_tpm_protein_coding<-data_expr_tpm[row.names(data_expr_tpm)%in%protein_coding_gene$gene_name,]
data_expr_tpm_lncRNA<-data_expr_tpm[row.names(data_expr_tpm)%in%lncRNA_gene$gene_name,]


write.table(data.frame(ID=rownames(data_expr_tpm_protein_coding),data_expr_tpm_protein_coding,check.names = F), file = paste0(opt$outdir,"/",opt$project,'_gene_expression_TPM_protein_coding.tsv'), sep="\t", row.names =F, quote = F)
write.table(data.frame(ID=rownames(data_expr_tpm_lncRNA),data_expr_tpm_lncRNA,check.names = F), file = paste0(opt$outdir,"/",opt$project,'_gene_expression_TPM_lncRNA.tsv'), row.names =F, sep="\t", quote = F)


cat("patient num:",length(unique(sample.info$patient)),"\n")

cat("sample num:",length(unique(sample.info$sample)),"\n")

cat("TP or NT sample num:\n")
table(sample.info$definition)
