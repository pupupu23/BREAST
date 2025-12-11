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

parser <- ArgumentParser(description='lasso cox regression analysis')

parser$add_argument("-i", "--data", type="character",required=T,
                    help="input data  file  path[required]",
                    metavar="data")

parser$add_argument( "-t", "--time", type="character",required=T,
                     help="set suvival time column name [required]",
                     metavar="time")
parser$add_argument( "-e", "--event", type="character",required=T, 
                     help="set event  column name [required]",
                     metavar="event")
parser$add_argument( "-v", "--variate", type="character",nargs='+',required=T, 
                     help=" variate  for cox analysis [required]",
                     metavar="variate")
parser$add_argument( "-s", "--seed", type="integer",required=F, default=2021,
                     help="set random seed [default 2021]",
                     metavar="seed")
parser$add_argument( "-l", "--lambda", type="integer",required=F, default=NULL,
                     help="set lambda cut off [default NULL]",
                     metavar="lambda")
parser$add_argument( "-P", "--predict.time", type="double",nargs='+',required=F, default=c(365,1095,1825),
                      help=" Time point of the ROC curve to select cutoff [default 365 1095 1825]",
                      metavar="predict.time")
# parser$add_argument( "--log2", action='store_true',
#                       help="whether do log2 transfrom for expression data [optional, default: False]")
parser$add_argument( "-o", "--outdir", type="character", default=getwd(),
                     help="output file directory [default cwd]",
                     metavar="outdir")
parser$add_argument("-p", "--prefix", type="character", default="lasso_cox",
                    help="out file name prefix [default lasso_cox]",
                    metavar="prefix")
# parser$add_argument( "-H", "--height", type="double", default=10,
#                      help="the height of pic   inches  [default 10]",
#                      metavar="height")
# parser$add_argument("-W", "--width", type="double", default=10,
#                     help="the width of pic   inches [default 10]",
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

package_list <- c("tidyverse","glmnet","survival","survminer","RColorBrewer","survivalROC","rms","GGally","pheatmap","plotROC")
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p,  warn.conflicts = FALSE)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

# package_list <- c("RegParallel")
# 
# for(p in package_list){
#   if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
#     if (!requireNamespace("BiocManager", quietly = TRUE))
#       install.packages("BiocManager")
#     BiocManager::install(p)
#     suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
#   }
# }

mycol=c(brewer.pal(6,"Set1"),brewer.pal(8,"Dark2"),brewer.pal(7,"Accent"))


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





# f=as.formula(paste0("Surv(",opt$time,",",opt$event,")~",paste(Variable,collapse="+")))
# f


coxdata<-read.table(opt$data,header=T,check.names = F,row.names=1,stringsAsFactors = F,sep="\t",comment.char = "")

Variable=opt$variate
x=as.matrix(coxdata[,Variable]) 
coxdata=coxdata[,c(opt$time,opt$event,Variable)]

y=data.matrix(Surv(coxdata[,opt$time],coxdata[,opt$event])) 


set.seed(opt$seed)
fit <- glmnet(x, y, family = "cox", maxit = 1000,alpha = 1)  #cox

print(fit)

pdf(file=paste(opt$outdir,"/",opt$prefix,".lasso_lambda.pdf",sep=""), height=5, width=7)

plot(fit,xvar = "lambda",label = TRUE)
dev.off()
png(filename=paste(opt$outdir,"/",opt$prefix,".lasso_lambda.png",sep=""), height=5*300, width=7*300, res=300, units="px")

plot(fit,xvar = "lambda",label = TRUE)
dev.off()



#Cross-validation
set.seed(opt$seed)
cvfit <- cv.glmnet(x, y, family = "cox", type.measure = "deviance", maxit = 1000,nfolds=10)

pdf(file=paste(opt$outdir,"/",opt$prefix,".lasso_cv.pdf",sep=""), height=5, width=7)
plot(cvfit) 
#其中两条虚线分别指示了两个特殊的λ值 
abline(v = log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")

dev.off()
png(filename=paste(opt$outdir,"/",opt$prefix,".lasso_cv.png",sep=""), height=5*300, width=7*300, res=300, units="px")
plot(cvfit) 
#其中两条虚线分别指示了两个特殊的λ值 
abline(v = log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")

dev.off()





lambda=NULL
lassoGene=NULL
actCoef=NULL
if(is.null(opt$lambda)){
    
  
  lambda=cvfit$lambda.min
  
  coef = coef(fit, s = lambda) 
  
  print("select the optimal OS-related signature genes with nonzero coefficients\n")
  print(coef)
  index = which(coef != 0) 
  actCoef = coef[index] 
  lassoGene = row.names(coef)[index] 
  geneCoef = cbind(Gene=lassoGene,Coef=actCoef) 
  geneCoef   
  write.table(geneCoef, file = paste(opt$outdir,"/",opt$prefix,".geneCoef.lambda.min.tsv",sep=""), row.names =F,quote = F,sep="\t")
  
  
  #write.table(geneCoef, file = paste(opt$outdir,"/",opt$prefix,".geneCoef.lambda.1se.tsv",sep=""), row.names =F,quote = F,sep="\t")
  
  
  
  #x=as.matrix(coxdata[,Variable]) 
  #y=data.matrix(Surv(coxdata[,opt$time],coxdata[,opt$event]))

}else{
  lambda=opt$lambda
  coef = coef(fit, s = lambda) 
  print("select the optimal OS-related GRGs with nonzero coefficients\n")
  index = which(coef != 0) 
  actCoef = coef[index] 
  lassoGene = row.names(coef)[index] 
  geneCoef = cbind(Gene=lassoGene,Coef=actCoef) 
  geneCoef
  
  write.table(geneCoef, file = paste(opt$outdir,"/",opt$prefix,".geneCoef.lambda.tsv",sep=""), row.names =F,quote = F,sep="\t")
  
}

lassoGene.add=gsub("(.+-.+)","`\\1`",lassoGene)

f=as.formula(paste0("Surv(",opt$time,",",opt$event,")~",paste(lassoGene.add,collapse="+")))
f

Multi_coxph_model <- coxph( f,data = coxdata ,init=actCoef,iter.max=0)
Multi_coxph_model



model=Multi_coxph_model
save(file=paste0(opt$outdir,"/","lassomodel.rda",sep="") ,model)

sink(file=paste0(opt$outdir,"/",opt$prefix,".lasso.model.txt",sep=""))
cat("f=")
f

cat("\n===========summary model=====================\n")
summary(Multi_coxph_model)


sink()

#########################Hazard ratio###########################################
Variable=gsub("`","",lassoGene)

pdf(file=paste(opt$outdir,"/",opt$prefix,".Hazard-Ratio.pdf",sep=""), height=10, width=10)
ggforest(Multi_coxph_model,
         data = coxdata,
         main = 'Hazard ratio',
         cpositions = c(0.05, 0.15, 0.35),
         fontsize = 1,
         noDigits = 3
)
dev.off()
png(filename=paste(opt$outdir,"/",opt$prefix,".Hazard-Ratio.png",sep=""), height=10*300, width=10*300, res=300, units="px")
ggforest(Multi_coxph_model,
         data = coxdata,
         main = 'Hazard ratio',
         cpositions = c(0.05, 0.15, 0.35),
         fontsize = 1,
         noDigits = 3
)
dev.off()




# 汇总多因素分析结果
Multi_sum <- summary(Multi_coxph_model)
MHR <- round(Multi_sum$coefficients[,2],2)
Multi_name<-rownames(Multi_sum$coefficients)
MPValue <- Multi_sum$coefficients[,5]
MCIL <- round(Multi_sum$conf.int[,3],2)
MCIU <- round(Multi_sum$conf.int[,4],2)
MCI <- paste0(MCIL,'-', MCIU)
Multi_name=gsub("`","",Multi_name)

Multi_vars <- data.frame('Characteristics' = Multi_name,
                         "coefficients"=Multi_coxph_model$coefficients,
                         'Hazard Ratio' = MHR,
                         'CI95' = MCI,
                         'P.value' = MPValue)


write.table(Multi_vars, file = paste(opt$outdir,"/",opt$prefix,".Hazard-Ratio.tsv",sep=""), row.names =T,quote = F,sep="\t")


##################################################
zphmodel<-cox.zph(Multi_coxph_model,transform="km", global=TRUE)

#zphmodel$table

#ggcoxzph(zphmodel)
#dev.off()

write.table(zphmodel$table, file = paste(opt$outdir,"/",opt$prefix,".PHtest.tsv",sep=""), row.names =T,quote = F,sep="\t")



#相关性分析
corp=ggpairs(coxdata[,Variable],axisLabels = "show") +theme_bw()+ theme(  
  panel.grid=element_blank(), 
  axis.text.x=element_text(colour="black"),
  axis.text.y=element_text(colour="black"),
  panel.border=element_rect(colour = "black"),
  legend.key = element_blank(),plot.title = element_text(hjust = 0.5))

pdf(file=paste(opt$outdir,"/",opt$prefix,".gene_expr_cor.pdf",sep=""), height=length(Variable), width=length(Variable))
print(corp)
dev.off()
png(filename=paste(opt$outdir,"/",opt$prefix,".gene_expr_cor.png",sep=""), height=length(Variable)*300, width=length(Variable)*300, res=300, units="px")
print(corp)
dev.off()

for (t in 1:length(opt$predict.time)){
  set.seed(opt$seed)
  res.cox1 <- cph(f, data = coxdata,surv=T,x=TRUE, y=TRUE,time.inc=opt$predict.time[t],init=actCoef,iter.max=0)
  
  #print(Variable,        lassoGene)
  #cal <- calibrate(res.cox1, cmethod='KM', method="boot",u=t,m=38,B=228)
  cal <- calibrate(res.cox1, cmethod='KM', method="boot",u=opt$predict.time[t],m=floor(nrow(coxdata)/3))
  
  
  png(filename=paste(opt$outdir,"/",opt$prefix,sprintf(".nomogram_calibrate%s.png",opt$predict.time[t]),sep=""), height=7*300, width=7*300, res=300, units="px")
  par(mar=c(5,4,4,2)+0.1)
  plot(cal,lwd=2,lty=1,cex.lab=1.2, cex.axis=1, cex.main=1.2, cex.sub=0.6, 
       errbar.col=mycol[t],sub=F,
       xlab=sprintf("Nomogram-Predicted Probability of %s  OS",opt$predict.time[t])   ,
       ylab=sprintf("Actual  %s DFS (proportion)" ,opt$predict.time[t]) ,
       col=mycol[t])
  dev.off()
  
  pdf(file=paste(opt$outdir,"/",opt$prefix,sprintf(".nomogram_calibrate%s.pdf",opt$predict.time[t]),sep=""), height=7, width=7)
  par(mar=c(5,4,4,2)+0.1)
  plot(cal,lwd=2,lty=1,cex.lab=1.2, cex.axis=1, cex.main=1.2, cex.sub=0.6, 
       errbar.col=mycol[t],sub=F,
       xlab=sprintf("Nomogram-Predicted Probability of %s  OS",opt$predict.time[t])   ,
       ylab=sprintf("Actual  %s DFS (proportion)" ,opt$predict.time[t]) ,
       col=mycol[t])
  dev.off()
  
  
}



multi_model_risk_values <- predict(Multi_coxph_model,coxdata[,lassoGene],type="lp")

mylabel=NULL
cut_value=NULL
AUC=NULL
sroclong=NULL
for(t in opt$predict.time){
  
  ROC_opt= survivalROC(Stime=coxdata[,opt$time],
                       status=coxdata[,opt$event],
                       marker = multi_model_risk_values,
                       predict.time = t,method="KM")
  
  #print(ROC_opt)
  
  min_index <- with(ROC_opt, which.min(1-TP+ FP))
  cut_value = c(cut_value,ROC_opt$cut.values[min_index])
  AUC=c(AUC,round(ROC_opt$AUC,3))
  
  mylabel=c(mylabel,paste("AUC at",t,"=",round(ROC_opt$AUC,3)))
  
  sroclong=rbind(sroclong,data.frame(TPF = ROC_opt[["TP"]], FPF = ROC_opt[["FP"]], 
                                     c = ROC_opt[["cut.values"]], 
                                     time = rep(ROC_opt[["predict.time"]], length(ROC_opt[["FP"]]))))
}

#sroclong <- do.call(rbind, sroc)
write.table(sroclong, file = paste(opt$outdir,"/",opt$prefix,".roc.res.tsv",sep=""), row.names =F,quote = F,sep="\t")

data.frame(time=opt$predict.time,auc=AUC,optimal_cut_score=cut_value)


basicplot=ggplot(sroclong, aes(x = FPF, y = TPF, label = c, color = as.factor(time))) + 
  geom_roc(labels = FALSE, stat = "identity",size=0.8,n.cuts = 0) + 
  geom_abline(intercept = 0, slope = 1,color="black",size = 0.8,linetype="dashed",alpha=0.8)+
  #style_roc()+
  #style_roc(guide=guides(pcol="black"),theme = theme_bw)+
  xlab("False positive rate")+ ylab("True positive rate")+
  scale_color_brewer(palette = "Set1",name="Time")+
  theme_bw()+ theme(  
    panel.grid=element_blank(), 
    axis.text.x=element_text(colour="black"),
    axis.text.y=element_text(colour="black"),
    panel.border=element_rect(colour = "black"),
    legend.key = element_blank(),plot.title = element_text(hjust = 0.5))

p=basicplot+annotate("text", x = 0.5, y = 0.25,hjust = 0,
                     label = paste(mylabel,collapse = "\n"))


# ROC
pdf(file = paste(opt$outdir,"/",opt$prefix,".ROC_result.pdf",sep="") ,height = 5,width = 5)
print(p)
dev.off()
png(filename  = paste(opt$outdir,"/",opt$prefix,".ROC_result.png",sep="") ,height = 5*300,width = 5*300,res=300,units = "px")
print(p)
dev.off()






############################################################
############################################################


mysurv <- Surv(coxdata[,opt$time], coxdata[,opt$event])
coxdata$riskScore <- predict(Multi_coxph_model,coxdata[,lassoGene],type="risk")

cut.score=NULL
if(is.null(opt$cut.score)){
  cut.score=median(coxdata$riskScore)
}else{
  cut.score=opt$cut.score
}

cat("cut.score:")
cut.score



group=ifelse(coxdata$riskScore>=cut.score,'high','low')

#group=factor(group,levels = c("low","high"),order=T)

coxdata$risk.group <- group
kmfit <- survfit(mysurv ~ coxdata$risk.group,data=coxdata)
#summary(kmfit)


p<-ggsurvplot(kmfit,data=coxdata,
              pval = TRUE, conf.int = FALSE,
              risk.table = TRUE, 
              risk.table.col = "strata", 
              linetype = "strata", 
              surv.median.line = "hv", 
              ggtheme = theme_classic(), 
              palette =  "Set1",
              ncensor.plot = TRUE
)
#p$table <- p$table + theme(axis.line = element_blank())
p$plot <- p$plot + labs(title = "Survival Curves")



pdf(file=paste(opt$outdir,"/",opt$prefix,".km_plot.pdf",sep=""), height=10, width=6)
print(p)
dev.off()
png(filename=paste(opt$outdir,"/",opt$prefix,".km_plot.png",sep=""), height=10*300, width=6*300, res=300, units="px")
print(p)
dev.off()


# Risk Score
coxdata <- coxdata[order(coxdata$riskScore),]


p=ggplot(coxdata,aes(x=c(1:nrow(coxdata)), y=riskScore,color=risk.group))+ 
  geom_point()+geom_vline(xintercept = table(coxdata$risk.group)["low"],linetype =2)+
  scale_color_manual(values = c("high"="#e41a1c","low"="#377eb8"))+
  theme_bw()+ylab("Risk Score")+xlab("Index")+ theme(  
    panel.grid=element_blank(), 
    axis.text.x=element_text(colour="black"),
    axis.text.y=element_text(colour="black"),
    panel.border=element_rect(colour = "black"),
    legend.key = element_blank(),plot.title = element_text(hjust = 0.5))



pdf(file=paste(opt$outdir,"/",opt$prefix,".Risk_Score.pdf",sep=""), height=4, width=7)
print(p)
dev.off()
png(filename=paste(opt$outdir,"/",opt$prefix,".Risk_Score.png",sep=""), height=4*300, width=7*300, res=300, units="px")
print(p)
dev.off()


STATUS<-ifelse(coxdata[,opt$event]==1,"Dead","Alive")

p=ggplot(coxdata, aes(x=c(1:nrow(coxdata)), y=coxdata[,opt$time], colour=as.factor(STATUS))) + 
  geom_point()+scale_color_manual(values = c("Dead"="#e41a1c","Alive"="#377eb8"),name="")+
  theme_bw()+ylab("Survival TIME")+xlab("Index")+theme(  
    panel.grid=element_blank(), 
    axis.text.x=element_text(colour="black"),
    axis.text.y=element_text(colour="black"),
    panel.border=element_rect(colour = "black"),
    legend.key = element_blank(),plot.title = element_text(hjust = 0.5))

pdf(file=paste(opt$outdir,"/",opt$prefix,".survival_time_scatter.pdf",sep=""), height=4, width=7)
print(p)
dev.off()
png(filename=paste(opt$outdir,"/",opt$prefix,".survival_time_scatter.png",sep=""), height=4*300, width=7*300, res=300, units="px")
print(p)
dev.off()

annotation_col = data.frame(risk.group = factor(coxdata$risk.group))  

rownames(annotation_col)= rownames(coxdata)

gene_expr <- t(coxdata[,Variable])

#colnames(gene_expr) <- rownames(coxdata)





ann_colors = list(
  risk.group = c("high"="#e41a1c","low"="#377eb8")
  
)




pdf(file=paste(opt$outdir,"/",opt$prefix,".expr_heatmap.pdf",sep=""), height=5, width=9)
pheatmap(gene_expr,annotation_col=annotation_col,
         cluster_cols =  FALSE,scale = "row", 
         color = colorRampPalette(c("blue", "white", "red"))(100),
         show_colnames = F,border=FALSE,
         annotation_colors = ann_colors)

dev.off()
png(filename=paste(opt$outdir,"/",opt$prefix,".expr_heatmap.png",sep=""), height=5*300, width=9*300, res=300, units="px")
pheatmap(gene_expr,annotation_col=annotation_col,
         cluster_cols =  FALSE,scale = "row", 
         color = colorRampPalette(c("blue", "white", "red"))(100),
         show_colnames = F,border=FALSE,
         annotation_colors = ann_colors)

dev.off()


no.lassoGene=setdiff(opt$variate,lassoGene)

write.table(data.frame(ID=row.names(coxdata),coxdata[,!colnames(coxdata)%in%no.lassoGene],check.names = F), file = paste(opt$outdir,"/",opt$prefix,".Risk_Score.tsv",sep=""), row.names =F,quote = F,sep="\t")












