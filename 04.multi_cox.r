
local({r <- getOption("repos")  
r["CRAN"] <- "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"   
options(repos=r)}) 
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")


p="argparse"

if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
  install.packages(p,  warn.conflicts = FALSE)
  suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
}

parser <- ArgumentParser(description='cox regression analysis gene expression')


parser$add_argument("-i", "--data", type="character",required=T,
                    help="input data  file  path[required]",
                    metavar="data")

parser$add_argument( "-t", "--time", type="character",required=T,
                     help="set suvival time column name [required]",
                     metavar="time")
parser$add_argument( "-e", "--event", type="character",required=T, 
                     help="set event  column name must 0 or 1 code format [required]",
                     metavar="event")
parser$add_argument( "-v", "--variate", type="character",nargs='+',required=T, 
                     help=" variate  for cox analysis [required]",
                     metavar="variate")
parser$add_argument( "-P", "--predict.time", type="double",nargs='+',required=F, default=c(365,1095,1825),
                     help=" Time point to draw the ROC curve  [default 365 1095 1825]",
                     metavar="predict.time")
parser$add_argument( "-c", "--cut.score", type="double",required=F, default=NULL,
                     help="set cut score value to divide high and low risk groups  [default median]",
                     metavar="cut.score")
parser$add_argument( "-s", "--seed", type="integer",required=F, default=2021,
                     help="set random seed [default 2021]",
                     metavar="seed")
parser$add_argument( "-o", "--outdir", type="character", default=getwd(),
                     help="output file directory [default cwd]",
                     metavar="outdir")
parser$add_argument("-p", "--prefix", type="character", default="cox",
                    help="out file name prefix [default cox]",
                    metavar="prefix")
opt <- parser$parse_args()


if( !file.exists(opt$outdir) ){
  if( !dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE) ){
    stop(paste("dir.create failed: outdir=",opt$outdir,sep=""))
  }
}



package_list <- c("tidyverse","survival","survminer","survivalROC","rms","GGally","pheatmap","plotROC","RColorBrewer")
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p,  warn.conflicts = FALSE)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}



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


coxdata<-read.table(opt$data,header=T,check.names = F,row.names=1,stringsAsFactors = F,sep="\t",comment.char = "")


Variable=opt$variate

c(opt$time,opt$event,Variable)
c(opt$time,opt$event,Variable) %in% colnames(coxdata)
coxdata=coxdata[,c(opt$time,opt$event,Variable)]

Variable.add=gsub("(.+-.+)","`\\1`",Variable)

f=as.formula(paste0("Surv(",opt$time,",",opt$event,")~",paste(Variable.add,collapse="+")))
f


Multi_coxph_model <- coxph( f,data = coxdata )
#Multi_coxph_model <- cph( f,data = coxdata ,x=T,y=T)

model<-Multi_coxph_model
save(file=paste0(opt$outdir,"/","coxmodel.rda",sep="") ,model)
sink(file=paste0(opt$outdir,"/",opt$prefix,".model.txt",sep=""))
cat("f=")
f
cat("\n===========summary model=====================\n")
summary(Multi_coxph_model)

sink()

#########################Hazard ratio 绘图###########################################
Variable=gsub("`","",Variable)

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


Multi_sum <- summary(Multi_coxph_model)
MHR <- round(Multi_sum$coefficients[,2],2)
Multi_name<-rownames(Multi_sum$coefficients)
MPValue <- Multi_sum$coefficients[,5]
MCIL <- round(Multi_sum$conf.int[,3],2)
MCIU <- round(Multi_sum$conf.int[,4],2)
MCI <- paste0(MCIL,'-', MCIU)

Multi_vars <- data.frame('Characteristics' = Multi_name,
                           "coefficients"=Multi_coxph_model$coefficients,
                           'Hazard Ratio' = MHR,
                           'CI95' = MCI,
                           'P.value' = MPValue)


write.table(Multi_vars, file = paste(opt$outdir,"/",opt$prefix,".Hazard-Ratio.tsv",sep=""), row.names =F,quote = F,sep="\t")

zphmodel<-cox.zph(Multi_coxph_model,transform="km", global=TRUE)

#zphmodel$table

#ggcoxzph(zphmodel)
#dev.off()

write.table(zphmodel$table, file = paste(opt$outdir,"/",opt$prefix,".PHtest.tsv",sep=""), row.names =T,quote = F,sep="\t")

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

vif<-rms::vif(Multi_coxph_model)
cat("vif:","\n")
vif


cat(" stepwise elimination method:","\n",paste0(opt$outdir,"/",opt$prefix,".StepwiseAIC.txt\n",sep=""))
sink(file=paste0(opt$outdir,"/",opt$prefix,".StepwiseAIC.txt",sep=""))
stats::step(Multi_coxph_model)
sink()

# ########################################
# # ##############################nomogram
# ############################################
# coxm <- cph( f,data = coxdata ,x=T,y=T,surv=T)
# survx <- Survival(coxm)
# funs<-lapply(opt$predict.time, function(t){eval(parse(text=paste0("function(x) survx(",t,",x)")))})
# funlabel=paste(opt$predict.time,"Pr(Survival)")
# #funlabel
# #funs
# 
# 
# pdf(file=paste(opt$outdir,"/",opt$prefix,".nomogram.pdf",sep=""), height=10, width=12)
# dd = datadist(coxdata)
# options(datadist='dd')
# plot(nomogram(coxm, fun=funs,
#               funlabel=funlabel))
# dev.off()
# png(filename=paste(opt$outdir,"/",opt$prefix,".nomogram.png",sep=""), height=10*300, width=12*300, res=300, units="px")
# 
# plot(nomogram(coxm, fun=funs,
#               funlabel=funlabel))
# dev.off()
# 



# pdf(file=paste(opt$outdir,"/",opt$prefix,".nomogram1.pdf",sep=""), height=10, width=12)
# 
# regplot(Multi_coxph_model, clickable=F, pdf=T,
#               points=TRUE, rank="sd",failtime = opt$predict.time,droplines=T,prfail = T,
#               other=(list(bvcol="red",sq="green",obscol="blue")))
# dev.off()


################################################################################
#units(coxdata$TIME) <- "day"

for (t in 1:length(opt$predict.time)){
  set.seed(opt$seed)
  res.cox1 <- cph(f, data = coxdata,surv=T,x=TRUE, y=TRUE,time.inc=opt$predict.time[t])
  
  #cal <- calibrate(res.cox1, cmethod='KM', method="boot",u=t,m=38,B=228)
  cal <- calibrate(res.cox1, cmethod='KM', method="boot",u=opt$predict.time[t],m=floor(nrow(coxdata)/3),B=1000)
  
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







############################################################
############################################################
multi_model_risk_values <- predict(Multi_coxph_model,coxdata[,Variable],type="lp")  # log(h(t,x))-log(h0(t)) = b


#predict.time =1000
#plot(mysurv)
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


#ROC
pdf(file = paste(opt$outdir,"/",opt$prefix,".ROC_result.pdf",sep="") ,height = 5,width = 5)
print(p)
dev.off()
png(filename  = paste(opt$outdir,"/",opt$prefix,".ROC_result.png",sep="") ,height = 5*300,width = 5*300,res=300,units = "px")
print(p)
dev.off()






############################################################
############################################################


#method: Kaplan-Meier survival curve was analyzed by “survival” R package and evaluated by log-rank test.

mysurv <- Surv(coxdata[,opt$time], coxdata[,opt$event])
coxdata$riskScore <- multi_model_risk_values

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
#coxdata$index <- c(1:nrow(coxdata))

c(1:nrow(coxdata))

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


####################signature heatmap  biomarker

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



write.table(data.frame(ID=row.names(coxdata),coxdata,check.names = F), file = paste(opt$outdir,"/",opt$prefix,".Risk_Score.tsv",sep=""), row.names =F,quote = F,sep="\t")






