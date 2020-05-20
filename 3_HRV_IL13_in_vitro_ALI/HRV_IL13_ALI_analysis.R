#Load libraries
library(DESeq2)
library(edgeR)
library(openxlsx)
library(heatmap3)
library(plotrix)
library(venn)
library(grid)

######################## PRELIMINARIES
#Bring in data
endo_dir<-"/Seibold/proj/Endotype/190426_A00405_0096_BHK5YCDSXX/count_matrix/count.matrix.hgnc.ENS84.txt"
endo<-read.table(endo_dir,header=T)

#Keep only endotype columns
endo<-endo[,c(25:64)]

#Order columns to match order of the design matrix (imported below)
endo<-endo[,order(as.numeric(sapply(strsplit(colnames(endo),"_"),function(x)x[3])))]

#Read in design table
endo_design<-read.table("endo_design.txt",header=T,stringsAsFactors=F)

#Sort design more intuitively
endo_design<-endo_design[with(endo_design,order(treat1,treat2,donor)),]

#Sort expression matrix by sample ID
endo<-endo[,with(endo,endo_design$sample)]

#Rename samples according to treatments
colnames(endo)<-apply(endo_design[,c(3,4,5)],1,function(x)paste(x,collapse="_"))

#Toss ribosomal/mtDNA-pseudogenes (1,878 gone)
endo<-endo[-grep("^MT-|^MTAT|^MTCO|^MTCY|^MTND|^MTRN|^MRPL|^MRPS|^RPL|^RPS",rownames(endo)),]

#Toss genes that don't have at least a count of 1 in 2 or more samples
endo<-endo[-which(rowSums(endo != 0) <= 1),]



######################## PLOTING COVID-19 GENES (HRV)
#subset data
endo_hrv_bsa<-endo[,c(1:5,11:15)]

#Specify design
subject <- sapply(strsplit(colnames(endo_hrv_bsa), split ="_"),function(x)x[1]) 
treatment <- sapply(strsplit(colnames(endo_hrv_bsa), split ="_"),function(x)x[2]) 
enviro<- sapply(strsplit(colnames(endo_hrv_bsa), split ="_"),function(x)x[3])  
design <- data.frame(row.names = colnames(endo_hrv_bsa), subject = subject, treatment = treatment)
#Do blind analysis
dds_hrv_bsa <- DESeqDataSetFromMatrix(
  countData = endo_hrv_bsa,
  colData = design,
  design = ~1)
dds_hrv_bsa_norm<-estimateSizeFactors(dds_hrv_bsa)
dds_hrv_bsa_norm<-counts(dds_hrv_bsa_norm, normalized=T)
dds_hrv_bsa_norm_log<-log2(dds_hrv_bsa_norm + 1)

#pdf("Coronovirus_genes_HRVpilot.pdf",width=4,height=5)
#BSA plots
dev.new(width=4,height=5)
par(mfrow=c(1,2),bty="l")
boxplot(dds_hrv_bsa_norm["ACE2",]~design$treatment,las=1,main="ACE2",
	col=c("grey","tomato"),names=F,ylab="Normalized expression")	
points(dds_hrv_bsa_norm["ACE2",]~design$treatment,pch=16,col=c("black"))
segments(x0=rep(1,5),y0=dds_hrv_bsa_norm["ACE2",which(design$treatment == "ctrl")],
	x1=rep(2,5),y1=dds_hrv_bsa_norm["ACE2",which(design$treatment == "hrv")])
#mtext(paste("p = ",formatC(t.test(mean_termGenes_il13_up[which(design$enviro=="bsa" & design$treatment=="hrv"),i],
#	mean_termGenes_il13_up[which(design$enviro=="bsa" & design$treatment=="ctrl"),i],alternative="greater",paired=T)$p.value,
#	format="e",digits=2),sep=""),side=3,line=0.2,at=2.5,adj=1,cex=1)
text(x=c(1,2),y=par()$usr[3]-0.20*(par()$usr[4]-par()$usr[3]),labels=c("CTRL","HRV"),srt=45,adj=0.8,xpd=T)
#1.29E-51 - up

boxplot(dds_hrv_bsa_norm["TMPRSS2",]~design$treatment,las=1,main="TMPRSS2",
	col=c("grey","tomato"),names=F,ylab="Normalized expression")	
points(dds_hrv_bsa_norm["TMPRSS2",]~design$treatment,pch=16,col=c("black"))
segments(x0=rep(1,5),y0=dds_hrv_bsa_norm["TMPRSS2",which(design$treatment == "ctrl")],
	x1=rep(2,5),y1=dds_hrv_bsa_norm["TMPRSS2",which(design$treatment == "hrv")])
#mtext(paste("p = ",formatC(t.test(mean_termGenes_il13_up[which(design$enviro=="bsa" & design$treatment=="hrv"),i],
#	mean_termGenes_il13_up[which(design$enviro=="bsa" & design$treatment=="ctrl"),i],alternative="greater",paired=T)$p.value,
#	format="e",digits=2),sep=""),side=3,line=0.2,at=2.5,adj=1,cex=1)
text(x=c(1,2),y=par()$usr[3]-0.20*(par()$usr[4]-par()$usr[3]),labels=c("CTRL","HRV"),srt=45,adj=0.8,xpd=T)
#0.08249 - down



######################## PLOTING COVID-19 GENES (IL-13)
#subset data
endo_ctrl_il13<-endo[,c(1:10)]

#Specify design
subject <- sapply(strsplit(colnames(endo_ctrl_il13), split ="_"),function(x)x[1]) 
treatment <- sapply(strsplit(colnames(endo_ctrl_il13), split ="_"),function(x)x[2]) 
enviro<- sapply(strsplit(colnames(endo_ctrl_il13), split ="_"),function(x)x[3]) 
design <- data.frame(row.names = colnames(endo_ctrl_il13), subject = subject, treatment = enviro)
#Do blind analysis
dds_ctrl_il13 <- DESeqDataSetFromMatrix(
  countData = endo_ctrl_il13,
  colData = design,
  design = ~1)
dds_ctrl_il13_norm<-estimateSizeFactors(dds_ctrl_il13)
dds_ctrl_il13_norm<-counts(dds_ctrl_il13_norm, normalized=T)
dds_ctrl_il13_norm_log<-log2(dds_ctrl_il13_norm + 1)


#pdf("Coronovirus_genes_IL13pilot.pdf",width=4,height=5)
#BSA plots
dev.new(width=4,height=5)
par(mfrow=c(1,2),bty="l")
boxplot(dds_ctrl_il13_norm["ACE2",]~design$treatment,las=1,main="ACE2",
	col=c("grey","tomato"),names=F,ylab="Normalized expression")	
points(dds_ctrl_il13_norm["ACE2",]~design$treatment,pch=16,col=c("black"))
segments(x0=rep(1,5),y0=dds_ctrl_il13_norm["ACE2",which(design$treatment == "bsa")],
	x1=rep(2,5),y1=dds_ctrl_il13_norm["ACE2",which(design$treatment == "il13")])
#mtext(paste("p = ",formatC(t.test(mean_termGenes_il13_up[which(design$enviro=="bsa" & design$treatment=="il13"),i],
#	mean_termGenes_il13_up[which(design$enviro=="bsa" & design$treatment=="bsa"),i],alternative="greater",paired=T)$p.value,
#	format="e",digits=2),sep=""),side=3,line=0.2,at=2.5,adj=1,cex=1)
text(x=c(1,2),y=par()$usr[3]-0.20*(par()$usr[4]-par()$usr[3]),labels=c("CTRL","IL-13"),srt=45,adj=0.8,xpd=T)
#0.0048 - down

boxplot(dds_ctrl_il13_norm["TMPRSS2",]~design$treatment,las=1,main="TMPRSS2",
	col=c("grey","tomato"),names=F,ylab="Normalized expression")	
points(dds_ctrl_il13_norm["TMPRSS2",]~design$treatment,pch=16,col=c("black"))
segments(x0=rep(1,5),y0=dds_ctrl_il13_norm["TMPRSS2",which(design$treatment == "bsa")],
	x1=rep(2,5),y1=dds_ctrl_il13_norm["TMPRSS2",which(design$treatment == "il13")])
#mtext(paste("p = ",formatC(t.test(mean_termGenes_il13_up[which(design$enviro=="bsa" & design$treatment=="il13"),i],
#	mean_termGenes_il13_up[which(design$enviro=="bsa" & design$treatment=="bsa"),i],alternative="greater",paired=T)$p.value,
#	format="e",digits=2),sep=""),side=3,line=0.2,at=2.5,adj=1,cex=1)
text(x=c(1,2),y=par()$usr[3]-0.20*(par()$usr[4]-par()$usr[3]),labels=c("CTRL","IL-13"),srt=45,adj=0.8,xpd=T)
#5.17e-09 - up

