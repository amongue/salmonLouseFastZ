library(data.table)
setwd("C:/Users/andre/Documents/Documents/Documents/Salmon louse fast Z/salmon_louse_data/")
lp<-fread("Outgroup_filterstats_relevant.txt",header=T,stringsAsFactors=F)
lstart<-fread("L_salmonis_synnoncount_v0.5.csv")
 lstart <- lstart[order(lstart$GeneId , -abs(lstart$Syn.sites)), ] #sort by id and reverse of abs(value)
#Then: Remove items after the first within id groups
 lmid<-lstart[!duplicated(lstart$GeneId), ]     
 
 lhack<-merge(lp,lstart,by="GeneId") 
dn<-lhack$variants_effect_missense_variant/lhack$Non.sites
ds<-(lhack$variants_effect_stop_retained_variant + lhack$variants_effect_synonymous_variant)/lhack$Syn.sites

dnds<-dn/ds
boxplot(dnds,outline=F)

#let's finish setting up the salmon louse for analysis

#1. which chromosome is the Z?
mcov<-as.data.frame(fread("DNA_reseq1_100kb_cov_means.txt",header=T),stringsAsFactors=F)
fcov<-as.data.frame(fread("female_DNA_100kb_cov_means.txt",header=T),stringsAsFactors=F)
colnames(mcov)[3]<-"M"
colnames(fcov)[3]<-"F"
allcov<-merge(mcov,fcov,by=c("id","pos"))
M2F<-log2(allcov$M/allcov$F)
ac2<-as.data.frame(cbind(allcov,M2F),stringAsFactors=F)
#write.csv(ac2,"cov_calcs.csv")
ac2<-as.data.frame(fread("cov_calcs.csv"))
hist(ac2$M2F,breaks=50,xlim=c(-1,1.5))

#head(ac2[order(-M2F),])

boxplot(ac2$M2F~ac2$id,outline=F,xlim=c(70,85),cex.axis=0.5)
#HG994594.1 is the main Z


boxplot(ac2$M2F~ac2$id,outline=F)
myFun <- function(x) {
  c( mean = mean(x), median = median(x), 
    std = sd(x))
}

tapply(ac2$M2F,ac2$id, myFun)
#CAJNVT010000004.1
#CAJNVT010000004.1
#CAJNVT010000006.1
#CAJNVT010000008.1
#CAJNVT010000010.1
#CAJNVT010000011.1



#2. what scaffold does each gene belong to?
lstart<-fread("L_salmonis_synnoncount_v0.5.csv")
colnames(lstart)[1]<-"Gene"
ros<-fread("gene_to_scaffold_rosetta.txt")
ls2<-merge(ros,lstart,by="Gene",all.x=T)
Linkage<-rep("A",nrow(ls2))
ls3<-as.data.frame(cbind(ls2,Linkage))

ls3[which(ls3$Scaffold=="HG994594.1"),5]<-"Z"

ann<-as.data.frame(fread("DNAreseq_combined_filterstats_relevant.txt"),stringsAsFactors=F)
Pn<-ann$variants_effect_missense_variant+ann$variants_effect_stop_retained_variant
Ps<-ann$variants_effect_synonymous_variant
poly<-as.data.frame(cbind(ann$GeneId,Pn,Ps),stringsAsFactors=F)
colnames(poly)[1]<-"Gene"

ls4<-merge(ls3,poly,by="Gene",all.x=T)
#write.csv(ls4,"L_salmonis_mat_v1.csv")
ls4<-as.data.frame(fread("L_salmonis_mat_v1.csv"))

boxplot((ls4$Ps/ls4$Syn.sites)~ls4$Linkage,notch=T,outline=F)
boxplot((ls4$Pn/ls4$Non.sites)/(ls4$Ps/ls4$Syn.sites)~ls4$Linkage,notch=T,outline=F,ylab="Pn/Ps")
abline(h=0.2,lty=3)

#3. what's the sex-biased composition of the Z vs autos?
library("DESeq2")

#whether we do SPM or DEseq, we need to combine individual expression files into one
#the convention is Female1 = counts, F1 = TPM, etc
f1<-fread("female_RNAseq_1.tsv",header=T)
f2<-fread("female_RNAseq_2.tsv",header=T)
f3<-fread("female_RNAseq_3.tsv",header=T)
f4<-fread("female_RNAseq_4.tsv",header=T)
f5<-fread("female_RNAseq_5.tsv",header=T)
f6<-fread("female_RNAseq_6.tsv",header=T)

m1<-fread("male_RNAseq_1.tsv",header=T)
m2<-fread("male_RNAseq_2.tsv",header=T)
m3<-fread("male_RNAseq_3.tsv",header=T)
m4<-fread("male_RNAseq_4.tsv",header=T)
m5<-fread("male_RNAseq_5.tsv",header=T)
m6<-fread("male_RNAseq_6.tsv",header=T)

f12<-merge(f1,f2,by=c("target_id","length"))
f34<-merge(f3,f4,by=c("target_id","length"))
f56<-merge(f5,f6,by=c("target_id","length"))
f1234<-merge(f12,f34,by=c("target_id","length"))
fall<-merge(f1234,f56,by=c("target_id","length"))


m12<-merge(m1,m2,by=c("target_id","length"))
m34<-merge(m3,m4,by=c("target_id","length"))
m56<-merge(m5,m6,by=c("target_id","length"))
m1234<-merge(m12,m34,by=c("target_id","length"))
mall<-merge(m1234,m56,by=c("target_id","length"))

Lsam<-merge(fall,mall,by=c("target_id","length"))
#write.csv(Lsam,"L_salmonis_combined_RNA.csv")


fpsm<-read.csv("L_salmonis_male_RNA_FPKM_for_dosage.csv",header=T,sep=",")
fpsf<-read.csv("L_salmonis_female_RNA_FPKM_for_dosage.csv",header=T,sep=",")
linkage<-rep("Auto",nrow(fpsm))
fpsff<-cbind(fpsf,linkage)
fpsff[which(fpsff$Scaffold=="HG994594.1"),9]<-"Z"
fpsff[which(fpsff$Scaffold=="CAJNVT010000004.1"),9]<-"Z"
fpsff[which(fpsff$Scaffold=="CAJNVT010000006.1"),9]<-"Z"
fpsff[which(fpsff$Scaffold=="CAJNVT010000008.1"),9]<-"Z"
fpsff[which(fpsff$Scaffold=="CAJNVT010000010.1"),9]<-"Z"
fpsff[which(fpsff$Scaffold=="CAJNVT010000011.1"),9]<-"Z"
fpfmean<-as.data.frame(cbind(fpsff$Trx,rowMeans(fpsff[,3:8]),fpsff$linkage))
colnames(fpfmean)<-c("Trx","FPKM_mean_F","Linkage")
#fpfmean<-fpfmean[which(fpfmean$FPKM_mean>0),]
write.csv(fpfmean,"female_mean_FPKM.csv",quote=F,row.names = F)
fpfmean<-as.data.frame(read.csv("female_mean_FPKM.csv",header=T,stringsAsFactors = F))
boxplot(fpfmean$FPKM_mean~fpfmean$Linkage,notch=T,outline=F)



fpsmm<-cbind(fpsm,linkage)
fpsmm[which(fpsmm$Scaffold=="HG994594.1"),9]<-"Z"
fpsmm[which(fpsmm$Scaffold=="CAJNVT010000004.1"),9]<-"Z"
fpsmm[which(fpsmm$Scaffold=="CAJNVT010000006.1"),9]<-"Z"
fpsmm[which(fpsmm$Scaffold=="CAJNVT010000008.1"),9]<-"Z"
fpsmm[which(fpsmm$Scaffold=="CAJNVT010000010.1"),9]<-"Z"
fpsmm[which(fpsmm$Scaffold=="CAJNVT010000011.1"),9]<-"Z"
fpmmean<-as.data.frame(cbind(fpsmm$Trx,rowMeans(fpsmm[,3:8]),fpsmm$linkage))
colnames(fpmmean)<-c("Trx","FPKM_mean_M","Linkage")
#fpmmean<-fpmmean[which(fpmmean$FPKM_mean>0),]
write.csv(fpmmean,"male_mean_FPKM.csv",quote=F,row.names = F)
fpmmean<-as.data.frame(read.csv("male_mean_FPKM.csv",header=T,stringsAsFactors = F))

dosage<-merge(fpmmean,fpfmean,by=c("Trx","Linkage"))
dosa<-dosage[which(dosage$FPKM_mean_F>0 & dosage$FPKM_mean_M>0),]

par(mfrow=c(1,2))
boxplot(log2(dosa$FPKM_mean_F)~dosa$Linkage,notch=T,outline=F,ylab="",las=1,xlab="",cex.axis=1.7)
mtext("Female",side = 3, adj = 0.5, cex=2, line=1)
mtext(expression(log[2]~(FPKM)),side=2, adj=0.5, line=2, cex=1.7)
text(1.5,11.9,"***",cex=3)
boxplot(log2(dosa$FPKM_mean_M)~dosa$Linkage,notch=T,outline=F,ylab="",las=1,xlab="",cex.axis=1.7)
text(1.5,12,"**",cex=3)
mtext("Male",side = 3, adj = 0.5, cex=2, line=1)


wilcox.test(log2(dosa$FPKM_mean_F)~dosa$Linkage)
#W = 3101777, p-value = 0.01799
#Z < AA 
wilcox.test(log2(dosa$FPKM_mean_M)~dosa$Linkage)
#W = 3132811, p-value = 0.006427
#ZZ < AA

dosaa<-dosa[which((dosa$Linkage=="Auto")),]
dosaz<-dosa[which((dosa$Linkage=="Z")),]
wilcox.test(log2(dosaz$FPKM_mean_M),log2(dosaz$FPKM_mean_F))
#W = 77631, p-value = 0.2543
#Z=ZZ


###let's deconstruct to individual tissues
ldosf<-as.data.table(fread('L_salmonis_female_RNA_FPKM_for_dosage.csv', header = TRUE))
ldosfz<-ldosf[which(ldosf$Scaffold=="HG994594.1"),]
ldosfa<-ldosf[which(ldosf$Scaffold!="HG994594.1"),]

ldosm<-as.data.table(fread('L_salmonis_male_RNA_FPKM_for_dosage.csv', header = TRUE))
ldosmz<-ldosm[which(ldosm$Scaffold=="HG994594.1"),]
ldosma<-ldosm[which(ldosm$Scaffold!="HG994594.1"),]
#we did a poor job number and naming...but we want

#we want M4 vs F1 for antennae
#but we for individual differences we don't want dosage driven by just which genes are not expressed (=0 FPKM)
par(mfrow=c(1,2))
ldosfa1<-ldosfa[which(ldosfa$F1>0),]
ldosfz1<-ldosfz[which(ldosfz$F1>0),]
boxplot(log2(ldosfa1$F1),xlim=c(0.3,3),notch=T,outline=F,ylab="",las=1,xlab="",cex.axis=1.7)
boxplot(log2(ldosfz1$F1),add=T,at=2,notch=T,outline=F,ylab="",las=1,xlab="",cex.axis=1.7)
mtext("Female antennae",side = 3, adj = 0.5, cex=2, line=1)
mtext(expression(log[2]~(FPKM)),side=2, adj=0.5, line=2, cex=1.7)
wilcox.test(ldosfa1$F1,ldosfz1$F1)
#W = 1817837, p-value = 0.0006374
text(2,9,"**",cex=2)

ldosma4<-ldosma[which(ldosma$M4>0),]
ldosmz4<-ldosmz[which(ldosmz$M4>0),]
boxplot(log2(ldosma4$M4),xlim=c(0.3,3),notch=T,outline=F,ylab="",las=1,xlab="",cex.axis=1.7)
boxplot(log2(ldosmz4$M4),add=T,at=2,notch=T,outline=F,ylab="",las=1,xlab="",cex.axis=1.7)
mtext("Male antennae",side = 3, adj = 0.5, cex=2, line=1)
mtext(expression(log[2]~(FPKM)),side=2, adj=0.5, line=2, cex=1.7)
wilcox.test(ldosma4$M4,ldosmz4$M4)
#W = 2175984, p-value = 0.0241
text(2,11.1,"*",cex=2)


#M6 vs F6 for whole adult
par(mfrow=c(1,2))
ldosfa6<-ldosfa[which(ldosfa$F6>0),]
ldosfz6<-ldosfz[which(ldosfz$F6>0),]
boxplot(log2(ldosfa6$F6),xlim=c(0.3,3),notch=T,outline=F,ylab="",las=1,xlab="",cex.axis=1.7)
boxplot(log2(ldosfz6$F6),add=T,at=2,notch=T,outline=F,ylab="",las=1,xlab="",cex.axis=1.7)
mtext("whole Adult F",side = 3, adj = 0.5, cex=2, line=1)
mtext(expression(log[2]~(FPKM)),side=2, adj=0.5, line=2, cex=1.7)
wilcox.test(ldosfa6$F6,ldosfz6$F6)
#W = 2314237, p-value = 0.0658
text(2,11.1,".",cex=4)

ldosma6<-ldosma[which(ldosma$M6>0),]
ldosmz6<-ldosmz[which(ldosmz$M6>0),]
boxplot(log2(ldosma6$M6),xlim=c(0.3,3),notch=T,outline=F,ylab="",las=1,xlab="",cex.axis=1.7)
boxplot(log2(ldosmz6$M6),add=T,at=2,notch=T,outline=F,ylab="",las=1,xlab="",cex.axis=1.7)
mtext("Whole Adult M",side = 3, adj = 0.5, cex=2, line=1)
mtext(expression(log[2]~(FPKM)),side=2, adj=0.5, line=2, cex=1.7)
wilcox.test(ldosma6$M6,ldosmz6$M6)
#W = 2082727, p-value = 0.1186



#M5 vs F2 for testes vs ovaries
par(mfrow=c(1,2))
ldosfa2<-ldosfa[which(ldosfa$F2>0),]
ldosfz2<-ldosfz[which(ldosfz$F2>0),]
boxplot(log2(ldosfa2$F2),xlim=c(0.3,3),notch=T,outline=F,ylab="",las=1,xlab="",cex.axis=1.7)
boxplot(log2(ldosfz2$F2),add=T,at=2,notch=T,outline=F,ylab="",las=1,xlab="",cex.axis=1.7)
mtext("Ovaries",side = 3, adj = 0.5, cex=2, line=1)
mtext(expression(log[2]~(FPKM)),side=2, adj=0.5, line=2, cex=1.7)
wilcox.test(ldosfa2$F2,ldosfz2$F2)
#W = 1790610, p-value = 0.05704
text(2,11.6,".",cex=4)

ldosma5<-ldosma[which(ldosma$M5>0),]
ldosmz5<-ldosmz[which(ldosmz$M5>0),]
boxplot(log2(ldosma5$M5),xlim=c(0.3,3),notch=T,outline=F,ylab="",las=1,xlab="",cex.axis=1.7)
boxplot(log2(ldosmz5$M5),add=T,at=2,notch=T,outline=F,ylab="",las=1,xlab="",cex.axis=1.7)
mtext("Testes",side = 3, adj = 0.5, cex=2, line=1)
mtext(expression(log[2]~(FPKM)),side=2, adj=0.5, line=2, cex=1.7)
wilcox.test(ldosma5$M5,ldosmz5$M5)
#W = 2170698, p-value = 0.8946


#M1,2,3 vs F3,4,5 for preadults
#lazy subsetting for the last one
ldosfa<-ldosfa[which(ldosfa$F3>0),]
ldosfa<-ldosfa[which(ldosfa$F4>0),]
ldosfa<-ldosfa[which(ldosfa$F5>0),]
ldosfz<-ldosfz[which(ldosfz$F3>0),]
ldosfz<-ldosfz[which(ldosfz$F4>0),]
ldosfz<-ldosfz[which(ldosfz$F5>0),]
ldosma<-ldosma[which(ldosma$M1>0),]
ldosma<-ldosma[which(ldosma$M2>0),]
ldosma<-ldosma[which(ldosma$M3>0),]
ldosmz<-ldosmz[which(ldosmz$M1>0),]
ldosmz<-ldosmz[which(ldosmz$M2>0),]
ldosmz<-ldosmz[which(ldosmz$M3>0),]

par(mfrow=c(1,2))
boxplot(log2(rowSums(ldosfa[,5:7]/3)),xlim=c(0.3,3),notch=T,outline=F,ylab="",las=1,xlab="",cex.axis=1.7)
boxplot(log2(rowSums(ldosfz[,5:7]/3)),add=T,at=2,notch=T,outline=F,ylab="",las=1,xlab="",cex.axis=1.7)
mtext("Mean Preadult F",side = 3, adj = 0.5, cex=2, line=1)
mtext(expression(log[2]~(FPKM)),side=2, adj=0.5, line=2, cex=1.7)
wilcox.test(log2(rowSums(ldosfa[,5:7]/3)),log2(rowSums(ldosz[,5:7]/3)))
#W = 2220382, p-value = 5.148e-08
text(2,9.5,"***",cex=2)

boxplot(log2(rowSums(ldosma[,3:5]/3)),xlim=c(0.3,3),notch=T,outline=F,ylab="",las=1,xlab="",cex.axis=1.7)
boxplot(log2(rowSums(ldosmz[,3:5]/3)),add=T,at=2,notch=T,outline=F,ylab="",las=1,xlab="",cex.axis=1.7)
mtext("Mean Preadult M",side = 3, adj = 0.5, cex=2, line=1)
mtext(expression(log[2]~(FPKM)),side=2, adj=0.5, line=2, cex=1.7)
wilcox.test(log2(rowSums(ldosma[,3:5]/3)),log2(rowSums(ldosmz[,3:5]/3)))
text(2,9,"***",cex=2)
#W = 2205199, p-value = 2.093e-05



#let's try SPM
lspm<-fread('L_salmonis_combined_RNA_SPM2.0.csv', header = TRUE)
colnames(lspm)[1]<-"Trx"
#the problem is that we have alternative transcripts for some genes. 
#We don't have a great way to deal with this other than to pick one
#we'll go with longest = canonical
gtl<-fread("Gene_to_locus_rosetta.txt")
df<-as.data.frame(merge(gtl,lspm,by="Trx"),stringsAsFactors=F)

df<-df[order(df[,'Gene'],-df[,'length']),]
lspm <- df[!duplicated(df$Gene),]

f_all<-rowSums(lspm[,4:9])
m_all<-rowSums(lspm[,10:15])
spm_f<-(f_all^2/(f_all^2+m_all^2))
hist(spm_f)
SPM.F<-as.data.frame(cbind(lspm$Gene,spm_f),stringsAsFactors=F)
colnames(SPM.F)<-c("Gene","SPM.F")


ls5<-merge(ls4,SPM.F,by="Gene")
SPM.Bias<-rep("Unbiased",nrow(ls5))
ls6[which(ls6$SPM.F<0.2),9]<-"Male"
ls6[which(ls6$SPM.F>0.8),9]<-"Female"
ls6<-as.data.frame(cbind(ls5,SPM.Bias),stringsAsFactors=F)
#write.csv(ls6,"L_salmonis_mat_v2.1.csv",row.names=F)
#time to get to dN/dS
ln<-fread("filterstats.genes.L.norm.txt",header=T,stringsAsFactors=F)
colnames(ln)
Dn<-ln$variants_effect_missense_variant+ln$variants_effect_stop_retained_variant
Ds<-ln$variants_effect_synonymous_variant
div<-as.data.frame(cbind(ln$GeneId,Dn,Ds),stringsAsFactors=F)
colnames(div)[1]<-"Gene"
ls7<-merge(ls6,div,by="Gene",all.x=T)
ls8<-ls7[!duplicated(ls7$Gene), ]  
#write.csv(ls8,"L_salmonis_mat_v2.5.csv",row.names=F)
#let's install a switch to bring in the more curated SNV dataset
ls0<-as.data.frame(fread("L_salmonis_mat_v2.5_varless.csv",stringsAsFactors = F))
ls0.5<-as.data.frame(fread("L.salmonis_curatedSNVs_and_freqs.csv",stringsAsFactors = F))
ls6<-as.data.frame(merge(ls0,ls0.5,by="Gene",all.x=T))
ls6<-as.data.frame(fread("L_salmonis_mat_v2.5.csv"))
hist(ls6$SPM.F)
abline(v=c(0.2,0.8))


#expanding the Z franchise to unplaced scaffolds
#CAJNVT010000004.1
#CAJNVT010000006.1
#CAJNVT010000008.1
#CAJNVT010000010.1
#CAJNVT010000011.1

ls6[which(ls6$Scaffold=="CAJNVT010000004.1"),5]<-"Z"
ls6[which(ls6$Scaffold=="CAJNVT010000006.1"),5]<-"Z"
ls6[which(ls6$Scaffold=="CAJNVT010000008.1"),5]<-"Z"
ls6[which(ls6$Scaffold=="CAJNVT010000010.1"),5]<-"Z"
ls6[which(ls6$Scaffold=="CAJNVT010000011.1"),5]<-"Z"



table(ls6$Linkage,ls6$SPM.Bias)
table(ls6$Linkage,ls6$SPM.Bias)/rowSums(table(ls6$Linkage,ls6$SPM.Bias))


bias_test<-table(ls6$Linkage,ls6$SPM.Bias)
allchi<-chisq.test(bias_test)
allchi
round(allchi$residuals,3)
#data:  bias_test
#X-squared = 21.625, df = 2, p-value = 2.015e-05
#    Female   Male Unbiased
#A -0.350 -0.508    0.442
#Z  2.113  3.071   -2.675


boxplot(ls6$Ps/ls6$Syn.sites~ls6$Linkage,outline=F,notch=T)
wilcox.test(ls6$Ps/ls6$Syn.sites~ls6$Linkage)
#W = 2511553, p-value = 0.9527
#synonymous variation overall NOT lower on the Z, as expected if it has roughly equal effective pop size

par(mai=c(1,1,1,0.03))
boxplot((ls6$Pn/ls6$Non.sites)/(ls6$Ps/ls6$Syn.sites)~ls6$Linkage,outline=F,notch=T,
        las=1, ylab="",xlab="Linkage",col=c("grey67","goldenrod"),cex.axis=1.4,names=c("Autosomal","Z"),cex.lab=2)
mtext(text="pN/pS",side=2,line=3,adj=0.5,cex=2)
mtext(text="Linkage",side=1,line=3,adj=0.5,cex=2)
mtext(text="Salmon louse faster Z", side=3, line=2, adj=0.5,cex=2.5)
wilcox.test((ls6$Pn/ls6$Non.sites)/(ls6$Ps/ls6$Syn.sites)~ls6$Linkage)
#W = 2282453, p-value = 0.07516
#very marginally higher pN/pS on the Z


boxplot(ls6$Pn/ls6$Non.sites~ls6$Linkage,outline=F,notch=T)
wilcox.test(ls6$Pn/ls6$Non.sites~ls6$Linkage)
#W = 2384750, p-value = 0.1108 however, pN is not...

boxplot(ls6$Ds/ls6$Syn.sites~ls6$Linkage,outline=F,notch=T)
wilcox.test(ls6$Ds/ls6$Syn.sites~ls6$Linkage)
#W = 3009395, p-value = 1.082e-11 is lower on the Z

boxplot(ls6$Dn/ls6$Non.sites~ls6$Linkage,outline=F,notch=T)
wilcox.test(ls6$Dn/ls6$Non.sites~ls6$Linkage)
#as is dN, W = 2949314, p-value = 1.969e-09



par(mai=c(1,1,1,0.03))
boxplot((ls6$Dn/ls6$Non.sites)/(ls6$Ds/ls6$Syn.sites)~ls6$Linkage,outline=F,notch=T,
        las=1, ylab="",xlab="Linkage",col=c("grey67","goldenrod"),cex.axis=1.4,names=c("Autosomal","Z"),cex.lab=2)
mtext(text="dN/dS",side=2,line=3.5,adj=0.5,cex=2)
mtext(text="Linkage",side=1,line=3,adj=0.5,cex=2)
mtext(text="Salmon louse faster Z", side=3, line=2, adj=0.5,cex=2.5)
text(x=1.5,y=0.4,"**",cex=4)

wilcox.test((ls6$Dn/ls6$Non.sites)/(ls6$Ds/ls6$Syn.sites)~ls6$Linkage)
#W = 588745, p-value = 0.002827
#a faster Z!

la<-ls6[which(ls6$Linkage=="A"),]
la$SPM.Bias<-factor(la$SPM.Bias,levels=c("Male","Unbiased","Female"))
lz<-ls6[which(ls6$Linkage=="Z"),]
lz$SPM.Bias<-factor(lz$SPM.Bias,levels=c("Male","Unbiased","Female"))


kruskal.test((lz$Dn/lz$Non.sites)/(lz$Ds/lz$Syn.sites)~lz$SPM.Bias)
#Kruskall-Wallis chi-squared = 0.071851, df = 2, p-value = 0.9647

kruskal.test((la$Dn/la$Non.sites)/(la$Ds/la$Syn.sites)~la$SPM.Bias)
#Kruskal-Wallis chi-squared = 22.222, df = 2, p-value = 1.495e-05
dt<-dunnTest((la$Dn/la$Non.sites)/(la$Ds/la$Syn.sites)~la$SPM.Bias)
#Comparison          Z      P.unadj        P.adj
#Female - Male   3.601789   3.160350e-04 6.320700e-04
#Female - Unbiased 4.652518 3.279054e-06 9.837161e-06
#Male - Unbiased 1.270899   2.037645e-01 2.037645e-01
#in other words F > UB = M


#Figure 1
layout(matrix(c(1,2,2,1,3,3),nrow=2,byrow=T),widths=c(2.3,3))
par(mai=c(1,1,0.8,0.03))
boxplot((ls6$Dn/ls6$Non.sites)/(ls6$Ds/ls6$Syn.sites)~ls6$Linkage,outline=F,notch=T,
        las=1, ylab="",xlab="",col=c("grey67","goldenrod"),cex.axis=1.9,names=c("Autosomal","Z"),cex.lab=2,ylim=c(0,0.4))
mtext(text="dN/dS",side=2,line=4,adj=0.5,cex=2)
mtext(text="Linkage",side=1,line=4,adj=0.5,cex=2)
mtext(text="Salmon louse faster Z", side=3, line=2, adj=0.5,cex=1.8)
text(x=2,y=0.37,"**",cex=4.3)
par(mai=c(0.05,1,0.8,0.04))
boxplot((lz$Dn/lz$Non.sites)/(lz$Ds/lz$Syn.sites)~lz$SPM.Bias,outline=F,
        notch=T,las=1,col=c("dodgerblue","grey32","red3"),ylab="",xlab="",names=c('','',''),
        xaxt='n',cex.axis=2.2,ylim=c(0,0.52))
text(x=1.5,y=0.48,"Z chromosome",cex=3.8)
mtext("Sex-biased gene evolution",side=3, adj=0.5,line=2,cex=2)
par(mai=c(0.7,1,0.3,0.04))
boxplot((la$Dn/la$Non.sites)/(la$Ds/la$Syn.sites)~la$SPM.Bias,outline=F,notch=T,
        las=1,col=c("dodgerblue","grey32","red3"),ylab="",xlab="",cex.axis=2.2,ylim=c(0,0.45))
text(x=3,y=0.4,"**",cex=4.3)
text(x=1.5,y=0.4,"Autosomes",cex=3.8)
mtext("Sex-bias class", side=1, adj=0.5, line=3.4,cex=2)



wilcox.test((ls6$Pn/ls6$Non.sites)/(ls6$Ps/ls6$Syn.sites)~ls6$Linkage)
#W = 2282453, p-value = 0.07516

kruskal.test((lz$Pn/lz$Non.sites)/(lz$Ps/lz$Syn.sites)~lz$SPM.Bias)
#Kruskal-Wallis chi-squared = 7.0173, df = 2, p-value = 0.02994
library("FSA")
dunnTest((lz$Pn/lz$Non.sites)/(lz$Ps/lz$Syn.sites)~lz$SPM.Bias)
#Comparison                 Z      P.unadj       P.adj
#    Female - Male  2.6397741  0.008296131  0.02488839
#Female - Unbiased  2.1329251  0.032930874  0.06586175
#  Male - Unbiased -0.9784869  0.327833543  0.32783354
# so we have a case where for Pn/Ps on the Z F > M = UB

#what about the autos?

kruskal.test((la$Pn/la$Non.sites)/(la$Ps/la$Syn.sites)~la$SPM.Bias)
#Kruskal-Wallis chi-squared = 147.58, df = 2, p-value < 2.2e-16
dunnTest((la$Pn/la$Non.sites)/(la$Ps/la$Syn.sites)~la$SPM.Bias)
#         Comparison         Z      P.unadj        P.adj
#     Female - Male  7.572658 3.656635e-14 7.313271e-14
# Female - Unbiased 11.899228 1.194476e-32 3.583428e-32
#   Male - Unbiased  4.576156 4.735978e-06 4.735978e-06
#on the autos F > M > UB

##Figure 2
layout(matrix(c(1,2,2,1,3,3),nrow=2,byrow=T),widths=c(2.3,3))
par(mai=c(1,1,1,0.03))
boxplot((ls6$Pn/ls6$Non.sites)/(ls6$Ps/ls6$Syn.sites)~ls6$Linkage,outline=F,notch=T,
        las=1, ylab="",xlab="",col=c("grey67","goldenrod"),cex.axis=1.9,names=c("Autosomal","Z"),cex.lab=2)
mtext(text="pN/pS",side=2,line=3.8,adj=0.5,cex=2)
mtext(text="Linkage",side=1,line=4,adj=0.5,cex=2)
mtext(text="Within-species variation", side=3, line=2, adj=0.5,cex=1.8)

par(mai=c(0.05,1,0.8,0.04))
boxplot((lz$Pn/lz$Non.sites)/(lz$Ps/lz$Syn.sites)~lz$SPM.Bias,outline=F,
        notch=T,las=1,col=c("dodgerblue","grey32","red3"),ylab="",xlab="",names=c('','',''),
        xaxt='n',cex.axis=2.2)
text(x=1.5,y=1.2,"Z chromosome",cex=3.8)
text(x=3,y=0.46,"a",cex=4,col="white")
text(x=1,y=0.33,"b",cex=4,col="white")
text(x=2,y=0.33,"b",cex=4,col="white")
mtext("Sex-biased gene evolution",side=3, adj=0.5,line=2,cex=2)
par(mai=c(0.7,1,0.3,0.04))
boxplot((la$Pn/la$Non.sites)/(la$Ps/la$Syn.sites)~la$SPM.Bias,outline=F,notch=T,
        las=1,col=c("dodgerblue","grey32","red3"),ylab="",xlab="",cex.axis=2.2)
text(x=3,y=0.4,"a",cex=4,col="white")
text(x=1,y=0.35,"b",cex=4,col="white")
text(x=2,y=0.3,"c",cex=4,col="white")
text(x=1.5,y=1.4,"Autosomes",cex=3.8)
mtext("Sex-bias class", side=1, adj=0.5, line=3.4,cex=2)



NItgCalc<-function(dn,ds,pn,ps)
{
  unbiasedNI<-sum((ds*pn)/(ps+ds))/sum((ps*dn)/(ps+ds))
  return(unbiasedNI)
}

ls6a<-ls6[which(ls6$Linkage=="A"),]
ls6z<-ls6[which(ls6$Linkage=="Z"),]

lza<-ls6z[which((ls6z$Ds + ls6z$Ps)!=0),]
laa<-ls6a[which((ls6a$Ds + ls6a$Ps)!=0),]
za<-1-NItgCalc(lza$Dn,lza$Ds,lza$Pn,lza$Ps)
aa<-1-NItgCalc(laa$Dn,laa$Ds,laa$Pn,laa$Ps)
#overall alpha
#Z:     -0.3015723
#autos: -0.3389717
lzma<-lza[which(lza$SPM.Bias=="Male"),]
lzua<-lza[which(lza$SPM.Bias=="Unbiased"),]
lzfa<-lza[which(lza$SPM.Bias=="Female"),]

lama<-laa[which(laa$SPM.Bias=="Male"),]
laua<-laa[which(laa$SPM.Bias=="Unbiased"),]
lafa<-laa[which(laa$SPM.Bias=="Female"),]

zma<-1-NItgCalc(lzma$Dn,lzma$Ds,lzma$Pn,lzma$Ps)
zua<-1-NItgCalc(lzua$Dn,lzua$Ds,lzua$Pn,lzua$Ps)
zfa<-1-NItgCalc(lzfa$Dn,lzfa$Ds,lzfa$Pn,lzfa$Ps)

ama<-1-NItgCalc(lama$Dn,lama$Ds,lama$Pn,lama$Ps)
aua<-1-NItgCalc(laua$Dn,laua$Ds,laua$Pn,laua$Ps)
afa<-1-NItgCalc(lafa$Dn,lafa$Ds,lafa$Pn,lafa$Ps)

#split alphas
#       male     unbiased      female
#Z     -0.300    -0.281        -0.582
#autos -0.276    -0.358        -0.332



testor<-abs((za-aa))
size<-10000
bsstata<-rep(0,size)
mwhole<-rbind(lza,laa)
#here's the switch for old vs new poly-set
mwhole$pN<-mwhole$Pn
mwhole$pS<-mwhole$Ps
mwhole$dN<-mwhole$Dn
mwhole$dS<-mwhole$Ds
pa<-nrow(lza)
pb<-nrow(mwhole)
for(i in 1:size)
{
  #next we'll shuffle our combined geneset
  bsap<-mwhole[sample(nrow(mwhole),replace=F),]
  bsstata[i]<-abs((1-NItgCalc(bsap$dN[(pa+1):pb],bsap$dS[(pa+1):pb],bsap$pN[(pa+1):pb],bsap$pS[(pa+1):pb]))-(1-NItgCalc(bsap$dN[1:pa],bsap$dS[1:pa],bsap$pN[1:pa],bsap$pS[1:pa])))
}
bsdist<-sort(bsstata)
p_val<-1-min(which(bsdist > testor))/length(bsdist)
p_val
#0.6677


#autosomes:
testor<-abs((ama-aua))
size<-10000
bsstata<-rep(0,size)
mwhole<-rbind(lama,laua)
#here's the switch for old vs new poly-set
mwhole$pN<-mwhole$Pn
mwhole$pS<-mwhole$Ps
mwhole$dN<-mwhole$Dn
mwhole$dS<-mwhole$Ds
pa<-nrow(lama)
pb<-nrow(mwhole)
for(i in 1:size)
{
  #next we'll shuffle our combined geneset
  bsap<-mwhole[sample(nrow(mwhole),replace=F),]
  bsstata[i]<-abs((1-NItgCalc(bsap$dN[(pa+1):pb],bsap$dS[(pa+1):pb],bsap$pN[(pa+1):pb],bsap$pS[(pa+1):pb]))-(1-NItgCalc(bsap$dN[1:pa],bsap$dS[1:pa],bsap$pN[1:pa],bsap$pS[1:pa])))
}
bsdist<-sort(bsstata)
p_val<-1-min(which(bsdist > testor))/length(bsdist)
p_val
#0.0219 M > UB 

testor<-abs((ama-afa))
size<-10000
bsstata<-rep(0,size)
mwhole<-rbind(lama,lafa)
#here's the switch for old vs new poly-set
mwhole$pN<-mwhole$Pn
mwhole$pS<-mwhole$Ps
mwhole$dN<-mwhole$Dn
mwhole$dS<-mwhole$Ds
pa<-nrow(lama)
pb<-nrow(mwhole)
for(i in 1:size)
{
  #next we'll shuffle our combined geneset
  bsap<-mwhole[sample(nrow(mwhole),replace=F),]
  bsstata[i]<-abs((1-NItgCalc(bsap$dN[(pa+1):pb],bsap$dS[(pa+1):pb],bsap$pN[(pa+1):pb],bsap$pS[(pa+1):pb]))-(1-NItgCalc(bsap$dN[1:pa],bsap$dS[1:pa],bsap$pN[1:pa],bsap$pS[1:pa])))
}
bsdist<-sort(bsstata)
p_val<-1-min(which(bsdist > testor))/length(bsdist)
p_val
#0.3146

testor<-abs((aua-afa))
size<-10000
bsstata<-rep(0,size)
mwhole<-rbind(laua,lafa)
#here's the switch for old vs new poly-set
mwhole$pN<-mwhole$Pn
mwhole$pS<-mwhole$Ps
mwhole$dN<-mwhole$Dn
mwhole$dS<-mwhole$Ds
pa<-nrow(laua)
pb<-nrow(mwhole)
for(i in 1:size)
{
  #next we'll shuffle our combined geneset
  bsap<-mwhole[sample(nrow(mwhole),replace=F),]
  bsstata[i]<-abs((1-NItgCalc(bsap$dN[(pa+1):pb],bsap$dS[(pa+1):pb],bsap$pN[(pa+1):pb],bsap$pS[(pa+1):pb]))-(1-NItgCalc(bsap$dN[1:pa],bsap$dS[1:pa],bsap$pN[1:pa],bsap$pS[1:pa])))
}
bsdist<-sort(bsstata)
p_val<-1-min(which(bsdist > testor))/length(bsdist)
p_val
#0.5999

#and now the Z

testor<-abs((zma-zua))
size<-10000
bsstata<-rep(0,size)
mwhole<-rbind(lzma,lzua)
#here's the switch for old vs new poly-set
mwhole$pN<-mwhole$Pn
mwhole$pS<-mwhole$Ps
mwhole$dN<-mwhole$Dn
mwhole$dS<-mwhole$Ds
pa<-nrow(lzma)
pb<-nrow(mwhole)
for(i in 1:size)
{
  #next we'll shuffle our combined geneset
  bsap<-mwhole[sample(nrow(mwhole),replace=F),]
  bsstata[i]<-abs((1-NItgCalc(bsap$dN[(pa+1):pb],bsap$dS[(pa+1):pb],bsap$pN[(pa+1):pb],bsap$pS[(pa+1):pb]))-(1-NItgCalc(bsap$dN[1:pa],bsap$dS[1:pa],bsap$pN[1:pa],bsap$pS[1:pa])))
}
bsdist<-sort(bsstata)
p_val<-1-min(which(bsdist > testor))/length(bsdist)
p_val
#0.9374 

testor<-abs((zma-zfa))
size<-10000
bsstata<-rep(0,size)
mwhole<-rbind(lzma,lzfa)
#here's the switch for old vs new poly-set
mwhole$pN<-mwhole$Pn
mwhole$pS<-mwhole$Ps
mwhole$dN<-mwhole$Dn
mwhole$dS<-mwhole$Ds
pa<-nrow(lzma)
pb<-nrow(mwhole)
for(i in 1:size)
{
  #next we'll shuffle our combined geneset
  bsap<-mwhole[sample(nrow(mwhole),replace=F),]
  bsstata[i]<-abs((1-NItgCalc(bsap$dN[(pa+1):pb],bsap$dS[(pa+1):pb],bsap$pN[(pa+1):pb],bsap$pS[(pa+1):pb]))-(1-NItgCalc(bsap$dN[1:pa],bsap$dS[1:pa],bsap$pN[1:pa],bsap$pS[1:pa])))
}
bsdist<-sort(bsstata)
p_val<-1-min(which(bsdist > testor))/length(bsdist)
p_val
#0.4141

testor<-abs((zua-zfa))
size<-10000
bsstata<-rep(0,size)
mwhole<-rbind(lzua,lzfa)
#here's the switch for old vs new poly-set
mwhole$pN<-mwhole$Pn
mwhole$pS<-mwhole$Ps
mwhole$dN<-mwhole$Dn
mwhole$dS<-mwhole$Ds
pa<-nrow(lzua)
pb<-nrow(mwhole)
for(i in 1:size)
{
  #next we'll shuffle our combined geneset
  bsap<-mwhole[sample(nrow(mwhole),replace=F),]
  bsstata[i]<-abs((1-NItgCalc(bsap$dN[(pa+1):pb],bsap$dS[(pa+1):pb],bsap$pN[(pa+1):pb],bsap$pS[(pa+1):pb]))-(1-NItgCalc(bsap$dN[1:pa],bsap$dS[1:pa],bsap$pN[1:pa],bsap$pS[1:pa])))
}
bsdist<-sort(bsstata)
p_val<-1-min(which(bsdist > testor))/length(bsdist)
p_val
#0.4098

#finally just to be sure zfa vs afa
testor<-abs((afa-zfa))
size<-10000
bsstata<-rep(0,size)
mwhole<-rbind(lzfa,lzfa)
#here's the switch for old vs new poly-set
mwhole$pN<-mwhole$Pn
mwhole$pS<-mwhole$Ps
mwhole$dN<-mwhole$Dn
mwhole$dS<-mwhole$Ds
pa<-nrow(lzfa)
pb<-nrow(mwhole)
for(i in 1:size)
{
  #next we'll shuffle our combined geneset
  bsap<-mwhole[sample(nrow(mwhole),replace=F),]
  bsstata[i]<-abs((1-NItgCalc(bsap$dN[(pa+1):pb],bsap$dS[(pa+1):pb],bsap$pN[(pa+1):pb],bsap$pS[(pa+1):pb]))-(1-NItgCalc(bsap$dN[1:pa],bsap$dS[1:pa],bsap$pN[1:pa],bsap$pS[1:pa])))
}
bsdist<-sort(bsstata)
p_val<-1-min(which(bsdist > testor))/length(bsdist)
p_val


testor<-abs((aua-zua))
size<-10000
bsstata<-rep(0,size)
mwhole<-rbind(laua,lzua)
#here's the switch for old vs new poly-set
mwhole$pN<-mwhole$Pn
mwhole$pS<-mwhole$Ps
mwhole$dN<-mwhole$Dn
mwhole$dS<-mwhole$Ds
pa<-nrow(laua)
pb<-nrow(mwhole)
for(i in 1:size)
{
  #next we'll shuffle our combined geneset
  bsap<-mwhole[sample(nrow(mwhole),replace=F),]
  bsstata[i]<-abs((1-NItgCalc(bsap$dN[(pa+1):pb],bsap$dS[(pa+1):pb],bsap$pN[(pa+1):pb],bsap$pS[(pa+1):pb]))-(1-NItgCalc(bsap$dN[1:pa],bsap$dS[1:pa],bsap$pN[1:pa],bsap$pS[1:pa])))
}
bsdist<-sort(bsstata)
p_val<-1-min(which(bsdist > testor))/length(bsdist)
p_val
#0.5053


##let's try to calculate Ne
thetas<-as.data.frame(fread("Lsal.angsd.thetas.idx .pestPG.txt",header=T))

the<-thetas$tW/thetas$nSites
#HG994594.1 is the main Z
#expanding the Z franchise to unplaced scaffolds
#CAJNVT010000004.1
#CAJNVT010000006.1
#CAJNVT010000008.1
#CAJNVT010000010.1
#CAJNVT010000011.1

thetasz<-thetas[which(thetas$Chr=="HG994594.1"|thetas$Chr=="CAJNVT010000004.1"| thetas$Chr=="CAJNVT010000006.1"|thetas$Chr=="CAJNVT010000008.1"|thetas$Chr=="CAJNVT010000010.1"|thetas$Chr=="CAJNVT010000011.1"),]
thetasa<-thetas[which(thetas$Chr!="HG994594.1"& thetas$Chr!="CAJNVT010000004.1"& thetas$Chr!="CAJNVT010000006.1"& thetas$Chr!="CAJNVT010000008.1"& thetas$Chr!="CAJNVT010000010.1"&thetas$Chr!="CAJNVT010000011.1"),]

mean(thetasz$tW/thetasz$nSites)/mean(thetasa$tW/thetasa$nSites)
#0.9543701



#####updating graphs for the supplement
#separating out dN, dS, pN and pS
par(mfrow=c(2,2))
boxplot(ls6$Dn/ls6$Non.sites~ls6$Linkage,outline=F,notch=T,
        las=1, ylab="",xlab="",col=c("grey67","goldenrod"),cex.axis=1.9,names=c("Autosomal","Z"),cex.lab=2)
text(2,0.06,"dN",cex=3)
text(2,0.036,"*",cex=3.5)
text(2,0.042,"*",cex=3.5)
text(2,0.048,"*",cex=3.5)
boxplot(ls6$Ds/ls6$Syn.sites~ls6$Linkage,outline=F,notch=T,
        las=1, ylab="",xlab="",col=c("grey67","goldenrod"),cex.axis=1.9,names=c("Autosomal","Z"),cex.lab=2)
text(2,0.5,"dS",cex=3)
text(2,0.28,"*",cex=3.5)
text(2,0.34,"*",cex=3.5)
text(2,0.4,"*",cex=3.5)
boxplot(ls6$Pn/ls6$Non.sites~ls6$Linkage,outline=F,notch=T,
        las=1, ylab="",xlab="",col=c("grey67","goldenrod"),cex.axis=1.9,names=c("Autosomal","Z"),cex.lab=2)
text(1,0.025,"pN",cex=3)
boxplot(ls6$Ps/ls6$Syn.sites~ls6$Linkage,outline=F,notch=T,
        las=1, ylab="",xlab="",col=c("grey67","goldenrod"),cex.axis=1.9,names=c("Autosomal","Z"),cex.lab=2)
text(2,0.08,"pS",cex=3)



wilcox.test(ls6$Dn/ls6$Non.sites~ls6$Linkage)
#as is dN, W = 2949314, p-value = 1.969e-09

wilcox.test(ls6$Ds/ls6$Syn.sites~ls6$Linkage)
#W = 3009395, p-value = 1.082e-11 is lower on the Z


wilcox.test(ls6$Pn/ls6$Non.sites~ls6$Linkage)
#W = 2384750, p-value = 0.1108 however, pN is not...


wilcox.test(ls6$Ps/ls6$Syn.sites~ls6$Linkage)
#W = 2511553, p-value = 0.9527
#synonymous variation overall NOT lower on the Z, as expected if it has roughly equal effective pop size


#some basic coverage stats for the outgroup
covs<-as.data.table(fread("LnorHiFi.Lsal.sorted.cov",stringsAsFactors = F))
mean(covs$meandepth[1:15])
#3.3x for chromosomal scaffolds
sum(covs[which(covs$meandepth>3),3])/sum(covs$endpos)
#63%
sum(covs[which(covs$meandepth>1),3])/sum(covs$endpos)
#97%


#let's explore the relationship between alpha and allele freq on the Z and A




#let's play with this a bit


#split by Z or A
ls6a<-ls6[which(ls6$Linkage=="A"),]
ls6z<-ls6[which(ls6$Linkage=="Z"),]

#get median of per-gene alphas using different frequency cutoffs
alova<-1-((ls6a$Pn/ls6a$Ps)/(ls6a$Dn/ls6a$Ds))
ala<-median(alov[is.finite(alov)])
a1<-1-((ls6a$Pn.1/ls6a$Ps.1)/(ls6a$Dn/ls6a$Ds))
a1a<-median(a1[is.finite(a1)])
a2<-1-((ls6a$Pn.2/ls6a$Ps.2)/(ls6a$Dn/ls6a$Ds))
a2a<-median(a2[is.finite(a2)])
a3<-1-((ls6a$Pn.3/ls6a$Ps.3)/(ls6a$Dn/ls6a$Ds))
a3a<-median(a3[is.finite(a3)])
a4<-1-((ls6a$Pn.4/ls6a$Ps.4)/(ls6a$Dn/ls6a$Ds))
a4a<-median(a4[is.finite(a4)])
a5<-1-((ls6a$Pn.5/ls6a$Ps.5)/(ls6a$Dn/ls6a$Ds))
a5a<-median(a5[is.finite(a5)])
a6<-1-((ls6a$Pn.6/ls6a$Ps.6)/(ls6a$Dn/ls6a$Ds))
a6a<-median(a6[is.finite(a6)])
a7<-1-((ls6a$Pn.7/ls6a$Ps.7)/(ls6a$Dn/ls6a$Ds))
a7a<-median(a7[is.finite(a7)])
a8<-1-((ls6a$Pn.8/ls6a$Ps.8)/(ls6a$Dn/ls6a$Ds))
a8a<-median(a8[is.finite(a8)])
a9<-1-((ls6a$Pn.9/ls6a$Ps.9)/(ls6a$Dn/ls6a$Ds))
a9a<-median(a9[is.finite(a9)])

par(mfrow=c(1,1))

alova<-1-((ls6z$Pn/ls6z$Ps)/(ls6z$Dn/ls6z$Ds))
alz<-median(alov[is.finite(alov)])
a1<-1-((ls6z$Pn.1/ls6z$Ps.1)/(ls6z$Dn/ls6z$Ds))
a1z<-median(a1[is.finite(a1)])
a2<-1-((ls6z$Pn.2/ls6z$Ps.2)/(ls6z$Dn/ls6z$Ds))
a2z<-median(a2[is.finite(a2)])
a3<-1-((ls6z$Pn.3/ls6z$Ps.3)/(ls6z$Dn/ls6z$Ds))
a3z<-median(a3[is.finite(a3)])
a4<-1-((ls6z$Pn.4/ls6z$Ps.4)/(ls6z$Dn/ls6z$Ds))
a4z<-median(a4[is.finite(a4)])
a5<-1-((ls6z$Pn.5/ls6z$Ps.5)/(ls6z$Dn/ls6z$Ds))
a5z<-median(a5[is.finite(a5)])
a6<-1-((ls6z$Pn.6/ls6z$Ps.6)/(ls6z$Dn/ls6z$Ds))
a6z<-median(a6[is.finite(a6)])
a7<-1-((ls6z$Pn.7/ls6z$Ps.7)/(ls6z$Dn/ls6z$Ds))
a7z<-median(a7[is.finite(a7)])
a8<-1-((ls6z$Pn.8/ls6z$Ps.8)/(ls6z$Dn/ls6z$Ds))
a8z<-median(a8[is.finite(a8)])
a9<-1-((ls6z$Pn.9/ls6z$Ps.9)/(ls6z$Dn/ls6z$Ds))
a9z<-median(a9[is.finite(a9)])

par(mfrow=c(1,1))
plot(c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),c(ala,a1a,a2a,a3a,a4a,a5a,a6a,a7a,a8a,a9a),
     ylim=c(-0.1,1),pch=15, col="grey67",cex=3)
segments(x=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),x1=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),
         y0=c(ala+sd(),a1a,a2a,a3a,a4a,a5a,a6a,a7a,a8a,a9a),
         points(c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),c(alz,a1z,a2z,a3z,a4z,a5z,a6z,a7z,a8z,a9z)
                ,pch=19, col="goldenrod",cex=3)
         abline(h=0,lty=3)
         
         
         
         
         
         #what if we just boxplot this?
         alova<-1-((ls6a$Pn/ls6a$Ps)/(ls6a$Dn/ls6a$Ds))
         ala<-(alova[is.finite(alova)])
         a1<-1-((ls6a$Pn.1/ls6a$Ps.1)/(ls6a$Dn/ls6a$Ds))
         a1a<-(a1[is.finite(a1)])
         a2<-1-((ls6a$Pn.2/ls6a$Ps.2)/(ls6a$Dn/ls6a$Ds))
         a2a<-(a2[is.finite(a2)])
         a3<-1-((ls6a$Pn.3/ls6a$Ps.3)/(ls6a$Dn/ls6a$Ds))
         a3a<-(a3[is.finite(a3)])
         a4<-1-((ls6a$Pn.4/ls6a$Ps.4)/(ls6a$Dn/ls6a$Ds))
         a4a<-(a4[is.finite(a4)])
         a5<-1-((ls6a$Pn.5/ls6a$Ps.5)/(ls6a$Dn/ls6a$Ds))
         a5a<-(a5[is.finite(a5)])
         a6<-1-((ls6a$Pn.6/ls6a$Ps.6)/(ls6a$Dn/ls6a$Ds))
         a6a<-(a6[is.finite(a6)])
         a7<-1-((ls6a$Pn.7/ls6a$Ps.7)/(ls6a$Dn/ls6a$Ds))
         a7a<-(a7[is.finite(a7)])
         a8<-1-((ls6a$Pn.8/ls6a$Ps.8)/(ls6a$Dn/ls6a$Ds))
         a8a<-(a8[is.finite(a8)])
         a9<-1-((ls6a$Pn.9/ls6a$Ps.9)/(ls6a$Dn/ls6a$Ds))
         a9a<-(a9[is.finite(a9)])
         
         #let's put some p-vals on the differences
         median(a1a)
         #-0.179
         median(a9a)
         #0.208
         wilcox.test(a1a,a9a)
         #W = 14136569, p-value < 2.2e-16
         #alpha increases with increasing stringency 
         
         
         par(mfrow=c(2,1))
         par(mai=c(0.02,1.1,1,0.01))
         boxplot(ala,notch=T,xlim=c(0,10),ylim=c(-2,1),las=1)
         mtext(text="Autosomes",side =3, line =1.5, adj=0.5,cex=2.2)
         mtext(text=expression(alpha),side=2,line=3,cex=3)
         boxplot(a1a,notch=T,add=T,at=2,las=1)
         boxplot(a2a,notch=T,add=T,at=3,las=1)
         boxplot(a3a,notch=T,add=T,at=4,las=1)
         boxplot(a4a,notch=T,add=T,at=5,las=1)
         boxplot(a5a,notch=T,add=T,at=6,las=1)
         boxplot(a6a,notch=T,add=T,at=7,las=1)
         boxplot(a7a,notch=T,add=T,at=8,las=1)
         boxplot(a8a,notch=T,add=T,at=9,las=1)
         boxplot(a9a,notch=T,add=T,at=10,las=1)
         abline(h=0,lty=3)
         alov<-1-((ls6z$Pn/ls6z$Ps)/(ls6z$Dn/ls6z$Ds))
         alz<-(alov[is.finite(alov)])
         a1<-1-((ls6z$Pn.1/ls6z$Ps.1)/(ls6z$Dn/ls6z$Ds))
         a1z<-(a1[is.finite(a1)])
         a2<-1-((ls6z$Pn.2/ls6z$Ps.2)/(ls6z$Dn/ls6z$Ds))
         a2z<-(a2[is.finite(a2)])
         a3<-1-((ls6z$Pn.3/ls6z$Ps.3)/(ls6z$Dn/ls6z$Ds))
         a3z<-(a3[is.finite(a3)])
         a4<-1-((ls6z$Pn.4/ls6z$Ps.4)/(ls6z$Dn/ls6z$Ds))
         a4z<-(a4[is.finite(a4)])
         a5<-1-((ls6z$Pn.5/ls6z$Ps.5)/(ls6z$Dn/ls6z$Ds))
         a5z<-(a5[is.finite(a5)])
         a6<-1-((ls6z$Pn.6/ls6z$Ps.6)/(ls6z$Dn/ls6z$Ds))
         a6z<-(a6[is.finite(a6)])
         a7<-1-((ls6z$Pn.7/ls6z$Ps.7)/(ls6z$Dn/ls6z$Ds))
         a7z<-(a7[is.finite(a7)])
         a8<-1-((ls6z$Pn.8/ls6z$Ps.8)/(ls6z$Dn/ls6z$Ds))
         a8z<-(a8[is.finite(a8)])
         a9<-1-((ls6z$Pn.9/ls6z$Ps.9)/(ls6z$Dn/ls6z$Ds))
         a9z<-(a9[is.finite(a9)])
         #let's put some p-vals on the differences
         median(a1z)
         #-0.198
         median(a9z)
         #-0.429
         wilcox.test(a1z,a9z)
         #W = 7071, p-value = 0.5728
         wilcox.test(a1a,a1z)
         #W = 566850, p-value = 0.8179
         wilcox.test(a9a,a9z)
         #W = 242396, p-value = 0.001378
         
         #out of curiosity
         wilcox.test(a1z,a2z)
         #still no difference here
         #W = 10667, p-value = 0.7166
         
         par(mai=c(1,1.1,0.8,0.01))
         boxplot(alz,notch=T,xlim=c(0,10),ylim=c(-2,1),col="goldenrod",las=1)
         mtext(text="Z chromosome",side =3, line =1.5, adj=0.5,cex=2.2)
         boxplot(a1z,notch=T,add=T,at=2,col="goldenrod",las=1)
         boxplot(a2z,notch=T,add=T,at=3,col="goldenrod",las=1)
         boxplot(a3z,notch=T,add=T,at=4,col="goldenrod",las=1)
         boxplot(a4z,notch=T,add=T,at=5,col="goldenrod",las=1)
         boxplot(a5z,notch=T,add=T,at=6,col="goldenrod",las=1)
         boxplot(a6z,notch=T,add=T,at=7,col="goldenrod",las=1)
         boxplot(a7z,notch=T,add=T,at=8,col="goldenrod",las=1)
         boxplot(a8z,notch=T,add=T,at=9,col="goldenrod",las=1)
         boxplot(a9z,notch=T,add=T,at=10,col="goldenrod",las=1)
         abline(h=0,lty=3)
         mtext(text=expression(alpha),side=2,line=3,cex=3)
         mtext(text="0.0",side=1,adj=0.1,las=0,cex=1.5,line=1)
         mtext(text="0.1",side=1,adj=0.2,las=0,cex=1.5,line=1)
         mtext(text="0.2",side=1,adj=0.3,las=0,cex=1.5,line=1)
         mtext(text="0.3",side=1,adj=0.4,las=0,cex=1.5,line=1)
         mtext(text="0.4",side=1,adj=0.5,las=0,cex=1.5,line=1)
         mtext(text="0.5",side=1,adj=0.6,las=0,cex=1.5,line=1)
         mtext(text="0.6",side=1,adj=0.7,las=0,cex=1.5,line=1)
         mtext(text="0.7",side=1,adj=0.8,las=0,cex=1.5,line=1)
         mtext(text="0.8",side=1,adj=0.9,las=0,cex=1.5,line=1)
         mtext(text="0.9",side=1,adj=1,las=0,cex=1.5,line=1)
         mtext("Polymorphism frequency >X",side=1,adj=0.5,line=3.4,cex=2.1)
         
         
         
         #####now impMKT-like
         #it can get more complicated, but as a first pass, what if we filter Pn but not Ps freqs
         #what if we just boxplot this?
         alova<-1-((ls6a$Pn/ls6a$Ps)/(ls6a$Dn/ls6a$Ds))
         ala<-(alova[is.finite(alova)])
         a1<-1-((ls6a$Pn.1/ls6a$Ps)/(ls6a$Dn/ls6a$Ds))
         a1a<-(a1[is.finite(a1)])
         a2<-1-((ls6a$Pn.2/ls6a$Ps)/(ls6a$Dn/ls6a$Ds))
         a2a<-(a2[is.finite(a2)])
         a3<-1-((ls6a$Pn.3/ls6a$Ps)/(ls6a$Dn/ls6a$Ds))
         a3a<-(a3[is.finite(a3)])
         a4<-1-((ls6a$Pn.4/ls6a$Ps)/(ls6a$Dn/ls6a$Ds))
         a4a<-(a4[is.finite(a4)])
         a5<-1-((ls6a$Pn.5/ls6a$Ps)/(ls6a$Dn/ls6a$Ds))
         a5a<-(a5[is.finite(a5)])
         a6<-1-((ls6a$Pn.6/ls6a$Ps)/(ls6a$Dn/ls6a$Ds))
         a6a<-(a6[is.finite(a6)])
         a7<-1-((ls6a$Pn.7/ls6a$Ps)/(ls6a$Dn/ls6a$Ds))
         a7a<-(a7[is.finite(a7)])
         a8<-1-((ls6a$Pn.8/ls6a$Ps)/(ls6a$Dn/ls6a$Ds))
         a8a<-(a8[is.finite(a8)])
         a9<-1-((ls6a$Pn.9/ls6a$Ps)/(ls6a$Dn/ls6a$Ds))
         a9a<-(a9[is.finite(a9)])
         
         #let's put some p-vals on the differences
         median(a1a)
         #-0.145
         median(a9a)
         #0.946
         wilcox.test(a1a,a9a)
         #W = 11843464, p-value < 2.2e-16
         #alpha increases with increasing stringency 
         
         
         par(mfrow=c(2,1))
         par(mai=c(0.02,1.1,1,0.01))
         boxplot(ala,notch=T,xlim=c(0,10),ylim=c(-2,1),las=1)
         mtext(text="Autosomes",side =3, line =1.5, adj=0.5,cex=2.2)
         mtext(text=expression(alpha),side=2,line=3,cex=3)
         boxplot(a1a,notch=T,add=T,at=2,las=1)
         boxplot(a2a,notch=T,add=T,at=3,las=1)
         boxplot(a3a,notch=T,add=T,at=4,las=1)
         boxplot(a4a,notch=T,add=T,at=5,las=1)
         boxplot(a5a,notch=T,add=T,at=6,las=1)
         boxplot(a6a,notch=T,add=T,at=7,las=1)
         boxplot(a7a,notch=T,add=T,at=8,las=1)
         boxplot(a8a,notch=T,add=T,at=9,las=1)
         boxplot(a9a,notch=T,add=T,at=10,las=1)
         abline(h=0,lty=3)
         alov<-1-((ls6z$Pn/ls6z$Ps)/(ls6z$Dn/ls6z$Ds))
         alz<-(alov[is.finite(alov)])
         a1<-1-((ls6z$Pn.1/ls6z$Ps)/(ls6z$Dn/ls6z$Ds))
         a1z<-(a1[is.finite(a1)])
         a2<-1-((ls6z$Pn.2/ls6z$Ps)/(ls6z$Dn/ls6z$Ds))
         a2z<-(a2[is.finite(a2)])
         a3<-1-((ls6z$Pn.3/ls6z$Ps)/(ls6z$Dn/ls6z$Ds))
         a3z<-(a3[is.finite(a3)])
         a4<-1-((ls6z$Pn.4/ls6z$Ps)/(ls6z$Dn/ls6z$Ds))
         a4z<-(a4[is.finite(a4)])
         a5<-1-((ls6z$Pn.5/ls6z$Ps)/(ls6z$Dn/ls6z$Ds))
         a5z<-(a5[is.finite(a5)])
         a6<-1-((ls6z$Pn.6/ls6z$Ps)/(ls6z$Dn/ls6z$Ds))
         a6z<-(a6[is.finite(a6)])
         a7<-1-((ls6z$Pn.7/ls6z$Ps)/(ls6z$Dn/ls6z$Ds))
         a7z<-(a7[is.finite(a7)])
         a8<-1-((ls6z$Pn.8/ls6z$Ps)/(ls6z$Dn/ls6z$Ds))
         a8z<-(a8[is.finite(a8)])
         a9<-1-((ls6z$Pn.9/ls6z$Ps)/(ls6z$Dn/ls6z$Ds))
         a9z<-(a9[is.finite(a9)])
         #let's put some p-vals on the differences
         median(a1z)
         #-0.1666667
         median(a9z)
         #0.7826797
         wilcox.test(a1z,a9z)
         #W = 6202, p-value = 3.052e-13
         wilcox.test(a1a,a1z)
         #W = 570584, p-value = 0.7522
         wilcox.test(a9a,a9z)
         #W = 636300, p-value = 0.002814
         
         #out of curiosity
         wilcox.test(a2a,a2z)
         #still no difference here
         #W = 603072, p-value = 0.1185
         median(a2a)
         #0.4368
         median(a2z)
         #0.3531
         
         par(mai=c(1,1.1,0.8,0.01))
         boxplot(alz,notch=T,xlim=c(0,10),ylim=c(-2,1),col="goldenrod",las=1)
         mtext(text="Z chromosome",side =3, line =1.5, adj=0.5,cex=2.2)
         boxplot(a1z,notch=T,add=T,at=2,col="goldenrod",las=1)
         boxplot(a2z,notch=T,add=T,at=3,col="goldenrod",las=1)
         boxplot(a3z,notch=T,add=T,at=4,col="goldenrod",las=1)
         boxplot(a4z,notch=T,add=T,at=5,col="goldenrod",las=1)
         boxplot(a5z,notch=T,add=T,at=6,col="goldenrod",las=1)
         boxplot(a6z,notch=T,add=T,at=7,col="goldenrod",las=1)
         boxplot(a7z,notch=T,add=T,at=8,col="goldenrod",las=1)
         boxplot(a8z,notch=T,add=T,at=9,col="goldenrod",las=1)
         boxplot(a9z,notch=T,add=T,at=10,col="goldenrod",las=1)
         abline(h=0,lty=3)
         mtext(text=expression(alpha),side=2,line=3,cex=3)
         mtext(text="0.0",side=1,adj=0.1,las=0,cex=1.5,line=1)
         mtext(text="0.1",side=1,adj=0.2,las=0,cex=1.5,line=1)
         mtext(text="0.2",side=1,adj=0.3,las=0,cex=1.5,line=1)
         mtext(text="0.3",side=1,adj=0.4,las=0,cex=1.5,line=1)
         mtext(text="0.4",side=1,adj=0.5,las=0,cex=1.5,line=1)
         mtext(text="0.5",side=1,adj=0.6,las=0,cex=1.5,line=1)
         mtext(text="0.6",side=1,adj=0.7,las=0,cex=1.5,line=1)
         mtext(text="0.7",side=1,adj=0.8,las=0,cex=1.5,line=1)
         mtext(text="0.8",side=1,adj=0.9,las=0,cex=1.5,line=1)
         mtext(text="0.9",side=1,adj=1,las=0,cex=1.5,line=1)
         mtext("Polymorphism frequency >X",side=1,adj=0.5,line=3.4,cex=2.1)
         
         #### so let's take our lesson from the impMKT-like and try to split genes by class
         ###what about per-class alphas now
         ls6a<-ls6[which(ls6$Linkage=="A"),]
         ls6am<-ls6a[which(ls6a$SPM.Bias=="Male"),]
         ls6au<-ls6a[which(ls6a$SPM.Bias=="Unbiased"),]
         ls6af<-ls6a[which(ls6a$SPM.Bias=="Female"),]
         
         ls6z<-ls6[which(ls6$Linkage=="Z"),]
         ls6zm<-ls6z[which(ls6z$SPM.Bias=="Male"),]
         ls6zu<-ls6z[which(ls6z$SPM.Bias=="Unbiased"),]
         ls6zf<-ls6z[which(ls6z$SPM.Bias=="Female"),]
         
         ##############################Female-biased#########################################
         alova<-1-((ls6af$Pn/ls6af$Ps)/(ls6af$Dn/ls6af$Ds))
         ala<-(alova[is.finite(alova)])
         a1<-1-((ls6af$Pn.1/ls6af$Ps)/(ls6af$Dn/ls6af$Ds))
         a1a<-(a1[is.finite(a1)])
         a2<-1-((ls6af$Pn.2/ls6af$Ps)/(ls6af$Dn/ls6af$Ds))
         a2a<-(a2[is.finite(a2)])
         a3<-1-((ls6af$Pn.3/ls6af$Ps)/(ls6af$Dn/ls6af$Ds))
         a3a<-(a3[is.finite(a3)])
         a4<-1-((ls6af$Pn.4/ls6af$Ps)/(ls6af$Dn/ls6af$Ds))
         a4a<-(a4[is.finite(a4)])
         a5<-1-((ls6af$Pn.5/ls6af$Ps)/(ls6af$Dn/ls6af$Ds))
         a5a<-(a5[is.finite(a5)])
         a6<-1-((ls6af$Pn.6/ls6af$Ps)/(ls6af$Dn/ls6af$Ds))
         a6a<-(a6[is.finite(a6)])
         a7<-1-((ls6af$Pn.7/ls6af$Ps)/(ls6af$Dn/ls6af$Ds))
         a7a<-(a7[is.finite(a7)])
         a8<-1-((ls6af$Pn.8/ls6af$Ps)/(ls6af$Dn/ls6af$Ds))
         a8a<-(a8[is.finite(a8)])
         a9<-1-((ls6af$Pn.9/ls6af$Ps)/(ls6af$Dn/ls6af$Ds))
         a9a<-(a9[is.finite(a9)])
         
         
         #let's put some p-vals on the differences
         median(a1a)
         #-0.1393
         median(a9a)
         #0.8796
         wilcox.test(a1a,a9a)
         #W = 29361, p-value < 2.2e-16
         #alpha increases with increasing stringency 
         
         
         par(mfrow=c(2,1))
         par(mai=c(0.02,1.1,1,0.01))
         boxplot(ala,notch=T,xlim=c(0,10),ylim=c(-2,1),las=1)
         mtext(text="Female-biased Autosomal",side =3, line =1.5, adj=0.5,cex=2.2)
         mtext(text=expression(alpha),side=2,line=3,cex=3)
         boxplot(a1a,notch=T,add=T,at=2,las=1)
         boxplot(a2a,notch=T,add=T,at=3,las=1)
         boxplot(a3a,notch=T,add=T,at=4,las=1)
         boxplot(a4a,notch=T,add=T,at=5,las=1)
         boxplot(a5a,notch=T,add=T,at=6,las=1)
         boxplot(a6a,notch=T,add=T,at=7,las=1)
         boxplot(a7a,notch=T,add=T,at=8,las=1)
         boxplot(a8a,notch=T,add=T,at=9,las=1)
         boxplot(a9a,notch=T,add=T,at=10,las=1)
         abline(h=0,lty=3)
         
         alov<-1-((ls6zf$Pn/ls6zf$Ps)/(ls6zf$Dn/ls6zf$Ds))
         alz<-(alov[is.finite(alov)])
         a1<-1-((ls6zf$Pn.1/ls6zf$Ps)/(ls6zf$Dn/ls6zf$Ds))
         a1z<-(a1[is.finite(a1)])
         a2<-1-((ls6zf$Pn.2/ls6zf$Ps)/(ls6zf$Dn/ls6zf$Ds))
         a2z<-(a2[is.finite(a2)])
         a3<-1-((ls6zf$Pn.3/ls6zf$Ps)/(ls6zf$Dn/ls6zf$Ds))
         a3z<-(a3[is.finite(a3)])
         a4<-1-((ls6zf$Pn.4/ls6zf$Ps)/(ls6zf$Dn/ls6zf$Ds))
         a4z<-(a4[is.finite(a4)])
         a5<-1-((ls6zf$Pn.5/ls6zf$Ps)/(ls6zf$Dn/ls6zf$Ds))
         a5z<-(a5[is.finite(a5)])
         a6<-1-((ls6zf$Pn.6/ls6zf$Ps)/(ls6zf$Dn/ls6zf$Ds))
         a6z<-(a6[is.finite(a6)])
         a7<-1-((ls6zf$Pn.7/ls6zf$Ps)/(ls6zf$Dn/ls6zf$Ds))
         a7z<-(a7[is.finite(a7)])
         a8<-1-((ls6zf$Pn.8/ls6zf$Ps)/(ls6zf$Dn/ls6zf$Ds))
         a8z<-(a8[is.finite(a8)])
         a9<-1-((ls6zf$Pn.9/ls6zf$Ps)/(ls6zf$Dn/ls6zf$Ds))
         a9z<-(a9[is.finite(a9)])
         
         #let's put some p-vals on the differences
         median(a1z)
         #-0.53846
         median(a9z)
         #0.75
         wilcox.test(a1z,a9z)
         #W = 6202, p-value = 3.052e-13
         wilcox.test(a1a,a1z)
         #W = 1937, p-value = 0.4023
         wilcox.test(a9a,a9z)
         #W = 636300, p-value = 0.002814
         
         #out of curiosity
         wilcox.test(a2a,a2z)
         #still no difference here
         #W = 1652, p-value = 0.9691
         median(a2a)
         #0.391
         median(a2z)
         #0.3333
         
         par(mai=c(1,1.1,0.8,0.01))
         boxplot(alz,notch=T,xlim=c(0,10),ylim=c(-2,1),col="goldenrod",las=1)
         mtext(text="Female-biased Z-linked",side =3, line =1.5, adj=0.5,cex=2.2)
         boxplot(a1z,notch=T,add=T,at=2,col="goldenrod",las=1)
         boxplot(a2z,notch=T,add=T,at=3,col="goldenrod",las=1)
         boxplot(a3z,notch=T,add=T,at=4,col="goldenrod",las=1)
         boxplot(a4z,notch=T,add=T,at=5,col="goldenrod",las=1)
         boxplot(a5z,notch=T,add=T,at=6,col="goldenrod",las=1)
         boxplot(a6z,notch=T,add=T,at=7,col="goldenrod",las=1)
         boxplot(a7z,notch=T,add=T,at=8,col="goldenrod",las=1)
         boxplot(a8z,notch=T,add=T,at=9,col="goldenrod",las=1)
         boxplot(a9z,notch=T,add=T,at=10,col="goldenrod",las=1)
         abline(h=0,lty=3)
         mtext(text=expression(alpha),side=2,line=3,cex=3)
         mtext(text="0.0",side=1,adj=0.1,las=0,cex=1.5,line=1)
         mtext(text="0.1",side=1,adj=0.2,las=0,cex=1.5,line=1)
         mtext(text="0.2",side=1,adj=0.3,las=0,cex=1.5,line=1)
         mtext(text="0.3",side=1,adj=0.4,las=0,cex=1.5,line=1)
         mtext(text="0.4",side=1,adj=0.5,las=0,cex=1.5,line=1)
         mtext(text="0.5",side=1,adj=0.6,las=0,cex=1.5,line=1)
         mtext(text="0.6",side=1,adj=0.7,las=0,cex=1.5,line=1)
         mtext(text="0.7",side=1,adj=0.8,las=0,cex=1.5,line=1)
         mtext(text="0.8",side=1,adj=0.9,las=0,cex=1.5,line=1)
         mtext(text="0.9",side=1,adj=1,las=0,cex=1.5,line=1)
         mtext("Polymorphism frequency >X",side=1,adj=0.5,line=3.4,cex=2.1)
         
         #####now unbiased ##############################
         alova<-1-((ls6au$Pn/ls6au$Ps)/(ls6au$Dn/ls6au$Ds))
         ala<-(alova[is.finite(alova)])
         a1<-1-((ls6au$Pn.1/ls6au$Ps)/(ls6au$Dn/ls6au$Ds))
         a1a<-(a1[is.finite(a1)])
         a2<-1-((ls6au$Pn.2/ls6au$Ps)/(ls6au$Dn/ls6au$Ds))
         a2a<-(a2[is.finite(a2)])
         a3<-1-((ls6au$Pn.3/ls6au$Ps)/(ls6au$Dn/ls6au$Ds))
         a3a<-(a3[is.finite(a3)])
         a4<-1-((ls6au$Pn.4/ls6au$Ps)/(ls6au$Dn/ls6au$Ds))
         a4a<-(a4[is.finite(a4)])
         a5<-1-((ls6au$Pn.5/ls6au$Ps)/(ls6au$Dn/ls6au$Ds))
         a5a<-(a5[is.finite(a5)])
         a6<-1-((ls6au$Pn.6/ls6au$Ps)/(ls6au$Dn/ls6au$Ds))
         a6a<-(a6[is.finite(a6)])
         a7<-1-((ls6au$Pn.7/ls6au$Ps)/(ls6au$Dn/ls6au$Ds))
         a7a<-(a7[is.finite(a7)])
         a8<-1-((ls6au$Pn.8/ls6au$Ps)/(ls6au$Dn/ls6au$Ds))
         a8a<-(a8[is.finite(a8)])
         a9<-1-((ls6au$Pn.9/ls6au$Ps)/(ls6au$Dn/ls6au$Ds))
         a9a<-(a9[is.finite(a9)])
         
         
         #let's put some p-vals on the differences
         median(a1a)
         #-0.1559379
         median(a9a)
         #1
         wilcox.test(a1a,a9a)
         #W = 7002109, p-value < 2.2e-16
         #alpha increases with increasing stringency 
         
         
         par(mfrow=c(2,1))
         par(mai=c(0.02,1.1,1,0.01))
         boxplot(ala,notch=T,xlim=c(0,10),ylim=c(-2,1),las=1)
         mtext(text="Unbiased Autosomal",side =3, line =1.5, adj=0.5,cex=2.2)
         mtext(text=expression(alpha),side=2,line=3,cex=3)
         boxplot(a1a,notch=T,add=T,at=2,las=1)
         boxplot(a2a,notch=T,add=T,at=3,las=1)
         boxplot(a3a,notch=T,add=T,at=4,las=1)
         boxplot(a4a,notch=T,add=T,at=5,las=1)
         boxplot(a5a,notch=T,add=T,at=6,las=1)
         boxplot(a6a,notch=T,add=T,at=7,las=1)
         boxplot(a7a,notch=T,add=T,at=8,las=1)
         boxplot(a8a,notch=T,add=T,at=9,las=1)
         boxplot(a9a,notch=T,add=T,at=10,las=1)
         abline(h=0,lty=3)
         
         alov<-1-((ls6zu$Pn/ls6zu$Ps)/(ls6zu$Dn/ls6zu$Ds))
         alz<-(alov[is.finite(alov)])
         a1<-1-((ls6zu$Pn.1/ls6zu$Ps)/(ls6zu$Dn/ls6zu$Ds))
         a1z<-(a1[is.finite(a1)])
         a2<-1-((ls6zu$Pn.2/ls6zu$Ps)/(ls6zu$Dn/ls6zu$Ds))
         a2z<-(a2[is.finite(a2)])
         a3<-1-((ls6zu$Pn.3/ls6zu$Ps)/(ls6zu$Dn/ls6zu$Ds))
         a3z<-(a3[is.finite(a3)])
         a4<-1-((ls6zu$Pn.4/ls6zu$Ps)/(ls6zu$Dn/ls6zu$Ds))
         a4z<-(a4[is.finite(a4)])
         a5<-1-((ls6zu$Pn.5/ls6zu$Ps)/(ls6zu$Dn/ls6zu$Ds))
         a5z<-(a5[is.finite(a5)])
         a6<-1-((ls6zu$Pn.6/ls6zu$Ps)/(ls6zu$Dn/ls6zu$Ds))
         a6z<-(a6[is.finite(a6)])
         a7<-1-((ls6zu$Pn.7/ls6zu$Ps)/(ls6zu$Dn/ls6zu$Ds))
         a7z<-(a7[is.finite(a7)])
         a8<-1-((ls6zu$Pn.8/ls6zu$Ps)/(ls6zu$Dn/ls6zu$Ds))
         a8z<-(a8[is.finite(a8)])
         a9<-1-((ls6zu$Pn.9/ls6zu$Ps)/(ls6zu$Dn/ls6zu$Ds))
         a9z<-(a9[is.finite(a9)])
         
         #let's put some p-vals on the differences
         median(a1z)
         #-0.04772727
         median(a9z)
         #0.7660494
         wilcox.test(a1z,a9z)
         #W = 2274, p-value = 4.786e-08
         wilcox.test(a1a,a1z)
         #W = 251365, p-value = 0.7946
         wilcox.test(a9a,a9z)
         #W = 288194, p-value = 0.02366
         
         #out of curiosity
         wilcox.test(a2a,a2z)
         #still no difference here
         #W = 269179, p-value = 0.3681
         median(a2a)
         #0.4337598
         median(a2z)
         #0.4020219
         
         par(mai=c(1,1.1,0.8,0.01))
         boxplot(alz,notch=T,xlim=c(0,10),ylim=c(-2,1),col="goldenrod",las=1)
         mtext(text="Unbiased Z-linked",side =3, line =1.5, adj=0.5,cex=2.2)
         boxplot(a1z,notch=T,add=T,at=2,col="goldenrod",las=1)
         boxplot(a2z,notch=T,add=T,at=3,col="goldenrod",las=1)
         boxplot(a3z,notch=T,add=T,at=4,col="goldenrod",las=1)
         boxplot(a4z,notch=T,add=T,at=5,col="goldenrod",las=1)
         boxplot(a5z,notch=T,add=T,at=6,col="goldenrod",las=1)
         boxplot(a6z,notch=T,add=T,at=7,col="goldenrod",las=1)
         boxplot(a7z,notch=T,add=T,at=8,col="goldenrod",las=1)
         boxplot(a8z,notch=T,add=T,at=9,col="goldenrod",las=1)
         boxplot(a9z,notch=T,add=T,at=10,col="goldenrod",las=1)
         abline(h=0,lty=3)
         mtext(text=expression(alpha),side=2,line=3,cex=3)
         mtext(text="0.0",side=1,adj=0.1,las=0,cex=1.5,line=1)
         mtext(text="0.1",side=1,adj=0.2,las=0,cex=1.5,line=1)
         mtext(text="0.2",side=1,adj=0.3,las=0,cex=1.5,line=1)
         mtext(text="0.3",side=1,adj=0.4,las=0,cex=1.5,line=1)
         mtext(text="0.4",side=1,adj=0.5,las=0,cex=1.5,line=1)
         mtext(text="0.5",side=1,adj=0.6,las=0,cex=1.5,line=1)
         mtext(text="0.6",side=1,adj=0.7,las=0,cex=1.5,line=1)
         mtext(text="0.7",side=1,adj=0.8,las=0,cex=1.5,line=1)
         mtext(text="0.8",side=1,adj=0.9,las=0,cex=1.5,line=1)
         mtext(text="0.9",side=1,adj=1,las=0,cex=1.5,line=1)
         mtext("Polymorphism frequency >X",side=1,adj=0.5,line=3.4,cex=2.1)
         
         ########################Male-biased#############################
         alova<-1-((ls6am$Pn/ls6am$Ps)/(ls6am$Dn/ls6am$Ds))
         ala<-(alova[is.finite(alova)])
         a1<-1-((ls6am$Pn.1/ls6am$Ps)/(ls6am$Dn/ls6am$Ds))
         a1a<-(a1[is.finite(a1)])
         a2<-1-((ls6am$Pn.2/ls6am$Ps)/(ls6am$Dn/ls6am$Ds))
         a2a<-(a2[is.finite(a2)])
         a3<-1-((ls6am$Pn.3/ls6am$Ps)/(ls6am$Dn/ls6am$Ds))
         a3a<-(a3[is.finite(a3)])
         a4<-1-((ls6am$Pn.4/ls6am$Ps)/(ls6am$Dn/ls6am$Ds))
         a4a<-(a4[is.finite(a4)])
         a5<-1-((ls6am$Pn.5/ls6am$Ps)/(ls6am$Dn/ls6am$Ds))
         a5a<-(a5[is.finite(a5)])
         a6<-1-((ls6am$Pn.6/ls6am$Ps)/(ls6am$Dn/ls6am$Ds))
         a6a<-(a6[is.finite(a6)])
         a7<-1-((ls6am$Pn.7/ls6am$Ps)/(ls6am$Dn/ls6am$Ds))
         a7a<-(a7[is.finite(a7)])
         a8<-1-((ls6am$Pn.8/ls6am$Ps)/(ls6am$Dn/ls6am$Ds))
         a8a<-(a8[is.finite(a8)])
         a9<-1-((ls6am$Pn.9/ls6am$Ps)/(ls6am$Dn/ls6am$Ds))
         a9a<-(a9[is.finite(a9)])
         
         
         #let's put some p-vals on the differences
         median(a1a)
         #-0.124841
         median(a9a)
         #0.9213788
         wilcox.test(a1a,a9a)
         #W = 386155, p-value < 2.2e-16
         #alpha increases with increasing stringency 
         
         
         par(mfrow=c(2,1))
         par(mai=c(0.02,1.1,1,0.01))
         boxplot(ala,notch=T,xlim=c(0,10),ylim=c(-2,1),las=1)
         mtext(text="Male-biased Autosomal",side =3, line =1.5, adj=0.5,cex=2.2)
         mtext(text=expression(alpha),side=2,line=3,cex=3)
         boxplot(a1a,notch=T,add=T,at=2,las=1)
         boxplot(a2a,notch=T,add=T,at=3,las=1)
         boxplot(a3a,notch=T,add=T,at=4,las=1)
         boxplot(a4a,notch=T,add=T,at=5,las=1)
         boxplot(a5a,notch=T,add=T,at=6,las=1)
         boxplot(a6a,notch=T,add=T,at=7,las=1)
         boxplot(a7a,notch=T,add=T,at=8,las=1)
         boxplot(a8a,notch=T,add=T,at=9,las=1)
         boxplot(a9a,notch=T,add=T,at=10,las=1)
         abline(h=0,lty=3)
         
         alov<-1-((ls6zm$Pn/ls6zm$Ps)/(ls6zm$Dn/ls6zm$Ds))
         alz<-(alov[is.finite(alov)])
         a1<-1-((ls6zm$Pn.1/ls6zm$Ps)/(ls6zm$Dn/ls6zm$Ds))
         a1z<-(a1[is.finite(a1)])
         a2<-1-((ls6zm$Pn.2/ls6zm$Ps)/(ls6zm$Dn/ls6zm$Ds))
         a2z<-(a2[is.finite(a2)])
         a3<-1-((ls6zm$Pn.3/ls6zm$Ps)/(ls6zm$Dn/ls6zm$Ds))
         a3z<-(a3[is.finite(a3)])
         a4<-1-((ls6zm$Pn.4/ls6zm$Ps)/(ls6zm$Dn/ls6zm$Ds))
         a4z<-(a4[is.finite(a4)])
         a5<-1-((ls6zm$Pn.5/ls6zm$Ps)/(ls6zm$Dn/ls6zm$Ds))
         a5z<-(a5[is.finite(a5)])
         a6<-1-((ls6zm$Pn.6/ls6zm$Ps)/(ls6zm$Dn/ls6zm$Ds))
         a6z<-(a6[is.finite(a6)])
         a7<-1-((ls6zm$Pn.7/ls6zm$Ps)/(ls6zm$Dn/ls6zm$Ds))
         a7z<-(a7[is.finite(a7)])
         a8<-1-((ls6zm$Pn.8/ls6zm$Ps)/(ls6zm$Dn/ls6zm$Ds))
         a8z<-(a8[is.finite(a8)])
         a9<-1-((ls6zm$Pn.9/ls6zm$Ps)/(ls6zm$Dn/ls6zm$Ds))
         a9z<-(a9[is.finite(a9)])
         
         #let's put some p-vals on the differences
         median(a1z)
         #-0.21212
         median(a9z)
         #0.789
         wilcox.test(a1z,a9z)
         #W = 738.5, p-value = 2.271e-05
         wilcox.test(a1a,a1z)
         #W = 38351, p-value = 0.5465
         wilcox.test(a9a,a9z)
         #W = 42366, p-value = 0.03764
         
         #out of curiosity
         wilcox.test(a2a,a2z)
         #still no difference here
         #W = 40946, p-value = 0.1356
         median(a2a)
         #0.464
         median(a2z)
         #0.2333
         
         par(mai=c(1,1.1,0.8,0.01))
         boxplot(alz,notch=T,xlim=c(0,10),ylim=c(-2,1),col="goldenrod",las=1)
         mtext(text="Male-biased Z-linked",side =3, line =1.5, adj=0.5,cex=2.2)
         boxplot(a1z,notch=T,add=T,at=2,col="goldenrod",las=1)
         boxplot(a2z,notch=T,add=T,at=3,col="goldenrod",las=1)
         boxplot(a3z,notch=T,add=T,at=4,col="goldenrod",las=1)
         boxplot(a4z,notch=T,add=T,at=5,col="goldenrod",las=1)
         boxplot(a5z,notch=T,add=T,at=6,col="goldenrod",las=1)
         boxplot(a6z,notch=T,add=T,at=7,col="goldenrod",las=1)
         boxplot(a7z,notch=T,add=T,at=8,col="goldenrod",las=1)
         boxplot(a8z,notch=T,add=T,at=9,col="goldenrod",las=1)
         boxplot(a9z,notch=T,add=T,at=10,col="goldenrod",las=1)
         abline(h=0,lty=3)
         mtext(text=expression(alpha),side=2,line=3,cex=3)
         mtext(text="0.0",side=1,adj=0.1,las=0,cex=1.5,line=1)
         mtext(text="0.1",side=1,adj=0.2,las=0,cex=1.5,line=1)
         mtext(text="0.2",side=1,adj=0.3,las=0,cex=1.5,line=1)
         mtext(text="0.3",side=1,adj=0.4,las=0,cex=1.5,line=1)
         mtext(text="0.4",side=1,adj=0.5,las=0,cex=1.5,line=1)
         mtext(text="0.5",side=1,adj=0.6,las=0,cex=1.5,line=1)
         mtext(text="0.6",side=1,adj=0.7,las=0,cex=1.5,line=1)
         mtext(text="0.7",side=1,adj=0.8,las=0,cex=1.5,line=1)
         mtext(text="0.8",side=1,adj=0.9,las=0,cex=1.5,line=1)
         mtext(text="0.9",side=1,adj=1,las=0,cex=1.5,line=1)
         mtext("Polymorphism frequency >X",side=1,adj=0.5,line=3.4,cex=2.1)
         
         
         ##what if we did something like the Lep faster-Z alpha figure with our impMKT-likes
         #so
         #Overall Autos vs Z
         # then auto vs Z for male, unbiased, female
         
         
         #we'll do across the board
         #all genes
         a<-1-((ls6a$Pn/ls6a$Ps)/(ls6a$Dn/ls6a$Ds))
         aa_full<-(a[is.finite(a2)])
         a<-1-((ls6z$Pn/ls6z$Ps)/(ls6z$Dn/ls6z$Ds))
         az_full<-(a[is.finite(a)])
         
         #we'll do the >0.2 cutoff Pn across the board
         #all genes
         a2<-1-((ls6a$Pn.2/ls6a$Ps)/(ls6a$Dn/ls6a$Ds))
         a2a_full<-(a2[is.finite(a2)])
         a2<-1-((ls6z$Pn.2/ls6z$Ps)/(ls6z$Dn/ls6z$Ds))
         a2z_full<-(a2[is.finite(a2)])
         
         #we'll do the >0.8 cutoff Pn across the board
         #all genes
         a8<-1-((ls6a$Pn.8/ls6a$Ps)/(ls6a$Dn/ls6a$Ds))
         a8a_full<-(a8[is.finite(a8)])
         a8<-1-((ls6z$Pn.8/ls6z$Ps)/(ls6z$Dn/ls6z$Ds))
         a8z_full<-(a8[is.finite(a8)])
         
         #we'll do the >0.9 cutoff Pn across the board
         #all genes
         a9<-1-((ls6a$Pn.9/ls6a$Ps)/(ls6a$Dn/ls6a$Ds))
         a9a_full<-(a9[is.finite(a9)])
         a9<-1-((ls6z$Pn.9/ls6z$Ps)/(ls6z$Dn/ls6z$Ds))
         a9z_full<-(a9[is.finite(a9)])
         
         #by sex-bias
         
         #simple alpha
         a<-1-((ls6am$Pn/ls6am$Ps)/(ls6am$Dn/ls6am$Ds))
         aa_mal<-(a[is.finite(a)])
         a<-1-((ls6zm$Pn/ls6zm$Ps)/(ls6zm$Dn/ls6zm$Ds))
         az_mal<-(a[is.finite(a)])
         
         a<-1-((ls6au$Pn/ls6au$Ps)/(ls6au$Dn/ls6au$Ds))
         aa_unb<-(a[is.finite(a)])
         a<-1-((ls6zu$Pn/ls6zu$Ps)/(ls6zu$Dn/ls6zu$Ds))
         az_unb<-(a[is.finite(a)])
         
         a<-1-((ls6af$Pn/ls6af$Ps)/(ls6af$Dn/ls6af$Ds))
         aa_fem<-(a[is.finite(a)])
         a<-1-((ls6zf$Pn/ls6zf$Ps)/(ls6zf$Dn/ls6zf$Ds))
         az_fem<-(a[is.finite(a)])
         
         a2<-1-((ls6af$Pn.2/ls6af$Ps)/(ls6af$Dn/ls6af$Ds))
         a2a_fem<-(a2[is.finite(a2)])
         a2<-1-((ls6zf$Pn.2/ls6zf$Ps)/(ls6zf$Dn/ls6zf$Ds))
         a2z_fem<-(a2[is.finite(a2)])
         
         
         #impMKT-like Pn>0.2
         a2<-1-((ls6am$Pn.2/ls6am$Ps)/(ls6am$Dn/ls6am$Ds))
         a2a_mal<-(a2[is.finite(a2)])
         a2<-1-((ls6zm$Pn.2/ls6zm$Ps)/(ls6zm$Dn/ls6zm$Ds))
         a2z_mal<-(a2[is.finite(a2)])
         
         a2<-1-((ls6au$Pn.2/ls6au$Ps)/(ls6au$Dn/ls6au$Ds))
         a2a_unb<-(a2[is.finite(a2)])
         a2<-1-((ls6zu$Pn.2/ls6zu$Ps)/(ls6zu$Dn/ls6zu$Ds))
         a2z_unb<-(a2[is.finite(a2)])
         
         a2<-1-((ls6af$Pn.2/ls6af$Ps)/(ls6af$Dn/ls6af$Ds))
         a2a_fem<-(a2[is.finite(a2)])
         a2<-1-((ls6zf$Pn.2/ls6zf$Ps)/(ls6zf$Dn/ls6zf$Ds))
         a2z_fem<-(a2[is.finite(a2)])
         
         a2<-1-((ls6af$Pn.2/ls6af$Ps)/(ls6af$Dn/ls6af$Ds))
         a2a_fem<-(a2[is.finite(a2)])
         a2<-1-((ls6zf$Pn.2/ls6zf$Ps)/(ls6zf$Dn/ls6zf$Ds))
         a2z_fem<-(a2[is.finite(a2)])
         
         
         #impMKT-like Pn>0.8
         a8<-1-((ls6am$Pn.8/ls6am$Ps)/(ls6am$Dn/ls6am$Ds))
         a8a_mal<-(a8[is.finite(a8)])
         a8<-1-((ls6zm$Pn.8/ls6zm$Ps)/(ls6zm$Dn/ls6zm$Ds))
         a8z_mal<-(a8[is.finite(a8)])
         
         a8<-1-((ls6au$Pn.8/ls6au$Ps)/(ls6au$Dn/ls6au$Ds))
         a8a_unb<-(a8[is.finite(a8)])
         a8<-1-((ls6zu$Pn.8/ls6zu$Ps)/(ls6zu$Dn/ls6zu$Ds))
         a8z_unb<-(a8[is.finite(a8)])
         
         a8<-1-((ls6af$Pn.8/ls6af$Ps)/(ls6af$Dn/ls6af$Ds))
         a8a_fem<-(a8[is.finite(a8)])
         a8<-1-((ls6zf$Pn.8/ls6zf$Ps)/(ls6zf$Dn/ls6zf$Ds))
         a8z_fem<-(a8[is.finite(a8)])
         
         a8<-1-((ls6af$Pn.8/ls6af$Ps)/(ls6af$Dn/ls6af$Ds))
         a8a_fem<-(a8[is.finite(a8)])
         a8<-1-((ls6zf$Pn.8/ls6zf$Ps)/(ls6zf$Dn/ls6zf$Ds))
         a8z_fem<-(a8[is.finite(a8)])
         
         
         #what about >0.9
         a9<-1-((ls6am$Pn.9/ls6am$Ps)/(ls6am$Dn/ls6am$Ds))
         a9a_mal<-(a9[is.finite(a9)])
         a9<-1-((ls6zm$Pn.9/ls6zm$Ps)/(ls6zm$Dn/ls6zm$Ds))
         a9z_mal<-(a9[is.finite(a9)])
         
         a9<-1-((ls6au$Pn.9/ls6au$Ps)/(ls6au$Dn/ls6au$Ds))
         a9a_unb<-(a9[is.finite(a9)])
         a9<-1-((ls6zu$Pn.9/ls6zu$Ps)/(ls6zu$Dn/ls6zu$Ds))
         a9z_unb<-(a9[is.finite(a9)])
         
         a9<-1-((ls6af$Pn.9/ls6af$Ps)/(ls6af$Dn/ls6af$Ds))
         a9a_fem<-(a9[is.finite(a9)])
         a9<-1-((ls6zf$Pn.9/ls6zf$Ps)/(ls6zf$Dn/ls6zf$Ds))
         a9z_fem<-(a9[is.finite(a9)])
         
         a9<-1-((ls6af$Pn.9/ls6af$Ps)/(ls6af$Dn/ls6af$Ds))
         a9a_fem<-(a9[is.finite(a9)])
         a9<-1-((ls6zf$Pn.9/ls6zf$Ps)/(ls6zf$Dn/ls6zf$Ds))
         a9z_fem<-(a9[is.finite(a9)])
         
         
         
         par(mfrow=c(3,1))
         par(mai=c(0.03,1,0.8,0.01))
         boxplot(aa_full,notch=T,xlim=c(0,9.2),ylim=c(-2,1.08),col="grey67",las=1, width=1.1)
         boxplot(az_full,notch=T,add=T,at=2,col="goldenrod",las=1,width=1.1)
         abline(h=0,lty=3)
         abline(v=3)
         boxplot(aa_mal,notch=T,add=T,at=3.7,col="grey67",las=1)
         boxplot(az_mal,notch=T,add=T,at=4.5,col="goldenrod",las=1)
         boxplot(aa_unb,notch=T,add=T,at=5.7,col="grey67",las=1)
         boxplot(az_unb,notch=T,add=T,at=6.5,col="goldenrod",las=1)
         boxplot(aa_fem,notch=T,add=T,at=7.7,col="grey67",las=1)
         boxplot(az_fem,notch=T,add=T,at=8.5,col="goldenrod",las=1)
         mtext(text=expression(alpha),side=2,line=3,cex=3,las=1)
         #mtext(text="Overall",side=1, line =1, adj=0.15, cex=2)
         #mtext(text="Male-biased",side=1, line =1, adj=0.45, cex=1.4)
         #mtext(text="Unbiased",side=1, line =1, adj=0.7, cex=1.4)
         #mtext(text="Female-biased",side=1, line =1, adj=0.96, cex=1.4)
         mtext(text="Considering all Non-synonymous polymorphisms (Pn)",side=3, line =1, adj=0.5, cex=2)
         
         
         par(mai=c(0.03,1,0.6,0.01))
         boxplot(a2a_full,notch=T,xlim=c(0,9.2),ylim=c(-1,1.08),col="grey67",las=1, width=1.1)
         boxplot(a2z_full,notch=T,add=T,at=2,col="goldenrod",las=1,width=1.1)
         abline(h=0,lty=3)
         abline(v=3)
         boxplot(a2a_mal,notch=T,add=T,at=3.7,col="grey67",las=1)
         boxplot(a2z_mal,notch=T,add=T,at=4.5,col="goldenrod",las=1)
         boxplot(a2a_unb,notch=T,add=T,at=5.7,col="grey67",las=1)
         boxplot(a2z_unb,notch=T,add=T,at=6.5,col="goldenrod",las=1)
         boxplot(a2a_fem,notch=T,add=T,at=7.7,col="grey67",las=1)
         boxplot(a2z_fem,notch=T,add=T,at=8.5,col="goldenrod",las=1)
         mtext(text=expression(alpha),side=2,line=3,cex=3,las=1)
         mtext(text="Considering Pn Frequency > 0.2",side=3, line =1, adj=0.5, cex=2)
         
         par(mai=c(0.38,1,0.6,0.01))
         boxplot(a9a_full,notch=T,xlim=c(0,9.2),ylim=c(-1,1.14),col="grey67",las=1, width=1.1)
         boxplot(a9z_full,notch=T,add=T,at=2,col="goldenrod",las=1,width=1.1)
         text(1,1.12,"**",cex=2.8)
         abline(h=0,lty=3)
         abline(v=3)
         boxplot(a9a_mal,notch=T,add=T,at=3.7,col="grey67",las=1)
         boxplot(a9z_mal,notch=T,add=T,at=4.5,col="goldenrod",las=1)
         text(3.7,1.12,"*",cex=2.8)
         boxplot(a9a_unb,notch=T,add=T,at=5.7,col="grey67",las=1)
         boxplot(a9z_unb,notch=T,add=T,at=6.5,col="goldenrod",las=1)
         text(5.7,1.12,"*",cex=2.8)
         boxplot(a9a_fem,notch=T,add=T,at=7.7,col="grey67",las=1)
         boxplot(a9z_fem,notch=T,add=T,at=8.5,col="goldenrod",las=1)
         mtext(text=expression(alpha),side=2,line=3,cex=3,las=1)
         mtext(text="Overall",side=1, line =1, adj=0.15, cex=2)
         mtext(text="Male-biased",side=1, line =1, adj=0.45, cex=1.4)
         mtext(text="Unbiased",side=1, line =1, adj=0.7, cex=1.4)
         mtext(text="Female-biased",side=1, line =1, adj=0.96, cex=1.4)
         mtext(text="Considering Pn Frequency > 0.9",side=3, line =1, adj=0.5, cex=2)
         
         
         
         #ancillary tests
         wilcox.test(aa_full,az_full)
         #W = 86262, p-value = 0.2172
         wilcox.test(aa_mal,az_mal)
         #W = 38563, p-value = 0.5
         wilcox.test(aa_unb,az_unb)
         #W = 250392, p-value = 0.7465
         wilcox.test(aa_fem,az_fem)
         #W = 1978.5, p-value = 0.3344
         
         #ancillary tests
         wilcox.test(a2a_full,a2z_full)
         #W = 603072, p-value = 0.1185
         wilcox.test(a2a_mal,a2z_mal)
         #W = 40946, p-value = 0.1356
         wilcox.test(a2a_unb,a2z_unb)
         #W = 269179, p-value = 0.3681
         wilcox.test(a2a_fem,a2z_fem)
         #W = 1652, p-value = 0.9691
         
         #ancillary tests
         wilcox.test(a9a_full,a9z_full)
         #W = 636300, p-value = 0.002814
         wilcox.test(a9a_mal,a9z_mal)
         #W = 42366, p-value = 0.03764
         wilcox.test(a9a_unb,a9z_unb)
         #W = 288194, p-value = 0.02366
         wilcox.test(a9a_fem,a9z_fem)
         #W = 1735.5, p-value = 0.8206