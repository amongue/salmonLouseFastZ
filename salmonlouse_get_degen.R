setwd("/Users/Andrew/Documents/Salmon_louse/")

library(GenomicFeatures)

library(Biostrings)

library(BSgenome)

library(parallel)


#first create the genome "package", only need to do once, check the github page for details
#forgeBSgenomeDataPkg("LsalmonisSeedFile.txt") 


#after that is complete, proceed below:

cores.avail <- detectCores()

cores <- ifelse(cores.avail > 1, detectCores() - 1, 1 ) # setting the number of cores to use for parallel processing: all of them but 1. Set cores<-1 for serial execution.




#cores <- detectCores() - 1 

library("LsalmonisBSgenome")  # Load library of genome sequence

#GTFs can SUCK IT
taxdb <- makeTxDbFromGFF(file = "GCA_905330665.1_LSAL_IOA00_R2_genomic.gff", format="gff")  # Load library for TxDb


cdsby <- cdsBy(taxdb, by="tx", use.names=T)  # extract CDS into a GRangesList

cds_seqs <- extractTranscriptSeqs(LsalmonisBSgenome, cdsby) # use cdsby GrangesList to extract every CDS into a BiostringsList

#let's not go through that multiple times in debugging:
#saveRDS(cds_seqs, "salmonis_cds_seqs.rds")
#cds_seqs<-readRDS("salmonis_cds_seqs.rds")
#saveRDS(cdsby, "salmonis_cdsby.rds")
#cdsby<-readRDS("salmonis_cdsby.rds")

#############

# Filtering! currently the script does no filtering of transcripts, so it will output data for every transcript.

# This will produce substantial redundancies when there are multiple transcripts per gene. Fortunately, for the case

# of monarch OGS2, there are not multiple transcripts per locus for any gene (that I'm aware). So I've ignored this

# step at this point. One obvious solution is to select the longest transcript, and then randomly sample among equally

# long transcripts.

############



# import functions for parsing codon degeneracy.

source(file = "/Users/Andrew/Documents/map_coords_functions.R")



system.time(my.coords <- mclapply(1:length(cdsby), map.coords, mc.cores = cores) )

my.coords <- lapply(1:length(cdsby), map.coords)
#this isn't working...gotta decompose it to debug it

trxp <- names(cdsby)
	degenx <- cds.degen(my.cds = cds_seqs[1])
	scaffx <- get.coords(cds.gr = cdsby[[1]])
test<-cbind(trxp[2],degenx,scaffx)


#okay now we're getting somewhere, the error comes from when nrow(degenx) =/= nrow(scaffx)
#this is a problem in the very first element, but not in 2, 3, or 4....so how common is it?


degenx<- vector(mode = "list", length = length(cdsby))
scaffx<- vector(mode = "list", length = length(cdsby))
for(i in 1:length(cdsby))
{
	scaffx[[i]] <- get.coords(cds.gr = cdsby[[i]])
	degenx[[i]]<-cds.degen(my.cds = cds_seqs[i])
	}
	
plot(lengths(scaffx)/3,lengths(degenx)/7)
abline(a=0,b=1)


trxp <- names(cdsby)

mct<-cbind(scaffx,degenx)

len_check<-rep(-9,nrow(mct))
for(i in 1:nrow(mct))
{
	len_check[i]<-nrow(scaffx[[i]])==nrow(degenx[[i]])
}

#there are 32 problem genes
#let's split on this
len_check<-as.logical(len_check)

behaved<-mct[len_check,]
problem_genes<-mct[!len_check,]

bhv.table<-mapply(cbind, behaved[,1], behaved[,2], SIMPLIFY=F)

#saveRDS(bhv.table, "bhv.table.rds")
#FUCKING SUCCESSS!!!


#saveRDS(problem_genes, "problem_genes.rds")
#cds_seqs<-readRDS("problem_genes.rds")
#now let's reckon with the problem genes
len_geno<-sapply(problem_genes[,1], NROW)
len_gff<-sapply(problem_genes[,2], NROW)

diff_vec<-len_geno-len_gff

hist(diff_vec,xlim=c(-4,4),breaks=8)
#update 1





mctul<-do.call(rbind, degenx)

#let's table the full degen table, and hack it to get get #syn and nonsyn sites for L salmonis
map.coords  <- function(N) {
	trxp <- names(cdsby)
	degenx <- cds.degen(my.cds = cds_seqs[N])
	return(cbind(trxp,degenx)) # inlcuding 'N' for trouble shooting purposes
	}
	
my.coords<-cbind(trxp,degenx)
#saveRDS(my.coords, "citri_my.coords2.rds")
#my.coords<-readRDS("citri_my.coords2.rds")
#cds_seqs<-readRDS("citri_cds_seqs.rds")

my.coords.table <- do.call(rbind, degenx) # sequentialy rbind every element in the list.

#my.coords.table <- do.call(rbind, my.coords) # sequentialy rbind every element in the list.



colnames(my.coords.table)  <- c("coord.cds", "transcript", "codon", "codon.base", "AminoAcid", "codon.posxn", "degeneracy")  # update column names to be informatively human-readable




colnames(my.coords.table)  <- c("trxp","coord.cds",  "transcript","codon", "codon.base", "AA", "codon.pos", "degeneracy")  # update column names to be informatively human-readable
mct<-as.data.frame(my.coords.table)


save(mct, file= "Lsalmonis_0x4x_hack.Rdata")


dum<-rep(0,length(mct$degeneracy))
sumstart<-as.data.frame(cbind(mct$transcript,as.numeric(as.character(mct$degeneracy)),dum,dum),stringsAsFactors=F)
colnames(sumstart)<-c("Gene","Degeneracy","Nonsyn","Syn")
#sumstart<-sumstart[which(nchar(as.character(sumstart$Degeneracy))<2),]
#sumstart$Degeneracy<-as.numeric(as.character(sumstart$Degeneracy))
#now we're thinkin with vectors!
case1<-sumstart$Degeneracy=="0"
case2<-sumstart$Degeneracy=="2"
case3<-sumstart$Degeneracy=="3"
case4<-sumstart$Degeneracy=="4"
sumstart$Nonsyn[case1]<-1
sumstart$Nonsyn[case2]<-2/3
sumstart$Nonsyn[case3]<-1/3
#sumstart$Nonsyn[case4]<-0 this is unneccessary 
sumstart$Syn[case4]<-1
sumstart$Syn[case2]<-1/3
sumstart$Syn[case3]<-2/3
sumstart$Syn<-as.numeric(as.character(sumstart$Syn))
sumstart$Nonsyn<-as.numeric(as.character(sumstart$Nonsyn))
#now we need to sum by gene

genesyns<-aggregate(x = sumstart$Syn, by = list(sumstart$Gene), FUN = sum)
colnames(genesyns)<-c("Gene","Syn.sites")
genenons<-aggregate(x = sumstart$Nonsyn, by = list(sumstart$Gene), FUN = sum)
colnames(genenons)<-c("Gene","Non.sites")
wholecount<-merge(genesyns,genenons,by="Gene")
finalcount<-cbind(trxp,wholecount)
#write.csv(finalcount,"L_salmonis_synnoncount_v0.5.csv",row.names=F,quote=F)



# Use gzip for output text file, which tends to be quite big otherwise.

gz.out <- gzfile("Dplex_0x4x.txt.gz")

write.table(my.coords.table, file = gz.out, row.names = F, quote = F, sep = "\t")
