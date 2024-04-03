#Let's curate the single nucleotide variants to filter out false divergences and infer ancestral vs derived allele in the focal group
#this improves the accuracy of dN/dS and facilitates calculation of something the asymptotic MK
#it will require 2 things:
#1. scaffold, position, and base information for SNPs within the focal group
#2. scaffold, position, and base information for SNVs in the outgroup

#I've done a version of this analysis previously that parsed ancestral and derived polymorphisms across the genome
#both coding and intergenic, and updated the reference fasta for it
#that would still be useful for say demography, but given our focus on coding sequence variation
#we can focus on the subset (which is much smaller and thus easier to work with locally)

library("data.table")
setwd("C:/Users/andre/Documents/Documents/Documents/Salmon louse fast Z/salmon_louse_data")
polys<-as.data.frame(fread("Lsal.DNA_reseq.cohort.af.genic.txt",fill=T,header=T),stringsAsFactors=F)
#185,547 coding sequence polymorphisms
####
###
##
#PART I: SNP CURATION
#
##
###
####
#I'm pre-splitting HOM and HET calls in the outgroup because we do different things with these 
#starting with homozygous divergence calls. 
#Anytime we have a hom-div that does NOT coincide with a focal poly, it is a valid divergence
#the trickier case is when they coincide
#REF  ALT.Poly  ALT.Div
#  A       T       T   not a divergence
#  A       T       C   is a divergence

homdiv<-as.data.frame(fread("Lnor.hom.genic.txt",fill=T),stringsAsFactors=F)
colnames(homdiv)<-c("CHROM","POS","REF","ALT","Gene","Eff1hom","Eff2hom","Eff3hom")
#793,857 coding sequence variants

#homozygous calls should be easy, is the Div allele in the polys at all?
homcheck<-merge(polys,homdiv,by=c("CHROM","POS","REF","Gene"),all.y=T)
#there are 21,814 overlapping sites
colnames(homcheck)[5]<-"ALT.P"
colnames(homcheck)[10]<-"ALT.D"
homfail<-homcheck[which(homcheck$ALT.P==homcheck$ALT.D),]
nrow(homfail)
# 15,639 shared ancestral polymorphisms, discard
hompass<-homcheck[which(homcheck$ALT.P!=homcheck$ALT.D),]
#nrow(hompass)
#6,175 true divergences at polymorphic sites, keep

#now we just need to keep the non-coincident hom divs as well as the hompass verified divergences
homfilt<-(row.names(homdiv)%in%row.names(homfail))
#the above gives us matches to homfail, we want to reverse that
homfilt2<-!(homfilt)

homcure<-homdiv[homfilt2,]
#778,218 passing divergences

#let's doublecheck the math
nrow(homdiv)-nrow(homfail)
#778218
#it worked!

#now to deal with heterozygous calls

hetdiv<-as.data.frame(fread("Lnor.het.genic.txt",fill=T),stringsAsFactors=F)
colnames(hetdiv)<-c("CHROM","POS","REF","ALT","Gene","Eff1het","Eff2het")
#23,842

#then we have the following patterns to consider
##REF  ALT.Poly  ALT.Div
#  A         T        T
#  A         T        C
#  A         T        C,G


hetcheck<-merge(polys,hetdiv,by=c("CHROM","POS","REF","Gene"),all.y=T)
colnames(hetcheck)[5]<-"ALT.P"
colnames(hetcheck)[10]<-"ALT.D"

hetonly<-hetcheck[which(is.na(hetcheck$ALT.P)),]
#most of the het Div calls do not occur with polys
nrow(hetonly)
#23,349
#if a het div does not coincide with a poly, then the only options are 
#REF  ALT.Div
# A     (A,)T   This is an outgroup polymorphism and should be removed
# A       C,T   This should be very rare but can be counted as a divergence
#so right away we can subset out when ALT.D only lists one allele (since the other is the implied REF allele)
hetonly2<-hetonly[which(nchar(hetonly$ALT.D)>1),]
#there are 4,175 possible true het-only divergences
#of these, many include a deletion + a variant allele...
#I think we should exclude these because indels create frameshifts and which have outsized effects
#in other words does one syn/non-syn div matter to selection as much as the frameshifted partner allele?
hetonly2.5<-hetonly2[which(as.character(substr(hetonly2$ALT.D,1,1))!='*'),]
#that gets rid of when the deletion is the first scored allele
#for reasons beyond me, the second allele could instead be the deletion
#e.g. *,C is separate from C,*
hetonly3<-hetonly2.5[which(as.character(substr(hetonly2.5$ALT.D,3,3))!='*'),]
nrow(hetonly3)
#2427
#this leaves us with only true het allele calls (e.g. A,G)
#as a practical matter, we can't tell which of the two alleles is fixed and which is polymorphic 
#so we'll further subset to only consider cases where ALT1 and ALT2 have the same effect
#both syn or non-syn
#so that we aren't guessing as to the effect at that site
hetonly4<-hetonly3[which(hetonly3$Eff1het==hetonly3$Eff2het),]
nrow(hetonly4)
#that leaves 2252 het-only divergences that pass curation
hetfilt1<-(row.names(hetdiv)%in%row.names(hetonly4))
#the above gives us positive matches so we should be able to pass it directly
hetpass1<-hetdiv[hetfilt1,]
  
  
  
#now we need to turn to cases of coincident polymorphism and divergence
hetco<-hetcheck[which(!(is.na(hetcheck$ALT.P))),]
#there are only 493 coincident poly + hetdivs
#so what do we do with these?
#there are 3 patterns:
#REF  ALT.Poly  ALT.Div
#  A       T       T   not a divergence
#  A       A       T  [not captured in the merge but not a divergence]
#  A       T      C,G is a divergence
#  A       C      C,G is a divergence and a poly

#let's start easy
hetco2<-hetco[which(hetco$ALT.P!=hetco$ALT.D),]
nrow(hetco2)
#251
#so we are removing 242 ancestral polymorphisms
#many of the remainder have only a single called div allele, meaning they're heterozygous for the reference
hetco3<-hetco2[which(nchar(hetco2$ALT.D)>1),]
nrow(hetco3)
#142
#that removes another 109
#as with the hetonlys above, we will remove deletions as they create outsized fitness effects
hetco3.5<-hetco3[which(as.character(substr(hetco3$ALT.D,1,1))!='*'),]
hetco4<-hetco3.5[which(as.character(substr(hetco3.5$ALT.D,3,3))!='*'),]
nrow(hetco4)
#93...so we're dealing with very few variants here, but we're not even done yet, let's consider the following:
#     CHROM      POS   REF    Gene  ALT.P    AF              Eff1 Eff2 Eff3 ALT.D            
# HG994580.1  2855047   C  LSAA_547     G  0.75 synonymous_variant             A,G
#in the above the focal pop has C,G  and the div sample has A,G....I think parsimony says that 
#the G is ancestral, sorts into both lineages, then C poly in the focal and A poly in the outgroup 
#so this should not be a divergence
#I think we only want the furtive and elusive quadrallelic site focal: A,T out: C,G
#do we even have any of those?
hetco4.5<-hetco4[which(as.character(substr(hetco4$ALT.D,1,1))!=as.character(hetco4$ALT.P)),]
hetco5<-hetco4.5[which(as.character(substr(hetco4.5$ALT.D,3,3))!=as.character(hetco4.5$ALT.P)),]
#and that leaves us with
nrow(hetco5)
#21
#only 21 sites in the genome represented by all four bases in the focal and outgroup together
#but that's enough to have to check if any has differing effects for ALT 1 and 2
#so we'll have to remove these
hetco6<-hetco5[which(hetco5$Eff1het==hetco5$Eff2het),]
#that removes...1, leaving us with 20
hetfilt2<-(row.names(hetdiv)%in%row.names(hetco6))
#the above gives us positive matches so we should be able to pass it directly
hetpass2<-hetdiv[hetfilt2,]
hetcure<-rbind(hetpass1,hetpass2)

#now we can finally combine them
colnames(homcure)[6:8]<-c("Eff1","Eff2","Eff3")
colnames(hetcure)[6:7]<-c("Eff1","Eff2")
allpass<-rbind(homcure[,1:7],hetcure)

#some finally summary stats:
nrow(homcure)/nrow(homdiv)
#98.0% of homozygous outgroup variants pass
nrow(hetcure)/nrow(hetdiv)
#9.5% of heterozygous outgroup variants pass
nrow(allpass)/(nrow(homdiv)+nrow(hetdiv))
#95.4% of putative divergences pass overall

####
###
##
#PART II: counting Dn and Ds by gene
##
###
####
#as when using SNPeff's more convenient table, we want to subset out big effect mutations like start loss, stop gain, etc
mispass<-allpass[which(allpass$Eff1=="missense_variant"),]
#225,802
mistable<-as.data.frame(table(mispass$Gene),stringAsFactors=F)
colnames(mistable)<-c("Gene","Dn")
synpass<-allpass[which(allpass$Eff1=="synonymous_variant"),]
#538,300
syntable<-as.data.frame(table(synpass$Gene),stringAsFactors=F)
colnames(syntable)<-c("Gene","Ds")
curatedDn.Ds<-merge(mistable,syntable,by="Gene",all.x=T,all.y=T)
#but some genes will have zeros for one variant type or another, producing NAs
#we need to turn those into zeros
curatedDn.Ds[is.na(curatedDn.Ds)]<-0
#write.csv(curatedDn.Ds,"curatedDn.Ds.csv",row.names = F,quote=F)
#done
####
###
##
#PART III: subsetting low-frequency polymorphisms 
##
###
####
#now we can use similar logic to make our own asymptotic MK test inputs
#we have allele frequencies counted as AF currently, AF is the frequency of the ALT (non-ref) alleles 
#what we want for the purposes of an asymptotic MK is the frequencies of the derived alleles
#when the reference is ancestral AF = derived freq
#but the frequency should be 1 - AF if the reference allele is derived
#we can gain some insight via parsimony logic and the outgroup variants
#when we have a polymorphism call in the focal group but not outgroup,
#that means that the outgroup has the reference allele 
#from that we infer that the reference is ancestral and does not need a frequency update
#when these variants coincide in the focal and outgroup we have to evaluate them
#re-instantiating hom and het div to start fresh
homdiva<-as.data.frame(fread("Lnor.hom.genic.txt",fill=T),stringsAsFactors=F)
colnames(homdiva)<-c("CHROM","POS","REF","ALT","Gene","Eff1","Eff2","Eff3")
hetdiva<-as.data.frame(fread("Lnor.het.genic.txt",fill=T),stringsAsFactors=F)
colnames(hetdiva)<-c("CHROM","POS","REF","ALT","Gene","Eff1","Eff2")
#we can't combine hom and het div calls because a REF:A ALT:T ALT.D:T has different implications 
#for hom divs, that's an informative site, for het, it's not
#then use the trick above to preserve row indices while we make our subset filter
anscheck<-merge(polys,homdiva,by=c("CHROM","POS","REF","Gene"),all.x=T)
colnames(anscheck)[5]<-"ALT.P"
colnames(anscheck)[10]<-"ALT.D"
ansshare<-anscheck[which(!(is.na(anscheck$ALT.D))),]
#as math above showed there are 21,814 sites that could potentially require updating on ancestry
#probably the easiest case is this one:
#CAJNVT010000001.1 344454   A  LSAA_6     T  0.75       synonymous_variant         T
#where the ALT.P and ALT.D have the same state, suggesting that the reference allele is derived
ansflip1<-ansshare[which(ansshare$ALT.P==ansshare$ALT.D),]
nrow(ansflip1)
#15,639
anssimp<-ansflip1[which(nchar(ansflip1$AF)<6),]
#the above is for multi-allelic polys in this case, but we would remove them it we did
#we don't have any in this case but this step is here to make the script more general
flipfilt<-(row.names(polys)%in%row.names(anssimp))
#the above gives us matches to ansflip
#i'm just creating a checkpoint dummy variable below
polytest<-polys
#we take the complimentary freq
AFreplace<-(1-as.numeric(polys$AF[flipfilt]))
polytest$AF[flipfilt]<-AFreplace
#finally we can start to aggregate
mispoly<-polytest[which(polytest$Eff1=="missense_variant"),]
#71,596 missense polymorphisms
mistablep<-as.data.frame(table(mispoly$Gene),stringAsFactors=F)
colnames(mistablep)<-c("Gene","Pn")
synpoly<-polytest[which(polytest$Eff1=="synonymous_variant"),]
#88,059
syntablep<-as.data.frame(table(synpoly$Gene),stringAsFactors=F)
colnames(syntablep)<-c("Gene","Ps")
Pn.Ps<-merge(mistablep,syntablep,by="Gene",all.x=T,all.y=T)
#that's out baseline...can we subset on frequency now?
#for a start, there are 3,314 SNPs that have a frequency of zero after the flip
#these are places the reference genome diverged from the wild population
#this represents only 1.8% of coding sequences SNP calls (citrus mealybug was way worse)
polytests<-polytest[which(nchar(polytest$AF)<6),]
#we're going to throw out the tri-allelic sites because they contain commas in the AF
#which makes them strings and not numbers
#doing so loses 3319 SNP sites, but we have plenty still
polytests$AF<-as.numeric(polytests$AF)
sort(unique(polytests$AF))
#0.000 0.125 0.167 0.250 0.333 0.375 0.500 0.625 0.667 0.750 0.833 0.875 1.000
#above is the range of allele frequencies
hist(polytests[which(polytests$Eff1=="synonymous_variant"),5])
hist(polytests[which(polytests$Eff1=="missense_variant"),5])
#the above code plots AF for non-syn and syn stites and is a crude version of the unfolded site frequency spectrum
#I say crude because there's some weird stratification and nearly empty bins
#I believe that happens when sample missingness creates an uncertain call in 1 sample, creating a different denominator
#for our purposes it should be good enough though
#note the missense histogram has bigger low-frequency hump (<0.2) and lower peaks at higher frequencies than the syn
#just as expected from a real SFS
#now we can just create poly counts for things with frequency cutoffs
#we've already done the full count
#how about >0.1, 0.2, 0.3, 0.4, 0.5, etc.?
#0.1
polytests1<-polytests[which(polytests$AF>0.1),]
#178,914
mispoly1<-polytests1[which(polytests1$Eff1=="missense_variant"),]
#68,895 missense polymorphisms
mistablep1<-as.data.frame(table(mispoly1$Gene),stringAsFactors=F)
colnames(mistablep1)<-c("Gene","Pn.1")
synpoly1<-polytests1[which(polytests1$Eff1=="synonymous_variant"),]
#85,575
syntablep1<-as.data.frame(table(synpoly1$Gene),stringAsFactors=F)
colnames(syntablep1)<-c("Gene","Ps.1")
Pn.Ps.1<-merge(mistablep1,syntablep1,by="Gene",all.x=T,all.y=T)
#0.2
polytests2<-polytests[which(polytests$AF>0.2),]
#105,245
mispoly2<-polytests2[which(polytests2$Eff1=="missense_variant"),]
# 40,099 missense polymorphisms
mistablep2<-as.data.frame(table(mispoly2$Gene),stringAsFactors=F)
colnames(mistablep2)<-c("Gene","Pn.2")
synpoly2<-polytests2[which(polytests2$Eff1=="synonymous_variant"),]
#50,737
syntablep2<-as.data.frame(table(synpoly2$Gene),stringAsFactors=F)
colnames(syntablep2)<-c("Gene","Ps.2")
Pn.Ps.2<-merge(mistablep2,syntablep2,by="Gene",all.x=T,all.y=T)
#0.3
polytests3<-polytests[which(polytests$AF>0.3),]
#84,904
mispoly3<-polytests3[which(polytests3$Eff1=="missense_variant"),]
#32,839 missense polymorphisms
mistablep3<-as.data.frame(table(mispoly3$Gene),stringAsFactors=F)
colnames(mistablep3)<-c("Gene","Pn.3")
synpoly3<-polytests3[which(polytests3$Eff1=="synonymous_variant"),]
#40,788
syntablep3<-as.data.frame(table(synpoly3$Gene),stringAsFactors=F)
colnames(syntablep3)<-c("Gene","Ps.3")
Pn.Ps.3<-merge(mistablep3,syntablep3,by="Gene",all.x=T,all.y=T)
#0.4
polytests4<-polytests[which(polytests$AF>0.4),]
#73,281
mispoly4<-polytests4[which(polytests4$Eff1=="missense_variant"),]
#28714 missense polymorphisms
mistablep4<-as.data.frame(table(mispoly4$Gene),stringAsFactors=F)
colnames(mistablep4)<-c("Gene","Pn.4")
synpoly4<-polytests4[which(polytests4$Eff1=="synonymous_variant"),]
#35132
syntablep4<-as.data.frame(table(synpoly4$Gene),stringAsFactors=F)
colnames(syntablep4)<-c("Gene","Ps.4")
Pn.Ps.4<-merge(mistablep4,syntablep4,by="Gene",all.x=T,all.y=T)
#0.5
polytests5<-polytests[which(polytests$AF>0.5),]
#62,839
mispoly5<-polytests5[which(polytests5$Eff1=="missense_variant"),]
#24,549 missense polymorphisms
mistablep5<-as.data.frame(table(mispoly5$Gene),stringAsFactors=F)
colnames(mistablep5)<-c("Gene","Pn.5")
synpoly5<-polytests5[which(polytests5$Eff1=="synonymous_variant"),]
#30,234
syntablep5<-as.data.frame(table(synpoly5$Gene),stringAsFactors=F)
colnames(syntablep5)<-c("Gene","Ps.5")
Pn.Ps.5<-merge(mistablep5,syntablep5,by="Gene",all.x=T,all.y=T)
#0.6
polytests6<-polytests[which(polytests$AF>0.6),]
#62839
mispoly6<-polytests6[which(polytests6$Eff1=="missense_variant"),]
#24549 missense polymorphisms
mistablep6<-as.data.frame(table(mispoly6$Gene),stringAsFactors=F)
colnames(mistablep6)<-c("Gene","Pn.6")
synpoly6<-polytests6[which(polytests6$Eff1=="synonymous_variant"),]
#30234
syntablep6<-as.data.frame(table(synpoly6$Gene),stringAsFactors=F)
colnames(syntablep6)<-c("Gene","Ps.6")
Pn.Ps.6<-merge(mistablep6,syntablep6,by="Gene",all.x=T,all.y=T)
#0.7
polytests7<-polytests[which(polytests$AF>0.7),]
#55858
mispoly7<-polytests7[which(polytests7$Eff1=="missense_variant"),]
#22070 missense polymorphisms
mistablep7<-as.data.frame(table(mispoly7$Gene),stringAsFactors=F)
colnames(mistablep7)<-c("Gene","Pn.7")
synpoly7<-polytests7[which(polytests7$Eff1=="synonymous_variant"),]
#26784
syntablep7<-as.data.frame(table(synpoly7$Gene),stringAsFactors=F)
colnames(syntablep7)<-c("Gene","Ps.7")
Pn.Ps.7<-merge(mistablep7,syntablep7,by="Gene",all.x=T,all.y=T)
#0.8
polytests8<-polytests[which(polytests$AF>0.8),]
#48534
mispoly8<-polytests8[which(polytests8$Eff1=="missense_variant"),]
#19562 missense polymorphisms
mistablep8<-as.data.frame(table(mispoly8$Gene),stringAsFactors=F)
colnames(mistablep8)<-c("Gene","Pn.8")
synpoly8<-polytests8[which(polytests8$Eff1=="synonymous_variant"),]
#23021
syntablep8<-as.data.frame(table(synpoly8$Gene),stringAsFactors=F)
colnames(syntablep8)<-c("Gene","Ps.8")
Pn.Ps.8<-merge(mistablep8,syntablep8,by="Gene",all.x=T,all.y=T)
#0.9
polytests9<-polytests[which(polytests$AF>0.9),]
#35796
mispoly9<-polytests9[which(polytests9$Eff1=="missense_variant"),]
#14880 missense polymorphisms
mistablep9<-as.data.frame(table(mispoly9$Gene),stringAsFactors=F)
colnames(mistablep9)<-c("Gene","Pn.9")
synpoly9<-polytests9[which(polytests9$Eff1=="synonymous_variant"),]
#16717
syntablep9<-as.data.frame(table(synpoly9$Gene),stringAsFactors=F)
colnames(syntablep9)<-c("Gene","Ps.9")
Pn.Ps.9<-merge(mistablep9,syntablep9,by="Gene",all.x=T,all.y=T)
####now we need to merge everything into one dataframe
base<-merge(curatedDn.Ds,Pn.Ps,by="Gene",all.x=T,all.y=T)
base1<-merge(base,Pn.Ps.1,by="Gene",all.x=T)
base2<-merge(base1,Pn.Ps.2,by="Gene",all.x=T)
base3<-merge(base2,Pn.Ps.3,by="Gene",all.x=T)
base4<-merge(base3,Pn.Ps.4,by="Gene",all.x=T)
base5<-merge(base4,Pn.Ps.5,by="Gene",all.x=T)
base6<-merge(base5,Pn.Ps.6,by="Gene",all.x=T)
base7<-merge(base6,Pn.Ps.7,by="Gene",all.x=T)
base8<-merge(base7,Pn.Ps.8,by="Gene",all.x=T)
base9<-merge(base8,Pn.Ps.9,by="Gene",all.x=T)
base9[is.na(base9)]<-0
#write.csv(base9,"L.salmonis_curatedSNVs_and_freqs.csv",row.names=F,quote=F)

