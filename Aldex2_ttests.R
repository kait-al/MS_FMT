# BDV Aldex2 script 2020

#ALDEx2 is a Bioconductor package available here:
#http://www.bioconductor.org/packages/release/bioc/html/ALDEx2.html

# Please see the Bioconductor page for the manual and install instructions

library(ALDEx2)

################################################
#Taxonomy
################################################

#read in a table of counts data
d<-read.table("SV_table_01filt.txt", sep="\t", quote="", check.names=F, header=T, row.names=1, comment.char="")

#subsetting a table based on sample names (column names)
#donor samples vs. patient samples (at baseline)
#baseline is sample right before FMT 
#i.e. either at 0 months or 6 months, depending on Early or Late intervention respectively
D<-c("Donor.1.01.01.18","Donor.1.01.20.18","Donor.1.03.28.18","Donor.1.21.01.18","Donor.1.22.01.18","Donor.2.01.05.18","Donor.2.14.07.18","Donor.2.16.07.18","Donor.2.23.01.18","Donor.2.23.05.18")
P<-c("MK.FMT.001.0m","MK.FMT.002.6m","MK.FMT.003.6m","MK.FMT.004.0m","MK.FMT.007.6m","MK.FMT.008.6m","MK.FMT.010.0m","MK.FMT.011.0m","MK.FMT.012.6m")

#this will retain the same order as the lists above
# NOTE: ALDEx input must be a DATA FRAME *not* a matrix
aldex.in<-d[,c(D, P)]

#Make a vector of conditions. This must be in the same order and the same number as the columns (samples) of the input table (aldex.in)
conds<-c(rep("D", length(D)), rep("P", length(P)))

#get the clr values
#this is the main ALDEx function for all downstream analyses
#mc.samples=128 is often sufficient
x <- aldex.clr(aldex.in, conds, mc.samples=128, verbose=TRUE)

#perform t-test (both Welches and Wilcoxon, plus a Benjamini-Hochberg multiple test correction)
x.tt <- aldex.ttest(x, conds, paired.test=FALSE)

#estimate effect size and the within and between condition values
#include indiv. samples or not
x.effect <- aldex.effect(x, conds)

#merge the data
x.all <- data.frame(x.tt, x.effect)

#significant?
sig <- which(x.tt$we.eBH < 0.05)
#sig
// integer(0)
#therefore nothing significant between donors and patients at baseline by Welch's
#But many with large effect size (Bacteroides more prevalent in P than D, Prevotella more prevalent in D than P)

#write a .txt with results
write.table(x.all, file="aldex/aldex_DvP_tax.txt", sep="\t", quote=F, col.names=NA)

#read in a table of counts data
d<-read.table("SV_table_01filt.txt", sep="\t", quote="", check.names=F, header=T, row.names=1, comment.char="")

################################################
#additional comparisons? 

#subsetting a table based on sample names (column names)
#Patients receiving Donor1 baseline vs. post final FMT
#baseline is sample right before first FMT 
#i.e. either at 0 months or 6 months, depending on Early or Late intervention respectively
###IMPORTANT, set paired.test=TRUE!! 

B1<-c("MK.FMT.001.0m","MK.FMT.003.6m","MK.FMT.004.0m","MK.FMT.007.6m","MK.FMT.010.0m")
P1<-c("MK.FMT.001.5m","MK.FMT.003.11m","MK.FMT.004.5m","MK.FMT.007.10m","MK.FMT.010.5m")
aldex.in<-d[,c(B1, P1)]
conds<-c(rep("B1", length(B1)), rep("P1", length(P1)))
x <- aldex.clr(aldex.in, conds, mc.samples=128, verbose=TRUE)
x.tt <- aldex.ttest(x, conds, paired.test=TRUE)
x.effect <- aldex.effect(x, conds)
x.all <- data.frame(x.tt, x.effect)
write.table(x.all, file="aldex/aldex_B1vP1_tax.txt", sep="\t", quote=F, col.names=NA)


#Patients receiving Donor2 baseline vs. post final FMT
B2<-c("MK.FMT.002.6m","MK.FMT.008.6m","MK.FMT.011.0m","MK.FMT.012.6m")
P2<-c("MK.FMT.002.12m","MK.FMT.008.7m","MK.FMT.011.5m","MK.FMT.012.7m")
aldex.in<-d[,c(B2, P2)]
conds<-c(rep("B2", length(B2)), rep("P2", length(P2)))
x <- aldex.clr(aldex.in, conds, mc.samples=128, verbose=TRUE)
x.tt <- aldex.ttest(x, conds, paired.test=TRUE)
x.effect <- aldex.effect(x, conds)
x.all <- data.frame(x.tt, x.effect)
write.table(x.all, file="aldex/aldex_B2vP2_tax.txt", sep="\t", quote=F, col.names=NA)

#repeat with picrust2outputs
#write.table(x.all, file="aldex/aldex_B1vP1_EC.txt", sep="\t", quote=F, col.names=NA)
#write.table(x.all, file="aldex/aldex_B1vP1_KO.txt", sep="\t", quote=F, col.names=NA)
#write.table(x.all, file="aldex/aldex_B1vP1_paths.txt", sep="\t", quote=F, col.names=NA)
#write.table(x.all, file="aldex/aldex_B2vP2_EC.txt", sep="\t", quote=F, col.names=NA)
#write.table(x.all, file="aldex/aldex_B2vP2_KO.txt", sep="\t", quote=F, col.names=NA)
#write.table(x.all, file="aldex/aldex_B2vP2_paths.txt", sep="\t", quote=F, col.names=NA)


################################################
#EC Enzymes
################################################

d<-read.table("picrust2/output/EC_metagenome_out/pred_metagenome_unstrat_rounded.txt", sep="\t", quote="", check.names=F, header=T, row.names=1, comment.char="")
d.1 <- data.frame(d)

#same for everything else from aldex.in... to

#significant?
sig <- which(x.tt$we.eBH < 0.05)
#sig
// integer(0)
#therefore nothing significant between donors and patients at baseline by Welch's
#

#write a .txt with results
write.table(x.all, file="aldex/aldex_DvP_EC.txt", sep="\t", quote=F, col.names=NA)

################################################
#KO KEGG orthologs
################################################

d<-read.table("picrust2/output/KO_metagenome_out//pred_metagenome_unstrat_rounded.txt", sep="\t", quote="", check.names=F, header=T, row.names=1, comment.char="")
d.1 <- data.frame(d)

#same for everything else from aldex.in... to

#significant?
sig <- which(x.tt$we.eBH < 0.05)
#sig
// integer(0)
#therefore nothing significant between donors and patients at baseline by Welch's, but lots of high effect sizes
#

#write a .txt with results
write.table(x.all, file="aldex/aldex_DvP_KO.txt", sep="\t", quote=F, col.names=NA)

################################################
#Picrust2 Pathways
################################################

d<-read.table("picrust2/output/pathways_out/path_abun_unstrat_rounded.txt", sep="\t", quote="", check.names=F, header=T, row.names=1, comment.char="")
d.1 <- data.frame(d)

#same for everything else from aldex.in... to

#significant?
sig <- which(x.tt$we.eBH < 0.05)
#sig
// integer(0)
#therefore nothing significant between donors and patients at baseline by Welch's, lots of high effect sizes
#

#write a .txt with results
write.table(x.all, file="aldex/aldex_DvP_paths.txt", sep="\t", quote=F, col.names=NA)
