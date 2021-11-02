#R
library(zCompositions)
library(vegan)
library(ALDEx2)

################################################
#Taxonomy
################################################
d<-read.table("SV_table_01filt.txt", sep="\t", quote="", check.names=F, header=T, row.names=1, comment.char="")
d.1 <- data.frame(d)
#use only counts table, or remove taxonomy column

# Can't have zeroes. Bayesian-Multiplicative replacement of count zeros
d.czm <- cmultRepl(t(d.1),  label=0, method="CZM")

#Centre-log-ratio transform the data
d.clr<-apply(d.czm, 2, function(x){log2(x) - mean(log2(x))})

#calculate principal components
#features are COLUMNS
#samples as ROWS
d.pcx <- prcomp(d.clr)

m <- read.table("metadata.txt", header=T, sep='\t', comment.char = "", row.names=1)
m<-data.frame(m)

fit1 <- envfit(d.pcx, data.frame(m[, c("Patient","FMT","FMT_Donor","EarlyLate","Time"), drop = FALSE]), permutations = 999, na.rm = TRUE)
fit1
// 
// ***VECTORS
// 
//                PC1      PC2     r2 Pr(>r)    
// Patient    0.90698 -0.42117 0.4356  0.001 ***
// FMT_Donor  0.98115  0.19325 0.2506  0.001 ***
// Time      -0.95290  0.30328 0.0454  0.120    
// ---
// Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
// Permutation: free
// Number of permutations: 999
// 
// ***FACTORS:
// 
// Centroids:
//                PC1     PC2
// FMTFMT     -0.6247  2.5407
// FMTNo      -0.7065  0.9828
// EarlyLateE -2.1788 -8.8185
// EarlyLateL  0.3494  8.8227
// 
// Goodness of fit:
//               r2 Pr(>r)    
// FMT       0.0012  0.904    
// EarlyLate 0.1516  0.001 ***
// ---
// Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
// Permutation: free
// Number of permutations: 999
// 
// 10 observations deleted due to missingness
######Same but with patient only counts and metadata tables (no donors)
// ***VECTORS
// 
//                PC1      PC2     r2 Pr(>r)    
// Patient    0.82096  0.57098 0.4445  0.001 ***
// FMT_Donor  0.99600 -0.08938 0.2500  0.001 ***
// Time      -0.88574 -0.46419 0.0515  0.104    
// ---
// Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
// Permutation: free
// Number of permutations: 999
// 
// ***FACTORS:
// 
// Centroids:
//                PC1     PC2
// FMTFMT      0.1809 -0.8731
// FMTNo      -0.1658  0.8004
// EarlyLateE -3.3100  9.8137
// EarlyLateL  2.2267 -6.6020
// 
// Goodness of fit:
//               r2 Pr(>r)    
// FMT       0.0014  0.883    
// EarlyLate 0.1379  0.001 ***
// ---
// Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
// Permutation: free
// Number of permutations: 999
// 

#using patient only data (no donors)

mm <- model.matrix(~Time+FMT_Donor, m)

#This replaces rows that had NAs so aldex.glm doesn't produce an error in later steps
mm <- mm[match(rownames(m),rownames(mm)),]
rownames(mm) <- rownames(mm)

#calculate aldex CLR 
clr <- aldex.clr(d.1, mm, mc.samples = 128)

#Aldex GLM functions
glm_p <- aldex.glm(clr, mm)
glm_effect <- aldex.glm.effect(clr, useMC=TRUE)
glm_output <- data.frame(glm_p,glm_effect)
write.table(glm_output, file="glm_tax_TimeDonor.txt", sep="\t", quote=F, col.names=NA)

################################################
#EC enzymes
################################################
d<-read.table("picrust2/output/EC_metagenome_out/pred_metagenome_unstrat.tsv", sep="\t", quote="", check.names=F, header=T, row.names=1, comment.char="")

#picrust output countains decimals, round to 0 digits
d.1<-round(d, digits=0)
write.table(d.1, file="picrust2/output/EC_metagenome_out/pred_metagenome_unstrat_rounded.txt", sep="\t", quote=F, col.names=NA)

#remove unnecessary samples from #008 and #013 
#read in the new table
d<-read.table("picrust2/output/EC_metagenome_out/pred_metagenome_unstrat_rounded.txt", sep="\t", quote="", check.names=F, header=T, row.names=1, comment.char="")
d.1 <- data.frame(d)

# Can't have zeroes. Bayesian-Multiplicative replacement of count zeros
d.czm <- cmultRepl(t(d.1),  label=0, method="CZM")

#Centre-log-ratio transform the data
d.clr<-apply(d.czm, 2, function(x){log2(x) - mean(log2(x))})

#calculate principal components
#features are COLUMNS
#samples as ROWS
d.pcx <- prcomp(d.clr)

fit1 <- envfit(d.pcx, data.frame(m[, c("Patient","FMT","FMT_Donor","EarlyLate","Time"), drop = FALSE]), permutations = 999, na.rm = TRUE)
fit1

// ***VECTORS
// 
//                PC1      PC2     r2 Pr(>r)    
// Patient   -0.77965 -0.62621 0.2971  0.001 ***
// FMT_Donor -0.34653 -0.93804 0.1455  0.001 ***
// Time       0.92526  0.37933 0.0165  0.509    
// ---
// Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
// Permutation: free
// Number of permutations: 999
// 
// ***FACTORS:
// 
// Centroids:
//                PC1     PC2
// FMTFMT      3.1308 -0.5531
// FMTNo       1.1305  0.0258
// EarlyLateE  0.8372  5.2178
// EarlyLateL  2.9281 -3.9302
// 
// Goodness of fit:
//               r2 Pr(>r)
// FMT       0.0009  0.913
// EarlyLate 0.0172  0.225
// Permutation: free
// Number of permutations: 999
// 
// 10 observations deleted due to missingness

################################################
#KO KEGG orthologs
################################################
d<-read.table("picrust2/output/KO_metagenome_out/pred_metagenome_unstrat.tsv", sep="\t", quote="", check.names=F, header=T, row.names=1, comment.char="")

#picrust output countains decimals, round to 0 digits
d.1<-round(d, digits=0)
write.table(d.1, file="picrust2/output/KO_metagenome_out/pred_metagenome_unstrat_rounded.txt", sep="\t", quote=F, col.names=NA)

#remove unnecessary samples from #008 and #013 
#read in the new table
d<-read.table("picrust2/output/KO_metagenome_out//pred_metagenome_unstrat_rounded.txt", sep="\t", quote="", check.names=F, header=T, row.names=1, comment.char="")
d.1 <- data.frame(d)

# Can't have zeroes. Bayesian-Multiplicative replacement of count zeros
d.czm <- cmultRepl(t(d.1),  label=0, method="CZM")

#Centre-log-ratio transform the data
d.clr<-apply(d.czm, 2, function(x){log2(x) - mean(log2(x))})

#calculate principal components
#features are COLUMNS
#samples as ROWS
d.pcx <- prcomp(d.clr)

fit1 <- envfit(d.pcx, data.frame(m[, c("Patient","FMT","FMT_Donor","EarlyLate","Time"), drop = FALSE]), permutations = 999, na.rm = TRUE)
fit1
// 
// ***VECTORS
// 
//                PC1      PC2     r2 Pr(>r)    
// Patient   -0.15973  0.98716 0.3419  0.001 ***
// FMT_Donor  0.22222  0.97500 0.2115  0.001 ***
// Time       0.09436 -0.99554 0.0189  0.429    
// ---
// Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
// Permutation: free
// Number of permutations: 999
// 
// ***FACTORS:
// 
// Centroids:
//                PC1     PC2
// FMTFMT      2.2793 -6.3074
// FMTNo       3.4653  1.4293
// EarlyLateE  0.5079 -4.5975
// EarlyLateL  4.5060 -0.7057
// 
// Goodness of fit:
//               r2 Pr(>r)
// FMT       0.0031  0.749
// EarlyLate 0.0015  0.863
// Permutation: free
// Number of permutations: 999
// 
// 10 observations deleted due to missingness

################################################
#Picrust2 Pathways
################################################
d<-read.table("picrust2/output/pathways_out/path_abun_unstrat.tsv", sep="\t", quote="", check.names=F, header=T, row.names=1, comment.char="")

#picrust output countains decimals, round to 0 digits
d.1<-round(d, digits=0)
write.table(d.1, file="picrust2/output/pathways_out/path_abun_unstrat_rounded.txt", sep="\t", quote=F, col.names=NA)

#remove unnecessary samples from #008 and #013 
#read in the new table
d<-read.table("picrust2/output/pathways_out/path_abun_unstrat_rounded.txt", sep="\t", quote="", check.names=F, header=T, row.names=1, comment.char="")
d.1 <- data.frame(d)

# Can't have zeroes. Bayesian-Multiplicative replacement of count zeros
d.czm <- cmultRepl(t(d.1),  label=0, method="CZM")

#Centre-log-ratio transform the data
d.clr<-apply(d.czm, 2, function(x){log2(x) - mean(log2(x))})

#calculate principal components
#features are COLUMNS
#samples as ROWS
d.pcx <- prcomp(d.clr)

fit1 <- envfit(d.pcx, data.frame(m[, c("Patient","FMT","FMT_Donor","EarlyLate","Time"), drop = FALSE]), permutations = 999, na.rm = TRUE)
fit1
// 
// ***VECTORS
// 
//                PC1      PC2     r2 Pr(>r)    
// Patient   -0.37643 -0.92644 0.2672  0.001 ***
// FMT_Donor -0.13229 -0.99121 0.2094  0.001 ***
// Time       0.59378  0.80463 0.0273  0.295    
// ---
// Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
// Permutation: free
// Number of permutations: 999
// 
// ***FACTORS:
// 
// Centroids:
//                PC1     PC2
// FMTFMT      1.8831  2.0869
// FMTNo      -0.5039 -1.6573
// EarlyLateE -1.0288  1.3888
// EarlyLateL  1.7588 -0.7111
// 
// Goodness of fit:
//               r2 Pr(>r)
// FMT       0.0215  0.143
// EarlyLate 0.0128  0.282
// Permutation: free
// Number of permutations: 999
// 
// 10 observations deleted due to missingness
