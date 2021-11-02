# R
# Maaslin2 determines associations
# use significant covariates from envfit results

################################################
#Taxonomy
################################################
library(zCompositions)
library(Maaslin2)

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
mp <- read.table("metadata_patient.txt", header=T, sep='\t', comment.char = "", row.names=1)
m1 <- read.table("metadata_donor1.txt", header=T, sep='\t', comment.char = "", row.names=1)
m2 <- read.table("metadata_donor2.txt", header=T, sep='\t', comment.char = "", row.names=1)

m<-data.frame(m)
mp<-data.frame(mp)
m1<-data.frame(m1)
m2<-data.frame(m2)

# make the name of the new folder it will with all your data in it

fit_data <- Maaslin2(d.czm, m, 'maaslin/tax_all_samples', transform = "none", normalization = "CLR", analysis_method = "LM", 
			fixed_effects = c("Time","FMT_Donor","EarlyLate"),
			random_effects = c("Patient"),
			standardize = FALSE)

#repeat with subgroups and change fixed effects to significant features from envfit

fit_data <- Maaslin2(d.czm, mp, 'maaslin/tax_patients_only', transform = "none", normalization = "CLR", analysis_method = "LM", 
			fixed_effects = c("Time","FMT_Donor","EarlyLate"),
			random_effects = c("Patient"),
			standardize = FALSE)
fit_data <- Maaslin2(d.czm, m1, 'maaslin/tax_donor1_recipients', transform = "none", normalization = "CLR", analysis_method = "LM", 
			fixed_effects = c("Time"),
			random_effects = c("Patient"),
			standardize = FALSE)
fit_data <- Maaslin2(d.czm, m2, 'maaslin/tax_donor2_recipients', transform = "none", normalization = "CLR", analysis_method = "LM", 
			fixed_effects = c("Time"),
			random_effects = c("Patient"),
			standardize = FALSE)

################################################
#EC Enzymes
################################################

d<-read.table("picrust2/output/EC_metagenome_out/pred_metagenome_unstrat_rounded.txt", sep="\t", quote="", check.names=F, header=T, row.names=1, comment.char="")
d.1 <- data.frame(d)

# Can't have zeroes. Bayesian-Multiplicative replacement of count zeros
d.czm <- cmultRepl(t(d.1),  label=0, method="CZM")

# make the name of the new folder it will with all your data in it

fit_data <- Maaslin2(d.czm, m, 'maaslin/EC_all_samples', transform = "none", normalization = "CLR", analysis_method = "LM", 
			fixed_effects = c("FMT_Donor","Time"),
			random_effects = c("Patient"),
			standardize = FALSE)
fit_data <- Maaslin2(d.czm, mp, 'maaslin/EC_patients_only', transform = "none", normalization = "CLR", analysis_method = "LM", 
			fixed_effects = c("FMT_Donor"),
			random_effects = c("Patient"),
			standardize = FALSE)
fit_data <- Maaslin2(d.czm, m1, 'maaslin/EC_donor1_recipients', transform = "none", normalization = "CLR", analysis_method = "LM", 
			fixed_effects = c("Time"),
			random_effects = c("Patient"),
			standardize = FALSE)
fit_data <- Maaslin2(d.czm, m2, 'maaslin/EC_donor2_recipients', transform = "none", normalization = "CLR", analysis_method = "LM", 
			fixed_effects = c("Time"),
			random_effects = c("Patient"),
			standardize = FALSE)

################################################
#KO KEGG orthologs
################################################

d<-read.table("picrust2/output/KO_metagenome_out//pred_metagenome_unstrat_rounded.txt", sep="\t", quote="", check.names=F, header=T, row.names=1, comment.char="")
d.1 <- data.frame(d)

# Can't have zeroes. Bayesian-Multiplicative replacement of count zeros
d.czm <- cmultRepl(t(d.1),  label=0, method="CZM")

# make the name of the new folder it will with all your data in it

fit_data <- Maaslin2(d.czm, m, 'maaslin/KO_all_samples', transform = "none", normalization = "CLR", analysis_method = "LM", 
			fixed_effects = c("FMT_Donor"),
			random_effects = c("Patient"),
			standardize = FALSE)
fit_data <- Maaslin2(d.czm, mp, 'maaslin/KO_patients_only', transform = "none", normalization = "CLR", analysis_method = "LM", 
			fixed_effects = c("FMT_Donor"),
			random_effects = c("Patient"),
			standardize = FALSE)
fit_data <- Maaslin2(d.czm, m1, 'maaslin/KO_donor1_recipients', transform = "none", normalization = "CLR", analysis_method = "LM", 
			fixed_effects = c("Time"),
			random_effects = c("Patient"),
			standardize = FALSE)
fit_data <- Maaslin2(d.czm, m2, 'maaslin/KO_donor2_recipients', transform = "none", normalization = "CLR", analysis_method = "LM", 
			fixed_effects = c("Time"),
			random_effects = c("Patient"),
			standardize = FALSE)

################################################
#Picrust2 Pathways
################################################

d<-read.table("picrust2/output/pathways_out/path_abun_unstrat_rounded.txt", sep="\t", quote="", check.names=F, header=T, row.names=1, comment.char="")
d.1 <- data.frame(d)

# Can't have zeroes. Bayesian-Multiplicative replacement of count zeros
d.czm <- cmultRepl(t(d.1),  label=0, method="CZM")

# make the name of the new folder it will with all your data in it

fit_data <- Maaslin2(d.czm, m, 'maaslin/paths_all_samples', transform = "none", normalization = "CLR", analysis_method = "LM", 
			fixed_effects = c("FMT_Donor"),
			random_effects = c("Patient"),
			standardize = FALSE)
fit_data <- Maaslin2(d.czm, mp, 'maaslin/paths_patients_only', transform = "none", normalization = "CLR", analysis_method = "LM", 
			fixed_effects = c("FMT_Donor"),
			random_effects = c("Patient"),
			standardize = FALSE)
fit_data <- Maaslin2(d.czm, m1, 'maaslin/paths_donor1_recipients', transform = "none", normalization = "CLR", analysis_method = "LM", 
			fixed_effects = c("Time"),
			random_effects = c("Patient"),
			standardize = FALSE)
fit_data <- Maaslin2(d.czm, m2, 'maaslin/paths_donor2_recipients', transform = "none", normalization = "CLR", analysis_method = "LM", 
			fixed_effects = c("Time"),
			random_effects = c("Patient"),
			standardize = FALSE)
