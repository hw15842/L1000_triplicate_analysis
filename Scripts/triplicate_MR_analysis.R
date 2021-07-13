
##############################
## triplicate_MR_analysis.R ##
##############################



args  <-  commandArgs(trailingOnly=TRUE)
results_location <- toString(args[1])
data_location <- toString(args[2])
prot_or_expression <- toString(args[3])  ## whether we are reading in the protein pairs or the expression pairs 


library("data.table")
library("dplyr")
library("plyr")
library("TwoSampleMR")


PE_pairs <- get(load(paste0(data_location, prot_or_expression, "pairs_for_triplicate_analysis.rdata"))

## add in a column that is the IEUgwas IDs for the proteins or exposures in the correct orientation

PE_pairs$ID_PE1 <- ifelse(PE_pairs$steiger_dir == "TRUE", PE_pairs$id.exposure, PE_pairs$id.outcome)

PE_pairs$ID_PE2 <- ifelse(PE_pairs$steiger_dir == "TRUE", PE_pairs$id.outcome, PE_pairs$id.exposure)



## PE will be used to refer to either the protein or the expression throughout this script 
# i.e. PE1 is either protein 1 or expression 1 depending on which data has been read in 

### We need:

# betaPE = Verified effect estimate between PE1 and PE2 - is this the beta value already in the dataset or should we re do it with all instruments? 
# betaT = measured effect of PE2 on T 
# betaF = measured effect of PE1 on T 



## beta of PE1 on PE2 = b_PE1_PE2        (verified direction from L1000)
## beta of PE1 on trait = b_PE1_T
## beta of PE2 on trait = b_PE2_t




## read in the trait list 

trait_list <- get(load(paste0(data_location, "......................")))




# perform the MR of each PE on all the traits
# we will use all the instruments for both PE1 and PE2 not just the pQTL/eQTL that means they are associated 


MR_func <- function(trait){

	MR_PE1_on_trait <- mr(make_dat(PE_pairs$ID_PE1, trait), method_list = "mr_ivw")
	MR_PE1_on_trait <- MR_PE1_on_trait[,c(1,2,6:9)]
	colnames(MR_PE1_on_trait) <- paste(colnames(MR_PE1_on_trait), "PE1_T", sep = "_")

	MR_PE2_on_trait <- mr(make_dat(PE_pairs$ID_PE2, trait), method_list = "mr_ivw")
	MR_PE2_on_trait <- MR_PE2_on_trait[,c(1,2,6:9)]
	colnames(MR_PE2_on_trait) <- paste(colnames(MR_PE2_on_trait), "PE2_T", sep = "_")

	
	res <- cbind(MR_PE1_on_trait, MR_PE2_on_trait)
	res <- res[,c(1,7,2:6,9:12)]
	colnames(res)[1:3] <- c("id_PE1", "id_PE2", "id_trait")

	MR_PE1_on_PE2 <- mr(make_dat(res$id_PE1, res$id_PE2), method_list = "mr_ivw")
	MR_PE1_on_PE2 <- MR_PE1_on_PE2[,c(1,2,6:9)]
	MR_PE1_on_PE2 <- subset(MR_PE1_on_PE2, MR_PE1_on_PE2$id.exposure %in% res$id_PE1 & MR_PE1_on_PE2$id.outcome %in% res$id_PE2)

	
}





















