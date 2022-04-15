#Import libraries
library(dplyr)
library(tidyr)
library(data.table)
library(tibble)
library(ggplot2)
library(stringr)
#Set working directory
setwd("/home/ionadmin/myeloseq_report")
QC_file.list = list.files(path = getwd(), "*_filtered_variants.tsv", full.names = T)
run_id <- read.csv(QC_file.list[1],sep = "\t")[1,1]

process_SC_control <- function(x) {
  ######## process sensitive control data ######################
  myelodata<-read.csv(x,sep = "\t")
  if (grepl( "m2", myelodata[1,2], fixed = TRUE)) {
    chip = 2
  } else { chip = 1 }
  myelodata$locus<-paste0(myelodata$Genes,":", myelodata$Coding)
  sort(myelodata$locus)
  DNAqcdata<-subset(myelodata,  locus  %in% c (
  "ABL1:c.944C>T",
  "BRAF:c.1799T>A",
  "CALR:c.1099_1150delCTTAAGGAGGAGGAAGAAGACAAGAAACGCAAAGAGGAGGAGGAGGCAGAGG",
  "CBL:c.1139T>C",
  "CBL:c.1259G>A",
  "CEBPA:c.937_939dup",
  "CSF3R:c.1853C>T",
  "FLT3:c.1759_1800dup",
  "FLT3:c.2503G>T",
  "FLT3:c.1806_1807insGGGGCTTTCAGAGAATATGAATATGATCTCAAA",
  "IDH1:c.394C>T",
  "JAK2:c.1624_1629delAATGAA",
  "JAK2:c.1849G>T",
  "MPL:c.1544G>T",
  "NPM1:c.863_864insTCTG",
  "SF3B1:c.2098A>G",
  "SF3B1:c.1998G>T",
  "MIR636,SRSF2:c.284_307delCCCCGGACTCACACCACAGCCGCC",
  "U2AF1,U2AF1L5:c.-8623C>T, c.101C>T")  )
  DNAqcdata$locus<- gsub("CALR:c.1099_1150delCTTAAGGAGGAGGAAGAAGACAAGAAACGCAAAGAGGAGGAGGAGGCAGAGG","CALR:c.1099_1150del", DNAqcdata$locus)
  DNAqcdata$locus<- gsub("FLT3:c.1806_1807insGGGGCTTTCAGAGAATATGAATATGATCTCAAA","FLT3:c.1806_1807ins33", DNAqcdata$locus)
  DNAqcdata$locus<- gsub("FLT3:c.1799_1800insTAATGAGTACTTCTACGTTGATTTCAGAGAATATGAATATGA","FLT3:c.1799_1800ins42", DNAqcdata$locus)
  DNAqcdata$locus<- gsub("JAK2:c.1624_1629delAATGAA","JAK2:c.1624_1629del", DNAqcdata$locus)
  DNAqcdata$locus<- gsub("SRSF2\\|MIR636\\|MFSD11:c.284_307delCCCCGGACTCACACCACAGCCGCC.*1500CGGCGGCTGTGGTGTGAGTCCGGGG>C","SRSF2:c.284_307del", DNAqcdata$locus)
  DNAqcdata$locus<- gsub("U2AF1\\|U2AF1L5:c.101C>T\\|c.-8623C>T","U2AF1:c.101C>T", DNAqcdata$locus)
  run_id <- DNAqcdata[1,1]
  DNAqcdata$AF <- as.numeric(as.character(DNAqcdata$X..Frequency))
  DNAQC <- ggplot(DNAqcdata, aes(x= locus, y=AF,color= locus)) +
    geom_point(size=6)+
    ggtitle(paste("Myeloid SC-Seraseq_",run_id,"DNA QC", paste("(chip#",chip,")",sep=""),":19 Sites Expected",sep = " ")) +
    xlab("Genomic Location") + ylab("Allele Frequency") +
    geom_hline(yintercept=3, linetype="dashed",  color = "red", size=1)+
    geom_hline(yintercept=15, linetype="dashed",  color = "blue", size=1)+
    theme(axis.text.x = element_text(size = 7, angle = 45), legend.title = element_text(size=8),legend.text = element_text(size = 7), legend.key.size = unit(0.3, "cm")  )+
    scale_color_manual(values=c("goldenrod3", "brown2", "deepskyblue1", "plum1", "turquoise2", "orchid", "greenyellow","seagreen3","hotpink", "lightskyblue2", "mediumpurple1", "firebrick1","olivedrab", "blue",  "orange1", "darkolivegreen3", "tomato3", "tan4", "red4", "purple1", "gray4", "grey66","goldenrod2" ))+
    guides(col = guide_legend(ncol = 1)) + theme(plot.title = element_text(hjust = 0.5))
  fusions_df <- subset(myelodata, Type == "FUSION" | Type == "RNAExonVariant")
  #fusions_df$Fusions <- paste(fusions_df$Genes,"(",fusions_df$Exon,")",sep="")
  RNAqcdata<-subset(fusions_df, Genes  %in% c (
    "KAT6A(17) - CREBBP(2)",
    "ETV6(4) - ABL1(2)",
    "ETV6(5) - ABL1(2)",
    "PCM1(23) - JAK2(12)",
    "FIP1L1(11) - PDGFRA(12)",
    "BCR(14) - ABL1(2)",
    "RUNX1(3) - RUNX1T1(3)",
    "PML(6) - RARA(3)"  )   )
  RNAqcdata$Read.Counts[is.na(RNAqcdata$Read.Counts)] <- 1
  RNAQC <- ggplot(RNAqcdata, aes(x=Genes, y=log2(as.numeric(as.character(Read.Counts))),color= Genes)) +
    geom_point(size=6)+ xlab("Fusions(Exome)") + ylab("log2 Read Counts") +
    ggtitle(paste("Myeloid SC_", run_id, "RNA QC",paste("(chip#",chip,")",sep=""),": 8 Events Expected",sep = " ")) +
    geom_hline(yintercept=log2(100), linetype="dashed",  color = "red", size=1)+
    theme(axis.text.x = element_text(size = 7, angle = 45), legend.title = element_text(size=8),legend.text =  element_text(size = 7), legend.key.size = unit(0.3, "cm")  )+
    scale_color_manual(values=c("seagreen3","hotpink", "lightskyblue2", "mediumpurple1", "firebrick1","olivedrab", "blue",  "orange1", "darkolivegreen3", "tomato3", "tan4", "red4", "purple1", "gray4", "grey66","goldenrod2", "brown", "deepskyblue3", "plum3", "turquoise", "orchid2", "greenyellow","salmon2" ))+
    guides(col = guide_legend(ncol = 1)) + theme(plot.title = element_text(hjust = 0.5))
  list(DNAQC, RNAQC)
}

QCplot = lapply(QC_file.list, process_SC_control) #list(DNAQC, RNAQC)
pdf(paste("Myeloid-SC-QC-",run_id,".pdf",sep=""))
QCplot
dev.off()
####### process SNP QC ##############
#read unfilted TSV files
files_to_read <- list.files(path=paste(getwd(),"Variants",sep="/"), pattern = "_Non-Filtered.*-oncomine\\.tsv$", recursive = TRUE)
#Read all tsv files in directory
all_files <- lapply(files_to_read,function(x) {
  read.csv(file = paste(getwd(),"Variants",x,sep="/"), quote="", sep = '\t', header = TRUE, skip=2)
})
#Combine file content list and file name list
all_lists <- mapply(c, files_to_read,all_files, SIMPLIFY = FALSE)
#Unlist all lists and add file name as new column
all_result <- rbindlist(all_lists, fill = T)
#Rename new column
names(all_result)[1] <- "Sample"
#Shorten Sample name by removing all characters after the first underscore
all_result$Sample<-gsub("_.*","",all_result$Sample)
#List of locus
snps <- read.csv("/home/ionadmin/myeloseq_report/QC_data/Myeloid.SNP_locus.csv")
snps <- snps$Locus
#add the 2 CNV number replacing the VAF
combined_final <- distinct(all_result)
combined_final <- as.data.frame(sapply(all_result, function(x) gsub("\"", "", x)))
data0<-combined_final[,c("Sample"  , "X.CHROM." ,"X.POS." , "X.INFO.A.AF.")]
data0$X..Frequency<-data0$X.INFO.A.AF.
data0$Locus <- paste(data0$X.CHROM.,data0$X.POS.,sep = ":")
data0$X..Frequency<-gsub("NA","",data0$X..Frequency)
data0 <- data0[,c("Sample","Locus","X..Frequency")]
data0 <- distinct(data0)
res.targeloci<-subset(data0, Locus  %in% snps)
data.e <- pivot_wider(
  res.targeloci,
  names_from = "Sample",
  values_from = "X..Frequency",
  values_fn = list(X..Frequency= list)) %>%
  unchop(everything())
#data.e<-spread(res.targeloci, Sample, X..Frequency)
data.e[is.na(data.e)] <- "0"
write.csv(data.e, file=paste("QC-SNPs-final-",run_id,".csv",sep=""), row.names=FALSE)
write.csv(combined_final, file=paste("combined_final_",run_id,".csv",sep=""), row.names=FALSE)
#imbalancefusion<-subset(combined_final,  (grepl(glob2rx("5p3pAssays") , X.rowtype.) )
#& (as.numeric(X.INFO.1.5P_3P_ASSAYS.) > 0)
#& !(X.call. == "NEG") & !(X.call. == "NOCALL"))
#write.csv(imbalancefusion, file=paste("imbalancefusion_final_",run_id,".csv",sep=""), row.names=FALSE)
# expcontrol QC
expcontrol<-subset(combined_final,  (grepl(glob2rx("ExprControl") , X.rowtype.) ))
write.csv(expcontrol, file=paste("expcontrol_",run_id,".csv",sep=""), row.names=FALSE)
#Read.Counts
expQC<- ggplot(expcontrol, aes(y=log2(as.numeric(as.character(X.INFO...READ_COUNT.))), x=Sample, color=Sample)) +
  geom_point(size=3)+ ylab("log2 Read Counts") +
  theme(axis.text.x = element_text(size = 6, angle = 65))+
  facet_wrap( ~X.INFO.1.GENE_NAME., ncol =2)+
  geom_hline(yintercept=log2(1800), linetype="dashed",  color = "red", size=0.5)+
  geom_hline(yintercept=log2(700), linetype="dashed",  color = "blue", size=0.5)+
  ggtitle(paste("Case Expression Control Myeloid QC_",run_id,sep="")) +
  guides(col = guide_legend ( ncol=1 )) +
  theme(legend.text=element_text(size=7))
pdf(paste("Case_Expression_Control_Myeloid_QC_",run_id,".pdf",sep=""))
expQC
dev.off()
## process filtered full TSV file #######
tsv_files_to_read <- list.files(path=paste(getwd(),"Variants",sep="/"), pattern = "-full.tsv$", recursive = TRUE)
#Read all tsv files in directory
tsv_files <- lapply(tsv_files_to_read,function(x) {
  read.csv(file = paste(getwd(),"Variants",x,sep="/"), quote="", sep = '\t', header = TRUE, skip=2)
})
#Combine file content list and file name list
content_lists <- mapply(c, tsv_files_to_read,tsv_files, SIMPLIFY = FALSE)
#Unlist all lists and add file name as new column
filtered_result <- rbindlist(content_lists, fill = T)
#Rename new column
names(filtered_result)[1] <- "Sample"
#Shorten Sample name by removing all characters after the first underscore
filtered_result$Sample<-gsub("_.*","",filtered_result$Sample)
#List of locus
#add the 2 CNV number replacing the VAF
combined_tsv <- distinct(filtered_result)
combined_tsv <- as.data.frame(sapply(filtered_result, function(x) gsub("\"", "", x)))
biallelicdata<-subset(combined_tsv, (grepl("|", protein, fixed = TRUE))
                      & (filter == "PASS") & (type %in% c("MNV","SNV")))
write.csv(biallelicdata, file=paste("biallelicdata-final-",run_id,".csv",sep=""), row.names=FALSE)