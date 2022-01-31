#Import libraries
library(dplyr)
library(tidyr)
library(data.table)
library(tibble)
library(ggplot2)

#Set working directory
setwd("/home/ionadmin/ion_report")

######## process sensitive control data ######################
oncominedata<-read.csv("sc_filtered_variants.tsv",sep = "\t")
oncominedata$locus<-paste0(oncominedata$Genes,":", oncominedata$Coding)

DNAqcdata<-subset(oncominedata,  locus  %in% c (
  "NRAS:c.182A>G",
  "IDH1:c.394C>T",
  "PIK3CA:c.1633G>A",
  "KIT:c.2447A>T",
  "EGFR:c.2236_2250delGAATTAAGAGAAGCA",
  "EGFR:c.2573T>G",
  "EGFR:c.2369C>T",
  "BRAF:c.1799T>A",
  "GNAQ:c.626A>C",
  "KRAS:c.35G>A",
  "ERBB2:c.2313_2324dup",
  "GNA11:c.626A>T" )  )

run_id <- DNAqcdata[1,1]
DNAqcdata$AF <- as.numeric(as.character(DNAqcdata$X..Frequency))

DNAQC <- ggplot(DNAqcdata, aes(x= locus, y=AF,color= locus)) +
  geom_point(size=6)+
  ggtitle(paste("Oncomine Solid",run_id,"DNA QC: 12 Sites Expected",sep = " ")) +
  xlab("Genomic Location") + ylab("Allele Frequency") +
  geom_hline(yintercept=3, linetype="dashed",  color = "red", size=1)+
  geom_hline(yintercept=15, linetype="dashed",  color = "blue", size=1)+
  theme(axis.text.x = element_text(size = 7, angle = 45), legend.title = element_text(size=8),legend.text = element_text(size = 7), legend.key.size = unit(0.3, "cm")  )+
  scale_color_manual(values=c("goldenrod3", "brown2", "deepskyblue1", "plum1", "turquoise2", "orchid", "greenyellow","seagreen3","hotpink", "lightskyblue2", "mediumpurple1", "firebrick1","olivedrab", "blue",  "orange1", "darkolivegreen3", "tomato3", "tan4", "red4", "purple1", "gray4", "grey66","goldenrod2" ))+
  guides(col = guide_legend(ncol = 1)) + theme(plot.title = element_text(hjust = 0.5))

fusions_df <- subset(oncominedata, Type == "FUSION" | Type == "RNAExonVariant")
fusions_df$Fusions <- paste(fusions_df$Genes,"(",fusions_df$Exon,")",sep="")
RNAqcdata<-subset(fusions_df, Fusions  %in% c (
  "EML4|ALK(13|20)",
  "KIF5B|RET(24|11)",
  "NCOA4|RET(7|12)",
  "CD74|ROS1(6|34)",
  "SLC34A2|ROS1(4|34)",
  "TPM3|NTRK1(7|10)",
  "FGFR3|BAIAP2L1(17|2)",
  "PAX8|PPARG(9|2)",
  "FGFR3|TACC3(17|11)",
  "ETV6|NTRK3(5|15)",
  "LMNA|NTRK1(2|11)",
  "SLC45A3|BRAF(1|8)",
  "TMPRSS2|ERG(1|2)",
  "MET|MET(13|15)"  )  )

RNAqcdata$Read.Counts[is.na(RNAqcdata$Read.Counts)] <- 1

RNAQC <- ggplot(RNAqcdata, aes(x=Fusions, y=log2(as.numeric(as.character(Read.Counts))),color= Fusions)) +
  geom_point(size=6)+ xlab("Fusions(Exome)") + ylab("log2 Read Counts") +
  ggtitle(paste("Oncomine Solid SC-QC", run_id, "RNA QC: 14 Events Expected",sep = " ")) +
  geom_hline(yintercept=log2(100), linetype="dashed",  color = "red", size=1)+
  theme(axis.text.x = element_text(size = 7, angle = 45), legend.title = element_text(size=8),legend.text =  element_text(size = 7), legend.key.size = unit(0.3, "cm")  )+
  scale_color_manual(values=c("seagreen3","hotpink", "lightskyblue2", "mediumpurple1", "firebrick1","olivedrab", "blue",  "orange1", "darkolivegreen3", "tomato3", "tan4", "red4", "purple1", "gray4", "grey66","goldenrod2", "brown", "deepskyblue3", "plum3", "turquoise", "orchid2", "greenyellow","salmon2" ))+
  guides(col = guide_legend(ncol = 1)) + theme(plot.title = element_text(hjust = 0.5))

QCplot = list(DNAQC, RNAQC)
pdf(paste("Oncomine-Solid-SC-QC-",run_id,".pdf",sep=""))
QCplot
dev.off()

####### process SNP QC ##############
#read unfilted TSV files
files_to_read <- list.files(path="/home/ionadmin/ion_report/Variants", pattern = "_Non-Filtered.*-oncomine\\.tsv$", recursive = TRUE)

#Read all tsv files in directory
all_files <- lapply(files_to_read,function(x) {
  read.csv(file = paste("/home/ionadmin/ion_report/Variants",x,sep="/"), quote="", sep = '\t', header = TRUE, skip=4)
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
snps <- read.csv("/home/ionadmin/ion_report/data/ONCsolid.SNP_locus.csv")
snps <- snps$Locus

#add the 2 CNV number replacing the VAF
combined_final <- distinct(all_result)
combined_final <- as.data.frame(sapply(all_result, function(x) gsub("\"", "", x)))

data0<-combined_final[,c("Sample"  , "X.CHROM." ,"X.POS." , "X.INFO.A.AF.", "X.INFO...CI.")]
data0$X..Frequency<-paste0(data0$X.INFO.A.AF., data0$X.INFO...CI.)
data0$Locus <- paste(data0$X.CHROM.,data0$X.POS.,sep = ":")
data0$X..Frequency<-gsub("NA","",data0$X..Frequency)
data0 <- data0[,c("Sample","Locus","X..Frequency")]
data0 = data0[!duplicated(data0$Locus),]
res.targeloci<-subset(distinct(data0), Locus  %in% snps)
data.e<-spread(res.targeloci, Sample, X..Frequency)
data.e[is.na(data.e)] <- 0

write.csv(data.e, file=paste("QC-SNPs-final-",run_id,".csv",sep=""), row.names=FALSE)
write.csv(combined_final, file=paste("combined_final_",run_id,".csv",sep=""), row.names=FALSE)

imbalancefusion<-subset(combined_final,  (grepl(glob2rx("5p3pAssays") , X.rowtype.) )
                        & (as.numeric(X.INFO.1.5P_3P_ASSAYS.) > 0)
                        & !(X.call. == "NEG") & !(X.call. == "NOCALL"))
write.csv(imbalancefusion, file=paste("imbalancefusion_final_",run_id,".csv",sep=""), row.names=FALSE)

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
  ggtitle("Case Expression Control ONC.Solid QC") +
  guides(col = guide_legend ( ncol=1 )) +
  theme(legend.text=element_text(size=7))

pdf(paste("Case_Expression_Control_ONC_Solid_QC_",run_id,".pdf",sep=""))
expQC
dev.off()

## process filtered full TSV file #######
tsv_files_to_read <- list.files(path="/home/ionadmin/ion_report/Variants", pattern = "-full.tsv$", recursive = TRUE)

#Read all tsv files in directory
tsv_files <- lapply(tsv_files_to_read,function(x) {
  read.csv(file = paste("/home/ionadmin/ion_report/Variants",x,sep="/"), quote="", sep = '\t', header = TRUE, skip=2)
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
snps <- read.csv("/home/ionadmin/ion_report/data/ONCsolid.SNP_locus.csv")
snps <- snps$Locus

#add the 2 CNV number replacing the VAF
combined_tsv <- distinct(filtered_result)
combined_tsv <- as.data.frame(sapply(filtered_result, function(x) gsub("\"", "", x)))
biallelicdata<-subset(combined_tsv, (grepl("?|p", protein, fixed = TRUE))
                      & (filter == "PASS") & (type %in% c("MNV","SNV","INDEL,MNV","INDEL","INDEL,SNV")))
write.csv(biallelicdata, file=paste("biallelicdata-final-",run_id,".csv",sep=""), row.names=FALSE)