
#Data cleaning of eQTLGen####
setwd()   # Setting working firectionary
library(data.table)
library(tidyverse)      #Packages loading and library
library(dplyr)
library(readxl)
eQTLgen <- fread(file = "2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt",
                 header = T)

druggable_gene <- read_xlsx("Druggable gene.xlsx")  #Drugable genes 
colnames(druggable_gene) <- "GeneSymbol"
head(eQTLgen)
eQTLgen <- inner_join(eQTLgen, druggable_gene, by="GeneSymbol")   #Intersection 
length(unique(eQTLgen$GeneSymbol))

max(eQTLgen$Pvalue)
eQTLgen_sig <- subset(eQTLgen, eQTLgen$Pvalue<5e-08)     #P-value < 5e-08

Freq<-fread(input = "2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt", header = T )
Freq <- Freq[ ,c("SNP","AlleleB_all")]
eQTLgen_sig<-merge(eQTLgen_sig, Freq, by="SNP")
head(eQTLgen_sig)
eQTLgen_sig$Beta <- eQTLgen_sig$Zscore/sqrt(2*eQTLgen_sig$AlleleB_all*(1-eQTLgen_sig$AlleleB_all)*(eQTLgen_sig$NrSamples+eQTLgen_sig$Zscore^2))
eQTLgen_sig$SE <- abs(eQTLgen_sig$Beta/eQTLgen_sig$Zscore)   #Gaining beta and SE of each gene

library(TwoSampleMR)
library(ieugwasr)
Sys.setenv(OPENGWAS_JWT= "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJzb25namlhaGFvMjU2N0AxNjMuY29tIiwiaWF0IjoxNzQyMDU0NTY2LCJleHAiOjE3NDMyNjQxNjZ9.phLwIW-htGqSWwtvgQU6PqQMZkVpkP3Ody6xDp2EIW28yooRrTMQVDQ-gy8B1zaOrfg7PdI9-9Y-LJrvTBMxG0z7Xe8Hho6p1c7QTd2rNb_Q1Ricl8LAxJQZjjYCW6LuYL-llBtDKRmbkePYSHcVU7sOmPOki__5XLp7a9a8T9lLtL_0CL2wml9UC8_9ayvQgzXYNAYuPzfTTMuVqfM2XnZDuyahlqVbq4QlrA9xcmJribTkxVaLFo_2XUPc96u6TTR9Qf_TwRzYUh6esJa0POAuPRvnCv8Vnsv06FDFnjNYCGkiSTYCV-cgTWISiqTGOMM0hU74o43V9KKwSVNYkA")
eQTLgen_exp <- read_exposure_data(filename="eQTLgen_sig.csv",  
                                  clump = FALSE,
                                  sep = ",",
                                  phenotype_col = "GeneSymbol",
                                  snp_col = "SNP",
                                  beta_col = "Beta",
                                  se_col = "SE",
                                  eaf_col = "AlleleB_all",
                                  effect_allele_col = "AssessedAllele",
                                  other_allele_col = "OtherAllele",
                                  pval_col = "Pvalue",
                                  units_col = "units",
                                  ncase_col = "ncase",
                                  ncontrol_col = "ncontrol",
                                  samplesize_col = "NrSamples",
                                  gene_col = "GeneSymbol",
                                  id_col = "GeneSymbol",
                                  min_pval = 1e-200,
                                  log_pval = FALSE,
                                  chr_col = "SNPChr",
                                  pos_col = "SNPPos")
eQTLgen_clump <- clump_data(dat=eQTLgen_exp,             #Clumping to remove those with LD
                            clump_kb = 1000,
                            clump_r2 = 0.1,
                            clump_p1 = 1,
                            clump_p2 = 1,
                            pop = "EUR")


eQTLgen_clump$statistics <- (eQTLgen_clump$beta.exposure/eQTLgen_clump$se.exposure)^2  #F-statistics calculation

#PsychENCODE data cleaning####
setwd()
library(tidyverse)
library(dplyr)
library(TwoSampleMR)
library(data.table)

psych <- fread(file = "DER-08a_hg19_eQTL.significant.txt", header = T)
psych$gene_id_clean <- sub("\\..*", "", psych$gene_id)
head(psych)
AF <- fread(file = "SNP_Information_Table_with_Alleles.txt", header = T)    #Allele information
colnames(AF)[colnames(AF) == "PEC_id"] <- "SNP_id"
length(unique(AF$SNP_id))
head(AF)
psych_AF <- inner_join(psych, AF, by="SNP_id")
head(psych_AF)
psych_AF <- psych_AF %>% dplyr::select(Rsid, everything())
max(psych_AF$SNP_distance_to_TSS)
min(psych_AF$SNP_distance_to_TSS)      #Less than 1 Mb
psych_AF$SE <- sqrt(((psych_AF$regression_slope)^2)/qchisq(psych_AF$nominal_pval, 1, lower.tail = F))
library(biomaRt)
hsmart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")
mapping <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = psych_AF$gene_id_clean,
  mart = hsmart)
druggable_gene <- read_xlsx("Druggable gene.xlsx")
druggable_gene <- druggable_gene[,1]
head(mapping)
head(merge)
merge <- merge[,1]
merge <- as.data.frame(merge)
colnames(merge) <- "gene_id_clean"
merge<-merge[!duplicated(merge$gene_id_clean),]
merge<-merge[!duplicated(merge$Symbol),]
head(psych_AF_merge)
colnames(mapping)[colnames(mapping) == "hgnc_symbol"] <- "Symbol"
merge <- inner_join(mapping, druggable_gene, by="Symbol")
colnames(merge)[colnames(merge) == "ensembl_gene_id"] <- "gene_id_clean"
psych_AF_merge <- merge(merge, psych_AF, by="gene_id_clean")
length(unique(psych_AF_merge$gene_id_clean))
length(unique(psych_AF_merge$Symbol))
length(unique(merge$gene_id_clean))
length(unique(mapping$ensembl_gene_id))
length(unique(psych_AF$gene_id_clean))

psych_AF_merge <- left_join(psych_AF_merge, merge, by="gene_id_clean")
psych_AF_merge_sig <- subset(psych_AF_merge, psych_AF_merge$nominal_pval<5e-08)   #P-value < 5e-08
length(unique(psych_AF_merge_sig$gene_id_clean))
length(unique(psych_AF_merge_sig$Symbol.y))
identical(psych_AF_merge_sig$Symbol.x, psych_AF_merge_sig$Symbol.y)
write.csv(psych_AF_merge_sig, file = "psych_AF_merge_sig.csv")
write.csv(psych_AF_merge, file = "psych_AF_merge.csv")

psych_AF_merge_sig <- read.csv(file = "psych_AF_merge_sig.csv", header = T)
length(unique(psych_AF_merge_sig$gene_id_clean))
head(psych_AF_merge_sig)
psych_AF_merge <- read.csv(file = "psych_AF_merge.csv", header = T)
length(unique(psych_AF_merge$gene_id_clean))
psych_AF_exp <- read_exposure_data(filename="psych_AF_merge_sig.csv",
                                  clump = FALSE,
                                  sep = ",",
                                  phenotype_col = "Symbol.x",
                                  snp_col = "Rsid",
                                  beta_col = "regression_slope",
                                  se_col = "SE",
                                  eaf_col = "",
                                  effect_allele_col = "ALT",
                                  other_allele_col = "REF",
                                  pval_col = "nominal_pval",
                                  units_col = "units",
                                  ncase_col = "ncase",
                                  ncontrol_col = "ncontrol",
                                  samplesize_col = "NrSamples",
                                  gene_col = "",
                                  id_col = "Symbol.x",
                                  min_pval = 1e-200,
                                  log_pval = FALSE,
                                  chr_col = "chr",
                                  pos_col = "position")
write.csv(psych_AF_exp, file = "psych_AF_exp.csv")
psych_AF_clump <-  clump_data(dat=psych_AF_exp,       #Clumping procedures
                              clump_kb = 1000,
                              clump_r2 = 0.1,
                              clump_p1 = 1,
                              clump_p2 = 1,
                              pop = "EUR")
length(unique(psych_AF_clump$exposure))
write.csv(psych_AF_clump, file = "psych_AF_clump.csv")
psych_AF_clump <- read.csv(file = "psych_AF_clump.csv")
psych_AF_clump$samplesize.exposure <- 1387
psych_AF_clump$F <- (psych_AF_clump$beta.exposure/psych_AF_clump$se.exposure)^2
class(psych_AF_clump$samplesize.exposure)
#Outcome data reading####
setwd()
library(data.table)
library(tidyverse)
library(TwoSampleMR)
library(readxl)
CVT <- fread(file = "finngen_R12_I9_THROMBICV.txt", header = T)
head(CVT)
eQTLgen_clump <- read.csv(file = "eQTLgen_clump.csv", header = T)

psych_AF_clump <- read.csv(file = "psych_AF_clump.csv", header = T)

CVT_outcome <- read_outcome_data( filename="finngen_R12_I9_THROMBICV.txt",
                                  snps = eQTLgen_clump$SNP,
                                  sep = "\t",
                                  phenotype_col = "Phenotype",
                                  snp_col = "rsids",
                                  beta_col = "beta",
                                  se_col = "sebeta",
                                  eaf_col = "af_alt",
                                  effect_allele_col = "alt",
                                  other_allele_col = "ref",
                                  pval_col = "pval",
                                  units_col = "",
                                  ncase_col = "ncase",
                                  ncontrol_col = "ncontrol",
                                  samplesize_col = "samplesize",
                                  gene_col = "gene",
                                  id_col = "id",
                                  min_pval = 1e-200,
                                  log_pval = FALSE,
                                  chr_col = "chrom",
                                  pos_col = "pos")
CVT_outcome$outcome <- "CVT"
CVT_outcome$samplesize.outcome <- 454144
class(CVT_outcome$samplesize.outcome)
eQTLgen_CVT <- harmonise_data(eQTLgen_clump, CVT_outcome)
Results_eQTLgen_CVT <- mr(eQTLgen_CVT)


#Forest Plot####
colnames(Results_IVW_sig)
forest_data <- Results_IVW_sig[ ,c(4, 6, 7, 8, 9, 10, 11)]
colnames(forest_data)[1] <- "Gene"
colnames(forest_data)[2] <- "Method"
colnames(forest_data)[3] <- "nSNPs"
colnames(forest_data)[6] <- "P-value"
colnames(forest_data)[7] <- "Adjusted P-value"
forest_data$Method <- "IVW"
forest_data$OR <- exp(forest_data$b)
forest_data$OR_lci <- exp(forest_data$b - 1.96*forest_data$se)
forest_data$OR_uci <- exp(forest_data$b + 1.96*forest_data$se)
forest_data$`P-value` <- format(forest_data$`P-value`, scientific = TRUE, digits = 3)
forest_data$`Adjusted P-value` <- format(forest_data$`Adjusted P-value`, scientific = TRUE, digits = 3)
range(forest_data$OR_lci)
range(forest_data$OR_uci)
library(dplyr)
library(magrittr)
forest_data$OR <- as.numeric(forest_data$OR)
forest_data$OR_lci <- as.numeric(forest_data$OR_lci)
forest_data$OR_uci <- as.numeric(forest_data$OR_uci)
forest_data <- forest_data[ ,-11]
forest_data$'Forest plot'<-paste(rep(" ",20), collapse=" ")
forest_data$'OR (95%CI)'<-ifelse(is.na(forest_data$OR),"",sprintf('%.3f (%.3f-%.3f)',
                                                               forest_data$OR, forest_data$OR_lci,
                                                               forest_data$OR_uci))

library(readxl)
library(forestploter)
library(grid)

tm<-forest_theme(data_size=10,
                 ci_pch=20,
                 ci_col = "#4575b4",
                 ci_lty = 1,
                 ci_lwd = 2.3,
                 ci_Theight = 0,
                 refline_lwd = 1.5,
                 refline_lty = "dashed",
                 refline_col = "red",
                 summary_fill = "#4575b4",
                 summary_col = "#4575b4",
                 footnote_cex = 1.1,
                 footnote_fontface = "italic",
                 footnote_col = "blue")

p_F<-forest(forest_data[,c(1:3, 11:12, 6:7)],est = forest_data$OR,lower = forest_data$OR_lci,
           upper = forest_data$OR_uci, sizes = 0.6, ci_column = 4,
           ref_line =1, xlim = c(0.2, 4), ticks_at = c(0.2, 1, 2, 3, 4),
           arrow_lab = c('protective factor','risk factor'),
           footnote = 'p<0.05 was considered statistically significant',
           theme=tm)
p_F

#Manhattan plot####
install.packages("CMplot")
install.packages("qqman")
library(CMplot)
library(qqman)
Results_IVW <- read.csv(file = "Results_eQTLgen_CVT_IVW.csv")
eQTLgen <- fread(file = "2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt",
                 header = T)
colnames(eQTLgen)
colnames(Results_IVW)
eQTLgen <- eQTLgen[ ,c(9, 10, 11)] 
Results_IVW <- Results_IVW[ ,c(5, 10, 11)]
colnames(eQTLgen)[1] <- "exposure"
library(dplyr)
eQTLgen <- distinct(eQTLgen, GeneSymbol, .keep_all = T)
Mydata <- inner_join(Results_IVW, eQTLgen, by = "exposure")
colnames(Mydata)
head(Mydata)
Mydata <- Mydata[,c(1, 4, 5, 2, 3)]
Mydata <- Mydata[,-5]
colnames(Mydata) <- c("SNP", "CHR", "BP", "P")
mycol <- c("#3399ff", "#99cc33", "#ffcc00")
p <- manhattan(Mydata, chr = "CHR", bp = "BP", p = "P", snp = "SNP", 
               logp = T, highlight = T, ylim = c(0, 10), tck = -0.02, 
               cex = 0.7, cex.axis = 0.5, cex.lab = 0.7, col = mycol,
               suggestiveline = F, genomewideline = F, chrlabs = c(1: 20, "P", "Q"))
manhattan(Mydata, highlight = snpsOfInterest)
manhattan(Mydata, annotatePval = 3.78e-04)
Mydata$Gene <- Mydata[, -5]
CMplot(Mydata, plot.type = "m", threshold = 3.78e-04, threshold.col=c('grey','black'),       
       threshold.lty = c(1,2), threshold.lwd = c(1,1), amplify = T, signal.cex = c(1,1), 
       signal.pch = c(20,20),signal.col = c("red","orange"), highlight = text$SNP,
       highlight.text = text$Gene, file.name = "re1_plus_gene", file="pdf")
dev.off

#Volcano plot####
if(!require("cowplot")) install.packages("cowplot")
library("cowplot")
color_scale <- scale_colour_gradientn(
  colours = c("#1681D8", "green","yellow", "#CE1256"),
  values = c(0, 0.4,0.7, 1),
  guide = guide_colorbar(title = "logp"))
max_size <- 4
min_size <- 1
Vol_mydata$Size <- sqrt(min_size + (max_size - min_size) * (Vol_mydata$logp - min(Vol_mydata$logp)) / (max(Vol_mydata$logp) - min(Vol_mydata$logp)))
Vol_mydata$gene_name <- Vol_mydata$exposure
Vol_mydata <- Vol_mydata %>% mutate(label = ifelse(gene_name %in% c("BMPR2","IL18","COMT", 
                                                                    "C4B","SERPINC1","ST14","ACP2","GRIK4", 
                                                                    "TNXB","RPN1","APOBEC3A","THBS1", 
                                                                    "PTGDR", "GYPC", "H1F0", "RELT", "CXCR6", 
                                                                    "KLHL36", "FCRL5"), gene_name, ""))


p <- ggplot(Vol_mydata, aes(x = b, y = logp)) +
  geom_point(aes(color = logp), size = 3) +
  color_scale +
  labs(x = "Beta", y = "-log(P-vale)") +
  theme(legend.position = c(0.02, 0.88))+
  theme_cowplot()+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1), 
        legend.position = c(0.02, 0.6)) + 
  geom_hline(yintercept = -log10(3.78e-04), linetype = "dotted", color = "black")+
  geom_text(aes(x = 2, y = 3.3, label = "FDR=0.05"), vjust = 1, hjust = -0.2,size = 4.5)+
  guides(size = "none") + 
  geom_label_repel(data = Vol_mydata, aes(label = label), color = "black", 
                   size = 3) 


p


#Brain-derived genes on risk of CVT####
psych_AF_clump <- read.csv(file = "psych_AF_clump.csv")
length(unique(psych_AF_clump$exposure))
psych_CVT_outcome <- read_outcome_data( filename="finngen_R12_I9_THROMBICV.txt",
                                  snps = psych_AF_clump$SNP,
                                  sep = "\t",
                                  phenotype_col = "Phenotype",
                                  snp_col = "rsids",
                                  beta_col = "beta",
                                  se_col = "sebeta",
                                  eaf_col = "af_alt",
                                  effect_allele_col = "alt",
                                  other_allele_col = "ref",
                                  pval_col = "pval",
                                  units_col = "",
                                  ncase_col = "ncase",
                                  ncontrol_col = "ncontrol",
                                  samplesize_col = "samplesize",
                                  gene_col = "gene",
                                  id_col = "id",
                                  min_pval = 1e-200,
                                  log_pval = FALSE,
                                  chr_col = "chrom",
                                  pos_col = "pos")
psych_CVT_outcome$outcome <- "CVT"
psych_CVT_outcome$samplesize.outcome <- 454144
write.csv(psych_CVT_outcome, file = "psych_CVT_outcome.csv")
psych_CVT <- harmonise_data(psych_AF_clump, psych_CVT_outcome)
write.csv(psych_CVT, file = "psych_CVT_harm.csv")
Results_psych_CVT <- mr(psych_CVT)
write.csv(Results_psych_CVT, file = "Results_psych_CVT.csv")
Psych_selec <- Results_psych_CVT[Results_psych_CVT$method %in% c("Inverse variance weighted", "Wald ratio"), ]
write.csv(Psych_selec, file = "Psych_selec.csv")
Psych_selec$Ajusted_P <- p.adjust(Psych_selec$pval, method = "BH")
write.csv(Psych_selec, file = "Psych_selec.csv")

Psych_eQTL_CVT_IVW <- read.csv(file = "Psych_selec.csv")
head(Psych_eQTL_CVT_IVW)
Psych_eQTL_CVT_IVW$Ajusted_P <- p.adjust(Psych_eQTL_CVT_IVW$pval, method = "BH")
write.csv(Psych_eQTL_CVT_IVW, file = "Psych_eQTL_CVT_IVW.csv")

#Manhattan Plot####
Psych_eQTL <- fread(file = "DER-08a_hg19_eQTL.significant.txt")
Psych_eQTL_symbol <- fread(file = "psych_AF_merge.csv")
colnames(Psych_eQTL)
colnames(Psych_eQTL_symbol)
head(Psych_eQTL_symbol)
identical(Psych_eQTL_symbol$Symbol.x, Psych_eQTL_symbol$Symbol.y)
Psych_eQTL_symbol <- Psych_eQTL_symbol[ ,c(3, 6, 7)]
Psych_eQTL_symbol <- distinct(Psych_eQTL_symbol, Symbol.x, .keep_all = T)
colnames(Psych_eQTL_symbol)[1] <- "exposure"
colnames(Psych_eQTL_CVT_IVW)
Results_psy_IVW <- Psych_eQTL_CVT_IVW[ ,c(5, 8, 10)]
Mydata_psy <- inner_join(Results_psy_IVW, Psych_eQTL_symbol, by = "exposure")
head(Mydata_psy)
colnames(Mydata_psy) <- c("SNP", "Beta", "Pvalue", "CHR", "BP")
Mydata_psy <- Mydata_psy[ ,c(1, 4, 5, 2, 3)]
Mydata_psy$CHR <- gsub("chr", "", Mydata_psy$CHR)
Mydata_psy <- Mydata_psy[ , -4]
text <- data_frame(SNP = "ZP3",
                   Gene = "ZP3")
CMplot(Mydata_psy, plot.type = "m", threshold = 4.39e-05, threshold.col=c('grey','black'),       
       threshold.lty = c(1,2), threshold.lwd = c(1,1), amplify = T, signal.cex = c(1,1), 
       signal.pch = c(20,20),signal.col = c("red","orange"), highlight = text$Gene,
       highlight.text = text$Gene, file.name = "Psy_mhantan", file="pdf")

colnames(Psych_eQTL_CVT_IVW)
Vol_Mydata_Psy <- Psych_eQTL_CVT_IVW[ ,c(5, 8, 10)]
Vol_Mydata_Psy$logp <-  -log10(Vol_Mydata_Psy$pval)
range(Vol_Mydata_Psy$b)



color_scale <- scale_colour_gradientn(
  colours = c("#1681D8", "green","yellow", "#CE1256"),
  values = c(0, 0.4,0.7, 1),
  guide = guide_colorbar(title = "logp"))
max_size <- 4
min_size <- 1

Vol_Mydata_Psy$gene_name <- Vol_Mydata_Psy$exposure
Vol_Mydata_Psy <- Vol_Mydata_Psy %>% mutate(label = ifelse(gene_name %in% c("ZP3"), gene_name, ""))

p <- ggplot(Vol_Mydata_Psy, aes(x = b, y = logp)) +
  geom_point(aes(color = logp), size = 3) +
  color_scale +
  labs(x = "Beta", y = "-log(P-vale)") +
  theme(legend.position = c(0.02, 0.88))+
  theme_cowplot()+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1), 
        legend.position = c(0.02, 0.6)) + 
  geom_hline(yintercept = -log10(4.39e-05), linetype = "dotted", color = "black")+
  geom_text(aes(x = 1, y = 4.2, label = "FDR=0.05"), vjust = 1, hjust = -0.2,size = 4.5)+
  guides(size = "none") + 
  geom_label_repel(data = Vol_Mydata_Psy, aes(label = label), color = "black", 
                   size = 3) 

p



#colocalization#####
setwd("D:/CVT药靶/CVTdrug/Outcome")

library(data.table)
CVT<-fread(file = "finngen_R12_I9_THROMBICV.txt")
head(CVT)


CVT<-CVT[, c("rsids", "alt", "ref", "af_alt", "beta", "sebeta", 
          "pval")]
CVT$sample.size <- 454144
colnames(CVT)<-c("SNP", "A1", "A2", "freq", "b", "se", "p", "n")
str(CVT)
f$SNP<-as.numeric(f$SNP)
f$A1<-as.numeric(f$A1)
f$A2<-as.numeric(f$A2)
write.table(f, file = "IA.txt", sep =" ", row.names = FALSE, quote = FALSE)
str(f)
dele<-which(CVT$SNP=="")
CVT<-CVT[-dele, ]
grepl(",", df$SNP)
grepl(",", CVT$SNP)
CVT<-CVT[!grepl(",", CVT$SNP), ]
CVT$SNP
duplicated(CVT$SNP)
CVT<-CVT[!duplicated(CVT$SNP), ]
write.table(CVT, file = "CVT.txt",row.names = F,sep = "\t",quote = F)
write.table(f, gzfile("IA.txt.gz"), sep = "\t",row.names = F,quote = F)
save.image()

setwd()
CVT <- fread(file = "finngen_R12_I9_THROMBICV.txt")
install.packages("vcfR")
library(vcfR)
SERPINC1 <- read.vcfR("eqtl-a-ENSG00000093010.vcf.gz")
library(tidyverse)
gt <- SERPINC1@gt
gt <- as.data.frame(gt)
gt <- separate(gt, col = "eqtl-a-ENSG00000093010", into=
                 c("ES", "SE", "LP", "AF", "SS", "ID"), 
               sep = ":")
gt <- na.omit(gt)
colnames(gt) <- c("format", "beta", "se", "logpvalue", 
                  "eaf", "samplesize", "snp")

gt$beta <- as.numeric(gt$beta)
gt$se <- as.numeric(gt$se)
gt$logpvalue <- as.numeric(gt$logpvalue)
gt$eaf <- as.numeric(gt$eaf)
gt$samplesize <- as.numeric(gt$samplesize)
gt$format <- NULL

fix <- SERPINC1@fix
fix <- as.data.frame(fix)
colnames(fix) <- c("chr", "pos", "snp", "ref", "alt")
fix <- fix[,c(1:5)]

SERPINC1_eqtl <- left_join(fix, gt, by="snp")
SERPINC1_eqtl <- na.omit(SERPINC1_eqtl)
SERPINC1_eqtl$maf <- ifelse(SERPINC1_eqtl$eaf < 0.5, 
                            SERPINC1_eqtl$eaf, 
                            1-SERPINC1_eqtl$eaf)
SERPINC1_eqtl$eaf <- NULL

#rs6435149
#chr: pos=22: 19921978
SERPINC1_eqtl <- SERPINC1_eqtl[SERPINC1_eqtl$chr==22, ]
SERPINC1_eqtl$logpvalue <- as.numeric(SERPINC1_eqtl$logpvalue)
SERPINC1_eqtl$p_value <- 10^(-SERPINC1_eqtl$logpvalue)
SERPINC1_eqtl$pos <- as.numeric(SERPINC1_eqtl$pos)
SERPINC1_eqtl <- SERPINC1_eqtl[SERPINC1_eqtl$pos > 19921978-100000, ]
SERPINC1_eqtl <- SERPINC1_eqtl[SERPINC1_eqtl$pos < 19921978+100000, ]
my_eqtl <- SERPINC1_eqtl[,c("snp", "p_value", "maf")]
colnames(my_eqtl) <- c("snp", "pvalues", "MAF")
my_eqtl <- na.omit(my_eqtl)
my_eqtl <- my_eqtl[my_eqtl$MAF>0, ]

library(TwoSampleMR)
CVT_outcome <- read_outcome_data( filename="finngen_R12_I9_THROMBICV.txt",
                                  snps = my_eqtl$snp,
                                  sep = "\t",
                                  phenotype_col = "Phenotype",
                                  snp_col = "rsids",
                                  beta_col = "beta",
                                  se_col = "sebeta",
                                  eaf_col = "af_alt",
                                  effect_allele_col = "alt",
                                  other_allele_col = "ref",
                                  pval_col = "pval",
                                  units_col = "",
                                  ncase_col = "ncase",
                                  ncontrol_col = "ncontrol",
                                  samplesize_col = "samplesize",
                                  gene_col = "gene",
                                  id_col = "id",
                                  min_pval = 1e-200,
                                  log_pval = FALSE,
                                  chr_col = "chrom",
                                  pos_col = "pos")
CVT_outcome$beta.outcome <- as.numeric(CVT_outcome$beta.outcome)
CVT_outcome$se.outcome <- as.numeric(CVT_outcome$se.outcome)
head(CVT_outcome)
CVT_outcome$varbate <- CVT_outcome$se.outcome^2
CVT_outcome <- CVT_outcome[ ,c("SNP", "pval.outcome",
                               "beta.outcome", "varbate")]
colnames(CVT_outcome) <- c("snp", "pvalues", "beta", "varbeta")
CVT_outcome <- CVT_outcome[!duplicated(CVT_outcome$snp), ]
my_eqtl <- my_eqtl[my_eqtl$snp %in% CVT_outcome$snp, ]
my_eqtl=as.list(my_eqtl)
my_eqtl[["type"]]="quant"
my_eqtl[["N"]]=31036

CVT_outcome=as.list(CVT_outcome)
CVT_outcome[["type"]]="cc"
CVT_outcome[["N"]]=454144

library(coloc)
coloc.abf(dataset1 = CVT_outcome, dataset2 = my_eqtl)
list(rm=ls())
#Visualization of colocalization####
library(dplyr)
gwas_fn <- CVT_outcome[ ,c("snp", "pvalues")] %>% 
  rename(rsid=snp, pval=pvalues)
eqtl_fn <- my_eqtl[ ,c("snp", "pvalues")] %>% 
  rename(rsid=snp, pval=pvalues)
library(locuscomparer)
print(locuscompare(in_fn1 = gwas_fn,
                   in_fn2 = eqtl_fn, 
                   title1 = "GWAS", 
                   title2 = "eQTL"))
