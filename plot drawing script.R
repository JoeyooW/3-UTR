#Reformatting for the raw data
#compare the tool efficiency for detection in 3'UTR region
#version 1 - 23.01.2025

#set up my environment
library("tidyverse")
library('ggplot2')
library('dtplyr')
library('dplyr')
library(tidyr)
library(ggpubr)
#library(pROC)
library(broom)

#------------get the VCF format for clinic and control dataset-------------------

#load my data
setwd('L:/Lab_MandyS/Joey Wu/data')
origin_VCF = read.delim('/Users/JOEY/Desktop/UQ/Research/data/240215_cisreg_3p_FILTER.txt', header = TRUE, sep = "\t")


#step by step plan of what plan to do. 
#add a '.' column
origin_VCF$dbSNP <- c(".")
summary(as.vector(origin_VCF))


#select the information we need
summary(origin_VCF)
selected_VCF <- origin_VCF %>% select('CHROM','POS','contsetid','REF','ALT') 
#colomn check and chromosome check
colnames(selected_VCF)
summary(as.factor(selected_VCF$CHROM))
summary(as.factor(origin_VCF$CHROM))
colnames(selected_VCF) = c("CHROM", "POS", "ID", "REF", "ALT")
#add the rest colomn'QUAl''FILTER''INFO'
TEST_VCF = selected_VCF %>% 
  add_column(QUAL = ".") %>%
  add_column(FILTER = ".") %>%
  add_column(INFO = ".")
# %>% means and then, pipeline


#load data
write.table(TEST_VCF, "TEST_VCF Both.txt", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
# sep means the field separator string

##-----------------box and whisker plot for CADD score-------------
#load the data
setwd('L:/Lab_MandyS/Joey Wu/data')
CADD_VCF = read.delim('GRCh38-v1.7_a834707af8ad410e4f350dd6f12db931-CADD.txt', header = TRUE, sep = "\t")
colnames(CADD_VCF) <- c('CHROM','POS','REF','ALT','RawScore','PHRED')
summary(as.vector(CADD_VCF))

#matching the score with ID and Group(C/D)
##systhesis the unique identify ID
CADD_VCF$CHROM <- gsub("chr", "", CADD_VCF$CHROM)
CADD_VCF$UIcode = gsub(" ", "", paste(CADD_VCF$CHROM,"-",CADD_VCF$POS,"-",CADD_VCF$REF,"-",CADD_VCF$ALT))
summary(as.vector(CADD_VCF$UIcode))

#left join for the table
#Score_VCF %>% left_join(Score_VCF,origin_VCF,by = c("UIcode", "variant38"))
#Score_VCF %>% left_join(origin_VCF, by = c("UIcode" = "variant38"))

Pic_VCF_CADD <- CADD_VCF %>%
  left_join(origin_VCF %>% select(UIcode = variant38, INFO, contsetid), 
            by = c("UIcode" = "UIcode"))

write.table(Pic_VCF_CADD, file = "CADD_Likelihood.txt", sep = "\t", row.names = FALSE, quote = FALSE)


#box plot
ggplot(data = Pic_VCF_CADD, mapping = aes(x = INFO, y = PHRED)) +
  geom_boxplot() 

ggplot(data = Pic_VCF_CADD, mapping = aes(x = INFO, y = PHRED)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.4) +
  theme_minimal() +
  labs(title = "Disease Impact Score by Variant Group",
       x = "Variant Group", y = "CADD Score") +
  theme(legend.position = "none")

#point plot
ggplot(Pic_VCF_CADD, aes(CHROM, PHRED)) +  geom_point() 

#violin plot + distribution + point +box
ggplot(Pic_VCF_CADD, aes(x = INFO, y = PHRED)) + geom_violin(aes(fill = INFO), size = 0.3, bw = .15, color = NA) + 
  scale_fill_manual(values = c("lightcyan2", "mistyrose1")) +
  geom_boxplot(
    width = 0.1, fill = "white",
    size = 0.3, outlier.shape = NA
  ) +
  ggdist::stat_halfeye(
    adjust = .33, 
    width = .55, 
    color = NA, 
    position = position_nudge(x = 0.1)  # Adjust nudging value for proper positioning
  ) +
  geom_point(
    position = position_nudge(x = -0.1),  # Adjust nudging value for proper positioning
    shape = 95, size = 4, alpha = 0.25  
  ) +  theme_minimal() 
  

  #density plot
  #mean
summary(as.factor(Pic_VCF_CADD$INFO))
Pic_VCF_CADD %>%  group_by(INFO) %>% summarise_at(vars(PHRED), list(name = mean))
#  INFO   name
#  <chr> <dbl>
#  1 Cont   5.85
#  2 Dis   11.3 


# Change fill color by sex and add mean line
# Use semi-transparent fill: alpha = 0.4
ggplot(Pic_VCF_CADD, aes(x=PHRED, color=INFO)) + geom_density(aes(fill = INFO), alpha = 0.4) +
  geom_vline(xintercept = 5.85, color = "#868686FF", linetype = "dashed") +
  geom_vline(xintercept = 11.3, color = "#868686FF", linetype = "dashed") +
  scale_fill_manual(values = c("cyan", "red"),na.value = "darkslategray") +
  theme(legend.position="none", text = element_text(size= 12))


ggplot(Pic_VCF_CADD, aes(x=PHRED, color=INFO)) + 
  geom_density(alpha = 0.4) +
  theme_minimal() +
  labs(title = "Score Distribution by Group", x = "Score", y = "Density")


#-----------------box and whisker plot for REMM pre_score----------


#get sorced VCF for the 3' UTR database(contain both clinical and common variant) and make box plot
#version 1 - 24.01.2025

#set the environment


#load the data
setwd('L:/Lab_MandyS/Joey Wu/data')
Score_VCF = read.delim('v0.4.hg38_3f9525378ef8bd631df26ac045d65d76-REMM-computed scores.tsv', header = FALSE, sep = "\t")
colnames(Score_VCF) <- c('CHROM','POS','REF','ALT','REMM_pre')
summary(as.vector(Score_VCF))

#matching the score with ID and Group(C/D)
##systhesis the unique identify ID
Score_VCF$CHROM <- gsub("chr", "", Score_VCF$CHROM)
Score_VCF$UIcode = gsub(" ", "", paste(Score_VCF$CHROM,"-",Score_VCF$POS,"-",Score_VCF$REF,"-",Score_VCF$ALT))
summary(as.vector(Score_VCF$UIcode))

#left join for the table
#Score_VCF %>% left_join(Score_VCF,origin_VCF,by = c("UIcode", "variant38"))
#Score_VCF %>% left_join(origin_VCF, by = c("UIcode" = "variant38"))

Pic_VCF_REMM <- Score_VCF %>%
  left_join(origin_VCF %>% select(UIcode = variant38, INFO, contsetid), 
            by = c("UIcode" = "UIcode"))
write.table(Pic_VCF_REMM, file = "REMM_Likelihood.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#box plot
ggplot(data = Pic_VCF_REMM, mapping = aes(x = INFO, y = REMM_pre)) +
  geom_boxplot() 


ggplot(data = Pic_VCF_REMM, mapping = aes(x = INFO, y = REMM_pre)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.4) +
  theme_minimal() +
  labs(title = "Disease Impact Score by Variant Group",
       x = "Variant Group", y = "REMM Score") +
  theme(legend.position = "none")

#point plot
ggplot(Pic_VCF_REMM, aes(CHROM, REMM_pre)) +  geom_point() 

#violin plot + distribution + point +box
ggplot(Pic_VCF_REMM, aes(x = INFO, y = REMM_pre)) + geom_violin(aes(fill = INFO), size = 0.3, bw = .15, color = NA) + 
  scale_fill_manual(values = c("lightcyan2", "mistyrose1")) +
  geom_boxplot(
    width = 0.1, fill = "white",
    size = 0.3, outlier.shape = NA
  ) +
  ggdist::stat_halfeye(
    adjust = .33, 
    width = .55, 
    color = NA, 
    position = position_nudge(x = 0.1)  # Adjust nudging value for proper positioning
  ) +
  geom_point(
    position = position_nudge(x = -0.1),  # Adjust nudging value for proper positioning
    shape = 95, size = 4, alpha = 0.25  
  ) +  theme_minimal() 
  

  #density plot
ggplot(Pic_VCF_REMM, aes(x=REMM_pre, color=INFO)) +
  geom_density() + geom_vline(aes(xintercept = grp.mean, color = INFO), linetype = "dashed") +scale_fill_grey()

#density plot
#mean
summary(as.factor(Pic_VCF_REMM$INFO))
Pic_VCF_REMM %>%  group_by(INFO) %>% summarise_at(vars(REMM_pre), list(name = mean))
#  INFO   name
#  <chr> <dbl>
#  1 Cont 0.528
#  2 Dis  0.766


# Change fill color by sex and add mean line
# Use semi-transparent fill: alpha = 0.4
ggplot(Pic_VCF_REMM, aes(x=REMM_pre, color=INFO)) + geom_density(aes(fill = INFO), alpha = 0.4) +
  geom_vline(xintercept = 0.528, color = "#868686FF", linetype = "dashed") +
  geom_vline(xintercept = 0.766, color = "#868686FF", linetype = "dashed") +
  scale_fill_manual(values = c("cyan", "red"),na.value = "darkslategray") +
  theme(legend.position="none", text = element_text(size= 12))


ggplot(Pic_VCF_REMM, aes(x=REMM_pre, color=INFO)) + 
  geom_density(alpha = 0.4) +
  theme_minimal() +
  labs(title = "Score Distribution by Group", x = "Score", y = "Density")

                                                               
              
#-----------------get the RegVar format file------------

#import VEP file
setwd('L:/Lab_MandyS/Joey Wu/data')
VEP_annotation = read.delim('VEP annotation.txt', header = TRUE, sep = "\t")
summary(as.vector(VEP_annotation))

VEP_annotation$Uploaded_variation <- VEP_annotation$X.Uploaded_variation

#check the distinct(not necessary)#
VEP_annotation %>%
  group_by(Uploaded_variation) %>%
  summarise(
    check = ifelse(n_distinct(Existing_variation) == 1, "Yes", "No")
  ) %>%  filter(check == "No") %>%
  pull(Uploaded_variation) %>%
  print()

# ^^^ output YES, for a variant just a dbSNP id annotation

#extract the dbSNP id
VEP_annotation_NoDup <- VEP_annotation
VEP_annotation_NoDup$dbSNP <- str_extract(VEP_annotation_NoDup$Existing_variation, "^rs[^,]+")





#match with Marrisa file and find the corresponding gene
gene_VCF = read.delim('L:/Lab_MandyS/Joey Wu/data/3prime_variants_test for VCF forming.csv', header = TRUE, sep = ",")
summary(gene_VCF)
summary(as.factor(origin_VCF$INFO))
summary(as.factor(gene_VCF$Identifier))
summary(as.factor(VEP_annotation_NoDup$X.Uploaded_variation))
duplicated(gene_VCF$Identifier)
duplicated(gene_VCF$GRCh38_POS)

VEP_annotation_NoDup$MC_ID <- gsub("^Mcont", "", VEP_annotation_NoDup$Uploaded_variation)

VEP_annotation_NoDup$Associated_gene <- gene_VCF$Associated_gene[match(VEP_annotation_NoDup$MC_ID, gene_VCF$MC_ID)]

#remove annotation records with unmatched gene for clincic set
VEP_annotation_NoDup$yes_or_no <- ifelse(VEP_annotation_NoDup$SYMBOL == VEP_annotation_NoDup$Associated_gene, "yes", "no")

VEP_annotation_NoDup <- VEP_annotation_NoDup %>%
  filter(yes_or_no == "yes" | is.na(yes_or_no))


#^^^------already get 168 annotated clinic variant(same number with 240215_cisreg_3p_FILTER)

#select the control set annotation
VEP_annotation_NoDup <- VEP_annotation_NoDup %>%
  mutate(
    yes_or_no_control = ifelse(is.na(yes_or_no) & BIOTYPE == "protein_coding" & Consequence == "3_prime_UTR_variant", "yes", NA)
  )

VEP_annotation_NoDup <- VEP_annotation_NoDup %>%
  filter(
    !is.na(yes_or_no) | !is.na(yes_or_no_control))


write.table(VEP_annotation_NoDup, file = "3.txt", sep = "\t",row.names = FALSE, quote = FALSE)


#not use # remove duplicate rows(need to be after the matching, need to remove whole duplicate row)
VEP_annotation_NoDup <- VEP_annotation_NoDup %>%
  distinct(Uploaded_variation, Existing_variation,SYMBOL, .keep_all = TRUE)
#keep unique dbSNP id and Mcontid combination, and with multi gene name after annotation
summary(VEP_annotation_NoDup)




#format
TEST_VCF_RegVar <- VEP_annotation_NoDup[,c("dbSNP","SYMBOL")]

write.table(TEST_VCF_RegVar, file = "TEST_VCF_RegVar.txt", sep = "\t", row.names = FALSE, quote = FALSE)


#-----------------box and whisker plot for RegVar score for liver tissue---------
#load the data
setwd('L:/Lab_MandyS/Joey Wu/data')
RegVar_VCF = read.delim('/Users/JOEY/Desktop/UQ/Research/data/RegVar/regvar-result-1e0bdfaf8af54c5eb59a28a1a24f2570.csv', header = TRUE, sep = ",")
summary(as.vector(RegVar_VCF))

#retrieve the regvar data to   the clinic and control set
#choose the shortest distance


RegVar_VCF <- RegVar_VCF %>%
  group_by(SNP) %>%
  slice_min(Distance, n = 1, with_ties = FALSE) %>%  
  ungroup()
write.table(RegVar_VCF, file = "RegVar_VCF.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# matching gene and snp and get ID 
RegVar_VCF <- RegVar_VCF %>%
  mutate(key = paste(SNP, Gene, sep = "_"))

VEP_annotation_NoDup <- VEP_annotation_NoDup %>%
  mutate(key = paste(dbSNP, SYMBOL, sep = "_"))
colnames(VEP_annotation_NoDup)

result <- RegVar_VCF %>%
  left_join(VEP_annotation_NoDup %>% select(key, Uploaded_variation), by = "key")

print(result)


write.table(result, file = "result.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#GRCH version check: regvar outcome:Grch37, others:Grch38
#dup macthing:ID109(17-8173550-C-G)/110(17-8173550-C-A),ID112(17-8173581-C-G)/327(17-8173581-C-T),ID24(17-8173581-C-T)/25(11-5225488-A-G)
#above dup matching is all missense
#some variant cannot retrieve the ID(all because regvar allocate it with different  gene), already search and add manually 



# the dup variant can be fill in manually after checking the indel or missense
# the NA variant can be fill in maually after checking the version of grch, and 
# also need to check the origin version and vep version 
# the third step should be checking the regvar score, select the online server
# or GUI script to validate
# Finally try to use the different tissue to see will the outcome be different or not
# can try whole blood first

Pic_VCF_REGVAR = read.delim('REGVAR result.csv', header = TRUE, sep = ",")
#box plot
ggplot(data = Pic_VCF_REGVAR, mapping = aes(x = INFO, y = SCORE)) +
  geom_boxplot() 


ggplot(data = Pic_VCF_REGVAR, mapping = aes(x = INFO, y = SCORE, fill = INFO)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.4) +
  theme_minimal() +
  labs(title = "Disease Impact Score by Variant Group",
       x = "Variant Group", y = "REGVAR Score") +
  theme(legend.position = "none")

for #density plot
#mean
summary(as.factor(Pic_VCF_REGVAR$INFO))
Pic_VCF_REGVAR %>%  group_by(INFO) %>% summarise_at(vars(SCORE), list(name = mean))
#  INFO   name
#  <chr> <dbl>
#  1 Cont  0.165
#  2 Dis   0.257

ggplot(Pic_VCF_REGVAR, aes(x = SCORE)) +
  geom_histogram(bins = 30)




# Change fill color by sex and add mean line
# Use semi-transparent fill: alpha = 0.4
ggplot(Pic_VCF_REGVAR, aes(x=log1p(SCORE), color=INFO)) + 
  geom_density(aes(fill = INFO), alpha = 0.4) +
  geom_vline(xintercept = 0.165, color = "#868686FF", linetype = "dashed") +
  geom_vline(xintercept = 0.257, color = "#868686FF", linetype = "dashed") +
  scale_fill_manual(values = c("cyan", "red"), na.value = "darkslategray") +
  theme(legend.position="none", text = element_text(size= 5))

ggplot(Pic_VCF_REGVAR, aes(x=log1p(SCORE), fill=INFO)) + 
  geom_density(alpha = 0.4) +
  theme_minimal() +
  labs(title = "Score Distribution by Group", x = "Score", y = "Density")



#-----------------box and whisker plot for RegVar score for whole blood---------
#load the data
setwd('L:/Lab_MandyS/Joey Wu/data')
RegVar_VCF_wblood = read.delim('/Users/JOEY/Desktop/UQ/Research/data/RegVar/regvar-result-1e0bdfaf8af54c5eb59a28a1a24f2570.csv', header = TRUE, sep = ",")
summary(as.vector(RegVar_VCF_wblood))
colnames(RegVar_VCF_wblood)[1] <- "SNP"
colnames(RegVar_VCF_wblood)[3] <- "Gene"

#retrieve the regvar data to   the clinic and control set
#choose the shortest distance
RegVar_VCF_wblood <- RegVar_VCF_wblood %>%
  group_by(SNP) %>%
  slice_min(Distance, n = 1, with_ties = FALSE) %>%  
  ungroup()


# matching gene and snp and get ID 
RegVar_VCF_wblood <- RegVar_VCF_wblood %>%
  mutate(key = paste(SNP, Gene, sep = "_"))

VEP_annotation_NoDup <- VEP_annotation_NoDup %>%
  mutate(key = paste(dbSNP, SYMBOL, sep = "_"))

result <- RegVar_VCF_wblood %>%
  left_join(VEP_annotation_NoDup %>% select(key, Uploaded_variation), by = "key")

print(result)


write.table(result, file = "result2.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#GRCH version check: regvar outcome:Grch37, others:Grch38
#dup macthing:ID109(17-8173550-C-G)/110(17-8173550-C-A),ID112(17-8173581-C-G)/327(17-8173581-C-T),ID24(17-8173581-C-T)/25(11-5225488-A-G)
#above dup matching is all missense
#some variant cannot retrieve the ID(all because regvar allocate it with different  gene), already search and add manually 



Pic_VCF_REGVA_wblood = read.delim('result2.csv', header = TRUE, sep = ",")
colnames(Pic_VCF_REGVA_wblood)[7] <- 'SCORE'
#box plot
ggplot(data = Pic_VCF_REGVA_wblood, mapping = aes(x = INFO, y = SCORE)) +
  geom_boxplot() 


for #density plot
#mean
summary(as.factor(Pic_VCF_REGVA_wblood$INFO))
Pic_VCF_REGVA_wblood %>%  group_by(INFO) %>% summarise_at(vars(SCORE), list(name = mean))
#  INFO   name
#  <chr> <dbl>
#  1 Cont  0.227
#  2 Dis   0.271

ggplot(Pic_VCF_REGVA_wblood, aes(x = SCORE)) +
  geom_histogram(bins = 30)


# Change fill color by sex and add mean line
# Use semi-transparent fill: alpha = 0.4
ggplot(Pic_VCF_REGVA_wblood, aes(x=log1p(SCORE), color=INFO)) + 
  geom_density(aes(fill = INFO), alpha = 0.4) +
  geom_vline(xintercept = 0.227, color = "#868686FF", linetype = "dashed") +
  geom_vline(xintercept = 0.271, color = "#868686FF", linetype = "dashed") +
  scale_fill_manual(values = c("cyan", "red"), na.value = "darkslategray") +
  theme(legend.position="none", text = element_text(size= 5))









  
  






#-----------------components analysis for CADD--------
setwd("/Users/JOEY/Downloads/CADD  components table plot")
CADD_elements = read_tsv('/Users/JOEY/Desktop/UQ/Research/data/CADD  components table plot/GRCh38-v1.7_anno_a834707af8ad410e4f350dd6f12db931.tsv')

# select 3' UTR only
CADD_elements <- CADD_elements[CADD_elements$Consequence == "3PRIME_UTR", ]

colnames(CADD_elements)[1:4] <- c('CHROM','POS','REF','ALT')

#matching the score with ID and Group(C/D)
##systhesis the unique identify ID
CADD_elements$CHROM <- gsub("chr", "", CADD_elements$CHROM)
CADD_elements$UIcode = gsub(" ", "", paste(CADD_elements$CHROM,"-",CADD_elements$POS,"-",CADD_elements$REF,"-",CADD_elements$ALT))
summary(as.vector(CADD_elements$UIcode))

#left join for the table
#Score_VCF %>% left_join(Score_VCF,origin_VCF,by = c("UIcode", "variant38"))
#Score_VCF %>% left_join(origin_VCF, by = c("UIcode" = "variant38"))
CADD_elements <-CADD_elements %>%
  left_join(origin_VCF %>% select(UIcode = variant38, INFO, contsetid), 
            by = c("UIcode" = "UIcode"))
summary(CADD_elements)
colnames(CADD_elements)[ncol(CADD_elements)] <- "ID"
summary(CADD_elements$ID)
summary(CADD_elements$INFO)

write.csv(CADD_elements, "CADD_scores.csv", row.names = FALSE)


--------------------####---------
CADD_com_analysis <- CADD_elements

# find na column
na_columns <- colnames(CADD_com_analysis)[sapply(CADD_com_analysis, function(col) all(is.na(col)))]

# print NA columns
if (length(na_columns) > 0) {
  print("[below column values are NA, will be removed:")
  print(na_columns)
  
  # remove NA column
  CADD_com_analysis <- CADD_com_analysis[, !colnames(CADD_com_analysis) %in% na_columns]
}

# setting components group
group_col <- "INFO"  # INFO group
components <- setdiff(colnames(CADD_com_analysis), c("Chrom","Pos","Ref","Alt","Type","Length", "AnnoType", "Consequence", "INFO", "RawScore", "PHRED", "ID"))  # components
components <- components[sapply(CADD_com_analysis[components], is.numeric)]

# set a table to storage
results <- data.frame(Component = character(),
                      Mean_Cont = numeric(),
                      Mean_Dis = numeric(),
                      SD_Cont = numeric(),
                      SD_Dis = numeric(),
                      Normality_P_Value_Cont = numeric(),
                      Normality_P_Value_Dis = numeric(),
                      Wilcoxon_P_Value = numeric(),
                      stringsAsFactors = FALSE)

# calculate the INFO group
for (component in components) {
  # remove NA and filter
  cont_group <- CADD_com_analysis[CADD_com_analysis[[group_col]] == "Cont", component, drop = FALSE]
  dis_group <- CADD_com_analysis[CADD_com_analysis[[group_col]] == "Dis", component, drop = FALSE]
  
  # Check if both groups have enough data points (more than 1)
  if (length(cont_group[[component]]) > 1 & length(dis_group[[component]]) > 1) {
    # mean and sd calculation
    mean_cont <- mean(cont_group[[component]], na.rm = TRUE)
    mean_dis <- mean(dis_group[[component]], na.rm = TRUE)
    sd_cont <- sd(cont_group[[component]], na.rm = TRUE)
    sd_dis <- sd(dis_group[[component]], na.rm = TRUE)
    
    # Store results
    results <- rbind(results, data.frame(Component = component, 
                                         Mean_Cont = mean_cont, 
                                         Mean_Dis = mean_dis, 
                                         SD_Cont = sd_cont, 
                                         SD_Dis = sd_dis))
  } else {
    # Not enough data for mean or sd, set to NA
    results <- rbind(results, data.frame(Component = component, 
                                         Mean_Cont = NA, 
                                         Mean_Dis = NA, 
                                         SD_Cont = NA, 
                                         SD_Dis = NA))
  }
}

# normality
for (component in components) {
  cont_group <- CADD_com_analysis[CADD_com_analysis[[group_col]] == "Cont", component, drop = FALSE]
  dis_group <- CADD_com_analysis[CADD_com_analysis[[group_col]] == "Dis", component, drop = FALSE]
  
  # variability
  if (length(unique(cont_group[[component]])) < 2) {
    results[component, "Normality_P_Value_Cont"] <- NA  # 设置为 NA
    message(paste("Skipping normality test for Cont group - component:", component, "because all values are identical"))
  } else {
    # normality check
    normality_cont <- shapiro.test(cont_group[[component]])
    results[component, "Normality_P_Value_Cont"] <- normality_cont$p.value
  }
  
  # variability
  if (length(unique(dis_group[[component]])) < 2) {
    results[component, "Normality_P_Value_Dis"] <- NA  # 设置为 NA，表示没有进行检验
    message(paste("Skipping normality test for Dis group - component:", component, "because all values are identical"))
  } else {
    # normality check
    normality_dis <- shapiro.test(dis_group[[component]])
    results[component, "Normality_P_Value_Dis"] <- normality_dis$p.value
  }
}

# Wilcoxon test
for (i in seq_along(components)) {
  component <- components[i]
  
  cont_group <- CADD_com_analysis[CADD_com_analysis[[group_col]] == "Cont", component, drop = FALSE]
  dis_group <- CADD_com_analysis[CADD_com_analysis[[group_col]] == "Dis", component, drop = FALSE]
  
  if (length(cont_group[[component]]) > 1 & length(dis_group[[component]]) > 1) {
    # Wilcoxon test
    wilcox_test <- wilcox.test(cont_group[[component]], dis_group[[component]], exact = FALSE)
    results[i, "Wilcoxon_P_Value"] <- wilcox_test$p.value
  } else {
    # Not enough data for Wilcoxon, set to NA
    results[i, "Wilcoxon_P_Value"] <- NA
  }
}

# Final results
head(results)


------------###second try for anotation analysis and boxplot------------

# annotation selection
annotations <- c(
  "ConsScore", "GC", "CpG", "minDistTSS", "minDistTSE", "priPhCons", "mamPhCons", "verPhCons",
  "priPhyloP", "mamPhyloP", "verPhyloP", "bStatistic", "mirSVR-Score", "mirSVR-E", "mirSVR-Aln",
  "GerpRS", "GerpRSpval", "GerpN", "GerpS", "EncodeH3K4me1-sum", "EncodeH3K4me1-max",
  "EncodeH3K4me2-sum", "EncodeH3K4me2-max", "EncodeH3K4me3-sum", "EncodeH3K4me3-max",
  "EncodeH3K9ac-sum", "EncodeH3K9ac-max", "EncodeH3K9me3-sum", "EncodeH3K9me3-max",
  "EncodeH3K27ac-sum", "EncodeH3K27ac-max", "EncodeH3K27me3-sum", "EncodeH3K27me3-max",
  "EncodeH3K36me3-sum", "EncodeH3K36me3-max", "EncodeH3K79me2-sum", "EncodeH3K79me2-max",
  "EncodeH4K20me1-sum", "EncodeH4K20me1-max", "EncodeH2AFZ-sum", "EncodeH2AFZ-max",
  "EncodeDNase-sum", "EncodeDNase-max", "EncodetotalRNA-sum", "EncodetotalRNA-max",
  "Dist2Mutation", "Freq100bp", "Rare100bp", "Sngl100bp", "Freq1000bp", "Rare1000bp",
  "Sngl1000bp", "Freq10000bp", "Rare10000bp", "Sngl10000bp", "RemapOverlapTF", "RemapOverlapCL",
  "RegSeq0", "RegSeq1", "RegSeq2", "RegSeq3", "RegSeq4", "RegSeq5", "RegSeq6", "RegSeq7",
  "ZooPriPhyloP", "ZooVerPhyloP", "RawScore", "PHRED"
)

# create a new df
results <- lapply(annotations, function(col) {
  x <- CADD_elements[[col]]
  group <- CADD_elements$INFO
  
  
  # group
  x_cont <- x[group == "Cont"]
  x_dis  <- x[group == "Dis"]
  
  # culculate 
  mean_score <- mean(x, na.rm = TRUE)
  mean_cont <- mean(x_cont, na.rm = TRUE)
  mean_dis  <- mean(x_dis, na.rm = TRUE)
  sd_cont   <- sd(x_cont, na.rm = TRUE)
  sd_dis    <- sd(x_dis, na.rm = TRUE)
  num_cont  <- sum(!is.na(x_cont))
  num_dis   <- sum(!is.na(x_dis))
  
  # Wilcoxon test
  p_wilcox <- tryCatch({
    wilcox.test(x_dis, x_cont)$p.value
  }, error = function(e) NA)
  
  # normality
  p_norm <- tryCatch({
    if (length(na.omit(x)) >= 3) shapiro.test(x)$p.value else NA
  }, error = function(e) NA)
  
  # rebind the df
  data.frame(
    annotation = col,
    mean_score = mean_score,
    mean_cont = mean_cont,
    mean_dis = mean_dis,
    sd_cont = sd_cont,
    sd_dis = sd_dis,
    wilcoxon_p_value = p_wilcox,
    normal_distribution_p = p_norm,
    num_cont = num_cont,
    num_dis = num_dis
  )
})


# output
annotations <- do.call(rbind, results)


head(annotations)

#write.csv(annotations, "CADD_annotation_summary.csv", row.names = FALSE)

str(annotations)

# select the top different
annotations <- annotations %>%
  mutate(abs_diff = abs(mean_dis - mean_cont)) %>%
  arrange(wilcoxon_p_value, desc(abs_diff))

top_annotations <- annotations %>%
  slice_head(n = 6) %>%
  pull(annotation)

cat("Top annotations selected:\n")
print(top_annotations)

# get selected data
top_annotations_data <- CADD_elements %>%
  select(INFO, all_of(top_annotations)) %>%
  pivot_longer(
    cols = -INFO, 
    names_to = "annotation", 
    values_to = "score"  
  )


# boxplot

ggplot(top_annotations_data, aes(x = INFO, y = score, fill = INFO)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 1, color = "black") +
  facet_wrap(~annotation, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("Cont" = "lightblue", "Dis" = "salmon")) +
  labs(
    title = "Top CADD annotation features by group",
    x = "Group",
    y = "Annotation score"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey80", fill = NA)
  )

##---------------每个图单独输出保存------

library(stringr)

# FDR 校正
annotations <- annotations %>%
  mutate(FDR = p.adjust(wilcoxon_p_value, method = "BH")) %>%
  arrange(FDR, wilcoxon_p_value, desc(abs_diff))

# 显著项（FDR < 0.05）
sig_annos <- annotations %>% filter(!is.na(FDR), FDR < 0.05) %>% pull(annotation)
message("Significant components (FDR<0.05): ", length(sig_annos))

# 输出目录
dir_sig <- "figures/components_sig"
dir.create(dir_sig, showWarnings = FALSE, recursive = TRUE)

# 单图绘制并保存（boxplot）
plot_one <- function(feature, df = CADD_elements, outdir = dir_sig) {
  d <- df |> select(INFO, !!sym(feature)) |> rename(score = !!sym(feature)) |> tidyr::drop_na(score)
  if (length(unique(d$INFO)) < 2) return(invisible(NULL))
  
  pval <- tryCatch(wilcox.test(score ~ INFO, data = d, exact = FALSE)$p.value, error = function(e) NA_real_)
  subt <- ifelse(is.na(pval), "Wilcoxon p = NA", paste0("Wilcoxon p = ", signif(pval, 3)))
  
  p <- ggplot(d, aes(INFO, score, fill = INFO)) +
    geom_boxplot(outlier.shape = 21, outlier.size = 1, color = "black") +
    scale_fill_manual(values = c(Cont="lightblue", Dis="salmon")) +
    labs(title = paste0("Component: ", feature), subtitle = subt, x = "Group", y = "Score") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none", panel.grid.minor = element_blank(),
          panel.border = element_rect(color="grey80", fill=NA))
  
  fn <- paste0(str_replace_all(feature, "[^A-Za-z0-9_]+", "_"), ".png")
  ggsave(file.path(outdir, fn), p, width = 6, height = 4, dpi = 300)
}

invisible(lapply(sig_annos, plot_one))


#------ annotation里面的图-------