# This script is to draw ROC curve with number of combination with/without physical interaction 
# version considering number of cell line without RBP expression
# 2023/06/27 made

# make new directory
setwd("C:/Rdata")
dir.create("20230627_ROC_curve_of_pvalue_after_cutoff_using_TPM_data")

# activate package to draw ROC curve
library(Epi)

# import result of correlation analysis between expression level of 410 RBPs and 377 miRNAs using TPM data
# this result is located at "https://github.com/Ryosuke-Hirota/20230517_CCLE_correlation_between_residual_and_RBPDB_RBP_RSEM_TPM"
setwd("C:/Rdata/20230517_CCLE_correlation_between_residual_and_RBPDB_RBP_RSEM_TPM")
c.result <-read.table("summary_of_correlation_between_residual_and_RBPDB_RBP_RSEM_TPM.txt",sep="\t",header = T,stringsAsFactors = F)

# edit list and convert p.value (because I couldn't draw correct ROC curve without converting p.value)
c.result <-subset(c.result,!is.na(c.result[,5]))

# extract the name of primary miRNA, primary transcript and mature miRNA 
pri.transcript <-c.result[,c(1:3)]
pri.transcript <-pri.transcript[!duplicated(pri.transcript),]
colnames(pri.transcript)[1:2] <-c("primary","mature")

c.result[,1] <-paste0(c.result[,1],"_",c.result[,2],"_",c.result[,3],"_vs_",c.result[,4])
c.result[,6] <-log10(c.result[,6])*-1
c.result <-c.result[,c(-2,-3,-4)]
colnames(c.result)[1] <-"combination"

# Import list of physical interaction
# This list is located at "https://github.com/Ryosuke-Hirota/20221122_ROC_curve_of_pvalue_after_cutoff"
# Attention : this list is made by pri-miRNA. Thus, when you correspond this list to the correlation analysis using mature miRNA, duplicated rows appear. 
setwd("C:/Rdata/20221017_ROC_curve_for_cutoff_with_functional_interactions")
phy.list <-read.table("list_of_treiber_physical_interaction_between_RBP_and_miRNA.txt",sep="\t",header = T,stringsAsFactors = F)

# edit list of physical interaction 
phy.list[,1] <-paste0(phy.list[,1],"_",phy.list[,2],"_",phy.list[,3],"_vs_",phy.list[,4])
colnames(phy.list)[1] <-"combination"
phy.list <-phy.list[,c(-2,-3,-4)]

# merge correlation analysis result and list of physical interaction
merge.df <-merge(c.result,phy.list,by="combination",all=T)
merge.df <-subset(merge.df,!is.na(merge.df[,2]))

# annotate whether treiber's combinations match with CCLE correlation analysis result
merge.df[!is.na(merge.df[,5]),6] <-"match"
merge.df[is.na(merge.df[,5]),6] <-"not_match"
colnames(merge.df)[6] <-"match_with_CCLE"
merge.df[,5] <-ifelse(is.na(merge.df[,5]),0,merge.df[,5])

# check number of combinations with/without physical interaction
# 354(unique:431) treiber's combinations with physical interaction(score>3) are matched with CCLE correlation analysis result
# 4098(unique:4420) treiber's combinations without physical interaction(score<=3) are matched with CCLE correlation analysis result
phy <-merge.df[merge.df[,5]>=3&merge.df[,6]=="match",1]
length(unique(phy))

no.phy <-merge.df[merge.df[,5]<3&merge.df[,6]=='match',1]
length(unique(no.phy))

# set value of cutoff
cutoff <-seq(50,900,50)

setwd("C:/Rdata/20230627_ROC_curve_of_pvalue_after_cutoff_using_TPM_data")

# draw ROC curve including all combinations
for (i in 1:length(cutoff)) {
  
  result <-merge.df
  
  # annotate positive and negative by physical interaction
  result[result[,5]>=3,7] <-1
  result[result[,5]<3,7] <-0
  colnames(result)[7] <-"physical_interaction"
  
  # cutoff number of cell line
  result <-result[result[,4]>=cutoff[i],]
  
  # count number of combinations with/without physical interaction
  # caution : this number is duplicated. If you wanna know non-duplicated number, write additional script.
  count <-as.data.frame(table(result[,7]))
  p <-count[2,2]
  np <-count[1,2]
  
  # draw ROC curve
  pdf(paste0("ROC_curve_pvalue_cutoff_",cutoff[i],"_including_all_combinations.pdf"))
  ROC(test=result$p.value, stat=result$physical_interaction, plot="ROC")
  mtext(text = paste0("physical interaction : ",p," , no physical interaction : ",np))
  dev.off()
}

# draw ROC curve not including all combinations
for (i in 1:length(cutoff)) {
  
  result <-merge.df
  result <-result[result[,6]=="match",]
  
  # annotate positive and negative by physical interaction
  result[result[,5]>=3,7] <-1
  result[result[,5]<3,7] <-0
  colnames(result)[7] <-"physical_interaction"
  
  # cutoff number of cell line
  result <-result[result[,4]>=cutoff[i],]
  
  # count number of combinations with/without physical interaction
  # caution : this number is duplicated. If you wanna know non-duplicated number, write additional script.
  count <-as.data.frame(table(result[,7]))
  p <-count[2,2]
  np <-count[1,2]
  
  # draw ROC curve
  pdf(paste0("ROC_curve_pvalue_cutoff_",cutoff[i],"_not_including_all_combinations.pdf"))
  ROC(test=result$p.value, stat=result$physical_interaction, plot="ROC")
  mtext(text = paste0("physical interaction : ",p," , no physical interaction : ",np))
  dev.off()
}
