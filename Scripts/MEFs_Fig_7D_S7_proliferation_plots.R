#!/storage/Software/packages/anaconda3/bin/Rscript
rm(list=ls())
gc()


##############################################################################################################
#                                                                                                            #
#   Project: MBU_spb54_005_MouseEmbryo                                                                       #
#   Malwina Prater (mn367@cam.ac.uk), 2022                                                                   #
#   MRC MBU, University of Cambridge                                                                         #
#   Script: bulk RNA-seq in mouse MEFs clones - proliferation plots                                          # 
#                                                                                                            #
##############################################################################################################


message("+--- Loading in the libraries (start up messages supressed) ---+")
suppressPackageStartupMessages({
  library("ggplot2")
  library("ggrepel")
  library("cowplot")
  library("dplyr")
  theme_set(theme_cowplot())
  library(gridExtra)
  library(ggchromatic) 
  library(ggpubr)
})
  
  



remove_outliers <- TRUE
experiment_no <- "expt1_2" #  "expt1" "expt2" "expt1_2"


Project        <- "MBU_spb54_005_MEF_5024_siRNA_plots_"
baseDir        <- "/Users/xxx/Documents/xxx/xxx/xxx" # replace with your path
setwd(baseDir)


Heteroplasmy_meta <- data.frame(clone=c("clone_62","clone_59","clone_44","clone_83","clone_26","clone_114",
                                        "clone_101","clone_109","clone_17","clone_12","clone_33","clone_48"), 
                                Heteroplasmy=c(88,89,91,90,89,88,   19,31,27,34,23,30))



if(experiment_no=="expt1"){
  # expt 1 :::
  siRNA_df <- read.csv("Input/siRNA_IncuCyte_Raw_Data.csv")
  colnames(siRNA_df)[colnames(siRNA_df) == "Clone.62.Mock"] <- "Mean_Confluency"
} else if(experiment_no=="expt2"){
  # expt 2 :::
  siRNA_df <- read.csv("Input/siRNA_IncuCyte_Raw_Data_Expt2.csv")
  siRNA_df$Clone <- paste0("clone_", siRNA_df$Clone)
  colnames(siRNA_df)[colnames(siRNA_df) == "Mean.Confluency"] <- "Mean_Confluency"
} else if(experiment_no=="expt1_2"){
  siRNA_df1 <- read.csv("Input/siRNA_IncuCyte_Raw_Data.csv")
  siRNA_df2 <- read.csv("Input/siRNA_IncuCyte_Raw_Data_Expt2.csv")
  siRNA_df2$Clone <- paste0("clone_", siRNA_df2$Clone)
  siRNA_df1$expt <- 1
  siRNA_df2$expt <- 2
  colnames(siRNA_df1)[colnames(siRNA_df1) == "Clone.62.Mock"] <- "Mean_Confluency"
  colnames(siRNA_df2)[colnames(siRNA_df2) == "Mean.Confluency"] <- "Mean_Confluency"
  siRNA_df <- rbind(siRNA_df1, siRNA_df2)
} else{ print("please prive expriment no!!!")}


if(remove_outliers==TRUE){
  siRNA_df <- siRNA_df[siRNA_df$Clone != "clone_17",]
  siRNA_df <- siRNA_df[siRNA_df$Clone != "clone_114",]
}




siRNA_df$Heteroplasmy <- Heteroplasmy_meta[match(siRNA_df$Clone , Heteroplasmy_meta$clone),]$Heteroplasmy
head(siRNA_df)

siRNA_df$Mean_Confluency <- as.numeric(as.character(siRNA_df$Mean_Confluency))
siRNA_df$SEM <- as.numeric(as.character(siRNA_df$SEM))




message("+-------------------------------------------------------------------------------")
message("+                     plot siRNA line plots                                     ")
message("+-------------------------------------------------------------------------------")

stderror <- function(x) sd(x)/sqrt(length(x))

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      se = sd(x[[col]], na.rm=TRUE) / sqrt(sum(!is.na(x[[col]]))))  
  }
  
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}




plot_per_sirna <- function(siRNA_target){
  
  df <- siRNA_df[siRNA_df$siRNA %in% c(siRNA_target, "Mock", "NT"), c("Clone","Heteroplasmy","Heteroplasmy_bin", "siRNA","well","Elapsed","Mean_Confluency", "SEM"  )]
  
  df$group <- paste(df$Heteroplasmy_bin , df$siRNA)
  if(experiment_no == "expt1_2"){
    df2 <- data_summary(data=df, varname="Mean_Confluency", groupnames=c("Heteroplasmy_bin","siRNA","Elapsed", "Clone"))
    df_summarised <- data_summary(data=df2, varname="Mean_Confluency", groupnames=c("Heteroplasmy_bin","siRNA","Elapsed"))
  } else{
    df_summarised <- data_summary(data=df, varname="Mean_Confluency", groupnames=c("Heteroplasmy_bin","siRNA","Elapsed"))
  }
  df_summarised$group <- paste(df_summarised$Heteroplasmy_bin , df_summarised$siRNA)
  df_summarised$group <- factor(  df_summarised$group, levels= c("low_het Mock", "high_het Mock", "low_het NT", "high_het NT", paste0("low_het ",siRNA_target), paste0("high_het ",siRNA_target)))
  df_summarised$se[is.na(df_summarised$se)] <- 0
  df_summarised$sd[is.na(df_summarised$sd)] <- 0
  
  plt2 <- ggplot(data=df_summarised, aes(x=Elapsed, y=Mean_Confluency, group=group, color=group, shape= Heteroplasmy_bin)) +
    geom_line()+ geom_errorbar(aes(ymin=Mean_Confluency-se, ymax=Mean_Confluency+se), width=.1) + geom_point() + ggtitle(siRNA_target)+ labs(x="Time (h)")   + scale_color_manual(values=c("thistle",'#999999', "seagreen3","green3", "indianred", "red"))  
  plt2
  
  return(plt2)
}

unique(siRNA_df$siRNA)

plt_Gapdh <- plot_per_sirna("Gapdh")
plt_E2f3 <- plot_per_sirna("E2f3")
plt_Taf1 <- plot_per_sirna("Taf1")
plt_Klf12 <- plot_per_sirna("Klf12")
plt_Cnot3 <- plot_per_sirna("Cnot3")
plt_Maz <- plot_per_sirna("Maz")
plt_Bclaf1 <- plot_per_sirna("Bclaf1")

common_legend <- get_legend(plt_Gapdh)


pdf(paste(Project, "siRNA",  "lineplots", clones_to_plot, experiment_no, "RMoutliers" ,".pdf", sep="_"), onefile=FALSE, width=12, height=6)
par(bg=NA)
ggarrange(plt_Gapdh,  plt_E2f3, plt_Taf1,
          plt_Klf12, plt_Cnot3, plt_Maz,
          plt_Bclaf1, ncol=4, nrow=2,  align="hv" , common.legend = TRUE,legend="right" )# ,  legend.grob = common_legend
dev.off()




message("+-------------------------------------------------------------------------------")
message("+                     stats for line plots                                      ")
message("+-------------------------------------------------------------------------------")

timepoint <- 60


df <- siRNA_df[siRNA_df$Elapsed ==timepoint, c("Clone","Heteroplasmy", "siRNA","Mean_Confluency", "Heteroplasmy_bin"  )]
df$Mean_Confluency <- as.numeric(as.character(df$Mean_Confluency))


# continuous correlation Mean_Confluency ~ Heteroplasmy 

for(siRNA in unique(df$siRNA)){
  df_sub <- df[df$siRNA == siRNA,]
  print(c("Spearman and Kendall correlation for ", siRNA))
  cor_res <- cor.test(log2(df_sub$Heteroplasmy), df_sub$Mean_Confluency, method = "spearman") # SIG
  print(cor_res$p.value)
  
  # mixed linear model 
  print("mixed linear model:: Time_s ~ log2(Burden)")
  obj.lme=nlme::lme(Mean_Confluency~log2(Heteroplasmy), data= df_sub, random = ~ 1|Clone)
  summary(obj.lme)
  print(anova(obj.lme))
}




# t-test  Mean_Confluency in high vs low Heteroplasmy 

t_test_list <- list()
count <- 1
for(siRNA in unique(df$siRNA)){
  df_sub <- df[df$siRNA == siRNA,]
  print(c("ANOVA for ", siRNA))
  ttest_res <- t.test(Mean_Confluency ~ Heteroplasmy_bin, data=df_sub)
  t_test_list[[count]] <- ttest_res$p.value
  count <- count +1
}
names(t_test_list) <-unique(df$siRNA)
t_test_df <- data.frame(t_test_pval=unlist(t_test_list))

t_test_df$siRNA <- rownames(t_test_df)
t_test_df <- t_test_df[t_test_df$siRNA != "Cnot3",]

t_test_df$padj_FDR <- p.adjust(t_test_df$t_test_pval, method = "fdr")
t_test_df$padj_bonferroni <- p.adjust(t_test_df$t_test_pval, method = "bonferroni")

#write.csv(t_test_df, paste0(Project,"_", timepoint,"_t_test.csv"))





# anova on AUC score of Mean_Confluency in high vs low Heteroplasmy

library(bayestestR)
area_under_curve(x=siRNA_df$Mean_Confluency, y=siRNA_df$Elapsed, method = "trapezoid")


# AUC for individual clones

AUC_list <- list()
Clone_list <- list()
count <- 1
for(siRNA in unique(siRNA_df$siRNA)){
  df_sub <- siRNA_df[siRNA_df$siRNA == siRNA,]
  count2 <- 1
  for(Clone in unique(siRNA_df$Clone)){
    Clone_list[[count2]] <- area_under_curve(x=df_sub[df_sub$Clone == Clone,]$Mean_Confluency, y=df_sub[df_sub$Clone == Clone,]$Elapsed, method = "trapezoid")
    count2 <- count2 + 1
  }
  names(Clone_list) <- unique(siRNA_df$Clone)
  AUC_list[[count]] <- Clone_list
  count <- count + 1
}
names(AUC_list) <- unique(siRNA_df$siRNA)
AUC_list

AUC_df <- data.frame(AUC=unlist(AUC_list))
AUC_df$clone <- gsub( ".*clone", "clone", rownames(AUC_df))
AUC_df$siRNA <- gsub( ".clone.*", "", rownames(AUC_df))
AUC_df$Heteroplasmy <- Heteroplasmy_meta[match( AUC_df$clone, Heteroplasmy_meta$clone),]$Heteroplasmy
AUC_df$Heteroplasmy_bin <- ifelse(AUC_df$Heteroplasmy < 50, "LowHet", "HighHet")


t_test_res_list <- list()
count <-1
for(siRNA in unique(AUC_df$siRNA)){
  print(siRNA)
  print(test_ANOVA(AUC_df[AUC_df$siRNA == siRNA,], "AUC", groups = "Heteroplasmy_bin"))
  t_test_res <- test_ANOVA(AUC_df[AUC_df$siRNA == siRNA,], "AUC", groups = "Heteroplasmy_bin")
  t_test_res_list[[count]] <- t_test_res$ANOVA_pval
  count <- count + 1
}
names(t_test_res_list) <- unique(AUC_df$siRNA)

AUC_t_test <- data.frame(pval= unlist(t_test_res_list))
AUC_t_test$siRNA <- rownames(AUC_t_test)
AUC_t_test <- AUC_t_test[AUC_t_test$siRNA != "Cnot3",]

AUC_t_test$padj_FDR <- p.adjust(AUC_t_test$pval, method = "fdr")
AUC_t_test$padj_bonferroni <- p.adjust(AUC_t_test$pval, method = "bonferroni")

#write.csv(AUC_t_test, paste0(Project, "_AUC_t_test.csv"))

















