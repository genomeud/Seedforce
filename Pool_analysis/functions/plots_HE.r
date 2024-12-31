
################
# # ALL SNPs # #
################

# Boxplot HE

pools <- c('ACI','ALI','AVE','BSA','CAE','CAM','CPU','CSA','ECA','GIL','GLI','GPA','HAD','KOS','LIS','LPS','MAQ','RSA','SHI','WOR')

for (spp in pools)
{
setwd(paste0("/projects/marroni/seedforce/pool/de_novo_pipeline_catalog_spp_pops/",spp,""))
ALL <- read.table(paste0("",spp,"_allfreq_HE_name.tbl"), header=T)
full <- as.data.frame(ALL[,c(11,12)])
colnames(full) <- NULL
colnames(full)=c('EXP_HET', 'NAME')
library("gplots")
fileOut=(paste0("boxplot_exp_het_pool_",spp,".pdf"))
pdf(fileOut, paper="A4")
par(cex.axis = 1, cex.lab = 1.2)
par(cex.main = 1.6)
par(plt = c(0.1, 0.5, 0.1, 0.5), mar = c(8, 6, 6, 4) + 0.1)
colors <- c("seagreen","tomato2","hotpink4","steelblue","lightgoldenrod","grey","cyan3","orange3","snow2")
boxplot2(full$EXP_HET ~ full$NAME, 
    xlab = "Category", 
    ylab = "Expected heterozygosity",
    col = colors,
    range = 0.1,
    main = paste0("Expected heterozygosity Pool ",spp,""),frame = FALSE, top = TRUE)
dev.off() 
}

# mean exp heterozygosity with error bars. first plot with standard deviation and second plot with standard error
# remove standard deviation in negative part
library(tidyverse)

#For all the SNPs belonging to one group compute average heterozygosity (N_HET_EXP and P_HET_EXP). 

pools <- c('ACI','ALI','AVE','BSA','CAE','CAM','CPU','CSA','ECA','GIL','GLI','GPA','HAD','KOS','LIS','LPS','MAQ','RSA','SHI','WOR')

for (spp in pools)
{
setwd(paste0("/projects/marroni/seedforce/pool/de_novo_pipeline_catalog_spp_pops/",spp,""))
ALL <- read.table(paste0("",spp,"_allfreq_HE_name.tbl"), header=T)
full <- as.data.frame(ALL[,c(11,12)])
# calculate mean expected heterozygosity
mean_P_HET_EXP <- (ALL %>% 
   select(POP_NAME, P_HET_EXP) %>% 
   group_by(POP_NAME) %>% 
   summarise(P_HET_EXP_mean = mean(P_HET_EXP)))
# calculate standard deviation
sd <- aggregate(ALL$P_HET_EXP, list(ALL$POP_NAME), FUN=sd) 
colnames(sd) <- NULL
colnames(sd)=c('POP_NAME', 'SD')
# combine both tables
ALL <- cbind(mean_P_HET_EXP,sd$SD)
colnames(ALL) <- NULL
colnames(ALL)=c('POP_NAME','mean_P_HET_EXP','SD')
out_color=paste0("mean_exp_het_",spp,"_sd.jpeg")
jpeg(out_color,width=15,height=16,units="cm",res=300,type="cairo")
# fileOut=(paste0("mean_exp_het_",spp,".pdf"))
# pdf(fileOut, paper="A4")
library(ggplot2)
# function for plot size
# width <- "7"
# height <- "7"
# fig <- function(width, heigth){
     # options(repr.plot.width = width, repr.plot.height = heigth)
# }
# Default bar plot
print(ggplot(ALL, aes(x=POP_NAME, y=mean_P_HET_EXP, fill=POP_NAME)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean_P_HET_EXP, ymax=mean_P_HET_EXP+SD), width=.2,
                 position=position_dodge(.9)) +
labs(title=paste0("Expected heterozygosity ",spp,""), x="", y = "Mean expected heterozygosity", size=22)+
   theme_classic()+
scale_y_continuous(limit = c(0, 0.5)) +
labs(fill = "Populations", size=16) +
theme(plot.title = element_text(hjust = 0.5, size=17), axis.text.x=element_text(angle=90, hjust=1, size=16),axis.text=element_text(size=14), axis.text.y=element_text(hjust=1, size=17), panel.spacing.x=unit(0.5, "lines")))
dev.off()
# calculate standard error
SE <- (HE %>%
  group_by(POP_NAME) %>%
  summarize(SE=plotrix::std.error(P_HET_EXP)))
colnames(sd) <- NULL
colnames(sd)=c('POP_NAME', 'SE')
# combine both tables
ALL <- cbind(mean_P_HET_EXP,SE$SE)
colnames(ALL) <- NULL
colnames(ALL)=c('POP_NAME','mean_P_HET_EXP','SE')
ALL$spp=spp
write.table(ALL,"exp_heterozygosity_and_se_per_pop.tbl",sep="\t", quote = FALSE, row.names=FALSE)


# fileOut=(paste0("mean_exp_het_",spp,"_SErr.pdf"))
# pdf(fileOut, paper="A4")
library(ggplot2)
# Default bar plot
out_color=paste0("mean_exp_het_",spp,"_se.jpeg")
jpeg(out_color,width=15,height=16,units="cm",res=300,type="cairo")
print(ggplot(ALL, aes(x=POP_NAME, y=mean_P_HET_EXP, fill=POP_NAME)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean_P_HET_EXP-SE, ymax=mean_P_HET_EXP+SE), width=.2,
                 position=position_dodge(.9)) +
labs(title=paste0("Expected heterozygosity ",spp,""), x="", y = "Mean expected heterozygosity", size=22)+
   theme_classic()+
scale_y_continuous(limit = c(0, 0.5)) +
labs(fill = "Populations", size=16) +
theme(plot.title = element_text(hjust = 0.5, size=17), axis.text.x=element_text(angle=90, hjust=1, size=16),axis.text=element_text(size=14), axis.text.y=element_text(hjust=1, size=17), panel.spacing.x=unit(0.8, "lines")))
dev.off()
}


# wilcoxon test and t test

# WILCOXON TEST
pp<-pairwise.wilcox.test(full[,1], full[,2], p.adjust.method = "none", paired = FALSE)
mymat<-tri.to.squ(pp$p.value)
myletters<-multcompLetters(mymat,compare="<=",threshold=0.05,Letters=letters)
# Convert myletters to a dataframe so it can be used in ggplot
letters_df <- data.frame(POP_NAME = names(myletters$Letters), Letters = myletters$Letters)

ALL <- merge(ALL, letters_df, by = "POP_NAME")

library(ggplot2)
# Default bar plot
out_color=paste0("mean_exp_het_",spp,"_se_wilcx.jpeg")
jpeg(out_color,width=15,height=15,units="cm",res=300,type="cairo")
print(ggplot(ALL, aes(x=POP_NAME, y=mean_P_HET_EXP, fill=POP_NAME)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean_P_HET_EXP-SE, ymax=mean_P_HET_EXP+SE), width=.2,
                 position=position_dodge(.9)) +
labs(title=paste0("Expected heterozygosity ",spp,""), x="", y = "Mean expected heterozygosity", size=18)+
   theme_classic()+
scale_y_continuous(limit = c(0, 0.3)) +
labs(fill = "Populations", size=15) +
theme(plot.title = element_text(hjust = 0.5, size=17), axis.text.x=element_text(angle=90, hjust=1, size=15),axis.text=element_text(size=14), axis.text.y=element_text(hjust=1, size=15), panel.spacing.x=unit(0.8, "lines")) +
geom_text(aes(x = factor(POP_NAME), y = max(mean_P_HET_EXP) + 0.12, label = Letters), 
              color = "black", size = 6))
dev.off()


# T TEST PAIRWISE COMPARISON
pp <- pairwise.t.test(full[,1], full[,2], p.adjust.method = "none", paired = FALSE)
mymat <- tri.to.squ(pp$p.value)
myletters <- multcompLetters(mymat, compare = "<=", threshold = 0.05, Letters = letters)
letters_pw <- data.frame(POP_NAME = names(myletters$Letters), Letters = myletters$Letters)
# calculate standard error
SE <- (HE %>%
  group_by(POP_NAME) %>%
  summarize(SE=plotrix::std.error(P_HET_EXP)))
colnames(sd) <- NULL
colnames(sd)=c('POP_NAME', 'SE')
# combine both tables
ALL <- cbind(mean_P_HET_EXP,SE$SE)
colnames(ALL) <- NULL
colnames(ALL)=c('POP_NAME','mean_P_HET_EXP','SE')
ALL <- merge(ALL, letters_pw, by = "POP_NAME")
out_color=paste0("mean_exp_het_",spp,"_se_ttest.jpeg")
jpeg(out_color,width=16,height=15,units="cm",res=300,type="cairo")
print(ggplot(ALL, aes(x=POP_NAME, y=mean_P_HET_EXP, fill=POP_NAME)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean_P_HET_EXP-SE, ymax=mean_P_HET_EXP+SE), width=.2,
                 position=position_dodge(.9)) +
labs(title=paste0("Expected heterozygosity ",spp,""), x="", y = "Mean expected heterozygosity", size=18)+
   theme_classic()+
scale_y_continuous(limit = c(0, 0.3)) +
labs(fill = "Populations", size=15) +
theme(plot.title = element_text(hjust = 0.5, size=17), axis.text.x=element_text(angle=90, hjust=1, size=15),axis.text=element_text(size=14), axis.text.y=element_text(hjust=1, size=15), panel.spacing.x=unit(0.8, "lines")) +
geom_text(aes(x = factor(POP_NAME), y = max(mean_P_HET_EXP) + 0.12, label = Letters),color = "black", size = 6))
dev.off()



#####################
# # FILTERED SNPs # #
#####################

pools <- c('ACI','ALI','AVE','BSA','CAE','CAM','CPU','CSA','ECA','GIL','GLI','GPA','GPA2','HAD','KOS','LIS','LPS','MAQ','RSA','SHI','WOR')

for (spp in pools)
{
    for (reads in c(5))
    {
    setwd(paste0("/projects/marroni/seedforce/pool/de_novo_pipeline_catalog_spp_pops/",spp,""))
    ALL <- read.table(paste0("",spp,"_allfreq_HE_name_minorALL",reads,"reads.tbl"), header=T)
    full <- as.data.frame(ALL[,c(11,12)])
    # calculate mean expected heterozygosity
    mean_P_HET_EXP <- (ALL %>% 
       select(POP_NAME, P_HET_EXP) %>% 
       group_by(POP_NAME) %>% 
       summarise(P_HET_EXP_mean = mean(P_HET_EXP)))
    # calculate standard deviation
    sd <- aggregate(ALL$P_HET_EXP, list(ALL$POP_NAME), FUN=sd) 
    colnames(sd) <- NULL
    colnames(sd)=c('POP_NAME', 'SD')
    # combine both tables
    ALL <- cbind(mean_P_HET_EXP,sd$SD)
    colnames(ALL) <- NULL
    colnames(ALL)=c('POP_NAME','mean_P_HET_EXP','SD')
    ALL$spp=spp
    write.table(ALL,paste0("exp_heterozygosity_per_pop_minorALL",reads,"reads.tbl"),sep="\t", quote = FALSE, row.names=FALSE)
    # outpath <- paste0("/projects/marroni/seedforce/pool/delivery_20231108-1426/de_novo_pipeline_catalog_spp_pops/",spp,"")

    # out=paste(outpath,"mean_exp_het_",spp,".jpeg",sep="")
    fileOut=(paste0("mean_exp_het_",spp,"_filtered_",reads,"_reads.pdf"))
    pdf(fileOut, paper="A4")
    # jpeg(out,width=16,height=8,units="cm",res=300,type="cairo")
    library(ggplot2)
    # Default bar plot
    print(ggplot(ALL, aes(x=POP_NAME, y=mean_P_HET_EXP, fill=POP_NAME)) + 
      geom_bar(stat="identity", color="black", 
               position=position_dodge()) +
      geom_errorbar(aes(ymin=mean_P_HET_EXP, ymax=mean_P_HET_EXP+SD), width=.2,
                     position=position_dodge(.9)) +
    labs(title=paste0("Expected heterozygosity ",spp,""), x="Population", y = "Mean expected heterozygosity", size=18)+
       theme_classic()+
    scale_y_continuous(limit = c(0, 0.5)) +
    labs(fill = "Populations", ) +
    theme(plot.title = element_text(hjust = 0.5, size=18), axis.text.x=element_text(angle=90, hjust=1, size=12), axis.text.y=element_text(hjust=1, size=12), panel.spacing.x=unit(0.5, "lines")))
    # calculate standard error
    SE <- (full %>%
      group_by(POP_NAME) %>%
      summarize(SE=plotrix::std.error(P_HET_EXP)))
    colnames(sd) <- NULL
    colnames(sd)=c('POP_NAME', 'SE')
    # combine both tables
    ALL <- cbind(mean_P_HET_EXP,SE$SE)
    colnames(ALL) <- NULL
    colnames(ALL)=c('POP_NAME','mean_P_HET_EXP','SE')
    ALL$spp=spp
    write.table(ALL,paste0("exp_heterozygosity_and_se_per_pop_minorALL",reads,"reads.tbl"),sep="\t", quote = FALSE, row.names=FALSE)
    library(ggplot2)
    # Default bar plot
    print(ggplot(ALL, aes(x=POP_NAME, y=mean_P_HET_EXP, fill=POP_NAME)) + 
      geom_bar(stat="identity", color="black", 
               position=position_dodge()) +
      geom_errorbar(aes(ymin=mean_P_HET_EXP-SE, ymax=mean_P_HET_EXP+SE), width=.2,
                     position=position_dodge(.9)) +
    labs(title=paste0("Expected heterozygosity ",spp,""), x="Population", y = "Mean expected heterozygosity", size=18)+
       theme_classic()+
    scale_y_continuous(limit = c(0, 0.5)) +
    labs(fill = "Populations", ) +
    theme(plot.title = element_text(hjust = 0.5, size=18), axis.text.x=element_text(angle=90, hjust=1, size=12), axis.text.y=element_text(hjust=1, size=12), panel.spacing.x=unit(0.5, "lines")))
    dev.off()
}
}
    # wilcoxon test and t test

# WILCOXON TEST
library(multcompView)
pools <- c('ACI','ALI','BSA','CAE','CAM','CPU','CSA','ECA','GLI','GPA2','HAD','KOS','LPS','MAQ','SHI','WOR')
# pools <- c('ACI','ALI','AVE','BSA','CAE','CAM','CPU','CSA','ECA','GIL','GLI','GPA','HAD','KOS','LIS','LPS','MAQ','RSA','SHI','WOR')

for (spp in pools)
{
    for (reads in c(5))
    {
        setwd(paste0("/projects/marroni/seedforce/pool/de_novo_pipeline_catalog_spp_pops/",spp,""))
    ALL <- read.table(paste0("",spp,"_allfreq_HE_name_minorALL",reads,"reads.tbl"), header=T)
    full <- as.data.frame(ALL[,c(11,12)])
    # calculate mean expected heterozygosity
    mean_P_HET_EXP <- (full %>% 
       select(POP_NAME, P_HET_EXP) %>% 
       group_by(POP_NAME) %>% 
       summarise(P_HET_EXP_mean = mean(P_HET_EXP)))
       SE <- (ALL %>%
      group_by(POP_NAME) %>%
      summarize(SE=plotrix::std.error(P_HET_EXP)))
    colnames(SE) <- NULL
    colnames(SE)=c('POP_NAME', 'SE')
    # combine both tables
    ALL <- cbind(mean_P_HET_EXP,SE$SE)
    colnames(ALL) <- NULL
    colnames(ALL)=c('POP_NAME','mean_P_HET_EXP','SE')
    ALL$spp=spp
tri.to.squ<-function(x)
{
rn<-row.names(x)
cn<-colnames(x)
an<-unique(c(cn,rn))
myval<-x[!is.na(x)]
mymat<-matrix(1,nrow=length(an),ncol=length(an),dimnames=list(an,an))
for(ext in 1:length(cn))
{
 for(int in 1:length(rn))
 {
 if(is.na(x[row.names(x)==rn[int],colnames(x)==cn[ext]])) next
 mymat[row.names(mymat)==rn[int],colnames(mymat)==cn[ext]]<-x[row.names(x)==rn[int],colnames(x)==cn[ext]]
 mymat[row.names(mymat)==cn[ext],colnames(mymat)==rn[int]]<-x[row.names(x)==rn[int],colnames(x)==cn[ext]]
 }
}
return(mymat)
}
pp<-pairwise.wilcox.test(full[,1], full[,2], p.adjust.method = "none", paired = FALSE)
mymat<-tri.to.squ(pp$p.value)
myletters<-multcompLetters(mymat,compare="<=",threshold=0.05,Letters=letters)
# Convert myletters to a dataframe so it can be used in ggplot
letters_df <- data.frame(POP_NAME = names(myletters$Letters), Letters = myletters$Letters)

ALL <- merge(ALL, letters_df, by = "POP_NAME")

library(ggplot2)
    fileOut=(paste0("mean_exp_het_",spp,"_filtered_",reads,"_se_wilcx.jpeg"))
jpeg(fileOut,width=15,height=15,units="cm",res=300,type="cairo")
print(ggplot(ALL, aes(x=POP_NAME, y=mean_P_HET_EXP, fill=POP_NAME)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean_P_HET_EXP-SE, ymax=mean_P_HET_EXP+SE), width=.2,
                 position=position_dodge(.9)) +
labs(title=paste0("Expected heterozygosity ",spp,""), x="", y = "Mean expected heterozygosity", size=18)+
   theme_classic()+
scale_y_continuous(limit = c(0, 0.5)) +
labs(fill = "Populations", size=15) +
theme(plot.title = element_text(hjust = 0.5, size=17), axis.text.x=element_text(angle=90, hjust=1, size=15),axis.text=element_text(size=14), axis.text.y=element_text(hjust=1, size=15), panel.spacing.x=unit(0.8, "lines")) +
geom_text(aes(x = factor(POP_NAME), y = max(mean_P_HET_EXP) + 0.12, label = Letters), color = "black", size = 6))
dev.off()
}
}



# T TEST PAIRWISE COMPARISON
library(multcompView)
library(ggplot2)
library(tidyverse)
pools <- c('ACI','ALI','BSA','CAE','CAM','CPU','CSA','ECA','GLI','GPA2','HAD','KOS','LPS','MAQ','SHI','WOR')
# 'AVE','GIL','GPA','LIS','RSA'
for (spp in pools)
{
    for (reads in c(5))
    {
        setwd(paste0("/projects/marroni/seedforce/pool/de_novo_pipeline_catalog_spp_pops/",spp,""))
    ALL <- read.table(paste0("",spp,"_allfreq_HE_name_minorALL",reads,"reads.tbl"), header=T)
    full <- as.data.frame(ALL[,c(11,12)])
    # calculate mean expected heterozygosity
    mean_P_HET_EXP <- (full %>% 
       select(POP_NAME, P_HET_EXP) %>% 
       group_by(POP_NAME) %>% 
       summarise(P_HET_EXP_mean = mean(P_HET_EXP)))
       SE <- (ALL %>%
      group_by(POP_NAME) %>%
      summarize(SE=plotrix::std.error(P_HET_EXP)))
    colnames(SE) <- NULL
    colnames(SE)=c('POP_NAME', 'SE')
    # combine both tables
    ALL <- cbind(mean_P_HET_EXP,SE$SE)
    colnames(ALL) <- NULL
    colnames(ALL)=c('POP_NAME','mean_P_HET_EXP','SE')
    ALL$spp=spp
tri.to.squ<-function(x)
{
rn<-row.names(x)
cn<-colnames(x)
an<-unique(c(cn,rn))
myval<-x[!is.na(x)]
mymat<-matrix(1,nrow=length(an),ncol=length(an),dimnames=list(an,an))
for(ext in 1:length(cn))
{
 for(int in 1:length(rn))
 {
 if(is.na(x[row.names(x)==rn[int],colnames(x)==cn[ext]])) next
 mymat[row.names(mymat)==rn[int],colnames(mymat)==cn[ext]]<-x[row.names(x)==rn[int],colnames(x)==cn[ext]]
 mymat[row.names(mymat)==cn[ext],colnames(mymat)==rn[int]]<-x[row.names(x)==rn[int],colnames(x)==cn[ext]]
 }
}
return(mymat)
}
pp <- pairwise.t.test(full[,1], full[,2], p.adjust.method = "none", paired = FALSE)
mymat <- tri.to.squ(pp$p.value)
myletters <- multcompLetters(mymat, compare = "<=", threshold = 0.05, Letters = letters)
letters_pw <- data.frame(POP_NAME = names(myletters$Letters), Letters = myletters$Letters)
# calculate standard error
SE <- (full %>%
  group_by(POP_NAME) %>%
  summarize(SE=plotrix::std.error(P_HET_EXP)))
colnames(SE) <- NULL
colnames(SE)=c('POP_NAME', 'SE')
# combine both tables
ALL <- cbind(mean_P_HET_EXP,SE$SE)
colnames(ALL) <- NULL
colnames(ALL)=c('POP_NAME','mean_P_HET_EXP','SE')
ALL <- merge(ALL, letters_pw, by = "POP_NAME")
out_color=paste0("mean_exp_het_",spp,"_",reads,"reads_se_ttest.jpeg")
jpeg(out_color,width=7,height=7,units="cm",res=300,type="cairo")
print(ggplot(ALL, aes(x=POP_NAME, y=mean_P_HET_EXP, fill=POP_NAME)) + 
  geom_bar(stat="identity", color="black", show.legend = FALSE,
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean_P_HET_EXP-SE, ymax=mean_P_HET_EXP+SE), width=.1,
                 position=position_dodge(.9)) +
labs(title=paste0("",spp,""), x="", y = "Expected Heterozygosity", size=8)+
   theme_classic()+
scale_y_continuous(limit = c(0, 0.5)) +
labs(fill = "Populations", size=6) +
theme(plot.title = element_text(hjust = 0.5, size=8), axis.text.x=element_text(angle=90, hjust=1, size=8),axis.text=element_text(size=8), axis.text.y=element_text(hjust=1, size=6), panel.spacing.x=unit(0.8, "lines")) +
geom_text(aes(x = factor(POP_NAME), y = 0.45, label = Letters),color = "black", size = 3))
dev.off()
    }
}



# ----- #

######################################
# Command for group SNPs by distance #
######################################

# this comand will create a column at the end with the SNP group were they belong (if SNPs are near to each other they will belong to the same group

setwd("/projects/marroni/pool_GBS/paloma/seedforce/05_SNPs_de_novo")
ind=read.table("08_mixed_01_ind_de_novo_2.txt", header=T)
pool=read.table("08_mixed_01_pool_de_novo_2.txt", header=T)

# REMOVE NAs
ind <- na.omit(ind)
pool <- na.omit(pool)

# create new column SNP_group with group number depending on position, near positions will have same SNP_group number
ind$SNP_group<-ind$CHROM
pool$SNP_group<-pool$CHROM

# create a table with last column including the SNP_group, one for ind and other for pool
write.table(ind, file = "08_mixed_01_ind_SNP_group_de_novo_2.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(pool, file = "08_mixed_01_pool_SNP_group_de_novo_2.txt", quote=FALSE, row.names=FALSE, sep="\t")

library(tidyverse)
head(ind)

#For all the SNPs belonging to one group compute average heterozygosity (N_HET_EXP and P_HET_EXP). 
# IND
mean_het_exp_N <- (ind %>% 
   select(SNP_group, N_HET_EXP) %>% 
   group_by(SNP_group) %>% 
   summarise(N_HET_EXP_mean = mean(N_HET_EXP)))
# POOL


mean_het_exp_P <- (pool %>% 
   select(SNP_group, P_HET_EXP) %>% 
   group_by(SNP_group) %>% 
   summarise(P_HET_EXP_mean = mean(P_HET_EXP)))

browser()

#If you manage to do that, plot the scatterplot of the new measures
# # SCATTERPLOT SAME AS BEFORE WITH NEW VALUES
#Regression line
z <- lm(mean_het_exp_N$N_HET_EXP_mean ~ mean_het_exp_P$P_HET_EXP_mean)
y <- lm(mean_het_exp_P$P_HET_EXP_mean ~ mean_het_exp_N$N_HET_EXP_mean)
# rsquare
rsquare<-round(summary(z)$adj.r.squared,3)
# calculate p-value
pvalue<-cor.test(mean_het_exp_P$P_HET_EXP_mean,mean_het_exp_N$N_HET_EXP_mean)$p.value
# create pdf output file
fileOut=("scatterplot_MEAN_exp_het_ind_vs_pool_reg_clean_de_novo_2.pdf")
# open pdf output file
pdf(fileOut, paper="A4")
# create the plot inside the opened pdf
# scatterplot mean heterozygosity POOL vs IND
plot(mean_het_exp_P$P_HET_EXP_mean,mean_het_exp_N$N_HET_EXP_mean,xlab="Mean Expected Heterozygosity POOL",ylab="Mean Expected Heterozygosity IND",pch=19,las=1,cex=0.3,main="Mean Expected Heterozygosity POOL vs IND de novo",xlim=c(0,0.5), ylim=c(0,0.5))
text(x=0.05,y=0.475,labels=bquote(R^2==~.(rsquare)))
text(x=0.05,y=0.45,labels=bquote(p==~.(pvalue)))
abline(z,col="green")
# scatterplot IND vs POOL
plot(mean_het_exp_N$N_HET_EXP_mean,mean_het_exp_P$P_HET_EXP_mean,ylab="Mean Expected Heterozygosity POOL",xlab="Expected Heterozygosity IND",pch=19,las=1,cex=0.3,main="Mean Expected Heterozygosity POOL vs IND de novo",xlim=c(0,0.5), ylim=c(0,0.5))
text(x=0.05,y=0.475,labels=bquote(R^2==~.(rsquare)))
text(x=0.05,y=0.45,labels=bquote(p==~.(pvalue)))
abline(y,col="blue")
# close the pdf file
dev.off()



