args <- commandArgs(trailingOnly = TRUE)

spp=args[1]
pop=args[2]

####################
# one plot per spp #
####################

# # # # ALLELE FREQUENCY AND HETEROZYGOSITY IN POOL WO REFERENCE# # # #
# Hardy-Weinberg-Equilibrium (HWE)
# p^2 +2pq + q^2 = 1
# p=frequency of allele homozygous dominant
# q=frequency of allele homozygous recesive (p-1)
# 2pq=frequency of heterozygous

# # merged pool

# setwd(paste0("/projects/marroni/seedforce/pool/de_novo_pipeline_catalog_spp_pops/",spp,"/"))
# tbl = read.table(paste0("",spp,"_FINAL.tbl"), header=F)
# names(tbl)=c("CHROM", "POS","REF", "ALT", "P_REF", "P_ALT", "Min_all")
# # calculate P_ALT_FREQ
# # P_ALT_FREQ=P_ALT/(P_REF+P_ALT)
# P_ALT_FREQ=tbl$P_ALT/(tbl$P_REF+tbl$P_ALT)
# # paste new column in the previous table and create a new table
# tbl2 <- cbind(tbl, P_ALT_FREQ)
# # calculate P_REF_FREQ
# P_REF_FREQ=(1-tbl2$P_ALT_FREQ)
# # paste new column in the previous table and create a new table
# tbl3 <- cbind(tbl2, P_REF_FREQ)
# # MAF
# MAF = tbl$Min_all/(tbl$P_REF+tbl$P_ALT)
# tbl4 <- cbind(tbl3, MAF)
# # calculate expected heterozygosity. 2pq
# P_HET_EXP = (2 * tbl4$P_ALT_FREQ * tbl4$P_REF_FREQ)
# tbl5 <- cbind(tbl4, P_HET_EXP)
# write.table(tbl5, file = paste0("",spp,"_allfreq_HE.tbl"), quote=FALSE, row.names=FALSE,sep="\t")


# histogram of the ALT allele frequency #

# create pdf output file. 
fileOut=(paste0("histogram_exp_het_pool_",spp,".pdf"))
pdf(fileOut, paper="A4")
# histogram with added parameters
hist(tbl5$MAF, breaks=10,cex.lab=1.3,cex.axis=1.2,cex.main=1.5,
main=paste0("MAF POOL ",spp,""),
xlab="MAF POOL",
ylab="Total SNPs",
xlim=c(0,0.5),
col="seagreen4",
# freq=FALSE
)
dev.off()

###############
# one per pop #
###############
# each pop separated
setwd(paste0("/projects/marroni/seedforce/pool/de_novo_pipeline_catalog_spp_pops/",spp,"/"))
tbl = read.table(paste0("",spp,"_",pop,"_FINAL.tbl"), header=F)
names(tbl)=c("CHROM", "POS","REF", "ALT", "P_REF", "P_ALT", "Min_all")
# calculate P_ALT_FREQ
# P_ALT_FREQ=P_ALT/(P_REF+P_ALT)
P_ALT_FREQ=tbl$P_ALT/(tbl$P_REF+tbl$P_ALT)
# paste new column in the previous table and create a new table
tbl2 <- cbind(tbl, P_ALT_FREQ)
# calculate P_REF_FREQ
P_REF_FREQ=(1-tbl2$P_ALT_FREQ)
# paste new column in the previous table and create a new table
tbl3 <- cbind(tbl2, P_REF_FREQ)
# MAF
MAF = tbl$Min_all/(tbl$P_REF+tbl$P_ALT)
tbl4 <- cbind(tbl3, MAF)
# calculate expected heterozygosity. 2pq
P_HET_EXP = (2 * tbl4$P_ALT_FREQ * tbl4$P_REF_FREQ)
tbl5 <- cbind(tbl4, P_HET_EXP)
write.table(tbl5, file = paste0("",spp,"_",pop,"_allfreq_HE.tbl"), quote=FALSE, row.names=FALSE,sep="\t")


# create pdf output file. 
fileOut=(paste0("histogram_MAF_pool_",spp,"_",pop,".pdf"))
pdf(fileOut, paper="A4")
# histogram with added parameters
hist(tbl5$MAF, breaks=10,cex.lab=1.3,cex.axis=1.2,cex.main=1.5,
main=paste0("MAF POOL ",spp," ",pop,""),
xlab="MAF POOL",
ylab="Total SNPs",
xlim=c(0,0.5),
col="coral2",
# freq=FALSE
)
dev.off()

