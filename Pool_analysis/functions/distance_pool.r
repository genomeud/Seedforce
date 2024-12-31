suppressPackageStartupMessages({
  library(optparse)
})
# different options need to be adapt manually depending on the samples we want to analyse
option_list = list(
# here include the sample files to be analysed within the vcf
    make_option(c("-s", "--spp"), type="character", default="GPA", 
              help="name of the studied spp [default= %default]", metavar="character"),
    make_option(c("-I", "--indir"), type="character", default="/projects/marroni/seedforce/pool/de_novo_pipeline_catalog_spp_pops/GPA", 
               help="output file name [default= %default]", metavar="character"),
    make_option(c("-o", "--outfile"), type="character", default="/projects/marroni/seedforce/pool/de_novo_pipeline_catalog_spp_pops/GPA/PCA.pdf", 
              help="output file name [default= %default]", metavar="character"),
    make_option(c("-O", "--outpath"), type="character",
    default="/projects/marroni/seedforce/pool/de_novo_pipeline_catalog_spp_pops/GPA",help="Output directory [default= %default]", metavar="character"),
    make_option(c("-S","--step"), action="store",dest="to_do",default="distance", type='character', help="Step of the pipeline (distance,plot) [%default]")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

  if (is.null(opt$spp)) {
  stop("WARNING: No spp specified with '-s' flag.")
} else {  cat ("spp is ", opt$spp, "\n")
  spp <- opt$spp  
  }
  
  if (is.null(opt$indir)) {
  stop("WARNING: No indir specified with '-I' flag.")
} else {  cat ("indir is ", opt$indir, "\n")
  indir <- opt$indir  
  }

  if (is.null(opt$outfile)) {
  stop("WARNING: No outfile specified with '-o' flag.")
} else {  cat ("outfile is ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }
  
if (is.null(opt$outpath)) {
  stop("WARNING: No outpath specified with '-O' flag.")
} else {  cat ("outpath ", opt$outpath, "\n")
  outpath <- opt$outpath  
  }
 
  if (is.null(opt$to_do)) {
  stop("WARNING: No step specified with '-S' flag.")
} else {  cat ("Step to run is", opt$to_do, "\n")
  to_do <- opt$to_do  
  } 

run_distance<-function(spp)

{
library(data.table)
library(ggplot2)
library(dplyr)
library("ggrepel")
#We have a parameter so thar you do not hardcode the file name in the function
# remove all lines that contains ##
cat("Reading vcf...\n")
# spp=ACI
# browser()
invcf<-dir(indir,pattern=".vcf",full.names=T)
innames<-basename(invcf)
allvcf<-fread(cmd=paste("grep -v ^##", invcf[1]),data.table=F)
#Remove indels and triallelic SNPs in one shot
allvcf<-allvcf[nchar(allvcf$ALT)==1,]
names(allvcf)[10]<-innames[1]
ciccio<-unlist(lapply(strsplit(allvcf[,innames[1]],":"),"[",3))
ciccio1<-as.numeric(unlist(lapply(strsplit(ciccio,","),"[",1)))
ciccio2<-as.numeric(unlist(lapply(strsplit(ciccio,","),"[",2)))
allvcf[paste0("ALT_",innames[1])]<-ciccio2/(ciccio1+ciccio2)
allvcf<-allvcf[,!names(allvcf)%in%c("ID","QUAL","FILTER","INFO")]
for(aaa in 2:length(invcf))
{
tvcf<-fread(cmd=paste("grep -v ^##", invcf[aaa]),data.table=F)
#Remove indels and triallelic SNPs in one shot
tvcf<-tvcf[nchar(tvcf$ALT)==1,]
names(tvcf)[10]<-innames[aaa]
ciccio<-unlist(lapply(strsplit(tvcf[,innames[aaa]],":"),"[",3))
ciccio1<-as.numeric(unlist(lapply(strsplit(ciccio,","),"[",1)))
ciccio2<-as.numeric(unlist(lapply(strsplit(ciccio,","),"[",2)))
tvcf[paste0("ALT_",innames[aaa])]<-ciccio2/(ciccio1+ciccio2)
tvcf<-tvcf[,!names(tvcf)%in%c("ID","QUAL","FILTER","INFO")]
allvcf<-merge(allvcf,tvcf,by=c("#CHROM","POS","REF","ALT","FORMAT"),all=F,sort=F)
}
# browser()
mymat<-select(allvcf,starts_with("ALT_"))
write.table(mymat,paste(outpath,"PCA.all_pops.",spp,".txt",sep=""),sep="\t",quote=F,row.names=F,col.names=T)
}


run_plot<-function(spp)

{
library(data.table)
library(ggplot2)
library(dplyr)
library("ggrepel")
# # plot
# LPS

setwd(paste0("/projects/marroni/seedforce/pool/de_novo_pipeline_catalog_spp_pops/",spp,""))
mymat<-read.table(paste0("/projects/marroni/seedforce/pool/de_novo_pipeline_catalog_spp_pops/",spp,"/PCA.all_pops.",spp,"_names.txt"), header=T)
# colnames(mymat)<-c("GOZ1","GOZ2","GOZ3","LIN","LPS")
# mymat<-as.data.frame(mymat[,c(1,2,3,4)])
mymat<-t(mymat)
mypca<-prcomp(mymat)
#scores <- mypca$x[,1:3]            # scores for first three PC's
pca_data <- as.data.frame(mypca$x)
pop<-row.names(pca_data)
#scores$myvar<-row.names(scores)
pdf(outfile)
pippo <- ggplot(pca_data, aes(x = PC1, y = PC2, color = pop), size = 5) +
    geom_point(size = 3) +
    geom_text_repel(aes(label = pop), size = 3.5) +
    coord_fixed() +
    theme(legend.text = element_text(size = 9))
    # ggtitle(paste0("PCA ",spp,"")) +
    # plot.title = element_text(size=18)
print(pippo)
dev.off()
}


if (to_do=="distance")
{
    run_distance(spp=spp)
}   else if (to_do=="plot")
{
    run_plot(spp=spp)
}



