# Copyright:	    Fabio Marroni & Paloma Perez 2024
# Aim:              Create table with allele frequencies
# To add:		     
# Suggestions: 
# Fixes:  


suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
# here include the sample IDs to be analysed within the vcf
  make_option(c("-S", "--samples"), type="character", default="10-GAL-BCO1-UNIPA-2-B2,11-GAL-BCO1-UNIPA-3-C2,12-GAL-BCO1-UNIPA-4-D2,13-GAL-BCO1-UNIPA-5-E2,14-GAL-BCO1-UNIPA-6-F2,15-GAL-BCO1-UNIPA-7-G2,16-GAL-BCO1-UNIPA-8-H2,17-GAL-BCO2-UNIPA-1-A3,18-GAL-BCO2-UNIPA-2-B3,19-GAL-BCO2-UNIPA-3-C3,1-GAL-TR1-UNIPA-1-A1,20-GAL-BCO2-UNIPA-4-D3,21-GAL-BCO2-UNIPA-5-E3,22-GAL-BCO2-UNIPA-6-F3,23-GAL-BCO2-UNIPA-7-G3,24-GAL-BCO2-UNIPA-8-H3,25-GAL-SOS-UNIPA-1-A4,26-GAL-SOS-UNIPA-2-B4,27-GAL-SOS-UNIPA-3-C4,28-GAL-SOS-UNIPA-4-D4,29-GAL-SOS-UNIPA-5-E4,2-GAL-TR1-UNIPA-2-B1,30-GAL-SOS-UNIPA-6-F4,31-GAL-SOS-UNIPA-7-G4,32-GAL-SOS-UNIPA-8-H4,33-GAL-BIA-UNIPA-1-A5,34-GAL-BIA-UNIPA-2-B5,35-GAL-BIA-UNIPA-3-C5,36-GAL-BIA-UNIPA-4-D5,37-GAL-BIA-UNIPA-5-E5,38-GAL-BIA-UNIPA-6-F5,39-GAL-BIA-UNIPA-7-G5,3-GAL-TR1-UNIPA-3-C1,40-GAL-BIA-UNIPA-8-H5,41-GAL-BA1-UNIPA-1-A6,42-GAL-BA1-UNIPA-2-B6,43-GAL-BA1-UNIPA-3-C6,44-GAL-BA1-UNIPA-4-D6,45-GAL-BA1-UNIPA-5-E6,46-GAL-BA1-UNIPA-6-F6,47-GAL-BA1-UNIPA-7-G6,48-GAL-BA1-UNIPA-8-H6,49-GAL-BA4-UNIPA-1-A7,4-GAL-TR1-UNIPA-4-D1,50-GAL-BA4-UNIPA-2-B7,51-GAL-BA4-UNIPA-3-C7,52-GAL-BA4-UNIPA-4-D7,53-GAL-BA4-UNIPA-5-E7,54-GAL-BA4-UNIPA-6-F7,55-GAL-BA4-UNIPA-7-G7,56-GAL-BA4-UNIPA-8-H7,57-GAL-BA5-UNIPA-1-A8,58-GAL-BA5-UNIPA-2-B8,59-GAL-BA5-UNIPA-3-C8,5-GAL-TR1-UNIPA-5-E1,60-GAL-BA5-UNIPA-4-D8,61-GAL-BA5-UNIPA-5-E8,62-GAL-BA5-UNIPA-6-F8,63-GAL-BA5-UNIPA-7-G8,64-GAL-BA5-UNIPA-8-H8,65-GAL-CUS-UNIPA-1-A9,66-GAL-CUS-UNIPA-2-B9,67-GAL-CUS-UNIPA-3-C9,68-GAL-CUS-UNIPA-4-D9,69-GAL-CUS-UNIPA-5-E9,6-GAL-TR1-UNIPA-6-F1,70-GAL-CUS-UNIPA-6-F9,71-GAL-CUS-UNIPA-7-G9,72-GAL-CUS-UNIPA-8-H9,73-GAL-3F-UNIPA-1-A10,74-GAL-3F-UNIPA-2-B10,75-GAL-3F-UNIPA-3-C10,76-GAL-3F-UNIPA-4-D10,77-GAL-3F-UNIPA-5-E10,78-GAL-3F-UNIPA-6-F10,79-GAL-3F-UNIPA-7-G10,7-GAL-TR1-UNIPA-7-G1,80-GAL-3F-UNIPA-8-H10,81-GAL-BA2-UNIPA-1-A11,82-GAL-BA2-UNIPA-2-B11,83-GAL-BA2-UNIPA-3-C11,84-GAL-BA2-UNIPA-4-D11,85-GAL-BA2-UNIPA-5-E11,86-GAL-BA2-UNIPA-6-F11,87-GAL-BA6-UNIPA-1-G11,88-GAL-BA6-UNIPA-2-H11,89-GAL-BA6-UNIPA-3-A12,8-GAL-TR1-UNIPA-8-H1,90-GAL-BA6-UNIPA-4-B12,91-GAL-BA6-UNIPA-5-C12,92-GAL-BA6-UNIPA-6-D12,93-GAL-BCO2-UNIPA-9-E12,94-GAL-BIA-UNIPA-9-F12,95-GAL-BA1-UNIPA-9-G12,96-GAL-CUS-UNIPA-9-H12,9-GAL-BCO1-UNIPA-1-A2", 
              help="Comma separated list of samples to keep for analysis [default= %default]", metavar="character"),
# put the path and file with vcf for individuals approach
  make_option(c("-V", "--vcffile"), type="character", default="populations.snps.filtered.recode.vcf", 
              help="path to vcf input file [default= %default]", metavar="character"),
# put path and name for the output file
  make_option(c("-O", "--outfile"), type="character", default="allele_frequency.tbl", 
              help="output file name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


  if (is.null(opt$samples)) {
  stop("WARNING: No samples specified with '-S' flag.")
} else {  cat ("samples is ", opt$samples, "\n")
  samples <- opt$samples  
  }

  if (is.null(opt$vcffile)) {
  stop("WARNING: No vcffile specified with '-V' flag.")
} else {  cat ("vcffile is ", opt$vcffile, "\n")
  vcffile <- opt$vcffile  
  }

  if (is.null(opt$outfile)) {
  stop("WARNING: No outfile specified with '-O' flag.")
} else {  cat ("outfile is ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }
  

filter_vcf<-function()
{
library(data.table)
#We have a parameter so that you do not hardcode the file name in the function
# remove all lines that contains ##
vcf<-fread(cmd=paste("grep -v ^##", vcffile),data.table=F)
cat("Selecting samples...\n")
# each sample is sepparated by commas
mysample<-unlist(strsplit(samples,","))
# first column header
# firstcol<-c("#CHROM","POS","REF","ALT","QUAL","INFO","FORMAT")
firstcol<-c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT")

allcol<-c(firstcol,mysample)
myvcf<-vcf[,names(vcf)%in%allcol]
# orig<-myvcf
# myvcf<-myvcf[1:100,]
# create 3 columns one for ind REF, one for ind_ALT and one for ind missing data 
myvcf$N_REF<-0
myvcf$N_ALT<-0
myvcf$N_MIS<-0


# for each sample in my sample do the following loop
for(aaa in 1:length(mysample))
{
# get the haplotype info which could be 0/0 when both haplotypes are like ref; 1/1 when both are like alt and 0/1 when we have 1 ref and 1 alt. Then we will sum how many times we have one or the other (ref or alt) and how many times is missing ./.
tmp<-unlist(lapply(strsplit(myvcf[,mysample[aaa]],":"),"[",1))
all1<-unlist(lapply(strsplit(tmp,"/"),"[",1))
all2<-unlist(lapply(strsplit(tmp,"/"),"[",2))
myvcf$N_REF<-apply(cbind(myvcf$N_REF,all1=="0"),1,sum)
myvcf$N_REF<-apply(cbind(myvcf$N_REF,all2=="0"),1,sum)
myvcf$N_ALT<-apply(cbind(myvcf$N_ALT,all1=="1"),1,sum)
myvcf$N_ALT<-apply(cbind(myvcf$N_ALT,all2=="1"),1,sum)
myvcf$N_MIS<-apply(cbind(myvcf$N_MIS,all1=="."),1,sum)
myvcf$N_MIS<-apply(cbind(myvcf$N_MIS,all2=="."),1,sum)
}
myvcf<-myvcf[,!names(myvcf)%in%mysample]

write.table(myvcf,outfile,quote=F,sep="\t",row.names=F)
cat("Done\n")
}

filter_vcf()
