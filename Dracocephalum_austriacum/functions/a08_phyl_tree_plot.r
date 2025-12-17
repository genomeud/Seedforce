# Copyright:	    Paloma perez & Fabio Marroni 2024
# Aim:              Phylogenetic trees based on neighbor joining distance method. Unrooted and rooted phylogenetic trees
# To add:		     
# Suggestions: 
# Fixes:  


suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-V", "--vcffile"), type="character",
  default="/projects/marroni/seedforce/dracocephalum_austriacum/github/tests/populations.snps.filtered.recode_MIS_filt_header.vcf",help="vcf file with header [default= %default]", metavar="character"),
  make_option(c("-M", "--metadata"), type="character",
  default="/projects/marroni/seedforce/dracocephalum_austriacum/raw_reads/sample_ID_and_populations.txt",help="metadata with individual IDs and pop names [default= %default]", metavar="character"),
  make_option(c("-O", "--outdir"), type="character",
  default="phyl_tree_NJ/",help="Output directory for phylogenetic trees [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$vcffile)) {
  stop("WARNING: No vcf file specified with '-V' flag.")
} else {  cat ("vcf file ", opt$vcffile, "\n")
  vcffile <- opt$vcffile  
  }
  
if (is.null(opt$metadata)) {
  stop("WARNING: No metadata file specified with '-M' flag.")
} else {  cat ("metadata file ", opt$metadata, "\n")
  metadata <- opt$metadata  
  }
  
if (is.null(opt$outdir)) {
  stop("WARNING: No output directory specified with '-O' flag.")
} else {  cat ("output directory ", opt$outdir, "\n")
  outdir <- opt$outdir  
  }

library(vcfR)
library(ape)

dir.create(outdir,showWarnings=F,recursive=TRUE)

vcf <- read.vcfR(vcffile)
# extract genotype information for each individual
genotypes <- extract.gt(vcf, element = "GT", as.numeric = TRUE)
# convert genotype matrix to distant matrix
dist_matrix <- dist(t(genotypes))
# construct phylogenetic tree with Neighbor-Joining
nj_tree <- nj(dist_matrix)
tree <- nj(dist_matrix)
# load metadata with group info.
metadata <- read.table(metadata, header=TRUE)
colnames(metadata)<-NULL
colnames(metadata)<-c("sample_id","group")
# put pop names in alphabetic order! Use same colors that were previously defined for other analyses (PCA, IBD, HWE...)
group_colors <- c("CH1" = "darkblue","CH2" = "darkmagenta","BZ1" = "maroon", "BZ2" = "palevioletred2","SO1" = "slateblue1","COR" = "peachpuff","MAL" = "skyblue1","CN1" = "mediumaquamarine", "ALL"="white")
colors=c("darkblue","darkmagenta","maroon","palevioletred2","slateblue1","peachpuff","skyblue1","mediumaquamarine")
pops <- c("CH1","CH2","BZ1","BZ2","SO1","COR","MAL","CN1")
metadata$group <- factor(metadata$group, levels = pops)

# ROOTED TREE
# out_root=file.path(outdir,"phylogenetic_tree_with_legend_rooted.pdf")
# pdf(out_root, width=4.5, height=7.5)
tip_colors <- group_colors[metadata$group[match(tree$tip.label, metadata$sample_id)]]
plot.phylo(tree, tip.color=tip_colors, cex=0.45, label.offset=0.05, main="Phylogenetic Tree")
legend("bottomright", legend=names(group_colors), fill=group_colors, cex=0.7, border=FALSE, bty="n")
dev.off()

# UNROOTED TREE
outdir <- "/projects/marroni/seedforce/dracocephalum_austriacum/github/tests"
out_unroot=file.path(outdir,"phylogenetic_tree_with_legend_unrooted.jpeg")
jpeg(out_unroot,width=14,height=12,units="cm",res=1200, type="cairo")
par(mar=c(1.6,3,1.2,0.1), mgp=c(1,0.2,0),tck=-0.03,oma=c(1.2,0.1,0.3,0.1))
plot.phylo(tree, tip.color=tip_colors, type = "unrooted", cex=0.45, label.offset=0.05)
legend("bottomright", legend=names(group_colors), fill=group_colors, cex=0.7, border=FALSE, bty="n")
# add the panel letter in the format "(b)"
pos=2
cex_text=1
mtext(paste0("(", tolower(LETTERS[pos]), ")"), 2, line=1, at=0, padj=-22.5, adj=1, cex=cex_text, las=2, font=2)
dev.off()


