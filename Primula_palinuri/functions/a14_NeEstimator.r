# Copyright:	    Paloma Perez & Fabio Marroni 2025
# Aim:              Prepare input files for NeEstimator. From vcf format to GENEPOP format
# To add:		     
# Suggestions: 
# Fixes:  

R

library(vcfR)
library(adegenet)

# load dartR functions
# Path to the folder containing your .R files
folder_path <- "/projects/marroni/seedforce/softwares/dartR_R_library/dartR/R"
# List all .R files in the folder
r_files <- list.files(path = folder_path, pattern = "\\.r$", full.names = TRUE)
R_files <- list.files(path = folder_path, pattern = "\\.R$", full.names = TRUE)
# Source each file with the functions
lapply(r_files, source)
lapply(R_files, source)

# Read the VCF file
# # # # LAM
pop <- c("LAM")
input.vcf <- paste0("filtered_",pop,".recode_inALLinds.vcf")


vcf <- read.vcfR(input.vcf)

# Convert VCF to a genlight object (adegenet)
x <- vcfR2genlight(vcf)
pop(x) <- c("LAM","LAM","LAM","LAM","LAM","LAM","LAM","LAM")

ploidy(x) <- 2

x <- gl.filter.allna(x)

# genelight to genepop
gl2genepop(
  x,
  outfile = paste0("",pop,".gen"),
  outpath = "Ne_effective_pop_size/NeEstimator_V2/pops_sep/",
  pop_order = "alphabetic",
  output_format = "2_digits",
  verbose = NULL
)

# # # # FIU
pop <- c("FIU")
input.vcf <- paste0("filtered_",pop,".recode_inALLinds.vcf")
vcf <- read.vcfR(input.vcf)

# Convert VCF to a genlight object (adegenet)
x <- vcfR2genlight(vcf)
pop(x) <- c("FIU","FIU","FIU","FIU","FIU","FIU","FIU","FIU","FIU","FIU","FIU","FIU")

ploidy(x) <- 2 

x <- gl.filter.allna(x)

# genelight to genepop
gl2genepop(
  x,
  outfile = paste0("",pop,".gen"),
  outpath = "Ne_effective_pop_size/NeEstimator_V2/pops_sep/",
  pop_order = "alphabetic",
  output_format = "2_digits",
  verbose = NULL
)

# # # # ID
pop <- c("ID")

pop <- c("ID")
input.vcf <- paste0("filtered_",pop,".recode_inALLinds.vcf")
vcf <- read.vcfR(input.vcf)

# Convert VCF to a genlight object (adegenet)
x <- vcfR2genlight(vcf)
pop(x) <- c("ID","ID","ID","ID","ID","ID","ID","ID","ID")

ploidy(x) <- 2 

x <- gl.filter.allna(x)

# genelight to genepop
gl2genepop(
  x,
  outfile = paste0("",pop,".gen"),
  outpath = "Ne_effective_pop_size/NeEstimator_V2/pops_sep/",
  pop_order = "alphabetic",
  output_format = "2_digits",
  verbose = NULL
)

# # # # PC
pop <- c("PC")

pop <- c("PC")
input.vcf <- paste0("filtered_",pop,".recode_inALLinds.vcf")
vcf <- read.vcfR(input.vcf)

# Convert VCF to a genlight object (adegenet)
x <- vcfR2genlight(vcf)
pop(x) <- c("PC","PC","PC","PC","PC","PC","PC","PC","PC","PC","PC","PC")

ploidy(x) <- 2 

x <- gl.filter.allna(x)

# genelight to genepop
gl2genepop(
  x,
  outfile = paste0("",pop,".gen"),
  outpath = "Ne_effective_pop_size/NeEstimator_V2/pops_sep/",
  pop_order = "alphabetic",
  output_format = "2_digits",
  verbose = NULL
)

# # # # SGP
pop <- c("SGP")
input.vcf <- paste0("filtered_",pop,".recode_inALLinds.vcf")
vcf <- read.vcfR(input.vcf)

# Convert VCF to a genlight object (adegenet)
x <- vcfR2genlight(vcf)
pop(x) <- c("SGP","SGP","SGP","SGP","SGP","SGP","SGP","SGP","SGP","SGP","SGP","SGP")

ploidy(x) <- 2 

x <- gl.filter.allna(x)

# genelight to genepop
gl2genepop(
  x,
  outfile = paste0("",pop,".gen"),
  outpath = "Ne_effective_pop_size/NeEstimator_V2/pops_sep/",
  pop_order = "alphabetic",
  output_format = "2_digits",
  verbose = NULL
)

# # # # CIM
pop <- c("CIM")
input.vcf <- paste0("filtered_",pop,".recode_inALLinds.vcf")
vcf <- read.vcfR(input.vcf)

# Convert VCF to a genlight object (adegenet)
x <- vcfR2genlight(vcf)
pop(x) <- c("CIM","CIM","CIM","CIM","CIM","CIM","CIM","CIM","CIM","CIM","CIM","CIM","CIM","CIM","CIM","CIM","CIM","CIM","CIM","CIM")

ploidy(x) <- 2 

x <- gl.filter.allna(x)

# genelight to genepop
gl2genepop(
  x,
  outfile = paste0("",pop,".gen"),
  outpath = "Ne_effective_pop_size/NeEstimator_V2/pops_sep/",
  pop_order = "alphabetic",
  output_format = "2_digits",
  verbose = NULL
)

# # # # PPA
pop <- c("PPA")
input.vcf <- paste0("filtered_",pop,".recode_inALLinds.vcf")
vcf <- read.vcfR(input.vcf)

# Convert VCF to a genlight object (adegenet)
x <- vcfR2genlight(vcf)
pop(x) <- c("PPA","PPA","PPA","PPA","PPA","PPA","PPA","PPA","PPA","PPA","PPA","PPA","PPA","PPA","PPA","PPA","PPA","PPA")

ploidy(x) <- 2 

x <- gl.filter.allna(x)

# genelight to genepop
gl2genepop(
  x,
  outfile = paste0("",pop,".gen"),
  outpath = "Ne_effective_pop_size/NeEstimator_V2/pops_sep/",
  pop_order = "alphabetic",
  output_format = "2_digits",
  verbose = NULL
)
