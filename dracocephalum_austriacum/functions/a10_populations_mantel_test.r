# Copyright:	    Paloma Perez & Fabio Marroni 2024
# Aim:              mantel test between genetic distances (measured by fst) and geographic coordinates and create scatterplot with pairwise comparisons
# To add:		     
# Suggestions: 
# Fixes: 

# Run with --help or -h flag for help.

suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-F", "--genfile"), type="character",
  default="populations.fst_summary_for_heatmap.tsv",help="FST matrix with FST between populations [default= %default]", metavar="character"),
  make_option(c("-G", "--geofile"), type="character",
  default="pop_geographic_coordinates.txt",help="Geographic coordinates for each population[default= %default]", metavar="character"),
  make_option(c("-O", "--outpath"), type="character",
  default="NULL",help="Output directory [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$genfile)) {
  stop("WARNING: No FST file specified with '-F' flag.")
} else {  cat ("FST matrix file ", opt$genfile, "\n")
  genfile <- opt$genfile  
  }

if (is.null(opt$geofile)) {
  stop("WARNING: No geofile with geographic coordinates fopr each population specified with '-G' flag. Desired input file contains Latitude\tLongitude\tpopulation")
} else {  cat ("geofile ", opt$geofile, "\n")
  geofile <- opt$geofile  
  }

if (is.null(opt$outpath)) {
  stop("WARNING: No outpath specified with '-O' flag.")
} else {  cat ("outpath ", opt$outpath, "\n")
  outpath <- opt$outpath  
  }
  
library(vegan)
library(geosphere)

fst_data  <- read.table(genfile, header=T, row.names=1, sep="\t")
coords  <- read.table(geofile, header=T, sep="\t")

fst_matrix <- as.matrix(fst_data)

# Copy upper triangle values to the lower triangle to symmetrize the matrix
fst_matrix[lower.tri(fst_matrix)] <- t(fst_matrix)[lower.tri(fst_matrix)]

# The diagonal should remain zeros, so we do not need to change it
diag(fst_matrix) <- 0

# Ensure that the populations are in the same order in both the FST matrix and geographic coordinates
coords <- coords[match(rownames(fst_matrix), coords$population), ]

# Calculate geographic distance matrix
geo_dist_matrix <- distm(coords[, c("Longitude", "Latitude")], fun = distHaversine)

rownames(geo_dist_matrix) <- coords$population
colnames(geo_dist_matrix) <- coords$population

# mantel_test <- mantel(fst_matrix, geo_dist_matrix, method="spearman")
mantel_test <- mantel(fst_matrix, geo_dist_matrix, method = "pearson", permutations = 9999)
print(mantel_test)

# scatterplot

# Flatten the distance matrices into vectors
genetic_dist <- as.vector(fst_matrix)
geographic_dist <- as.vector(geo_dist_matrix/1000)
# Filter out the self-comparison points (where both distances are 0)
valid_points <- geographic_dist != 0 & genetic_dist != 0

# Create scatterplot
fileOut=("scatterplot_FST_geographic_populations_dracocephalum.png")
png(fileOut,width=13,height=14,units="cm",res=1600, type="cairo")
par(cex.axis = 0.7, cex.lab = 1)
par(cex.main = 1.1)
par(mar = c(3, 3, 2, 2),oma=c(1.5,0.3,0.3,0.1), mgp=c(1.7,0.7,0))
plot(geographic_dist[valid_points], genetic_dist[valid_points],
     xlab = "Geographic Distance (Km)",
     ylab = expression("Genetic Distance (F" [ST] * ")"),
     # main = "Genetic Distance vs Geographic Distance",
     main = c(expression(paste(italic("D. austriacum")))),
     # sub = expression(italic("P. palinuri")),
     pch = 2, cex.main=1.2, col = "firebrick", cex=1.4)

# mtext(expression(bold(italic("D. austriacum"))), side = 3, line = 0.1, cex = 1.2)

# Add a regression line
abline(lm(genetic_dist[valid_points] ~ geographic_dist[valid_points]), col = "black")

# Add Mantel statistic and significance to the plot
text(x = 347, y = 0.06,
     labels = paste("r =", round(mantel_test$statistic, 4), "\np =", round(mantel_test$signif, 3)),
     cex = 0.9, col = "black")
# pos=3
# cex_text=1.2
# mtext(paste0("", toupper(LETTERS[pos]), ")"), 2, line=1, at=0, padj=-40.5, adj=1, cex=cex_text, las=2, font=2)
dev.off()
     
# # # # # FILTER CN1
# Load data
fst_data  <- read.table("populations.fst_summary_for_heatmap.tsv", header=TRUE, row.names=1, sep="\t")
coords  <- read.table("pop_geographic_coordinates.txt", header=TRUE, sep="\t")

# Remove CN1 from both datasets
fst_data <- fst_data[!rownames(fst_data) %in% "CN1", !(colnames(fst_data) %in% "CN1")]
coords   <- coords[!coords$population %in% "CN1", ]

# Convert to matrix
fst_matrix <- as.matrix(fst_data)

# Symmetrize the matrix
fst_matrix[lower.tri(fst_matrix)] <- t(fst_matrix)[lower.tri(fst_matrix)]
diag(fst_matrix) <- 0

# Ensure same order in both datasets
coords <- coords[match(rownames(fst_matrix), coords$population), ]

# Geographic distance matrix
geo_dist_matrix <- distm(coords[, c("Longitude", "Latitude")], fun = distHaversine)

rownames(geo_dist_matrix) <- coords$population
colnames(geo_dist_matrix) <- coords$population

# mantel_test <- mantel(fst_matrix, geo_dist_matrix, method="spearman")
mantel_test <- mantel(fst_matrix, geo_dist_matrix, method = "pearson", permutations = 9999)
print(mantel_test)

# scatterplot

# Flatten the distance matrices into vectors
genetic_dist <- as.vector(fst_matrix)
geographic_dist <- as.vector(geo_dist_matrix/1000)
# Filter out the self-comparison points (where both distances are 0)
valid_points <- geographic_dist != 0 & genetic_dist != 0

# Create scatterplot
fileOut=("scatterplot_FST_geographic_populations_dracocephalum_noCN1.png")
png(fileOut,width=13,height=14,units="cm",res=1600, type="cairo")
par(cex.axis = 0.7, cex.lab = 1)
par(cex.main = 1.1)
par(mar = c(3, 3, 2, 2),oma=c(1.5,0.3,0.3,0.1), mgp=c(1.7,0.7,0))
plot(geographic_dist[valid_points], genetic_dist[valid_points],
     xlab = "Geographic Distance (Km)",
     ylab = expression("Genetic Distance (F" [ST] * ")"),
     # main = "Genetic Distance vs Geographic Distance",
     main = c(expression(paste(italic("D. austriacum")))),
     # sub = expression(italic("P. palinuri")),
     pch = 2, cex.main=1.2, col = "firebrick", cex=1.4)

# mtext(expression(bold(italic("D. austriacum"))), side = 3, line = 0.1, cex = 1.2)

# Add a regression line
abline(lm(genetic_dist[valid_points] ~ geographic_dist[valid_points]), col = "black")

# Add Mantel statistic and significance to the plot
text(x = 80, y = 0.16,
     labels = paste("r =", round(mantel_test$statistic, 4), "\np =", round(mantel_test$signif, 3)),
     cex = 0.9, col = "black")
pos=3
cex_text=1.2
# mtext(paste0("", toupper(LETTERS[pos]), ")"), 2, line=1, at=0, padj=-40.5, adj=1, cex=cex_text, las=2, font=2)
dev.off()
