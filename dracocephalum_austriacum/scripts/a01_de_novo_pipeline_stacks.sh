# Copyright:	    Paloma Perez & Fabio Marroni 2025
# Aim:              De novo pipeline for building loci from ddRAD paired end experiments using STACKS + filtering of genotyping data to ensure high quality + genetic analyses including MAF, mean expected and observed heterozygosity, FST, IBD, PCA
# To add:		     
# Suggestions: 
# Fixes:  

spp=dracocephalum_austriacum

# create output folders
mkdir -p ${INPUT_DIR}/${spp}/de_novo_pipeline/populations_stacks

out1=${INPUT_DIR}/${spp}/de_novo_pipeline
out2=${INPUT_DIR}/${spp}/de_novo_pipeline/populations_stacks
# folder containing fastq sequences
raw_files=${INPUT_DIR}/${spp}/raw_reads
# Folder containing the functions
FUNC_DIR=functions

# Required population file inside raw_files folder 
# population_names.txt example:
# ID1 pop1
# ID2 pop1
# ID3 pop2
# ID4 pop2

###############################
# Run STACKS de novo pipeline #
###############################

perl ${STACKS_path}/denovo_map.pl -M 2 -T 24 -o ${out1} --popmap ${raw_files}/population_names.txt --samples ${raw_files} --paired &> ${out1}/de_novo_pipeline.log

# Run populations. At least 75% of the individuals with SNP information. hwe for hardy weinberg equilibrium estimates and fst for F statistics
${STACKS_path}/populations -P ${out1} -O ${out2} -M ${raw_files}/population_names.txt -t 24 --fstats --hwe --vcf -R 0.75 --max-obs-het 0.8 &> ${out2}/populations_stacks_pipeline.log

# Filter vcf file for individuals with at least 5 reads
vcftools --vcf ${out2}/populations.snps.vcf --minDP 5 --recode --out ${out2}/populations.snps.filtered &> ${out2}/vcf_filter.log

# Get allele count table from filtered vcf considering all the individuals
Rscript ${FUNC_DIR}/a01_allele_frequencies_individuals.r \
-S "10-DRA-CH2-MUSE-2-B2,11-DRA-CH2-MUSE-3-C2,12-DRA-CH2-MUSE-4-D2,13-DRA-CH2-MUSE-5-E2,14-DRA-CH2-MUSE-6-F2,15-DRA-CH2-MUSE-7-G2,16-DRA-CH2-MUSE-8-H2,17-DRA-SO1-MUSE-1-A3,18-DRA-SO1-MUSE-2-B3,19-DRA-SO1-MUSE-3-C3,1-DRA-CH1-MUSE-1-A1,20-DRA-SO1-MUSE-4-D3,21-DRA-SO1-MUSE-5-E3,22-DRA-SO1-MUSE-6-F3,23-DRA-SO1-MUSE-7-G3,24-DRA-SO1-MUSE-8-H3,25-DRA-CN1-MUSE-1-A4,26-DRA-CN1-MUSE-2-B4,27-DRA-CN1-MUSE-3-C4,28-DRA-CN1-MUSE-4-D4,29-DRA-CN1-MUSE-5-E4,2-DRA-CH1-MUSE-2-B1,30-DRA-CN1-MUSE-6-F4,31-DRA-CN1-MUSE-7-G4,32-DRA-CN1-MUSE-8-H4,33-DRA-MAL-MUSE-1-A5,34-DRA-MAL-MUSE-2-B5,35-DRA-MAL-MUSE-3-C5,36-DRA-MAL-MUSE-4-D5,37-DRA-MAL-MUSE-5-E5,38-DRA-MAL-MUSE-6-F5,39-DRA-MAL-MUSE-7-G5,3-DRA-CH1-MUSE-3-C1,40-DRA-MAL-MUSE-8-H5,41-DRA-COR-MUSE-1-A6,42-DRA-COR-MUSE-2-B6,43-DRA-COR-MUSE-3-C6,44-DRA-COR-MUSE-4-D6,45-DRA-COR-MUSE-5-E6,46-DRA-COR-MUSE-6-F6,47-DRA-COR-MUSE-7-G6,48-DRA-COR-MUSE-8-H6,49-DRA-SO1-MUSE-16-A7,4-DRA-CH1-MUSE-4-D1,50-DRA-SO1-MUSE-15-B7,51-DRA-SO1-MUSE-13-C7,52-DRA-SO1-MUSE-12-D7,53-DRA-SO1-MUSE-10-E7,54-DRA-COR-MUSE-15-F7,55-DRA-COR-MUSE-14-G7,56-DRA-COR-MUSE-13-H7,57-DRA-COR-MUSE-12-A8,58-DRA-COR-MUSE-11-B8,59-DRA-COR-MUSE-10-C8,5-DRA-CH1-MUSE-5-E1,60-DRA-COR-MUSE-9-D8,61-DRA-MAL-MUSE-15-E8,62-DRA-MAL-MUSE-14-F8,63-DRA-MAL-MUSE-12-G8,64-DRA-MAL-MUSE-11-H8,65-DRA-MAL-MUSE-10-A9,66-DRA-MAL-MUSE-9-B9,67-DRA-BZ1-MUSE-1-C9,68-DRA-BZ1-MUSE-2-D9,69-DRA-BZ1-MUSE-3-E9,6-DRA-CH1-MUSE-6-F1,70-DRA-BZ1-MUSE-4-F9,71-DRA-BZ1-MUSE-5-G9,72-DRA-BZ1-MUSE-6-H9,73-DRA-BZ2-MUSE-1-A10,74-DRA-BZ2-MUSE-2-B10,75-DRA-BZ2-MUSE-3-C10,77-DRA-BZ2-MUSE-5-E10,78-DRA-BZ2-MUSE-6-F10,79-DRA-BZ2-MUSE-7-G10,7-DRA-CH1-MUSE-7-G1,80-DRA-CH1-MUSE-9-H10,81-DRA-CH1-MUSE-10-A11,82-DRA-CH1-MUSE-11-B11,83-DRA-CH1-MUSE-12-C11,84-DRA-CH1-MUSE-13-D11,85-DRA-CH1-MUSE-14-E11,86-DRA-CH1-MUSE-15-F11,87-DRA-CH1-MUSE-16-G11,88-DRA-CH1-MUSE-17-H11,89-DRA-CH1-MUSE-18-A12,8-DRA-CH1-MUSE-8-H1,90-DRA-CH1-MUSE-19-B12,91-DRA-CH1-MUSE-20-C12,9-DRA-CH2-MUSE-1-A2" \
-V ${out2}/populations.snps.filtered.recode.vcf \
-O ${out2}/allele_info.tbl &> ${out2}/allele_info.log

# Create column with minor allele considering all individuals together
sed 1d ${out2}/allele_info.tbl | \
awk 'OFS="\t" {if($NF<($(NF-2)+$(NF-1))) print $0}' | awk 'OFS="\t" {if($10<$11)print $0,$10;if($11<$10)print $0,$11;if($10==$11)print $0,$10}' | cat <(paste <(head -1 ${out2}/allele_info.tbl) <( echo -e "MinAll")) - > ${out2}/allele_info_MinAll.tbl

# Calculate minor allele frequency, expected heterozygosity and histogram with minor allele frequency
Rscript ${FUNC_DIR}/a02_HWE_analysis.r \
-F ${out2}/allele_info_MinAll.tbl \
-O ${out2}/allele_freq_and_exp_het.tbl \
-P ${out2}/histogram_MAF.jpeg


# Create column with pop name = ALL. For some plots we will include all the individuals together as a unique sample
sed 1d ${out2}/allele_freq_and_exp_het.tbl | \
awk 'OFS="\t" {print $0,"ALL"}' | \
cat <(paste <(head -1 ${out2}/allele_freq_and_exp_het.tbl) <(echo -e "population")) - > ${out2}/allele_freq_and_exp_het_All.tbl

# Filter positions prior to divide by populations. We only want positions that were genotyped with at least 5 reads in 50% of the individuals. file without vcf header "##"
join -a 1 -1 1 -2 1 <(awk 'OFS="\t" {print $1"|"$2}' ${out2}/allele_freq_and_exp_het.tbl | sed 1d | sort -k1,1 ) <(grep -v "##" ${out2}/populations.snps.filtered.recode.vcf | sed 1d | awk 'OFS="\t" {print $1"|"$2,$0}' | sort -k1,1) -t $'\t' | cut -f1 --complement | sort -n -k1,2 | \
cat <(grep -v "##" ${out2}/populations.snps.filtered.recode.vcf | head -1) - > ${out2}/populations.snps.filtered.recode_MIS_filt.vcf

# keep normal vcf format for other analyses
cat <(grep "##" ${out2}/populations.snps.filtered.recode.vcf) ${out2}/populations.snps.filtered.recode_MIS_filt.vcf > ${out2}/populations.snps.filtered.recode_MIS_filt_header.vcf

# Divide the vcf by populations. This will generate a separate file for each population + histogram with minor allele frequency for every population.
# output file for each population: the number of reference alleles (POP1_REF); the number of alternative alleles (POP1_ALT); the number of missing alleles in the population (POP1_MIS); the total number of alleles for that population (POP1_REF + POP1_ALT + POP1_MIS) (POP1_N); expected heterozygosity (POP1_He); count for the minor allele (POP1_minALL) and the minor allele frequency (POP1_minALLfreq); population ID (Popname)
# otput format example
# #CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT          POP1_REF  POP1_ALT  POP1_MIS  POP1_N  POP1_He  POP1_minALL  POP1_minALLfreq     Popname
# 2       215  .   C    T    .     PASS    .     GT:DP:AD:GQ:GL  23        15        2         40      0.43125  15           0.394736842105263   POP1
# 2       30   .   G    T    .     PASS    .     GT:DP:AD:GQ:GL  38        0         2         40      0        0            0                   POP1
# 2       96   .   C    T    .     PASS    .     GT:DP:AD:GQ:GL  38        0         2         40      0        0            0                   POP1

cd ${out2}
Rscript ${FUNC_DIR}/a03_divide_populations_and_HWE_analysis.r \
-S "10-DRA-CH2-MUSE-2-B2,11-DRA-CH2-MUSE-3-C2,12-DRA-CH2-MUSE-4-D2,13-DRA-CH2-MUSE-5-E2,14-DRA-CH2-MUSE-6-F2,15-DRA-CH2-MUSE-7-G2,16-DRA-CH2-MUSE-8-H2,17-DRA-SO1-MUSE-1-A3,18-DRA-SO1-MUSE-2-B3,19-DRA-SO1-MUSE-3-C3,1-DRA-CH1-MUSE-1-A1,20-DRA-SO1-MUSE-4-D3,21-DRA-SO1-MUSE-5-E3,22-DRA-SO1-MUSE-6-F3,23-DRA-SO1-MUSE-7-G3,24-DRA-SO1-MUSE-8-H3,25-DRA-CN1-MUSE-1-A4,26-DRA-CN1-MUSE-2-B4,27-DRA-CN1-MUSE-3-C4,28-DRA-CN1-MUSE-4-D4,29-DRA-CN1-MUSE-5-E4,2-DRA-CH1-MUSE-2-B1,30-DRA-CN1-MUSE-6-F4,31-DRA-CN1-MUSE-7-G4,32-DRA-CN1-MUSE-8-H4,33-DRA-MAL-MUSE-1-A5,34-DRA-MAL-MUSE-2-B5,35-DRA-MAL-MUSE-3-C5,36-DRA-MAL-MUSE-4-D5,37-DRA-MAL-MUSE-5-E5,38-DRA-MAL-MUSE-6-F5,39-DRA-MAL-MUSE-7-G5,3-DRA-CH1-MUSE-3-C1,40-DRA-MAL-MUSE-8-H5,41-DRA-COR-MUSE-1-A6,42-DRA-COR-MUSE-2-B6,43-DRA-COR-MUSE-3-C6,44-DRA-COR-MUSE-4-D6,45-DRA-COR-MUSE-5-E6,46-DRA-COR-MUSE-6-F6,47-DRA-COR-MUSE-7-G6,48-DRA-COR-MUSE-8-H6,49-DRA-SO1-MUSE-16-A7,4-DRA-CH1-MUSE-4-D1,50-DRA-SO1-MUSE-15-B7,51-DRA-SO1-MUSE-13-C7,52-DRA-SO1-MUSE-12-D7,53-DRA-SO1-MUSE-10-E7,54-DRA-COR-MUSE-15-F7,55-DRA-COR-MUSE-14-G7,56-DRA-COR-MUSE-13-H7,57-DRA-COR-MUSE-12-A8,58-DRA-COR-MUSE-11-B8,59-DRA-COR-MUSE-10-C8,5-DRA-CH1-MUSE-5-E1,60-DRA-COR-MUSE-9-D8,61-DRA-MAL-MUSE-15-E8,62-DRA-MAL-MUSE-14-F8,63-DRA-MAL-MUSE-12-G8,64-DRA-MAL-MUSE-11-H8,65-DRA-MAL-MUSE-10-A9,66-DRA-MAL-MUSE-9-B9,67-DRA-BZ1-MUSE-1-C9,68-DRA-BZ1-MUSE-2-D9,69-DRA-BZ1-MUSE-3-E9,6-DRA-CH1-MUSE-6-F1,70-DRA-BZ1-MUSE-4-F9,71-DRA-BZ1-MUSE-5-G9,72-DRA-BZ1-MUSE-6-H9,73-DRA-BZ2-MUSE-1-A10,74-DRA-BZ2-MUSE-2-B10,75-DRA-BZ2-MUSE-3-C10,77-DRA-BZ2-MUSE-5-E10,78-DRA-BZ2-MUSE-6-F10,79-DRA-BZ2-MUSE-7-G10,7-DRA-CH1-MUSE-7-G1,80-DRA-CH1-MUSE-9-H10,81-DRA-CH1-MUSE-10-A11,82-DRA-CH1-MUSE-11-B11,83-DRA-CH1-MUSE-12-C11,84-DRA-CH1-MUSE-13-D11,85-DRA-CH1-MUSE-14-E11,86-DRA-CH1-MUSE-15-F11,87-DRA-CH1-MUSE-16-G11,88-DRA-CH1-MUSE-17-H11,89-DRA-CH1-MUSE-18-A12,8-DRA-CH1-MUSE-8-H1,90-DRA-CH1-MUSE-19-B12,91-DRA-CH1-MUSE-20-C12,9-DRA-CH2-MUSE-1-A2" \
-V ${out2}/populations.snps.filtered.recode_MIS_filt.vcf \
-F ${out2}/histogram_MAF_pops_ALL.pdf &> ${out2}/divide_populations_and_HWE.log &


# perform PCA analysis and IBD with SNPRelate
for step in pca_analysis ibd_analysis 
do
    cd ${out2}
    Rscript ${FUNC_DIR}/a04_SNPRelate_PCA_HO_KIN_IBD_stripchart.r \
    -p ${spp} \
    -V ${out2}/populations.snps.filtered.recode_MIS_filt.vcf \
    -P ${raw_files}/sample_ID_and_populations.txt \
    -O ${out2}/SNPRelate/ \
    -s ${step} &> ${out2}/${step}.log
done

# FST heatmap. We need to adapt the file that was generated with stacks populations.fst_summary.tsv for the heatmap.
# # example input table
# ISS  RMA       CRE       LSO       MOL       TON       RIC       PLA
# ISS  0.578954  0.56816   0.572382  0.380488  0.662919  0.590242  0.673136
# RMA  0.58424   0.442205  0.453758  0.615775  0.518718  0.626302
# CRE  0.546331  0.443563  0.700695  0.639801  0.72085
# LSO  0.43062   0.640166  0.571839  0.654934
# MOL  0.553193  0.509198  0.558448
# TON  0.593288  0.711866
# RIC  0.308493


# # example desired output table
# 0    ISS  RMA       CRE      LSO       MOL       TON       RIC       PLA
# ISS  0    0.578954  0.56816  0.572382  0.380488  0.662919  0.590242  0.673136
# RMA  0    0         0.58424  0.442205  0.453758  0.615775  0.518718  0.626302
# CRE  0    0         0        0.546331  0.443563  0.700695  0.639801  0.72085
# LSO  0    0         0        0         0.43062   0.640166  0.571839  0.654934
# MOL  0    0         0        0         0         0.553193  0.509198  0.558448
# TON  0    0         0        0         0         0         0.593288  0.711866
# RIC  0    0         0        0         0         0         0         0.308493
# PLA  0    0         0        0         0         0         0         0


# prepare input table
awk -F'\t' -v OFS='\t' '{ for (i = 1; i <= NF; ++i) sub(/^$/, 0, $i) } 1' ${out2}/populations.fst_summary.tsv | cat - <(echo -e "BZ2\t0\t0\t0\t0\t0\t0\t0\t0") > ${out2}/populations.fst_summary_for_heatmap.tsv

# FST heatmap
Rscript ${FUNC_DIR}/a05_FST_heatmap.r \
-I ${out2}/populations.fst_summary_for_heatmap.tsv \
-O ${out2}/ &> ${out2}/FST.log

# plots
# boxplot with expected and observed heterozygosity and histograms with mean exp and obs heterozygosity
Rscript ${FUNC_DIR}/a06_plots_heterozygosity.r \
-I ${out2}/populations.sumstats.tsv

# Count heterozygous SNPs + scatterplot
Rscript ${FUNC_DIR}/a07_count_het_SNPs_per_ind.r \
-V ${out2}/populations.snps.filtered.recode_MIS_filt_header.vcf \
-O ${out2}/genotype_counts.tbl

# Create phylogenetic trees
Rscript ${FUNC_DIR}/a08_phyl_tree_plot.r \
-V ${out2}/populations.snps.filtered.recode_MIS_filt_header.vcf \
-M ${raw_files}/sample_ID_and_populations.txt \
-O ${out2}/phyl_tree_NJ

# mantel test between genetic distances (measured by fst) and geographic coordinates
# geographic coords for each population as follows:
# Latitude   Longitude    population
# 40.032339  15.272951    POP1
# 40.032329  15.272896    POP2
# 40.032261  15.272726    POP3
# 40.032179  15.272605    POP4

Rscript ${FUNC_DIR}/a10_populations_mantel_test.r \
-F ${out2}/populations.fst_summary_for_heatmap.tsv \
-G ${out2}/pop_geographic_coordinates.txt \
-O ${out2}/

# count SNP distribution and transition/transversion from SNPs identified
awk 'BEGIN {OFS="\t"} !/^#/ {changes[$4"/"$5]++} END {for (c in changes) print c, changes[c]}' ${out2}/populations.snps.filtered.recode_MIS_filt_header.vcf | \
awk 'OFS="\t" {if($1=="C/T" || $1=="T/C" ||$1=="A/G" ||$1=="G/A")print $0,"transition"; else print $0,"transversion"}' | sort -k3,3 > ${out2}/SNP_type.txt
cat <(echo -e "type\tSNPs") <(awk 'OFS="\t" {if($NF=="transition") print $0}' ${out2}/SNP_type.txt | cut -f3 --complement) <(awk 'OFS="\t" {if($NF=="transition") print $0}' ${out2}/SNP_type.txt | awk '{ sum+=$2} END {print "transitions",sum}') \
<(awk 'OFS="\t" {if($NF=="transversion") print $0}' ${out2}/SNP_type.txt | cut -f3 --complement) \
<(awk 'OFS="\t" {if($NF=="transversion") print $0}' ${out2}/SNP_type.txt | awk '{ sum+=$2} END {print "transversions",sum}') > ${out2}/SNP_type_final.txt

# histogram with transition/transposition
Rscript ${FUNC_DIR}/a11_histogram_transition_transversion.r \
-I ${out2}/SNP_type_final.txt \
-O ${out2}/

# gene flow Nm
# Wright method for gene flow estimation
# Nm=(1-Fst)/4Fst
# Wright believed that if the Fst value of the populations was between 0 and 0.05, it will indicate that there is no genetic differentiation among populations, if Fst value is between 0.05 and 0.15, it is moderately differentiated, if Fst value is between 0.15 and 0.25, then it is highly differentiated
mkdir -p ${out1}/gene_flow_Nm

awk 'NR==1{print; next} {
  printf $1 "\t"
  for(i=2; i<=NF; i++){
    fst = $i
    if(fst == 0){
      nm = "NA"   # avoid division by zero
    } else {
      nm = (1 - fst) / (4 * fst)
    }
    printf nm (i<NF ? "\t" : "\n")
  }
}' ${out1}/populations.fst_summary_for_heatmap.tsv | sed 's/NA/0/g' > ${out1}/gene_flow_Nm/Nm_${spp}_sq_matrix.tbl

Rscript ${FUNC_DIR}/a12_gene_flow_Nm.r --infile ${out1}/gene_flow_Nm/Nm_${spp}_sq_matrix.tbl --outpath ${out1}/gene_flow_Nm/
