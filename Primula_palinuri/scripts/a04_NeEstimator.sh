# Copyright:	    Paloma Perez & Fabio Marroni 2025
# Aim:              calculate Effective population size using NeEstimator
# To add:		     
# Suggestions: 
# Fixes:  

spp=primula_palinuri

out1=${INPUT_DIR}/${spp}/de_novo_pipeline
out2=${INPUT_DIR}/${spp}/de_novo_pipeline/populations_stacks
# folder containing fastq sequences
raw_files=${INPUT_DIR}/${spp}/raw_reads
# Folder containing the functions
FUNC_DIR=functions

# https://neestimator.software.informer.com/download/

##################
# POPS SEPARATED #
##################
mkdir -p ${out2}/Ne_effective_pop_size/NeEstimator_V2/pops_sep
path=${out2}/Ne_effective_pop_size/NeEstimator_V2/

# get one vcf for each population separated
for pop in LAM FIU ID PC SGP CIM
do
    cd ${path}/pops_sep
    sed 1d /projects/marroni/seedforce/primula_palinuri/raw_reads/sample_ID_and_populations.txt | cut -f1 | grep ${pop} > ${path}/pops_sep/${pop}.txt
    vcftools --vcf ${path}/populations.snps.filtered.recode_MIS_filt_PLINK_header.vcf --keep ${path}/pops_sep/${pop}.txt --recode --out filtered_${pop}
    grep -v "\./\.:" ${path}/pops_sep/filtered_${pop}.recode.vcf > ${path}/pops_sep/filtered_${pop}.recode_inALLinds.vcf
done

# PPA pop is different because is also the spp name
pop=PPA
cd ${path}/pops_sep
paste <(sed 1d /projects/marroni/seedforce/primula_palinuri/raw_reads/sample_ID_and_populations.txt | cut -f1 | sed 's/-/\t/g' | cut -f3) <(sed 1d /projects/marroni/seedforce/primula_palinuri/raw_reads/sample_ID_and_populations.txt | cut -f1) | awk -v pop=${pop} 'OFS="\t" {if($1==pop) print $2}' > ${path}/pops_sep/${pop}.txt
vcftools --vcf ${path}/populations.snps.filtered.recode_MIS_filt_PLINK_header.vcf --keep ${path}/pops_sep/${pop}.txt --recode --out filtered_${pop}
grep -v "\./\.:" ${path}/pops_sep/filtered_${pop}.recode.vcf > ${path}/pops_sep/filtered_${pop}.recode_inALLinds.vcf


# change ind ID for pop ID
cat <(grep "##" ${path}/populations.snps.filtered.recode_MIS_filt_PLINK_header.vcf ) <(grep -v "##" ${path}/populations.snps.filtered.recode_MIS_filt_PLINK_header.vcf | head -1 | tr '\t' '\n' | sed -E 's/^[0-9]+-PPA-([^-]+)(-[^-]+)?(-[0-9]+-[A-Z0-9]+)?/\1/'| awk '{printf "%s\t", $0} END {print ""}') <(grep -v "#" ${path}/populations.snps.filtered.recode_MIS_filt_PLINK_header.vcf) > ${path}/pops_sep/populations.snps.filtered.recode_MIS_filt_PLINK_header_names.vcf

# change vcf format to GENPOP neccessary for NeEstimator
Rscript ${FUNC_DIR}/a14_NeEstimator.r

# Run NeEstimator for LD random model
conda deactivate
conda activate /iga/scripts/dev_modules/miniconda3/envs/java
cd /iga/scripts/dev_modules/miniconda3/envs/java
java -jar /projects/marroni/seedforce/softwares/NeEstimator/test/NeEstimator2x1.jar
