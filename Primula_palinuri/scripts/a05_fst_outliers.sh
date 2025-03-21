# Copyright:	    Paloma Perez & Fabio Marroni 2024
# Aim:              FST outlier analyses to identify potential adaptive loci
# To add:		     
# Suggestions: 
# Fixes:  

spp=primula_palinuri

# create output folders
out1=${INPUT_DIR}/${spp}/de_novo_pipeline
out2=${INPUT_DIR}/${spp}/de_novo_pipeline/populations_stacks
# folder containing fastq sequences
raw_files=${INPUT_DIR}/${spp}/raw_reads
# Folder containing the functions
FUNC_DIR=functions

# With every file including the comparisons. For each SNP between all comparisons do the media for the FST value. We will have Fst value for every SNP. First merge all tables using the column with SNP position
# Classify SNPs according to FST value. Those SNPs with higher Fst, we see if they are present in stacks with a lot of SNPs. Take a subset and do blastx

# first calculate SNP frequency spectrum (SFS)
# FILE WITH SFS WITH MORE INFO
mkdir -p ${out2}/fst_outlier_analysis
cd ${out2}/fst_outlier_analysis

for file in populations.fst_CIM-FIU populations.fst_CIM-SGP populations.fst_ID-SGP populations.fst_PPA-FIU populations.fst_SGP-LAM populations.fst_CIM-ID populations.fst_FIU-LAM populations.fst_PC-FIU  populations.fst_PPA-ID populations.fst_CIM-LAM populations.fst_FIU-SGP populations.fst_PC-ID populations.fst_PPA-LAM populations.fst_CIM-PC populations.fst_ID-FIU populations.fst_PC-LAM populations.fst_PPA-PC populations.fst_CIM-PPA populations.fst_ID-LAM populations.fst_PC-SGP populations.fst_PPA-SGP
do
    join -a 1 -e 'NA' -o 'auto' -1 1 -2 1 <(grep -v "#" ${out2}/populations.snps.filtered.recode_MIS_filt_header.vcf | awk 'OFS="\t" {print $1"_"$2}' | sort -k1,1) <(sed 1d ${out2}/${file}.tsv |awk 'OFS="\t" {print $1"_"$6-1,$0}'|  sort -k1,1 ) -t $'\t' | \
        # at the end I removed the header  but just to know the proper format
    cat <(echo -e "Locus_ID\tLocus\tPop_1_ID\tPop_2_ID\tChr\tBP\tColumn\tOverall_Pi\tAMOVA_Fst\tFisher's_P\tOdds_Ratio\tCI_Low\tCI_High\tLOD\tCorrected_AMOVA_Fst\tSmoothed_AMOVA_Fst\t     Smoothed_AMOVA_Fst_P-value\tWindow_SNP_Count") - | \
    awk 'OFS="\t" {print $3"-"$4,$9}' | sed 1d > ${out2}/fst_outlier_analysis/${file}_fst_popnames.tbl &
done

paste ${out2}/fst_outlier_analysis/*_fst_popnames.tbl | paste <(grep -v "#" ${out2}/populations.snps.filtered.recode_MIS_filt_header.vcf | awk 'OFS="\t" {print $1"_"$2}' | sort -k1,1) - | \
# remove lines with only NAs
awk '{ for(i=2; i<=NF; i++) if ($i != "NA" && $i != "NA-NA") { print; next } }' > ${out2}/fst_outlier_analysis/ALL_fst_popnames.tbl

# FILE WITH SFS SUMMARY
for file in populations.fst_CIM-FIU populations.fst_CIM-SGP populations.fst_ID-SGP populations.fst_PPA-FIU populations.fst_SGP-LAM populations.fst_CIM-ID populations.fst_FIU-LAM populations.fst_PC-FIU  populations.fst_PPA-ID populations.fst_CIM-LAM populations.fst_FIU-SGP populations.fst_PC-ID populations.fst_PPA-LAM populations.fst_CIM-PC populations.fst_ID-FIU populations.fst_PC-LAM populations.fst_PPA-PC populations.fst_CIM-PPA populations.fst_ID-LAM populations.fst_PC-SGP populations.fst_PPA-SGP
do
    join -a 1 -e 'NA' -o 'auto' -1 1 -2 1 <(grep -v "#" ${out2}/populations.snps.filtered.recode_MIS_filt_header.vcf | awk 'OFS="\t" {print $1"_"$2}' | sort -k1,1) <(sed 1d ${out2}/${file}.tsv |awk 'OFS="\t" {print $1"_"$6-1,$0}'|  sort -k1,1 ) -t $'\t' | \
    # at the end I removed the header  but just to know the proper format
    cat <(echo -e "Locus_ID\tLocus\tPop_1_ID\tPop_2_ID\tChr\tBP\tColumn\tOverall_Pi\tAMOVA_Fst\tFisher's_P\tOdds_Ratio\tCI_Low\tCI_High\tLOD\tCorrected_AMOVA_Fst\tSmoothed_AMOVA_Fst\t     Smoothed_AMOVA_Fst_P-value\tWindow_SNP_Count") - | \
    awk 'OFS="\t" {print $3"-"$4,$9}' | sed 1d | cut -f2 > ${out2}/fst_outlier_analysis/${file}_fst.tbl &
done

paste ${out2}/fst_outlier_analysis/*_fst.tbl | paste <(grep -v "#" ${out2}/populations.snps.filtered.recode_MIS_filt_header.vcf | awk 'OFS="\t" {print $1"_"$2}' | sort -k1,1) - | \
# remove lines with only NAs
awk '{ for(i=2; i<=NF; i++) if ($i != "NA") { print; next } }' > ${out2}/fst_outlier_analysis/ALLfst.tbl

# header
for file in populations.fst_CIM-FIU populations.fst_CIM-SGP populations.fst_ID-SGP populations.fst_PPA-FIU populations.fst_SGP-LAM populations.fst_CIM-ID populations.fst_FIU-LAM populations.fst_PC-FIU  populations.fst_PPA-ID populations.fst_CIM-LAM populations.fst_FIU-SGP populations.fst_PC-ID populations.fst_PPA-LAM populations.fst_CIM-PC populations.fst_ID-FIU populations.fst_PC-LAM populations.fst_PPA-PC populations.fst_CIM-PPA populations.fst_ID-LAM populations.fst_PC-SGP populations.fst_PPA-SGP
do
    join -a 1 -e 'NA' -o 'auto' -1 1 -2 1 <(grep -v "#" ${out2}/populations.snps.filtered.recode_MIS_filt_header.vcf | awk 'OFS="\t" {print $1"_"$2}' | sort -k1,1) <(sed 1d ${out2}/${file}.tsv |awk 'OFS="\t" {print $1"_"$6-1,$0}'|  sort -k1,1 ) -t $'\t' | \
    # at the end I removed the header  but just to know the proper format
    cat <(echo -e "Locus_ID\tLocus\tPop_1_ID\tPop_2_ID\tChr\tBP\tColumn\tOverall_Pi\tAMOVA_Fst\tFisher's_P\tOdds_Ratio\tCI_Low\tCI_High\tLOD\tCorrected_AMOVA_Fst\tSmoothed_AMOVA_Fst\t     Smoothed_AMOVA_Fst_P-value\tWindow_SNP_Count") - | \
    awk 'OFS="\t" {print $3"-"$4,$9}' | sed 1d | cut -f1 > ${out2}/fst_outlier_analysis/${file}_fst_header.tbl &
done

paste ${out2}/fst_outlier_analysis/*_fst_header.tbl | sort -u | grep -v "NA" > ${out2}/fst_outlier_analysis/fst_header.txt
rm *_fst_header.tbl


paste <( cut -f1 ${out2}/fst_outlier_analysis/ALL_fst_popnames.tbl) ${out2}/fst_outlier_analysis/ALLfst.tbl | \
cat <(paste <(echo -e "Chr|pos") ${out2}/fst_outlier_analysis/fst_header.txt) - > ${out2}/fst_outlier_analysis/ALLfst_pops.tbl

# mean fst per snp across all comparisons
paste <(awk 'OFS="\t" {out=$0; for(i=2;i<=NF;i++) t+=$i; print out,t; t=0}' ${out2}/fst_outlier_analysis/ALLfst_pops.tbl) \
<(awk -F'\t' 'OFS="\t"{count=0; for(i=2; i<=NF; i++) if($i!="NA") count++; print count}' ${out2}/fst_outlier_analysis/ALLfst_pops.tbl) | \
awk 'OFS="\t" {print $0,$(NF-1)/$NF}' | \
sed 1d | \
sort -nr -k25,25 | \
cat <(paste <(head -1 ${out2}/fst_outlier_analysis/ALLfst_pops.tbl) <(echo -e "sum_FST\tcount_FST\tmean_FST")) - > ${out2}/fst_outlier_analysis/ALLfst_FINAL.tbl

# Classify SNPs according to FST value. 5%, 1% and 0.1% top SNPs with higher Fst
# top 5%, 1%, 0.1%
conda deactivate 
conda activate /iga/scripts/dev_modules/miniconda3/envs/seedforce
for top in 5 1 0.1
do
    cd ${out2}/fst_outlier_analysis/
    # calculate 0.1, 1 and 5% of total Fst outliers
    var=$(awk -v top=${top} 'OFS="\t" {print int(NR*top/100)}' ${out2}/fst_outlier_analysis/stacks_fst_sorted.txt | tail -1)
    head -${var} ${out2}/fst_outlier_analysis/stacks_fst_sorted.txt | xargs samtools faidx ${out2}/../catalog.fa > ${out2}/fst_outlier_analysis/stacks_fst_top${top}perc.fa
done
# ----

# Classify SNPs according to FST value. Those SNPs with higher Fst, we see if they are present in stacks with a lot of SNPs. Take a subset and do blast (check top 20,30 ,50, 100 etc) blastx only

# get fasta sequences for interesting stacks
cut -f1 ${out2}/fst_outlier_analysis/ALLfst_FINAL.tbl | sed 1d | \
sed 's/_/\t/g' | cut -f1 | uniq > ${out2}/fst_outlier_analysis/stacks_fst_sorted.txt

conda deactivate 
conda activate /iga/scripts/dev_modules/miniconda3/envs/seedforce

samtools faidx ${out2}/../catalog.fa 
cd ${out2}/fst_outlier_analysis/
cat ${out2}/fst_outlier_analysis/stacks_fst_sorted.txt | xargs samtools faidx ${out2}/../catalog.fa > ${out2}/fst_outlier_analysis/stacks_fst.fa

conda deactivate
conda activate /projects/novabreed/share/anaconda-3/envs/blast_2.10.1
# blastx All snps
blastx -db /projects/marroni/seedforce/linaria_flava/de_novo_pipeline/blastx_stacks/database/swissprot -query ${out2}/fst_outlier_analysis/stacks_fst.fa -num_threads 20 -outfmt 6 -evalue 1E-5 -out ${out2}/fst_outlier_analysis/blastx_nr_ALLfst.tab -num_alignments 5


# get accession descriptions

cut -f2 ${out2}/fst_outlier_analysis/blastx_nr_ALLfst.tab > ${out2}/fst_outlier_analysis/accessions.txt




