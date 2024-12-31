#!/bin/bash

mkdir -p /projects/marroni/seedforce/pool/de_novo_pipeline_catalog_spp_pops
out1=/projects/marroni/seedforce/pool/de_novo_pipeline_catalog_spp_pops
raw_files=/projects/marroni/seedforce/pool/raw_reads

t=20
id=1

# Twenty-one species (Astragalus verrucosus (AVE), Bassia Saxicola (BSA), Campanula sabatia (CSA), Cytisus aeolicus )(CAE), Limonium strictissimum (LIS), Linum muelleri (Failed due to low number of reads), Ribes sardoum (RSA), Silene hicesiae (SHI), Adenophora liliifolia (ALI), Centranthus amazonum (Failed due to low number of raw reads), Crepis pusilla (CPU), Eleocharis carniolica (ECA), Gentiana ligustica (GLI), Gladiolus palustris (GPA), Himantoglossum adriaticum (HAD), Kosteletzkya pentacarpos (KOS), Leucojum nicaeense (ACI), Linaria pseudolaxiflora (LPS), Marsilea quadrifolia (MAQ), Woodwardia radicans (WOR) and Gladiolus illyricus (GIL), this last one not included in the original plan) were sequenced following the procedure for the pool approach (table 1) and DNA was extracted from 613 individuals that were included in 77 populations. 
# Spp are listed in ${raw_files}/population_files_per_spp/spp.txt. Example of header below:
# ALI
# AVE
# BSA

# required modules
module load it/rad/stacks/2.4
module load it/aligners/bwa/0.7.17
module load it/tools/samtools/1.7
module load sw/bio/picard-tools/1.88
/projects/marroni/software/freebayes/freebayes-1.3.6
module load lang/java/jdk1.7.0_07
module load lang/r/2.15.1
# create catalog for each spp and combine with all pops per spp together

while read spp
do
    mkdir -p ${out1}/${spp}
    echo -e "I started de_novo_pipeline stacks $spp\n"
    # skip rows starting with #
    [[ "$spp" =~ \#.* ]]&& continue
    module purge
    module load it/rad/stacks/2.4
    perl /iga/scripts/packages/stacks-2.4/bin/denovo_map.pl -M 2 -T ${t} -o ${out1}/${spp} --popmap ${raw_files}/population_files_per_spp/population_${spp}.txt --samples ${raw_files} --paired &> ${out1}/${spp}/de_novo_pipeline.log
    echo -e "I finished de_novo_pipeline stacks $spp\n"
    # create index for the "reference" catalog
    echo -e "creating index for catalog $spp\n"
    module load it/aligners/bwa/0.7.17
    bwa index ${out1}/${spp}/catalog.fa.gz
done < ${raw_files}/population_files_per_spp/spp.txt 

while read spp
do
    # # skip rows starting with #
    [[ "$spp" =~ \#.* ]]&& continue
        while read spp pop ID
        do
        # skip rows starting with #
        [[ "$pop" =~ \#.* ]]&& continue
        # map with BWA mem using the stacks catalog
        echo -e "I started BWA mem $spp ${pop}\n"
        module purge
        module load it/aligners/bwa/0.7.17
        bwa mem ${out1}/${spp}/catalog.fa.gz ${raw_files}/*${spp}*${pop}*.1.fq.gz ${raw_files}/*${spp}*${pop}*.2.fq.gz -t ${t} > ${out1}/${spp}/${spp}_${pop}.sam
        echo -e "I finished bwa mem $spp ${pop}\n"
        module load it/tools/samtools/1.7
        samtools view -b ${out1}/${spp}/${spp}_${pop}.sam -o ${out1}/${spp}/${spp}_${pop}.bam
        rm ${out1}/${spp}/${spp}_${pop}.sam
        #Sort BAM files by coordinate
        echo -e "I started samtools sort $spp $pop $ID\n"
        samtools sort -O BAM -o ${out1}/${spp}/${spp}_${pop}_sorted.bam ${out1}/${spp}/${spp}_${pop}.bam -@ ${t}
        echo -e "I finished samtools sort $spp $pop $ID\n"
        done < /projects/marroni/seedforce/pool/raw_reads/population_files_per_spp/${spp}/population_spp_pops_IDs.txt
done < ${raw_files}/population_files_per_spp/spp.txt

# ALL pops together merged in a unique file        
while read spp
do
    # # skip rows starting with #
    [[ "$spp" =~ \#.* ]]&& continue       
    samtools merge -f ${out1}/${spp}/${spp}_sorted_merged.bam ${out1}/${spp}/${spp}_*_sorted.bam -@ ${t}        
    echo -e "I started samtools quality filter mapping quality=5 $spp $pop $ID\n"
    samtools view -bq 5 -o ${out1}/${spp}/${spp}_sorted_merged_q5.bam ${out1}/${spp}/${spp}_sorted_merged.bam -@ ${t}
    echo -e "I finished samtools quality filter mapping quality=5 $spp $pop $ID\n"
    #Add read group with picardtools
    module load sw/bio/picard-tools/1.88
    tmp_dir=/projects/novabreed/share/epaparelli/tmp
    PICARD_PATH=/iga/stratocluster/packages/sw/bio/picard-tools/1.88
    BASEDIR=${out1}/${spp}/
    ALDIR=${out1}/${spp}/
    for INFILE in ${ALDIR}/*_sorted_merged_q5.bam
    do
    SAMPLE=$(basename $INFILE)
    SAMPLE=${SAMPLE/_HQ.bam/}
    # THIS STEP IS IMPORTANT FOR CHANGING THE NAMES AND MAKE THE FILES SUITABLE FOR SNP CALLING
    # better to asign the sample name as this MYRGPU
    MYRGPU=$(samtools view $INFILE | head -n1 | cut -d":" -f3)
    MYRGPU=$SAMPLE
    # create random ID. IS IMPORTANT BECAUSE WE WANT TO HAVE A DIFFERENT ID FOR EACH FILE.
    RANDOMID=$(head -5 /dev/urandom | tr -cd '[:alnum:]' | cut -c -10)
    RGLINE="RGID=$RANDOMID RGPL="ILLUMINA" RGCN="Istituto_di_Genomica_Applicata" RGPU=$MYRGPU_name RGLB=${SAMPLE}_library RGSM=${SAMPLE}"
    #RGDS=nome/numero della corsa+indice (specifico per corsa e campione)   RGLB=nome libreria  RGSM=nome campione
    module load lang/java/jdk1.7.0_07
    module load lang/r/2.15.1
    java -Xmx4g -Djava.io.tmpdir=${tmp_dir} \
    -jar ${PICARD_PATH}/AddOrReplaceReadGroups.jar INPUT=${INFILE} \
    OUTPUT=${ALDIR}/${SAMPLE}_RG.bam ${RGLINE} RGLB=${SAMPLE}\
    VALIDATION_STRINGENCY=LENIENT TMP_DIR=${tmp_dir}
    done
    #prepare catalog to run freebayes to call SNPs
    echo -e "preparing catalog for freebayes to call SNPs $spp\n"
    gunzip ${out1}/${spp}/catalog.fa.gz
    GENOME=${out1}/${spp}/catalog.fa
    samtools faidx ${out1}/${spp}/catalog.fa
    #We need freebayes to call SNPs on "polyploids"
    #First attempt: we use continuous pool model (we pretend we do not know pool size).
    echo -e "I started freebayes $spp\n"
    GENOME=${out1}/${spp}/catalog.fa
    ALDIR=${out1}/${spp}
    /projects/marroni/software/freebayes/freebayes-1.3.6 -f ${GENOME} \
    --use-best-n-alleles 4 --pooled-continuous -F 0.01 -C 2 --report-genotype-likelihood-max --min-mapping-quality 5 --min-base-quality 20 --min-coverage 40 \
    ${ALDIR}/${spp}_sorted_merged_q5.bam_RG.bam > ${ALDIR}/${spp}.vcf
    echo -e "I finished freebayes $spp\n"        
    # filter indels and prepare vcf file
    echo -e "I started filtering indels and preparing vcf file $spp\n"
    grep -v "##" ${out1}/${spp}/${spp}.vcf | awk 'OFS="\t" {print $1,$2,$4,$5,$NF}' | awk 'BEGIN { FS = OFS = "\t" } { gsub(",", "\t", $NF) } 1' | \
    sed 's/:/\t/g' | cut -f1,2,3,4,7,8 | sed 1d | \
    # remove indels in REF or ALT column
    awk -F'\t' 'length($3) <= 1 && length($4) <= 1' | \
    # minor allele column
    awk 'OFS="\t" {if($(NF-1)<$NF)print $0,$(NF-1); else print $0,$NF}' | \
    cat <(echo -e "#CHROM\tPOS\tREF\tALT\tP_REF\tP_ALT\tMin_all") - > ${out1}/${spp}/${spp}_FINAL.tbl
    # # MAF plot and exp HET is done in the following comands
    echo -e "I finished POOL de novo pipeline $spp\n"
done < ${raw_files}/population_files_per_spp/spp.txt

# get SNPs per each pop separately
while read spp
do
    # skip rows starting with #
    [[ "$spp" =~ \#.* ]]&& continue
    while read spp pop ID
    do
    t=10
    # skip rows starting with #
    [[ "$pop" =~ \#.* ]]&& continue
    samtools view -bq 5 -o ${out1}/${spp}/${spp}_${pop}_sorted_q5.bam ${out1}/${spp}/${spp}_${pop}_sorted.bam -@ ${t}
    sleep 2
    echo -e "I finished samtools quality filter for 5 reads $spp $pop $ID\n"
    #Add read group with picardtools
    module load sw/bio/picard-tools/1.88
    tmp_dir=/projects/novabreed/share/epaparelli/tmp
    PICARD_PATH=/iga/stratocluster/packages/sw/bio/picard-tools/1.88
    BASEDIR=${out1}/${spp}/
    ALDIR=${out1}/${spp}/
    for INFILE in ${ALDIR}/*_sorted_q5.bam
    do
    SAMPLE=$(basename $INFILE)
    SAMPLE=${SAMPLE/_HQ.bam/}
    # THIS STEP IS IMPORTANT FOR CHANGING THE NAMES AND MAKE THE FILES SUITABLE FOR SNP CALLING
    # better to asign the sample name as this MYRGPU
    MYRGPU=$(samtools view $INFILE | head -n1 | cut -d":" -f3)
    MYRGPU=$SAMPLE
    # create random ID. IS IMPORTANT BECAUSE WE WANT TO HAVE A DIFFERENT ID FOR EACH FILE.
    RANDOMID=$(head -5 /dev/urandom | tr -cd '[:alnum:]' | cut -c -10)
    RGLINE="RGID=$RANDOMID RGPL="ILLUMINA" RGCN="Istituto_di_Genomica_Applicata" RGPU=$MYRGPU_name RGLB=${SAMPLE}_library RGSM=${SAMPLE}"
    #RGDS=nome/numero della corsa+indice (specifico per corsa e campione)   RGLB=nome libreria  RGSM=nome campione
    module load lang/java/jdk1.7.0_07
    module load lang/r/2.15.1
    java -Xmx4g -Djava.io.tmpdir=${tmp_dir} \
    -jar ${PICARD_PATH}/AddOrReplaceReadGroups.jar INPUT=${INFILE} \
    OUTPUT=${ALDIR}/${SAMPLE}_RG.bam ${RGLINE} RGLB=${SAMPLE}\
    VALIDATION_STRINGENCY=LENIENT TMP_DIR=${tmp_dir}
    done
    #We need freebayes to call SNPs on "polyploids"
    #First attempt: we use continuous pool model (we pretend we do not know pool size).
    echo -e "I started freebayes $spp ${pop}\n"
    GENOME=${out1}/${spp}/catalog.fa
    ALDIR=${out1}/${spp}
    /projects/marroni/software/freebayes/freebayes-1.3.6 -f ${GENOME} \
    --use-best-n-alleles 4 --pooled-continuous -F 0.01 -C 2 --report-genotype-likelihood-max --min-mapping-quality 5 --min-base-quality 20 --min-coverage 40 \
    ${ALDIR}/${spp}_${pop}_sorted_q5.bam > ${ALDIR}/${spp}_${pop}.vcf &
    # how can i create a log file with freebayes? the following comand is not working
    # /projects/marroni/software/freebayes/freebayes-1.3.6 -f ${GENOME} \
    # --use-best-n-alleles 4 --pooled-continuous -F 0.01 -C 2 --report-genotype-likelihood-max --min-mapping-quality 5 --min-base-quality 20 --min-coverage 40 \
    # ${ALDIR}/${spp}_${pop}_sorted_q5.bam > ${ALDIR}/${spp}_${pop}.vcf &> ${ALDIR}/freebayes_${pop}.log &
    done < /projects/marroni/seedforce/pool/raw_reads/population_files_per_spp/${spp}/population_spp_pops_IDs.txt
done < ${raw_files}/population_files_per_spp/spp.txt    

while read spp
do
    # skip rows starting with #
    [[ "$spp" =~ \#.* ]]&& continue
    while read spp pop ID
    do
    sleep 1
    # skip rows starting with #
    [[ "$pop" =~ \#.* ]]&& continue   
    echo -e "I finished freebayes $spp ${pop}\n"        
    # filter indels and prepare vcf file
    echo -e "I started filtering indels and preparing vcf file $spp\n"
    grep -v "##" ${out1}/${spp}/${spp}_${pop}.vcf | awk 'OFS="\t" {print $1,$2,$4,$5,$NF}' | awk 'BEGIN { FS = OFS = "\t" } { gsub(",", "\t", $NF) } 1' | \
    sed 's/:/\t/g' | cut -f1,2,3,4,7,8 | sed 1d | \
    # remove indels in REF or ALT column
    awk -F'\t' 'length($3) <= 1 && length($4) <= 1' | \
    # minor allele column
    awk 'OFS="\t" {if($(NF-1)<$NF)print $0,$(NF-1); else print $0,$NF}' | \
    cat <(echo -e "#CHROM\tPOS\tREF\tALT\tP_REF\tP_ALT\tMin_all") - > ${out1}/${spp}/${spp}_${pop}_FINAL.tbl
    done < /projects/marroni/seedforce/pool/raw_reads/population_files_per_spp/${spp}/population_spp_pops_IDs.txt &
done < ${raw_files}/population_files_per_spp/spp.txt

while read spp
do
# skip rows starting with #
[[ "$spp" =~ \#.* ]]&& continue
    while read spp pop ID
    do
    # skip rows starting with #
    [[ "$pop" =~ \#.* ]]&& continue 
    # MAF plot and exp HET
    echo -e "I started Rscript $spp $pop\n"
    Rscript /projects/marroni/seedforce/scripts/HWE_analysis_MAF_plot.r ${spp} ${pop}
    echo -e "I finished Rscript $spp $pop\n"
    echo -e "I finished POOL MAF plot and He $spp ${pop}\n"
    done < ${raw_files}/population_files_per_spp/${spp}/population_spp_pops_IDs.txt
done < ${raw_files}/population_files_per_spp/spp.txt  

# calculate mean coverage
rm ${out1}/mean_coverage_all.txt
while read spp
do
    while read spp pop ID
    do
    samtools depth -aa ${out1}/${spp}/${spp}_${pop}_sorted.bam  | awk -v spp=${spp} -v pop=${pop} '{sum+=$3} END { print "Mean_coverage "spp,pop,sum/NR}' 
    done < ${raw_files}/population_files_per_spp/${spp}/population_spp_pops_IDs.txt
done < ${raw_files}/population_files_per_spp/spp.txt >> ${out1}/mean_coverage_all.txt &

# The output for samtools depth is a simple tab-separated table with three columns: reference name, position, and coverage depth

# calculate total mapped reads per pop
rm ${out1}/mapped_reads_all.txt
while read spp
do
    while read spp pop ID
    do
        paste <(echo -e "mapped_reads\t${spp}\t${pop}") <(samtools view -c -F 260 ${out1}/${spp}/${spp}_${pop}_sorted.bam)
    done < ${raw_files}/population_files_per_spp/${spp}/population_spp_pops_IDs.txt
done < ${raw_files}/population_files_per_spp/spp.txt >> ${out1}/mapped_reads_all.txt

# calculate total stacks per catalog per spp
rm ${out1}/catalog_all.txt
while read spp
do
    paste <(echo -e "catalog\t${spp}\t") <(grep ">" ${out1}/${spp}/catalog.fa | wc -l)
done < ${raw_files}/population_files_per_spp/spp.txt >> ${out1}/catalog_all.txt

# total SNPS with indels per pop
rm ${out1}/SNPS_with_indels_per_pop.txt
while read spp
do
    while read spp pop ID
    do
        paste <(echo -e "SNPS_with_indels\t${spp}\t${pop}") <(grep -v "##" ${out1}/${spp}/${spp}_${pop}.vcf | sed 1d | wc -l)
    done < ${raw_files}/population_files_per_spp/${spp}/population_spp_pops_IDs.txt
done < ${raw_files}/population_files_per_spp/spp.txt >> ${out1}/SNPS_with_indels_per_pop.txt

# total SNPS per pop
rm ${out1}/SNPS_per_pop.txt
while read spp
do
    while read spp pop ID
    do
        paste <(echo -e "SNPS\t${spp}\t${pop}") <(grep -v "##" ${out1}/${spp}/${spp}_${pop}_FINAL.tbl | sed 1d | wc -l)
    done < ${raw_files}/population_files_per_spp/${spp}/population_spp_pops_IDs.txt
done < ${raw_files}/population_files_per_spp/spp.txt >> ${out1}/SNPS_per_pop.txt

# total SNPS with indels per spp
rm ${out1}/SNPS_with_indels_per_spp.txt
while read spp
do
    paste <(echo -e "SNPS_with_indels\t${spp}\t") <(grep -v "##" ${out1}/${spp}/${spp}.vcf | sed 1d | wc -l)
done < ${raw_files}/population_files_per_spp/spp.txt >> ${out1}/SNPS_with_indels_per_spp.txt


# total SNPS per spp
rm ${out1}/SNPS_per_spp.txt
while read spp
do
    paste <(echo -e "SNPS\t${spp}\t") <(grep -v "##" ${out1}/${spp}/${spp}_FINAL.tbl | sed 1d | wc -l)
done < ${raw_files}/population_files_per_spp/spp.txt >> ${out1}/SNPS_per_spp.txt


# calculate total raw reads per pop
 
rm ${raw_files}/raw_reads_per_pop.txt
while read spp
do
    while read spp pop ID
    do
    zcat ${raw_files}/${ID}.1.fq.gz | grep "@" | wc -l | paste <(echo -e "${pop}") -
    done < ${raw_files}/population_files_per_spp/${spp}/population_spp_pops_IDs.txt
done < ${raw_files}/population_files_per_spp/spp.txt >> ${raw_files}/raw_reads_per_pop.txt


# Rscript to create PCA file
while read spp
do
    sleep 1
    Rscript /projects/marroni/seedforce/scripts/distance_pool.r -s ${spp} -I ${out1}/${spp}/ -o ${out1}/${spp}/PCA.pdf -O ${out1}/${spp}/ -S distance &
done < ${raw_files}/population_files_per_spp/spp.txt

# correct names and do PCA plot
while read spp
do
    sed -e 's/.vcf//g' ${out1}/${spp}/PCA.all_pops.${spp}.txt | sed -e 's/ALT_//g' | awk 'BEGIN{FS=OFS="\t"}{NF--; print}'  > ${out1}/${spp}/PCA.all_pops.${spp}_names.txt
done < ${raw_files}/population_files_per_spp/spp.txt

while read spp
do
    sleep 1
    Rscript /projects/marroni/seedforce/scripts/distance_pool.r -s ${spp} -I ${out1}/${spp}/ -o ${out1}/${spp}/PCA.pdf -O ${out1}/${spp}/ -S plot &
done < ${raw_files}/population_files_per_spp/spp.txt

# expected heterozygosity and standard error
rm ${out1}/exp_heterozygosity_and_se_ALL.tbl
while read spp
do
    sed 1d ${out1}/${spp}/exp_heterozygosity_and_se_per_pop.tbl | sed 's/\./,/g' >> ${out1}/exp_heterozygosity_and_se_ALL.tbl
done < ${raw_files}/population_files_per_spp/spp.txt

# add header
cat <(head -1 ${out1}/ACI/exp_heterozygosity_and_se_per_pop.tbl) ${out1}/exp_heterozygosity_and_se_ALL.tbl > ${out1}/exp_heterozygosity_and_se_ALL_final.tbl

# expected heterozygosity and standard error filtered for at least 5 or 10 reads MAF
while read spp
do
    for reads in 5 10
    do
    rm ${out1}/exp_heterozygosity_and_se_ALL_minorALL${reads}reads.tbl
    sed 1d ${out1}/${spp}/exp_heterozygosity_per_pop_minorALL${reads}reads.tbl | sed 's/\./,/g' >> ${out1}/exp_heterozygosity_and_se_ALL_minorALL${reads}reads.tbl
    done
done < ${raw_files}/population_files_per_spp/spp.txt

# expected heterozygosity per pool. combine all pops inside a spp
while read spp
do
    while read spp pop ID
    do
    sed 1d ${out1}/${spp}/${spp}_${pop}_allfreq_HE.tbl | awk -v pop=${pop} 'OFS="\t" {print $0,pop}' > ${out1}/${spp}/${spp}_${pop}_allfreq_HE_name.tbl
    done < ${raw_files}/population_files_per_spp/${spp}/population_spp_pops_IDs.txt
done < ${raw_files}/population_files_per_spp/spp.txt 

# filter by 5 or 10 reads in Minor allele
while read spp
do
    cat ${out1}/${spp}/${spp}_*_allfreq_HE_name.tbl | cat <(paste <(head -1 ${out1}/ACI/ACI_VIN_allfreq_HE.tbl) <(echo -e "POP_NAME")) - > ${out1}/${spp}/${spp}_allfreq_HE_name.tbl
    # filter by MAF at least 5 reads or 10 reads
    awk 'OFS="\t" {if($7>=5) print $0}' ${out1}/${spp}/${spp}_allfreq_HE_name.tbl > ${out1}/${spp}/${spp}_allfreq_HE_name_minorALL5reads.tbl
    awk 'OFS="\t" {if($7>=10) print $0}' ${out1}/${spp}/${spp}_allfreq_HE_name.tbl > ${out1}/${spp}/${spp}_allfreq_HE_name_minorALL10reads.tbl
done < ${raw_files}/population_files_per_spp/spp.txt 

# expected heterozygosity (He) for each pop. all SNPs and filtered for 5 and 10 reads in MAF
Rscript /projects/marroni/seedforce/scripts/plots_HE.r




