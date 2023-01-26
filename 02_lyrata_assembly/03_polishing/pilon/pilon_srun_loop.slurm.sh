#!/bin/bash

set -o errexit # exit on errors

module purge
module load SAMtools/1.9-GCC-8.2.0-2.31.1
module load Bowtie2/2.3.5.1-GCC-8.2.0-2.31.1
module load Pilon/1.23-Java-11
module list

THREADS="$SLURM_CPUS_PER_TASK"
ASSEMBLY=$(ls data/*.fasta)     # name of the (link to the) fasta
R1=$(ls data/*_R1_*.fastq.gz)
R2=$(ls data/*_R2_*.fastq.gz)
PREFIX=$(basename "$ASSEMBLY" .fasta)

PILON="java -Xmx16G -jar $EBROOTPILON/pilon.jar"
echo "Pilon command is: ${PILON}"
# Here are some interesting pilon options
OPTIONS="--diploid"  # assume diploid
OPTIONS="--vcf"       # output a VCF
OPTIONS="--variant"   # heuristic for variants not assembly
OPTIONS="--fix all"   # correct bases, indels, and local misassembly and write a FASTA file
OPTIONS="--changes"
OPTIONS="--threads 14"   # show changes to the FASTA
# Use this set of options to generate the homozygous genome informed by reads.
OPTIONS="--fix all --changes" 
echo "Pilon options set to: ${OPTIONS}"
    
run_bowtie () {
    echo START BOWTIE
    rm -rf bowtie_round${ROUND}
    mkdir bowtie_round${ROUND}
    cd bowtie_round${ROUND}
    echo "bowtie2-build --threads ${THREADS} ../${ASSEMBLY} ${PREFIX}_round${ROUND}"
    bowtie2-build --threads ${THREADS} ../${ASSEMBLY} ${PREFIX}_round${ROUND}
    echo "bowtie2 --threads ${THREADS} --no-unal --no-mixed --no-discordant --sensitive --end-to-end --score-min 'C,0,-1' -x ${PREFIX}_round${ROUND} -1 ../${R1} -2 ../${R2} -S ${PREFIX}_round${ROUND}.sam"
    bowtie2 --threads ${THREADS} --no-unal --no-mixed --no-discordant --sensitive --end-to-end --score-min 'C,0,-1' -x ${PREFIX}_round${ROUND} -1 ../${R1} -2 ../${R2} -S ${PREFIX}_round${ROUND}.sam
    echo "samtools sort -@ ${THREADS} -T tmp -O BAM -o ${PREFIX}_round${ROUND}.sorted.bam ${PREFIX}_round${ROUND}.sam"
    samtools sort -@ ${THREADS} -T tmp -O BAM -o ${PREFIX}_round${ROUND}.sorted.bam ${PREFIX}_round${ROUND}.sam
    echo "samtools index ${PREFIX}_round${ROUND}.sorted.bam"
    samtools index ${PREFIX}_round${ROUND}.sorted.bam
    echo "samtools flagstat ${PREFIX}_round${ROUND}.sorted.bam > ${PREFIX}_round${ROUND}.sorted.bam.flagstat"
    samtools flagstat ${PREFIX}_round${ROUND}.sorted.bam > ${PREFIX}_round${ROUND}.sorted.bam.flagstat
    echo "samtools stats -@ ${THREADS} ${PREFIX}_round${ROUND}.sorted.bam | grep '^SN' | cut -f 2- > ${PREFIX}_round${ROUND}.sorted.bam.stats.SN"
    samtools stats -@ ${THREADS} ${PREFIX}_round${ROUND}.sorted.bam | grep '^SN' | cut -f 2- > ${PREFIX}_round${ROUND}.sorted.bam.stats.SN
    echo DONE BOWTIE
    date
    cd ..
}

run_pilon() {
    echo START PILON
    rm -rf pilon_round${ROUND}
    mkdir pilon_round${ROUND}
    cd pilon_round${ROUND}
    BAMFILE="../bowtie_round${ROUND}/${PREFIX}_round${ROUND}.sorted.bam"
    OUTPUTS=${PREFIX}_pilon_round${ROUND}
    echo genome ../${ASSEMBLY}
    echo BAMFILE $BAMFILE
    echo output ${OUTPUTS}.fasta
    echo "${PILON} --genome ../${ASSEMBLY} --frags ${BAMFILE} ${OPTIONS} --output ${OUTPUTS} | tee round${ROUND}.pilon"
    ${PILON} --genome ../${ASSEMBLY} --frags ${BAMFILE} ${OPTIONS} --output ${OUTPUTS} | tee round${ROUND}.pilon
    echo -n $?; echo " exit status"
    echo DONE PILON
    date
    cd ..
}

for ROUND in {1}
do
	echo "ROUND: ${ROUND}"
	run_bowtie
	#echo "ROUND: ${ROUND}"
	#run_pilon
	#ASSEMBLY=pilon_round${ROUND}/${OUTPUTS}.fasta
	#echo "next ASSEMBLY: $ASSEMBLY"
done

ls -lr

