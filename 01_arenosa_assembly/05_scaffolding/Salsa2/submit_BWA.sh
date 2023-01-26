#!/bin/bash

#if [ -z $1 ]; then
#    echo "Usage: ./submit_BWA.sh <prefix.fasta> <R1> <R2> "
#    echo -e "\t Reference fasta and reads"
#    exit -1
#fi

for fasta in $(ls */*.fasta); do  
echo $fasta
R1=$(ls reads/*R1*.fastq.gz)
R2=$(ls reads/*R2*.fastq.gz)
echo $R1 $R2

#fasta=$1
#R1=$2
#R2=$3

PREFIX=$(basename $fasta .fasta)

if [[ ! -f $fasta || ! -f $R1 || ! -f $R2 ]] ; then
    echo "Reference fasta or reads: $fasta $R1 $R2 are missing. Exit."
    exit -1
fi

mkdir -p logs
mkdir -p bwa

account=NN9525K
name=BWA_$PREFIX
wait_for=""
cpus=16
mem=10G
partition=bigmem
#test="--qos=devel"
#walltime=00:30:00
test=""
walltime=10:00:00
log=logs/$name.%j.log
script=./BWA_align.sh
args="$fasta $R1 $R2"
mail="--mail-user=jonathan.bramsiepe@ibv.uio.no --mail-type=ALL"


# BWA_align.sh
#align_reads assembly.fasta read1.fastq.gz read2.fastq.gz
echo -e "\
sbatch --account=$account --job-name=$name $wait_for --cpus-per-task=$cpus  --mem-per-cpu=$mem --partition=$partition $test --time=$walltime --error=$log --output=$log $mail $script $args"
sbatch --account=$account --job-name=$name $wait_for --cpus-per-task=$cpus  --mem-per-cpu=$mem --partition=$partition $test --time=$walltime --error=$log --output=$log $mail $script $args > ${name}_jid
wait_for=`cat ${name}_jid |awk '{print "afterok:"$NF}' |tr '\n' ',' |awk '{print substr($0, 1, length($0)-1)}'`
echo "wait_for $wait_for"
done
