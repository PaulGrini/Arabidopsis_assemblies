#!/bin/bash

module purge
module load Anaconda3/2019.03
source activate busco_env

i=$SLURM_ARRAY_TASK_ID

if [[ ! -f fasta.fofn ]] ; then
    echo "fasta.fofn missing. Exit."
    exit -1
fi

FASTA=$(sed -n ${i}p fasta.fofn)

cd BUSCO

if [[ ! -f $FASTA ]] ; then
    echo "Contigs: $FASTA missing. Exit."
    exit -1
fi

PREFIX=$(basename "$FASTA" .fasta)
echo "output is BUSCO/run_busco_embryophyta_${PREFIX}/"

run_BUSCO.py -i $FASTA -c $SLURM_CPUS_PER_TASK -o ${PREFIX} -l /cluster/projects/nn9525k/databases/odb9/embryophyta_odb9/ -m genome
