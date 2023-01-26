#!/bin/bash
#SBATCH --job-name=getChloro
#SBATCH --output=slurm-%j.base
#SBATCH --account=NN9525K
##SBATCH --qos=devel
##SBATCH --time=00:30:00
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2G
#SBATCH --mail-user=jonathan.bramsiepe@ibv.uio.no
#SBATCH --mail-type=END
##SBATCH --partition=bigmem


module load GetOrganelle/1.7.1-foss-2019b-Python-3.7.4
module load pigz/2.4-GCCcore-8.3.0 

if [ -z $1 ]; then
    echo "Usage: getChloro.slurm <_R1_rep.fastq.gz> <_R2_rep.fastq.gz>"
    echo -e "\t seed:  seeds/ChrC_Athaliana_447_TAIR10.fa"
    echo -e "\t ouput: $PREFIX_plastome_output" 
    exit -1
fi

R1=$1
R2=$2
PREFIX=$(basename $R1 _R1_rep.fastq.gz)
OUT=${PREFIX}_plastome_output

if [[ ! -f $R1 || ! -f $R2 ]] ; then
    echo "reads: $R1 and $R2 are missing. Exit."
    exit -1
fi

#R1=reads/illumina_2Az1_9/Sample_06-F-2Az1-9_R1_rep.fastq.gz
#R2=reads/illumina_2Az1_9/Sample_06-F-2Az1-9_R2_rep.fastq.gz
SEED=seeds/ChrC_Athaliana_447_TAIR10.fa
#OUT=Sample_06-F-2Az1-9_plastome_output

get_organelle_from_reads.py -1 $R1 -2 $R2 -o $OUT -R 15 -k 21,45,65,85,105 -F embplant_pt -s $SEED -t 8