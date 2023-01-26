#!/bin/bash
#SBATCH --job-name=canu
#SBATCH --account=nn9525k
#SBATCH --output=slurm-%j.base
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=10G
#SBATCH --partition=bigmem



module purge
module load canu/1.8-GCCcore-8.2.0-Perl-5.28.1 

#STR=/cluster/work/users/anderkkr/01_arenosa_genome/01_aenosa_pacbio/fasta/arenosa_lordec_corrected.fasta.gz
STR=/cluster/work/users/anderkkr/01_arenosa_genome/01_aenosa_pacbio/fastq/reads.fastq.gz
canu -d canu_out -pacbio-raw "$STR" -p arenosa_pb_uncorr_assembly genomeSize=200m useGrid=false 

