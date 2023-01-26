#!/bin/bash
#SBATCH --job-name=quastlg
#SBATCH --output=slurm-%j.base
#SBATCH --account=NN9525K
#SBATCH --qos=devel
#SBATCH --time=00:30:00
##SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-user=jonathan.bramsiepe@ibv.uio.no
#SBATCH --mail-type=ALL

module --quiet purge  
module load QUAST/5.0.2-foss-2018b-Python-2.7.15

PRJDIR=$(pwd)

quast.py --large -t $SLURM_CPUS_PER_TASK ${PRJDIR}/contigs/*contigs.fasta
