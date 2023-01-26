#!/bin/bash
#SBATCH --job-name=quastlg
#SBATCH --output=slurm-%j.base
#SBATCH --account=NN9525K
#SBATCH --qos=devel
#SBATCH --time=00:30:00
##SBATCH --time=10:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-user=jonathan.bramsiepe@ibv.uio.no
#SBATCH --mail-type=ALL

set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module --quiet purge  
module load QUAST/5.0.2-foss-2018b-Python-2.7.15
module load Anaconda3/2019.03

export PS1=\$
source ${EBROOTANACONDA3}/etc/profile.d/conda.sh
conda activate /cluster/projects/nn9525k/Programs/conda_envs/circos

PRJDIR=$(pwd)
#REFERENCE="$PRJDIR"/reference/Alyrata_384_v1.fasta
SIZE=179629584 

#quast.py --large -t $SLURM_CPUS_PER_TASK -r ${REFERENCE} --k-mer-stats --circos ${PRJDIR}/assemblies/*.fasta
#quast.py --large -t $SLURM_CPUS_PER_TASK --est-ref-size $SIZE *.fasta
quast.py --large -t $SLURM_CPUS_PER_TASK --est-ref-size $SIZE */*.fasta
