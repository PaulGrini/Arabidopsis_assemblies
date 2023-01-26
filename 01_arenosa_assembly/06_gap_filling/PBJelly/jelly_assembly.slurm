#!/bin/bash
#SBATCH --job-name=assembly
#SBATCH --output=slurm-%j.base
#SBATCH --account=NN9525K
##SBATCH --qos=devel
##SBATCH --time=00:30:00
#SBATCH --time=40:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=10G
#SBATCH --mail-user=jonathan.bramsiepe@ibv.uio.no
#SBATCH --mail-type=END
#SBATCH --partition=bigmem

module --quiet purge
module load Anaconda3/2019.03

export PS1=\$
source ${EBROOTANACONDA3}/etc/profile.d/conda.sh
conda activate /cluster/projects/nn9525k/Programs/conda_envs/PBJelly

Jelly.py assembly Protocol.xml -x "--nproc=16"
