#!/bin/bash
#SBATCH --job-name=support
#SBATCH --output=slurm-%j.base
#SBATCH --account=NN9525K
#SBATCH --qos=devel
#SBATCH --time=00:30:00
##SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-user=jonathan.bramsiepe@ibv.uio.no
#SBATCH --mail-type=ALL
##SBATCH --partition=bigmem

module --quiet purge
module load Anaconda3/2019.03

export PS1=\$
source ${EBROOTANACONDA3}/etc/profile.d/conda.sh
conda activate /cluster/projects/nn9525k/Programs/conda_envs/PBJelly

Jelly.py support Protocol.xml

Jelly.py extraction Protocol.xml