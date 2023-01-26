#!/bin/bash
#SBATCH --job-name=orthofinder
#SBATCH --account=nn9525k
#SBATCH --output=slurm-%j.base
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=6G
#SBATCH --partition=normal
#SBATCH --mail-user=jonathan.bramsiepe@ibv.uio.no
#SBATCH --mail-type=END

set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module --quiet purge  # Reset the modules to the system default
module load Anaconda3/2019.03
module list    # For easier debugging

#source of the Anaconda environment
export PS1=\$
source ${EBROOTANACONDA3}/etc/profile.d/conda.sh
conda activate /cluster/projects/nn9525k/Programs/conda_envs/orthofinder

echo "My current directory  " $PWD

#Run orthofinder (add -t 16 to run with 16 threads)
orthofinder -t $SLURM_CPUS_PER_TASK -f proteomes