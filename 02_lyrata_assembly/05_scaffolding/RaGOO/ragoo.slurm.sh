#!/bin/bash
#SBATCH --job-name=RaGOO
#SBATCH --output=slurm-%j.base
#SBATCH --account=NN9525K
#SBATCH --qos=devel
#SBATCH --time=00:30:00
##SBATCH --time=10:00:00
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-user=jonathan.bramsiepe@ibv.uio.no
#SBATCH --mail-type=END

set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error

module --quiet purge  # Reset the modules to the system default
module load Anaconda3/2019.03

#(needed in the source of the Anaconda environment)
export PS1=\$
source ${EBROOTANACONDA3}/etc/profile.d/conda.sh
conda activate /cluster/projects/nn9525k/Programs/conda_envs/RaGOO

ragoo.py -R PB_corr_reads.fasta -T corr -t 6 -g 100 -C petraea_pilon9.fasta Alyrata_384_v1.fasta
echo "output "$(pwd)"/ragoo_output"