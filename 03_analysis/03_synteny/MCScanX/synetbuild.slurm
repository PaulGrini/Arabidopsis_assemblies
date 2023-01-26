
#!/bin/bash
#SBATCH --job-name=Synet
#SBATCH --output=slurm-%j.base
#SBATCH --account=NN9525K
#SBATCH --qos=devel
#SBATCH --time=00:30:00
##SBATCH --time=40:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-user=jonathan.bramsiepe@ibv.uio.no
#SBATCH --mail-type=END
##SBATCH --partition=bigmem

module --quiet purge
module load DIAMOND/0.9.30-iccifort-2019.5.281

PATH=$PATH:/cluster/projects/nn9525k/Programs/MCScanX/

./SynetBuild-X.sh 6 5 25 10