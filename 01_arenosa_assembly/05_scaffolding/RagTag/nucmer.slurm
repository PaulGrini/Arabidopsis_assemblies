#!/bin/bash
#SBATCH --job-name=nucmer
#SBATCH --output=slurm-%j.base
#SBATCH --account=NN9525K
#SBATCH --qos=devel
#SBATCH --time=00:30:00
##SBATCH --time=10:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
##SBATCH --partition=bigmem
#SBATCH --mail-user=jonathan.bramsiepe@ibv.uio.no
#SBATCH --mail-type=END

if [ -z $1 ]; then
	echo "Usage: sbatch nucmer.slrum <reference file> <query file> [fasta]"
	echo ""
	exit -1
fi

REF=$1
ASS=$2
if [[ ! -e $REF || ! -e $ASS ]]; then
    echo "reference=$REF or query=$ASS are missing"
    exit -1
fi

module load MUMmer/4.0.0beta2-foss-2018b

PREFIX=$(basename "$ASS" .fasta)
echo "\n
nucmer --maxmatch -l 150 -c 500 -t "$SLURM_CPUS_PER_TASK" -p "$PREFIX" "$REF" "$ASS""
nucmer --maxmatch -l 150 -c 500 -t "$SLURM_CPUS_PER_TASK" -p "$PREFIX" "$REF" "$ASS"
#nucmer --maxmatch -l 100 -c 500 -t "$SLURM_CPUS_PER_TASK" -p "$PREFIX" "$REF" "$ASS"
echo "\n
dnadiff -p "$PREFIX" -d "$PREFIX".delta"
dnadiff -p "$PREFIX" -d "$PREFIX".delta

echo "goodbye"
