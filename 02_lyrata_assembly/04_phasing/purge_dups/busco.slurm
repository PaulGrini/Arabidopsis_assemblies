#!/bin/bash
#SBATCH --job-name=busco
#SBATCH --account=nn9525k
#SBATCH --output=slurm-%j.base
#SBATCH --time=16:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-user=jonathan.bramsiepe@ibv.uio.no
#SBATCH --mail-type=END

module purge

module load Anaconda3/2019.03

source activate busco_env

for f in petraea_pilon9/petraea_pilon9.fasta purge_work/purged.fa purge_work/hap.fa
do
 echo "f is $f"
 prefix=$(basename $f)
 echo "output is busco_embryophyta_"${prefix}""
 run_BUSCO.py -i "$f" -c 8 -o busco_embryophyta_"${prefix}" -l /cluster/projects/nn9525k/databases/odb9/embryophyta_odb9/ -m genome
done

mkdir busco_summary
cd busco_summary/
cp ../run_busco_embryophyta_*/short_summary_busco_embryophyta_*.txt .
cd ..
/cluster/projects/nn9525k/Programs/conda_envs/busco_env/bin/generate_plot -wd busco_summary

for f in $(find . -name 'run_busco*' -type d);
do back=$(pwd) ;  cd $f ;
tar --remove-files -czvf hmmer_output.tar.gz hmmer_output/ ;
tar --remove-files -czvf blast_output.tar.gz blast_output/;
tar --remove-files -czvf augustus_output.tar.gz augustus_output/;
tar --remove-files -czvf single_copy_busco_sequences.tar.gz single_copy_busco_sequences/;
cd $back ;done