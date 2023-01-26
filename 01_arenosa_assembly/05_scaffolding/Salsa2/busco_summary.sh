#!/bin/bash

module --quiet purge
module load Anaconda3/2019.03
export PS1=\$
source ${EBROOTANACONDA3}/etc/profile.d/conda.sh
conda activate /cluster/projects/nn9525k/Programs/conda_envs/busco_env

mkdir -p BUSCO/busco_summary

if [[ ! -f fasta.fofn ]] ; then
    echo "fasta.fofn missing. Exit."
    exit -1
fi

echo "copy files"
cp BUSCO/run_*/short_summary_* BUSCO/busco_summary/
#for FASTA in fasta.fofn; do
#	PREFIX=$(basename "$FASTA" .fasta)
	#cp BUSCO/run_busco_embryophyta_${PREFIX}/short_summary_${PREFIX}.txt BUSCO/busco_summary/short_summary_${PREFIX}.txt
#done

echo "plot"
/cluster/projects/nn9525k/Programs/conda_envs/busco_env/bin/generate_plot -wd BUSCO/busco_summary

echo "tar"

for d in $(find . -name 'run_*' -type d);
do back=$(pwd) ;  cd $d ;
tar --remove-files -czvf hmmer_output.tar.gz hmmer_output/ ;
tar --remove-files -czvf blast_output.tar.gz blast_output/;
tar --remove-files -czvf augustus_output.tar.gz augustus_output/;
tar --remove-files -czvf single_copy_busco_sequences.tar.gz single_copy_busco_sequences/;
cd $back ;done

echo "finished"
