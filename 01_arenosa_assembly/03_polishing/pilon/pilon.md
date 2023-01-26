# copy data

``` bash
cd pilon_pb_uncorr

mkdir -p data

ln -s /cluster/work/users/jonathbr/arenosa_assembly/01_Data/contigs/arenosa_pb_uncorr_assembly.contigs.fasta data/arenosa_pb_uncorr.fasta

ln -s /cluster/work/users/jonathbr/arenosa_assembly/01_Data/reads/*.fastq.gz data
```

# run pilon loops

``` bash
sbatch pilon_loop_1.slurm
sbatch pilon_loop_2.slurm
sbatch pilon_loop_3.slurm
sbatch pilon_loop_4.slurm
```

# fast commands for pilon check

``` bash
wc -l */*.changes

grep "overall alignment rate" slurm*
```

# Rename header, save and continue with

``` bash
cd ..

cat pilon_pb_uncorr/pilon_round12/arenosa_pb_uncorr_pilon_round12.fasta | sed 's/_pilon_pilon_pilon_pilon_pilon_pilon_pilon_pilon_pilon_pilon_pilon_pilon/_arenosa_pilon12/' > arenosa_pilon12.fasta

wc arenosa_pilon12.fasta #3196391   3196391 258828483 arenosa_pilon12.fasta
```

# tar and remove leftovers

``` bash
cd pilon_pb_uncorr

mv */*.changes .
mv */*.flagstat .
mv */*.SN .

rm -r data
rm -r pilon_round*
rm -r bowtie_round*

cd ..

tar -czvf pilon_pb_uncor.tar.gz pilon_pb_uncorr
rm -r pilon_pb_uncorr

gzip arenosa_pilon12.fasta
```
