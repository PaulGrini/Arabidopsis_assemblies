RaGOO on pilon9
================
2020-06-10

  - [setup dir](#setup-dir)
  - [links](#links)
  - [rename pilon fasta](#rename-pilon-fasta)
  - [RaGOO](#ragoo)
  - [tar copy and save](#tar-copy-and-save)
  - [continue with](#continue-with)

## setup dir

``` bash
cd /cluster/work/users/jonathbr/
mkdir ragoo_petraea_pilon9
cd ragoo_petraea_pilon9
```

## links

``` bash
## Reference
ln -s /cluster/projects/nn9525k/lyrata_genome_phytomzome/Alyrata_384_v1.fa Alyrata_384_v1.fasta

## Reads 
ln -s /cluster/projects/nn9525k/hybrids/jasonrm/GenomeAsssembly/A.lyrata.ver1/lyrata.correctedReads.fasta PB_corr_reads.fasta

ln -s /cluster/projects/nn9525k/hybrids/molbar_illumina_DNA/trimmed_reads/Sample_02-B-2Lz3-4/Sample_02-B-2Lz3-4_R?_rep.fastq.gz .

# Assembly
ln -s /cluster/projects/nn9525k/jonathan/pilon_loop/pilon_round9/petraea.contigs_pilon_round9.fasta . #compare purged haploid vs diploid
cp /cluster/projects/nn9525k/jonathan/purge_dups_petraea_pilon9/petraea_canu_pilon9_purged.fasta .
```

## rename pilon fasta

``` bash
cat petraea.contigs_pilon_round9.fasta | sed 's/_pilon_pilon_pilon_pilon_pilon_pilon_pilon_pilon_pilon/_petraea_pilon9/' > petraea_pilon9.fasta
```

## RaGOO

``` bash
sbatch ragoo.slurm # petraea_pilon9.fasta diploid 
#job job 856477 01:09:08   00:12:33
cd ragoo_output
mv ragoo.fasta petraea_canu_pilon9_ragoo_dontuse.fasta
cd ..
mv ragoo_output ragoo_output_diploid

sbatch ragoo_purged.slurm #job 870317 00:32:43   00:06:25 2367204K
```

## tar copy and save

``` bash
ln -s ragoo_output/ragoo.fasta petraea_canu_pilon9_purged_ragoo.fasta
wc !$ #110       110 184504108 petraea_canu_pilon9_purged_ragoo.fasta

for f in $(find . -name 'ragoo_output*' -type d);
do back=$(pwd) ;  cd $f ;
tar --remove-files -czvf groupings.tar.gz groupings/ ;
tar --remove-files -czvf orderings.tar.gz orderings/;
cd $back ;done

rm Sample_02-B-2Lz3-4_R?_rep.fastq.gz
rm PB_corr_reads.fasta
rm Alyrata_384_v1.fasta
du -sh . #2.9G

cp -rv /cluster/work/users/jonathbr/ragoo_petraea_pilon9 /cluster/projects/nn9525k/jonathan/
```

## continue with

``` bash
/cluster/projects/nn9525k/jonathan/ragoo_petraea_pilon9/petraea_canu_pilon9_purged_ragoo.fasta
wc #110       110 184504108
```
