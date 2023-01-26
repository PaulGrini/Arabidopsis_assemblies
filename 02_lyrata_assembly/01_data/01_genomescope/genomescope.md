# genomescope-GS estimation

2021-05-19

-   [Setup directory](#setup-directory)
-   [input.fofn](#input.fofn)
-   [Get the right k size](#get-the-right-k-size)
-   [Build k-mer dbs with meryl](#build-k-mer-dbs-with-meryl)
-   [srun genomescope](#srun-genomescope)

## Links

https://github.com/marbl/merqury

https://github.com/tbenavi1/genomescope2.0


## Setup directory {#setup-directory}

``` bash
cd /cluster/work/users/jonathbr/merqury
ln -s $MERQURY/merqury.sh . 

mkdir reads
cp /cluster/projects/nn9525k/hybrids/jasonrm/GenomeAsssembly/A.lyrata.ver1/lyrata.correctedReads.fasta reads/

cp /cluster/projects/nn9525k/hybrids/molbar_illumina_DNA/trimmed_reads/Sample_02-B-2Lz3-4/Sample_02-B-2Lz3-4_R?_rep.fastq.gz reads/
```

## input.fofn {#input.fofn}

``` bash
ls reads/*.fastq.gz > input.fofn
```

## Get the right k size {#get-the-right-k-size}

``` bash
sh $MERQURY/best_k.sh 180000000
#k=19
```

/cluster/projects/nn9525k/Programs/merqury/\_submit_build.sh has been changed for use on saga

## Build k-mer dbs with meryl {#build-k-mer-dbs-with-meryl}

``` bash
$MERQURY/_submit_build.sh 19 input.fofn petrea_illumina
```

``` bash
cat petrea_illumina.k19.hist | sed 's/\t/ /' > petraea_illumina.k19_space.hist
#rm arenosa_illumina.k19.hist
```

## srun genomescope {#srun-genomescope}

``` bash
srun --cpus-per-task=2 --mem-per-cpu=2G --time=00:30:00 --account=nn9525k --x11 --pty bash -i
cd /cluster/work/users/jonathbr/petraea_genomescope/merqury
module purge
module load R/3.6.2-foss-2019b
Rscript /cluster/projects/nn9525k/Programs/genomescope2.0/genomescope.R -i petraea_illumina.k19_space.hist -o genomescope -n petraea_illumina -k 19 -p 2
```

| petraea_illumina_log               | petraea_illumina_linear               |
|------------------------------------|---------------------------------------|
| ![](petraea_illumina_log_plot.png) | ![](petraea_illumina_linear_plot.png) |

length: 201,144,702 bp
