install and run PBJelly
================
2020-06-20

  - [Install](#install)
      - [Requirements](#requirements)
  - [Run PBJelly](#run-pbjelly)
  - [add fake Phred score to fasta](#add-fake-phred-score-to-fasta)
  - [Protocol.xml](#protocol.xml)
  - [Setup your files](#setup-your-files)
  - [Mapping your data](#mapping-your-data)
  - [Support The Gaps](#support-the-gaps)
  - [Extract Useful Reads](#extract-useful-reads)
  - [Support The Gaps & Extract Useful Reads as slurm
    script](#support-the-gaps-extract-useful-reads-as-slurm-script)
  - [Assemble The Gaps](#assemble-the-gaps)
  - [Output Your Results](#output-your-results)
  - [tar copy and save](#tar-copy-and-save)
  - [continue with](#continue-with)
  - [Extras](#extras)

## Install

some help:

<https://github.com/esrice/PBJelly>

<https://github.com/cgjosephlee/PBJelly>

<https://github.com/alvaralmstedt/Tutorials/wiki/Gap-closing-with-PBJelly>

<https://github.com/PacificBiosciences/blasr/wiki>

**Reminder:**

Sometimes it fails without returning error code, go check each xxx.err.

### Requirements

``` bash
#Blasr 
#Python 2.7 
#Networkx 2.2 
#biopython for fasta_to_fastq.py (from Bio import SeqIO)

conda create -n PBJelly python=2.7 blasr networkx biopython
conda activate /cluster/projects/nn9525k/Programs/conda_envs/PBJelly

cd /cluster/projects/nn9525k/Programs/conda_envs/PBJelly
git clone https://github.com/esrice/PBJelly.git

#This is the path where you've install the suite.
export SWEETPATH=/cluster/projects/nn9525k/Programs/conda_envs/PBJelly/PBJelly \
#for python modules
export PYTHONPATH=$PYTHONPATH:$SWEETPATH \
#for executables
export PATH=$PATH:$SWEETPATH/bin/ \

cd $CONDA_PREFIX
mkdir -p ./etc/conda/activate.d
#mkdir -p ./etc/conda/deactivate.d
touch ./etc/conda/activate.d/env_vars.sh
#touch ./etc/conda/deactivate.d/env_vars.sh

echo "
#!/bin/sh
#This is the path where you've install the suite. 
export SWEETPATH=/cluster/projects/nn9525k/Programs/conda_envs/PBJelly/PBJelly 
#for python modules 
export PYTHONPATH=$PYTHONPATH:$SWEETPATH 
#for executables 
export PATH=$PATH:$SWEETPATH/bin/ 
"  >> /cluster/projects/nn9525k/Programs/conda_envs/PBJelly/etc/conda/activate.d/env_vars.sh

# test
Jelly.py -h
```

## Run PBJelly

``` bash
cd /cluster/work/users/jonathbr
mkdir PBJelly
cd PBJelly
```

``` bash
# reads 
cp /cluster/projects/nn9525k/hybrids/jasonrm/GenomeAsssembly/A.lyrata.ver1/lyrata.correctedReads.fasta .

# Assembly
cp /cluster/projects/nn9525k/jonathan/ragoo_petraea_pilon9/petraea_canu_pilon9_purged_ragoo.fasta .
wc #110       110 184504108
```

## add fake Phred score to fasta

“Fasta files (reference or reads) must have an associated qual file with
the same name sitting in the same directory beside the .fasta file.”

“Qual files should contain the Phred Scores of bases (0-93) and should
not be encoded (i.e. no Sanger/Solexa, only the number for the score) If
you do not have qualities for your sequences, use fakeQuals.py –help”

``` bash
#fakeQuals.py lyrata.correctedReads.fasta petraea.pbcorr.fastq
#./fasta_to_fastq.py lyrata.correctedReads.fasta petraea.pbcorr.fastq 

#module load BBMap/38.50b-GCC-8.2.0-2.31.1
#reformat.sh in=test.fa out=fake.fq qfake=35
#reformat.sh in=lyrata.correctedReads.fasta out=petraea.pbcorr.fastq qfake=40

####fakeQuals.py data/reference/lambda.fasta data/reference/lambda.qual
fakeQuals.py petraea_canu_pilon9_purged_ragoo.fasta petraea_canu_pilon9_purged_ragoo.qual
fakeQuals.py reads/lyrata.correctedReads.fasta reads/lyrata.correctedReads.qual


# fix read names # Sequence names CANNOT have spaces.
cat lyrata.correctedReads.fasta.original | sed 's/ id=.*$//' > lyrata.correctedReads.fasta
fakeQuals.py lyrata.correctedReads.fasta lyrata.correctedReads.qual
```

``` bash
summarizeAssembly.py petraea_canu_pilon9_purged_ragoo.fasta > sum_petra_canu.txt
```

## Protocol.xml

``` bash
<jellyProtocol>
    <reference>/cluster/work/users/jonathbr/PBJelly/petraea_canu_pilon9_purged_ragoo.fasta</reference>
    <outputDir>/cluster/work/users/jonathbr/PBJelly/jelly_out/</outputDir>
    <blasr>--minMatch 8 --maxMatch 20 --sdpTupleSize 8 --minPctIdentity 75 --bestn 2 --nCandidates 10 --maxScore -500 --nproc 20 --noSplitSubreads</blasr>
    <input baseDir="/cluster/work/users/jonathbr/PBJelly/reads/">
        <job>lyrata.correctedReads.fasta</job>
    </input>
</jellyProtocol>
```

## Setup your files

``` bash
Jelly.py setup Protocol.xml
```

## Mapping your data

``` bash
#Jelly.py mapping Protocol.xml
#change Protocol.xml to
#--nproc 20 #20 threads
sbatch jelly_mapping.slurm #job 873975 09:56:03   00:30:15  9 420 884K TIME LIMIT (20cpus)
sbatch jelly_mapping.slurm #job 874058 879442  19-04:59:+   11:41:35 41 748 156K (40cpus)
```

## Support The Gaps

``` bash
#login-1
screen -S jelly
srun --cpus-per-task=4 --mem-per-cpu=4G --time=05:00:00 --account=nn9525k --pty bash -i
conda activate /cluster/projects/nn9525k/Programs/conda_envs/PBJelly
cd /cluster/work/users/jonathbr/PBJelly3

Jelly.py support Protocol.xml #very fast/short
#[WARNING] Read ['m54047_191216_052037/20644150/1048_36444'] self-extends Node ref0000003e5 Possible Evidence of Tandem Repeat on Singleton
```

## Extract Useful Reads

``` bash
Jelly.py extraction Protocol.xml #2:20min
```

## Support The Gaps & Extract Useful Reads as slurm script

``` bash
# Jelly.py support Protocol.xml #and 
# Jelly.py extraction Protocol.xml #included 
sbatch jelly_support.slurm #job 894161  00:06:09   00:08:17 63920952K
tail jelly_out/assembly/extraction.err
2020-06-18 14:59:13,862 [INFO] Parsed 15926 Reads
2020-06-18 15:00:03,452 [INFO] Finished
```

## Assemble The Gaps

``` bash
#Jelly.py assembly Protocol.xml -x "--nproc=16"
sbatch jelly_assembly.slurm #job 904286 4-08:00:46   13:24:59 105 315 400K
```

## Output Your Results

``` bash
#Jelly.py output Protocol.xml
sbatch jelly_output.slurm #job 905313 00:02:14   00:04:18 2879.50M
```

``` bash

grep -c "^>" jelly.out.fasta
#55
grep -Ho N jelly.out.fasta | uniq -c
# 38610 jelly.out.fasta:N

summarizeAssembly.py jelly_out/jelly.out.fasta > sum_petra_canu_jelly.txt
```

## tar copy and save

``` bash
du -sh .
#89G

rm -r reads 
rm petraea.pbcorr.fastq

tar --remove-files -czvf slurm.base.tar.gz slurm-*
tar --remove-files -czvf petraea_canu_input.tar.gz petraea_canu_pilon9_purged_ragoo.*

cp jelly_out/jelly.out.fasta petraea_canu_jelly.fasta

tar --remove-files -czvf jelly_out.tar.gz jelly_out


du -sh . #1.2G

cp -rv /cluster/work/users/jonathbr/PBJelly2 /cluster/projects/nn9525k/jonathan/PBJelly_petraea
```

## continue with

``` bash
/cluster/projects/nn9525k/jonathan/PBJelly_petraea/petraea_canu_jelly.fasta
wc  #110       110 187826863 /cluster/projects/nn9525k/jonathan/PBJelly_petraea/petraea_canu_jelly.fasta
```

## Extras

cp (<https://github.com/esrice/PBJelly>) .

``` bash
blasrToBed.py
This script will convert blasr's .m4 or .m5 format into a
BED Format file ( http://genome.ucsc.edu/FAQ/FAQformat.html#format1 )
If you would like to visualize the alignments, I
reccommend using IGB ( http://bioviz.org/igb/index.html ).

bedToCoverageWig.py
Turn a bed file with alignments into a depth of coverage
WIG Format file ( http://genome.ucsc.edu/FAQ/FAQformat.html#format6 ).
```
