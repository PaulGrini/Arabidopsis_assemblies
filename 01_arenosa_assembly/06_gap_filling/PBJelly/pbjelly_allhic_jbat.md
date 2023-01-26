PBJelly on allhic\_jbat\_v4
================
2020-12-01

  - [add fake Phred score to fasta](#add-fake-phred-score-to-fasta)
  - [Jelly.py setup](#jelly.py-setup)
  - [Jelly.py mapping](#jelly.py-mapping)
  - [Support The Gaps & Extract Useful Reads as slurm
    script](#support-the-gaps-extract-useful-reads-as-slurm-script)
  - [Assemble The Gaps](#assemble-the-gaps)
  - [Output Results](#output-results)
  - [tar copy and save](#tar-copy-and-save)
  - [continue with](#continue-with)
  - [Extras](#extras)

``` bash
cd $USERWORK

mkdir pbjelly_jbat
cd pbjelly_jbat

## Reads 
zcat /cluster/projects/nn9525k/arenosa_genome_pb_uncorr/canu_out/arenosa_pp_uncorr_assembly.correctedReads.fasta.gz > pbReads.fasta

## Assembly
cp /cluster/projects/nn9525k/arenosa_genome/02_scaffolds/07_jbat/12_arenosa_purged_sr_groups_13_groups_jbat_juicer_v4/12_arenosa_purged_sr_groups_13_groups.asm.FINAL.fasta arenosa_jbat.fasta.bak

cp /cluster/projects/nn9525k/jonathan/arenosa_assembly/06_gap_filling/PBJelly_v3/simplifyFastaHeaders.pl .

./simplifyFastaHeaders.pl arenosa_jbat.fasta.bak arenosa arenosa_jbat.fasta header.map

fakeQuals.py arenosa_jbat.fasta arenosa_jbat.qual

mkdir jellyout
```

``` bash

# fix read names # Sequence names CANNOT have spaces.
cat pbReads.fasta | sed 's/ id=.*$//' > PBreads.fasta # can be skipped and done later

split -l 30000 PBreads.fasta pbReads.split.

for f in pbReads.split.* ; do  mv $f ${f}.fasta ; done
for f in pbReads.split.* ; do  echo "<job>${f}</job>" ; done
```

``` bash
------- Protocol.xml
<jellyProtocol>
    <reference>/cluster/work/users/jonathbr/pbjelly_jbat/arenosa_jbat.fasta</reference>
    <outputDir>/cluster/work/users/jonathbr/pbjelly_jbat/jellyout/</outputDir>
    <blasr>--minMatch 8 --sdpTupleSize 8 --minPctSimilarity 75 --bestn 1 --nCandidates 10 --maxScore -500 --noSplitSubreads</blasr>
    <input baseDir="/cluster/work/users/jonathbr/pbjelly_jbat/">
    <job>pbReads.split.aa.fasta</job>
<job>pbReads.split.ab.fasta</job>
<job>pbReads.split.ac.fasta</job>
<job>pbReads.split.ad.fasta</job>
<job>pbReads.split.ae.fasta</job>
<job>pbReads.split.af.fasta</job>
<job>pbReads.split.ag.fasta</job>
<job>pbReads.split.ah.fasta</job>
<job>pbReads.split.ai.fasta</job>
<job>pbReads.split.aj.fasta</job>
<job>pbReads.split.ak.fasta</job>
<job>pbReads.split.al.fasta</job>
<job>pbReads.split.am.fasta</job>
<job>pbReads.split.an.fasta</job>
<job>pbReads.split.ao.fasta</job>
<job>pbReads.split.ap.fasta</job>
<job>pbReads.split.aq.fasta</job>
<job>pbReads.split.ar.fasta</job>
<job>pbReads.split.as.fasta</job>
<job>pbReads.split.at.fasta</job>
<job>pbReads.split.au.fasta</job>
<job>pbReads.split.av.fasta</job>
    </input>
</jellyProtocol>
----------
```

## add fake Phred score to fasta

``` bash
conda activate /cluster/projects/nn9525k/Programs/conda_envs/PBJelly

##fakeQuals.py data/reference/lambda.fasta data/reference/lambda.qual
fakeQuals.py arenosa_jbat.fasta arenosa_jbat.qual

for FASTA in pbReads.split.*.fasta;
do 
 PREFIX=$(basename $FASTA .fasta)
 mv $FASTA ${FASTA}.original
 echo -e "\n
 cat ${FASTA}.original | sed 's/ id=.*$//' > ${FASTA}"
 cat ${FASTA}.original | sed 's/ id=.*$//' > ${FASTA}
 echo -e "\n
 fakeQuals.py $FASTA ${PREFIX}.qual"
 fakeQuals.py $FASTA ${PREFIX}.qual
done
```

``` bash
summarizeAssembly.py arenosa_jbat.fasta > sum_arenosa.txt
```

## Jelly.py setup

``` bash
Jelly.py setup Protocol.xml
```

## Jelly.py mapping

``` bash
sbatch jelly_mapping.slurm #job 1588445
#27-07:15:+ 1-11:4 13 395 344K
```

## Support The Gaps & Extract Useful Reads as slurm script

``` bash
# Jelly.py support Protocol.xml #and 
# Jelly.py extraction Protocol.xml #included 
sbatch jelly_support.slurm #8G job 1594842 00:05:11   00:14:00  3 185 248K
```

## Assemble The Gaps

``` bash
#Jelly.py assembly Protocol.xml -x "--nproc=16"
sbatch jelly_assembly.slurm  #job 1595029 6-16:16:03   15:57:32 101 215 408K
```

## Output Results

``` bash
conda activate /cluster/projects/nn9525k/Programs/conda_envs/PBJelly

Jelly.py output Protocol.xml

grep -c "^>" jellyout/jelly.out.fasta
# 264
grep -c "^>" arenosa_jbat.fasta.bak
# 266

grep -Ho N jellyout/jelly.out.fasta | uniq -c
# 15205 jellyout/jelly.out.fasta:N 
grep -Ho N arenosa_jbat.fasta.bak | uniq -c
# 16352 arenosa_jbat.fasta.bak:N

grep -c "filled" jellyout/gap_fill_status.txt
#45
grep -c "extend" jellyout/gap_fill_status.txt
#78

wc -l jellyout/gap_fill_status.txt
#298 jellyout/gap_fill_status.txt

summarizeAssembly.py jellyout/jelly.out.fasta > sum_arenosa_jelly.txt
```

## tar copy and save

``` bash
du -sh .
# 127G

#remove input 
#rm *.fasta*
#rm *.qual
mv arenosa_jbat.fasta.original arenosa_jbat.original
gzip arenosa_jbat.original

rm arenosa_jbat.fasta*
rm arenosa_jbat.qual

#remove reads
rm pbReads*
rm PBreads.fasta

rm core.*

tar --remove-files -czvf slurm.base.tar.gz slurm-*

cp jellyout/jelly.out.fasta arenosa_jbat_jelly.fasta
gzip arenosa_jbat_jelly.fasta

tar --remove-files -czvf jelly_out.tar.gz jellyout

du -sh . #156M

mkdir /cluster/projects/nn9525k/arenosa_genome/02_scaffolds/08_pbjelly

cp -rv /cluster/work/users/jonathbr/pbjelly_jbat/* /cluster/projects/nn9525k/arenosa_genome/02_scaffolds/08_pbjelly
```

## continue with

``` bash
wc /cluster/projects/nn9525k/arenosa_genome/02_scaffolds/08_pbjelly/arenosa_jbat_jelly.fasta.gz
#158631   864616 41256291 /cluster/projects/nn9525k/arenosa_genome/02_scaffolds/08_pbjelly/arenosa_jbat_jelly.fasta.gz
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
