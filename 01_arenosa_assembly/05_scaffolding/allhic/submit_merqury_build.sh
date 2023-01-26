#!/bin/bash

PRJDIR=$(pwd)

mkdir -p MERQURY

ls ${PRJDIR}/reads/*.fastq.gz > input.fofn
export MERQURY=/cluster/projects/nn9525k/Programs/merqury

#$MERQURY/_submit_build.sh 19 input.fofn arenosa_illumina
$MERQURY/_submit_build_2h.sh 19 input.fofn arenosa_illumina

#cat arenosa_illumina.k19.hist | sed 's/\t/ /' > arenosa_illumina.k19_space.hist
#rm arenosa_illumina.k19.hist
