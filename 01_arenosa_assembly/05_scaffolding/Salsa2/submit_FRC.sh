#!/bin/bash

PRJDIR=$(pwd)

mkdir -p FRC_curve
mkdir -p logs


ls ${PRJDIR}/bwa/*.bam > bam.fofn
njobs=$(wc -l bam.fofn |awk '{print $1}')


account=NN9525K
name=FRCurve
array="--array=1-$njobs"
wait_for=""
cpus=8
mem=2G
partition=normal
#test="--qos=devel"
#walltime=00:30:00
test=""
walltime=02:00:00
log=logs/$name.%A_%a.log
script=./FRCurve.sh
#args="$BAM 176685550"
args="176685550"
mail="--mail-user=jonathan.bramsiepe@ibv.uio.no --mail-type=ALL"

# FRCurve.sh
#FRC --pe-sam "$ALIGNMENT" --genome-size "$GENOME_SIZE" --output "${PREFIX}"
echo -e "\
sbatch --account=$account --job-name=$name $array $wait_for --cpus-per-task=$cpus  --mem-per-cpu=$mem --partition=$partition $test --time=$walltime --error=$log --output=$log $mail $script $args"
sbatch --account=$account --job-name=$name $array $wait_for --cpus-per-task=$cpus  --mem-per-cpu=$mem --partition=$partition $test --time=$walltime --error=$log --output=$log $mail $script $args > ${name}_jid
wait_for=`cat ${name}_jid |awk '{print "--dependency=afterok:"$NF}' |tr '\n' ',' |awk '{print substr($0, 1, length($0)-1)}'`
echo "wait_for $wait_for"


name=FRC_plot
#wait_for=""
cpus=2
mem=1G
partition=normal
test="--qos=devel"
walltime=00:30:00
#test=""
#walltime=12:00:00
log=logs/$name.%j.log
script=./FRC_plot.sh
args=""
mail="--mail-user=jonathan.bramsiepe@ibv.uio.no --mail-type=ALL"

# FRC_plot.sh
#FRC --pe-sam "$ALIGNMENT" --genome-size "$GENOME_SIZE" --output "${PREFIX}"
echo -e "\
sbatch --account=$account --job-name=$name $wait_for --cpus-per-task=$cpus  --mem-per-cpu=$mem --partition=$partition $test --time=$walltime --error=$log --output=$log $mail $script $args"
sbatch --account=$account --job-name=$name $wait_for --cpus-per-task=$cpus  --mem-per-cpu=$mem --partition=$partition $test --time=$walltime --error=$log --output=$log $mail $script $args > ${name}_jid
wait_for=`cat ${name}_jid |awk '{print "--dependency=afterok:"$NF}' |tr '\n' ',' |awk '{print substr($0, 1, length($0)-1)}'`
echo "wait_for $wait_for"

