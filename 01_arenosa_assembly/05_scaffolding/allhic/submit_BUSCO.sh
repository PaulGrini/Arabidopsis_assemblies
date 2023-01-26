#!/bin/bash


PRJDIR=$(pwd)

mkdir -p BUSCO
mkdir -p logs


ls ${PRJDIR}/*/*.fasta > fasta.fofn
njobs=$(wc -l fasta.fofn |awk '{print $1}')


account=NN9525K
name=busco
array="--array=1-$njobs"
wait_for=""
cpus=8
mem=4G
partition=normal
#test="--qos=devel"
#walltime=00:30:00
test=""
walltime=12:00:00
log=logs/$name.%A_%a.log
script=./busco.sh
args=""
mail="--mail-user=jonathan.bramsiepe@ibv.uio.no --mail-type=ALL"

# busco.sh
echo -e "\
sbatch --account=$account --job-name=$name $array $wait_for --cpus-per-task=$cpus  --mem-per-cpu=$mem --partition=$partition $test --time=$walltime --error=$log --output=$log $mail $script $args"
sbatch --account=$account --job-name=$name $array $wait_for --cpus-per-task=$cpus  --mem-per-cpu=$mem --partition=$partition $test --time=$walltime --error=$log --output=$log $mail $script $args > ${name}_jid
wait_for=`cat ${name}_jid |awk '{print "--dependency=afterok:"$NF}' |tr '\n' ',' |awk '{print substr($0, 1, length($0)-1)}'`
echo "wait_for: $wait_for"


name=busco_summary
#wait_for=""
cpus=2
mem=1G
partition=normal
test="--qos=devel"
walltime=00:30:00
#test=""
#walltime=12:00:00
log=logs/$name.%j.log
script=./busco_summary.sh
args=""
mail="--mail-user=jonathan.bramsiepe@ibv.uio.no --mail-type=ALL"

# busco_summary.sh
echo -e "\
sbatch --account=$account --job-name=$name $wait_for --cpus-per-task=$cpus  --mem-per-cpu=$mem --partition=$partition $test --time=$walltime --error=$log --output=$log $mail $script $args"
sbatch --account=$account --job-name=$name $wait_for --cpus-per-task=$cpus  --mem-per-cpu=$mem --partition=$partition $test --time=$walltime --error=$log --output=$log $mail $script $args > ${name}_jid
wait_for=`cat ${name}_jid |awk '{print "--dependency=afterok:"$NF}' |tr '\n' ',' |awk '{print substr($0, 1, length($0)-1)}'`
echo "done $wait_for"

