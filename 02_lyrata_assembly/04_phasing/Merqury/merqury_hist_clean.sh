#!/bin/bash

#cat arenosa_illumina.k19.hist | sed 's/\t/ /' > arenosa_illumina.k19_space.hist

#rm arenosa_illumina.k19.hist
rm meryl_union_sum.jid
rm meryl_count.meryl.list
rm meryl_count.jid
rm input.fofn

for f in $(find . -name '*.meryl' -type d);
do back=$(pwd) ;
 dir=$(dirname $f);
 base=$(basename $f)
 cd $dir;
  tar --remove-files -czvf ${base}.tar.gz ${base} ;
 cd $back ;
done
