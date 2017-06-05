#!/bin/bash

# short script to take pairwise output from synbuilder and converts into an R usable format.
# synbuilder url: http://bioinfo.konkuk.ac.kr/synteny_portal/htdocs/synteny_builder.php

# $REF is our reference genome
# $QUE is our query genome
# $FILE is our input 

#invoke with:
# QUE=<queGenome> REF=<refGenome> FILE=<synBuilderOutput> bash synBuildConvert.bash > outfile

grep $REF $FILE > ref.temp
grep $QUE $FILE > que.temp

paste ref.temp que.temp > ref_que.temp

rm ref.temp que.temp

sed -i .rm '1d' ref_que.temp

rm *.rm

perl -p -e 's/ -/ c/g' ref_que.temp | perl -p -e 's/:/\t/g' | perl -p -e 's/ /\t/g' | perl -p -e 's/\./\t/g' | perl -p -e 's/-/\t/g' | perl -p -e 's/c\n/-\n/g' 

rm ref_que.temp