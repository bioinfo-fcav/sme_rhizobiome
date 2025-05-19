#!/bin/bash

in=$1
out=$2

for fq in $(ls ${in}/*fastq) ; do

        bn=$(basename ${fq} .fastq)

        echo "[TRIMSIZE] ${bn}"

        usearch11 \
                -fastq_eestats2  ${in}/${bn}.fastq \
                -output ${out}/${bn}_trimSizes.txt \
                -length_cutoffs 240,260,2 \
                -ee_cutoffs 1,2,3,4 #> /dev/null 2> /dev/null

done
