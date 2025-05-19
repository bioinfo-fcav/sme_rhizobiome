#!/bin/bash

in=$1
out=$2

for fq in $(ls ${in}/*fastq) ; do

        bn=$(basename ${fq} .fastq)

        echo "[FastQC] ${bn}"

        fastqc \
                ${in}/${bn}.fastq \
                --outdir ${out} \
		--nogroup > /dev/null 2> /dev/null

done

cd ${out} 

multiqc ./ -n $(basename $in)

rm *.zip
