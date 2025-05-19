#!/bin/bash

in=$1
out=$2

for fq in $(ls ${in}/*_clean_R1.fastq) ; do

        bn=$(basename ${fq} _clean_R1.fastq)

        echo "[FLASH] ${bn}"

        flash \
                ${in}/${bn}_clean_R1.fastq \
                ${in}/${bn}_clean_R2.fastq \
                --min-overlap 4 \
                --max-overlap 250 \
                --output-prefix ${out}/${bn} #> /dev/null 2> /dev/null

        rm -rf ${out}/*.hist*
        rm -rf ${out}/*.notCombined_*
        mv ${out}/${bn}.extendedFrags.fastq ${out}/${bn}_joined.fastq

done
