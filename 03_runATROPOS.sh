#!/bin/bash

in=$1
out=$2

for fq in $(ls ${in}/*_R1.fastq) ; do

        bn=$(basename ${fq} _R1.fastq)

        echo "[ATROPOS] ${bn}"

        atropos \
                -g GTGYCAGCMGCCGCGGTAA \
                -G GACTACHVGGGTATCTAATCC \
                --minimum-length=120 \
                --pair-filter=any \
                --discard-untrimmed \
                --match-read-wildcards \
		--threads 12 \
                -o ${out}/${bn}_noadapt_R1.fastq \
                -p ${out}/${bn}_noadapt_R2.fastq \
                -pe1 ${in}/${bn}_R1.fastq \
                -pe2 ${in}/${bn}_R2.fastq > /dev/null 2> /dev/null

done
