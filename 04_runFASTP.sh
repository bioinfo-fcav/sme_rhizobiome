#!/bin/bash

in=$1
out=$2

for fq in $(ls ${in}/*_noadapt_R1.fastq) ; do

        bn=$(basename ${fq} _noadapt_R1.fastq)

        echo "[FASTP] ${bn}"

        fastp \
                -i ${in}/${bn}_noadapt_R1.fastq \
		-I ${in}/${bn}_noadapt_R2.fastq \
                -o ${out}/${bn}_clean_R1.fastq \
		-O ${out}/${bn}_clean_R2.fastq \
                --length_required 120 \
		--average_qual 25 \
		--thread 12 \
                --html /dev/null \
                --json /dev/null > /dev/null 2> /dev/null

done
