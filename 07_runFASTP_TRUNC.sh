#!/bin/bash

in=$1
out=$2

for fq in $(ls ${in}/*_joined.fastq) ; do

        bn=$(basename ${fq} _joined.fastq)

        echo "[FASTP] ${bn}"

        fastp \
                -i ${in}/${bn}_joined.fastq \
                -o ${out}/${bn}_final.fastq \
                --disable_adapter_trimming \
                --length_required 248 \
		--length_limit 256 \
                --html /dev/null \
                --json /dev/null #> /dev/null 2> /dev/null

done
