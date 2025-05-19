#!/bin/bash

in=$1
out=$2

for fq in $(ls ${in}/*fastq) ; do

        bn=$(basename ${fq} .fastq)

        echo "[FINDADAPT] ${bn}"

        head -n 2000 ${fq} > ${in}/tmp.fq

        usearch11 \
                -search_oligodb ${in}/tmp.fq \
                -db ../misc/primers.fa \
                -maxdiffs 2 \
                -strand both \
                -userout ${out}/${bn}_primers.txt \
                -userfields query+target+qstrand+diffs+tlo+thi+qlor+qhir+trowdots > /dev/null 2> /dev/null

        rm -rf ${in}/tmp.fq

done

