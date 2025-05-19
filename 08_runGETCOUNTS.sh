#!/bin/bash

indir=$1

if [ ! ${indir} ]; then
	echo "Missing input base directory with ./data and ./analyis subdirectories"
	exit
fi

outdir=$2

if [ ! ${outdir} ]; then
	echo "Missing output base directory for ./misc"
	exit
fi
if [ ! -d ${outdir} ]; then
	echo "Wrong output base directory for ./misc (${outdir})"
	exit
fi
	
mkdir -p ${outdir}/misc

tbl=${outdir}/misc/counts.tsv

if [ -f ${tbl} ]; then
	rm -f ${tbl}
fi

echo -e "ID\tRaw\twAdapt\tQC\tMerged\tTrunc" > ${tbl}

for fq in $(ls ${indir}/data/*_R1.fastq) ; do

        bn=$(basename ${fq} _R1.fastq)

        raw=$(cat ${indir}/data/${bn}_R1.fastq | wc -l | sed 's/ /\t/g' | awk '{print $1/4}')
        apt=$(cat ${indir}/analysis/atropos/${bn}_noadapt_R1.fastq | wc -l | sed 's/ /\t/g' | awk '{print $1/4}')
        qlc=$(cat ${indir}/analysis/fastp/${bn}_clean_R1.fastq | wc -l | sed 's/ /\t/g' | awk '{print $1/4}')
        joi=$(cat ${indir}/analysis/flash/${bn}_joined.fastq | wc -l | sed 's/ /\t/g' | awk '{print $1/4}')
        trc=$(cat ${indir}/analysis/final/${bn}_final.fastq | wc -l | sed 's/ /\t/g' | awk '{print $1/4}')

        echo -e "${bn}\t${raw}\t${apt}\t${qlc}\t${joi}\t${trc}" >> ${tbl}

done
