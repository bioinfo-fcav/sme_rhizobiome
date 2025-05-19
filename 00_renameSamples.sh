#!/bin/bash

in_dir=$1
out_dir=$2
meta_table=$3

for i in $(ls ${in_dir}/*_R1_*.fastq) ; do

	bn=$(basename ${i})

        n_old=$(grep -w ${bn} ${meta_table} | cut -f 3)
        n_new=$(grep -w ${bn} ${meta_table} | cut -f 2)

	if [ ! ${n_old} == "" ]; then
		
		ln -s ${in_dir}/${n_old} ${out_dir}/${n_new}_R1.fastq
		ln -s ${in_dir}/$(echo ${n_old} | sed 's/_R1_/_R2_/g') ${out_dir}/${n_new}_R2.fastq
	fi

done
