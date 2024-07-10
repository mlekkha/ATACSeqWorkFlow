#!/bin/bash

### vars
samplesheet=$1

### Download data from SRA
while read -r condition seqnum srr
do

        prefetch ${srr}
        fasterq-dump ${srr} --skip-technical --split-files --outfile ${condition}_${seqnum}_temp.fastq

done < ${samplesheet}

### Merge the files based on condition
for i in $(ls ./*.fastq | awk -F "/" {'print $NF'} | cut -f1 -d "_" | sort | uniq)
do

    echo merging ${i} sample files
    cat ${i}*_1.fastq > ${i}_R1.fastq
    cat ${i}*_2.fastq > ${i}_R2.fastq

done

### clean up the unneccessary fastq files
rm *temp*.fastq
rm -R SRR*

### zip the data
gzip *.fastq
