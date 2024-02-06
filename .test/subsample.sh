#!/usr/bin/env bash

mkdir -p subsample

FASTQ=$(ls *.fastq.gz)

for file in $FASTQ;
do
	zcat $file | head -1000000 | pigz > "subsample/${file}" 
done
