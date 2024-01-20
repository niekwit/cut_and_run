#!/usr/bin/env bash
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR166/018/SRR16694118/SRR16694118_1.fastq.gz -o IgG_Rabbit_07_R1_001.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR166/018/SRR16694118/SRR16694118_2.fastq.gz -o IgG_Rabbit_07_R2_001.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR166/019/SRR16694119/SRR16694119_1.fastq.gz -o control_06_R1_001.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR166/019/SRR16694119/SRR16694119_2.fastq.gz -o temp/control_06_R2_001.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR166/020/SRR16694120/SRR16694120_1.fastq.gz -o temp/control_04_R1_001.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR166/020/SRR16694120/SRR16694120_2.fastq.gz -o temp/control_04_R2_001.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR166/014/SRR16694114/SRR16694114_1.fastq.gz -o temp/setdb1_R1_001.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR166/014/SRR16694114/SRR16694114_2.fastq.gz -o temp/setdb1_R2_001.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR166/013/SRR16694113/SRR16694113_1.fastq.gz -o temp/setdb1_R1_001.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR166/013/SRR16694113/SRR16694113_2.fastq.gz -o temp/setdb1_R2_001.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR166/017/SRR16694117/SRR16694117_1.fastq.gz -o temp/IgG_Rabbit_R1_001.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR166/017/SRR16694117/SRR16694117_2.fastq.gz -o temp/IgG_Rabbit_R2_001.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR166/021/SRR16694121/SRR16694121_1.fastq.gz -o temp/control_03_R1_001.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR166/021/SRR16694121/SRR16694121_2.fastq.gz -o temp/control_03_R2_001.fastq.gz

# only keep 2.5M reads
for i in $(ls temp/*.fastq.gz)
do
	zcat $i | head -n 5000000 | gzip > $(basename $i)
done




