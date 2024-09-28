#!/bin/bash -V
source /home/yh/3-tool-ref/anaconda3/etc/profile.d/conda.sh
conda activate rna-seq
cd /media/yh/Nano2/YL/m6a/bamraw
bampath=/media/yh/Nano2/YL/m6a/bamnpt
bed=/media/yh/Nano2/YL/m6a/Arabidopsis_thaliana.include_regions.bed
for i in `ls *.bam`;
do
	bedtools intersect -v -abam $i -b /media/yh/Nano2/YL/m6a/rRNA.bed > ${bampath}/filtered.bam
	samtools index ${bampath}/filtered.bam
	samtools view -h -q 60 -L $bed ${bampath}/filtered.bam | grep -P "(NH:i:1|^@)" | samtools view -Sb - > ${bampath}/${i}
	samtools index ${bampath}/${i}
	samtools view -h ${bampath}/${i} | awk 'substr($0,1,1)=="@" || ($9>-20000 && $9<20000)' | samtools view -Sb > ${bampath}/${i/bam/u20k.bam}
	samtools index ${bampath}/${i/bam/u20k.bam}
	rm ${bampath}/${i} ${bampath}/${i}.bai ${bampath}/filtered.bam ${bampath}/filtered.bam.bai
done
