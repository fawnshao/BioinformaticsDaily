#!/bin/bash
mkdir FastqcRes
fastqc -t 48 -o FastqcRes/ ./*fastq.gz
cd FastqcRes/
for f in *.zip; do unzip $f & done
wait
grep -A 2 "Overrepresented sequences" *_R1_001_fastqc/fastqc_data.txt | sed -n '3~4p'
