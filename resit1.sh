#!/bin/bash

# specify the directory where the input files are located
input_dir=/localdisk/data/BPSM/ICA1/fastq

# specify the output directory where the fastqc reports will be saved
output_dir=./

# loop through the input files and run fastqc
for i in $(seq 5053 6950); do
    fastqc ${input_dir}/Tco-${i}_1.fq.gz -o ${output_dir}
    fastqc ${input_dir}/Tco-${i}_2.fq.gz -o ${output_dir}
done
