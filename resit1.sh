#!/bin/bash

# Prompt the user for the input directory
echo "Please enter the path to the directory containing the input files:F.e./localdisk/data/BPSM/ICA1/fastq/"
read INPUT_DIR

# Copy the input files to the current directory
cp ${INPUT_DIR}/Tco-{5053..6950}_{1,2}.fq.gz .

# loop through the input files and run fastqc
for file in Tco-*.fq.gz; do
    fastqc "${file}" 
done
echo "FastQC reports are saved in the current directory."


# Prompt the user for the path to the reference genome FASTA file
echo "Please enter the path to the reference genome FASTA file. F.e: /localdisk/home/data/BPSM/ICA1/Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz"
read REF_GENOME

# Prompt the user for the desired name for the Bowtie2 index
echo "Please enter a name for the Bowtie2 index:F.e.TcongolenseIL3000_index"
read INDEX_NAME

# Build the Bowtie2 index using the reference genome FASTA file
bowtie2-build $REF_GENOME $INDEX_NAME

# Loop through all paired-end read files in the directory and align them to the reference genome
for file in Tco-*_1.fq.gz
do
  base=$(basename $file _1.fq.gz)
  bowtie2 -x $INDEX_NAME -1 ${base}_1.fq.gz -2 ${base}_2.fq.gz -S ${base}.sam
done



