#!/bin/bash

# Prompt the user for the input directory
echo "Please enter the path to the directory containing the input files:F.e.   /localdisk/data/BPSM/ICA1/fastq/"
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
echo "Please enter a name for the Bowtie2 index:F.e.  TcongolenseIL3000_index"
read INDEX_NAME

# Build the Bowtie2 index using the reference genome FASTA file
bowtie2-build $REF_GENOME $INDEX_NAME
echo "Reference index is created."
# Loop through all paired-end read files in the directory and align them to the reference genome
for file in Tco-*_1.fq.gz
do
  base=$(basename $file _1.fq.gz)
  bowtie2 -x $INDEX_NAME -1 ${base}_1.fq.gz -2 ${base}_2.fq.gz -S ${base}.sam
done
#!/bin/bash
echo "Fastq files are aligned to $INDEX_NAME and *.sam files generated for each sample"
for file in *.sam
do
  samtools sort $file -o $(basename $file .sam)_sorted.bam
done
echo "Sam files sorted and saved into current directory as  *_sorted.bam"


echo "Please enter the full path+filename to  bedfile: F.e. /localdisk/data/BPSM/ICA1/TriTrypDB-46_TcongolenseIL3000_2019.bed "
read BED_DIR

# Loop over all sorted BAM files in the current directory
for file in *_sorted.bam
do
  # Get the sample name
  sample=$(basename $file .sorted.bam)

  # Use bedtools multicov to generate counts data
  bedtools multicov -bams $file -bed $BED_DIR -names $sample > $sample.counts.txt
done





















