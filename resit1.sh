#!/bin/bash"
# Prompt the user for the input directory
##echo "Please enter the path to the directory containing the input files:F.e.   /localdisk/data/BPSM/ICA1/fastq/"
##read INPUT_DIR
##echo "Dear USER, this directory contains some other files, to which we dont have access. Therefore we will copy the ones we need into current directory, it will be a bit messy now, but will help us to reduce other messines further. Thanks for understanding and collaboration"
# Copy the input files to the current directory
##cp ${INPUT_DIR}/Tco-{5053..6950}_{1,2}.fq.gz .

# loop through the input files and run fastqc
#for file in Tco-*.fq.gz; do
#    fastqc "${file}" 
#done
#echo "FastQC reports are saved in the current directory."
#echo "Please analyse them before further processing"

# Prompt the user for the path to the reference genome FASTA file
#echo "Please enter the path to the reference genome FASTA file. F.e: /localdisk/home/data/BPSM/ICA1/Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz"
#read REF_GENOME

# Prompt the user for the desired name for the Bowtie2 index
#echo "Please enter a name for the Bowtie2 index:F.e.  TcongolenseIL3000_index"
#read INDEX_NAME

# Build the Bowtie2 index using the reference genome FASTA file
#bowtie2-build $REF_GENOME $INDEX_NAME
#echo "Reference index is created."
# Loop through all paired-end read files in the directory and align them to the reference genome
#for file in Tco-*_1.fq.gz
#do
  # Extract the common part of the filename
 # base=$(basename $file _1.fq.gz)

  # Align the paired-end reads to the reference genome using Bowtie2
 # bowtie2 -x $INDEX_NAME -1 ${base}_1.fq.gz -2 ${base}_2.fq.gz -S ${base}.sam
  
  # Convert the SAM file to BAM file and sort it
 # samtools view -bS ${base}.sam | samtools sort -o ${base}.sorted.bam
  
  # Index the sorted BAM file
  #samtools index ${base}.sorted.bam
  
  # Remove the intermediate SAM file
  #rm ${base}.sam
#done
echo "Fastq files are aligned to ${INDEX_NAME} and *.sam files generated for each sam"
echo "Sam files sorted, indexed and saved into current directory as  *_sorted.bam"


echo "Please enter the full path+filename to  bedfile: F.e. /localdisk/data/BPSM/ICA1/TriTrypDB-46_TcongolenseIL3000_2019.bed "
read BED_DIR
# Create an array of sorted BAM file paths in sampleID order
bam_files=($(ls -1v *.sorted.bam))

# Loop over the sorted BAM files in sampleID order
for file in "${bam_files[@]}"
do
  # Skip any files that end in ".sorted.bam.bai"
  if [[ $file == *.sorted.bam.bai ]]; then
    continue
  fi

  # Get the sample name
  sample=$(basename "$file" .sorted.bam)

  # Use bedtools intersect to identify overlaps with the genomic regions
  bedtools intersect -a "$BED_DIR" -b "$file" -bed > "$sample.overlaps.bed"

  # Use bedtools multicov to generate counts data
  bedtools multicov -bams "$file" -bed "$BED_DIR" > "$sample.counts.txt"
done

# Concatenate all the count data files into a common count data file with all the samples ordered
cat *.counts.txt | sort -k1,1 > all_counts.txt






































