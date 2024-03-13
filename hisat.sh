#!/bin/bash
#SBATCH --job-name=hisat
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 15
#SBATCH --mem=220G
#SBATCH --partition=himem
#SBATCH --qos=himem
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nidhi.vijayan@uconn.edu
#SBATCH -o hisat_%j.out
#SBATCH -e hisat_%j.err


module load hisat2
module load samtools/1.9

hisat2-build /labs/Nyholm/Es_V2_Genome/Lachesis_assembly.fasta /labs/Nyholm/HISAT2_Genome_V2/Es_genome

FASTQ_FILE=/labs/Nyholm/Oleg_project/HCVFYDSX3_1_R13227_20220410/demultiplexed/ANG_CC/TRIMMED/*.1


CLEANED_READS=/labs/Nyholm/Oleg_project/HCVFYDSX3_1_R13227_20220410/demultiplexed/ANG_CC/TRIMMED/

#BAM=/labs/Nyholm/Oleg_project/HCVFYDSX3_1_R13227_20220410/demultiplexed/ANG_CC/Star_BAM/


input=/labs/Nyholm/Oleg_project/HHH72DSX2_1_R11994_20210815/demultiplexed/ANG_CC/
out=/labs/Nyholm/Oleg_project/HHH72DSX2_1_R11994_20210815/demultiplexed/ANG_CC/hisat_V2/
genome=/labs/Nyholm/HISAT2_Genome_V2/Es_genome

for FILE in ${FASTQ_FILE}
	do
		#This sets prefixes for in/outs downstream# 
		prefix=`echo ${FILE} | cut -d "/" -f 9 | cut -d "." -f 1`

		hisat2 -p 8 -x ${genome} -1 ${CLEANED_READS}${prefix}_FP.fq.gz -2 ${CLEANED_READS}${prefix}_RP.fq.gz \
		-S ${out}${prefix}.sam --rg PL:ILLUMINA --rna-strandness FR
	
done

