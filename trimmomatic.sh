#!/bin/bash
#SBATCH --job-name=trim
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 20
#SBATCH --mem=350G
#SBATCH --partition=himem
#SBATCH --qos=himem
#SBATCH --mail-type=ALL
#SBATCH -o trim_%j.out
#SBATCH -e trim_%j.err

module load Trimmomatic/0.36

for FILE in ${input}
	do
		prefix=`echo ${FILE} | cut -d "/" -f 8 | cut -d "." -f 1`

		java -jar $Trimmomatic PE ${FILE}/${prefix}_L001_R1_001.fastq.gz ${FILE}/${prefix}_L001_R2_001.fastq.gz ${TRIMMED}${prefix}_FP.fq.gz ${TRIMMED}${prefix}_FuP.fq.gz ${TRIMMED}${prefix}_RP.fq.gz ${TRIMMED}${prefix}_RuP.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:149
done
