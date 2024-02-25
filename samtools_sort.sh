module load samtools/1.7

for FILE in ${input}
	do
		prefix=`echo ${FILE} | cut -d "/" -f 9 | cut -d "." -f 2`
samtools sort ${input}/167756_ANGS9Aligned.sortedByCoord.out.bam -o ${output}/167756_ANGS9.out.bam

done
