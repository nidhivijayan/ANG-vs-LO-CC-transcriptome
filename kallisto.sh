module load kallisto/0.44.0

#First index the fasta file
kallisto index -i trinity.fasta.index Trinity.fasta

#Then quantify with kallisto
kallisto quant -i trinity.fasta.index -o kallisto_quant_CCQ3 -t 8 ${input}/167753_CCQ3_FP.fastq.gz ${input}/167753_CCQ3_RP.fastq.gz 

for FILE in ${input}
	do
#		#This sets prefixes for in/outs downstream# 
	prefix=`echo ${FILE} | cut -d "/" -f 9 | cut -d "." -f 1`
	
		kallisto quant -i trinity.fasta.index -o kallisto_quant_${prefix} -t 8 ${FILE}/${prefix}_FP.fq.gz ${FILE}/${prefix}_RP.fq.gz 
done
