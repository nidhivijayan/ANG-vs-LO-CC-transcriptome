module load trinity/2.8.5
module load kallisto/0.44.0
module load R/3.6.1
Rscript -e "library(edgeR)"

/isg/shared/apps/trinity/2.8.5/util/abundance_estimates_to_matrix.pl --est_method kallisto \
	--gene_trans_map trinity.fasta.gene_trans_map \
	--out_prefix kallisto \
	--name_sample_by_basedir \
#  *output of kallisto/abundance.tsv of all samples
