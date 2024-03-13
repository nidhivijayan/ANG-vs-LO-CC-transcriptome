#!/bin/bash
#SBATCH --job-name=GO
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=200G
#SBATCH --partition=himem
#SBATCH --qos=himem
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nidhi.vijayan@uconn.edu
#SBATCH -o go_%j.out
#SBATCH -e go_%j.err

module load trinity/2.8.5
module load R/4.2.0
Rscript -e "library(edgeR)"
Rscript -e "library(limma)"
Rscript -e "library(DESeq2)"
#Rscript -e "library(ctc)"
Rscript -e "library(Biobase)"
Rscript -e "library(gplots)"
#Rscript -e "library(ape)"
#Rscript -e "library(goseq)"

GO_path=/home/FCAM/nvijayan/transcriptome/ANG_CC/bam/merged_GG_trinity/trinotate/

/isg/shared/apps/trinity/2.8.5/Analysis/DifferentialExpression/run_DE_analysis.pl \
	--matrix ../kallisto_ed/kallisto_ed.gene.counts.matrix \
	--method edgeR \
	--samples_file samples_sym.txt \
	--min_reps_min_cpm 15,10 \
	--output DE_edgR_sym \
	--contrasts contrasts_2.txt

/isg/shared/apps/trinity/2.8.5/Analysis/DifferentialExpression/analyze_diff_expr.pl \
	--matrix ../../kallisto_ed/kallisto_ed.gene.TMM.EXPR.matrix --samples samples_sym.txt\
	--examine_GO_enrichment --GO_annots ${GO_path}/go_annotations.txt --gene_lengths Trinity.gene_lengths

#/isg/shared/apps/trinity/2.8.5/Analysis/DifferentialExpression/run_GOseq.pl \
#	--genes_single_factor kallisto_ed.gene.counts.matrix.ANG_vs_CC.edgeR.DE_results.P0.001_C2.ANG-UP.subset \
#	--GO_assignments ${GO_path}/go_annotations.txt \
#	--lengths Trinity.gene_lengths.txt 