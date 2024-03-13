#!/bin/bash
#SBATCH --job-name=Ex90
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 10
#SBATCH --mem=100G
#SBATCH --partition=himem
#SBATCH --qos=himem
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nidhi.vijayan@uconn.edu
#SBATCH -o ex90_%j.out
#SBATCH -e ex90_%j.err

module load trinity/2.8.5
module load kallisto/0.44.0

#/isg/shared/apps/trinity/2.8.5/util/misc/contig_ExN50_statistic.pl \
#kallisto/kallisto2.isoform.TMM.EXPR.matrix centroids.fasta | tee ExN50.stats

/isg/shared/apps/trinity/2.8.5/util/misc/plot_ExN50_statistic.Rscript  ExN50.stats \
xpdf ExN50.stats.plot.pdf