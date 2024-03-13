module load Trinotate/3.2.1
module load signalp/4.1

hmmscan --cpu 30 --domtblout pfam.domtblout /isg/shared/databases/Pfam/Pfam-A.hmm Cat_merged_all.fasta.transdecoder_dir/longest_orfs.pep

# BLASTp annotation
module load blast/2.7.1

db=~/transcriptome
input=~/TRANSDECODER

makeblastdb -in uniprot_sprot.pep -dbtype prot

blastp -query transcripts_AA.fasta -db ${db}/uniprot_sprot.pep \
-max_target_seqs 1 -outfmt 6 -evalue 1e-5 \
> swissprot.blastp.outfmt6

Trinotate Trinotate.sqlite init \
     --gene_trans_map ~/TRANSDECODER/centroids_w_gills_ed_gene.trans_map \
     --transcript_fasta ~/TRANSDECODER/centroids_w_gills_ed8.fasta \
     --transdecoder_pep ~/TRINITY_V2/TRANSDECODER/Trinity_combined_ang_cc_gills.fasta.transdecoder.pep

Trinotate Trinotate.sqlite LOAD_pfam ../TrinotatePFAM2.out

Trinotate Trinotate2.sqlite LOAD_swissprot_blastp ../BLAST/swissprot.blastp.outfmt6

Trinotate Trinotate.sqlite report > Trinotate_v2.xls

Trinotate --db <sqlite.db> --init \
           --gene_trans_map <file> \
           --transcript_fasta <file> \
           --transdecoder_pep <file>

