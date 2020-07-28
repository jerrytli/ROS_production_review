makeblastdb -in Housekeeping_gene_list.fasta -parse_seqids -dbtype 'nucl' -out Housekeeping_genes_blast_db

blastp -query NADPH-oxidase-cTerm-reference.fasta -db ./BLAST_database/NADPH-oxidase-homologs-db -evalue 1e-8 -outfmt "6 qseqid sseqid pident evalue len qstart qend sseq" -out file.txt
