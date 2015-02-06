#!/bin/bash

if [ -z $1 ] || ! [  -a $1 ]; then
	echo "Error: Please provide the input fasta file!"
	exit
fi

if ! [ -a $2 ]; then
	echo "Error: Please provide the trinotate database!"
	exit
fi
$TRINITY_HOME/trinity-plugins/TransDecoder_r20140704/TransDecoder -t $1 2>>trinotAdd.log

blastx -query $1 -db /home/GLBRCORG/omoskvin/dbs/uniprot_sprot.fasta -num_threads 8 -max_target_seqs 1 -outfmt 6 > blastx.outfmt6  2>>trinotAdd.log &

blastp -query "$1.transdecoder.pep" -db /home/GLBRCORG/omoskvin/dbs/uniprot_sprot.fasta -num_threads 8 -max_target_seqs 1 -outfmt 6 > blastp.outfmt6 2>>trinotAdd.log &

hmmscan --cpu 8 --domtblout NVS.unmapped.PFAM.out /home/GLBRCORG/omoskvin/dbs/Pfam-A.hmm "$1.transdecoder.pep" > pfam.log 2>>trinotAdd.log &

wait

#We are using our own mapping
#$TRINITY_HOME/util/support_scripts/get_Trinity_gene_to_trans_map.pl $1 >  Trinity.fasta.gene_trans_map 2>>trinotAdd.log

./fake_map.pl $1 >  Trinity.fasta.gene_trans_map

$TRINOTATE_HOME/Trinotate $2 init --gene_trans_map Trinity.fasta.gene_trans_map --transcript_fasta $1 --transdecoder_pep "$1.transdecoder.pep" 2>>trinotAdd.log

$TRINOTATE_HOME/Trinotate $2 LOAD_blastp blastp.outfmt6 2>>trinotAdd.log

$TRINOTATE_HOME/Trinotate $2 LOAD_blastx blastx.outfmt6 2>>trinotAdd.log

$TRINOTATE_HOME/Trinotate $2 report > trinotate_annotation_report.xls 2>>trinotAdd.log
