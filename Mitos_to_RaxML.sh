#!/usr/bin/env bash
#This takes the Mitos annotated fasta file and makes into individual aligned gene files and one concatenated file with the genes in the same order for phylogegentic reconstruction and the partition file for RaxML
# needs MAFFT in the PATH
#needs a file with the sample names without extension. Same name as in the Mitos output. (eg. if the Mitos output is Contig1.fas, it should read Contig1)  This file should be called "contig_list.txt"


#Preping the mitos output so the actual sequences are in one line
#this is the mitos annotated fasta. Only Mitos annotated fasta should be in the folder, and should end in .fas  
for input in `ls *fas`
        do
          	awk  '{if(/[A|T|C|G|a|t|c|g]$/){printf "%s",$0}else{print}}' ${input} | sed 's/>/\n>/g' > ${input}_seq.fasta

        for gene in `egrep '[+,-]\; (\w+.+)' ${input}_seq.fasta | awk '{print $4}'`
                do
                grep ${gene} -A 1 ${input}_seq.fasta >> ${gene}_all.fas
                mafft --auto ${gene}_all.fas > ${gene}_all_align.fas
        done
done

#making all the alignments sequential
for align in `ls *all_align.fas`
	do 
	awk  '{if(/[A|T|C|G|a|t|c|g]$/){printf "%s",$0}else{print}}' ${align} | sed 's/>/\n>/g' > ${align}_seq.fasta
done


#creating a single sequence for each mitocondrial genome, that will be aligned to the others mt genomes and will be in the order of the partition file. Needs a list with the sample names without extension. Will call it contig_list.txt.

for contig in $(cat contig_list.txt) 
do
for gene_align in `ls *fas_seq.fasta`
	do
	grep -A1  -e ${contig} ${gene_align} |grep -v '>' |awk  '{if(/[A|T|C|G|a|t|c|g]$/){printf "%s",$0}else{print}}' >> ${contig}_for_RaxMl.fasta
		done
	done


# giving the name to the sequence
for genome in `ls *_for_RaxMl.fasta`
do
 sed -i -e "1i >${genome}" ${genome}
 done

#creating the partition file
for contig in $(cat contig_list.txt) 
	do 
for gene_align in `ls *fas_seq.fasta`
	 do grep -A1 -e ${contig} ${gene_align} |grep -v '>'|tr "\n" "#" |grep -aob '#' |awk -F: '{print $1}'|awk 'p{printf "DNA, gene%d=%d-%d\n", ++r, p, $1-1} {p=$1}' |sort -V|uniq >> Partition.txt
	  done
	done



#creating the concatenated file with all aligned genomes

for genome in `ls *_for_RaxMl.fasta`
	do
	(cat "${genome}"; echo) >> all_genomes_for_RaxML.fasta
done

