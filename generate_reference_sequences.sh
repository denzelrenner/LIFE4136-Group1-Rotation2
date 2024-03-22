#!/bin/bash

# make directories to put the reference sequences into
mkdir -p ~/reference_sequences

# Using the genomic coordinates in the gff, get the coding sequence for our genes of interest
# g46214
# exon1
samtools faidx ~/Reference_Genome/C_excelsa_V5.fasta Cexcelsa_scaf_2:36768726-36768935 > \
~/reference_sequences/g46214_exon1.fasta

# exon2
samtools faidx ~/Reference_Genome/C_excelsa_V5.fasta Cexcelsa_scaf_2:36769070-36769459 > \
~/reference_sequences/g46214_exon2.fasta

# exon3
samtools faidx ~/Reference_Genome/C_excelsa_V5.fasta Cexcelsa_scaf_2:36769937-36770296 > \
~/reference_sequences/g46214_exon3.fasta

# concatenate the exons to make the final coding sequence
cat ~/reference_sequences/*_exon*.fasta > ~/reference_sequences/g46214_coding_sequence_preliminary.fasta

# Remove the fasta headers and get only the nucleotide information for the coding sequence
grep -v '>' ~/reference_sequences/g46214_coding_sequence_preliminary.fasta > \
~/reference_sequences/g46214_coding_sequence.fasta

# remove all the exons
rm ~/reference_sequences/*_exon*.fasta
rm ~/reference_sequences/g46214_coding_sequence_preliminary.fasta


# Using the genomic coordinates in the gff, get the coding sequence for our genes of interest
# g10577
# exon1
samtools faidx ~/Reference_Genome/C_excelsa_V5.fasta Cexcelsa_scaf_1:3432326-3436300 > \
~/reference_sequences/g10577_exon1.fasta

# exon2
samtools faidx ~/Reference_Genome/C_excelsa_V5.fasta Cexcelsa_scaf_1:3436373-3436406 > \
~/reference_sequences/g10577_exon2.fasta

# exon3
samtools faidx ~/Reference_Genome/C_excelsa_V5.fasta Cexcelsa_scaf_1:3436660-3436961 > \
~/reference_sequences/g10577_exon3.fasta

# concatenate the exons to make the final coding sequence
cat ~/reference_sequences/*_exon*.fasta > ~/reference_sequences/g10577_coding_sequence_preliminary.fasta

# Remove the fasta headers and get only the nucleotide information for the coding sequence
grep -v '>' ~/reference_sequences/g10577_coding_sequence_preliminary.fasta > \
~/reference_sequences/g10577_coding_sequence.fasta

# remove all the exons
rm ~/reference_sequences/*_exon*.fasta
rm ~/reference_sequences/g10577_coding_sequence_preliminary.fasta
