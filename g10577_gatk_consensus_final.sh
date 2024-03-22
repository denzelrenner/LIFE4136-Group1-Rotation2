#!/bin/bash 

## setup 
# make output directory 
mkdir -p ~/g10577_gatk_output/diploid_output
mkdir -p ~/g10577_gatk_output/tetraploid_output


## index vcf files
# for diploids
gatk IndexFeatureFile \
	-I ~/VCFs/UK_scan_dips.vcf

# for tetraploids
gatk IndexFeatureFile \
	-I ~/VCFs/UK_scan_tets.vcf

## prodce dictionary file for reference fasta
gatk CreateSequenceDictionary \
	-R ~/Reference_Genome/C_excelsa_V5.fasta

## filter the vcf to only contain biallelic variants we are interested in
# for diploids
gatk SelectVariants \
	-R ~/Reference_Genome/C_excelsa_V5.fasta \
	-V ~/VCFs/UK_scan_dips.vcf \
	-L Cexcelsa_scaf_1:3432326-3436961 \
	-select-type SNP \
	-select-type INDEL \
	--select "AF > 0.49" \
	--output ~/g10577_gatk_output/diploid_output/diploid_filtered.vcf \
	--restrict-alleles-to BIALLELIC

# for tetraploids
gatk SelectVariants \
	-R ~/Reference_Genome/C_excelsa_V5.fasta \
	-V ~/VCFs/UK_scan_tets.vcf \
	-L Cexcelsa_scaf_1:3432326-3436961 \
	--output ~/g10577_gatk_output/tetraploid_output/tetraploid_filtered.vcf \
	-select-type SNP \
	-select-type INDEL \
	--select "AF > 0.49" \
	--restrict-alleles-to BIALLELIC
	
## index the new filtered vcf files
# for diploids
gatk IndexFeatureFile \
	-I ~/g10577_gatk_output/diploid_output/diploid_filtered.vcf

# for tetraploids
gatk IndexFeatureFile \
	-I ~/g10577_gatk_output/tetraploid_output/tetraploid_filtered.vcf


## produce consensus sequences
# for diploids
# exon1
gatk FastaAlternateReferenceMaker \
	-R ~/Reference_Genome/C_excelsa_V5.fasta \
	-O ~/g10577_gatk_output/diploid_output/g10577_diploid_exon1.fasta \
	-V ~/g10577_gatk_output/diploid_output/diploid_filtered.vcf \
	-L Cexcelsa_scaf_1:3432326-3436300

# exon2
gatk FastaAlternateReferenceMaker \
        -R ~/Reference_Genome/C_excelsa_V5.fasta \
        -O ~/g10577_gatk_output/diploid_output/g10577_diploid_exon2.fasta \
        -V ~/g10577_gatk_output/diploid_output/diploid_filtered.vcf \
        -L Cexcelsa_scaf_1:3436373-3436406

# exon3
gatk FastaAlternateReferenceMaker \
        -R ~/Reference_Genome/C_excelsa_V5.fasta \
        -O ~/g10577_gatk_output/diploid_output/g10577_diploid_exon3.fasta \
        -V ~/g10577_gatk_output/diploid_output/diploid_filtered.vcf \
        -L Cexcelsa_scaf_1:3436660-3436961

# concatenate the exons to make the final coding sequence
cat ~/g10577_gatk_output/diploid_output/*diploid_exon*.fasta > ~/g10577_gatk_output/diploid_output/g10577_coding_sequence_preliminary.fasta

# Remove the fasta headers and get only the nucleotide information for the coding sequence
grep -v '>' ~/g10577_gatk_output/diploid_output/g10577_coding_sequence_preliminary.fasta > \
~/g10577_gatk_output/diploid_output/g10577_diploid_coding_sequence.fasta

# remove all the exons
rm ~/g10577_gatk_output/diploid_output/*diploid_exon*
rm ~/g10577_gatk_output/diploid_output/g10577_coding_sequence_preliminary.fasta


# for tetraploids
# exon1
gatk FastaAlternateReferenceMaker \
        -R ~/Reference_Genome/C_excelsa_V5.fasta \
        -O ~/g10577_gatk_output/tetraploid_output/g10577_tetraploid_exon1.fasta \
        -V ~/g10577_gatk_output/tetraploid_output/tetraploid_filtered.vcf \
        -L Cexcelsa_scaf_1:3432326-3436300

# exon2
gatk FastaAlternateReferenceMaker \
        -R ~/Reference_Genome/C_excelsa_V5.fasta \
        -O ~/g10577_gatk_output/tetraploid_output/g10577_tetraploid_exon2.fasta \
        -V ~/g10577_gatk_output/tetraploid_output/tetraploid_filtered.vcf \
        -L Cexcelsa_scaf_1:3436373-3436406

# exon3
gatk FastaAlternateReferenceMaker \
        -R ~/Reference_Genome/C_excelsa_V5.fasta \
        -O ~/g10577_gatk_output/tetraploid_output/g10577_tetraploid_exon3.fasta \
        -V ~/g10577_gatk_output/tetraploid_output/tetraploid_filtered.vcf \
        -L Cexcelsa_scaf_1:3436660-3436961

# concatenate the exons to make the final coding sequence
cat ~/g10577_gatk_output/tetraploid_output/*tetraploid_exon*.fasta > ~/g10577_gatk_output/tetraploid_output/g10577_coding_sequence_preliminary.fasta

# Remove the fasta headers and get only the nucleotide information for the coding sequence
grep -v '>' ~/g10577_gatk_output/tetraploid_output/g10577_coding_sequence_preliminary.fasta > \
~/g10577_gatk_output/tetraploid_output/g10577_tetraploid_coding_sequence.fasta

# remove all the exons
rm ~/g10577_gatk_output/tetraploid_output/*tetraploid_exon*      
rm ~/g10577_gatk_output/tetraploid_output/g10577_coding_sequence_preliminary.fasta


# Now produce a file with all the multiallelic sites so we can manually choose a site with high frequency
# for diploids
gatk SelectVariants \
        -R ~/Reference_Genome/C_excelsa_V5.fasta \
        -V ~/VCFs/UK_scan_dips.vcf \
        -L Cexcelsa_scaf_1:3432326-3436961 \
        --output ~/g10577_gatk_output/diploid_output/multiallelic_diploid_filtered.vcf \
        --restrict-alleles-to MULTIALLELIC


# for tetraploids
gatk SelectVariants \
        -R ~/Reference_Genome/C_excelsa_V5.fasta \
        -V ~/VCFs/UK_scan_tets.vcf \
        -L Cexcelsa_scaf_1:3432326-3436961 \
        --output ~/g10577_gatk_output/tetraploid_output/multiallelic_tetraploid_filtered.vcf \
        --restrict-alleles-to MULTIALLELIC


# The vcf with multiallelic sites is very messy and we only want to have to look through variants that have an AF recorded
# for diploid
grep 'AF' ~/g10577_gatk_output/diploid_output/multiallelic_diploid_filtered.vcf > \
~/g10577_gatk_output/diploid_output/final_multiallelic_diploid.vcf

rm ~/g10577_gatk_output/diploid_output/multiallelic_diploid_filtered.vcf

# for tetraploid
grep 'AF' ~/g10577_gatk_output/tetraploid_output/multiallelic_tetraploid_filtered.vcf > \
~/g10577_gatk_output/tetraploid_output/final_multiallelic_tetraploid.vcf

rm ~/g10577_gatk_output/tetraploid_output/multiallelic_tetraploid_filtered.vcf

