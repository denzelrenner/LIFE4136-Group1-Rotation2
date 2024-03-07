# LIFE4136_Group1_Rotation2
This repositry contained all the Scripts and data needed to reproduce results from group1 in rotation2 of the groupwork projects.

## What is the problem we have been presented with?
We have been provided with some genes in cochaleria 

## What are the aims of our analysis?

## What are the expected outcomes?
We expect to 

# Prerequisites
## Required files
All these files should be downloaded into the same directory before following the rest of this document.

For replicating results on the cloud :
1)VCFs directory containing {UK_scan_dips.vcf and UK_scan_tets.vcf}
2)Reference_Genome directory containing {C_excelsa_V5_braker2_wRseq.gff3, C_excelsa_V5.fasta and C_excelsa_V5.fasta.fai}

For replicating results from your local machine
1)

## Tool versions and links
These are all the tools that were used in our analysis with versions and links provided where applicable.

ncbi blast (https://blast.ncbi.nlm.nih.gov/)
orf finder (https://www.ncbi.nlm.nih.gov/orffinder/)
uniprot (https://www.uniprot.org/align), clustalO version 1.2.4 (clustal was ran through the align function on the uniprot website)
jalview version 2.11.3.2
alphafold collab (ColabFold v1.5.5: AlphaFold2 using MMseqs2
(https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb)
PyMOL version 2.5.8
Homebrew version 4.2.10 (Mac Intel i5)
ffmpeg version 6.1.1 
GATK verison 4.2.2.0, (HTSJDK version 2.24.1, Picard version 2.25.4)
samtools version 1.19.2, (htslib version 1.19.1)
RAxML-NG version 0.9.0

## Tool intallation 
These scripts/commands should be ran to install all the tools necessary to reproduce results

On the cloud:
raxml_install.sh


On your local machine:
brew install ffmpeg (Mac Intel i5)

# THE BEGINNING 

## Consensus Sequences (Nucleotide)
The first thing we will do is get consensus sequences for our diploids and tetraploids. We will filter the vcfs we have to only include biallelic variants with an allele frequency greater than 0.4. The resulting filtered vcfs will then be used, along with the reference fasta and gff, to produce our consensus sequences for the genes of interest only and not the entire scaffold. This is accomplished by running the script below

<insert_script>

## Protein Sequences 
Amongst the resulting files there should be fasta files containing the entire consensus gene sequence for g46214 and g10577 in our diploids and tetraploids. We then followed the link to the orf finder website and input the consensus sequences for each gene in the different contrasts.

1)Input nucleotide sequence from g42614 diploid into the query sequence field and submit the job. From the resulting output stitch the protein sequence for the three longest open reading frames to form the full length protein(ORF1+ORF3+ORF4) (results can be achieved if only the coding sequence is taken and not the entire gene. Methodology not shown) 

2)Input nucleotide sequence from g42614 tetraploid into the query sequence field and submit the job. From the resulting output stitch the protein sequence for the three longest open reading frames to form the full length protein(ORF1+ORF3+ORF4) (results can be achieved if only the coding sequence is taken and not the entire gene. Methodology not shown) 

3)Input nucleotide sequence from g10577 diploid into the query sequence field and submit the job. From the resulting output take the longest open reading frame which is the full length protein (ORF1). If you take only coding sequence is the tetraploid truncated?

4)Input nucleotide sequence from g10577 tetraploid into the query sequence field and submit the job. From the resulting output take the longest open reading frame which is the full length protein (ORF1).












# Tool References
These are references (where applicable) for all the tools used in our analysis
