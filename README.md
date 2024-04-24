# LIFE4136_Group1_Rotation2
This repository contains all the scripts and lists of data files needed to reproduce results from group1 in rotation2 of the groupwork projects.

Note: Unless on the cloud or stated otherwise, all command line code was ran on a Mac Intel i5 

## What is the problem we have been presented with?
We have been provided with some genes in Cochaleria that are under selection during the stabalisation of polyploidy. Polyploids are individuals that arise from a whole genome duplication (WGD) event, a process that results in multiple sets of chromosomes within an individual. These polyploids evolve mechanisms to handle the hardships that come with WGD, such as meiotic instability or regulation of gene expression, and also develop adaptations to tolerate the harsh environmental conditions that can sometimes lead to polyploidy. In this project we are focusing on two genes caught in the selection scan, g46214 and g10577. For each gene, we are to investigate the amino acid sequences and tetriary structures of their proteins in diploids and tetraploids to determine the consequence of any genetic variants.

## What are the aims of our analysis?
We were aiming to translate the g46214 and g10577 genes and investigate the primary sequences of their proteins in diploids and tetraploids to first identify important domains and predict the function of the protein. We also aimed to identify any mutations that might result in altered function of the protein in diploids or tetraploids by comparing the primary sequences with that of close homologs. Another objective of the analysis was to generate structural predictions for the diploid and tetraploid alleles for each of our candidate genes which would allow us to compare the tetriary structures of the proteins at the different ploidy levels, and observe the impact of the different mutations on protein folding and other protein properties like hydrophobicity or electrostatic charge. 

## What are the expected outcomes?
Depending on the gene and given our understanding of the system, we expect to see mutations and amino acid changes in the tetraploid g46214 and g10577 proteins, relative to the diploid g46214 and g10577 proteins, that help stabilise the polyploids and confer a fitness advantage through enhanced ecological tolerance for example. The altered function of the proteins could result from the mutation of a conserved, functional residue which enhances or reduces the interaction with other ions and proteins, whilst in some cases a mutation can cause a truncation and remove an entire functional domain.

# Prerequisites
## Required files and data
All these files should be downloaded into the same directory before following the rest of this document.

The vcf files contain variants at only the genomic regions we are interested in. The `UK_scan_dips.vcf` contains variants in diploids and `UK_scan_tets.vcf` contains variants in tetraploids.

The `C_excelsa_V5_braker2_wRseq.gff3` contains annotations for the reference genome (Cochlearia excelsa). The `C_excelsa_V5.fasta` file contains sequence information for the reference genome (Cochlearia excelsa) and has been indexed to produce a .fai file. The `C_excelsa_V5.fasta.fai` file is a fasta index file which allows us to find a particular nucelotide at specific genomic coordinates in the `C_excelsa_V5.fasta` file. You can read more about it on the [GATK website](https://gatk.broadinstitute.org/hc/en-us/articles/360035531652-FASTA-Reference-genome-format).

|Directory|Files|
|---------|-----|
|VCFs|UK_scan_dips.vcf<br>UK_scan_tets.vcf|
|Reference_Genome|C_excelsa_V5_braker2_wRseq.gff3<br>C_excelsa_V5.fasta<br>C_excelsa_V5.fasta.fai|


## Tool versions and links
These are all the tools that were used in our analysis with versions and links provided where applicable. Dependencies for certain packages, and their versions, are placed in parentheses. Some references were chosen based on what was recommended on the tool's online help page/documentation.

| Tool | Version | Reference(Harvard Style) |
|------|---------|-----------|
|[ncbi blast](https://blast.ncbi.nlm.nih.gov/)|NA| Johnson, M., Zaretskaya, I., Raytselis, Y., Merezhuk, Y., McGinnis, S. and Madden, T.L., 2008. NCBI BLAST: a better web interface. Nucleic acids research, 36(suppl_2), pp.W5-W9.|
|[TAIR BLAST](https://www.arabidopsis.org/Blast/index.jsp)|version 2.9.0+ (uses NCBI BLAST 2.9.0+)| Garcia-Hernandez, M., Berardini, T., Chen, G., Crist, D., Doyle, A., Huala, E., Knee, E., Lambrecht, M., Miller, N., Mueller, L.A. and Mundodi, S., 2002. TAIR: a resource for integrated Arabidopsis data. Functional & integrative genomics, 2, pp.239-253.|
|[orf finder](https://www.ncbi.nlm.nih.gov/orffinder/)|NA| NA|
|[uniprot](https://www.uniprot.org/align)|clustalO version 1.2.4 (clustal was ran through the align function on the uniprot website)|UniProt Consortium, 2019. UniProt: a worldwide hub of protein knowledge. Nucleic acids research, 47(D1), pp.D506-D515.|
|[alphafold collab](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb)|(ColabFold v1.5.5: AlphaFold2 using MMseqs2|Mirdita, M., Schütze, K., Moriwaki, Y., Heo, L., Ovchinnikov, S. and Steinegger, M., 2022. ColabFold: making protein folding accessible to all. Nature methods, 19(6), pp.679-682.|
|[InterPro](https://www.ebi.ac.uk/interpro/)|version 98.0|Hunter, S., Apweiler, R., Attwood, T.K., Bairoch, A., Bateman, A., Binns, D., Bork, P., Das, U., Daugherty, L., Duquenne, L. and Finn, R.D., 2009. InterPro: the integrative protein signature database. Nucleic acids research, 37(suppl_1), pp.D211-D215.|
|[ITOL:Interactive Tree Of Life](https://itol.embl.de/)|version 6.8.2|Letunic, I. and Bork, P., 2021. Interactive Tree Of Life (iTOL) v5: an online tool for phylogenetic tree display and annotation. Nucleic acids research, 49(W1), pp.W293-W296.|
|[MEGA](https://www.megasoftware.net/)|version 11|Tamura, K., Dudley, J., Nei, M. and Kumar, S., 2007. MEGA4: molecular evolutionary genetics analysis (MEGA) software version 4.0. Molecular biology and evolution, 24(8), pp.1596-1599.|
|[jalview](https://www.jalview.org/)|version 2.11.3.2|Waterhouse, A.M., Procter, J.B., Martin, D.M., Clamp, M. and Barton, G.J., 2009. Jalview Version 2—a multiple sequence alignment editor and analysis workbench. Bioinformatics, 25(9), pp.1189-1191.|
|[PyMOL](https://pymol.org/)|version 2.5.8|The PyMOL Molecular Graphics System, Version 2.0 Schrödinger, LLC.|
|[GATK](https://github.com/broadinstitute/gatk) (HTSJDK,Picard)|version 4.2.2.0 (version 2.24.1,version 2.25.4)|Van der Auwera GA & O'Connor BD. (2020). Genomics in the Cloud: Using Docker, GATK, and WDL in Terra (1st Edition). O'Reilly Media.|
|[samtools](https://github.com/samtools/samtools/blob/develop/README.md) (htslib)|version 1.19.2 (version 1.19.1)|Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., Durbin, R. and 1000 Genome Project Data Processing Subgroup, 2009. The sequence alignment/map format and SAMtools. bioinformatics, 25(16), pp.2078-2079.|
|[RAxML-NG](https://github.com/amkozlov/raxml-ng)|version 0.9.0|Kozlov, A.M., Darriba, D., Flouri, T., Morel, B. and Stamatakis, A., 2019. RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference. Bioinformatics, 35(21), pp.4453-4455.|
|[SWISS-MODEL](https://swissmodel.expasy.org/interactive)|NA|Waterhouse, A., Bertoni, M., Bienert, S., Studer, G., Tauriello, G., Gumienny, R., Heer, F.T., de Beer, T.A.P., Rempfer, C., Bordoli, L. and Lepore, R., 2018. SWISS-MODEL: homology modelling of protein structures and complexes. Nucleic acids research, 46(W1), pp.W296-W303.|
|[Conda](https://github.com/conda/conda)|version 23.5.2|NA|
|[Homebrew](https://brew.sh/)|version 4.2.10|NA|
|ffmpeg|version 6.1.1|NA|

## Tool intallation 
The guidance below outlines the scripts,steps or commands that have to be ran to install some of the tools needed to reproduce our results.

### On the cloud HPC:
To install raxml-ng run the script below in the command line.
```bash
bash ~/raxml_install.sh
```

To install [GATK](https://github.com/broadinstitute/gatk), [samtools](https://github.com/samtools/samtools/blob/develop/README.md), and [Conda](https://github.com/conda/conda) follow the guidance on their respective github pages as those tools were installed by a bioinformatics technician and code to install them was not generated from scratch in our analysis. 


### On your local machine:
To install the `homebrew` package manager follow the steps outlined on the [homebrew website](https://brew.sh/)

To install the video formatting tool `ffmpeg` input the command below into the command line once `homebrew` has been installed

```bash
brew install ffmpeg
``` 

To install the [Jalview](https://www.jalview.org/),[PyMOL](https://pymol.org/),and [MEGA](https://www.megasoftware.net/) softwares, navigate to their respective websites and follow the download guidance for your machine.

## Script description

| Script Name | Description |
|-------------|-------------|
| `raxml_install.sh` | installs raxml-ng |
| `phylogenetic_tree_g46214.sh` | runs raxml-ng to build a maximum likelihood tree using the allignment for g46214 proteins in tetraploid and diploid and their homologs |
| `generate_reference_sequences.sh` | produces a nucleotide sequence (in fasta format) for the g10577 and g46214 reference genes |
| `g46214_gatk_consensus_final.sh` | produces a consensus nucleotide sequence (in fasta format) for the g46214 diploid and tetraploid genes |
| `g10577_gatk_consensus_final.sh` | produces a consensus nucleotide sequence (in fasta format) for the g10577 diploid and tetraploid genes |
| `colorh.py` | colours proteins by hydrophobicity. This was retrieved from the [PyMOL wiki color h page](https://pymolwiki.org/index.php/Color_h) |
| `tetraploid_domain_highlight.py`| produces a series of images (rotated by 90 degrees) of the tetraploid g46214 protein with its domains coloured and highlighted, the electrostatic potential of the whole protein, and the hydrophobicity of the whole protein |
| `diploid_domain_highlight.py`| produces a series of images (rotated by 90 degrees) of the diploid g46214 protein with its domains coloured and highlighted, the electrostatic potential of the whole protein, and the hydrophobicity of the whole protein |
| `tetraploid_temporary_image_generation.py` | produces a series of images rotated by 1 degree across 360 degrees for the tetraploid g46214 protein |
| `diploid_temporary_image_generation.py` | prdocues a series of images rotated by 1 degree across 360 degrees for the diploid g46214 protein |
| `make_protein_rotation_movie.sh` | stitches together the 360 diploid and tetraploid g46214 protein images to make a short movie of the protein rotating and showing all domains/faces of the protein |


# THE ANALYSIS

## Part1 - Consensus Sequences (Nucleotide)
All code in this section needs to be ran on the cloud HPC. 

Running all the scripts in this `consensus sequence(Nucleotide)` section will also produce `C_excelsa_V5.dict` and `C_excelsa_V5.fasta.fai` files which are required by GATK tools to access specified regions of the reference fasta. The `.dict` file describes the contents of our fasta file, and as mentioned before the `.fai` file is a fasta index file which allows us to find a particular nucelotide at specific genomic coordinates in the FASTA file. You can read more about these file formats on the [GATK website](https://gatk.broadinstitute.org/hc/en-us/articles/360035531652-FASTA-Reference-genome-format)

The first thing we will do is get consensus sequences for the reference g46214 and g10577 genes so it can then be used in later stages of the analysis. We knew the genomic positions for the different genes of interest by opening the `C_excelsa_V5_braker2_wRseq.gff3` file in a text editor and manually searching for the gene names. To create the consensus sequences for the reference gene you should enter the command below into the command line.

```bash
conda activate /shared/apps/conda/bio2
bash ~/generate_reference_sequences.sh
```

This should produce a directory called `reference_sequences` and within that directory there should be a file called `g10577_coding_sequence.fasta` which has the g10577 reference consensus sequence and another file called `g46214_coding_sequence.fasta` which has the g46214 reference consensus sequence.

### g46214
With the `bio2` conda environemnt still active, we will now filter the vcfs we have to only include biallelic variants with an allele frequency greater than 0.49 (to include allele frequencies of 0.5). The resulting filtered vcfs will then be used, along with the reference fasta and gff, to produce our consensus sequences for the genes of interest in the diploid and tetraploid. This is accomplished by running the script below:

```bash
bash ~/g46214_gatk_consensus_final.sh
```

This script produces a directory called `g46214_gatk_output` with sub-directories for diploid output (`diploid_output`) and tetraploid output (`tetraploid_output`).

For g46214 diploids, the most important output files in the `diploid_output` directory are called `g46214_diploid_coding_sequence.fasta` which contains the nucleotide consensus sequence for the diploid, and `final_multiallelic_diploid.vcf` which contains multiallelic variants. You have to manually identify and extract alleles in the `final_multiallelic_diploid.vcf` file that have an allele frequency greater than 0.5, and then insert them into the `g46214_diploid_coding_sequence.fasta` file to get a final sequence for diploids. 

For g46214 tetraploids, the most important output files in the `tetraploid_output` directory are called `g46214_tetraploid_coding_sequence.fasta` which contains the nucleotide consensus sequence for the tetraploid, and `final_multiallelic_tetraploid.vcf` which contains multiallelic variants. You have to manually identify and extract alleles in the `final_multiallelic_tetraploid.vcf` file that have an allele frequency greater than 0.5, and then insert them into the `g46214_tetraploid_coding_sequence.fasta` file to get a final sequence for the tetraploids.

### g10577
With the `bio2` conda environemnt still active, we will filter the vcfs we have to only include biallelic variants with an allele frequency greater than 0.6. The resulting filtered vcfs will then be used, along with the reference fasta and gff, to produce our consensus sequences for the genes of interest in the tetraploid and diploid. This is accomplished by running the script below:

```bash
bash ~/g10577_gatk_consensus_final.sh
```
This script produces a directory called `g10577_gatk_output` with sub-directories for diploid output (`diploid_output`) and tetraploid output (`tetraploid_output`).

For g10577 diploids, the most important output files in the `diploid_output` directory are called `g10577_diploid_coding_sequence.fasta` which contains the nucleotide consensus sequence for the diploid, and `final_multiallelic_diploid.vcf` which contains multiallelic variants. You have to manually identify and extract alleles in the `final_multiallelic_diploid.vcf` file that have an allele frequency greater than 0.6, and then insert them into the `g10577_diploid_coding_sequence.fasta` file to get a final sequence for diploids.

For g10577 tetraploids, the most important output files in the `tetraploid_output` directory are called `g10577_tetraploid_coding_sequence.fasta` which contains the nucleotide consensus sequence for the diploid, and `final_multiallelic_tetraploid.vcf` which contains multiallelic variants. You have to manually identify and extract alleles in the `final_multiallelic_tetraploid.vcf` file that have an allele frequency greater than 0.6, and then insert them into the `g10577_tetraploid_coding_sequence.fasta` file to get a final sequence for tetraploids.

## Part2 - Protein Sequences 
Now that we have the nucleotide sequences for g46214 and g10577 in our diploids, tetraploids, and reference, we can translate them to get our protein sequences. This can be accomplished by following the steps below.

### g46214

 1. Follow the link to the [orf finder](https://www.ncbi.nlm.nih.gov/orffinder/) website.

 2. For the diploid protein, input the nucleotide seqeunce from the `g46214_diploid_coding_sequence.fasta` file into the query seqeunce field and submit the job. On the output page, copy and paste the amino acid seqeuence for ORF1 (the longest open reading frame) into a new fasta file with any identifiable headers you prefer (i.e `>diploid_g46214`).

 3. For the tetraploid protein, input the nucleotide seqeunce from the `g46214_tetraploid_coding_sequence.fasta` file into the query seqeunce field and submit the job. On the output page, copy and paste the amino acid seqeuence for ORF1 (the longest open reading frame) into a new fasta file with any identifiable headers you prefer (i.e `>tetraploid_g46214`).

 4. For the reference protein, input the nucleotide seqeunce from the `g46214_coding_sequence.fasta` file into the query seqeunce field and submit the job. On the output page, copy and paste the amino acid seqeuence for ORF1 (the longest open reading frame) into a new fasta file with any identifiable headers you prefer (i.e `>reference_g46214`).

### g10577

 1. Follow the link to the [orf finder](https://www.ncbi.nlm.nih.gov/orffinder/) website.

 2. For the diploid protein, input the nucleotide seqeunce from the `g10577_diploid_coding_sequence.fasta` file into the query seqeunce field and submit the job. On the output page, copy and paste the amino acid seqeuence for ORF1 (the longest open reading frame) into a new fasta file with any identifiable headers you prefer (i.e `>diploid_10577`).

 3. For the tetraploid protein, input the nucleotide seqeunce from the `g10577_tetraploid_coding_sequence.fasta` file into the query seqeunce field and submit the job. On the output page, copy and paste the amino acid seqeuence for ORF1 (the longest open reading frame) into a new fasta file with any identifiable headers you prefer (i.e `>tetraploid_g10577`).

 4. For the reference protein, input the nucleotide seqeunce from the `g10577_coding_sequence.fasta` file into the query seqeunce field and submit the job. On the output page, copy and paste the amino acid seqeuence for ORF1 (the longest open reading frame) into a new fasta file with any identifiable headers you prefer (i.e `>reference_g10577`).


## Part3 - Homolog identification
We have now retrieved the sequences for our proteins, so we want to figure out what they might actually be. One good way to do that is finding their closest homologs in other plant or animal species which can also give you an idea on how the protein may function.

### g46214
You should create a fasta file called `g46214_homologs.fasta` to store all homologs. Homologous proteins were selected based on having 100% query cover and >60% percentage identity to the reference protein sequence. 

To identify g46214 homologs:

 1. Follow the link to the [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/) webpage

 2. Select `Protein BLAST` 

 3. Input the g46214 reference protein sequence into the `Enter Query Sequence` field, and select `BLAST` at the bottom of the page

 4. Manually retrieve the protein sequence of the homologs from the blastp results page by selecting their ncbi dataset accession code (i.e `XP_018447019.1`).

 5. Select the `FASTA` option at the top of the page, and finally copy and paste the protein sequence into the `g46214_homologs.fasta` fasta file

 6. The protein sequences for the reference, diploid and tetraploid generated in Part2 of the analysis should also be added to the `g46214_homologs.fasta` file in fasta format.

The accession codes for the homologs identified in this step and used in subsequent steps of the analysis are outlined below:
| Homolog | Accession Code | 
|--------|------------------|
| Hirschfeldia incana | KAJ0255869.1 |
| Brassica rapa | XP_009106205.1 |
| Raphanus sativus | XP_018447019.1 |
| Brassica napus | XP_013648757.1 |
| Arabidopsis thaliana | NP_177686.1 |
| Arabidopsis suecica | KAG7659674.1 |
| Arabidopsis thaliana x Arabidopsis arenosa | KAG7651806.1 |
| Eutrema salsugineum | XP_006390294.1 |
| Capsella rubella | XP_006302550.1 |

The closest homologs of g46214 were bbx21 proteins in different plant species.

### g10577

For g10577, homologous proteins were selected based on >70% query cover and >40% percentage identity to the tetraploid protein sequence. You should create a fasta file called `g10577_homologs.fasta` to store all homologs. 

To identify g10577 homologs in other plant species we followed the steps below:

 1. Follow the link to the [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/) webpage

 2. Select `Protein BLAST` 

 3. Input the g10577 tetraploid protein sequence into the `Enter Query Sequence` field, and select `BLAST` at the bottom of the page

 4. Manually retrieve the protein sequence of the homologs from the blastp results page by selecting their NCBI dataset accession code (i.e `XP_018447019.1`) in the blastp results page

 5. Select the `FASTA` option at the top of the page, and finally copy and paste the protein sequence into the `g10577_homologs.fasta` fasta file

 6. The reference, diploid and tetraploid protein sequences generated in Part2 of the analysis should also be added to the `g10577_homologs.fasta` file in fasta format. Note that we also searched for homologs in the model species Arabidopsis thaliana on the TAIR website, but for g10577 the identified homolog did not make biological sense when investigated further through multiple sequence allignments and the literature cited in the `Domain Identification` section of this analysis so we did not include that protein as a homolog.

The accession codes for the homologs used in subsequent steps of the analysis are outlined below:
| Homolog | Accession Code | 
|--------|------------------|
|Microthlasmi erraticum|CAA7027704.1|
|A thaliana x A arenosa|KAG7542151.1|
|Trifolium repens|KAK2377062.1|
|Trifolium pratense|PNX95763.1|
|Vitis vinifera|RVW62496.1|
|Brassica rapa|XP_033145473.1|
|Lactuca sativa|XP_042756658.2|
|Manihot esculenta|XP_021629864.2|
|Brassica napa|XP_022570911.2|
|Hellanthus annus|XP_022003179.2|
|Sparassis crispa|XP_027612870.1|

The closest homologs for g10577 were retrotransposons.

## Part4 - Domain identification
We have discovered the closest homologs for our proteins in different species so we can begin investigating the structural domains in our proteins to get an idea of how they function.

### g46214

To verify the functional domains within the g46214 proteins we followed these steps:

 1. Follow the link to the [InterPro website](https://www.ebi.ac.uk/interpro/)

 2. Input the diploid protein sequence (in fasta format) for g46214 into the query field labelled `Enter your sequence` and choose `search`

 3. Manually insert the positions of the domains into a txt file called `g46214_diploid_domains.txt` so we have the exact coordinates of the different domains in the protein
 
 4. Follow the link in step 1 to return to the Interpro home page. Input the tetraploid protein sequence for g46214 into the query field labelled `Enter your sequence` and choose `search`

 5. Manually insert the positions of the domains into a txt file called `g46214_tetraploid_domains.txt` so we have the exact coordinates of the different domains in the protein

Note: Using the domain information we got from Interpro and reading this paper on bbox proteins (Crocco, C.D. and Botto, J.F., 2013. BBX proteins in green plants: insights into their evolution, structure, feature and functional diversification. Gene, 531(1), pp.44-52.), the domains and their positions in the tetraploid and diploid protein were adjusted. A final list of all the domains, motifs, and relevant mutations in our protein are highlighted in the table below.

| Domain | Diploid Position | Tetraploid Position |
|--------|------------------|---------------------|
| Bbox domains | 5-105 | 5-108 |
| Non Classical Nuclear Localisation Signal | 309-312 | 312-315 |
| C-Terminus VP pair | 300-301 | 303-304 |
| M6 Motif | 173-186 | 176-189 |
| M7 Motif | 100-111 | 103-114 |
| Mutation 1 | 150 | 153 |
| Mutation 2 | 204 | 207 |

### g10577

To verify the functional domains within the g10577 proteins we followed these steps:

 1. Follow the link to the [InterPro website](https://www.ebi.ac.uk/interpro/)

 2. Input the diploid protein sequences (in fasta format) for g10577 into the query field labelled `Enter your sequence` and choose `search`

 3. Manually insert the positions of the domains into a txt file `g10577_diploid_domains.txt` so we have the exact coordinates of the different domains in the protein

 4. Follow the link in step 1 to return to the Interpro home page. Input the tetraploid protein sequence (in fasta format) for g10577 into the query field labelled `Enter your sequence` and choose `search`

 5. Manually insert the positions of the domains into a txt file `g10577_tetraploid_domains.txt` so we have the exact coordinates of the different domains in the protein
    
Note that domains were introduced and domain positions were adjusted based on information in the research papers mentioned below. 

(Papolu, P.K., Ramakrishnan, M., Mullasseri, S., Kalendar, R., Wei, Q., Zou, L.H., Ahmad, Z., Vinod, K.K., Yang, P. and Zhou, M., 2022. Retrotransposons: How the continuous evolutionary front shapes plant genomes for response to heat stress. Frontiers in plant science, 13, p.1064847.). 

(Peterson-Burch, B.D. and Voytas, D.F., 2002. Genes of the Pseudoviridae (Ty1/copia retrotransposons). Molecular biology and evolution, 19(11), pp.1832-1845.).

(Systematic survey of plant LTR-retrotransposons elucidates phylogenetic relationships of their polyprotein domains and provides a reference for element classification). 

(Heslop-Harrison, J.S., Schwarzacher, T. and Liu, Q., 2023. Polyploidy: its consequences and enabling role in plant diversification and evolution. Annals of Botany, 131(1), pp.1-10.)

The final list of domain positions was not provided so only the final domain names are given in a table below. 

| Domain |
|--------|
| GAG domain|
| Integrase | 
| Reverse Transcriptase |
| RNaseH |
| Protease | 


## Part5 - Multiple Sequence Allignments and Phylogenetic Tree Building 

Code in this section should be ran on the cloud HPC.

This is the last step where we remain at the primary sequence level. We have sucessfully determined the important domains in our proteins, as well as their homlogs across different species. The next step is to identify conserved residues in our proteins and determine any important mutations between our diploid and tetraploid proteins by comparing the sequences of our proteins with their close homologs. We will also build phylogenetic trees as another form of visualising and representing how much the proteins in Cochalearia have diverged from their homologs, and also how much the tetraploid protein might have diverged from the diploid protein.

### g46214
 
By following the steps below, a multiple sequence allignment of the protein sequences was generated.

 1. Follow the link to the [uniprot](https://www.uniprot.org/align) website and navigate to the tab labelled `Align`

 2. Copy and paste all the sequences from the `g46214_homologs.fasta` you generated in Part3 of the anlysis into the query field then click the `Run align` button

 3. On the results page, click on the `download` option and choose `FASTA` format from the drop-down menu before downloading the allignment. When prompted to enter a file name, enter `g46214_allignment.fasta`

 4. To build a phylogenetic tree from this allignment we will use RAxML-NG by running the command below. This will produce a file with the extension `.besttree` which is the unrooted best maximum likelihood tree. To find out more about the `.besttree` extension, you can read through this paper (Kozlov, A.M., Darriba, D., Flouri, T., Morel, B. and Stamatakis, A., 2019. RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference. Bioinformatics, 35(21), pp.4453-4455.) or have a look on the [RAxML-NG github page](https://github.com/amkozlov/raxml-ng/wiki/Output:-files-and-settings#output-prefix)


```bash
sbatch ~/phylogenetic_tree_g46214.sh
```

We will view the `.besttree` file using the ITOL website. To do this the steps below were followed:

  1. Follow the link provided to the [ITOL website](https://itol.embl.de/)

  2. Click on the `Upload a tree` option

  3. Navigate to the `Tree file` field and select the `choose file` option. Search through your files to locate the raxml file with the `.besttree` extension and upload it.

  4. Replace the accession code names at the tip of the tree by clicking on the accession code (i.e XP_009106205.1), navigating to `label`, then select `edit label` and replace it with the correct genus and species. You can find the correct genus and species for the accession code within a table formed in Part3 of this analysis 

  5. Take a picture of the screen using CMD+SHIFT+4 and manually adjust the size to have the whole phylogenetic tree inside it.

Next, to be able to view and compare the primary sequence of the different homologs in the multiple sequence allignment, we used Jalview. The steps below should be followed:
  
  1. Open Jalview on your local machine. Select `file` in the toolbar, then select `Input Allignment` and finally choose `From File`. Navigate to your `g46214_allignment.fasta` and load it in.

  2. Using the positions of the mutations mentioned in the `Domain Identification` section of the analysis, manually search for the mutations in the allignment to determine if the residue is conserved across all homologs. To highlight residues in the allignment you should drag your cursor across that position in the allignment (i.e the residue at position 50 in all the different homlogs). Right click on the highlighted box. Choose `selection`, then `Edit New Group`, then `Group Colour`, and choose the colour `Hydrophobicity`.

  3. To make the allignment easier to understand when we create our images, we want to have the actual genus and species name for each homolog instead of the NCBI accession code. In Jalview, you have to right-click on the NCBI accession code and then select `Edit Name/Description` and use the table in the `Homolog Identification` section of the analysis to insert the correct genus and species name.

  4. Take a picture of the screen using CMD+SHIFT+4 and manually adjust the size to have the whole allignment in it.


### g10577

Note: This is written acknowledgement that these steps detailing the multiple sequence allignment and tree building for g10577 were  written and explained by a group member, Luke. File nomenclature in Steps 1 and 6 were edited for increased clarity.

The steps below will allow you to create a multiple sequence allignment for the g10577 proteins and their closest homologs, as well as neighbour joining tree showing genetic relationships .

1. Upload the g10577 homolog fasta file into MEGA11: Align -> Edit/Build Alignment -> Retrieve sequences from a file -> `g10577_homologs.fasta`

2. Produce a multiple sequence alignment by navigating to the tool bar and selecting Alignment: Alignment -> Align by ClustalW -> select all -> used default parameters except for Delay Divergent Cutoff (%) and selected 45%.

3. Next, selected Data -> Phylogenetic Analysis from the tool bar then produced a phylogenetic tree by selecting Phylogeny from the toolbar -> Construct/Test Neighbour-Joining Tree with the following parameters.

```  
1. Test of Phylogeny: Bootstrap method

2. No. Bootstrap Replications: 500

3. Substitutions type: Amino acid

4. Model/Method: Poisson model

5. Rates among Sites: Gamma Distributed (G)

6. Gamma Parameter: 1.00

7. Pattern among Lineages: Same (Homogeneous)

8. Gaps/Missing Data Treatment: Pairwise Deletion

9. Number of Threads: 3
```

5. This produced an NJ tree with default layout which was customised with the following

```
1.Taxon names -> Font -> Arial -> Bold Italic -> 10
    
2.Layout -> Toggle Scaling of the Tree + Auto-size Tree
    
3.Manually customise the length and the width of the tree to make it more aesthetic
```

6. To save a PNG file of the tree, select Image from the toolbar: Image -> Save a PNG file, then give the file path to your Desktop with an appropriate file name such as `g10577_phylogenetic_tree.png`.

## Part6 - Protein Structure Modelling

Note:All code in this section should be ran from the command line on your local machine. It is also important to note that alphafold collab has a limit to how much modelling you can do within a given time frame, so you may need multiple google accounts or a colleagues machine to complete modelling for all the domains and proteins in this step.  

We now want to visualise the three dimensional structure of our proteins. We will use protein visualisation tools and the files we are using are in the PDB format. PDB files contain information about the atoms in a protein and their coordinates. You can read more about this file type on the [rcsb website](https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/introduction#:~:text=A%20typical%20PDB%20formatted%20file,the%20atoms%20and%20their%20coordinates.)

### g46214

We will first create a directory in our home directory to host all the protein stuctures and any modelling related output by following the command below:

```bash
mkdir -p ~/g46214_modelling_output/tetraploid_g46214_protein_images/movie
mkdir -p ~/g46214_modelling_output/diploid_g46214_protein_images/movie
```

To obtain 3D structure models for our diploid and tetraploid proteins we followed the steps below:
 
 1. Follow the link to the [alphafold collab](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb) website 
 
 2. Input the diploid protein sequence into the `query sequence` field. The job name entered should be `diploid_g46214`.

 3. Navigate to the options at the top of the page, select `Runtime` and choose `Run all`

 4. When the modelling has been completed, on the Safari web browser (version 15.6) you will be prompted to allow the resulting file to be downloaded and selecting `allow` will download a zip file into your downloads folder (Mac). Note if you have selected a different directory as your default directory for downloads to be sent to, you will have to change it back to the `Downloads` folder for the purpose of following this analysis.

 5. Follow the link to the [alphafold collab](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb) website or refresh the page.
 
 6. Input the tetraploid protein sequence into the `query sequence` field, and the job name entered should be `tetraploid_g46214`.

 7. Repeat steps 3 and 4

 8. Move all the downloaded `.zip` files (including those from modelling the g46214 diploid protein) from your `Downloads` folder to the directory we created earlier for protein structures, and open the files following the commands below:


```bash
mv ~/Downloads/*g46214*.zip ~/g46214_modelling_output
for file in ~/g46214_modelling_output/*.zip; do unzip "$file"; done
```

By unzipping the files,there should be a directory for the g46214 diploid protein modelling output and the g46214 tetraploid protein modelling output. The name of the directory should start with `tetraploid_g46214` or `diploid_g46214` , but alphafold assigns a random string of characters to the end of the directory name (i.e erk453). Each directory will contain the top 5 models that alphafold generated in pdb format for that modelling job. 

### g10577

We will first create a directory to host all the protein stuctures and any modelling related output by following the command below:

```bash
mkdir ~/g10577_modelling_output
```

For g10577 we had to use a different approach to model the proteins due to limitations with alphafold's memory and being unable to model the whole 1000+ amino acid long protein. Instead of putting the whole tetraploid or diploid sequence into alphafold, the protein was fragmented into the longest possible segments that would avoid alphafold from running out of time whilst modelling. We decided on this method as opposed to using an alternative modelling software like Phyre2 because we did not get biologically sensible output when those were used to model our proteins. 

This approach to modelling also requires you to obtain a complete reference protein model using [SWISS-MODEL](https://swissmodel.expasy.org/interactive) so you can essentially 'map' the protein fragments onto the reference protein to re-build our diploid and tetraploid proteins. The protein we chose for our reference had the uniprot ID `Q9LPK1` and we settled on that as a reference because of its high coverage and sequence identity (74.62%) to our reference g10577 protein sequence. 

Note: For the approach used in this stage of the analysis to be reproduceable you will need the exact coordinates of the fragments we are taking from the diploid or tetraploid protein(i.e if you should take residue 1-400, followed by 401-600 etc), however these were not provided and so although the steps outlined below are what was done, they are missing that bit of information.

The analysis can be performed by following the steps below:
    
 1. We will first need to get the 3D structure for the reference protein. Follow the link to the [SWISS-MODEL](https://swissmodel.expasy.org/interactive) website, select `Start Modelling`, input the g10577 reference protein sequence into the `Target Sequence` field and select `Build Model`. Identify the model with template `Q9LPK1.1.A`, and download the model in PDB format to your downloads folder. The name of the reference model downloaded from SWISS-MODEL is `model_01.pdb` and that can be changed by navigating to the directory with the `model_01.pdb` and following the command below.
    
```bash
mv model_01.pdb reference_g10577.pdb
```
 2. Follow the link to the [alphafold collab](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb) website

 3. Input the first fragment for the diploid sequences into the `query sequence` field, and the job name entered should be `diploid_framgment1_g10577`. 

 4. Navigate to the options at the top of the page, select `Runtime` and choose `Run all`

 5. When the modelling has been completed, on the Safari web browser (version 15.6) you will be prompted to allow the resulting file to be downloaded and selecting `allow` will download a zip file into your downloads folder (Mac). Note if you have selected a different directory as your default directory for downloads to be sent to, you will have to change it back to the `Downloads` folder for the purpose of following this anaysis.

 6. Repeat steps 2-5 for the other fragments of the diploid g10577 protein. At step 3 you should change the job name to reflect which fragment is being modelled (i.e using `diploid_fragment2_g10577` when modelling the second fragment)

 7. Follow the link to the [alphafold collab](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb) website or refresh the page

 8. Input the first fragment for the tetraploid g10577 protein sequence into the `query sequence` field, and the job name entered should be `tetraploid_framgment1_g10577`.

 9. Repeat steps 2-5 for the other fragments of the tetraploid g10577 protein. At step 3 you should change the job name to reflect which fragment is being modelled (i.e using `tetraploid_fragment2_g10577` when modelling the second fragment)
 
 10. Manualy move the downloaded `.zip` files, and reference pdb model from your `Downloads` folder to the directory we created earlier for g10577 protein structures.

By unzipping the files,there should be a directory for each of the diploid g10577 fragment protein modelling output and the tetraploid g10577 fragment protein modelling output. The name of the directory should start with `tetraploid_fragment` or `diploid_fragment` , but alphafold assigns a random string of characters to the end of the directory name (i.e mhy429). Each directory will contain the top 5 models that alphafold generated in pdb format for that modelling job. 

## Part7 - Image Generation

All the code in this section of the analysis is entered into the PyMOL command line.

We have successfully modelled our proteins and now want to actually investigate the mutations in three dimensional space and create images to be used in our papers/presentations. We will load the different pdb files into PyMOL, highlight domains or motifs of interest, and take snapshots of our proteins. 

As mentioned before, the alphafold output directories are named with unique identifiers such that the directory name for the modelling job is always different every time you run it. As such the paths given in our commands are generic to not cause any confusion. The output directory for each modelling job should have a number of different models ranked from 001 to 005. We will always be using the rank 001 model for the purposes of investigating our protein structure and mutations, and we take it as our best estimate of the true protein structure.

### g46214 Tetraploids

To get images for the g46214 tetraploid protein, we first load the protein into PyMOL and colour it by hydrophibicity. To do this you need to enter the commands below into the PyMOL command line:

```bash
load ~/g46214_modelling_output/your_alphafold_tetraploid_output_directory/your_rank_001_tetraploid.pdb, tetraploid_g46214

select tetraploid_g46214_bbox_domains, tetraploid_g46214 and resi 5-108

run ~/g46214_modelling_output/colorh.py

color_h tetraploid_g46214
```
Now that the tetraploid g46214 protein has been loaded into PyMOL, using your mouse or trackpad, manually adjust the orientation of the tetraploid protein to a position you are happy with. Follow the next few steps:

 1. Navigate to the header of PyMOL and select the plugin tab, select APBS electrostatics.

 2. Select the drop down menu in the selection entry field (selection:[       ]) and select `polymer & tetraploid_g46214` depending on which protein you are investigating. This will produce an object in PyMOL showing the electrostatic potential across the whole protein. When that has completed close the pop-up that comes afterwards.

 3. Repeat steps 1 and 2 but select `polymer & tetraploid_g46214_bbox_domains` in the selection entry field. This will produce an object showing the electrostatic potential across the bbox domains only. Once this is completed you might have to manually zoom out so your whole protein is showing on the screen and the bbox domain is not being focused on. You should not select or deselect any of the objects named on the right hand side of the PyMOL window as this causes the scripts that generate images to not function correctly.

Now we have everything we need to produce our images and you can run the `tetraploid_domain_highlight.py` script to get your figures that show the tetraploid g46214 protein with domains highlighted, coloured by hydrophobicity, and coloured by electrostatic potential at 90 degree angles. The images will be located in the `~/g46214_modelling_output/tetraploid_g46214_protein_images` directory and are clearly named based on if that frame was coloured by hydrophobicity for example. In the PyMOL commmand line you can should enter the command below:

```bash
run ~/path/to/python/script/tetraploid_domain_highlight.py
```

NOTE: Regarding the images showing the protein coloured by electrostatic potential, if you are using too much RAM on your machine the whole protein may be coloured white or only a few areas of your protein may be coloured by electrostatic potential. If that occurs close all other running applications on your device and rereun the steps, or restart your device and rereun the steps. 

### g46214 Diploids

To get good images for the g46214 diploid protein, we first load the protein into PyMOL and colour it by hydrophibicity. To do this you need to enter the commands below into the PyMOL command line:

```bash
load ~/g46214_modelling_output/your_alphafold_diploid_output_directory/your_rank_001_diploid.pdb, diploid_g46214

select diploid_g46214_bbox_domains, diploid_g46214 and resi 5-105

run ~/g46214_modelling_output/colorh.py

color_h diploid_g46214
```

Now that the diploid g46214 protein has been loaded into PyMOL, Using your mouse or trackpad, manually adjust the orientation of the diploid protein to a position you are happy with. Follow the next few steps:

 1. Navigate to the header of PyMOL and select the plugin tab, select APBS electrostatics.

 2. Select the drop down menu in the selection entry field (selection:[       ]) and select `polymer & diploid_g46214` depending on which protein you are investigating. This will produce an object in PyMOL showing the electrostatic potential across the whole protein. When that has completed close the pop-up that comes afterwards.

 3. Repeat steps 1 and 2 but select `polymer & diploid_g46214_bbox_domains`. This will produce an object showing the electrostatic potential across the bbox domains only. Once this is completed you might have to zoom out so your whole protein is showing on the screen and the bbox domain is not being focused on. You should not select or deselect any of the objects named on the right hand side of the PyMOL window as this causes the scripts that generate images to not function correctly.

Now we have everything we need to produce our images and you can run the `diploid_domain_highlight.py` script to get your nice figures that show the diploid g46214 protein with domains highlighted, coloured by hydrophobicity, and coloured by electrostatic potential at 90 degree angles. The images will be located in the `~/g46214_modelling_output/diploid_g46214_protein_images` directory. In the PyMOL commmand line you should enter the command below:

```bash
run ~/path/to/python/script/diploid_domain_highlight.py
```

### g10577 Tetraploids

Due to the limitations with alphafold's memory we have modelled fragments of the tetraploid g10577 protein and not the whole protein.

We will first have to load in the reference protein which was retrieved from SWISS-MODEL (uniprot ID Q9LPK1), then we will load in all the different fragments of the protein and one-by-one these will be aligned to the reference protein to effectively stitch together the original protein. 

This can be achieved by following the steps outlined below:

 1. Open a new PyMOL window and load the reference SWISSMODEL protein into PyMOL using the commands below

 ```bash
 load /path/to/your/reference/protein.pdb, reference_g10577
 ```
 2. For each fragment you have modelled in the `Protein Structure Modelling` section of this analysis, we want to load that bit of protein sequence into PyMOL and give the object a name within PyMOL so that it is easily identifiable. The [load command](https://pymolwiki.org/index.php/Load) takes the pdb file you want to load, then the name you want to assign to that object once it is loaded into PyMOL. You will again be choosing the rank 001 model from the alphafold output because it gives us the best estimate at the actual protein structure. Because the exact fragments were not provided, only the code for fragment 1 is outlined below, but ideally you would repeat this step for each fragment that was modelled and change the name of the object to reflect the fragment that was loaded in.
   
```bash
load /path/to/your/tetraploid/fragment1/directory/your_rank_001_tetraploid_fragment1.pdb, tetraploid_framgment1
```
   
 3. Now we will align each fragment to the reference SWISS-MODEL protein to try and map the domains and recreate a complete protein. The [align](https://pymolwiki.org/index.php/Align) command in PyMOL first takes the name of the object you want to align, followed by what object you want to align it to. Like with the command above this one, because the exact fragments were not provided only the code for fragment 1 is outlined below, but ideally you would repeat this step for each fragment that was modelled and change the name of the object to reflect the fragment that was being aligned.

```bash
align tetraploid_fragment1, reference_g10577
```

 4. We want to create an object for the different domains within our newly assembled protein so they can be highlighted different colours in our images. The protein fragments that we have modelled contain the different protein domains within them so for each fragment we have to take into account the position of the domain within the full protein. For example if fragment 1 spans from residue 1-400, the integrase domain might only be from residue 200-330 of that fragment of the protein. The select command in PyMOL allows us to specify a specific set of residues to highlight. The [select](https://pymolwiki.org/index.php/Select) command in PyMOL first takes the name for the new object it will be creating (the name of the actual domain), then it takes the object you are selecting a given set of residues from (the fragment with the domain within it) and the residues you are selecting (the positions of the domain within the fragment). Following the code below would allow you to to select a given domain within a given fragment. x and y represent the start and stop position for the domain.

```bash
 select integrase, tetraploid_fragment1 and resi x-y`'
```
5. Colour the different domains by navigating to the object list at the right hand side of the window. Look for the name of the domain of interest, select the `C` and then select the colour based on the colour scheme outlined blow.

   GAG_domain -> sky blue
   
   integrase -> orange

   protease -> yellow

   reverse_transcriptase -> green

   RNaseH -> purple

6. We want to only highlight the functional domains we have just coloured so we will first hide all protein models that have been loaded in. On the right hand side of the PyMOL window look for an object called `all`. Select the `H` and then select `everything`. This should hide all objects you have loaded into PyMOL so far. Each of the domains mentioned in step 5 should have an object on the right hand side of the PyMOL window. For each domain, select the `S` next to the object's name and then select `cartoon`. The PyMOL window will now only contain the domains of the protein.

7. Manually rotate the protein based on personal preference and take a picture of the domains using the command below in the PyMOL command line.
```bash
png ~/path/to/ouput_image/directory/tetraploid_image.png, 3500, 3500, -1, ray=0, dpi=500
```
8. To get the electrostatic potential for all the domains navigate to the header of PyMOL and select the plugin tab, and then select APBS electrostatics. Select the drop down menu in the selection entry field (selection:[       ]) and select `polymer & <domain name>`. This will produce an object in PyMOL showing the electrostatic potential across the given domain. When that has completed close the pop-up that comes afterwards.

9. Repeat step 8 for all the domains we have outlined in step 5 to colour all the domains by electrostatic potential.

10. Manually rotate the protein (now coloured by electrostatic potential) based on personal preference and take a picture using the command below in the PyMOL command line
```bash
png ~/path/to/ouput_image/directory/tetraploid_electrostatic_image.png, 3500, 3500, -1, ray=0, dpi=500
```

   
### g10577 Diploids
Due to the limitations with alphafold's memory we have modelled fragments of the diploid g10577 protein and not the whole protein.

We will first have to load in the reference protein which was retrieved from SWISS-MODEL (uniprot ID Q9LPK1), then we will load in all the different fragments of the protein and one-by-one these will be aligned to the reference protein to effectively stitch together the original protein. 

This can be achieved by following the steps outlined below:

 1. Open a new PyMOL window and load the reference SWISSMODEL protein into PyMOL using the commands below

 ```bash
 load /path/to/your/reference/protein.pdb, reference_g10577
 ```
 2. For each fragment you have modelled in the `Protein Structure Modelling` section of this analysis, we want to load that bit of protein sequence into PyMOL and give the object a name within PyMOL so that it is easily identifiable. The load command takes the pdb file you want to load, then the name you want to assign to that object once its loaded into PyMOL. You will again be choosing the rank 001 model from the alphafold output because it gives us the best estimate at the actual protein structure modelled by alphafold. Because the exact fragments were not provided, only the code for fragment 1 is outlined below, but ideally you would repeat this step for each fragment that was modelled and change the name of the object to reflect the fragment that was loaded in.
   
```bash
load /path/to/your/diploid/fragment1/directory/your_rank_001_diploid_fragment1.pdb, diploid_framgment1
```
   
 3. Now we will align each fragment to the reference SWISS-MODEL protein to try and map the domains and recreate a complete protein. The [align](https://pymolwiki.org/index.php/Align) command in PyMOL first takes the name of the object you want to align, followed by what object you want to align it to. Like with the command above this one, because the exact fragments were not provided only the code for fragment 1 is outlined below, but ideally you would repeat this step for each fragment that was modelled and change the name of the object to reflect the fragment that was loaded in.

```bash
align diploid_fragment1, reference_g10577
```

 4. We want to create an object for the different domains within our newly assembled protein so they can be highlighted different colours in our images. The protein fragments that we have modelled contain the different protein domains within them so for each fragment we have to take into account the position of the domain within the full protein. For example if fragment 1 spans from residue 1-400 and the integrase domain would be from residue 200-330 of that fragment of the protein and the select command in PyMOL allows us to do that. The [select](https://pymolwiki.org/index.php/Select) command in PyMOL first takes the name for the new object it will be creating (the name of the actual domain), then it takes the object you are selecting a given set of residues from (the fragment with the domain) and the residues you are selecting (the positions of the domain within the fragment). Following the code below would allow you to to select a given domain within a given fragment. x and y represent the start and stop position for the domain.

```bash
 select integrase, diploid_fragment1 and resi x-y`'
```
5. Colour the different domains by navigating to the object list at the right hand side of the window. Look for the name of the domain of interest, select the `C` and then select the colour based on the colour scheme outlined blow.

   GAG_domain -> sky blue
   
   integrase -> orange

   protease -> yellow

   reverse_transcriptase -> green

   RNaseH -> purple

6. We want to only highlight the functional domains we have just coloured so we will first hide all protein models that have been loaded in. On the right hand side of the PyMOL window look for an object called `all`. Select the `H` and then select `everything`. This should hide all objects you have loaded into PyMOL so far. Each of the domains mentioned in step 5 should have an object on the right hand side of the PyMOL window. For each domain, select the `S` next to the object's name and then select `cartoon`. The PyMOL window will now only contain the domains of the protein.

7. Manually rotate the protein based on personal preference and take a picture of the domains using the command below in the PyMOL command line.
```bash
png ~/path/to/ouput_image/directory/diploid_image.png, 3500, 3500, -1, ray=0, dpi=500
```
8. To get the electrostatic potential for all the domains navigate to the header of PyMOL and select the plugin tab, and then select APBS electrostatics. Select the drop down menu in the selection entry field (selection:[       ]) and select `polymer & <domain name>`. This will produce an object in PyMOL showing the electrostatic potential across the given domain. When that has completed close the pop-up that comes afterwards.

9. Repeat step 8 for all the domains we have outlined in step 5 to colour all the domains by electrostatic potential.

10. Manually rotate the protein (now coloured by electrostatic potential) based on personal preference and take a picture using the command below in the PyMOL command line
```bash
png ~/path/to/ouput_image/directory/diploid_electrostatic_image.png, 3500, 3500, -1, ray=0, dpi=500
```
## Step7 - Movie Generation
We also want to have a short video showing some of our protein molecules rotating over time. The logic behind this is the same as with a flipbook. Over 360 degrees, we will rotate the protein molecule by 1 degree across a given axis, and take a picture after each rotation. Combining those individual images together will make a movie. This can be accomplished by running the scripts below: 

### g46214 Tetraploids

You will have to open a new PyMOL window and then run the commands below in the PyMOL command line
```bash
load ~/g46214_modelling_output/gatk_04AF_tetraploid_b50ff/gatk_04AF_tetraploid_b50ff_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb, tetraploid_g46214

select tetraploid_g46214_bbox_domains,tetraploid_g46214 and resi 5-108

run ~/g46214_modelling_output/colorh.py

color_h tetraploid_g46214
```

Now that the tetraploid g46214 protein has been loaded into pymol, Using your mouse or trackpad, manually adjust the orientation of the tetraploid protein to a position you are happy with. Follow the next few steps:

 1. Navigate to the header of pymol and select the plugin tab, select APBS electrostatics.

 2. Select the drop down menu in the selection entry field (selection:[       ]) and select `polymer & tetraploid_g46214` . This will produce an object in PyMOL showing the electrostatic potential across the whole protein. When that has completed close the pop-up that comes afterwards.

Now rotate molecule to choose a good starting position that will show everything you want when it rotates. It will be rotating up to down and this can be accomplished running the script below to produce a series of images after every rotation

In the pymol command line run the code below to produce a directory with a collection of tetraploid g46214 images.
```bash
run ~/tetraploid_temporary_image_generation.py 
```

### g46214 Diploids

You will have to open a new PyMOL window and then run the commands below in the PyMOL command line

```bash
load ~/g46214_modelling_output/gatk_04AF_dip_29144/gatk_04AF_dip_29144_unrelaxed_rank_001_alphafold2_ptm_model_2_seed_000.pdb, diploid_g46214

select diploid_g46214_bbox_domains,diploid_g46214 and resi 5-105

run ~/g46214_modelling_output/colorh.py

color_h diploid_g46214
```

Now that the diploid g46214 protein has been loaded into pymol, Using your mouse or trackpad, manually adjust the orientation of the diploid protein to a position you are happy with. Follow the next few steps:

 1. Navigate to the header of pymol and select the plugin tab, select APBS electrostatics.

 2. Select the drop down menu in the selection entry field (selection:[       ]) and select `polymer & diploid_g46214` . This will produce an object in PyMOL showing the electrostatic potential across the whole protein. When that has completed close the pop-up that comes afterwards.

Now rotate molecule to choose a good starting position that will show everything you want when it rotates. It will be rotating up to down and this can be accomplished running the script below to produce a series of images after every rotation

In the pymol command line run the code below to produce a directory with a collection of diploid g46214 images
```bash
run ~/diploid_temporary_image_generation.py 
```

The series of images produced for the diploid protein only and the tetraploid protein only will be stitched together to produce a movie file for the diploid protein rotating and a movie file for the tetraploid protein rotating. This can be accomplished by running the script below in the command line

```bash
bash make_protein_rotation_movie.sh
```

# Conclusion
By following the steps outlined in this analysis, we were able to identify g46214 as showing strong similarity to bbox21 proteins from the bbox family, and g10577 showed strong similarity to retrotransposons.

In our analysis, we found a truncation in the C-terminus of the tetraploid g10577 which might impair its' function. Given the genomic instability polyploids are faced with it might be better to not have this transposon function

With regards to g46214, we found mutations that increased the hydrophobicity of a motif that might be necessary for protein-protein interactions. These results still need to be experimentally validated.

