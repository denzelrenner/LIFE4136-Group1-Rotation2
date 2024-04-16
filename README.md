# LIFE4136_Group1_Rotation2
This repository contains all the scripts and lists of data files needed to reproduce results from group1 in rotation2 of the groupwork projects.

Note: Unless on the cloud or stated otherwise, all command line code was ran on a Mac Intel i5 or Lukes Mac

## What is the problem we have been presented with?
We have been provided with some genes in Cochaleria that are under selection during the stabalisation of polyploidy. Polyploids are individuals that arise from a whole genome duplication (WGD) event, a process that results in multiple sets of chromosomes within an individual. These polyploids evolve mechanisms to handle the hardships that come with WGD, such as meiotic instability and regulation of gene expression, as well as develop adaptations to tolerate the harsh environmental conditions that can sometimes lead to polyploidy. In this project we have been provided with two genes caught in the selection scan, g46214 and g10577, and we are to investigate each gene in diploids and tetraploids through studying the amino acid sequences and tetriary structures of their proteins to determine the consequence of any genetic variants.

## What are the aims of our analysis?
We were aiming to translate the genes and investigate the primary sequences of the proteins in diploids and tetraploids to first identify important domains and predict the function of the protein. We also aimed to identify any mutations that might result in altered function of the protein (or a specific domain within the protein) in diploids or tetraploids by comparing the primary sequences with that of close homologs. Another objective of the analysis was to generate structural predictions for the diploid and tetraploid alleles for each of our candidate genes which would allow us to compare the tetriary structures of the proteins at the different ploidy levels, and observe the impact of the different mutations on the tertiary structure and other protein properties like hydrophobicity or electric charge. 

## What are the expected outcomes?
Depending on the gene and given our understanding of the system, we expect to see mutations and amino acid changes in the tetraploid g46214 and g10577 proteins relative to the diploid proteins that help stabilise the polyploids by helping to manage some of the stresses that come with polyploidy, such as controlling gene expression, as well as mutations that confer enhanced ecological tolerance. The altered function of the proteins could result from the mutation of a conserved, functional residue which enhances or reduces the interaction with other ions and proteins, whilst in some cases a mutation can be deleterious and remove an entire functional domain.

# Prerequisites
## Required files and data
All these files should be downloaded into the same directory before following the rest of this document.
The vcf file contains variants. 
The gff contains annotations for the reference genome, and the fasta is of the reference and has been indexed to produce a .fai compatible with many tools like samtools and GATK

|Directory|Files|
|---------|-----|
|VCFs|UK_scan_dips.vcf<br>UK_scan_tets.vcf|
|Reference_Genome|C_excelsa_V5_braker2_wRseq.gff3<br>C_excelsa_V5.fasta<br>C_excelsa_V5.fasta.fai|

The `.fai` file is a fasta index file which allows us to find a particular nucelotide at specific genomic coordinates in the FASTA file. You can read more about it on the [GATK website](https://gatk.broadinstitute.org/hc/en-us/articles/360035531652-FASTA-Reference-genome-format)

## Tool versions and links
These are all the tools that were used in our analysis with versions and links provided where applicable. Dependencies for certain packages, and their versions, are placed in parentheses. Some references were chosen based on what was recommended on the tool's online help page.

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
|[Conda](https://github.com/conda/conda)|version 23.5.2|NA|
|[Homebrew](https://brew.sh/)|version 4.2.10|NA|
|ffmpeg|version 6.1.1|NA|

## Tool intallation 
The guidance below outlines the necessary scripts,steps or commands that need to be ran to install some of the tools necessary to reproduce our results.

### On the cloud HPC:
To install raxml-ng run the script below in the command line.
```bash
bash ~/raxml_install.sh
```

To install `GATK`, `samtools`, and `conda` follow the guidance on their respective github pages as those tools were installed by a bioinformatics technician and code to install them was not generated from scratch in our analysis. 


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
| `phylogenetic_tree_g46214.sh` | runs raxml-ng to build a maximum likelihood tree using the allignment for g46214 in tetraploid and diploid and the homologs |
| `generate_reference_sequences.sh` | produces a nucleotide sequence (in fasta format) for the g10577 and g46214 reference genes |
| `g46214_gatk_consensus_final.sh` | produces a consensus nucleotide sequence (in fasta format) for the g46214 diploid and tetraploid genes |
| `g10577_gatk_consensus_final.sh` | produces a consensus nucleotide sequence (in fasta format) for the g10577 diploid and tetraploid genes |
| `colorh.py` | colours proteins by hydrophobicity. You need to create a file on your local machine called colorh.py and copy and paste the code from the [PyMOL wiki color h page](https://pymolwiki.org/index.php/Color_h) into the file|
| `tetraploid_domain_highlight.py`| produces a series of images (rotated by 90 degrees) of the tetraploid g46214 protein with its domains coloured and highlighted, the electrostatic potential of the whole protein, and the hydrophobicity of the whole protein |
| `diploid_domain_highlight.py`| produces a series of images (rotated by 90 degrees) of the diploid g46214 protein with its domains coloured and highlighted, the electrostatic potential of the whole protein, and the hydrophobicity of the whole protein |
| `tetraploid_temporary_image_generation.py` | prdocues a series of images rotated by 1 degree across 360 degrees for the tetraploid g46214 protein |
| `diploid_temporary_image_generation.py` | prdocues a series of images rotated by 1 degree across 360 degrees for the diploid g46214 protein |
| `make_protein_rotation_movie.sh` | stitches together the 360 diploid and tetraploid g46214 protein images to make a short movie of the protein rotating and showing all domains/faces of the protein |


# THE ANALYSIS

## Consensus Sequences (Nucleotide)
The first thing we will do is get consensus sequences for our diploids and tetraploids. We will filter the vcfs we have to only include biallelic variants with an allele frequency greater than 0.49 (to include allele frequencies of 0.5) in the case of g46214, and greater than 0.6 in the case of g10577. The resulting filtered vcfs will then be used, along with the reference fasta and gff, to produce our consensus sequences for the genes of interest. This is accomplished by following the steps below:

```bash
conda activate /shared/apps/conda/bio2

bash ~/generate_reference_sequences.sh
bash ~/g46214_gatk_consensus_final.sh
bash ~/g10577_gatk_consensus_final.sh
```
Running these scripts will produce `C_excelsa_V5.dict` and `C_excelsa_V5.fasta.fai` files which are required by GATK tools to access specified regions of the reference fasta. The `.dict` file describes the contents of our fasta file, and as mentioned before the `.fai` file is a fasta index file which allows us to find a particular nucelotide at specific genomic coordinates in the FASTA file. You can read more about these file formats on the [GATK website](https://gatk.broadinstitute.org/hc/en-us/articles/360035531652-FASTA-Reference-genome-format)

## Protein Sequences 
Amongst the files produced from running the gatk consensus scripts there should be fasta files containing the entire consensus coding sequence for g46214 and g10577 in our diploids, tetraploids, and reference. Now that we have the nucleotide sequences, we can translate them to get our protein sequences. This can be accomplished by following the steps below.

 1. Follow the link to the [orf finder](https://www.ncbi.nlm.nih.gov/orffinder/) website.

 2. Input the nucleotide seqeunce into the query seqeunce field and submit the job. On the output page, copy and paste the amino acid seqeuence for ORF1 (the longest open reading frame) into a new fasta file with any identifiable headers you prefer (i.e `>diploid_g46214` or `>tetraploid_g10577`).


## Homolog identification

### g46214
Now that we have retrieved the sequences for our proteins, we want to figure out what they actually are, and what their closest homologs are in other plant or animal species. You should create a fasta file called `g46214_homologs.fasta` to store all homologs. Homologous proteins were selected based on having 100% query cover and >60% percentage identity to the reference protein sequence. 

To identify g46214 homologs:

 1. Follow the link to the [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/) webpage

 2. Select `Protein BLAST` 

 3. Input the protein sequence for the g46214 reference protein sequences into the `Enter Query Sequence` field, and select `BLAST` at the bottom of the page

 4. Manually retrieve the protein sequence of the homologs from the blastp results page by selecting their ncbi dataset accession code (i.e `XP_018447019.1`) in the blastp results page

 5. Select the `FASTA` option at the top of the page, and finally copy and paste the protein sequence into the `g46214_homologs.fasta` fasta file

 6. The consensus protein sequences for the reference, diploid and tetraploid protein sequences should also be added to the homolog files in fasta format.

The accession codes for the homologs used in subsequent steps of the analysis is outlined below:
| Homolog | Accession Code | 
|--------|------------------|
|  | KAJ0255869.1 |
| | XP_009106205.1 |
| | XP_018447019.1 |
|  | XP_013648757.1 |
|  | NP_177686.1 |
|  | KAG7659674.1 |
|  | KAG7651806.1 |
|  | XP_006390294.1 |
|  | XP_006302550.1 |

#### g10577

For g10577, homologous proteins were selected based on >70% query cover and >40% percentage identity to the reference protein sequence.

To identify g10577 homologs in other plant species we followed the steps below:

 1. Follow the link to the [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/) webpage

 2. Select `Protein BLAST` 

 3. Input the protein sequence for the g10577 reference proteins sequences into the `Enter Query Sequence` field, and select `BLAST` at the bottom of the page

 4. Manually retrieve the protein sequence of the homologs from the blastp results page by selecting their ncbi dataset accession code (i.e `XP_018447019.1`) in the blastp results page

 5. Select the `FASTA` option at the top of the page, and finally copy and paste the protein sequence into the `g10577_homologs.fasta` fasta file

 6. The consensus protein sequences for the reference, diploid and tetraploid protein sequences should also be added to the homolog files in fasta format.
Note that we also searched for homologs in the model species Arabidopsis thaliana on the TAIR website, but for g10577 the identified homolog did not make biological sense when investigated further through multiple sequence allignments and the literature cited in the `Domain Identification` section of this analysis so it was not we suspect it was not that.

## Domain identification
We have found out the closest homologs for our proteins in different species so we can begin to start investigating our proteins function, and structural domains.

### g46214

To verify the functional domains within the g46214 proteins we followed these steps:

 1. Follow the link to the [InterPro website](https://www.ebi.ac.uk/interpro/)

 2. Input the diploid and tetraploid protein sequences (in fasta format) for g46214 into the query field labelled `Enter your sequence` and choose `search`

 3. Manually insert the positions of the domains into a txt file called (i.e `g46214_tetraploid_domains.txt` and `g46214_diploid_domains.txt`) so we have the exact coordinates of the different domains in the protein

Note: Using the domain information from following the steps above and this paper on bbox proteins (Crocco, C.D. and Botto, J.F., 2013. BBX proteins in green plants: insights into their evolution, structure, feature and functional diversification. Gene, 531(1), pp.44-52.), the  domains and their position in the tetraploid/diploid protein found on interpro was adjusted. A final list of all the domains, motifs, and mutations in our protein are highlighted in the table below.

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

 2. Input the diploid and tetraploid protein sequences (in fasta format) for g10577 into the query field labelled `Enter your sequence` and choose `search`

 3. Manually insert the positions of the domains into a txt file (i.e `g10577_tetraploid_domains.txt` and `g10577_diploid_domains.txt`) so we have the exact coordinates of the different domains in the protein

Note that domains and domain positions were introduced and adjusted based on information in these papers.
(Papolu, P.K., Ramakrishnan, M., Mullasseri, S., Kalendar, R., Wei, Q., Zou, L.H., Ahmad, Z., Vinod, K.K., Yang, P. and Zhou, M., 2022. Retrotransposons: How the continuous evolutionary front shapes plant genomes for response to heat stress. Frontiers in plant science, 13, p.1064847.). 
(Peterson-Burch, B.D. and Voytas, D.F., 2002. Genes of the Pseudoviridae (Ty1/copia retrotransposons). Molecular biology and evolution, 19(11), pp.1832-1845.).
(Systematic survey of plant LTR-retrotransposons elucidates phylogenetic relationships of their polyprotein domains and provides a reference for element classification). 
(Heslop-Harrison, J.S., Schwarzacher, T. and Liu, Q., 2023. Polyploidy: its consequences and enabling role in plant diversification and evolution. Annals of Botany, 131(1), pp.1-10.)

| Domain |
|--------|
| GAG domain|
| Integrase | 
| Reverse Transcriptase |
| RNaseH |
| Protease | 



## Multiple Sequence Allignments and Phylogenetic Tree Building 

This is the last step where we remain at the primary sequence level. We have sucessfully determined the important domains in our proteins, as well as their homlogs across different species. The next step is to identify conserved residues in our proteins and determine any important mutations between our diploid, tetraploid,and reference proteins, by comparing the sequences of our proteins with their close homologs. We will also build phylogenetic trees as another form of visualising and representing how much the proteins in Cochalearia have diverged from their homologs, and also how much the tetraploid have diverged from the diploids.

### g46214
 
By following the steps below, a multiple sequence allignment of the protein sequences was generated.

 1. Follow the link to the uniprot website and navigate to the tab labelled `Align`

 2. Copy and paste all the sequences from the fasta file you generated above into the query field then click the `Run align` button

 3. On the results page, click on the `download` option and choose `FASTA` format from the drop-down menu before downloading the allignment. When prompted to enter a file name, enter `g46214_allignment.fasta`

 4. To build a phylogenetic tree from this allignment we will use raxml by running the command below. This will produce a file with the extension `.besttree` which is the unrooted best maximum likelihood tree.

```bash
sbatch ~/phylogenetic_tree_g46214.sh
```

We will view the `.besttree` file using the ITOL website. To do this the steps below were followed:

  1. Follow the link provided to the [ITOL website](https://itol.embl.de/)

  2. Click on the `Upload a tree` option

  3. Navigate to the `Tree file` field and select the `choose file` option. Search through your files to locate the raxml file with the `.besttree` extension and upload it.

  4. Replace the accession code names at the tip of the tree by clicking on the accession code (i.e XP_009106205.1), navigating to `label`, then select `edit label` and replace it with its Genus and Species

  5. Capture the screen using CMD+SHIFT+4 and adjusting the size to have the whole phylogenetic tree inside it.

Based on the results from the allignment and a [paper detailing different bbox proteins](Crocco, C.D. and Botto, J.F., 2013. BBX proteins in green plants: insights into their evolution, structure, feature and functional diversification. Gene, 531(1), pp.44-52.), new domains were documented for the bbox proteins, and the positions for some of the interpro domains were altered. The exact residues are detailed in scripts mentioned below and as such will not be repeated in this statement.

To find out more about the `.besttree` extension, you can read through this paper (Kozlov, A.M., Darriba, D., Flouri, T., Morel, B. and Stamatakis, A., 2019. RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference. Bioinformatics, 35(21), pp.4453-4455.) or the [RAxML-NG github page](https://github.com/amkozlov/raxml-ng/wiki/Output:-files-and-settings#output-prefix)

### g10577

Note: This step was written and explained by a group member, Luke.

1. Made a new fasta file of a multiple sequence alignment using the top BLASTp hits with the g10577_tetraploid_GATK_>0.75 amino acid sequence

2. Uploaded the fasta file into MEGA11: Align -> Edit/Build Alignment -> Retrieve sequences from a file -> `g10577_homologs.fasta`

3. Produced a multiple sequence alignment by navigating to the tool bar and selecting Alignment: Alignment -> Align by ClustalW -> select all -> used default parameters except for Delay Divergent Cutoff (%) and selected 45%.

4. Next, selected Data -> Phylogenetic Analysis from the tool bar then produced a phylogenetic tree by selecting Phylogeny from the toolbar -> Construct/Test Neighbour-Joining Tree with the following parameterss

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

5. This produced an NJ tree with default layout – I customised the tree with the following

```
1.Taxon names -> Font -> Arial -> Bold Italic -> 10
    
2.Layout -> Toggle Scaling of the Tree + Auto-size Tree
    
3.I customised the length and the width of the tree to make it more aesthetic
```

6. To save a PNG file of the tree I selected Image from the toolbar: Image -> Save a PNG file, then gave the file path to my Desktop with an appropriate file name.

Note, after performing multiple sequence allignments and reading through papers (which papers?) the positions of domains in the g10577 diploid and tetraploid were manually adjusted to better reflect the real positions of functional domains

## Protein Structure Modelling

We now want to visualise the three dimensional structure of our proteins. The proteins we will be loaded into visualisation tools and files are in the PDB format. PDB files contain information about the atoms in the protein and their coordinates. You can read more about this file type on the [rcsb website](https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/introduction#:~:text=A%20typical%20PDB%20formatted%20file,the%20atoms%20and%20their%20coordinates.)

### g46214

We will first create a directory in our home directory to host all the protein stuctures and any modelling related output by following the command below:

```bash
mkdir -p ~/g46214_modelling_output/tetraploid_g46214_protein_images/movie
mkdir -p ~/g46214_modelling_output/diploid_g46214_protein_images/movie
```

To obtain 3D structure models for our proteins we followed the steps below:
 
 1. Follow the link to the [alphafold collab](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb) website 
 
 2. Input the protein sequence for reference,diploid,and tetraploid sequences into the `query sequence` field, and for a given candidate gene , the job name should have some identifier they all share (i.e `diploid_g46214`,`tetraploid_g46214`). It is very important that common identifiers are given due to naming requirements in subsequent commands. Only a single protein sequence can be modelled at a time.

 3. Navigate to the options at the top of the page, select `Runtime` and choose `Run all`

 4. When the modelling has been completed, on the Safari web browser (version 15.6) you will be prompted to allow the resulting file to be downloaded and selecting 'allow' will download a zip file into your downloads folder (Mac). Note if you have selected a different directory as your default directory for downloads to be sent to, you will have to change it back to the `Downloads` folder for the purpose of following this anaysis.

 5. Move the zip files from your `Downloads` folder to the directory we created earlier for protein structures, and open the files following the commands below:

```bash
mv ~/Downloads/*_g46214.zip ~/g46214_modelling_output
for file in ~/g46214_modelling_output/*.zip; do unzip "$file"; done
```

### g10577

We now want to visualise the three dimensional structure of the g10577 protein in diploids and tetraploids. We will first create a directory in our home directory to host all the protein stuctures and any modelling related output by following the command below:

```bash
mkdir -p ~/g10577_modelling_output/tetraploid_g10577_protein_images/movie
mkdir -p ~/g10577_modelling_output/diploid_g10577_protein_images/movie
```

For g10577 we had to use a different approach to model the proteins due to limitations with alphafold's memory and being unable to model the whole 1000+ amino acid long protein. Instead of putting the whole sequence into alphafold, we used the domain positions we identified in the domain identification step above (adjusted based on the literature and multiple sequence allignments) to obtain 3D structure models for the different domains in our tetraploid and diploid proteins. We decided on this method as opposed to using an alternative modelling software because we did not get biologically sensible output using software like Phyre2. This approach to modelling also requires us to obtain a complete reference protein model using swissmodel so we can essentially 'map' the domains onto the reference protein to re-build our diploid and tetraploid proteins. The protein we chose for our reference had the uniprot ID Q9LPK1 and we settled on that as a reference because of its high coverage and sequence identity (74.62%) to our reference g10577 protein sequence. The analysis can be performed by following the steps below:

 1. Follow the link to the [alphafold collab](https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb) website

 2. Input the protein sequence for a given domain in the diploid and tetraploid sequences into the `query sequence` field, and the job name should have some identifier they all share (i.e `diploid_domain1_g10577`,`tetraploid_domain1_g10577`). It is very important that common identifiers are given due to naming requirements in subsequent commands. Only a single domain from a single protein can be modelled at a time, and alphafold collab can run out of memory, so you may need multiple google accounts or use a colleagues machine to get all the domains done.

 3. Repeat step 2 for all other domains in g10577 that we idenetified in the domain identification step.  

 4. Navigate to the options at the top of the page, select `Runtime` and choose `Run all`

 5. When the modelling has been completed, on the Safari web browser (version 15.6) you will be prompted to allow the resulting file to be downloaded and selecting 'allow' will download a zip file into your downloads folder (Mac). Note if you have selected a different directory as your default directory for downloads to be sent to, you will have to change it back to the `Downloads` folder for the purpose of following this anaysis.

 6. Follow the link to the SWISS-MODEL website, select `Start Modelling`, input the reference g10577 sequence into the `Target Sequence` field and select `Build Model`. Identify the model with template `Q9LPK1.1.A`, and download the model in PDB format.

 7. Move the zip files, and reference pdb model from your `Downloads` folder to the directory we created earlier for protein structures, and open the files following the commands below:

```bash
mv ~/Downloads/*_g10577.zip ~/g10577_modelling_output
mv ~/Downloads/model_01.pdb ~/g10577_modelling_output/reference_g10577.pdb
for file in ~/g10577_modelling_output/*.zip; do unzip "$file"; done
```
## Image Generation

We have successfully modelled our proteins and now want to actually investigate the mutations in three dimensional space and create good quality images to be used in our papers/presentations. We will open the different pdbs in pyMOL, highlight domains/motifs of interest, and take snapshots of our proteins. This can be acheived following the steps outlined below.

We only want to use the best model from the alphafold output, which should have rank_001 in its name. We will load that into pymol and change the name of the object following the command in the file named below: Say what a pdb is?

The alphafold output directory has a different suffix for the directory name (i.e eebdh) so use whichever youve been given and choose your rank001 model. For the sake of running the scripts below you can only use one protein at a time. Either you are investigating the tetraploid sturcture or the diploid structure.

The alphafold output directory is named with unique identifiers such that the directory name for the modelling job is always different every time you run it. As such the paths given in our commands is generic to not cause any confusion. The output directory for each modelling job should have a number of different models ranked from 001 to 005 based on SOMETHING IDK. We will alwyas be using the rank 001 model for the purposes of investigating our protein structure and mutations.

### g46214

### Tetraploids

```bash
load ~/g46214_modelling_output/your_alphafold_tetraploid_output_directory/your_rank_001_tetraploid.pdb, tetraploid_g46214

select tetraploid_g46214_bbox_domains,tetraploid_g46214 and resi 5-108

run ~/g46214_modelling_output/colorh.py

color_h tetraploid_g46214
```
Now that the tetraploid g46214 protein has been loaded into pymol, Using your mouse or trackpad, manually adjust the orientation of the tetraploid protein to a position you are happy with. Follow the next few steps:

 1. Navigate to the header of pymol and select the plugin tab, select APBS electrostatics.

 2. Select the drop down menu in the selection entry field (selection:[       ]) and select `polymer & tetraploid_g46214` depending on which protein you are investigating. This will produce an object in PyMOL showing the electrostatic potential across the whole protein. When that has completed close the pop-up that comes afterwards.

 3. Repeat steps 1 and 2 but select `polymer & tetraploid_g46214_bbox_domains`. This will produce an object showing the electrostatic potential across the bbox domains only. Once this is completed you might have to zoom out so your whole protein is showing on the screen and the bbox domain is not being focused on.
 
NOTE: If you are using too much RAM on your machine the whole protein may be coloured white or only a few areas of your protein may be coloured by electrostatic potential. If that occurs close all other running applications on your device and rereun the steps, or restart your device and rereun the ste[s.

Now we have everything we need to produce our images and you can run the `diploid_domain_highlight.py`  script to get your nice figures that show the protein and its electrostatic potential at 90 degree angles. In the pyMOL commmand line you can should enter the command below:

```bash
run ~/path/to/python/script/tetraploid_domain_highlight.py
```

### Diploids

```bash
load ~/g46214_modelling_output/your_alphafold_diploid_output_directory/your_rank_001_diploid.pdb, diploid_g46214

select diploid_g46214_bbox_domains,diploid_g46214 and resi 5-105

run ~/g46214_modelling_output/colorh.py

color_h diploid_g46214
```

Now that the diploid g46214 protein has been loaded into pymol, Using your mouse or trackpad, manually adjust the orientation of the diploid protein to a position you are happy with. Follow the next few steps:

 1. Navigate to the header of pymol and select the plugin tab, select APBS electrostatics.

 2. Select the drop down menu in the selection entry field (selection:[       ]) and select `polymer & diploid_g46214` depending on which protein you are investigating. This will produce an object in PyMOL showing the electrostatic potential across the whole protein. When that has completed close the pop-up that comes afterwards.

 3. Repeat steps 1 and 2 but select `polymer & diploid_g46214_bbox_domains` depending on which protein you are investigating. This will produce an object showing the electrostatic potential across the bbox domains only. Once this is completed you might have to zoom out so your whole protein is showing on the screen and the bbox domain is not being focused on.
 
NOTE: if you are using too much RAM on your machine the whole protein may be coloured white or only a few areas of your protein may be coloured by electrostatic potential. If that occurs close all other running applications on your device and rereun the steps, or restart your device and rereun the steps.

Now we have everything we need to produce our images and you can run the `diploid_domain_highlight.py`  script to get your nice figures that show the protein and its electrostatic potential at 90 degree angles. In the pyMOL commmand line you can should enter the command below:

```bash
run ~/path/to/python/script/diploid_domain_highlight.py
```

### g10577

Due to the limitations with alphafold's memory we have modelled domains of our g10577 protein and not the whole protein.

We will have to load in the reference protein which was retrieved from SWISS-MODEL (uniprot ID Q9LPK1), then we will load in all the different domains of the protein and one-by-one these will be alligned to the reference protein to effectivey stitch together our original protein. 

### Tetraploid

This can be achieved by following the steps outlined below:

 1. Open a new PyMOL window and load the reference SWISSMODEL protein into PyMOL using the commands below

    ```bash
    load /path/to/your/reference/protein.pdb, reference_g10577
    ```
2. For each domain you have modelled in the Protein Structure Modelling step run these commands in the PyMOL command line to load them into PyMOL. The domain object is the name you want to call the object in PyMOL. We recommend using names such as `RNaseH_domain` which reflect the underlying biology, rather than using numbered domains. You will again be choosing the rank 001 model from the alphafold output because it gives us the best estimate at the actual protein structure modelled by alphafold
   ```bash
   load /path/to/your/diploid/or/tetraploid/Gag_domain.pdb, GAG_domain
   load /path/to/your/diploid/or/tetraploid/integrase_domain.pdb, integrase
   load /path/to/your/diploid/or/tetraploid/protease_domain.pdb, protease
   load /path/to/your/diploid/or/tetraploid/reverse_transcriptase_domain.pdb, reverse_transcriptase
   load /path/to/your/diploid/or/tetraploid/rnaseh_domain.pdb, RNaseH
   ```
3. Now we will allign each domain to the reference SWISSMODEL protein to try and map the domains and recreate a complete protein. Repeat this step for each domain until they have all been alligned to the reference

   ```bash
   allign GAG_domain, reference_g10577
   allign integrase, reference_g10577
   allign protease, reference_g10577
   allign reverse_transcriptase, reference_g10577
   allign RNaseH, reference_g10577
   ```
4. Hide the reference protein manually by navigating to the object list at the right hand side of the window. Look for the name of the reference protein which should be `reference_g10577`, select the `H` and then select `everything` from the drop down menu.

5. Colour the different domains by navigating to the object list at the right hand side of the window. Look for the name of the domain of interest, select the `C` and then select the colour based on the colour scheme outlined blow.

   GAG_domain -> sky blue
   
   integrase -> orange

   protease -> yellow

   reverse_transcriptase -> green

   RNaseH -> purple

7. Manually rotate the protein based on personal preference and take a picture of the assembled protein using the command below in the PyMOL command line
   ```bash
   png ~/path/to/ouput_image/directory/tetraploid_image.png, 3500, 3500, -1, ray=0, dpi=500
   ```
8. To get the electrostatic potential for all the domains navigate to the header of PyMOL and select the plugin tab, and then select APBS electrostatics. Select the drop down menu in the selection entry field (selection:[       ]) and select `polymer & <domain name>`. This will produce an object in PyMOL showing the electrostatic potential across the given domain. When that has completed close the pop-up that comes afterwards.

9. Repeat step 7 for all the domains we have outlined in step 2.

10. Manually rotate the protein (now coloured by electrostatic potential) based on personal preference and take a picture using the command below in the PyMOL command line
   ```bash
   png ~/path/to/ouput_image/directory/tetraploid_electrostatic_image.png, 3500, 3500, -1, ray=0, dpi=500
   ```

   
### Diploid
Due to the limitations with alphafold's memory we have modelled domains/chunks of our g10577 protein and not the whole protein.

We will have to load in the reference protein which was retrieved from SWISS-MODEL (uniprot ID Q9LPK1), then we will load in all the different domains of the diploid protein and one-by-one these will be alligned to the reference protein to effectivey stitch together our original protein. 

This can be achieved by following the steps outlined below:

 1. Open a new PyMOL window and load the reference SWISSMODEL protein into PyMOL using the commands below

    ```bash
    load /path/to/your/reference/protein.pdb, reference_g10577
    ```
2. For each domain you have modelled in the `Protein Structure Modelling` step run these commands in the PyMOL command line to load them into PyMOL. The domain object is the name you want to call the object in PyMOL. We recommend using names such as `RNaseH_domain` which reflect the underlying biology, rather than using numbered domains. You will again be choosing the rank 001 model from the alphafold output because it gives us the best estimate at the actual protein structure modelled by alphafold
   ```bash
   load /path/to/your/diploid/Gag_domain.pdb, GAG_domain
   load /path/to/your/diploid/integrase_domain.pdb, integrase
   load /path/to/your/diploid/protease_domain.pdb, protease
   load /path/to/your/diploid/reverse_transcriptase_domain.pdb, reverse_transcriptase
   load /path/to/your/diploid/rnaseh_domain.pdb, RNaseH
   ```
3. Now we will allign each domain to the reference SWISSMODEL protein to try and map the domains and recreate a complete protein. Repeat this step for each domain until they have all been alligned to the reference

   ```bash
   allign GAG_domain, reference_g10577
   allign integrase, reference_g10577
   allign protease, reference_g10577
   allign reverse_transcriptase, reference_g10577
   allign RNaseH, reference_g10577
   ```
4. Hide the reference protein manually by navigating to the object list at the right hand side of the window. Look for the name of the reference protein which should be `reference_g10577`, select the `H` and then select `everything` from the drop down menu.
  
5. Colour the different domains by navigating to the object list at the right hand side of the window. Look for the name of the domain of interest, select the `C` and then select the colour based on the colour scheme outlined blow.

   GAG_domain -> sky blue
   integrase -> orange
   protease -> yellow
   reverse_transcriptase -> green
   RNaseH -> purple
   
6. When you have manually rotated the protein how you like you take a picture of the protein or domain using the command below in the PyMOL command line
   ```bash
   png ~/path/to/ouput_image/directory/diploid_image.png, 3500, 3500, -1, ray=0, dpi=500
   ```
   
7. To get the electrostatic potential for all the domains navigate to the header of PyMOL and select the plugin tab, and then select APBS electrostatics. Select the drop down menu in the selection entry field (selection:[       ]) and select `polymer & <domain name>`. This will produce an object in PyMOL showing the electrostatic potential across the given domain. When that has completed close the pop-up that comes afterwards.

8. Repeat step 7 for all the domains we have outlined

9. Manually rotate the protein how you like and take a picture of the electrostatic potential of the domains using the command below in the PyMOL command line
   ```bash
   png ~/path/to/ouput_image/directory/diploid_electrostatic_image.png, 3500, 3500, -1, ray=0, dpi=500
   ```
## Movie Generation
We also want to have a short video showing some of our protein molecules rotating over time. The logic behind this is the same as with a flipbook. Over 360 degrees, we will rotate the protein molecule by 1 degree across a given axis, and take a picture after each rotation. Combining those individuals together will make a movie. This can be accomplished by running the scripts below: 

### g46214
You will have to open a new PyMol window and then run the commands below

### Tetraploids
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

In the pymol command line run the code below to produce a directory with a collection of tetraploid g46214 images
```bash
run ~/tetraploid_temporary_image_generation.py 
```

You will have to open a new PyMol window again and run the commands below

### Diploids
```bash
load ~/g46214_modelling_output/gatk_04AF_dip_29144/gatk_04AF_dip_29144_unrelaxed_rank_001_alphafold2_ptm_model_2_seed_000.pdb, diploid_g46214

select diploid_g46214_bbox_domains,diploid_g46214 and resi 5-105

run ~/g46214_modelling_output/colorh.py

color_h diploid_g46214
```

Now that the diploid g46214 protein has been loaded into pymol, Using your mouse or trackpad, manually adjust the orientation of the tetraploid protein to a position you are happy with. Follow the next few steps:

 1. Navigate to the header of pymol and select the plugin tab, select APBS electrostatics.

 2. Select the drop down menu in the selection entry field (selection:[       ]) and select `polymer & diploid_g46214` . This will produce an object in PyMOL showing the electrostatic potential across the whole protein. When that has completed close the pop-up that comes afterwards.

Now rotate molecule to choose a good starting position that will show everything you want when it rotates. It will be rotating up to down and this can be accomplished running the script below to produce a series of images after every rotation

In the pymol command line run the code below to produce a directory with a collection of diploid g46214 images
```bash
run ~/diploid_temporary_image_generation.py 
```


The series of images produced for the diploid and the tetraploids will be stitched together individually for each ploidy. This can be accomplished by running the script below in the command line

```bash
bash make_protein_rotation_movie.sh
```

