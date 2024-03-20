# LIFE4136_Group1_Rotation2
This repositry contains all the ccripts and data needed to reproduce results from group1 in rotation2 of the groupwork projects.

Note: Unless on the cloud or stated otherwise, all command line code was ran on a Mac Intel i5 or Lukes Mac

## What is the problem we have been presented with?
We have been provided with some genes in Cochaleria that are under selection during the stabalisation of polyploidy where the plants are moving from having diploid genomes to tetraploid ones

## What are the aims of our analysis?
We were aiming to generate structural predictions for candidate genes (from the selection scan) in both diploids and tetraploids, then compare these structures as well as their primary sequences, to determine differences that might relate to differnet function in diploids and tetraploids.

## What are the expected outcomes?
Depending on the gene and given our understanding of the system, we expect to see mutations in the tetraploids that allow it to survive certain conditions better than the diploids 

# Prerequisites
## Required files and data
All these files should be downloaded into the same directory before following the rest of this document.

For replicating results on the cloud :

* `VCFs` directory containing `UK_scan_dips.vcf` and `UK_scan_tets.vcf`

* `Reference_Genome` directory containing `C_excelsa_V5_braker2_wRseq.gff3`, `C_excelsa_V5.fasta` and `C_excelsa_V5.fasta.fai`

For replicating results from your local machine

1)

Explain what fai is and put a link as well, and also .dict A more detailed explanation can be found here: (https://gatk.broadinstitute.org/hc/en-us/articles/360035531652-FASTA-Reference-genome-format)

## Tool versions and links
These are all the tools that were used in our analysis with versions and links provided where applicable.

ncbi blast (https://blast.ncbi.nlm.nih.gov/)

TAIR BLAST version 2.9.0+ (uses NCBI BLAST 2.9.0+) (https://www.arabidopsis.org/Blast/index.jsp)

orf finder (https://www.ncbi.nlm.nih.gov/orffinder/)

uniprot (https://www.uniprot.org/align), clustalO version 1.2.4 (clustal was ran through the align function on the uniprot website)

alphafold collab (ColabFold v1.5.5: AlphaFold2 using MMseqs2

(https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb)

InterPro version 98.0 (https://www.ebi.ac.uk/interpro/)

ITOL:Interactive Tree Of Life version 6.8.2 (https://itol.embl.de/)

MEGA version 11 (https://www.megasoftware.net/)

jalview version 2.11.3.2

PyMOL version 2.5.8

Homebrew version 4.2.10 (Mac Intel i5)

ffmpeg version 6.1.1 

GATK verison 4.2.2.0, (HTSJDK version 2.24.1, Picard version 2.25.4)

samtools version 1.19.2, (htslib version 1.19.1)

RAxML-NG version 0.9.0

Conda version 23.5.2



all performed on mac intel i5

## Tool intallation 
These scripts/commands should be ran to install all the tools necessary to reproduce results

On the cloud:
raxml_install.sh


On your local machine:

how to install homebrew

brew install ffmpeg 

Navigate to the Jalview,PyMOL websites and follow the download link for your machine

## Script description
`raxml_install.sh` - installs raxml-ng

`consenus.sh` - produces a consensus nucleotide sequence (in fasta format) for the reference, diploid, and tetraploids.


# THE ANALYSIS

## Consensus Sequences (Nucleotide)
The first thing we will do is get consensus sequences for our diploids and tetraploids. We will filter the vcfs we have to only include biallelic variants with an allele frequency greater than 0.49 (to include allele frequencies of 0.5) . The resulting filtered vcfs will then be used, along with the reference fasta and gff, to produce our consensus sequences for the genes of interest. This is accomplished by following the steps below:

```bash
conda activate /shared/apps/conda/bio2

bash `<insert_script>`
```
## Protein Sequences 
Amongst the files produced from running the <insert_script> there should be fasta files containing the entire consensus coding sequence for g46214 and g10577 in our diploids, tetraploids, and reference. Now that we have the nucleotide sequences, we can translate them to get our protein sequences. This can be accomplished by following the steps below.

  1.Follow the link to the orf finder website.

  2.Input the nucleotide seqeunce into the query seqeunce field and submit the job. On the output page, copy and paste the amino acid seqeuence for ORF1 (the longest open reading frame) into a new fasta file with any identifiable headers you prefer.


## Homolog identification

Now that we have retrieved the sequences for our proteins, we want to figure out what they actually are, and what their closest homologs are in other plant or animal species. You should create fasta files called `g10577_homologs.fasta` and `g46214_homologs.fasta` to store all homologs selected based on a given criteria outlined below.

For g46214, homologous proteins were selected based on having 100% query cover and >60% percentage identity to the reference protein sequence. 
For g10577, homologous proteins were selected based on >70% query cover and >40% percentage identity to the reference protein sequence

To identify g46214 and g10577 homologs in other plant species we followed the steps below:

 1. Follow the link to the NCBI BLAST webpage

 2. Select `Protein BLAST` 

 3. Input the protein sequence for the g10577 and g46214 reference proteins sequences into the `Enter Query Sequence` field, and select `BLAST` at the bottom of the page

 4. Manually retrieve the protein sequence of the homologs from the blastp results page by selecting their ncbi dataset accession code (i.e `XP_018447019.1`) in the blastp results page

 5. Select the `FASTA` option at the top of the page, and finally copy and paste the protein sequence into the `g10577_homologs.fasta` or `g46214_homologs.fasta` fasta files

 6. The consensus protein sequences for the reference, diploid and tetraploid protein sequences should also be added to the homolog files in fasta format.

Note that we also searched for homologs in the model species Arabidopsis thaliana on the TAIR website, but for g10577 the identified homolog did not make biological sense when investigated further through multiple sequence allignments and was not we suspect it was not that.

## Domain identification
We have found out the closest homologs for our proteins in different species so we can begin to start investigating our proteins function, and structural domains. We now want to verify the functional domains within our proteins. To this we followed these steps 

 1. Follow the link to the InterPro website

 2. Input the protein sequences (in fasta format) into the query field labelled `Enter your sequence` and choose `search`

 3. Manually insert the positions of the domains, and any other relevant information, into a txt file of your choosing so we have the exact coordinates of the different domains in the protein

## Multiple Sequence Allignments and Phylogenetic Tree Building 

This is the last step where we remain at the primary sequence level. We have sucessfully determined the important domains in our proteins, as well as their homlogs across different species. The next step is to identify conserved residues in our proteins and determine any important mutations between our diploid, tetraploids,and reference proteins, by comparing the sequences of our proteins with their close homologs. We will also build phylogenetic trees as another form of visualising and representing how much the proteins in Cochalearia have deviated from their homologs, and also how much the tetraploids have deviated from the diploids 

### g46214
 
By following the steps below, a multiple sequence allignment of the protein sequences was generated.

1.Follow the link to the uniprot website and navigate to the tab labelled `Align`

2.Copy and paste all the sequences from the fasta file you generated above into the query field then click the `Run align` button

3.On the results page, click on the `download` option and choose `FASTA` format from the drop-down menu before downloading the allignment. When prompted to enter a file name, enter `g46214_allignment.fasta`

To build a phylogenetic tree from this allignment we will use raxml by running the script below. This will produce a file with the extension `.besttree` which is the unrooted best maximum likelihood tree.  

`<insert_script>`

We will view the best tree file using the ITOL website. To do this the steps below were followed

  1.Follow the link provided to the ITOL website 

  2.Click on the `Upload a tree` option

  3.Navigate to the `Tree file` field and select the `choose file` option. Search through your files to locate the raxml file with the `.besttree` extension and upload it.

  4.Replace the accession code names at the tip of the tree by clicking on the accession code (i.e XP_009106205.1), navigating to `label`, then select `edit label` and replace it with its Genus and Species

### g10577

1.Made a new fasta file of a multiple sequence alignment using the top BLASTp hits with the g10577_tetraploid_GATK_>0.75 amino acid sequence

2.Uploaded the fasta file into MEGA11: Align -> Edit/Build Alignment -> Retrieve sequences from a file -> `g10577_homologs.fasta`

3.Produced a multiple sequence alignment by navigating to the tool bar and selecting Alignment: Alignment -> Align by ClustalW -> select all -> used default parameters except for Delay Divergent Cutoff (%) and selected 45%.

4.Next, selected Data -> Phylogenetic Analysis from the tool bar then produced a phylogenetic tree by selecting Phylogeny from the toolbar -> Construct/Test Neighbour-Joining Tree with the following parameterss

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

5.This produced an NJ tree with default layout – I customised the tree with the following

```
1.Taxon names -> Font -> Arial -> Bold Italic -> 10
    
2.Layout -> Toggle Scaling of the Tree + Auto-size Tree
    
3.I customised the length and the width of the tree to make it more aesthetic
```

6.To save a PNG file of the tree I selected Image from the toolbar: Image -> Save a PNG file  then gave the file path to my Desktop with an appropriate file name.


## Protein Structure Modelling

### g46214

We now want to visualise the three dimensional structure of our proteins. We will first create a directory in our home directory to host all the protein stuctures and any modelling related output by following the command below:

```bash
mkdir -p ~/g46214_modelling_output/tetraploid_g46214_protein_images/movie
mkdir -p ~/g46214_modelling_output/diploid_g46214_protein_images/movie

mkdir -p ~/g10577_modelling_output/tetraploid_g10577_protein_images/movie
mkdir -p ~/g10577_modelling_output/diploid_g10577_protein_images/movie
```

To obtain 3D structure models for our proteins we followed the steps below:

1.Follow the link to the alphafold collab website provided above

2.Input the protein sequence for reference,diploid,and tetraploid sequences into the `query sequence` field, and for a given candidate gene , the job name should have some identifier they all share (i.e `diploid_g46214`,`tetraploid_g46214`,`reference_g46214`) where the `*` is any extra information you want to add. Note it is very important that common identifiers are given due to naming requirements in subsequent scripts, and only a single protein sequence can be modelled at a time.

3.Navigate to the options at the top of the page, select `Runtime` and choose `Run all`

4.When the modelling has been completed, on the Safari web browser (version 15.6) you will be prompted to allow the resulting file to be downloaded and selecting 'allow' will download a zip file into your downloads folder (Mac). Note if you have selected a different directory as your default directory for downloads to be sent to, you will have to change it back to the `Downloads` folder for the purpose of following this anaysis.

5.Move the zip files from your `Downloads` folder to the directory we created earlier for protein structures, and open the files following the commands below

```bash
mv ~/Downloads/*_g46214.zip ~/g46214_modelling_output
unzip ~/g46214_modelling_output/*_g46214.zip
unzip ~/g46214_modelling_output/*tetraploid*.result.zip
unzip ~/g46214_modelling_output/*diploid*.result.zip
unzip ~/g46214_modelling_output/*reference*.result.zip
```

## Image and Movie generation
We have successfully modelled our proteins and now want to actually investigate the structure and create good quality images to be used in our papers/presentations. We will open the different pdbs in pyMOL, highlight domains/motifs of interest, and take snapshots of our proteins. This can be acheived following the steps outlined below.

We only want to use the best model from the alphafold output, it should have rank_001 in its name. We will load that into pymol and change the name of the object following the command in the file named below:

The alphafold tetraploid directory has a different suffix for the directory name (i.e eebdh) so use whichever youve been given and choose your rank001 model. You will also need to copy and paste or download the code for a script to colour proteins by hydrophobicity on the PyMOL wiki. Read more about it here (https://pymolwiki.org/index.php/Color_h)

```bash
load ~/g46214_modelling_output/your_alphafold_tetraploid_output_directory/your_rank_001_tetraploid.pdb, tetraploid_g46214

select tetraploid_g46214_bbox_domains,tetraploid_g46214 and resi 5-108

load ~/g46214_modelling_output/your_alphafold_diploid_output_directory/your_rank_001_diploid.pdb, diploid_g46214

select diploid_g46214_bbox_domains,diploid_g46214 and resi 5-105

run ~/g46214_modelling_output/colorh.py

color_h diploid_g46214

color_h tetraploid_g46214

```

Now the diploid and tetraploid g46214 proteins have been loaded into pymol. We alligned them so the rotations produce good output. Using your mouse or trackpad, manually adjust the orientation of either the diploid/tetraploid protein to a position you are happy with. Follow the next few steps:

 1. Navigate to the header of pymol and select the plugin tab, select APBS electrostatics.

 2. Select the drop down menu in the selection entry field (selection:[       ]) and select `polymer & tetraploid_g46214` to produce an object showing the electrostatic potential  across the whole tetraploid protein only. When that has completed a pop up will ask you to close something and you shoudl select yes.

 3. Repeat steps 1 and 2 but select {`polymer & tetraploid_g46214_bbox_domains`,`polymer & tetraploid_g46214`,`polymer & diploid_g46214_bbox_domains`} to be able to produce an   object showing the electrostatic potential across the bbox domains only.
 
NOTE: if you do not have sufficient RAM the electrostatic potential will only be coloured across portions of the protein. If that occurs close all other running applications on your device and rereun the steps.

Now we have everything we need to produce our images and you can run the `image_generation.py` script to get your nice figures that show the protein and its electrostatic potential at 90 degree angles. In the pyMOL commmand line you can should enter the command below"

```bash
run ~/path/to/python/script/image_generation.py
```


We also want to have a short video showing some of our protein molecules rotating over time. The logic behind this is the same as with a flipbook. Over 360 degrees, we will rotate the protein molecule by 1 degree across a given axis, and take a picture after each rotation. Combining those together will make a movie. This can be accomplished by running the scripts below:


```bash
load ~/g46214_modelling_output/gatk_04AF_tetraploid_b50ff/gatk_04AF_tetraploid_b50ff_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb, tetraploid_g46214
select tetraploid_g46214_bbox_domains,tetraploid_g46214 and resi 5-108

load ~/g46214_modelling_output/gatk_04AF_dip_29144/gatk_04AF_dip_29144_unrelaxed_rank_001_alphafold2_ptm_model_2_seed_000.pdb, diploid_g46214
select diploid_g46214_bbox_domains,diploid_g46214 and resi 5-105

run ~/g46214_modelling_output/colorh.py

colorh tetraploid_g46214
```

Now do the select thing

Now rotate molecule to choose a good starting position that will show everything you want when it rotates. It will be rotating up to down

This script will produce a series of images
<insert_script>

These series of images are then stitched together to produce a movie of a rotating protein
m













# Tool References
These are references (where applicable) for all the tools used in our analysis
