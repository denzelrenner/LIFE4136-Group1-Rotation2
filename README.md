# LIFE4136_Group1_Rotation2
This repositry contained all the Scripts and data needed to reproduce results from group1 in rotation2 of the groupwork projects.

## What is the problem we have been presented with?
We have been provided with some genes in cochaleria 

## What are the aims of our analysis?
We were aiming to look at the structure and function of different genes found in a selection scan

## What are the expected outcomes?
We expect to 

# Prerequisites
## Required files and data
All these files should be downloaded into the same directory before following the rest of this document.

For replicating results on the cloud :
1)VCFs directory containing {UK_scan_dips.vcf and UK_scan_tets.vcf}
2)Reference_Genome directory containing {C_excelsa_V5_braker2_wRseq.gff3, C_excelsa_V5.fasta and C_excelsa_V5.fasta.fai}

For replicating results from your local machine
1)

Explain what fai is and put a link as well

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

InterPro version 98.0 (https://www.ebi.ac.uk/interpro/)

ITOL:Interactive Tree Of Life version 6.8.2 (https://itol.embl.de/)

Conda version 23.5.2

all performed on mac intel i5

## Tool intallation 
These scripts/commands should be ran to install all the tools necessary to reproduce results

On the cloud:
raxml_install.sh


On your local machine:
brew install ffmpeg (Mac Intel i5)
Navigate to the Jalview,PyMOL websites and follow the download link for your machine

## Script description
`raxml_install.sh` - installs raxml-ng

`consenus.sh` - produces a consensus nucleotide sequence (in fasta format) for the reference, diploid, and tetraploids.


# THE ANALYSIS

## Consensus Sequences (Nucleotide)
The first thing we will do is get consensus sequences for our diploids and tetraploids. We will filter the vcfs we have to only include biallelic variants with an allele frequency greater than 0.4. The resulting filtered vcfs will then be used, along with the reference fasta and gff, to produce our consensus sequences for the genes of interest only and not the entire scaffold. This is accomplished by running the script below.

`<insert_script>`

## Protein Sequences 
Amongst the resulting files there should be fasta files containing the entire consensus gene sequence for g46214 and g10577 in our diploids and tetraploids. Follow the link to the orf finder website and input the consensus sequences for each gene in the different contrasts.

1.Input nucleotide sequence from g42614 diploid into the query sequence field and submit the job. From the resulting output stitch the protein sequence for the three longest open reading frames to form the full length protein(ORF1+ORF3+ORF4) (results can be achieved if only the coding sequence is taken and not the entire gene. Methodology not shown) 

2.Input nucleotide sequence from g42614 tetraploid into the query sequence field and submit the job. From the resulting output stitch the protein sequence for the three longest open reading frames to form the full length protein(ORF1+ORF3+ORF4) (results can be achieved if only the coding sequence is taken and not the entire gene. Methodology not shown) 

3.Input nucleotide sequence from g10577 diploid into the query sequence field and submit the job. From the resulting output take the longest open reading frame which is the full length protein (ORF1). If you take only coding sequence is the tetraploid truncated?

4.Input nucleotide sequence from g10577 tetraploid into the query sequence field and submit the job. From the resulting output take the longest open reading frame which is the full length protein (ORF1).

## Homolog identification
In order to find homologs blastp search on the tair database and on the ncbi 


## Domain identification
Now we have found out the closest homologs for our proteins. We now want to verify the functional domains within our proteins. To this we followed these steps 

1.Follow the link to the InterPro website

2.Input the protein sequences (in fasta format) into the query field labelled `Enter your sequence` and choose `search`

3)

## Multiple Sequence Allignments and Phylogenetic Tree Building 

### g46214
Homologous proteins to be used in pyhlogenetic tree building and mutliple sequence allignments were selected based on having 100% query cover and >40% percentage identity to our reference protein sequence. The protein sequences for the selected homologs were manually retrived by selecting their ncbi dataset accession code in the blastp results page, then selecting the 'FASTA' option at the top of the page, and finally copy and pasting the protein sequence into the fasta file. The consensus protein sequences for the reference, diploid and tetraploid protein sequences should also be added to these files in fasta format.

By following the steps below, a multiple sequence allignment of the protein sequences was generated.

1.Follow the link to the uniprot website and navigate to the tab labelled `Align`

2.Copy and paste all the sequences from the fasta file you generated above into the query field then click the `Run align` button

3.On the results page, click on the `download` option and choose `FASTA` format from the drop-down menu before downloading the allignment.

To build a phylogenetic tree from this allignment we will use raxml by running the script below. This will produce a file with the extension `.besttree` which is a file that . 

`<insert_script>`

We will view the best tree file using the ITOL website. To do this the steps below were followed

1)Follow the link provided to the ITOL website 

2)Click on the `Upload a tree` option

3)Navigate to the `Tree file` field and select the `choose file` option. Search through your files to locate the raxml file with the `.besttree` extension and upload it.

### g10577
1.Made a new fasta file of a multiple sequence alignment using the top BLASTp hits with the g10577_tetraploid_GATK_>0.75 amino acid sequence

2.Uploaded the fasta file into MEGA11: Align -> Edit/Build Alignment -> Retrieve sequences from a file -> `input.fasta`

3.Produced a multiple sequence alignment by navigating to the tool bar and selecting Alignment: Alignment -> Align by ClustalW -> select all -> used default parameters except for Delay Divergent Cutoff (%) and selected 45%.

4.Next, selected Data -> Phylogenetic Analysis from the tool bar then produced a phylogenetic tree by selecting Phylogeny from the toolbar -> Construct/Test Neighbour-Joining Tree with the following parameters
  
    1.Test of Phylogeny: Bootstrap method
    
    2.No. Bootstrap Replications: 500
    
    3.Substitutions type: Amino acid
    
    4.Model/Method: Poisson model
    
    5.Rates among Sites: Gamma Distributed (G)
    
    6.Gamma Parameter: 1.00
    
    7.Pattern among Lineages: Same (Homogeneous)
    
    8.Gaps/Missing Data Treatment: Pairwise Deletion
    
    9.Number of Threads: 3

5.This produced an NJ tree with default layout – I customised the tree with the following
  
    1.Taxon names -> Font -> Arial -> Bold Italic -> 10
    
    2.Layout -> Toggle Scaling of the Tree + Auto-size Tree
    
    3.I customised the length and the width of the tree to make it more aesthetic

6.To save a PNG file of the tree I selected Image from the toolbar: Image -> Save a PNG file  then gave the file path to my Desktop with an appropriate file name.


## Protein Structure Modelling
Say what you tried

Create a directory in your home directory (on your local machine) to host all the protein stuctures and any modelling related output by following the command below

```bash
mkdir ~/g46214_modelling_output
```



To obtain 3D strucutre models for our proteins we input the sequence .

1.Follow the link to the alphafold collab website provided above

2.Input the protein sequence for reference,diploid,and tetraploid sequences into the `query sequence` field, and enter reference_g46214,diplpid_g46214,and tetraploid_g46214 in the `job name` field respectively. Note it is very important these exact job names are given due to naming requirements in subsequent scripts, and only a single protein sequence can be modelled at a time.

3.Navigate to the options at the top of the page, select `Runtime` and choose `Run all`

4.When the modelling has been completed, on the Safari web browser (version 15.6) you will be prompted to allow the resulting file to be downloaded and selecting 'allow' will download a zip file into your downloads folder (Mac)

5.Move the zip files from your downloads folder to the directory we created earlier for protein structures, and open the files following the command below

```bash
mv ~/Downloads/*_g46214.zip ~/g46214_modelling_output
unzip *_g46214.zip
```

## Image and Movie generation
We have successfully modelled our proteins and now want to create good quality images to be used in our papers/presentations.

We also want to have a short video showing some of our protein molecules rotating over time. This can be accomplished by running the scripts below.

This script will produce a series of images
<insert_script>

These series of images are then stitched together to produce a movie of a rotating protein
m













# Tool References
These are references (where applicable) for all the tools used in our analysis
