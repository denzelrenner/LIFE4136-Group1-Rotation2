#!/bin/bash
#SBATCH --job-name=phylogenetic_tree
#SBATCH --partition=hpc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=30g
#SBATCH --time=24:00:00
#SBATCH --output=~/%x.out
#SBATCH --error=~/%x.err

# source home profile 
source $HOME/.bash_profile

# activate environment with all tools we need 
conda activate raxml_tree_build

# make directory for raxml output
mkdir -p ~/raxml_output

# move into the directory where we will put all the output
cd ~/raxml_output

# build our tree using LG as model of protein sequence evolution, and using Raphinus sativus as the outgroup
raxml-ng \
	--msa ~/g46214_allignment.fasta \
	--threads 8 \
	--prefix g46214_phylogeny_output \
	--model LG \
	--all \
	--tree rand{10},pars{10} \
	--bs-trees 1000 \
	--outgroup XP_018447019.1

# get the verison number of raxml that we are using 
raxml-ng --version

echo "The Job ID for this job is: $SLURM_JOB_ID"

