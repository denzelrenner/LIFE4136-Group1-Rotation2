#!/bin/bash
#SBATCH --job-name=downloading_raxml
#SBATCH --partition=hpc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=30g
#SBATCH --time=24:00:00
#SBATCH --output=/workhere/students_2023/Group_2_VDJ/Output_and_Error/%x.out
#SBATCH --error=/workhere/students_2023/Group_2_VDJ/Output_and_Error/%x.err

# source home profile 
source $HOME/.bash_profile

# create the environment with the name we want
conda create -n raxml_tree_build -c bioconda python=3.7 

# activate the environment 
conda activate raxml_tree_build

# install raxml
conda install -c bioconda raxml-ng

# deactivate env
conda deactivate 

# look at job id 
echo "The Job ID for this job is: $SLURM_JOB_ID"