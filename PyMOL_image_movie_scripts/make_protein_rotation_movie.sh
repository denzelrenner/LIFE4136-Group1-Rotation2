#!/bin/bash

# This script takes the series of png images created using the diploid and tetraploid python scripts and stitches them together to make a movie of the protein rotating

# Tetraploid Movie 

# Using the output from the temporary image generation python script, specify the series of pngs with -i. 
# %3d allows us to specify the pattern of the input images where they have the 001 numbering format 
# give the name and path to the output file in the last line of the command

ffmpeg -y \
    -i ~/g46214_modelling_output/tetraploid_g46214_protein_images/movie/frame_%3d.png \
    -c:v libx264 \
    -pix_fmt yuv420p \
    -r 60 \
    -vf "drawtext=text='Tetraploid g46214':x=0:y=5:fontsize=30:fontcolor=white" \
    ~/g46214_modelling_output/tetraploid_g46214_protein_images/movie/tetraploid_g46214_rotation_clip.mp4

# Remove all the individual frames from the tetraploid movie directory so it is not cluttered
rm ~/g46214_modelling_output/tetraploid_g46214_protein_images/movie/frame*.png

# Diploid Movie
ffmpeg -y \
    -i ~/g46214_modelling_output/diploid_g46214_protein_images/movie/frame_%3d.png \
    -c:v libx264 \
    -pix_fmt yuv420p \
    -r 60 \
    -vf "drawtext=text='Diploid g46214':x=0:y=5:fontsize=30:fontcolor=white" \
    ~/g46214_modelling_output/diploid_g46214_protein_images/movie/diploid_g46214_rotation_clip.mp4

# Remove all the individual frames from the diploid movie directory so it is not cluttered
rm ~/g46214_modelling_output/diploid_g46214_protein_images/movie/frame*.png

