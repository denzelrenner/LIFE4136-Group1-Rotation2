# This script allows snapshot images to be taken of the protein of interest. It takes a snapshot image every 90 degrees of rotation.
# Before running this script you will have to load your protein into pymol manually rotate your protein using the trackpad or mouse to get a
# good starting position you will be happy with 

# Example run command to be entered into the pymol command line
# run /path/to/python/file

# import required modules
from pymol import cmd
import os

# Highlight the mutations in the tetraploid
cmd.select('mutation1','tetraploid_g46214 and resi 207')
cmd.select('mutation2','tetraploid_g46214 and resi 153')

# cmd.show(representation='sticks',selection='mutations')
# Deselect everything so the red selection points dont show up in the output
cmd.deselect()

# Hide the cartoons showing just the normal protein
cmd.hide(representation='cartoon' ,selection='all')
# cmd.hide(representation='sticks',selection='mutations')

# Remove all objects from the electrostatic group. This step allows us to manipulate the output from running APBS electrostatics better
cmd.ungroup("run01")
cmd.ungroup("run02")

# change name of the object showing the electrostatic potential of the protein to something more identifiable
cmd.set_name('prepared01','tetraploid_g46214_electrostatics')
cmd.set_name('prepared02','tetraploid_g46214_bbox_domains')


# Set a new name for the images produced for showing the electrostatic potential across the protein
output_file_prefix1 = '~/g46214_modelling_output/tetraploid_g46214_protein_images/tetraploid_whole_protein_electrostatic'
    
# Do the same but for the visualisation of the electrostatic potential across the protein
for i in range(0,360,90):

    # Rotate the molecule by 90 degrees. replace object with the name of the protein in pymol 
    cmd.rotate(axis = 'y', angle = 90, object='tetraploid_g46214_electrostatics') 

    # Refresh the screen
    cmd.refresh()

    # Create an image of the protein
    cmd.png(f"{output_file_prefix1}_{i}.png",width=3500,height=3500,dpi=500,ray=0)


# Remove the bar at the bottom showing the electrostatically 
cmd.delete("apbs_ramp01")

# When we ungrouped the 
cmd.hide(representation='surface' ,selection='tetraploid_g46214_electrostatics')


# Show only the bboxes which are important to bbox protein transcriptional activity
output_file_prefix2 = '~/g46214_modelling_output/tetraploid_g46214_protein_images/tetraploid_bbox_domains_electrostatic'

# Go from 0 to 360 in increments of 90. This means the loop is only ran four times and so the protein is only rotated 4 times, and the output file
# has how much it was rotated by in the name 
for i in range(0,360,90):

    # Rotate the molecule by 90 degrees. replace object with the name of the protein in pymol 
    cmd.rotate(axis = 'y', angle = 90, object='tetraploid_g46214_bbox_domains') 

    # Refresh the screen
    cmd.refresh()

    # Create an image of the protein
    cmd.png(f"{output_file_prefix2}_{i}.png",width=3500,height=3500,dpi=500,ray=0)

# remove electrostatic potential thing from the bottom
cmd.delete("apbs_ramp02")

# Hide the bbox domain electrostatics
cmd.hide(representation='surface' ,selection='tetraploid_g46214_bbox_domains')

# Show only the hydrophobicity of the protein
cmd.show(representation='cartoon' ,selection='tetraploid_g46214')
cmd.show(representation='sticks',selection='mutation1')
cmd.show(representation='sticks',selection='mutation2')

# Store the path for the output file in a variable. For the output file name, do not give the .png file extension, only the prefix is needed
output_file_prefix3 = '~/g46214_modelling_output/tetraploid_g46214_protein_images/tetraploid_whole_protein_hydrophobicity'

# Go from 0 to 360 in increments of 90. This means the loop is only ran four times and so the protein is only rotated 4 times, and the output file
# has how much it was rotated by in the name 
for i in range(0,360,90):

    # Rotate the molecule by 90 degrees. replace object with the name of the protein in pymol 
    cmd.rotate(axis = 'y', angle = 90, object='tetraploid_g46214') 

    # Refresh the screen
    cmd.refresh()

    # Create an image of the protein
    cmd.png(f"{output_file_prefix3}_{i}.png",width=3500,height=3500,dpi=500,ray=0)

# Colour the protein by hydrophobicity 
# cmd.color_h('tetraploid_g46214')

# We wil begin by selecting domains in our protein to highlight and colour differently 

# Colour the whole protein a default colour
cmd.color('red', 'tetraploid_g46214')

# Highlight the m6 motif whose function is currently unknown 

cmd.select('tetraploid_m6_motif', 'tetraploid_g46214 and resi 176-189')

cmd.color('yellow', 'tetraploid_m6_motif')

# Highlight the m7 motif whose function is  currently unknown

cmd.select('tetraploid_m7_motif', 'tetraploid_g46214 and resi 103-114')

cmd.color('magenta', 'tetraploid_m7_motif')

# Highlight the VP pair near the C terminus, conserved in all bbox proteins and has been thought to be useful for their protein interactions (needs experimental validation)

cmd.select('tetraploid_VP_pair', 'tetraploid_g46214 and resi 303-304')

cmd.color('green', 'tetraploid_VP_pair')

# Highlight the non_classical nuclear localisation signal at the C-terminus. If it functions similar to bbx21 then it needs to be recognised by importins and sequestered to the nucleus to function (needs experimental validation to see if KKKT is a valid NLS)

cmd.select('non_classical_NLS', 'tetraploid_g46214 and resi 312-315')

cmd.color('cyan' , 'non_classical_NLS')

# select the first and second bbbox domain
# cmd.select('tetraploid_bbox1_domain','tetraploid_g46214 and resi 5-47')

# cmd.color('blue', 'tetraploid_bbox1_domain')

# # Select the second bbox domain
# cmd.select('tetraploid_bbox2_domain', 'tetraploid_g46214 and resi 63-105')

# cmd.color('white', 'tetraploid_bbox2_domain')

cmd.select('tetraploid_bbox1_and_bbox2_domain','tetraploid_g46214 and resi 5-108')

cmd.color('white', 'tetraploid_bbox1_and_bbox2_domain')

# Change the colour of the mutations

cmd.color('blue','mutation1')
cmd.color('chocolate','mutation2')

# cmd.show('sticks','mutations')

# Deselect everything so the red selection points dont show up in the output
cmd.deselect()

# Store the path to the output file in a variable
output_file_prefix4 = '~/g46214_modelling_output/tetraploid_g46214_protein_images/tetraploid_whole_protein'

for i in range(0,360,90):

    # Rotate the molecule by 90 degrees. replace object with the name of the protein in pymol 
    cmd.rotate(axis = 'y', angle = 90, object='tetraploid_g46214') 

    # Refresh the screen
    cmd.refresh()

    # Create an image of the protein
    cmd.png(f"{output_file_prefix4}_{i}.png",width=3500,height=3500,dpi=500,ray=0)









