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


# Deselect everything so the red selection points dont show up in the output
cmd.deselect()

# unpack all the objects from the group produced when you get the electrostatic image of something
cmd.ungroup("run01")

# change name of the object showing the electrostatic potential of the protein to something more identifiable
cmd.set_name('prepared01','electrostatic_tetraploid_g46214')

# get a clear canvas
cmd.hide(representation='cartoon' ,selection='all')
cmd.hide(representation='surface' ,selection='all')

# Determine the axis of rotation for the molecule of interest to be rotated on
axis_of_rotation = 'x'

# Decide how many degrees the molecule moves along the axis of rotation axis during each iteration of the loop
degrees_to_rotate_by = 1

# Stores the path to the directory where you want to store the images. 
tetraploid_directory_to_store_frames = '~/g46214_modelling_output/tetraploid_g46214_protein_images/movie'

# Electrostatic bbox images


# Tetraploid
# Show the electrostatic bbox2 domain of the tetraploid protein 
cmd.show(representation='surface', selection='electrostatic_tetraploid_g46214')

# Highlight the mutations in the tetraploid
# cmd.select('mutations','tetraploid_g46214 and resi 207+153')

# cmd.color('blue','mutations')

# cmd.show('sticks','mutations')

# Do electrostatic potential rotation first so that the bar at the bottom showing the elctrostatic potential is not in all clips
# Rotate the two bbox domains again, showing the elecostratic potential across them
for i in range(360):

    # The 
    cmd.rotate(axis = axis_of_rotation, angle = degrees_to_rotate_by, object='electrostatic_tetraploid_g46214') 

    # refresh the screen
    cmd.refresh()

    # Save each frame in each iteration of the loop into a directory of your choosing, the 03d essentially makes it so that numbers are written out as 001.png instead of just 1.png
    # It becomes useful later on when we are trying to stitch the images back together using the ffmpeg command line tool
    # Add 359 because rotating the first protein molecule ended at 359
    cmd.png("{}/frame_{:03d}.png".format(tetraploid_directory_to_store_frames,i+(359*2)),width=1000,height=1000,dpi=200)


# Hide the cartoon representation of the bbox domains
cmd.hide(representation='surface',selection='electrostatic_tetraploid_g46214')


# Delete the bar at the bottom showing the elctrostatic potential
cmd.delete("apbs_ramp01")

# Full Protein Molecule

# show the tetraploid protein molecule
cmd.show(representation='cartoon',selection='tetraploid_g46214')
cmd.show(representation='sticks',selection='mutation1')
cmd.show(representation='sticks',selection='mutation2')
# Go through a 360 rotation
for i in range(360):

    # Rotate the object of interest each iteration of the loop. The object is the protein/molecule you are interested in and want to rotate.
    cmd.rotate(axis = axis_of_rotation, angle = degrees_to_rotate_by, object='tetraploid_g46214') 

    # Refresh the screen
    cmd.refresh()

    # Save each frame in each iteration of the loop into a directory of your choosing, the 03d essentially makes it so that numbers are written out as 001.png instead of just 1.png
    # It becomes useful later on when we are trying to stitch the images back together using the ffmpeg command line tool
    cmd.png("{}/frame_{:03d}.png".format(tetraploid_directory_to_store_frames,i+359),width=1000,height=1000,dpi=200)

# Hide the tetraploid g46214 protein 
# cmd.hide(representation='cartoon',selection='tetraploid_g46214')


# TETRAPLOID PROTEIN SETUP
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
cmd.select('tetraploid_bbox1_and_bbox2_domain','tetraploid_g46214 and resi 5-108')

cmd.color('white', 'tetraploid_bbox1_and_bbox2_domain')

# Highlight the mutations in the tetraploid
# cmd.select('mutations','tetraploid_g46214 and resi 207+153')

cmd.color('blue','mutation1')
cmd.color('chocolate','mutation2')

# cmd.show('sticks','mutations')

# Deselect everything so the red selection points dont show up in the output
cmd.deselect()

## domains only with no electrostatics or hydrophobic
# tetraploid

# Show the  domains of the tetraploid g46214 protein only
cmd.show(representation='cartoon',selection='tetraploid_g46214')

# rotate the two bbox domains across the x axis.
for i in range(360):

    # the object is the protein/molecule you are interested in and want to rotate. you
    cmd.rotate(axis = axis_of_rotation, angle = degrees_to_rotate_by, object='tetraploid_g46214')

    # refresh the screen
    cmd.refresh()

    # Save each frame in each iteration of the loop into a directory of your choosing, the 03d essentially makes it so that numbers are written out as 001.png instead of just 1.png
    # It becomes useful later on when we are trying to stitch the images back together using the ffmpeg command line tool
    # Add 359 because rotating the first protein molecule ended at 359
    cmd.png("{}/frame_{:03d}.png".format(tetraploid_directory_to_store_frames,i),width=1000,height=1000,dpi=200)

# Hide the two bbox domains
# cmd.hide(representation='cartoon',selection='tetraploid_g46214_bbox_domains')


