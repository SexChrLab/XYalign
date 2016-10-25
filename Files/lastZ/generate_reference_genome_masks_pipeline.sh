#!/bin/bash

. /u/local/Modules/default/init/modules.sh
module load R
module load python
module load bedtools

cwd=$(pwd)

### To run, modify these paths
lastz_dir=$cwd/lastz-distrib/bin #This is the directory where to call lastZ
analysis_dir=$cwd/generate_chrY_masks #This is the directory where all the input and output files are stored. Assume that the fasta file for a specific chromosome is stored here. 
scripts_dir=$cwd/scripts #This is the directory where all of the scripts are stored
chr='chrY'
xKb_windowSize=10000

################################################################################################
# STEP 01: To run lastZ for self-self alignment
################################################################################################

$lastz_dir/lastz $analysis_dir/$chr".fa" --self --notransition --ambiguous=iupac --nogapped --nomirror --step=10 --exact=50 --notrivial --format=rdotplot > $analysis_dir/$chr"_"$chr".rdotplot"

################################################################################################
# STEP 02: Dot plot for visualization
################################################################################################

Rscript $scripts_dir/dotplot.R $analysis_dir/$chr"_"$chr".rdotplot" $analysis_dir/$chr"_"$chr".jpeg"

# Visually check the dotplot. If the dot plot looks good, proceed to STEP 03

################################################################################################
# STEP 03: Format the rdotplot output from lastZ into outputs to be used for python scripts.
# See the bash script "format_rdotplot.R" for more details
################################################################################################

$scripts_dir/format_rdotplot.sh $analysis_dir/$chr"_"$chr".rdotplot" $analysis_dir chrY

################################################################################################
# STEP 04: Generate Xkb non-overlapping windows (choose here 10kb)
################################################################################################
bedtools makewindows -g $analysis_dir/$chr"_length.txt" -w $xKb_windowSize > $analysis_dir/$chr"_"$xKb_windowSize"_nonOverlapping.txt"

###############################################################################################
# STEP 05: Run python script to obtain the masked regions
###############################################################################################
# There are 7 arguments for the Python script
arg_1=$analysis_dir/$chr"_"$xKb_windowSize"_nonOverlapping.txt"
arg_2=$analysis_dir/$chr"_target_clean_joined.txt"
arg_3=$analysis_dir/$chr"_query_clean_joined.txt"
arg_4=50000
arg_5=$analysis_dir/"target_out_bufferRegion"$arg_4".txt"
arg_6=$analysis_dir/"query_out_bufferRegion"$arg_4".txt"
arg_7=$analysis_dir/"window_query_out_bufferRegion"$arg_4".txt"

python $scripts_dir/generate_masks.py $arg_1 $arg_2 $arg_3 $arg_4 $arg_5 $arg_6 $arg_7

###############################################################################################
# STEP 07: Plotting for visualization
###############################################################################################
# A. Generate the rdotplot format for the "masked" regions
paste $arg_5 $arg_6 > $analysis_dir/"target_query_out_bufferRegion"$arg_4".rdotplot"
Rscript $scripts_dir/dotplot.R $analysis_dir/"target_query_out_bufferRegion"$arg_4".rdotplot" $analysis_dir/$chr"_"$chr"_masked_regions.jpeg"

# B. Overlaying
Rscript $scripts_dir/dotplot_overlap.R $analysis_dir/$chr"_"$chr".rdotplot" $analysis_dir/"target_query_out_bufferRegion"$arg_4".rdotplot" $analysis_dir/$chr"_"$chr"overlay_with_maskedRegions.jpeg"












