#!/bin/bash
##### TODO: Make lastZ command more flexible
. /u/local/Modules/default/init/modules.sh
module load R
module load python
module load bedtools

cwd=$(pwd)

### Ask the user what build (hg19 vs hg38) and what chromosome

echo -n "Enter the build, hg19 or hg38, and press [ENTER]: "
read build

echo -n "Enter the chromosome you want to create a mask of and press [ENTER]: "
read chr

echo -n "Enter the type of lastZ you are running (ie default) and press [ENTER] "
read lastZ

echo -n "Enter the size of the nonOverlapping window and press [ENTER]: "
read xKb_windowSize

### Create a directory to dump all output

mkdir ${cwd}/${build}_${chr}_lastZ${lastZ}_outDir

lastz_dir=$cwd/lastz-distrib/bin #This is the directory where to call lastZ
analysis_dir=${cwd}/${build}_${chr}_lastZ${lastZ}_outDir #This is the directory where all the input and output files are stored. 
scripts_dir=$cwd/scripts #This is the directory where all of the scripts are stored

################################################################################################
# STEP 01: To run lastZ for self-self alignment
################################################################################################

$lastz_dir/lastz ${cwd}/fastaData/ucsc_${build}/$chr".fa" --self --notransition --ambiguous=iupac --nogapped --nomirror --step=20 --notrivial --format=rdotplot > $analysis_dir/$chr"_"$chr".rdotplot"

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
arg_8=$analysis_dir/"toMaskRegions.bed" ### masked regions in bed format 
arg_9=$analysis_dir/"identityMapped.bed" ### regions that are on the identity line in bed format 
#python $scripts_dir/generate_masks.py $arg_1 $arg_2 $arg_3 $arg_4 $arg_5 $arg_6 $arg_7 $arg_8 $arg_9
python $scripts_dir/generate_masks.py $arg_1 $arg_2 $arg_3 $arg_4 $arg_5 $arg_6 $arg_7

###############################################################################################
# STEP 07: Plotting for visualization
###############################################################################################
# A. Generate the rdotplot format for the "masked" regions
paste $arg_5 $arg_6 > $analysis_dir/"target_query_out_bufferRegion"$arg_4".rdotplot"
Rscript $scripts_dir/dotplot.R $analysis_dir/"target_query_out_bufferRegion"$arg_4".rdotplot" $analysis_dir/$chr"_"$chr"_masked_regions.jpeg"

# B. Overlaying
Rscript $scripts_dir/dotplot_overlap.R $analysis_dir/$chr"_"$chr".rdotplot" $analysis_dir/"target_query_out_bufferRegion"$arg_4".rdotplot" $analysis_dir/$chr"_"$chr"overlay_with_maskedRegions.jpeg"

##############################################################################################
# STEP 08: Compute some statistics, ie. how many base pairs are uniquely mapped vs multimapped
##############################################################################################

# A. Sort the 2 bed files
sort -n -k 2,2 $analysis_dir/"toMaskRegions.bed" > $analysis_dir/"toMaskRegions_sorted.bed"
sort -n -k 2,2 $analysis_dir/"identityMapped.bed" > $analysis_dir/"identityMapped_sorted.bed"

# B. Use bedtools merge

bedtools merge -i $analysis_dir/"toMaskRegions_sorted.bed" > $analysis_dir/"toMaskRegions_sorted_merged.bed"
bedtools merge -i $analysis_dir/"identityMapped_sorted.bed" > $analysis_dir/"identityMapped_sorted_merged.bed"

# C. Compute the number of base pair within each region
awk '{ $4 = $3 - $2} 1' $analysis_dir/"toMaskRegions_sorted_merged.bed" > $analysis_dir/"toMaskRegions_sorted_merged_diff.bed"

awk '{ $4 = $3 - $2} 1' $analysis_dir/"identityMapped_sorted_merged.bed" > $analysis_dir/"identityMapped_sorted_merged_diff.bed"

toMaskSum=`awk '{ sum += $4 } END {print sum} ' $analysis_dir/"toMaskRegions_sorted_merged_diff.bed"`
identityMappedSum=`awk '{ sum += $4 } END {print sum} ' $analysis_dir/"identityMapped_sorted_merged_diff.bed"`

echo 'Total number of base pairs that are masked is' $toMaskSum
echo 'Total number of base pairs that are on the identity line is' $identityMappedSum








