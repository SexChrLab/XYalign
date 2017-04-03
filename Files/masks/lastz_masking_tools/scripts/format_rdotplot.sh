#!/bin/bash

#### Output from lastZ in rdotpot format (please see the documentation for lastZ for more details):
#chrY chrY
#13308   17343
#13375   17410
#NA      NA
#17578   17641
#17631   17694
#NA      NA

#### Want to format this rdotplot output such that the output will be 2 files. One file contains the coordinates for the target and one file contains the coordinates for the query. For consistency, I refer to the coordinates on the first column as the target, and the coordinates on the second column as query

lastZOutput=$1
dir=$2
chr=$3
### STEP 1: extract column 1 for target, extract column 2 for query
awk '{print$1}' $lastZOutput > $dir/$chr"_target.txt"
awk '{print$2}' $lastZOutput > $dir/$chr"_query.txt"

### STEP 2: Remove the first line and NA
grep -v chr $dir/$chr"_target.txt" | grep -v NA > $dir/$chr"_target_clean.txt"
grep -v chr $dir/$chr"_query.txt" | grep -v NA > $dir/$chr"_query_clean.txt"

### STEP 3: Join every two lines together so that the format isL: start	end
sed 'N;s/\n/\t/' $dir/$chr"_target_clean.txt" > $dir/$chr"_target_clean_joined.txt"
sed 'N;s/\n/\t/' $dir/$chr"_query_clean.txt" > $dir/$chr"_query_clean_joined.txt"
