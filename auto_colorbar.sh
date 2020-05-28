#! /bin/bash



# Copyright (C) 2013 by Arizona State University and Mark Robinson
# All rights reserved.
#
# $Author: mhenriksen $
#   $Date: 2020-05-27 17:00:55 -0700 (Wed, 27 May 2020) $
#    $Rev: 20778 $

#FOLLOWING ARE DEBUG INPUTS
#filename=cb_rumker.prn
#filename=inad_custom_cb7.lut
#globalMin=-2947
#globalMax=-1294

DEBUG=false #set to false/true to toggle temporarily generated files for debugging.

filename=$1
globalMin=$2
globalMax=$3

#FOLLOWING BLOCK OUTPUTS USAGE INFORMATION.
if [[ "$#" -eq 0 ]]; 
then 
     echo "
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
USAGE:
./auto_colobar.sh <file name> <global min elevationas integer> <global max elevation as integer> <-g>

NOTES:
<file name> can be either *.lut or *.prn.
<-g> is flag for doing continuous gradients between spacing for .lut files, this flag must be placed at end.

Example: ./auto_colorbar.sh cb_rumker.prn -2947 -1294
Example: ./auto_colorbar.sh inad_custom_cb7.lut -2947 -1294
Example: ./auto_colorbar.sh colorbrew_eq_cont.lut -2947 -1294 -g

OUTPUT:
The program output is the legend.png file.

TEMPORARY FILES:
Files generated and cleaned up during script are: 

elevationsList_debug.txt: Reads evelation segmentations from .lut or .prn file. 
colorsList_debug.txt: Tracks colors to be used as either fill or benchmarks for gradient setting.
temp_debug.lut: A version of the input .lut file without comments.
test_debug.prn: A file holding the rgb values for each section or section benchmark.

DEBUG:
Set DEBUG=true to preserve temporary files generated for data handling.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"
     exit 2
fi

checkFlag=false
#FOLLOWING BLOCK CHECKS NUMBER OF INPUTS FOR ERRORS
if [[ "$#" -eq 3 ]];then
    check1=true
elif [[ "$#" -eq 4 && $4 == "-g" ]]; then
    checkFlag=true
else
    echo "
ERROR: Incorrect number of inputs or flag, check usage by running ./auto_colorbar.sh with no inputs.
"
    exit 2
fi

#FOLLOWING BLOCK CHECKS INPUT FILETYPE FOR ERROR CHECKING
if [[ $filename == *.prn || $filename == *.lut ]]; then
    check2=true
else
    echo "
ERROR: Input file incorrect format
"
    exit 2
fi

#FOLLOWING BLOCK COMPARES MAX/MIN VALUES FOR ERROR CHECKING
if (( $(echo "$globalMax > $globalMin" |bc -l) )); then
    check3=true
else
    echo "
ERROR: Second input (global minimum) must be less than third input (gloal maximum).
"
    exit 2
fi

#SCRIPT CONTINUES IF ALL TESTS PASSED.

printf -v globalMin %.3f "$globalMin"  #Converts variable for easier mathematical handling
printf -v globalMax %.3f "$globalMax"  #Converts variable for easier mathematical handling

spread=$(echo  $globalMax - $globalMin | bc -l)  #spread used for calculating divisions later

convert -size 500x2550 xc:white legend.png #Overwrites any existing legend.png as a new blank, white legend.png of specified size

[ -e test_debug.prn ] && rm test_debug.prn #cleans old temp test_debug.prn file if it exists
[ -e temp_debug.lut ] && rm temp_debug.lut #cleans temp temp_debug.lut file
[ -e elevationsList_debug.txt ] && rm elevationsList_debug.txt #cleans temp colorsList_debug.txt file
[ -e colorsList_debug.txt ] && rm colorsList_debug.txt #cleans temp colorsList_debug.txt file

#BLOCK CREATES TEMPORARY FILE elevationsList_debug.txt TO TRACK DESIRED DEVISIONS BETWEEN MIN/MAX
#DIFFERENT EXTRACTION METHOD DEPENDING ON FILETYPE
if [[ $filename == *".lut"* ]];then
    sed 's|\(.*\)//.*|\1|' $filename > temp_debug.lut #Removes all comments from input .lut file.
    if [ "$checkFlag" == true ];then # Block for -g flag .lut option
	cat temp_debug.lut | while read -r percent1 r g b; do
	    #Loop reads through temp_debug.lut extracting colors and reading percent distributions into
	    #test_debug.prn as elevation segmentations,
	    if [ "$percent1" != "nv" ]; then
		deci1=`echo "${percent1::-1}/100.0" |bc -l`
		divider1=`echo "($deci1*$spread)+$globalMin"|bc -l`
		first=`printf %.0f $(echo $divider1)`
		echo "$first $r $g $b" >>test_debug.prn
		echo $first >>elevationsList_debug.txt
	    fi
	done
    else #Non -g flag .lut option
	cat temp_debug.lut | while read -r percent1 r g b; do
	    #Loop reads through temp_debug.lut extracting colors and reading percent distributions into
	    #test_debug.prn as elevation segmentations. It reads two lines at a time and skips 'nv' line.
	    if [ "$percent1" != "nv" ]; then
		deci1=`echo "${percent1::-1}/100.0" |bc -l`
		divider1=`echo "($deci1*$spread)+$globalMin"|bc -l`
		first=`printf %.0f $(echo $divider1)`
		read -r percent2 r g b 
		deci2=`echo "${percent2::-1}/100.0" |bc -l`
		divider2=`echo "($deci2*$spread)+$globalMin"|bc -l`
		second=`printf %.0f $(echo $divider2)`
		echo "$first $second $r $g $b" >>test_debug.prn
		echo $first >>elevationsList_debug.txt
		echo $second >>elevationsList_debug.txt
	    fi
	done
    fi
elif [[ $filename == *".prn"* ]];then #Following extracts elevations from .prn filetype.
    cat $filename | while read -r first second r g b; do 
	echo "$first $second $r $g $b" >>test_debug.prn
	echo $first >>elevationsList_debug.txt
	echo $second >>elevationsList_debug.txt
    done
fi

readarray -t arr < elevationsList_debug.txt #reads elevationsList_debug.txt into array variable
#[ -e temp_debug.lut ] && rm temp_debug.lut #cleans temp elevationsList_debug.txt file
#[ -e elevationsList_debug.txt ] && rm elevationsList_debug.txt #cleans temp elevationsList_debug.txt file
eval a=($(printf "%q\n" "${arr[@]}" | sort -un))

#INITIALIZE COORDINATES FOR BOUNDARIES OF COLORBAR
bottomEnd=$(( 1250*2 ))
topEnd=$(( 80*2 ))
range=$((bottomEnd-topEnd))
sizes=()

#########################################################################
#Algorithm to adjust the maximum bin size.
# If a percent is larger than $tolerance, this will shrink it to $tolerance
# and restribute the overage amount equally to the other bins.
# Process is iterative so that if a bin overflows $tolerance due to readjustment
# it will fix the new overflow.

tolerance=.30001 #USE THIS TO SET DESIRED MAX BIN SIZE, SUGGESTED USE IS 5 DECIMAL ACCURACY
#Example: Set tolerance=0.30001 for 30% bin size
#Example: Set tolerance=0.20001 for 20% bin size

max_value=1
arraylength="${#a[@]}"
percent_array=()

#Calculate the initial percent coverage of each bin
for (( i=1; i<${arraylength}; i++ ));
do
    percent=`echo "( ${a[$i]} - ${a[i-1]} )/( $globalMax - $globalMin )" | bc -l`
    percent_array[$i-1]=$percent
done

#While the max size in the array is above the tolerance, the following
#loop tries to redistribute the bin sizes to meet the tolerance.
while (( $(echo "$max_value > $tolerance" | bc -l) )); do
    over_total=0 # tracks total overage of 30% in all bins.
    index_array=() # tracks index of bins with > $tolerance coverage.
    under_count=0 # count number of bins with <= $tolerance.
    
    #For loop builds index_array as a zero array.
    for (( i=0; i<${arraylength}-1; i++ ));do
	index_array+=(0)
    done

    #Following checks for bins which violat the tolerance
    for (( i=1; i<${arraylength}; i++ ));
    do
	percent_check=${percent_array[$i-1]}
	if (( $(echo "$percent_check > $tolerance" | bc -l) ));then
	    over=`echo "($percent_check-$tolerance)" | bc -l` #overflow of current bin
	    percent_check=$tolerance #force bin size to $tolerance
	    over_total=`echo "($over_total+$over)" | bc -l`
	    index_array[$i-1]=1
	else
	    index_array[$i-1]=0
	fi
	percent=`printf %.10f $(echo $percent_check)`
	percent_array[$i-1]=$percent
    done
    
    #count number of non-overflow values in index_array
    for x in "${index_array[@]}";do
	if [ $x -eq 0 ];then
	    under_count=$((under_count+1))
	fi
    done

    #calculate amount to be added to each bin and add it
    each=`echo "( $over_total )/( $under_count )" | bc -l`
    for (( i=0; i<${arraylength}-1; i++ ));do
	if [ ${index_array[$i]} -eq 0 ];then
	    new_percent=`echo "( ${percent_array[$i]} + $each )" | bc -l`
	    percent_array[$i]=$new_percent
	fi
    done

    max_value=$tolerance #Setting this here will let loop finish if no offending values are found
    for value in "${percent_array[@]}";do
	if (( $(echo "$value > $tolerance" | bc -l) )); then
	    max_value=$value
	fi
    done
done

####################################################################################

#The following for-loop generates a list of the pixel sizes of each segment in the colorbar.
for (( i=1; i<${arraylength}; i++ ));
do
    b=`echo "${percent_array[$i-1]}9*$range" | bc -l`
    c=`printf %.0f $(echo $b)`
    sizes+=($c)
done

#The following for-loop builds the array coords to keep track of the pixel positions
#for each segment of the colorbar.
coords=($bottomEnd)
track=$bottomEnd
for x in "${sizes[@]}";do
    track=$(( track - x ))
    coords+=($track)
done
coords[-1]=$topEnd

#The following loop reads the rgb colors from test_debug.prn for each section
#and places them into the array 'arr'. For the -g flag test_debug.prn has one
#less column, so this if-statment accounts for that.
if [ "$checkFlag" == true ] ; then
    cat test_debug.prn | while read elevation r g b; do
	printf -v color "#%02X%02X%02X" $r $g $b
	echo "$color"
    done >> colorsList_debug.txt
else
    cat test_debug.prn | while read min max r g b; do
	printf -v color "#%02X%02X%02X" $r $g $b
	echo "$color"
    done >> colorsList_debug.txt
fi
readarray arr < colorsList_debug.txt


legendText=50
i=0
arrLength=`echo ${#arr[*]}`
if [ "$checkFlag" == true ];then # Block for -g flag option
    while [ $i -le $(( $arrLength-2 )) ]; do
	yShiftLine=$(( ${coords[i]} + 0 )) #Shift for drawing horizontal line on colorbar.
	ySize=$(( ${coords[i]} - ${coords[i+1]} )) #Calculates pixel size of each bin.
	sizeName="240x$ySize" #Makes string to pass into convert gradient command of -g flag.
	trueY=$(( ${coords[i]} - $ySize )) #The real pixel location where the gradient bin will be placed.
	yShiftText=$(( ${coords[i]} + 8 )) #Shift applied to drawing colorbar text.

	#Following line draws gradient bin
	convert legend.png \( -size "$sizeName" gradient:"${arr[$i+1]}-${arr[$i]}" \) -geometry +40+$trueY -compose over -composite legend.png
	#Following line writes elevation label
	convert legend.png -pointsize "$legendText" -fill black -draw "text 320,'$yShiftText' '${a[$i]}' " legend.png
	#Following line draws horizontal line when applicable.
	if [[ $i -gt 0 ]];then
	    convert legend.png -stroke black -strokewidth 5 -draw "line 280,'$yShiftLine' 310,'$yShiftLine'" legend.png 
	fi
	((i++))
    done
else #Block for non -g flag option
    for x in "${arr[@]}";do
	yShiftLine=$(( ${coords[i]} + 0 )) #Shift for drawing horizontal line on colorbar.
	yShiftText=$(( ${coords[i]} + 8 )) #Shift applied to drawing colorbar text.

	#Following line draws filled bin with solid color.
	convert legend.png -fill "$x" -draw "rectangle 40,'${coords[$i]}' 280,'${coords[$i+1]}'" legend.png
	#Following line writes elevation label
	convert legend.png -pointsize "$legendText" -fill black -draw "text 320,'$yShiftText' '${a[$i]}' " legend.png
	#Following line draws horizontal line when applicable.
	if [[ $i -gt 0 ]];then
	    convert legend.png -stroke black -strokewidth 5 -draw "line 280,'$yShiftLine' 310,'$yShiftLine'" legend.png 
	fi
	((i++))
    done
fi

yShift=$(( $topEnd + 8  ))
convert legend.png -pointsize "$legendText" -fill black -draw "text 320,'$yShift' '${a[$i]}' " legend.png #Writes final elevation label
convert legend.png -stroke black -strokewidth 5 -draw "line 38,'$bottomEnd' 310,'$bottomEnd'" legend.png #Draws bottom horizontal line
convert legend.png -stroke black -strokewidth 5 -draw "line 40,'$bottomEnd' 40,'$topEnd'" legend.png #Draws left colorbar border
convert legend.png -stroke black -strokewidth 5 -draw "line 280,'$bottomEnd' 280,'$topEnd'" legend.png #Draws right colorbar border
convert legend.png -stroke black -strokewidth 5 -draw "line 38,'$topEnd' 310,'$topEnd'" legend.png #Draws top colorbar border
convert legend.png -pointsize 70 -fill black -draw "text 33,80 'Elevation (m)'" legend.png #Writes colorbar title

#Following block will delete temporary files generated during script if DEBUG is set to false.
if [ "$DEBUG" == false ] ;then
    [ -e test_debug.prn ] && rm test_debug.prn #cleans old temp test_debug.prn file if it exists
    [ -e temp_debug.lut ] && rm temp_debug.lut #cleans temp temp_debug.lut file
    [ -e elevationsList_debug.txt ] && rm elevationsList_debug.txt #cleans temp elevationsList_debug.txt file
    [ -e colorsList_debug.txt ] && rm colorsList_debug.txt #cleans temp colorsList_debug.txt file
fi
