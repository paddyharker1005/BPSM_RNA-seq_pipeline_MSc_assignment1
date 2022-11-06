#!/bin/bash

#define fastq file information path as a variable
FASTQ_PATH=/localdisk/data/BPSM/AY21/fastq

#make a new directory to store mean counts data
mkdir mean_counts_data

#define a count variable and set it to zero
COUNT=0

#take all lines except the headers of the 100k.fqfiles sample information file and sort by all experimental variable fields before sorting by the replicate number field numerically
#pipe the output to a while loop that reads each line of the sorted information and... 
tail -n +2 $FASTQ_PATH/100k.fqfiles | sort -k4,4n -k2,2 -k5,5 -k3,3n | while read line;
do
#adds 1 to the count
    COUNT=$((COUNT+1))

#if the count is 1 then treat as the first replicate and assign the sample name to the variable REPLICATE1
    if test $COUNT -eq 1
    then
        REPLICATE1=$(echo $line | awk '{print $1}')

#if the count is 2 then treat as the second replicate and assign the sample name to the variable REPLICATE2
    elif test $COUNT -eq 2
    then 
        REPLICATE2=$(echo $line | awk '{print $1}')

#if the count is 3 then treat as the third replicate and assign the sample name to the variable REPLICATE3
#take the other experimental variable details for this sample and assign to the variable GROUP 
    else
        REPLICATE3=$(echo $line | awk '{print $1}')
        GROUP=$(echo $line | awk 'BEGIN {OFS = "_"} {print $2, $4 "h", $5}')

#use the replicate variables to input the replicate sample counts files into an awk command that extracts the the sixth (counts) column from each of the replicates counts files and calculates the mean counts.
#the gene name and description and the mean counts are stored in a tempory file.
	echo "Combining ${GROUP} replicates into mean counts data..."
        awk 'BEGIN{FS="\t"; OFS="\t"} ARGIND==1{f1[$4,$5]=$6;next} ARGIND==2{f2[$4,$5]=$6;next} ARGIND==3{sum = $6 + f1[$4,$5] + f2[$4,$5]; print $4,$5, sum/3}' counts_data/100k.${REPLICATE1}_counts.txt counts_data/100k.${REPLICATE2}_counts.txt counts_data/100k.${REPLICATE3}_counts.txt > tmp.txt
        
#echo the column headers and add the values from the temporary file to the end, saving the output in a new mean counts file with the sample group included in the file name.
	echo -e "Gene_name\tGene_description\tMean_counts" | cat - tmp.txt > mean_counts_data/${GROUP}_meancounts.txt

#remove temporary file and reset the counter to 0        
	rm tmp.txt
        COUNT=0
    fi
done

