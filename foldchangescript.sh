#!/bin/bash

#make new directory to store fold change data
mkdir fold_change_data

#define variables for all the mean counts files(groups), the total number of groups and a counter
GROUP_FILES=(mean_counts_data/*.txt)
TOTAL_GROUPS=$(ls mean_counts_data | wc -l)
COUNT=0

#for loop that loops through every group file and...
for GROUP in ${GROUP_FILES[@]}
do
#sets index to 0
    i=0
#while the index value is less than the total number of group files...
    while [ $i -lt $((TOTAL_GROUPS)) ]
    do
#if the counter does not equal the index (stops the group file being compared to itself) 
        if [ $COUNT -ne $i ]
        then
#assign the group name of the file in the current interation of the for loop to FILE1NAME and the group...
#name of the the group file at index position $i to FILENAME2
            FILE1NAME=$(basename $GROUP _meancounts.txt)
            FILE2NAME=$(basename ${GROUP_FILES[$i]} _meancounts.txt)

#awk command that takes the third field from each file if the value is not the header and the value from the second file does not equal 0 (dividing by zero would cause an error) and prints the fold change from file1 to file2.
#if the value in the second file does equal 0, print NA (fold change caluclation not possible).
#output is stored in a temporary file.
            echo "Creating fold change data for $FILE1NAME to $FILE2NAME"
	    awk 'BEGIN{FS="\t"; OFS="\t"} ARGIND==1 && NR!=1{a[$1,$2]=$3;next} ARGIND==2 {if($3!="Mean_counts" && a[$1,$2]!=0) print $1,$2,$3/a[$1,$2]; else if($3!="Mean_counts" && a[$1,$2]==0) print $1, $2, "NA"}' ${GROUP} ${GROUP_FILES[$i]} | sort -t$'\t' -k3,3nr > tmp.txt

#echo the header names and add the values in the temporary file to the end, storing the output in a new file named by the group name of the first file 'to' the group name of the second file.            
	    echo -e "Gene_name\tGene_description\tFold_change" | cat - tmp.txt > fold_change_data/${FILE1NAME}_to_${FILE2NAME}.txt

#remove the temporary file	    
	    rm tmp.txt
        fi

#add 1 to the index variable
        i=$(( i + 1 ))
    done

#add 1 to the counter
    COUNT=$((COUNT+1))
done
