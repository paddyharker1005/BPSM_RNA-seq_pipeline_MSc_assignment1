#!/bin/bash

#define the path of the raw sequence data as a variable
FASTQ_PATH=/localdisk/data/BPSM/AY21/fastq 

#make a directory to store the FASTQC results
mkdir fastqc_results

#for each sample file in the fastq directory with the extension .fq.gz perform a quality check using FASTQC and automatically extract the .zip output files   

fastqc -t 10 --extract ${FASTQ_PATH}/*.fq.gz -o fastqc_results

#concatenate all summary.txt files into one new temporary summary file
cat fastqc_results/*/summary.txt > fastqc_result_summaries.txt 

#wrangler using an awk command that stores all pass/warnings/fails for each sample and then prints them as individual lines
awk 'BEGIN{FS= "\t"; OFS= "\t"} {A[$3] = A[$3] FS $1} END{for(i in A) print i A[i]}' fastqc_result_summaries.txt > fastqc_wrangled.txt
#echos the ID and fastqc module headers and adds the values from the wranlged file to create a summary table 
echo -e "ID\tBasic Statistics\tPer base sequence quality\tPer sequence quality scores\tPer base sequence content\tPer sequence GC content\tPer base N content\tSequence Length Distribution\tSequence Duplication Levels\tOverrepresented sequences\tAdapter Content" | cat - fastqc_wrangled.txt > fastqc_results/fastqc_summary_table.txt
rm fastqc_wrangled.txt

#grep any fails and warnings and display the total counts to the user 
TOTAL_FAIL=$(grep -c "FAIL" fastqc_result_summaries.txt)
echo "${TOTAL_FAIL} FAILS detected:"

TOTAL_WARNING=$(grep -c "WARN" fastqc_result_summaries.txt)
echo "${TOTAL_WARNING} WARNINGS detected:"

#if the total number of fails or warnings is >1 ask if the user would like to see the details
#if the user replies y display the filenames and FASTQC modules of all warnings/fails
#if the use replies n then continue
#if the user replies with anything other than y/n, prompt them again
if [ $TOTAL_FAIL -gt 0 ] || [ $TOTAL_WARNING -gt 0 ] 
then
        while true; do
            read -p "View WARNINGS/FAILS? (y/n)" yn
            case $yn in
                [Yy]* ) echo "FAILS:" && cat fastqc_result_summaries.txt | grep "FAIL" | awk 'BEGIN {FS="\t"} ; {print $3, $2}'; echo "WARNINGS:" && cat fastqc_result_summaries.txt | grep "WARN" | awk 'BEGIN {FS="\t"}; {print $3, $2}'; break;;
                [Nn]* ) break;;
                * ) echo "Please answer either yes or no.";;
            esac
        done
fi

#remove the temporary summary file to tidy up
rm fastqc_result_summaries.txt
