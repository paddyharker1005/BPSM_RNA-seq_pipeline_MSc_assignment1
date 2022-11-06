#!/bin/bash

#inform user fastqc analysis is starting and execute the fastqc module of the pipeline
echo "Beginning fastqc quality check on raw data..."
./fastqcscript.sh &&
#prompt for user input based on fastqc summary - y to continue with alignment and quantification, n to exit pipeline
while true; do
    read -p "Do you wish to continue with alignment and quantification? (y/n)" yn
    case $yn in
        [Yy]* ) break;;
        [Nn]* ) exit;;
#if anything other than Y/y/N/n entered by the user, prompt again - keep looping until appropriate answer given
        * ) echo "Please answer either yes or no.";;
    esac
done &&
#inform user alignment and quantfication is starting and exute the alignment quanitification module of the pipeline
echo "Beginning alignment and quantification..."
./alignmentscript.sh &&
#inform user generation of mean counts data is starting and execute the mean counts module of the pipeline
echo "Generating mean counts data..."
./meancountsscript.sh &&
#inform user generation of fold change data is starting and execute the fold change module of the pipeline
echo "Generating fold change data..."
./foldchangescript.sh &&
#inform user that the pipeline analysis is complete and that all results are stored in sub-directories in the users current directory
echo "Analysis complete, results stored in subdirectories within current directory." 

