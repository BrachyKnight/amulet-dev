#!/bin/bash

RUN=${1?"insert run number as input"} 
PREPROCESSFOLDER=${2?"insert path to folder containing root files, example: ./logic_signal_analysis_entire_run.sh 1 ../DAQpreprocessed/"}
echo "Performing complete analysis on all the following files"
for filename in ${PREPROCESSFOLDER}/run${RUN}meas*_amulet.root; do
	echo "$filename"
	sleep 0.25
done

for filename in  ${PREPROCESSFOLDER}/run${RUN}meas*_amulet.root; do
	./../digitizer_scripts/logicSignalAnalysis "$filename"
done


