#!/bin/bash

RUN=${1?"insert run number as input"} 
DAQERMFOLDER=${2?"insert path to folder containing xml files"}
DAQPROCESSEDFOLDER=${3?"insert path to folder containing output processed files, example: ./preprocess_entire_run.sh 2 /run/media/anastasia/HDD/DAQERM/ ../DAQprocessed"}

echo "Performing complete analysis on all the following files"
for filename in ${DAQERMFOLDER}/run${RUN}meas*.xml; do
	echo "$filename"
	sleep 0.25
done

for filename in ${DAQERMFOLDER}/run${RUN}meas*.xml; do
	./xmltoTTreeRDF ${DAQPROCESSEDFOLDER} "$filename"
done

for filename in ${DAQPROCESSEDFOLDER}/run${RUN}meas*_amulet.root; do
	./logicSignalAnalysis "$filename"
done

./lifetime ${DAQPROCESSEDFOLDER} ${DAQPROCESSEDFOLDER} ${RUN}
