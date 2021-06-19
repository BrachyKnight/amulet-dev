#!/bin/bash

RUN=${1?"insert run number as input"} 
DAQERMFOLDER=${2?"insert path to folder containing xml files"}
echo "Performing complete analysis on all the following files"
for filename in ${DAQERMFOLDER}/run${RUN}meas*.xml; do
	echo "$filename"
	sleep 0.25
done

for filename in ${DAQERMFOLDER}/run${RUN}meas*.xml; do
	./../digitizer_scripts/xmltoTTreeRDF ../DAQpreprocessed/ "$filename"
done

for filename in ../DAQpreprocessed/run${RUN}meas*_amulet.root; do
	./../digitizer_scripts/logicSignalAnalysis ../DAQpreprocessed/ "$filename"
done

#hadd -f ../DAQpreprocessed/run${RUN}_allmeas_amulet.root ../DAQpreprocessed/run${RUN}meas*_amulet.root

