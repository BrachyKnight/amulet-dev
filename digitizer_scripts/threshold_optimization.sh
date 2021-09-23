#!/bin/bash

RUN=${1?"insert run number as input"} 
PREPROCESSFOLDER=${2?"insert path to folder containing root files, example: ./logic_signal_analysis_entire_run.sh 1 ../DAQpreprocessed/"}
PROCESSFOLDER=${3?"insert path to folder containing root files, example: ./logic_signal_analysis_entire_run.sh 1 ../DAQpreprocessed/ ../DAQprocessed/"}

for th in $(seq 1 2 50); do #FIRST INCREMENT LAST
	MINTH=$((th))
	MAXTH=$((100-th))

	echo "MinTh = " ${MINTH};
	echo "MaxTh = " ${MAXTH};

	echo "Performing complete analysis on all the following files"
	for filename in ${PREPROCESSFOLDER}/run${RUN}meas*_amulet.root; do
		echo "$filename"
		sleep 0.25
	done

	for filename in  ${PREPROCESSFOLDER}/run${RUN}meas*_amulet.root; do
		./../digitizer_scripts/logicSignalAnalysis "$filename" ${MINTH} ${MAXTH};
	done

	DIRTHR=${PROCESSFOLDER}/min${MINTH}max${MAXTH}/
	if [[ ! -d "$DIRTHR" ]]
	then
	    echo "$DIRTHR does not exists on your filesystem."
	    mkdir ${DIRTHR}
	fi
	./lifetime ${DIRTHR} ${PREPROCESSFOLDER} ${RUN};
done
