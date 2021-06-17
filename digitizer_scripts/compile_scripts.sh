#!/bin/bash
c++ -Wall -o xmltoTTreeRDF xmltoTTreeRDF.cpp  `root-config --cflags --glibs`
c++ -Wall -o logicSignalAnalysis logicSignalAnalysis.cpp `root-config --cflags --glibs`
c++ -Wall -o lifetime lifetime.cpp `root-config --cflags --glibs`

