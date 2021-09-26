#!/bin/bash
c++ -Wall -o xmltoTTreeRDF xmltoTTreeRDF.cpp  `root-config --cflags --glibs` -O3
c++ -Wall -o WVFdisplayer WVFdisplayer.cpp `root-config --cflags --glibs` -O3
c++ -Wall -o logicSignalAnalysis logicSignalAnalysis.cpp `root-config --cflags --glibs` -O3
c++ -Wall -o lifetime lifetime.cpp `root-config --cflags --glibs` -O3
c++ -Wall -o amulet amulet.cpp `root-config --cflags --glibs` -O3
