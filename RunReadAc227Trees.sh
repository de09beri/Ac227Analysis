#!/bin/bash

rm $AD_AC227ANALYSIS_DATA_PLOTS/*

root -l -b <<EOF

.L $AD_AC227ANALYSIS_SCRIPTS/ReadAc227Trees.C+
const int numCells = 154;
const int numTrees = 590;
const int PLOTFLAG = 1;
double PSDCUTLOW = 0.18;
double PSDCUTHIGH = 0.35; 
PlotResults(numCells,numTrees,PLOTFLAG,PSDCUTLOW,PSDCUTHIGH)

.q

EOF


