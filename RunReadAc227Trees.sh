#!/bin/bash

rm $AD_AC227ANALYSIS_DATA_PLOTS/*

root -l -b <<EOF

.L $AD_AC227ANALYSIS_SCRIPTS/ReadAc227Trees.C+
const int numCells = 154;
const int numTrees = 590/10;
const int PLOTFLAG = 2;
PlotResults(numCells,numTrees,PLOTFLAG)

.q

EOF


