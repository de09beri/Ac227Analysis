#!/bin/bash

rm $AD_AC227ANALYSIS_DATA_PLOTS/*

root -l -b <<EOF

.L ReadAc227Trees.C+
string treeDir = Form("%s",gSystem->Getenv("P50X_AC227ANALYSIS_TREES"));
const int numCells = 154;
const int numTrees = 645/5;
const int PLOTFLAG = 0;
PlotResults(treeDir,numCells,numTrees,PLOTFLAG)

.q

EOF


