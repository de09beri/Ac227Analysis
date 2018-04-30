#!/bin/bash

# This script will combine all Ac227 root trees created using the Ac227TreePlugin
# ./RunCombineTrees.sh <Data type> <Series> <First run> <Last run>
# Data type = 0, reactor on
# Data type = 1, reactor rampdown
# Data type = 2, reactor off

$AD_AC227ANALYSIS_SCRIPTS/RunCombineTrees.sh 0 9 1 33
$AD_AC227ANALYSIS_SCRIPTS/RunCombineTrees.sh 0 12 0 9
$AD_AC227ANALYSIS_SCRIPTS/RunCombineTrees.sh 0 13 1 10
$AD_AC227ANALYSIS_SCRIPTS/RunCombineTrees.sh 0 14 0 2
$AD_AC227ANALYSIS_SCRIPTS/RunCombineTrees.sh 0 15 0 20
$AD_AC227ANALYSIS_SCRIPTS/RunCombineTrees.sh 0 16 0 16
$AD_AC227ANALYSIS_SCRIPTS/RunCombineTrees.sh 0 17 1 5
$AD_AC227ANALYSIS_SCRIPTS/RunCombineTrees.sh 0 18 1 23
$AD_AC227ANALYSIS_SCRIPTS/RunCombineTrees.sh 0 19 0 24
$AD_AC227ANALYSIS_SCRIPTS/RunCombineTrees.sh 0 20 1 19
$AD_AC227ANALYSIS_SCRIPTS/RunCombineTrees.sh 0 21 0 4
$AD_AC227ANALYSIS_SCRIPTS/RunCombineTrees.sh 0 22 0 19
$AD_AC227ANALYSIS_SCRIPTS/RunCombineTrees.sh 0 23 0 112

$AD_AC227ANALYSIS_SCRIPTS/RunCombineTrees.sh 1 0 0 3

$AD_AC227ANALYSIS_SCRIPTS/RunCombineTrees.sh 2 0 0 22
$AD_AC227ANALYSIS_SCRIPTS/RunCombineTrees.sh 2 1 0 331

