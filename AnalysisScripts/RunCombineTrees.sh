#!/bin/bash

rm $AD_AC227ANALYSIS_SCRIPTS/filePath.txt

DATATYPE=$1

printf -v s "%03d" $2
RUN=$3
RUNEND=$4

echo "===================================================================="
echo "Going to combine trees of series ${s}, runs ${RUN} through ${RUNEND}"

echo "===================================================================="
echo "ADDING FILE PATHS TO TEXT FILE"

if [ $DATATYPE -eq 0 ]
then
	SERIESDIR=$(find $P2X_ANALYZED/WetCommissioning -type d -name "series${s}" -exec echo {} \;)
fi

if [ $DATATYPE -eq 1 ]
then
	SERIESDIR=$(find $P2X_ANALYZED/180316_Rampdown -type d -name "series${s}" -exec echo {} \;)
fi

if [ $DATATYPE -eq 2 ]
then
	SERIESDIR=$(find $P2X_ANALYZED/180316_Background -type d -name "series${s}" -exec echo {} \;)
fi

echo ${s} >> $AD_AC227ANALYSIS_SCRIPTS/filePath.txt

while [ $RUN -le $RUNEND ]; do

	printf -v r "%05d" $RUN

	RUNDIR=$(find $SERIESDIR -type d -name "s${s}_f${r}_ts*" -exec echo {} \;)
	echo $RUNDIR

	if [ ${RUNDIR} ]; then
	
		FILE=$(find $RUNDIR -type f -name "ADAnalyzer.root" -exec echo {} \;) 
		echo ${FILE} >> $AD_AC227ANALYSIS_SCRIPTS/filePath.txt	

	fi

	let RUN=RUN+1
done

echo "===================================================================="
echo "COMBINING FILES INTO ONE ROOT TREE"

root -l -b  <<EOF
int DataType = $DATATYPE 
.L $AD_AC227ANALYSIS_SCRIPTS/CombineTrees.C+
CombineTrees(DataType)
.q
EOF





