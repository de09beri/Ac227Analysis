#!/bin/bash

rm filePath.txt

printf -v s "%03d" $1
RUN=$2
RUNEND=$3

echo "===================================================================="
echo "Going to combine trees of series ${s}, runs ${RUN} through ${RUNEND}"

echo "===================================================================="
echo "ADDING FILE PATHS TO TEXT FILE"

SERIESDIR=$(find $P2X_ANALYZED/180316_Background -type d -name "series${s}" -exec echo {} \;)

echo ${s} >> filePath.txt

while [ $RUN -le $RUNEND ]; do

	printf -v r "%05d" $RUN

	RUNDIR=$(find $SERIESDIR -type d -name "s${s}_f${r}_ts*" -exec echo {} \;)
	echo $RUNDIR

	if [ ${RUNDIR} ]; then
	
		FILE=$(find $RUNDIR -type f -name "ADAnalyzer.root" -exec echo {} \;) 
		echo ${FILE} >> filePath.txt	

	fi

	let RUN=RUN+1
done

echo "===================================================================="
echo "COMBINING FILES INTO ONE ROOT TREE"

root -l -b  <<EOF
.L CombineTrees.C+
CombineTrees()
.q
EOF





