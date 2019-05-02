#!/bin/bash
paraFile=$1
source $paraFile

## ==========================
## print the parameter values
## ==========================
echo "projDir=${projDir}"
echo "fastqFile=${fastqFile}"
echo "ref=${ref}"
echo "refrag=${refrag}"
echo "outDir=${outDir}"
echo "simuName=${simuName}"
echo "summaryFile=${summaryFile}"

echo "bwaDir=${bwaDir}"
echo "samtoolsDir=${samtoolsDir}"
echo "bedtoolsDir=${bedtoolsDir}"

echo "train=${train}"
echo "postProcess=${postProcess}"
echo "coreN=$coreN"
echo "mismatchN=${mismatchN}"
echo "gapN=${gapN}"
echo "mismatchP=${mismatchP}"
echo "gapP=${gapP}"
echo "chimericP=${chimericP}"
echo "simuN=${simuN}"
echo "readLen=${readLen}"
echo "resolution=${resolution}"
echo "lowerBound=${lowerBound}"
echo "refragU=${refragU}"
echo "cutsite=${cutsite}"

## =================
## Start of FreeHi-C
## =================
bash freeHiC.sh "$projDir" "$fastqFile" "$outDir" "$simuName" "$ref" "$bwaDir" "$samtoolsDir" "$bedtoolsDir" "$train" "$postProcess" "$coreN" "$mismatchN" "$gapN" "$mismatchP" "$gapP" "$chimericP" "$simuN" "$readLen" "$resolution" "$lowerBound" "$refragU" "$cutsite" "$refrag" "$summaryFile" >$paraFile.log

