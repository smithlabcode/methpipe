#!/bin/sh

#run liftOver program on methcount file with all cytosines with context and strand information
#run-liftOver-allC.sh [Directory of liftOver tool] [xToX.over.chain.gz] [x.meth] [FILENAME-xToXindex] [FILENAME-unlifted]

LiftOver_DIR="$1" 
Chain="$2" 
Source="$3" 
IndexFile="$4" 
Unlift="$5" 

tmpfile=$(mktemp)
awk '{OFS="\t"; print $1,$2,$2+1,$1":"$2":"$2+1":"$4":"$3, 0,$3}' < ${Source} > ${tmpfile} 
${LiftOver_DIR}/liftOver ${tmpfile} ${Chain} ${IndexFile} ${Unlift} && rm ${tmpfile} 

