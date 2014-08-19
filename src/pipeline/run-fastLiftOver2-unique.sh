#!/bin/sh

#./run-fastLiftOver2-unique.sh [xToX.IndexFile] [x.meth] [x.meth.liftedToX] [x.meth.unmappedToX] [METHPIPE_DIR]

export LC_ALL=C

IndexFile="$1" 
FromFile="$2" 
ToFile="$3" 
Unmapped="$4" 
METHPIPE_DIR="$5" 

${METHPIPE_DIR}/fastLiftOver2 -i ${IndexFile} -f ${FromFile} -t ${ToFile} -v -u ${Unmapped}

tmpfile=$(mktemp)

####unique####
#Keep only 1-1 mapped positions

sort -k1,1 -k2,2n $ToFile \
|awk '{OFS="\t";}
NR == 1 {chr = $1; pos = $2; strand = $3; seq = $4; t = $6; m = $5 * $6; dup=0}
NR > 1 {
if ($1 == chr && $2 == pos && $3 == strand && $4 == seq)
{
t += $6; m += $5 * $6; dup =1;
}
else
{
if (t == 0) meth = 0.0; else meth = m / t;
if (dup == 0) print chr,pos,strand,seq,meth,t;
chr = $1; pos = $2; strand = $3; seq = $4; t = $6; m = $5 * $6; dup = 0;
}
}
END {
if (t == 0) meth = 0.0; else meth = m / t;
if (dup == 0) print chr,pos,strand,seq,meth,t;
}' > ${tmpfile} && mv ${tmpfile} ${ToFile}



