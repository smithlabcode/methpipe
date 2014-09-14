#!/bin/sh  
#./run-fastLiftOver2.sh -i indexfile -f fromfile -t tofile -M MethpipBinDir [-s -u]


usage()
{
cat << EOF
usage: $0 options
This script run the test1 or test2 over a machine.
OPTIONS:
   -h      Show this message
   -i      Index file
   -f      File to lift
   -t      Lifed file
   -M      Path to the methpipe root directory
   -s      Report sites on single strand (+ strand)
   -u      Report uniquely mapped sites (If not specified, multiple-to-1 mapping are collapsed)
EOF
}

IndexFile=
FromFile=
ToFile= 
METHPIPE_DIR= 
SINGLESTRAND= 
UNIQUE= 

while getopts "hi:f:t:M:su" OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         i)
             IndexFile=$OPTARG
             ;;
         f)
             FromFile=$OPTARG
             ;;
         t)
             ToFile=$OPTARG
             ;;
	 M)  
             METHPIPE_DIR=${OPTARG}
             ;;
         s)
             SINGLESTRAND=1
             ;;
	 u)
	     UNIQUE=1
             ;;
         ?)
             usage
             exit
             ;;
     esac
done

export LC_ALL=C

tmpfile=$(mktemp)
${METHPIPE_DIR}/fastLiftOver2 -i ${IndexFile} -f ${FromFile} -t ${ToFile} -v -u ${tmpfile}
rm ${tmpfile}

tmpfile=$(mktemp)
if [ "$SINGLESTRAND" ];
then
  awk '{OFS="\t"; 
  if($3=="+") print;
  if($3=="-") print $1,$2-1,"+",$4,$5,$6}' < $ToFile > ${tmpfile} && mv ${tmpfile} ${ToFile};
fi


if [ "$UNIQUE" ] ; 
then
  tmpfile=$(mktemp)
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
  }' > ${tmpfile} && mv ${tmpfile} ${ToFile};
else
  sort -k1,1 -k2,2n $ToFile \
  |awk '
  NR == 1 {chr = $1; pos = $2; strand = $3; seq = $4; t = $6; m = $5 * $6;}
  NR > 1 {
  if ($1 == chr && $2 == pos && $3 == strand && $4 == seq)
  {
  t += $6; m += $5 * $6;
  }
  else
  {
  if (t == 0) meth = 0.0; else meth = m / t;
  print chr,pos,strand,seq,meth,t;
  chr = $1; pos = $2; strand = $3; seq = $4; t = $6; m = $5 * $6;
  }
  }
  END {
  if (t == 0) meth = 0.0; else meth = m / t;
  print chr,pos,strand,seq,meth,t;
  }' > ${tmpfile} && mv ${tmpfile} ${ToFile};
fi



