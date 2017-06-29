#!/bin/bash

# This script will automatically run prediction softwares.
# Some require installation and are available at following links :
# Pyfasta (Python 2.7)
# ANCHOR : http://anchor.enzim.hu/Downloads.php
# SignalP : http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?signalp

# Set path variables yourself
DATA=~/Documents/Slavoff_SPEP # Path to your fasta file
ANCHOR=/usr/local/bin/anchor # Path to Anchor software
SIGNALP=~/Downloads/signalp-4.1 # Path to SignalP

# Do not touch this
CBS=http://www.cbs.dtu.dk
BIOMINE='http://biomine.cs.vcu.edu/servers/biomine.php'

################################################################################
#                     THE SCRIPT BEGIN HERE
################################################################################
# Script name variable
SCRIPT=`basename ${BASH_SOURCE[0]}`
# Help function
function HELP {
  echo "Help documentation for ${SCRIPT}\n"
  echo "Basic usage: $SCRIPT [-d][-h] file.ext\n"
  echo "Other command line options :"
  echo "-d  --Directory where you want the outputs stored (default directory : ./outputs)."
  echo "-h  --Displays this help message. No further functions are performed.\n"
  echo "Example: $SCRIPT file.ext\n"
  exit 1
}
# Set arguments
OUTDIR=$DATA/outputs;
NUMARGS=$#
if [ $NUMARGS -eq 0 ]; then
  HELP
fi
while getopts :d:h FLAG ; do
  case $FLAG in
    d) OUTDIR=$OPTARG;;
    h) HELP ;;
    \?)
      echo "Use $SCRIPT -h to see the help documentation.\n"
      exit 2 ;;
  esac
done
shift $((OPTIND-1))
if [[ $# -eq 0 ]]; then
  echo "No file provided.\nUse $SCRIPT -h to see the help documentation.\n"
  exit;
else FILE=$1
fi
if [ ! -e $FILE ]; then
  echo "$FILE does not exist. Please provide a valid file name.\n"
  exit;
fi
echo "Please enter an email address if you want to get some results by email (default : yu@yi.yo) :"
# read email
email='yu@yi.yo'
mkdir $OUTDIR
echo "\nThe script will automatically run several proteomics prediction softwares.
The job might take a few minutes, so just go take a coffee, relax and come back later.\n"

################################################################################
#                     LOCAL PROGRAMS
################################################################################

cat $1 | grep ">" | sed 's%>%%' | sed 's%.fasta%%' | sed 's%:%_%' | awk '{print $1}' > $OUTDIR/names
filename=$(basename $1)

# echo ...Splitting your big fasta file...
# # In case of big file split it into even smaller files
# if [ $(cat $OUTDIR/names | wc -l) -gt 99 ]; then
#   mkdir ${filename%.*}_split
#   split -l 98 $1
#   for f in $(ls x*); do mv $f ${filename%.*}_split/. ; done
# fi
# # Split multiline fasta file into single sequence files
FOLDER="${filename%.*}_single"
# rm -rf $FOLDER
# mkdir $FOLDER
# pyfasta split --header "$FOLDER/%(seqid)s.fasta" $FILE;
# rm -rf **/*.flat *.flat **/*.gdx *.gdx # Remove useless files
# echo "Results saved in $FOLDER directory."

# # Motif Search with Python script
# echo ...Performing ELM motif search...
# python $DATA/src/motifs.py $DATA/data/elm_classes.tsv $FILE $OUTDIR/Motifs.tsv
# echo "Results saved in $OUTDIR/Motifs.tsv."

# # ANCHOR binding regions prediction
# echo ...Running ANCHOR...
# rm -rf $OUTDIR/anchor_tmp
# mkdir $OUTDIR/anchor_tmp
# for spep in $(find ./$FOLDER -name '*.fasta');
# do
#   bname=$(basename $spep)
#   name=$( echo ${bname%.*} | sed 's%:%_%'  | awk '{print $1}')
#   # echo '# '$name >> $OUTDIR/Anchor.txt
#   # echo '# '$(basename $spep) >> $OUTDIR/Anchor.txt
#   perl $ANCHOR/anchor.pl -v $spep > $OUTDIR/anchor_tmp/$name
# done
# echo "Results temporarily saved in anchor_tmp directory"
#
# # SignalP prediction of signal peptide
# echo ...Running SignalP...
# $SIGNALP/signalp -t euk -f short $FILE > $OUTDIR/SignalP.txt;
# echo "Results saved in $OUTDIR/SignalP.txt."

# ################################################################################
# #                     SERVERS
# ################################################################################
#
# # TargetP
# echo ...Sending data to CBS TargetP server...
# curl -i -s -F SEQSUB=@$1 -F configfile="/usr/opt/www/pub/CBS/services/TargetP-1.1/TargetP.cf" $CBS/cgi-bin/webface2.fcgi > tmp
# url=$(cat tmp | grep -Eo 'cgi-bin/webface2.fcgi[a-zA-Z0-9\?\=]+')
# status=1
# while [ $status -ne 0 ]; do
#   sleep 5s
#   curl -s $CBS/$url > $OUTDIR/TargetP.txt
#   status=$(cat $OUTDIR/TargetP.txt | grep -Eo 'Your job [A-Z0-9]+ is' | wc -m)
# done
# echo Cleaning output file...
# grep -v -E "<|>" $OUTDIR/TargetP.txt > tmp && mv tmp $OUTDIR/TargetP.txt
# echo "Results saved in $OUTDIR/TargetP.txt"
#
# # TMHMM
# echo ...Sending data to CBS TMHMM server...
# curl -i -s -F seqfile=@$1 -F outform='-short' -F configfile="/usr/opt/www/pub/CBS/services/TMHMM-2.0/TMHMM2.cf" $CBS/cgi-bin/webface2.fcgi > tmp
# url=$(cat tmp | grep -Eo 'cgi-bin/webface2.fcgi[a-zA-Z0-9\?\=]+')
# status=1
# while [ $status -ne 0 ]; do
#   sleep 5s
#   curl -s $CBS/$url > $OUTDIR/TMHMM.txt
#   status=$(cat $OUTDIR/TMHMM.txt | grep -Eo 'Your job [A-Z0-9]+ is' | wc -m)
# done
# echo Cleaning output file...
# grep -v -E "<|>|with" $OUTDIR/TMHMM.txt > tmp && mv tmp $OUTDIR/TMHMM.txt
# echo "Results saved in $OUTDIR/TMHMM.txt"

# # NetSurfP
# echo ...Sending data to CBS NetSurfP server...
# curl -i -s -F uploadfile=@$1 -F configfile="/usr/opt/www/pub/CBS/services/NetSurfP-1.1/NetSurfP.cf" $CBS/cgi-bin/webface2.fcgi > tmp
# url=$(cat tmp | grep -Eo 'cgi-bin/webface2.fcgi[a-zA-Z0-9\?\=]+')
# status=1
# while [ $status -ne 0 ]; do
#   sleep 20s
#   curl -s $CBS/$url > $OUTDIR/NetSurfP.txt
#   status=$( cat $OUTDIR/NetSurfP.txt | grep -Eo 'Your job [A-Z0-9]+ is' | wc -m )
# done
# rm tmp
# echo Cleaning output file...
# grep -v -E "<|>" $OUTDIR/NetSurfP.txt > tmp && mv tmp $OUTDIR/NetSurfP.txt
# echo "Results saved in $OUTDIR/NetSurfP.txt"
#
# # DRNApred for DNA and RNA binding sites
# echo ...Sending data to Biomine DRNApred server...
# if [ $(cat $OUTDIR/names | wc -l) -gt 99 ]; then
#   rm -f $OUTDIR/DRNApred.txt
#   for f in $(ls ${filename%.*}_split/); do
#     echo $f
#     curl -s -F data=@${filename%.*}_split/$f -F email1=$email $BIOMINE'?name=DRNApred' > ${f}_tmp
#     url=$(cat ${f}_tmp | grep -Eo -m 1 $BIOMINE'\?id=[0-9]+' | awk '{print $1 }' )
#     echo $url
#     if [ $(echo $url | wc -m) -ne 0 ]; then
#       status=0
#       while [ $status -eq 0 ]; do
#         echo $status
#         sleep 10
#         curl -s $url > ${f}_tmp
#         status=$(cat ${f}_tmp | grep -Eo 'Your results are ready' | wc -m)
#       done
#       cat ${f}_tmp | grep -Eo -m 1 'href="http://biomine.cs.vcu.edu/webresults/DRNApred/[0-9]+'\
#       | sed 's%href="%%' | awk '{print $1 "/DRNA.pred" }' | xargs curl -s  >> $OUTDIR/DRNApred.txt
#       rm ${f}_tmp
#     fi
#   done
# else
#   curl -s -F data=@$1 -F email1=$email $BIOMINE'?name=DRNApred' > tmp
#   url=$(cat tmp | grep -Eo -m 1 $BIOMINE'\?id=[0-9]+' | awk '{print $1 }' )
#   if [ $(echo $url | wc -m) -ne 0 ]; then
#     status=0
#     while [ $status -eq 0 ]; do
#       sleep 10
#       curl -s $url > tmp
#       status=$(cat tmp | grep -Eo 'Your results are ready' | wc -m)
#     done
#     cat tmp | grep -Eo -m 1 'href="http://biomine.cs.vcu.edu/webresults/DRNApred/[0-9]+'\
#     | sed 's%href="%%' | awk '{print $1 "/DRNA.pred" }' | xargs curl -s  > $OUTDIR/DRNApred.txt
#     rm tmp
#   fi
# fi

# JPred secondary structure prediction
# This one is particularly long and won't work with very short sequences
echo ...Running JPred API...
csh bin/jpred.csh $FOLDER $OUTDIR bin/jpredapi;
echo "Results available in $FOLDER_output folder."

################################################################################
#                     RESULTS
################################################################################
#
mkdir cbs-anchor
echo "...Analysing results..."
rm -rf $OUTDIR/cbs-anchor
mkdir $OUTDIR/cbs-anchor
names=$(ls $OUTDIR/anchor_tmp | awk '{print $1}')
for name in $(cat $OUTDIR/names); do
  python src/genplot.py $name $OUTDIR > $OUTDIR/cbs-anchor/$name
done
echo "Plots generated."

echo "...Adding results to database and generating html file..."
python src/SQLscript.py $OUTDIR $(cat $OUTDIR/names | wc -l) > Slavoff_sEPs.html
echo "Done. You can find a html table with all results summarized :)"
# rm -rf $FOLDER
# rm -rf $OUTDIR/anchor_tmp
