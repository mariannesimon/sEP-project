#!/bin/bash

# This script will automatically run prediction softwares.
# Some require installation and are available at following links :
# Pyfasta (Python 2.7) : https://pypi.python.org/pypi/pyfasta
# Biopython : http://biopython.org/wiki/Download
# ANCHOR : http://anchor.enzim.hu/Downloads.php
# SignalP : http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?signalp

# Set your programs path variables
ANCHOR=/usr/local/bin/anchor # Path to Anchor software
SIGNALP=/usr/local/bin/signalp-4.1 # Path to SignalP

# Do not touch this
CBS=http://www.cbs.dtu.dk
BIOMINE='http://biomine.cs.vcu.edu/servers/biomine.php'

################################################################################
#                     THE SCRIPT BEGINS HERE
################################################################################

############# Script name variable #############
SCRIPT=`basename ${BASH_SOURCE[0]}`
############# Help function #############
function HELP {
  echo "Help documentation for ${SCRIPT}\n"
  echo "Basic usage: $SCRIPT [-d][-h] file.ext\n"
  echo "Other command line options :"
  echo "-d  --Directory where you want the outputs stored (default directory : ../outputs)."
  echo "-h  --Displays this help message. No further functions are performed.\n"
  echo "Example: $SCRIPT file.ext\n"
  exit 1
}
############ Set arguments #############
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
email='yu@yi.yo' # Can be changed to received some results by email
echo "\nThe script will automatically run several proteomics prediction softwares.
The job might take a few minutes (or a few hours), so have a coffee and come back later.\n"

################################################################################
#                     LOCAL PROGRAMS
################################################################################

filename=$(basename $1)
if [ ${#OUTDIR} -eq 0 ]; then OUTDIR=$(dirname $0)/../${filename%.*}_output; fi
mkdir $OUTDIR
cat $1 | grep ">" | sed 's%>%%' | sed 's%.fasta%%' | sed 's%:%_%' | awk '{print $1}' > $OUTDIR/names

############# Splits big files into smaller files #############

echo ...Splitting your big fasta file...
if [ $(cat $OUTDIR/names | wc -l) -gt 99 ]; then
  mkdir ${filename%.*}_split
  split -l 98 $1
  for f in $(ls x*); do mv $f ${filename%.*}_split/. ; done
fi

############# Split multiline fasta file into single sequence files #############

FOLDER="${filename%.*}_single"
rm -rf $FOLDER
mkdir $FOLDER
pyfasta split --header "$FOLDER/%(seqid)s.fasta" $FILE;
rm -rf **/*.flat *.flat **/*.gdx *.gdx # Remove useless files
echo "Results saved in $FOLDER directory."

# ############# Protein content analysis with ProtParam #############
#
# echo "...Analyzing your protein sequences..."
# python $(dirname $0)/protAna.py $OUTDIR $1
# echo "Results saved in $OUTDIR/protParams.tsv"
#
# ############# Linear Motif Search  #############
#
# echo ...Performing ELM motif search...
# python $(dirname $0)/motifs.py -m $(dirname $0)/../data/elm_classes.tsv -i $FILE -o $OUTDIR/Motifs.tsv
# echo "Results saved in $OUTDIR/Motifs.tsv."
#
# ############# ANCHOR binding regions prediction #############
#
# echo ...Running ANCHOR...
# echo '' > $OUTDIR/Anchor.txt
# # rm -rf $OUTDIR/anchor_tmp
# # mkdir $OUTDIR/anchor_tmp
# for spep in $(find ./$FOLDER -name '*.fasta');
# do
#   bname=$(basename $spep)
#   name=$( echo ${bname%.*} | sed 's%:%_%'  | awk '{print $1}')
#   # echo '>'$name >> $OUTDIR/Anchor.txt
#   # echo '# '$(basename $spep) >> $OUTDIR/Anchor.txt
#   # perl $ANCHOR/anchor.pl -v $spep > $OUTDIR/anchor_tmp/$name
#   perl $ANCHOR/anchor.pl -v $spep > $OUTDIR/anchor_tmp
#   cat $OUTDIR/anchor_tmp | grep -v '^#' | awk -v n=$name '{ print $0 FS n }' >> $OUTDIR/Anchor.txt
# done
# rm $OUTDIR/anchor_tmp
# echo "Results saved in $OUTDIR/Anchor.txt directory"
#
# ############# SignalP prediction of signal peptide #############
#
# echo ...Running SignalP...
# $SIGNALP/signalp -t euk -f short $FILE > $OUTDIR/SignalP.txt;
# echo "Results saved in $OUTDIR/SignalP.txt."
#
# ################################################################################
# #                     SERVERS
# ################################################################################
#
# ############# TargetP #############
#
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
# ############# TMHMM #############
#
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
#
# ############# NetSurfP #############
#
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
# ############ DRNApred for DNA and RNA binding sites #############
#
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
#
# # ############# JPred secondary structure prediction #############
#
# This one is particularly long and won't work with very short sequences
echo ...Running JPred API...
csh bin/jpred.csh $FOLDER $OUTDIR bin/jpredapi;
echo "Results available in $FOLDER_output folder."

# ################################################################################
# #                     RESULTS
# ################################################################################
#
# ############# Python plots #############
#
# echo "...Analysing results..."
# mkdir $OUTDIR/summary
# python $(dirname $0)/summary.py $OUTDIR

############# Update database and generate summary files #############

echo "...Adding results to SQLite database and generating html file..."
sql=$(dirname $0)/../sqlite
current=$(pwd)
cd $OUTDIR
ABSPATH=$(pwd)
cd $current
python $sql/main.py -d $ABSPATH -o $OUTDIR/sEP_table.csv -f csv -db $sql/spepDB.db -t spep
python $sql/main.py -o $OUTDIR/sEP_table.html -f html -db $sql/spepDB.db -t spep
echo "Done. You can find a html table with all results summarized :)"

############# R plots #############

mkdir $OUTDIR/Rplots
Rscript $(dirname $0)/resPlots.R $OUTDIR $OUTDIR/sEP_table.csv $(dirname $0)/../data/aaproperties.tsv
echo "Plots generated."

############# Deleting remaining files #############

rm -rf $FOLDER ${filename%.*}_split/
