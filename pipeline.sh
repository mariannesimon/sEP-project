#!/bin/bash

# This script will automatically run prediction softwares.
# Some require installation and are available at following links :
# ANCHOR : http://anchor.enzim.hu/Downloads.php
# JPred4 : http://www.compbio.dundee.ac.uk/jpred/api.shtml#download
# SignalP : http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?signalp
# Psipred : http://bioinfadmin.cs.ucl.ac.uk/downloads/psipred/

# Set path variables
DATA=$(pwd); # Path to your fasta file
ANCHOR=/usr/local/bin/anchor; # Path to Anchor software
JPRED=~/Downloads/jpredMassSubmitSchedule; # Path to JPred software
SIGNALP=~/Downloads/signalp-4.1; # Path to SignalP
PSIPRED=~/Downloads/psipred/BLAST+ # Path to Psipred

# Set Script name variable
SCRIPT=`basename ${BASH_SOURCE[0]}`
# Help function
function HELP {
  echo "Help documentation for ${SCRIPT}\n"
  echo "Basic usage: $SCRIPT [h] file.ext\n"
  echo "Other command line options :"
  # echo "-p  --Numbers of files to parse file into (in case of big file)."
  echo "-h  --Displays this help message. No further functions are performed.\n"
  echo "Example: $SCRIPT file.ext\n"
  exit 1
}
# Set arguments
NUMARGS=$#
if [ $NUMARGS -eq 0 ]; then
  HELP
fi
while getopts :h FLAG ; do
  case $FLAG in
    # p) PARSE=$OPTARG;;
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

echo "\nThe script will automatically run several proteomics prediction softwares.
Please be patient, the job might take a few minutes.\n"

# cat $1 | grep ">" | awk '{print $1}' > $1.names;
# cat $1 | grep -v "^>" | awk '{print $1}' > $1.seqs;

# Split multiline fasta file into several single sequence files
echo Splitting sPEP fasta file...
FOLDER="${FILE%.*}_single"
rm -rf $FOLDER
mkdir $FOLDER
pyfasta split --header "$FOLDER/%(seqid)s.fasta" $FILE;
rm -rf *.flat *.gdx # Remove useless files
echo "Results saved in $FOLDER directory."

# Motif Search with Python script
echo Performing simple motif search...
python $DATA/src/motifs.py $DATA/data/motifs_elm.txt $FILE /outputs/Motifs.txt
echo "Results saved in outputs/Motifs.txt."

# ANCHOR binding regions prediction
echo Running ANCHOR...
rm -f outputs/Anchor.txt
for spep in $(find ./$FOLDER -name '*.fasta');
do
  echo '>'$(basename $spep) >> outputs/Anchor.txt
  perl $ANCHOR/anchor.pl $spep >> outputs/Anchor.txt;
done
echo "Results saved in outputs/Anchor.txt."

# SignalP prediction of signal peptide
echo Running SignalP...
$SIGNALP/signalp -t euk -f short $FILE > outputs/SignalP.txt;
echo "Results saved in outputs/SignalP.txt."

# Psipred structure prediction
echo Running Psipred...
# $PSIPRED/runpsipredplus $FILE ;
echo "Results saved in outputs/Psipred.txt."

# JPred secondary structure prediction
echo Running JPred API...
# This one is particularly long and won't work with very short sequences
cd $JPRED;
./massSubmitScheduler.csh $DATA/$FOLDER;
cd $DATA;
echo "Results available online. Please find the links in sPEP_single_output folder."

# Clean data
rm -rf $FOLDER;

# Results summary : one single file
echo "Analysing results..."

echo "Results summary saved in outputs/Sum.txt"
