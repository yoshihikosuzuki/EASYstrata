#!/bin/bash
#PURPOSE: script to run braker
#AUTHOR: QR
#Date updated: 22-07-2025
#------------- EXTERNAL VARIABLE FROM CONFIG FILE -------------- #
source ../config/config
source ../config/cpu_mem
#------------- CONDA ACTIVATION  -------------- #
#eval "$(conda shell.bash hook)"
#conda activate superannot
#--- start of setting path ---- " 
CDB_PATH
TSEBR_PATH
PROTH_PATH
GMARK_PATH
AUGCO_PATH
AUGBI_PATH
AUGSC_PATH
#--- end of setting path ---- " 
############################################################
# ERROR TRACKING.                                          #
############################################################
set -eE -o functrace
failure() {
  local lineno=$1
  local msg=$2
   echo "command failed at line $lineno: $msg"
}
trap 'failure ${LINENO} "$BASH_COMMAND"' ERR
##################################################
#  ---- external data required arguments --------------- #
if (( $# < 3 )) ; then
    echo "USAGE: $0 <reference_genome> <species>  <fungus>(YES/NO)" 
    echo -e "Expecting the following parameters:\n
        1 - the reference genome\n
        2 - a species name\n
        3 - a string YES/NO stating wether data are from fungus or not \n
	4 - round the round ID\n" 
        exit 1
else
    genome=$1 #the reference genome here! "genome.wholemask.fa"
    species=$2
    fungus=$3
    round=$4
fi
NCPUS="$NCPUS_BRAKER"
TIME=$(date +%Y-%m-%d_%Hh%Mm%Ss)
#  --------- step 2 : BRAKER WITH REFERENCE DATABASE USING FIVE ROUNDS --------- ## 
FOLDER=06_braker/round"$round"_braker_on_refprot 
if [ ! -d $FOLDER ] ; then
    mkdir -p $FOLDER 
fi
echo -e "\n\n----------- round $round ------------\n\n" 

output="braker.gtf"

relatProt="relatProt.fa"

wd=${FOLDER}
if [ -f "$wd"/"$output" ]
then
    echo "file $output round $round already exist will skip the run"
else
   rm -rf "${wd:?}"/*
   if [[ $fungus = "YES" ]]
   then
       braker.pl --species="$species"_"$TIME"_round"$round"  --genome="$genome" --threads="$NCPUS" \
           --softmasking --prot_seq=$relatProt --workingdir=$wd --fungus  
   else
       braker.pl --species="$species"_"$TIME"_round"$round"  --genome="$genome" --threads="$NCPUS" \
           --softmasking --prot_seq=$relatProt --workingdir=$wd 
   fi
fi
