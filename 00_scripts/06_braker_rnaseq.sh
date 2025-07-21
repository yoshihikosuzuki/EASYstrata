#!/bin/bash
#PURPOSE: script to run braker
#AUTHOR: QR
#Date updated: 10-03-2023

#------------- EXTERNAL VARIABLE FROM CONFIG FILE -------------- #
source ../config/config
source ../config/cpu_mem
#------------- CONDA ACTIVATION  -------------- #
#eval "$(conda shell.bash hook)"
#conda activate superannot
#conda activate /scratch/qrougemont/new_superannot/
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

# keep track of the last executed command
#trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
#trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

#  ---- external data required arguments --------------- #
if (( $# < 4 )) ; then
    echo "USAGE: $0 <reference_genome> <species> <RNAseq>(YES/NO) <fungus>(YES/NO) <bamlist> (optional) <NCPUS>(optional)" 
    echo -e "Expecting the following parameters:\n
        1 - the reference genome\n
        2 - a species name\n
        3 - a string YES/NO for wether RNAseq should be considered or not\n\n
        4 - a string YES/NO stating wether data are from fungus or not \n\n" 
        exit 1
else
    genome=$1 #the reference genome here! "genome.wholemask.fa"
    species=$2
    RNAseq=$3 #YES/NO
    fungus=$4
    bamlist=$5
    NCPUS=$6
fi

NCPUS="$NCPUS_BRAKER"

TIME=$(date +%Y-%m-%d_%Hh%Mm%Ss)

#----------------- BAM data for RNAseq ------------------------ #
#previously we used a concatenated and filtered bam.
#however, providing a list of all bam gives similar results with a significant gain in time:
if [ -z "$bamlist" ] ; then
    alnBAM=$(echo 04_mapped/*sorted.bam |sed 's/ /,/g' )
else

    #assuming list of bam already exist !
    #these will be soft linked in 04_mapped
    if [ ! -d "04_mapped" ] ; then
           mkdir 04_mapped
    fi

    cd 04_mapped || exit 

    #for i in $(cat "$bamlist" ) ; do
    while read -r line 
    do 
        ln -s "$line" . 
    done < "$bamlist" 
    cd ../ ; 

    alnBAM=$(echo 04_mapped/*sorted.bam |sed 's/ /,/g' )
fi

## --------- step 1 : BRAKER WITH RNA SEQ  ---------  ##
## Note on braker: we found that running on RNAseq & proteins separately and combining afterwards 
## manually with TSEBRA give higher BUSCO scores. 
## in addition manually setting parameter for TSEBRA config file may prevent the removal of 
## biologically important genes, a case we have observed on all out dataset.
## Therefore we choose to only implement this approach
## for user having trouble installing braker and all it's dependencies it should still be possible
## to modify the code lightly and run through singularity.

#Run braker with RNAseq 

output="braker.gtf"

if [[ $RNAseq = "YES" ]]
then
    wd=06_braker/rnaseq
    if [ -f "$wd"/"$output" ]
    then
        echo -e "file $output RNAseq already exist will skip the run\n"
    else
        echo "starting annotation with RNAseq data only"
        mkdir -p $wd

        if [[ $fungus = "YES" ]]
        then
            echo -e "------ \n running braker on rnaseq data \n -------"
            echo -e "------ \n data are from fungus \n -------"
                braker.pl --species="$species"_"$TIME"_rnaseq --species="$species"_"$TIME"_rnaseq --fungus \
                    --genome="$genome" --threads="$NCPUS"  --softmasking --bam="$alnBAM" --workingdir=$wd 
        else
            echo -e "------ \n running braker on rnaseq data \n -------"
                braker.pl --species="$species"_"$TIME"_rnaseq --species="$species"_"$TIME"_rnaseq\
                    --genome="$genome" --threads="$NCPUS"  --softmasking --bam="$alnBAM" --workingdir=$wd 
        fi
    fi 
fi
