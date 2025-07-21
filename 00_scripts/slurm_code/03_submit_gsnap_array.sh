#!/bin/bash
#SBATCH --account=youraccount
#SBATCH --time=36:00:00
#SBATCH --job-name=gsnap
#SBATCH --output=log_gsnap-%J.out
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --array=1-24 #set size according to number of chromosome for instance

# Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

#activate conda

#activate the env

#give expected number of arguments from config file:
source config/config
source config/colors 

##  ------------------------ general parameters --------------------------------  ##
haplotype=$1 #the same name from the config user for haplotype1 or haplotype2
folder=$2    #either haplo1 or haplo2
listfile=$3  #list of files to be mapped 
if [ -z "$haplotype" ]; then
    echo "error no $haplotype name provided! "
    echo "this will be used as a short name for subsequent analyses"
    echo "please see readme.md"
    exit 1
fi
if [ -z "$folder" ]; then
    echo "error no $folder name provided! "
    echo "this will be the output folder in which the aanlysis will be located"
    echo "this should be either haplo1 or haplo2"
    echo "please see readme.md"
    exit 1
fi

echo rnaseqlist is "$RNAseqlist" 

if [ -n "${genome}" ] ; 
then
    b1=$(basename "${genome%.fa*}" )
    if [[ "$b1" =~ [^a-zA-Z0-9._] ]] ; 
    then 
        echo "error only alphanumeric character allowed in genomeIDs" 
        echo "see readme.md on github"
        #exit 
    else 
        echo "genome 1 is $genome" 
    fi
fi
if file --mime-type "$genome" | grep -q gzip$; then
   echo "$genome is gzipped"
   gunzip "$genome"
   genome=${genome%.gz}
   cd "$folder"/03_genome || exit 1 
   cp "$genome" "$haplotype".fa
   cd ../../
   #modifiy in config file as well:
   sed -i -E '/^genome/ s/.gz//g' config/config
else
   echo "$genome is not gzipped"
   genome=$genome
   cd "$folder"/03_genome || exit 1 
   cp "$genome" "$haplotype".fa
   cd ../../
fi

cd "$folder" || exit 1 

genome=03_genome/"$haplotype".fa

TIME=$(date +%Y-%m-%d_%Hh%Mm%Ss)
LOG_FOLDER="LOGS"
#create folder if not existent:
if [ ! -d "$LOG_FOLDER" ] ; then mkdir $LOG_FOLDER ; fi 

ncol=$(awk '{print NF}'  "$RNAseqlist" |uniq)

if [[ $ncol = 2 ]] ; then
    echo -e "\nrunning gsnap assuming reads are Paired-End\n" 
    #launch gsnap - samtools and read count:
    for read1 in ../02_trimmed/*R1.paired.fastq.gz  ; do
        [ -e "$read1" ] || continue 
        ../00_scripts/03_gsnap_PE.sh "$genome" "$read1" 2>&1 |tee "$LOG_FOLDER"/gsnap_"$(basename "$read1")"_"$TIME".log
        if [[  "${PIPESTATUS[0]}" -ne 0 ]]
        then
            echo "ERROR GSNAP FAILED"
            exit 1
        fi
    done 
else
    #assuming SE:
    #launch gsnap - samtools and read count:
    echo -e "\nrunning gsnap assuming reads are Single-End\n"

    for read1 in ../02_trimmed/*R1.fastq.gz ; do
        [ -e "$read1" ] || continue 
        ../00_scripts/03_gsnap_SE.sh "$genome" "$read1" 2>&1 |tee "$LOG_FOLDER"/gsnap_"$(basename "$read1")"_"$TIME".log
        if [[  "${PIPESTATUS[0]}" -ne 0 ]]
        then
            echo "ERROR GSNAP FAILED"
            exit 1
        fi
    done 
fi
