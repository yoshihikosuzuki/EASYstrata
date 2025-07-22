#!/bin/bash
##SBATCH --account=youraccount
#SBATCH --time=12:00:00
#SBATCH --job-name=gsnap
#SBATCH --output=log_gsnap-%J.out
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --array=1-3 #set according to number of RNAseq samples
# Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

#activate conda

#give expected number of arguments from config file:
source config/config
##  ------------------------ general parameters --------------------------------  ##
haplotype=$1 #the same name from the config user for haplotype1 or haplotype2
folder=$2    #either haplo1 or haplo2
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
cd "$folder" || exit 1 

genome=03_genome/"$haplotype".fa

TIME=$(date +%Y-%m-%d_%Hh%Mm%Ss)
LOG_FOLDER="LOGS"
#create folder if not existent:
if [ ! -d "$LOG_FOLDER" ] ; then mkdir $LOG_FOLDER ; fi 

ncol=$(awk '{print NF}'  "$RNAseqlist" |uniq)
awk '{print $1}'  "$RNAseqlist" > samples

wanted_samples=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples )

if [[ $ncol = 2 ]] ; then
    echo -e "\nrunning gsnap assuming reads are Paired-End\n" 
    #launch gsnap - samtools and read count:
        ../00_scripts/03_gsnap_PE_array.sh "$genome" "$wanted_samples" 2>&1 |tee "$LOG_FOLDER"/gsnap_"$(basename "$read1")"_"$TIME".log
        if [[  "${PIPESTATUS[0]}" -ne 0 ]]
        then
            echo "ERROR GSNAP FAILED"
            exit 1
        fi
else
    #assuming SE:
    #launch gsnap - samtools and read count:
    echo -e "\nrunning gsnap assuming reads are Single-End\n"
        [ -e "$read1" ] || continue 
        ../00_scripts/03_gsnap_SE_array.sh "$genome" "$wanted_samples" 2>&1 |tee "$LOG_FOLDER"/gsnap_"$(basename "$read1")"_"$TIME".log
        if [[  "${PIPESTATUS[0]}" -ne 0 ]]
        then
            echo "ERROR GSNAP FAILED"
            exit 1
        fi
fi

if [ -e samples ]; then rm samples ; fi
