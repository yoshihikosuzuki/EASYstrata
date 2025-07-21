#!/bin/bash
#SBATCH --account=youraccount
#SBATCH --time=03:00:00
#SBATCH --job-name=gmap
#SBATCH --output=log_gmap-%J.out
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

# Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

#activate conda

#give expected number of arguments from config file:
source config/config
source config/colors 

##  ------------------------ general parameters --------------------------------  ##
haplotype=$1 #the same name from the config user for haplotype1 or haplotype2
folder=$2    #either haplo1 or haplo2
genome=$3    #full path to the genome associated with haplo1 or haplo2

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
if [ -z "$genome" ]; then
    echo "error no $genome name provided! "
    echo "please see readme.md"
    exit 1
fi
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
   mkdir -p "$folder"/03_genome
   cd "$folder"/03_genome || exit 1 
   cp "$genome" "$haplotype".fa
   cd ../../
   #modifiy in config file as well:
   sed -i -E '/^genome/ s/.gz//g' config/config
else
   echo "$genome is not gzipped"
   genome=$genome
   mkdir -p "$folder"/03_genome
   cd "$folder"/03_genome || exit 1 
   cp "$genome" "$haplotype".fa
   cd ../../
fi

cd "$folder" || exit 1 

genome=03_genome/"$haplotype".fa

#launch gmap :
../00_scripts/02_gmap.sh "$genome" 
