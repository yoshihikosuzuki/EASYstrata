#!/bin/bash
##SBATCH --account=youraccount
#SBATCH --time=78:00:00
#SBATCH --job-name=repeatmodeler
#SBATCH --output=log_repeatmodeler-%J.out
#SBATCH --mem=40G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
# Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

#activate conda
#source /local/env/envconda3.sh 
#activate the env
#mamba activate repeatmodeler_env
#give expected number of arguments from config file:
source config/config
##  ------------------------ general parameters --------------------------------  ##
while [ $# -gt 0 ] ; do
  case $1 in
    -g | --genome) genome="$2" ;echo "the genome file  is: $genome" >&2;;
    -s | --haplotype) haplotype="$2" ;echo "the haplotype name will be $haplotype" >&2;;
    -f | --folder ) folder="$2"; echo "the data located in $folder will be processed">&2;;
    -h | --help) echo -e "Option required:
    -g/--genome \t the reference genome file 
    -s/--haplotype\t the haplotype name (used for database building and basename in several steps)
    -f/--folder\t the folder where the genome will be located (either haplo1 or haplo2). the genome will be copied here
    Optional:
    -m/--mask\t a string stating wether unknown TE should be removed for genome annotation (YES/NO) -- default: YES
    " >&2;exit 1;;
    esac
    shift
done
if [ -z "$genome" ] || [ -z "$haplotype" ] ; then
    echo -e  >&2 "Fatal error: Ref genome (-g), and haplotype name (-s) not defined\n
    see manual with -h or --help"
exit 2
fi
#optional parameters:
if [ -z "$Mask" ] ; then
    Mask=YES
fi
if [ -z "$folder" ]; then
    folder=haplo1
    mkdir "$folder"
    mkdir "$golder"/03_genome
    echo "warning, no folder path provided, the data will be located in haplo1"
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
   if [ ! -e $folder/03_genome ] ; then mkdir $folder/03_genome ; fi
   cd "$folder"/03_genome || exit 1 
   if [ ! -s  "$haplotype".fa ]
   then
       cp "$genome" "$haplotype".fa
   fi
   cd ../../
   #modifiy in config file as well:
   sed -i -E '/^genome/ s/.gz//g' config/config
else
   echo "$genome is not gzipped"
   genome=$genome
   if [ ! -e $folder/03_genome ] ; then mkdir $folder/03_genome ; fi
   cd "$folder"/03_genome || exit 1 
   if [ ! -s  "$haplotype".fa ]
   then
       cp "$genome" "$haplotype".fa
   fi
   cd ../../
fi
LOG_FOLDER="LOGS"
if [ ! -d "$LOG_FOLDER" ] ; then mkdir $LOG_FOLDER ; fi 

cd "$folder" || exit 1
# -------------------- run repeatmodeler ---------------------- #
#setting up path prior to running busco: 
rm_file="03_genome/genome.wholemask.fa"
if [  -e "$rm_file" ] ; then 
    echo -e "\n-----------------------------------------------------"
    echo -e "\trepeatmodeler output already exist\n\twill skip this step"; 
    echo -e "-----------------------------------------------------\n"
else 
    echo -e "\n-----------------------------------------------------"
    echo -e "no repeatmodeller output \n will launch repeatmodeller " ; 
    echo -e "-----------------------------------------------------\n"
    #remove any empty rm_file
    rm -rf $rm_file 

    ../00_scripts/05_repeatmodeler.sh 03_genome/"$haplotype".fa "$haplotype" "$Mask" 2>&1 |tee LOGS/log_rm
    if [[  "${PIPESTATUS[0]}" -ne 0 ]]
    then
        echo -e "${RED} ERROR: repeatmodeler failed.\n
        check the provided libraries and software dependancies  \n${NC}"   
        exit 
    else
        echo -e "\n${BLU}---- repeatmodeler run successfull ----\n${NC}"
        fi
fi
