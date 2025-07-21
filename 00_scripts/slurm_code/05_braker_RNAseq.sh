#!/bin/bash
##SBATCH --account=youraccount
#SBATCH --time=74:00:00
#SBATCH --job-name=braker
#SBATCH --output=log_braker-%J.out
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32

# Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

#activate conda
#source /local/env/envconda3.sh 
#activate the env
mamba activate superannot/

#give expected number of arguments from config file:
source config/config
source config/cpu_mem
##  ------------------------ general parameters --------------------------------  ##
while [ $# -gt 0 ] ; do
  case $1 in
    -g | --genome) genome="$2" ;echo "the genome file  is: $genome" >&2;;
    -s | --haplotype) haplotype="$2" ;echo "the haplotype name will be $haplotype" >&2;;
    -o | --folder ) folder="$2"; echo "the data located in $folder will be processed">&2;;
    -b | --bam )    bamlist="$2" ; echo "the optional bamlist of bam files will be $bamlist"  >&2;;
    -r | --rna )    RNAseq="$2" ; echo "Is RNAseq provided ? $RNAseq " >&2;; 
    -f | --fungus ) fungus="$2" ; echo "Do species belong to fungi? $fungus" >&2;;
    -h | --help) echo -e "Option required:
    -g/--genome \t the reference genome file 
    -s/--haplotype\t the haplotype name (used for database building and basename in several steps)
    -o/--output\t the folder where the genome will be located (either haplo1 or haplo2). the genome will be copied here
    Optional:
    -r/--rnaseq\t (a string YES/NO stating wether RNAseq should be used
    -b/--bam\t (a full path to a txt file specifiying a list of bam files, if available
        " >&2;exit 1;;
    esac
    shift
done
if [ -z "$genome" ] || [ -z "$haplotype" ] || [ -z "$folder" ] ; then
    echo -e  >&2 "Fatal error: Ref genome (-g), haplotype name (-s) and output folder (-o haplo1 or haplo2) not defined\n
    see manual with -h or --help"
exit 2
fi
if [ -z "$RNAseq" ] ; then 
  RNAseq=NO
fi
if [ -z "$fungus" ] ; then 
  fungus=NO
fi

# -------------------- run Braker  ---------------------------- #
#setting up path prior to running busco: 
tsebrapath=$(command -v tsebra.py |xargs dirname )
cdbpath=$(command -v cdbfasta |xargs dirname )
protpath=$(command -v prothint.py |xargs dirname)
gmarkpath=$(command -v gmes_petap.pl |xargs dirname)
augbin=$(command -v augustus |xargs dirname)
augscripts="${augbin//bin/scripts}"
augconf="${augbin//bin/config}"

#verify again that all path exist -----
[[ -z "$tsebrapath" ]] && { echo "Error: tsebra.py not found"; exit 1; }
[[ -z "$cdbpath" ]] && { echo "Error: cdbfasta not found"; exit 1; }
[[ -z "$protpath" ]] && { echo "Error: prothint not found"; exit 1; }
[[ -z "$gmarkpath" ]] && { echo "Error: genemark not found"; exit 1; }
[[ -z "$augbin" ]] && { echo "Error: Augustus binaries not found"; exit 1; }

# reshape braker code prior to run:
sed -i "s#CDB_PATH#export CDBTOOLS_PATH=$cdbpath#" 00_scripts/06_braker_rnaseq.sh
sed -i "s#TSEBR_PATH#export TSEBRA_PATH=$tsebrapath#" 00_scripts/06_braker_rnaseq.sh
sed -i "s#PROTH_PATH#export PROTHINT_PATH=$protpath#" 00_scripts/06_braker_rnaseq.sh
sed -i "s#GMARK_PATH#export GENEMARK_PATH=$gmarkpath#" 00_scripts/06_braker_rnaseq.sh
sed -i "s#AUGCO_PATH#export AUGUSTUS_CONFIG_PATH=$augconf#" 00_scripts/06_braker_rnaseq.sh
sed -i "s#AUGBI_PATH#export AUGUSTUS_BIN_PATH=$augbin#" 00_scripts/06_braker_rnaseq.sh
sed -i "s#AUGSC_PATH#export AUGUSTUS_SCRIPTS_PATH=$augscripts#" 00_scripts/06_braker_rnaseq.sh

echo -e "---- running braker now on $haplotype ----- " 
echo "see details in braker_log in case of bugs" 
if [ ! -d $folder ]
then
	echo "creating $folder folder"
	mkdir $folder
fi


cd "$folder" || exit 1 #either haplo1 or haplo2

if [ ! -d 03_genome ]
then
	echo "creating 03_genome folder"
	mkdir 03_genome
fi

if [[ $annotateTE = "NO" ]]
then
    echo -e "\n-----------------------------------------------------"
    echo -e "\tNO TE annotation requested\n\twill skip this step"; 
    echo -e "\tassuming genome already softmasked ";
    echo -e "-----------------------------------------------------\n"
    if [ -L 03_genome/genome.wholemask.fa ] ; then
            echo "symbolic link alreadry present"
    else
        if ! ln -s "$genome" 03_genome/genome.wholemask.fa     
        then
            echo "error could not copy genome!"
            exit 1
        fi
    fi
fi
LOG_FOLDER="LOGS"
if [ ! -d "$LOG_FOLDER" ] ; then mkdir $LOG_FOLDER ; fi 

../00_scripts/06_braker_rnaseq.sh 03_genome/genome.wholemask.fa \
    "$haplotype" \
    $RNAseq \
    "$fungus" \
    "$bamlist" 2>&1 |tee LOGS/log_braker  #NO for no rnaseq  
