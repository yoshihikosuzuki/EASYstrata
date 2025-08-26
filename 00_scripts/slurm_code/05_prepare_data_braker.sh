#!/bin/bash
##SBATCH --account=youraccount
#SBATCH --time=74:00:00
#SBATCH --job-name=braker
#SBATCH --output=log_braker_array-%J.out
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
##SBATCH --array=1-5
# Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

#activate conda
source /local/env/envconda3.sh 
#activate the env
conda activate /scratch/qrougemont/new_superannot/

#give expected number of arguments from config file:
source config/config
source config/cpu_mem
##  ------------------------ general parameters --------------------------------  ##
while [ $# -gt 0 ] ; do
  case $1 in
    -g | --genome) genome="$2" ;echo "the genome file  is: $genome" >&2;;
    -s | --haplotype) haplotype="$2" ;echo "the haplotype name will be $haplotype" >&2;;
    -o | --output ) folder="$2"; echo "the data located in $folder will be processed">&2;;
    -f | --fungus ) fungus="$2" ; echo "Do species belong to fungi? $fungus" >&2;;
    -h | --help) echo -e "Option required:
    -g/--genome \t the reference genome file 
    -s/--haplotype\t the haplotype name (used for database building and basename in several steps)
    -o/--output\t the folder where the genome will be located (either haplo1 or haplo2). the genome will be copied here
    -f/--fungus\t a string YES/NO to state whether the species belong to fungus or not
        " >&2;exit 1;;
    esac
    shift
done
if [ -z "$genome" ] || [ -z "$haplotype" ] || [ -z "$folder" ] ; then
    echo -e  >&2 "Fatal error: Ref genome (-g), haplotype name (-s) and output folder (-o haplo1 or haplo2) not defined\n
    see manual with -h or --help"
     exit 2
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
sed -i "s#CDB_PATH#export CDBTOOLS_PATH=$cdbpath#" 00_scripts/06_braker_db_array.sh
sed -i "s#TSEBR_PATH#export TSEBRA_PATH=$tsebrapath#" 00_scripts/06_braker_db_array.sh
sed -i "s#PROTH_PATH#export PROTHINT_PATH=$protpath#" 00_scripts/06_braker_db_array.sh
sed -i "s#GMARK_PATH#export GENEMARK_PATH=$gmarkpath#" 00_scripts/06_braker_db_array.sh
sed -i "s#AUGCO_PATH#export AUGUSTUS_CONFIG_PATH=$augconf#" 00_scripts/06_braker_db_array.sh
sed -i "s#AUGBI_PATH#export AUGUSTUS_BIN_PATH=$augbin#" 00_scripts/06_braker_db_array.sh
sed -i "s#AUGSC_PATH#export AUGUSTUS_SCRIPTS_PATH=$augscripts#" 00_scripts/06_braker_db_array.sh

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
#----------------OrthoDB and Other Protein data -------------- #
target="$orthoDBspecies"

if [ -z "${orthoDBspecies}" ]; 
then
    echo "no orthoDB species provided"
        if [ -z "${RelatedProt}" ] ; 
        then
            echo -e "\nno related protein"
            echo -e "WARNING - you have not provided any external data\n"
            echo -e "braker will try running with RNAseq only\n\n" 
        else
            echo -e "related protein is $RelatedProt\n\n"
            if file --mime-type "$RelatedProt" | grep -q gzip$; then
              echo "$RelatedProt is gzipped"
              gunzip "$RelatedProt"
              relatProt=${RelatedProt%.gz}
           else
              echo "$RelatedProt is not gzipped"
              relatProt=$RelatedProt
           fi
        fi
else
    echo "orthoDBspecies is $target"
    clades=("Metazoa" 
        "Vertebrata" 
        "Viridiplantae" 
        "Arthropoda" 
        "Eukaryota" 
        "Fungi" 
        "Alveolata" 
        "Stramenopiles")
    if [[ ${clades[*]} =~ $target ]]
    then
        if [ -d odb12 ] ; then
            echo "folder odb12 already present"
        else 
           mkdir odb12
        fi
        cd odb12 || exit
        for file in "$target".fa*
        do 
            if [ -f "$file" ] 
            then
                echo "warning file $target.fa already present "
		gunzip "${target}".fa.gz
                cd ..
                if [ -s relatProt.fa ] ; then
		       	relatProt="relatProt.fa" ; 
		else
	            if [ -z ${RelatedProt} ] ; then
                       echo "no related protein"
                       echo "copying file"
                       cp odb12/"${target}".fa  relatProt.fa
                       relatProt="relatProt.fa"
                       echo "compressing file"
                       gzip odb12/"${target}".fa
                    else
                       echo "combining $RelatedProt with odb12 data" 
                       cat "$RelatedProt"  odb12/"${target}".fa > relatProt.fa
                       relatProt="relatProt.fa"
		    fi
		fi
            else
                echo "download partionned odb12 for $target lineage"
                wget -q https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/"${target}".fa.gz
                gunzip "${target}".fa.gz
                cd ../ 
                if [ -z ${RelatedProt} ] ; then
                    echo "no related protein"
                    echo "copying file"
                    cp odb12/"${target}".fa  relatProt.fa
                    relatProt="relatProt.fa"
                    echo "compressing file"
                    gzip odb12/"${target}".fa
                else
                    echo "combining $RelatedProt odb12 data" 
                    cat "$RelatedProt"  odb12/"${target}".fa > relatProt.fa
                    relatProt="relatProt.fa"
                fi
            fi
        done 
    else
        echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        echo -e "error the clade name you provided is not in the orthoDB list !!\n"
        echo -e "please check the clade name"
        echo -e "this should be one among (case sensitive):
        [Metazoa Vertebrata Viridiplantae Arthropoda Eukaryota Fungi Alveolata]"
        echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
        exit 1
    fi
fi

