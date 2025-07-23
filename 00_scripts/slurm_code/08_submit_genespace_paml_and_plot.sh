#!/bin/bash
##SBATCH --account=youraccount
#SBATCH --time=09:00:00
#SBATCH --job-name=geneSpace_and_co
#SBATCH --output=log_geneSpace_and_co-%J.out
#SBATCH --mem=80G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
# Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

#activate conda
#eval "$(conda shell.bash hook)"
#mamba activate superannot
#give expected number of arguments from config file:
source config/config

if [ -z "$haplotype1" ] ; then
    echo "error haplotype1 name not provided check your config file"
    exit 1
fi
if [ -z "$haplotype2" ] ; then
    echo "error haplotype2 name not provided check your config file"
    exit 1
fi
if [ -z "$scaffold" ] ; then
    echo "error scaffold names are not provided check your config file"
    echo "see details in README.md and example in example_data"
    exit 1
fi
if [ -n "${ancestral_genome}" ] ; then
    echo "ancestral_species is $ancestral_genome "
fi
if [ -n "${ancestral_gtf}" ] ; then
    echo "ancestral_species is $ancestral_genome "
fi


##  ------------------------ general parameters --------------------------------  ##
while [ $# -gt 0 ] ; do
  case $1 in
    -o | --option ) opt="$2"; echo "the option $opt will be used">&2;;
    -h | --help) echo -e "Option required:
        -o/--option\t the option to be used for processing the data.
        Option between 3 and 8 must be chosen.
        see README.md
        " >&2;exit 1;;
    esac
    shift
done

if [ -z "$opt" ] ; then
    echo -e  >&2 "Fatal error: and option (between 3 to 8 must be provided - see details in README.md\n
    see manual with -h or --help"
exit 2
fi

LOG_FOLDER="LOGS"
if [ ! -d "$LOG_FOLDER" ] ; then mkdir $LOG_FOLDER ; fi 

if [ -n "${ancestral_genome}" ] ; then
     echo "running analyses with ancestral genome"
    ./00_scripts/11_run_GeneSpace_paml_ideogram.sh \
         -s1 "$haplotype1" \
         -s2 "$haplotype2" \
         -a "$ancestral_genome" \
         -g "$ancestral_gtf" \
         -c "$scaffold" \
         -o "$opt" 2>&1 |\
        tee LOGS/log_GeneSpace_and_Co
else
     echo "running analyses without ancestral genome"
    ./00_scripts/11_run_GeneSpace_paml_ideogram.sh \
         -s1 "$haplotype1" \
         -s2 "$haplotype2" \
         -c "$scaffold" \
         -o "$opt" 2>&1 |\
        tee LOGS/log_GeneSpace_and_Co
fi
