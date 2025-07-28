#!/bin/bash
##SBATCH --account=youraccount
#SBATCH --time=24:00:00
#SBATCH --job-name=chgpts
#SBATCH --output=log_chgpt-%J.out
#SBATCH --mem=90G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
# Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR
#source environnment 
#source /local/env/envconda3.sh 
mamba activate superannot/

#Purpose:
#script to perform changepoint analyses
#Author: QR

source config/config
#------------------------ step 8 -- model comparison -------------------------------------------------#
if [[ -d 02_results/modelcomp ]]
then
    echo -e "WARNING directory modelcomp already exists! check its content first"
    echo -e "will perform the next analysis from existing files"
    exit 0
else
    mkdir 02_results/modelcomp/
    if [ -n "${ancestral_genome}" ] ; then
        Rscript 00_scripts/Rscripts/06.MCP_model_comp.R YES || \
        { echo -e "ERROR! changepoint failed - check your data\n" ; exit 1 ; }
    else
        Rscript 00_scripts/Rscripts/06.MCP_model_comp.R NO || \
         { echo -e "ERROR! changepoint failed - check your data\n" ; exit 1 ; }
    fi
 fi 
