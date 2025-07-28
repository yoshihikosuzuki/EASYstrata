#!/bin/bash  
##SBATCH --account=youraccount
#SBATCH --time=2:00:00
#SBATCH --job-name=paml
#SBATCH --output=log_paml-%J.out
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
# Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR
#source environnment 
#source /local/env/envconda3.sh 
mamba activate superannot/

source config/colors
source config/config
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

#check if results already exists:
if [ "$ds_method" == "codeml" ] && [ -e 02_results/paml/results_codeml.txt ] && [ -e 02_results/paml/wanted_sequence ] ; then
    sizecodeml=$(wc -l 02_results/paml/results_codeml.txt |awk '{print $1}')
    sizeseq=$(wc -l 02_results/paml/wanted_sequence |awk '{print $1}')
    if [ "$sizeseq" == "$sizecodeml" ]; then 
        echo "warning results from codeml already exist"
        echo "not overwriting, please check the file" 
	    echo "will continue analysis from existing file"
	    #exit 0
    else
       cd 02_results/paml || exit 1
       cat wanted_sequence | parallel -j 32 ../../00_scripts/12_paml_parallel_version.sh {}
       cd ../../
    fi
elif [ "$ds_method" == "yn00" ] && [ -e 02_results/paml/results_YN.txt ] && [ -e 02_results/paml/wanted_sequence ] ; then
    sizecodeml=$(wc -l 02_results/paml/results_YN.txt |awk '{print $1}')
    sizeseq=$(wc -l 02_results/paml/wanted_sequence |awk '{print $1}')
    if [ "$sizeseq" == "$sizecodeml" ]; then 
        echo "warning results from codeml already exist"
        echo "not overwriting, please check the file" 
	    echo "will continue analysis from existing file"
     	#exit 0
    else
       cd 02_results/paml || exit 1
       cat wanted_sequence | parallel -j 32 ../../00_scripts/12_paml_parallel_version.sh {}
       cd ../../
    fi
fi

#===============================================================================
#then reshape a bit when all is good:
#===============================================================================
#we concatenate everyone to work with them in the next scripts:
cd 02_results/paml || exit 1

if [ -e results_YN.txt ] ; then rm results_YN.txt ; fi
cat sequence_files/tmp.*/resultat_Yang_Nielsen_2000_method.orthogp.txt >> results_YN.txt
cp results_YN.txt results_YN.txt.bkp

if [ -e results_codeml.txt ] ; then rm results_codeml.txt ; fi

cat sequence_files/tmp.*/resultat_codeml.txt >> results_codeml.txt
cp results_codeml.txt results_codeml.txt.bkp

if [ -e correspondance.table.hap1.txt ] && [ ! -e correspondance.table.hap2.txt ] ; then
   sed -i 's/>//g' correspondance.table.hap1.txt
   awk 'NR==FNR{a[$2]=$1;next}$5 in a{$5=a[$5]}1' correspondance.table.hap1.txt results_YN.txt > tmp
   awk 'NR==FNR{a[$2]=$1;next}$4 in a{$4=a[$4]}1' correspondance.table.hap1.txt results_codeml.txt > tmp_cdml
   mv tmp results_YN.txt
   mv tmp_cdml results_codeml.txt
fi

if [ -e correspondance.table.hap2.txt ] && [ -e correspondance.table.hap1.txt ] ; then
   sed -i 's/>//g' correspondance.table.hap1.txt
   awk 'NR==FNR{a[$2]=$1;next}$5 in a{$5=a[$5]}1' correspondance.table.hap1.txt results_YN.txt > tmp
   awk 'NR==FNR{a[$2]=$1;next}$4 in a{$4=a[$4]}1' correspondance.table.hap1.txt results_codeml.txt > tmp_cdml

    sed -i 's/>//g' correspondance.table.hap2.txt
    awk 'NR==FNR{a[$2]=$1;next}$6 in a{$6=a[$6]}1' correspondance.table.hap2.txt tmp > results_YN.txt
    awk 'NR==FNR{a[$2]=$1;next}$5 in a{$5=a[$5]}1' correspondance.table.hap2.txt tmp_cdml > results_codeml.txt

fi

if [ -e correspondance.table.hap2.txt ] && [ ! -e correspondance.table.hap1.txt ]  ; then
   sed -i 's/>//g' correspondance.table.hap2.txt
   awk 'NR==FNR{a[$2]=$1;next}$6 in a{$6=a[$6]}1' correspondance.table.hap2.txt results_YN.txt > tmp
   awk 'NR==FNR{a[$2]=$1;next}$5 in a{$5=a[$5]}1' correspondance.table.hap2.txt results_codeml.txt > tmp_cdml
   mv tmp results_YN.txt
   mv tmp_cdml results_codeml.txt
fi
cd ../../

#===============================================================================
#conting:
if [ "$ds_method" == "yn00" ] ; then
    pamlsize=$(wc -l 02_results/paml/results_YN.txt |awk '{print $1}' ) 
else
    pamlsize=$(wc -l 02_results/paml/results_codeml.txt |awk '{print $1}' )
fi
scpo=$(wc -l 02_results/paml/single.copy.orthologs |awk '{print $1}' )
echo -e "there is $pamlsize results for PAML \n"
echo -e "there is $scpo single copy orthologs \n" 

#===============================================================================
#---------------------------- step 7 plot paml results  -----------------------#
bedhaplo1="haplo1/08_best_run/"$haplotype1".v2.bed"
bedhaplo2="haplo2/08_best_run/"$haplotype2".v2.bed"
bedanc="ancestral_sp/ancestral_sp.bed"

if [ ! -d Rlogs ]; then mkdir Rlogs ; fi
    if [ -n "${ancestral_genome}" ]; then

        echo "using ancestral genome"
        if ! Rscript ./00_scripts/Rscripts/03.plot_paml.R "$ds_method" "$max_ds" \
        "$haplotype1" "$haplotype2" "$scaffold" \
        "$bedhaplo1" "$bedhaplo2" ancestral_sp "$bedanc" 2> Rlogs/Rlogs_plot_paml
        then
            echo -e "\nERROR plotting paml results failed\n"
            echo -e "\nplease check input files and logs!!\n\n"
            exit 1
        fi
    else
        if ! Rscript ./00_scripts/Rscripts/03.plot_paml.R "$ds_method" "$max_ds" "$haplotype1" "$haplotype2" \
        "$scaffold"  "$bedhaplo1" "$bedhaplo2" 2> Rlogs/Rlogs_plot_paml
        then
            echo -e "\nERROR plotting paml results failed\n"
            echo -e "\nplease check input files and logs!!\n\n"
            exit 1
        fi
    fi
#===============================================================================
#cleanup:
fasta1=sorted."$haplotype1".wanted_cds.fa  
fasta2=sorted."$haplotype2".wanted_cds.fa
f1=$(basename "$fasta1" )
f2=$(basename "$fasta2" )
newf1=${f1%.fa**}.nostopcodon.fasta
newf2=${f2%.fa**}.nostopcodon.fasta
rm "$f1" "$f2" "$newf1" "$newf2" "$haplotype2".linearised.cds "$haplotype1".linearised.cds ID1 ID2
