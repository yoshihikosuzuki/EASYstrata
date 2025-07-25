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

cd 02_results/paml || exit 1

cat wanted_sequence | parallel -j 32 ./12_paml_parallel_version.sh {}

#===============================================================================
#then reshape a bit when all is good:
#===============================================================================
#we concatenate everyone to work with them in the next scripts:
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
#cleanup:
fasta1=sorted."$haplotype1".wanted_cds.fa  
fasta2=sorted."$haplotype2".wanted_cds.fa
f1=$(basename "$fasta1" )
f2=$(basename "$fasta2" )
newf1=${f1%.fa**}.nostopcodon.fasta
newf2=${f2%.fa**}.nostopcodon.fasta
rm "$f1" "$f2" "$newf1" "$newf2" "$haplotype2".linearised.cds "$haplotype1".linearised.cds ID1 ID2
