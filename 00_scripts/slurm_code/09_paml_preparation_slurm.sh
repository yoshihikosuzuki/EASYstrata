#!/bin/bash
##SBATCH --account=youraccount
#SBATCH --time=8:00:00
#SBATCH --job-name=paml
#SBATCH --output=log_paml-%J.out
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
# Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR
#source environnment 
mamba activate superannot/

# -- some colors for warnings in the terminal  --:
source config/colors
source config/config
############################################################
# verify info from the config file                         #
if [ -z "${haplotype1}" ] || [ -z "${haplotype2}" ] ; then
    echo "error missing info in config file, please check README.md"
    exit 2
fi
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
if [ ! -d 02_results ] ; then mkdir 02_results ; fi #ignore if already existent
if [[ ! -d 02_resumts/paml/ ]] ;then  
    mkdir -p 02_results/paml 
fi
#check if results already exists:
if [ "$ds_method" == "codeml" ] && [ -e 02_results/paml/results_codeml.txt ] && [ -e 02_results/paml/wanted_sequence ] ; then
    sizecodeml=$(wc -l 02_results/paml/results_codeml.txt |awk '{print $1}')
    sizeseq=$(wc -l 02_results/paml/wanted_sequence |awk '{print $1}')
    if [ "$sizeseq" == "$sizecodeml" ]; then 
        echo "warning results from codeml already exist"
        echo "not overwriting, please check the file" 
        echo "will continue analysis using extisting file"
        #exit 1 #or simply continue analysis 
    fi
elif [ "$ds_method" == "yn00" ] && [ -e 02_results/paml/results_YN.txt ] && [ -e 02_results/paml/wanted_sequence ] ; then
    sizecodeml=$(wc -l 02_results/paml/results_YN.txt |awk '{print $1}')
    sizeseq=$(wc -l 02_results/paml/wanted_sequence |awk '{print $1}')
    if [ "$sizeseq" == "$sizecodeml" ]; then 
        echo "warning results from codeml already exist"
        echo "not overwriting, please check the file" 
        echo "will continue analysis using extisting file"
        #exit 1
    fi
fi
#------------------------------ step 1 prepare input files  -------------------------------------#
#cds file:
cdsfile1=haplo1/08_best_run/$haplotype1.spliced_cds.fa
cdsfile2=haplo2/08_best_run/$haplotype2.spliced_cds.fa

##---  linearise the cds file ---#
awk '$0~/^>/{if(NR>1){print sequence;sequence=""}print $0}$0!~/^>/{sequence=sequence""$0}END{print sequence}' "${cdsfile2}"  > 02_results/paml/"$haplotype2".linearised.cds
awk '$0~/^>/{if(NR>1){print sequence;sequence=""}print $0}$0!~/^>/{sequence=sequence""$0}END{print sequence}' "${cdsfile1}"  > 02_results/paml/"$haplotype1".linearised.cds

##---- recover the wanted sequences in the CDS file #
cd 02_results/paml || exit 1

#all is run from paml folder now
if [ -e sorted."$haplotype1".wanted_cds.fa ] ; then rm sorted."$haplotype1".wanted_cds.fa ; fi
if [ -e sorted."$haplotype2".wanted_cds.fa ] ; then rm sorted."$haplotype2".wanted_cds.fa ; fi

while read -r pattern ; 
do 
    grep -w -A1 "$pattern" "$haplotype2".linearised.cds  >> sorted."$haplotype2".wanted_cds.fa ;
done < sco."$haplotype2".txt 
#
while read -r pattern ; 
do 
    grep -w -A1 "$pattern" "$haplotype1".linearised.cds >> sorted."$haplotype1".wanted_cds.fa ; 
done < sco."$haplotype1".txt 

### ---- then run paml and dnds.sh  ----- #
#------ part specific to PAML --------------------- #
fasta1=sorted."$haplotype1".wanted_cds.fa  
fasta2=sorted."$haplotype2".wanted_cds.fa

f1=$(basename "$fasta1" )
f2=$(basename "$fasta2" )
newf1=${f1%.fa**}.nostopcodon.fasta
newf2=${f2%.fa**}.nostopcodon.fasta

##1 ----- remove stop codon -------
awk '{if($1 !~ /^>/) {if( substr(toupper($0), length($0)-2, length($0)) ~ /(TGA|TAG|TAA)/){ print substr($0, 1, length($0)-3)} else  {print $0 }}else{print $0}}' "${fasta1}" > "${newf1}" 
awk '{if($1 !~ /^>/) {if( substr(toupper($0), length($0)-2, length($0)) ~ /(TGA|TAG|TAA)/){ print substr($0, 1, length($0)-3)} else  {print $0 }}else{print $0}}' "${fasta2}" > "${newf2}"

##2 ------ split and cat pairwise sequence -------
#split 
if [ -d sequence_files/ ] ; then rm -rf sequence_files ; fi

mkdir sequence_files

grep ">"  "$newf1" > ID1
grep ">"  "$newf2" > ID2

#check that all gene names are below 32 characters otherwise paml will fail!
awk '{print $0"\t"length }' ID1 |awk '$2>32 {print $1}' > long_geneID.hap1
awk '{print $0"\t"length }' ID2 |awk '$2>32 {print $1}' > long_geneID.hap2

#if some gene have length above we rename them using a structure of the type: gene$id
#id is a seq from 1 to n with n the max number of gene to rename
if [ -s long_geneID.hap1 ] ; then
    #the file is not empty; so we will rename the "long" genes:
    #0 - print some important warning as this may affect the user expectation: 
    nb_genes=$(wc -l long_geneID.hap1 |awk '{print $1}' )
    echo -e "${RED}!!! warning !!!\n some gene names are too long! ${NC} \n 
    a total of $nb_genes genes names will be renamed\n
    you'll find their name in the file:\n
    correspondance.table.hap1.txt\n\n" 
    #1 - create a correspondance table:
    #j=0 ; for i in $(cat long_geneID.hap1) ; do j=$(( "$j" + 1)) ; echo -e "$i\t>gene.$j"  >> correspondance.table.hap1.txt; done
    j=0 ; 
    while read -r line ; 
    do
        j=$(( "$j" + 1)) ; 
        echo -e "$line\t>hap1_gene.$j"  >> correspondance.table.hap1.txt; 
    done < long_geneID.hap1

    #2 - keep a copy:
    oldf1=$newf1.original.genename.fa
    cp "$newf1" "$oldf1"

    #3 - rename the fasta with awk
    awk 'NR==1 { next } FNR==NR { a[$1]=$2; next } $1 in a { $1=a[$1] }1' \
        correspondance.table.hap1.txt "${oldf1}" > "${newf1}" 

    #4 - re-extract the corrected ID for paml to succeed!
    grep ">"  "$newf1" > ID1
fi

#do the same for haplotype2: 
if [ -s long_geneID.hap2 ] ; then
    #the file is not empty; so we will rename the "long" genes:
    #0 - print some important warning as this may affect the user expectation: 
    nb_genes=$(wc -l long_geneID.hap2 |awk '{print $1}' )
    echo -e "${RED}!!! warning !!!\n some gene names are too long! ${NC} \n 
    a total of $nb_genes genes names will be renamed\n
    you'll find their name in the file:\n
    correspondance.table.hap2.txt\n\n" 
    #1 - create a correspondance table:
    j=0 
    while read -r line ; 
    do
        j=$(( "$j" + 1)) ; 
        echo -e "$line\t>hap2_gene.$j"  >> correspondance.table.hap2.txt; 
    done < long_geneID.hap2
 
    #2 - keep a copy:
    oldf2=$newf2.original.genename.fa
    cp "$newf2" "$oldf2"

    #3 - rename the fasta with awk
    awk 'NR==1 { next } FNR==NR { a[$1]=$2; next } $1 in a { $1=a[$1] }1' \
        correspondance.table.hap2.txt "${oldf2}" > "${newf2}"

    #4 - re-extract the corrected ID for paml to succeed!
    grep ">"  "$newf2" > ID2
fi
paste ID1 ID2 |sed 's/>//g' > wanted_sequence
mkdir 02_results/paml/sequence_files
