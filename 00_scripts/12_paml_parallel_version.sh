#!/bin/bash
#Author: QR
#date: 2024
#purpose: compute dS, dN, and their SE based on paml using yn00 model
source ../../config/config
source ../../config/colors
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
##---- recover the wanted sequences in the CDS file #
wantedseq=$1 #wanted_sequence 
line0=$(echo "$wantedseq" | awk '{print $1}' )
line1=$(echo "$wantedseq" | awk '{print $2}' ) #"$wantedseq" )

fasta1=sorted."$haplotype1".wanted_cds.fa  
fasta2=sorted."$haplotype2".wanted_cds.fa
f1=$(basename "$fasta1" )
f2=$(basename "$fasta2" )
newf1=${f1%.fa**}.nostopcodon.fasta
newf2=${f2%.fa**}.nostopcodon.fasta

mkdir sequence_files/tmp."${line0}".vs."${line1}"
grep -w -A1 "${line0}"  "$newf1" > sequence_files/tmp."${line0}".vs."${line1}"/sequence.fasta
grep -w -A1 "${line1}"  "$newf2" >> sequence_files/tmp."${line0}".vs."${line1}"/sequence.fasta

cp ../../config/yn00_template.ctl sequence_files/tmp."${line0}".vs."${line1}"/
cp ../../config/codeml.ctl sequence_files/tmp."${line0}".vs."${line1}"/

translatorx_vLocal.pl -i sequence_files/tmp."${line0}".vs."${line1}"/sequence.fasta \
            -o sequence_files/tmp."${line0}".vs."${line1}"/results 2>&1 |tee log.translator
   
#sed -i 's/sequence_NT.fasta/results.nt_ali.fasta/g' sequence_files/tmp."${line0}".vs."${line1}"/yn00_template.ctl
#sed -i 's/sequence_NT.fasta/results.nt_ali.fasta/g' sequence_files/tmp."${line0}".vs."${line1}"/codeml.ctl

cmd=$(command -v macse_v2.07.jar )
java -jar $cmd -prog alignSequences -seq sequence_files/tmp."${line0}".vs."${line1}"/sequence.fasta
        sed -i 's/!/-/g' sequence_files/tmp."${line0}".vs."${line1}"/sequence_NT.fasta 

cd sequence_files/tmp."${line0}".vs."${line1}"/ || exit 1
    
#run paml :
yn00 yn00_template.ctl
echo -ne | codeml codeml.ctl
    
cd ../../

#extract the dS, dN and SE from the output: 
awk '/\+\-/ && !/(dS|SE)/ {split(FILENAME, a, "."); 
    print $(NF-2), $(NF), $(NF-5), $(NF-3),"'${line0}'","'${line1}'"}'  \
        sequence_files/tmp."${line0}".vs."${line1}"/out_yn00_orthogp \
        >  sequence_files/tmp."${line0}".vs."${line1}"/resultat_Yang_Nielsen_2000_method.orthogp.txt 

sed -E '/dS =/ s/dS =/dS = /' sequence_files/tmp."${line0}".vs."${line1}"/mlc |
awk '/dS =/ {split(FILENAME, a, "."); 
    print $(NF), $(NF-3), $(NF-6),"'${line0}'","'${line1}'"}'  \
        >  sequence_files/tmp."${line0}".vs."${line1}"/resultat_codeml.txt
