#!/bin/bash
##SBATCH --account=youraccount
#SBATCH --time=8:00:00
#SBATCH --job-name=after_changepoint
#SBATCH --output=log_after_changepoint-%J.out
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
# Move to directory where job was submitted
cd $SLURM_SUBMIT_DIR

source /local/env/envconda3.sh 
conda activate /scratch/qrougemont/new_superannot/

#Purpose:
#master script to prepare bed files, laucnch GeneSpace, run paml and launch downstream Rscript 
#will check existence of all dependencies
#Date: 2025
#Author: QR
# -- some colors for warnings in the terminal  --:
source config/colors
source config/config
############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo -e "master script to: \n 1 - create bed files, \n 2 - create circos, \n 3 - create ideogram)"
   echo " "
   echo "Usage: $0 [-s1|-s2|-f|-a|-g|-h|]"
   echo "options:"
   echo " -h|--help: Print this Help."
   echo " -o|--options : the type of analysis to be performed: either 'synteny_and_Ds' (GeneSpace+Minimap2+Ds+changepoint), 'synteny_only' (GeneSpace+Minimap2), 'Ds_only' (paml and changepoint)"
   echo " "
   echo "dependancies: R"
}

############################################################
# Process the input options.                               #
############################################################
while [ $# -gt 0 ] ; do
  case $1 in
    -o  | --options ) options="$2" ; 
        echo -e "options for computation are ***${options}*** \n" >&2 ;;
    -h  | --help ) Help ; exit 2 ;;
   esac
   shift
done 

if [ -z "${haplotype1}" ] || [ -z "${haplotype2}" ] || [ -z "${scaffold}" ] || [ -z "${options}" ]  ; then
    Help
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
if [ -z ${scaffold+x} ]; then
        echo "ERROR! scaffold file is unset";
        echo "please provide a file of target scaffold (e.g. X or ancestral state)"
        echo "see example in example_data/scaffold.txt"
        exit 1
else
        echo "scaffold file is set to '$scaffold'";
fi

bedanc="ancestral_sp/ancestral_sp.bed"
    
##################################################################################################
#VARIABLE DECLARATION
bedhaplo1="haplo1/08_best_run/$haplotype1.v2.bed"
bedhaplo2="haplo2/08_best_run/$haplotype2.v2.bed"

scafforientation="02_results/chromosomes_orientation.txt"
chromosomes="02_results/chromosomes.txt"
#check if a TEfile exist for genome1 and genome2:
if [ -s haplo1/03_genome/filtered."$haplotype1".TE.bed ] ;
then
    genome1TE="haplo1/03_genome/filtered.$haplotype1.TE.bed"
    annotateTE="YES"
elif [ -n "$TEgenome1" ] ;
then
    genome1TE="$TEgenome1"
    annotateTE="YES"
else
    annotateTE="NO"
fi

if [ -s haplo2/03_genome/filtered."$haplotype2".TE.bed ] ;
then
    genome2TE="haplo2/03_genome/filtered.$haplotype2.TE.bed"
    annotateTE="YES"
elif [ -n "$TEgenome2" ] ;
then
    genome2TE="$TEgenome2"
    annotateTE="YES"
else
    annotateTE="NO"
fi

echo -e "annotateTE is set to $annotateTE"

if [ -n "$ancestral_genome" ] ; then
    ancestral=$(head -n1 ancestral_sp/ancestral_sp.fa.fai \
    |cut -f1 \
    |awk '{gsub("_","\t",$0) ; print $1}')

fi
# --------------- step8 data processing for ideogram and circos-----------------#
if [[ $options = "synteny_and_Ds" ]]  || [[ $options = "Ds_only" ]] || [[ $options = "synteny_only" ]] ; then 
    if [ -n "${ancestral_genome}" ]
    then
        echo "inferring synteny with ancestral species: "
        join  -1 4 -2 4 <(sort -k4,4 02_results/orthologues)  \
                        <(sort -k4,4 "$bedanc" ) \
            | sed 's/ /\t/g' \
            | join -1 5 -2 4 <(sort -k5,5 -) \
                           <(sort -k4,4 "$bedhaplo1" )  \
            |awk 'NR==1 {print "HOG\tOG\tN0\tchrom1\tGene1\tstart1\tend1\tchrom2\tGene2\tstart2\tend2"}
                    {print $3"\t"$4"\t"$5"\t"$7"\t"$2"\t"$8"\t"$9"\t"$10"\t"$1"\t"$11"\t"$12}' \
            > 02_results/synteny_ancestral_sp_"$haplotype1".txt
       
        if [ -s 02_results/synteny_ancestral_sp_"$haplotype1".txt ]
        then
            size1=$(wc -l   02_results/synteny_ancestral_sp_"$haplotype1".txt |awk '{print $1}')
        else
            echo "synteny file between ancestral species and $haplotype1 is empty"
            echo "please check your data"
            exit 1
        fi
    
        join  -1 4 -2 4 <(sort -k4,4 02_results/orthologues)  \
                <(sort -k4,4 "$bedanc" ) \
                | sed 's/ /\t/g' \
                |join -1 6 -2 4 <(sort -k6,6 -) \
                                <(sort -k4,4 "$bedhaplo2" ) \
                |awk 'NR==1 {print "HOG\tOG\tN0\tchrom1\tGene1\tstart1\tend1\tchrom2\tGene2\tstart2\tend2"}
                    {print $3"\t"$4"\t"$5"\t"$7"\t"$2"\t"$8"\t"$9"\t"$10"\t"$1"\t"$11"\t"$12}' \
                >  02_results/synteny_ancestral_sp_"$haplotype2".txt
    
         if [ -s 02_results/synteny_ancestral_sp_"$haplotype2".txt ]
         then
             size2=$(wc -l   02_results/synteny_ancestral_sp_"$haplotype2".txt |awk '{print $1}')
         else
             echo "synteny file between ancestral species and $haplotype2 is empty"
             echo "please check your data"
             exit 1
         fi
       
           echo -e "number of lines in synteny file ancestral_sp vs $haplotype1 is $size1"
           echo -e "number of lines in synteny file ancestral_sp vs $haplotype2 is $size2"
     
        join  -1 5 -2 4 <(sort -k5,5 02_results/orthologues)  \
                        <(sort -k4,4 "$bedhaplo1" ) \
            | sed 's/ /\t/g' \
            | join -1 6 -2 4 <(sort -k6,6 -) \
                           <(sort -k4,4 "$bed$haplotype2" )  \
            |awk 'NR==1 {print "HOG\tOG\tN0\tchrom1\tGene1\tstart1\tend1\tchrom2\tGene2\tstart2\tend2"}
                    {print $3"\t"$4"\t"$5"\t"$7"\t"$2"\t"$8"\t"$9"\t"$10"\t"$1"\t"$11"\t"$12}' \
            > 02_results/synteny_"$haplotype1"_"$haplotype2".txt
    
         if [ -s 02_results/synteny_"$haplotype1"_"$haplotype2".txt ] ;
         then
             size3=$(wc -l   02_results/synteny_"$haplotype1"_"$haplotype2".txt |awk '{print $1}')
         else
             echo "synteny file between $haplotype1 and $haplotype2 is empty"
             echo "please check your data"
             exit 1
         fi
                   
         echo -e "number of lines in synteny file $haplotype1 vs $haplotype2 is $size3"

     else
        echo "no ancestral species assumed "
        echo "inferring synteny between $haplotype1 and $haplotype2"


       join  -1 4 -2 4 <(sort -k4,4 02_results/orthologues)  \
                       <(sort -k4,4 "$bedhaplo1" ) \
            | sed 's/ /\t/g' \
            | join -1 5 -2 4 <(sort -k5,5 -) \
                           <(sort -k4,4 "$bedhaplo2" )  \
            |awk 'NR==1 {print "HOG\tOG\tN0\tchrom1\tGene1\tstart1\tend1\tchrom2\tGene2\tstart2\tend2"}
                    {print $3"\t"$4"\t"$5"\t"$6"\t"$2"\t"$7"\t"$8"\t"$9"\t"$1"\t"$10"\t"$11}'   \
            > 02_results/synteny_"$haplotype1"_"$haplotype2".txt

         if [ -s 02_results/synteny_"$haplotype1"_"$haplotype2".txt ] ;
         then
             size3=$(wc -l 02_results/synteny_"$haplotype1"_"$haplotype2".txt |awk '{print $1}')
         else
             echo "synteny file between $haplotype1 and $haplotype2 is empty"
             echo "please check your data"
             exit 1
         fi
                   
         echo -e "number of lines in synteny file $haplotype1 vs $haplotype2 is $size3"

    fi

    #/!\ chromosomes should be reconstructed on the fly from the N0.tsv file
    if [ -n "${ancestral_genome}" ]
    then
        cat <(sed 1d 02_results/synteny_"$haplotype1"_"$haplotype2".txt \
             |awk -v var1="$haplotype1" -v var2="$haplotype2"  '{print var1"\t"$4"\n"var2"\t"$8}' \
             |sort |uniq -c |awk '$1>0 {print $2"\t"$3}'  ) \
             <(sed 1d 02_results/synteny_ancestral_sp_"$haplotype1".txt \
             |awk -v var1="$ancestral" -v var2="$haplotype1"  '{print var1"\t"$4"\n"var2"\t"$8}' \
             |sort |uniq -c |awk '$1>0 {print $2"\t"$3}' ) \

    else
           sed 1d 02_results/synteny_"$haplotype1"_"$haplotype2".txt \
             |awk -v var1="$haplotype1" -v var2="$haplotype2"  '{print var1"\t"$4"\n"var2"\t"$8}' \
             |sort |uniq -c |awk '$1>0 {print $2"\t"$3}' \
             |sort |uniq > 02_results/chromosomes.txt

    awk '{print $1"\t"$5"\t"$9}' 02_results/synteny_"$haplotype1"_"$haplotype2".txt |sed 1d > 02_results/sco


    fi
    chromosomes="02_results/chromosomes.txt"
    awk '{print $0"\tN"}' $chromosomes > 02_results/chromosomes_orientation.txt
    scafforientation="02_results/chromosomes_orientation.txt"

    #test if previous step was successfull else plot or exit with high levels of pain
    #take advantage of samtools to get length of genome
    samtools faidx haplo1/03_genome/"$haplotype1".fa 
    samtools faidx haplo2/03_genome/"$haplotype2".fa
fi
########################################################################################################
#	Data PLOTTING AFTER CHANGEPOINT
#######################################################################################################
if [[ $options = "synteny_and_Ds" ]]  || [[ $options = "Ds_only" ]] || [[ $options = "plots" ]] ; then
echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo -e "~ \tcreating ideogram colored by strata\t ~"
echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

for links in 02_results/modelcomp/noprior/classif.s*haplo1.haplo2 ; 
do 
    if [ -n "$scafforientation" ] ; then
        echo "particular orientation provided"
        Rscript ./00_scripts/Rscripts/04.ideogram.R \
                -c 02_results/sco \
                -i "$bedhaplo1" \
                -j "$bedhaplo2" \
                -f haplo1/03_genome/"$haplotype1".fa.fai \
                -g haplo2/03_genome/"$haplotype2".fa.fai \
                -l "$links" \
                -s "$scafforientation" 2> Rlogs/Rlogs_plot_ideogram_colored_"$(basename "$links")".txt 
    else
        Rscript ./00_scripts/Rscripts/04.ideogram.R \
                -c 02_results/sco \
                -i "$bedhaplo1" \
                -j "$bedhaplo2" \
                -f haplo1/03_genome/"$haplotype1".fa.fai \
                -g haplo2/03_genome/"$haplotype2".fa.fai \
                -l "$links" 2> Rlogs/Rlogs_plot_ideogram_colored_"$(basename "$links")".txt 
    fi
done

if [ -n "${ancestral_genome}" ] ; then
   for links in 02_results/modelcomp/noprior/classif.s*ancestral.haplo1 ; 
   do 
      if [ -n "$scafforientation" ] ; then
        echo "particular orientation provided"
        Rscript ./00_scripts/Rscripts/04.ideogram.R \
                -c 02_results/sco_anc \
                -i "$bedanc" \
                -j "$bedhaplo1"  \
                -f ancestral_sp/ancestral_sp.fa.fai \
                -g haplo1/03_genome/"$haplotype1".fa.fai \
                -l "$links" \
                -s "$scafforientation"  2> Rlogs/Rlogs_plot_ideogram_colored_"$(basename "$links")"_ancestral.txt
     else
         Rscript ./00_scripts/Rscripts/04.ideogram.R \
                -c 02_results/sco_anc \
                -i "$bedanc" \
                -j "$bedhaplo1"  \
                -f ancestral_sp/ancestral_sp.fa.fai \
                -g haplo1/03_genome/"$haplotype1".fa.fai \
                -l "$links" 2> Rlogs/Rlogs_plot_ideogram_colored_"$(basename "$links")".txt 
     fi
   done
fi
## for fun we now make circos plot with links consisting of the strata and colored by their ds Values:
#preparer des bed file pour faire des circos plots:
if [ ! -d 02_results/bed ] ; then mkdir 02_results/bed ; fi
if [ -n "${ancestral_genome}" ] ;
then
    cut -f 1-3,19 02_results/modelcomp/noprior/df.txt |sed 1d > 02_results/bed/ancestralspecies.3strata.bed
    cut -f 1-3,20 02_results/modelcomp/noprior/df.txt |sed 1d > 02_results/bed/ancestralspecies.4strata.bed
    cut -f 1-3,21 02_results/modelcomp/noprior/df.txt |sed 1d > 02_results/bed/ancestralspecies.5strata.bed
    cut -f 1-3,22 02_results/modelcomp/noprior/df.txt |sed 1d > 02_results/bed/ancestralspecies.6strata.bed
    cut -f 1-3,23 02_results/modelcomp/noprior/df.txt |sed 1d > 02_results/bed/ancestralspecies.7strata.bed
    cut -f 1-3,24 02_results/modelcomp/noprior/df.txt |sed 1d > 02_results/bed/ancestralspecies.8strata.bed
    cut -f 1-3,25 02_results/modelcomp/noprior/df.txt |sed 1d > 02_results/bed/ancestralspecies.9strata.bed

    #haplotype1 bed:
    join -1 1 -2 4 <(sort -k1,1 <(cut -f7,19 02_results/modelcomp/noprior/df.txt )) \
         <(sort -k4,4 "$bedhaplo1" )  \
         |awk '{print $3"\t"$4"\t"$4"\t"$2}'  > 02_results/bed/"$haplotype1".3strata.bed
    join -1 1 -2 4 <(sort -k1,1 <(cut -f7,20 02_results/modelcomp/noprior/df.txt )) \
        <(sort -k4,4 "$bedhaplo1" )  \
        |awk '{print $3"\t"$4"\t"$4"\t"$2}'  > 02_results/bed/"$haplotype1".4strata.bed
    join -1 1 -2 4 <(sort -k1,1 <(cut -f7,21 02_results/modelcomp/noprior/df.txt )) \
        <(sort -k4,4 "$bedhaplo1" )  \
        |awk '{print $3"\t"$4"\t"$4"\t"$2}'  > 02_results/bed/"$haplotype1".5strata.bed
    join -1 1 -2 4 <(sort -k1,1 <(cut -f7,22 02_results/modelcomp/noprior/df.txt )) \
        <(sort -k4,4 "$bedhaplo1" )  \
        |awk '{print $3"\t"$4"\t"$4"\t"$2}'  > 02_results/bed/"$haplotype1".6strata.bed
    join -1 1 -2 4 <(sort -k1,1 <(cut -f7,23 02_results/modelcomp/noprior/df.txt )) \
        <(sort -k4,4 "$bedhaplo1" )  \
        |awk '{print $3"\t"$4"\t"$4"\t"$2}'  > 02_results/bed/"$haplotype1".7strata.bed
    join -1 1 -2 4 <(sort -k1,1 <(cut -f7,24 02_results/modelcomp/noprior/df.txt )) \
        <(sort -k4,4 "$bedhaplo1" )  \
        |awk '{print $3"\t"$4"\t"$4"\t"$2}'  > 02_results/bed/"$haplotype1".8strata.bed
    join -1 1 -2 4 <(sort -k1,1 <(cut -f7,25 02_results/modelcomp/noprior/df.txt )) \
        <(sort -k4,4 "$bedhaplo1" )  \
        |awk '{print $3"\t"$4"\t"$4"\t"$2}'  > 02_results/bed/"$haplotype1".9strata.bed

    #haplotype2 bed:
    join -1 1 -2 4 <(sort -k1,1 <(cut -f12,19 02_results/modelcomp/noprior/df.txt )) \
        <(sort -k4,4 "$bedhaplo2" )  \
        |awk '{print $3"\t"$4"\t"$4"\t"$2}'  > 02_results/bed/"$haplotype2".3strata.bed
    join -1 1 -2 4 <(sort -k1,1 <(cut -f12,20 02_results/modelcomp/noprior/df.txt )) \
        <(sort -k4,4 "$bedhaplo2" )  \
        |awk '{print $3"\t"$4"\t"$4"\t"$2}'  > 02_results/bed/"$haplotype2".4strata.bed
    join -1 1 -2 4 <(sort -k1,1 <(cut -f12,21 02_results/modelcomp/noprior/df.txt )) \
        <(sort -k4,4 "$bedhaplo2" )  \
        |awk '{print $3"\t"$4"\t"$4"\t"$2}'  > 02_results/bed/"$haplotype2".5strata.bed
    join -1 1 -2 4 <(sort -k1,1 <(cut -f12,22 02_results/modelcomp/noprior/df.txt )) \
        <(sort -k4,4 "$bedhaplo2" )  \
        |awk '{print $3"\t"$4"\t"$4"\t"$2}'  > 02_results/bed/"$haplotype2".6strata.bed
    join -1 1 -2 4 <(sort -k1,1 <(cut -f12,23 02_results/modelcomp/noprior/df.txt )) \
        <(sort -k4,4 "$bedhaplo2" )  \
        |awk '{print $3"\t"$4"\t"$4"\t"$2}'  > 02_results/bed/"$haplotype2".7strata.bed
    join -1 1 -2 4 <(sort -k1,1 <(cut -f12,24 02_results/modelcomp/noprior/df.txt )) \
        <(sort -k4,4 "$bedhaplo2" )  \
        |awk '{print $3"\t"$4"\t"$4"\t"$2}'  > 02_results/bed/"$haplotype2".8strata.bed
    join -1 1 -2 4 <(sort -k1,1 <(cut -f12,25 02_results/modelcomp/noprior/df.txt )) \
        <(sort -k4,4 "$bedhaplo2" )  \
        |awk '{print $3"\t"$4"\t"$4"\t"$2}'  > 02_results/bed/"$haplotype2".9strata.bed

    cat 02_results/bed/ancestralspecies.3strata.bed 02_results/bed/"$haplotype1".3strata.bed > 02_results/bed/ancestralspecies_"$haplotype1".3strata.bed
    cat 02_results/bed/ancestralspecies.4strata.bed 02_results/bed/"$haplotype1".4strata.bed > 02_results/bed/ancestralspecies_"$haplotype1".4strata.bed
    cat 02_results/bed/ancestralspecies.5strata.bed 02_results/bed/"$haplotype1".5strata.bed > 02_results/bed/ancestralspecies_"$haplotype1".5strata.bed
    cat 02_results/bed/ancestralspecies.6strata.bed 02_results/bed/"$haplotype1".6strata.bed > 02_results/bed/ancestralspecies_"$haplotype1".6strata.bed
    cat 02_results/bed/ancestralspecies.7strata.bed 02_results/bed/"$haplotype1".7strata.bed > 02_results/bed/ancestralspecies_"$haplotype1".7strata.bed
    cat 02_results/bed/ancestralspecies.8strata.bed 02_results/bed/"$haplotype1".8strata.bed > 02_results/bed/ancestralspecies_"$haplotype1".8strata.bed
    cat 02_results/bed/ancestralspecies.9strata.bed 02_results/bed/"$haplotype1".9strata.bed > 02_results/bed/ancestralspecies_"$haplotype1".9strata.bed

    cat 02_results/bed/ancestralspecies.3strata.bed 02_results/bed/"$haplotype2".3strata.bed > 02_results/bed/ancestralspecies_"$haplotype2".3strata.bed
    cat 02_results/bed/ancestralspecies.4strata.bed 02_results/bed/"$haplotype2".4strata.bed > 02_results/bed/ancestralspecies_"$haplotype2".4strata.bed
    cat 02_results/bed/ancestralspecies.5strata.bed 02_results/bed/"$haplotype2".5strata.bed > 02_results/bed/ancestralspecies_"$haplotype2".5strata.bed
    cat 02_results/bed/ancestralspecies.6strata.bed 02_results/bed/"$haplotype2".6strata.bed > 02_results/bed/ancestralspecies_"$haplotype2".6strata.bed
    cat 02_results/bed/ancestralspecies.7strata.bed 02_results/bed/"$haplotype2".7strata.bed > 02_results/bed/ancestralspecies_"$haplotype2".7strata.bed
    cat 02_results/bed/ancestralspecies.8strata.bed 02_results/bed/"$haplotype2".8strata.bed > 02_results/bed/ancestralspecies_"$haplotype2".8strata.bed
    cat 02_results/bed/ancestralspecies.9strata.bed 02_results/bed/"$haplotype2".9strata.bed > 02_results/bed/ancestralspecies_"$haplotype2".9strata.bed
else
    #tester si pas d'ancestral
    #haplotype1 bed:
    cut -f1-3,15 02_results/modelcomp/noprior/df.txt |sed 1d > 02_results/bed/"$haplotype1".2strata.bed
    cut -f1-3,16 02_results/modelcomp/noprior/df.txt |sed 1d > 02_results/bed/"$haplotype1".3strata.bed
    cut -f1-3,17 02_results/modelcomp/noprior/df.txt |sed 1d > 02_results/bed/"$haplotype1".4strata.bed
    cut -f1-3,18 02_results/modelcomp/noprior/df.txt |sed 1d > 02_results/bed/"$haplotype1".5strata.bed
    cut -f1-3,19 02_results/modelcomp/noprior/df.txt |sed 1d > 02_results/bed/"$haplotype1".6strata.bed
    cut -f1-3,20 02_results/modelcomp/noprior/df.txt |sed 1d > 02_results/bed/"$haplotype1".7strata.bed
    cut -f1-3,21 02_results/modelcomp/noprior/df.txt |sed 1d > 02_results/bed/"$haplotype1".8strata.bed
    cut -f1-3,22 02_results/modelcomp/noprior/df.txt |sed 1d > 02_results/bed/"$haplotype1".9strata.bed
    cut -f1-3,23 02_results/modelcomp/noprior/df.txt |sed 1d > 02_results/bed/"$haplotype1".10strata.bed
   #haplotype2 bed:
    join -1 1 -2 4 <(sort -k1,1 <(cut -f10,15 02_results/modelcomp/noprior/df.txt )) \
       <(sort -k4,4 "$bedhaplo2" )  \
       |awk '{print $3"\t"$4"\t"$4"\t"$2}'  > 02_results/bed/"$haplotype2".2strata.bed
    join -1 1 -2 4 <(sort -k1,1 <(cut -f10,16 02_results/modelcomp/noprior/df.txt )) \
       <(sort -k4,4 "$bedhaplo2" )  \
       |awk '{print $3"\t"$4"\t"$4"\t"$2}'  > 02_results/bed/"$haplotype2".3strata.bed
    join -1 1 -2 4 <(sort -k1,1 <(cut -f10,17 02_results/modelcomp/noprior/df.txt )) \
       <(sort -k4,4 "$bedhaplo2" )  \
       |awk '{print $3"\t"$4"\t"$4"\t"$2}'  > 02_results/bed/"$haplotype2".4strata.bed
    join -1 1 -2 4 <(sort -k1,1 <(cut -f10,18 02_results/modelcomp/noprior/df.txt )) \
        <(sort -k4,4 "$bedhaplo2" )  \
        |awk '{print $3"\t"$4"\t"$4"\t"$2}'  > 02_results/bed/"$haplotype2".5strata.bed
    join -1 1 -2 4 <(sort -k1,1 <(cut -f10,19 02_results/modelcomp/noprior/df.txt )) \
        <(sort -k4,4 "$bedhaplo2" )  \
        |awk '{print $3"\t"$4"\t"$4"\t"$2}'  > 02_results/bed/"$haplotype2".6strata.bed
    join -1 1 -2 4 <(sort -k1,1 <(cut -f10,20 02_results/modelcomp/noprior/df.txt )) \
        <(sort -k4,4 "$bedhaplo2" )  \
        |awk '{print $3"\t"$4"\t"$4"\t"$2}'  > 02_results/bed/"$haplotype2".7strata.bed
    join -1 1 -2 4 <(sort -k1,1 <(cut -f10,21 02_results/modelcomp/noprior/df.txt )) \
        <(sort -k4,4 "$bedhaplo2" )  \
        |awk '{print $3"\t"$4"\t"$4"\t"$2}'  > 02_results/bed/"$haplotype2".8strata.bed
    join -1 1 -2 4 <(sort -k1,1 <(cut -f10,22 02_results/modelcomp/noprior/df.txt )) \
        <(sort -k4,4 "$bedhaplo2" )  \
        |awk '{print $3"\t"$4"\t"$4"\t"$2}'  > 02_results/bed/"$haplotype2".9strata.bed
    join -1 1 -2 4 <(sort -k1,1 <(cut -f10,23 02_results/modelcomp/noprior/df.txt )) \
        <(sort -k4,4 "$bedhaplo2" )  \
        |awk '{print $3"\t"$4"\t"$4"\t"$2}'  > 02_results/bed/"$haplotype2".10strata.bed
fi
cat 02_results/bed/"$haplotype1".2strata.bed 02_results/bed/"$haplotype2".2strata.bed> 02_results/bed/"$haplotype1"."$haplotype2".2strata.bed
cat 02_results/bed/"$haplotype1".3strata.bed 02_results/bed/"$haplotype2".3strata.bed> 02_results/bed/"$haplotype1"."$haplotype2".3strata.bed
cat 02_results/bed/"$haplotype1".4strata.bed 02_results/bed/"$haplotype2".4strata.bed> 02_results/bed/"$haplotype1"."$haplotype2".4strata.bed
cat 02_results/bed/"$haplotype1".5strata.bed 02_results/bed/"$haplotype2".5strata.bed> 02_results/bed/"$haplotype1"."$haplotype2".5strata.bed
cat 02_results/bed/"$haplotype1".6strata.bed 02_results/bed/"$haplotype2".6strata.bed> 02_results/bed/"$haplotype1"."$haplotype2".6strata.bed
cat 02_results/bed/"$haplotype1".7strata.bed 02_results/bed/"$haplotype2".7strata.bed> 02_results/bed/"$haplotype1"."$haplotype2".7strata.bed
cat 02_results/bed/"$haplotype1".8strata.bed 02_results/bed/"$haplotype2".8strata.bed> 02_results/bed/"$haplotype1"."$haplotype2".8strata.bed
cat 02_results/bed/"$haplotype1".9strata.bed 02_results/bed/"$haplotype2".9strata.bed> 02_results/bed/"$haplotype1"."$haplotype2".9strata.bed
cat 02_results/bed/"$haplotype1".10strata.bed 02_results/bed/"$haplotype2".10strata.bed> 02_results/bed/"$haplotype1"."$haplotype2".10strata.bed

echo -e "\n\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo -e "~ \tcreating circos colored by strata\t ~"
echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"

if [ -n "${ancestral_genome}" ] ; then
    echo -e "\n ancestral genome provided \n"
    for links in 02_results/bed/ancestralspecies"$haplotype1".*.strata.bed ; do
        if [[ $annotateTE = "YES" ]] ; then
            echo -e "\nTE bed provided\n"
            Rscript 00_scripts/Rscripts/05_plot_circos.R -s "$ancestral" -p "$haplotype1" \
               -c "$chromosomes" \
               -y 02_results/synteny_ancestral_sp_"$haplotype1".txt \
               -f ancestral_sp/ancestral_sp.fa.fai \
               -g haplo1/03_genome/"$haplotype1".fa.fai \
               -i "$bedanc"  \
               -j "$bedhaplo1"  \
               -t "$genome1TE" \
               -u "$genome2TE" \
               -l "$links" 2> Rlogs/Rlogs_plot_circos_TE_ancestral_sp_colored_"$(basename "$links")".txt
        else 
            echo -e "\nassuming noTE\n"
            Rscript 00_scripts/Rscripts/05_plot_circos.R -s "$ancestral" -p "$haplotype1" \
               -c "$chromosomes" \
               -y 02_results/synteny_ancestral_sp_"$haplotype1".txt \
               -f ancestral_sp/ancestral_sp.fa.fai \
               -g haplo1/03_genome/"$haplotype1".fa.fai \
               -i "$bedanc"  \
               -j "$bedhaplo1"  \
               -t "$genome1TE" \
               -u "$genome2TE" \
               -l "$links" 2> Rlogs/Rlogs_plot_circos_noTE_colored_"$(basename "$links")".txt
        fi
    done
else
    echo -e "\nno ancestral genome provided\n"
fi
#performing haplo1 vs haplo2 comparisons :
for links in 02_results/bed/"$haplotype1"."$haplotype2".*strata.bed ; do
    if [[ $annotateTE = "YES" ]] ; then
        echo -e "\nTE bed provided\n"02_results/dS.values.forchangepoint.txt
        echo -e "\runngin circos with links file $links\n\n"
        Rscript 00_scripts/Rscripts/05_plot_circos.R -s "$haplotype1" -p "$haplotype2" \
           -c "$chromosomes" \
           -y 02_results/synteny_"$haplotype1"_"$haplotype2".txt \
           -f haplo1/03_genome/"$haplotype1".fa.fai\
           -g haplo2/03_genome/"$haplotype2".fa.fai \
           -i "$bedhaplo1"  \
           -j "$bedhaplo2"  \
           -t "$genome1TE" \
           -u "$genome2TE" \
           -l "$links" 2> Rlogs/Rlogs_plot_circos_TE_colored_"$(basename "$links")".txt
    else 
        echo assuming noTE
        Rscript 00_scripts/Rscripts/05_plot_circos.R \
           -s "$haplotype1"  \
           -p "$haplotype2" \
           -c "$chromosomes" \
           -y 02_results/synteny_"$haplotype1"_"$haplotype2".txt \
           -f haplo1/03_genome/"$haplotype1".fa.fai\
           -g haplo2/03_genome/"$haplotype2".fa.fai \
           -i "$bedhaplo1"  \
           -j "$bedhaplo2"  \
           -l "$links" 2> Rlogs/Rlogs_plot_ciros_noTE_colored_"$(basename "$links")".txt
    fi
done

#now we will do the same discretisation of dS for plotting in ideogram: 
echo -e "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
echo -e "~ \tcreating ideogram colored by dS values\t ~"
echo -e "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"

if [  -n "${ancestral_genome}" ] ; then
    echo -e "\nancestral genome was provided for inference\n" 
    if [ -n "$scafforientation" ] ; then
        echo -e "\nparticular orientation will be used\n"
        #we will make an ideogram with it 
        if ! Rscript ./00_scripts/Rscripts/04.ideogram.R \
                -c 02_results/sco_anc \
                -i "$bedanc" \
                -j "$bedhaplo1"  \
                -d 02_results/dS.values.forchangepoint.txt \
                -f ancestral_sp/ancestral_sp.fa.fai \
                -g haplo1/03_genome/"$haplotype1".fa.fai \
                -s "$scafforientation" 2> Rlogs/Rlogs_plot_ideogram_colored_by_dsquantile.txt 
        then
             echo -e "\nERROR: ideograms failed /!\ \n
             please check logs and input data\n" 
             exit 1
        fi
    else
        if ! Rscript ./00_scripts/Rscripts/04.ideogram.R \
                -c 02_results/sco_anc \
                -i "$bedanc" \
                -j "$bedhaplo1"  \
                -d 02_results/dS.values.forchangepoint.txt \
                -f ancestral_sp/ancestral_sp.fa.fai \
                -g haplo1/03_genome/"$haplotype1".fa.fai \
                -s "$scafforientation" 2> Rlogs/Rlogs_plot_ideogram_colored_by_dsquantile.txt 
        then
             echo -e "\nERROR: ideograms failed /!\ \n
             please check logs and input data\n" 
             exit 1
        fi
    fi
fi

if [ -n "$scafforientation" ] ; then
    echo -e "particular orientation will be used"

    if ! Rscript ./00_scripts/Rscripts/04.ideogram.R \
        -c 02_results/sco \
        -i "$bedhaplo1" \
        -j "$bedhaplo2"  \
        -d 02_results/dS.values.forchangepoint.txt \
        -f haplo1/03_genome/"$haplotype1".fa.fai \
        -g haplo2/03_genome/"$haplotype2".fa.fai \
        -s "$scafforientation" 2> Rlogs/Rlogs_plot_ideogram_colored_by_dsquantile.txt 
    then
        echo -e "\nERROR: ideograms failed /!\ \n
        please check logs and input data\n" 
        exit 1
    fi
else
    if ! Rscript ./00_scripts/Rscripts/04.ideogram.R \
        -c 02_results/sco \
        -i "$bedhaplo1" \
        -j "$bedhaplo2"  \
        -d 02_results/dS.values.forchangepoint.txt \
        -f haplo1/03_genome/"$haplotype1".fa.fai \
        -g haplo2/03_genome/"$haplotype2".fa.fai 2> Rlogs/Rlogs_plot_ideogram_colored_by_dsquantile.txt
    then
        echo -e "\nERROR: ideograms failed /!\ \n
        please check logs and input data\n" 
        exit 1
    fi
fi
#finally create circos plot based on dS values quantiles:
if [[ $annotateTE = "YES" ]] ; then
                echo "assuming TE bed file exist"
                echo "assuming links"
                #bed file of TE should exist:
                if ! Rscript 00_scripts/Rscripts/05_plot_circos.R -s "$haplotype1" -p "$haplotype2" \
                    -c "$chromosomes" \
                    -y 02_results/synteny_"$haplotype1"_"$haplotype2".txt  \
                    -f haplo1/03_genome/"$haplotype1".fa.fai \
                    -g haplo2/03_genome/"$haplotype2".fa.fai \
                    -i "$bedhaplo1"  \
                    -j "$bedhaplo2"  \
                    -t "$genome1TE" \
                    -u "$genome2TE" \
                    -d 02_results/dS.values.forchangepoint.txt 2> Rlogs/Rlogs_plot_circos_color_by_dsquantile.txt
                then
                    echo -e "\nERROR: circos plots failed /!\ \n
                    please check logs and input data\n" 
                    exit 1
                fi
else #assume no TE: 
                echo "assuming no TE bed files"
                echo "assuming links"
                if ! Rscript 00_scripts/Rscripts/05_plot_circos.R -s "$haplotype1" -p "$haplotype2" \
                    -c "$chromosomes" \
                    -y 02_results/synteny_"$haplotype1"_"$haplotype2".txt  \
                    -f haplo1/03_genome/"$haplotype1".fa.fai \
                    -g haplo2/03_genome/"$haplotype2".fa.fai \
                    -i "$bedhaplo1"  \
                    -j "$bedhaplo2"  \
                    -d 02_results/dS.values.forchangepoint.txt 2> Rlogs/Rlogs_plot_circos_colored_by_dsquantile.txt
                then
                    echo -e "\nERROR: circos plots failed /!\ \n
                    please check logs and input data\n" 
                    exit 1
                fi
fi
fi 
echo "all analyses finished"
