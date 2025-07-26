#!/bin/bash
#Purpose:
#master script to prepare bed files, laucnch GeneSpace, run paml and launch downstream Rscript 
#will check existence of all dependencies
#Date: 2023
#Author: QR

# -- some colors for warnings in the terminal  --:
source config/colors
############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo -e "master script to: \n 1 - create bed files, \n 2 - launch GeneSpace, \n 3 - run paml and launch downstream Rscripts (Rideogram, plot paml, etc)"
   echo " "
   echo "Usage: $0 [-s1|-s2|-f|-a|-g|-h|]"
   echo "options:"
   echo " -h|--help: Print this Help."
   echo " -s1|--haplo1: the name of the first  focal haplotype "
   echo " -s2|--haplo2: the name of the second focal haplotype "
   echo " -a|--ancestral_genome: the name of the ancestral haplo to infer orthology and plot gene order"
   echo " -g|--ancestral_gtf: the name of the ancestral gtf associated with the ancestral genome"
   echo " -f|--folderpath: the path to the global folder containing haplo1 and haplo 2"
   echo " -c|--chromosome: a tab separated txt file listing the name of the reference species (e.g sp1), the corresponding set of chromosomes (e.g.: chrX , supergene, etc) and the orientation of the chromosome (N: Normal, R: Reverse) if their is more than one"
   echo " -o|--options : the type of analysis to be performed: either 'synteny_and_Ds' (GeneSpace+Minimap2+Ds+changepoint), 'synteny_only' (GeneSpace+Minimap2), 'Ds_only' (paml and changepoint)"
   echo " "
   echo "dependancies: orthofinder, mcscanx, GeneSpace, paml (yn00), Rideogram, translatorX minimap2"
}

############################################################
# Process the input options.                               #
############################################################
while [ $# -gt 0 ] ; do
  case $1 in
    -s1 | --haplo1) haplo1="$2" ; echo -e "haplotype 1 Name is ***${haplo1}*** \n" >&2;;
    -s2 | --haplo2) haplo2="$2" ; echo -e "haplotype 2 Name is ***${haplo2}*** \n" >&2;;
    -a  | --ancestral_genome) ancestral_genome="$2" ; 
        echo -e "ancestral haplo  Name is ***${ancestral_genome}*** \n" >&2;;
    -g  | --ancestral_gtf) ancestral_gtf="$2" ; 
        echo -e "ancestral gtf  Name is ***${ancestral_gtf}*** \n" >&2;;
    -f  | --folderpath  ) folderpath="$2"   ; 
        echo -e "global folder is  ${folderpath} \n" >&2;;
    -c  | --chromosome )  chromosome="$2"   ; 
        echo -e "target chromosome are ${chromosome} \n" >&2 ;; 
    -o  | --options ) options="$2" ; 
        echo -e "options for computation are ***${options}*** \n" >&2 ;;
    -h  | --help ) Help ; exit 2 ;;
   esac
   shift
done 

if [ -z "${haplo1}" ] || [ -z "${haplo2}" ] || [ -z "${chromosome}" ] || [ -z "${options}" ]  ; then
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

scaffold=$chromosome
if [ -z ${scaffold+x} ]; then
        echo "ERROR! scaffold file is unset";
        echo "please provide a file of target scaffold (e.g. X or ancestral state)"
        echo "see example in example_data/scaffold.txt"
        exit 1
else
        echo "scaffold file is set to '$scaffold'";
fi

#replace space (multiple or not) by a single tab:
sed -i 's/ \+/\t/g' "$scaffold"
if [ ! -d 02_results ] ; then mkdir 02_results ; fi #ignore if already existent

#make ancestral species optional
if [ -n "${ancestral_genome}" ] ; then
    echo "ancestral_species is $ancestral_genome "
    echo "will attempt to extract the CDS and PROT from it "
    if [ ! -d 02_results ] ; then mkdir ancestral_sp/03_genome ; fi #note: this folder exist already 
    # ----- check compression of fasta  ------ ##
    #check compression
    if file --mime-type "$ancestral_genome" | grep -q gzip$; then
       echo "$ancestral_genome is gzipped"
       gunzip "$ancestral_genome"
       ancestral_genome="${ancestral_genome%.gz}"
    else
       #trim any eventual .gz extension from file
       ancestral_genome="${ancestral_genome%.gz}"
       echo "$ancestral_genome is not gzipped"

    fi
   
    if file --mime-type "$ancestral_gtf" | grep -q gzip$; then
       echo "$ancestral_gtf is gzipped"
       gunzip "$ancestral_gtf"
       ancestral_gtf="${ancestral_gtf%.gz}"
    else
       echo "$ancestral_gtf is not gzipped"
       #trim any eventual .gz extension from file
       ancestral_gtf="${ancestral_gtf%.gz}"

    fi

    #check if gtf or gff:
    if file --mime-type "$ancestral_gtf" | grep -q gff; then
        echo file is gff ;
        echo converting into gtf
        gffread "$ancestral_gtf" -T -o tmp
        ancestral_gtf=${ancestral_gtf%.gff*}.gtf
        mv tmp "$ancestral_gtf" 
        sed -i -E '/ancestral_gtf/ s/.gff./.gtf/' config/config
    else
        echo " "
    fi
    if [ ! -d ancestral_sp ] ; then mkdir ancestral_sp ; fi
    cd ancestral_sp || exit 
    if [ -f ancestral_sp.fa ] ; then
        rm ancestral_sp.fa
    fi
    
    ln -s "${ancestral_genome}" ancestral_sp.fa ; 
    samtools faidx ancestral_sp.fa ; 
    cd ../
    gffread -g "${ancestral_genome}" -x ancestral_sp/ancestral_sp.spliced_cds.fa  \
    "${ancestral_gtf}" 
    transeq -sequence ancestral_sp/ancestral_sp.spliced_cds.fa \
        -outseq ancestral_sp/ancestral_sp_prot.fa
    awk '$3=="transcript" {print $1"\t"$4"\t"$5"\t"$10}' "$ancestral_gtf" |\
        sed 's/"//g' |sed 's/;//g'  > ancestral_sp/ancestral_sp.bed
    sed -i 's/_1//g'  ancestral_sp/ancestral_sp_prot.fa
    gene_ids_hap=$(cut -f 4 ancestral_sp/ancestral_sp.bed |awk -F "_" '{print NF}' |sort |uniq )
    if [ "$gene_ids_hap" = 3 ] ; then
        echo "gene_id hap1 is ok"
    else
       echo "error! gene structure should be of the type:"
       echo "[IndividualID]_[chromosomeID]_[geneID]" 
       echo -e "it should contain exactly two underscore to separate your individual/strain/species
             from the chromosome, and the chromosome from the geneID\n"
       echo "please reformat the ID in your gff/gtf"   
       exit 1
    fi
    bedanc="ancestral_sp/ancestral_sp.bed"
fi

#----------- step 1 prepare bed file for each haplo ----------------------------#
#
#test options :
if [[ $options = "Ds_only" ]] ;  
then
    mkdir -p 02_results/paml 
elif [[ $options = "plots" ]] ; 
then
    echo "only drawing plots" 
elif [[ $options = "synteny_and_Ds" ]] ; 
then
    if [ ! -d genespace ] ; then  
       mkdir -p genespace/bed genespace/peptide 
       mkdir -p  02_results/paml #02_results/plots 
    fi
elif [[ $options = "synteny_only" ]] ; 
then
    if [ ! -d genespace ] ; then  
       mkdir -p genespace/bed genespace/peptide 
       mkdir 02_results/paml
    fi
elif [[ $options == "changepoint" ]] ;
then
        echo "only changepoint will be performed"
fi

# create bed
gffread --bed -E haplo1/08_best_run/"$haplo1".final.gtf -o haplo1/08_best_run/"$haplo1".bed 
gffread --bed -E haplo2/08_best_run/"$haplo2".final.gtf -o haplo2/08_best_run/"$haplo2".bed 
cut -f 1-4  haplo1/08_best_run/"$haplo1".bed  > haplo1/08_best_run/"$haplo1".v2.bed
cut -f 1-4  haplo2/08_best_run/"$haplo2".bed  > haplo2/08_best_run/"$haplo2".v2.bed

bedhaplo1="haplo1/08_best_run/$haplo1.v2.bed"
bedhaplo2="haplo2/08_best_run/$haplo2.v2.bed"

#checking the architeture of gene id:
gene_ids_hap1=$(cut -f 4 haplo1/08_best_run/"$haplo1".v2.bed | awk -F "_" '{print NF}' |sort |uniq )
if [ "$gene_ids_hap1" = 3 ] ; then
    echo "gene_id hap1 is ok"
else
   echo "error! gene structure should be of the type:"
   echo "[IndividualID]_[chromosomeID]_[geneID]" 
   echo -e "it should contain exactly two underscore to separate your individual/strain/species
         from the chromosome, and the chromosome from the geneID\n"
   echo "please reformat the ID in your gff/gtf"   
   exit 1
fi

gene_ids_hap2=$(cut -f 4 haplo2/08_best_run/"$haplo2".v2.bed | awk -F "_" '{print NF}' |sort |uniq )
if [ "$gene_ids_hap2" = 3 ] ; then
    echo "gene_id hap1 is ok"
else
   echo "error! gene structure should be of the type:"
   echo "[IndividualID]_[chromosomeID]_[geneID]" 
   echo -e "it should contain exactly two underscore to separate your individual/strain/species
         from the chromosome, and the chromosome from the geneID\n"
   echo "please reformat the ID in your gff/gtf"   
   exit 1
fi

#------------------------------ step 2 prepare genespace data------------------#
#
if [[ $options = "synteny_and_Ds" ]] || [[ $options = "synteny_only" ]] ; 
then
    path1=$(readlink -f  haplo1/08_best_run/"$haplo1".v2.bed)
    path2=$(readlink -f haplo2/08_best_run/"$haplo2".v2.bed) 
    if [ -f genespace/bed/"$haplo1".bed ] ; then 
        rm genespace/bed/"$haplo1".bed ;
    fi
    if [ -f genespace/bed/"$haplo2".bed ] ; then 
        rm genespace/bed/"$haplo2".bed ;
    fi

    ln -s "$path1" genespace/bed/"$haplo1".bed
    ln -s "$path2" genespace/bed/"$haplo2".bed
    
    # simplify the protein file to match the bed (i.e. remove the _1 inserted by transeq and the CDS length info):
    #sed 's/_1.*$//g' haplo1/08_best_run/"$haplo1"_prot.final.clean.fa \
    #    > genespace/peptide/"$haplo1".fa
    #sed 's/_1.*$//g' haplo2/08_best_run/"$haplo2"_prot.final.clean.fa \
    #    > genespace/peptide/"$haplo2".fa
    if [ ! -s genespace/peptide/"$haplo1".fa ] ; then
         cp haplo1/08_best_run/"$haplo1"_prot.final.clean.fa \
         genespace/peptide/"$haplo1".fa
         cp haplo2/08_best_run/"$haplo2"_prot.final.clean.fa \
         genespace/peptide/"$haplo2".fa
    fi 
    #verify that IDs in bed and fasta file are matching - else exit  
    grep ">" genespace/peptide/"$haplo1".fa |sed 's/>//g' > tmp1
    grep ">" genespace/peptide/"$haplo2".fa |sed 's/>//g' > tmp2
    
    check1=$(grep -Ff tmp1 genespace/bed/"$haplo1".bed |wc -l )
    check2=$(grep -Ff tmp2 genespace/bed/"$haplo2".bed |wc -l )
    
    echo -e "number of protein in haplo1 present in bed is: $check2"
    echo -e "number of protein in haplo2 present in bed is:  $check1"
    
    bedsize1=$(wc -l genespace/bed/"$haplo1".bed |awk '{print $1}' )
    bedsize2=$(wc -l genespace/bed/"$haplo2".bed |awk '{print $1}' )
    
    echo -e "bedisze2  size is $bedsize2"
    echo -e "bedisze1  size is $bedsize1"
    
    #check that all is matching:
    if [ "$bedsize1" = "$check1" ]
    then
        echo "input1 passing quality check" 
        rm tmp1
    else
        echo "input1 NOT passing quality check"
        echo "check your data"
        if [ "$check1" -lt "$bedsize1" ] ;
        then
            #there is less protein than bed size, we subset the bed
            echo -e "\n\nattempting to subset the bed for haplo1\n"
            grep -Ff tmp1 genespace/bed/"$haplo1".bed > bed.hap1.tmp
            mv bed.hap1.tmp genespace/bed/"$haplo1".bed
        #exit 2
        fi
    fi
    
    if [ "$bedsize2" = "$check2" ]
    then
        echo "input2 passing quality check" 
        rm tmp2
    else
        echo "input2 NOT passing quality check"
        echo "check your data"
        if [ "$check2" -lt "$bedsize2" ] ;
        then
            #there is less protein than bed size, we subset the bed
            echo -e "\n\nattempting to subset the bed for haplo2\n"
            grep -Ff tmp2 genespace/bed/"$haplo2".bed > bed.hap2.tmp
            mv bed.hap2.tmp genespace/bed/"$haplo2".bed
        #exit 2
        fi
    fi
    # -- handling ancestral haplo ------
    # -- this part assumes that a bed and peptide file are existant for the 
    # -- ancestral haplo
    if [ -n "${ancestral_genome}" ] ; then
        cd genespace/bed/ || exit 1
        if [ ! -f ancestral_sp.bed ] ;then
            ln -s ../../ancestral_sp/ancestral_sp.bed . 
        fi 
        cd ../peptide || exit 1
        if [ ! -f ancestral_sp.fa ] ; then
            ln -s ../../ancestral_sp/ancestral_sp_prot.fa ancestral_sp.fa
        fi
        cd ../../
    fi
fi
#------------------- step 3 run GeneSpace -------------------------------------#
if [[ $options = "synteny_and_Ds" ]]  || [[ $options = "synteny_only" ]] ; then
    if [ -s genespace/results/gsParams.rda ] ; then
        echo -e "\ngenespace results already exist - skipping run\n"
    else
        cd genespace  || exit 1
        if [ -s orthofinder ] ; then
            #a previous unfinished run exist - remove its content:
            rm -rf genespace/dotplots/ genespace/orthofinder/ \
                genespace/pangenes/ genespace/results/ \
                genespace/riparian/ genespace/syntenicHits/ genespace/tmp
        fi
        MCScanpath=$(command -v MCScanX |xargs dirname )
        #just in case this is not already done:
        Rscript -e  'devtools::install_github("jtlovell/GENESPACE")'
        sed -i "s#mcpath#$MCScanpath#" ../00_scripts/Rscripts/01.run_geneSpace.R
        
        Rscript ../00_scripts/Rscripts/01.run_geneSpace.R || exit 1
        #plot genespace subspace of target chromosomes: 
        echo scaffold is "$scaffold"
        if [ ! -S scaffold.txt ] ; then
            ln -s "$scaffold" scaffold.txt
        fi
        echo -e "---------- making subplots using scaffold data ----------------"
        if [ -n "${ancestral_genome}" ] ; then
            Rscript ../00_scripts/Rscripts/02.plot_geneSpace.R ancestral_sp
        else
            Rscript ../00_scripts/Rscripts/02.plot_geneSpace.R "$haplo1"
        fi
        cd ../
        
        if ls genespace/*.pdf 1> /dev/null 2>&1 ; then
            cp genespace/*pdf 02_results/
        else
            echo "plot in genespace/*pdf does not exist"
            echo "check your scaffold file and whether genespace run properly"
            #exit 2
        fi
    fi 

echo -e "\n-- extract single copy orthologs from orthofinder for later use --\n"
        pathN0="genespace/orthofinder/Results_*/Phylogenetic_Hierarchical_Orthogroups/N0.tsv"
        sed -i -e "s/\r//g" $pathN0
        #check size 
        minsize=30  
        #minisize =minimum number of single copy orthologs for further analyses. 
        #This is really a low bound....

        haplo=$(head -n1 $pathN0 |awk '{print $7}')

        if [  -z "${ancestral_genome}" ] ; then 
            p1=$(awk -v hap="$haplo1" '{for(i=1;i<=NF;++i)if($i ~ hap )print $i}' <(grep -v "," $pathN0 |grep -Ff <(awk '{print $2}' $scaffold) | awk 'NF==5' ) )
            p2=$(awk -v hap="$haplo2" '{for(i=1;i<=NF;++i)if($i ~ hap )print $i}' <(grep -v "," $pathN0 |grep -Ff <(awk '{print $2}' $scaffold) | awk 'NF==5' ) )
            size1=$(paste <(echo "$p1") |wc -l )
            size2=$(paste <(echo "$p2") |wc -l )
            echo "there is $size1 single copy orthologs in $haplo1 data"
            echo "there is $size2 single copy orthologs in $haplo2 data"

            if [ "$size1" != "$size2" ] ; then 
                echo "error! number of single copy orthologs in $haplo1 and $haplo2 are not identical" ; 
                echo "please check your single copy orthologs file"
                exit 
            fi
            paste <(echo "$p1") <(echo "$p2") |\
                    awk '{print "OG\t"$0}' > 02_results/paml/single.copy.orthologs 
            
            filesize=$(wc -l 02_results/paml/single.copy.orthologs |awk '{print $1}' )
            if [ "$filesize" -lt "$minsize" ] ;
            then
                echo "error! file 02_results/paml/single.copy.orthologs is empty"
                echo "please check your single copy orthologs in N0.tsv"
                echo "please check your input data as well"
                exit 1
            fi 
            
        elif [ -n "${ancestral_genome}" ]
        then
            ancestral=$(head -n1 ancestral_sp/ancestral_sp.fa.fai \
                |cut -f1 \
                |awk '{gsub("_","\t",$0) ; print $1}')

            p1=$(awk -v hap="$haplo1" '{for(i=1;i<=NF;++i)if($i ~ hap )print $i}' <(grep -v "," $pathN0 |grep -Ff <(awk '{print $2}' $scaffold) | awk 'NF==6' ) )
            p2=$(awk -v hap="$haplo2" '{for(i=1;i<=NF;++i)if($i ~ hap )print $i}' <(grep -v "," $pathN0 |grep -Ff <(awk '{print $2}' $scaffold) | awk 'NF==6' ) )
            size1=$(paste <(echo "$p1") |wc -l )
            size2=$(paste <(echo "$p2") |wc -l )
            echo "there is $size1 single copy orthologs in $haplo1 data"
            echo "there is $size2 single copy orthologs in $haplo2 data"

            if [ "$size1" != "$size2" ] ; then 
                echo "error! number of single copy orthologs in $haplo1 and $haplo2 are not identical" ; 
                echo "please check your single copy orthologs file"
                exit 
            fi

            p_anc=$(awk -v hap="$ancestral" '{for(i=1;i<=NF;++i)if($i ~ hap )print $i}' <(grep -v "," $pathN0 |grep -Ff <(awk '{print $2}' $scaffold) | awk 'NF==6' ) ) 
            size_anc=$(paste <(echo "$p_anc") |wc -l )
            if [ "$size1" != "$size_anc" ] ; then 
                echo "error! number of single copy orthologs in ancestral species is different from haplo1"
                echo "there is $size1 single copy orthologs in $haplo1 data"
                echo "there is $size_anc single copy orthologs in ancestral data"
                exit
            fi
            paste <(echo "$p_anc" ) <(echo "$p1") <(echo "$p2") |\
                awk '{print "OG\t"$0}' > 02_results/paml/single.copy.orthologs 
            filesize=$(wc -l 02_results/paml/single.copy.orthologs |awk '{print $1}' )
            if [ "$filesize" -lt "$minsize" ] ;
            then
                echo "error! file 02_results/paml/single.copy.orthologs is empty"
                echo "please check your single copy orthologs in N0.tsv"
                echo "please check your input data as well"
                exit 1
            fi 

        fi  
        #create orthologues file:
        awk '{print "ortho1\tortho2\t"$0}' 02_results/paml/single.copy.orthologs > 02_results/orthologues
#------------------- step 4 run minimap2 --------------------------------------#
    echo -e "\n-----------------\n   perform whole genome synteny\n------------\n" 
    echo -e "\n-----------------\n   running minimap\n-------------------------\n" 
    
    if [ ! -d "02_results/minimap_alns" ] ; then mkdir 02_results/minimap_alns/ ; fi
    if [ ! -s 02_results/minimap_alns/aln."$haplo1"_"$haplo2".paf ] ;
    then
        echo "running $haplo1 vs $haplo2"
        #create paf file:
        minimap2 -cx asm5 \
            haplo1/03_genome/"$haplo1".fa \
            haplo2/03_genome/"$haplo2".fa \
            > 02_results/minimap_alns/aln."$haplo1"_"$haplo2".paf || \
            { echo -e "${RED} ERROR! minimap2 failed - check your data\n${NC} " ; exit 1 ; }

    fi 
    if [ -n "${ancestral_genome}" ] ; then
        echo -e "\n------- an ancestral genome was provided ------ "
        if [ ! -s 02_results/minimap_alns/aln.ancestral_sp_"$haplo2".paf ];
        then
            echo -e "running minimap for genome broad synteny plots  -------\n" 
            minimap2 -cx asm5 \
                ancestral_sp/ancestral_sp.fa \
                haplo2/03_genome/"$haplo2".fa \
                > 02_results/minimap_alns/aln.ancestral_sp_"$haplo2".paf || \
            { echo -e "${RED} ERROR! minimap2 faield - check your data\n${NC} " ; exit 1 ; }
        fi
        if [ ! -s 02_results/minimap_alns/aln.ancestral_sp_"$haplo1".paf ];
        then
            echo -e "running minimap for genome broad synteny plots  -------\n" 
            minimap2 -cx asm5 \
                ancestral_sp/ancestral_sp.fa \
                haplo1/03_genome/"$haplo1".fa \
                > 02_results/minimap_alns/aln.ancestral_sp_"$haplo1".paf  || \
            { echo -e "${RED} ERROR! minimap2 faield - check your data\n${NC} " ; exit 1 ; }
        fi 
        #preparing scaffold to highlight in minimap2 dotplot:
        #Note: this code assume the following structure of gene id: "haploXX_chrXX_geneXXX"
        #no other separator should be used 
        awk '{gsub("_","\t",$0) ; print $2"_"$3"\t"$5"_"$6}' 02_results/paml/single.copy.orthologs|\
            sort |\
            uniq -c|\
            awk '$1>10 '  > 02_results/scaff.anc.haplo1.txt
        if [ ! -s 02_results/scaff.anc.haplo1.txt ] ; then 
            echo "WARNING FILE 02_results/scaff.anc.haplo1.txt is empty" 
            echo "this may cause problem later "
        fi  
        awk '{gsub("_","\t",$0) ; print $2"_"$3"\t"$8"_"$9}' 02_results/paml/single.copy.orthologs|\
            sort |\
            uniq -c\
            |awk '$1>10 ' > 02_results/scaff.anc.haplo2.txt
        if [ ! -s 02_results/scaff.anc.haplo2.txt ] ; then 
            echo "WARNING FILE 02_results/scaff.anc.haplo2.txt is empty" 
            echo "this may cause problem later "
        fi  

        awk '{gsub("_","\t",$0) ; print $5"_"$6"\t"$8"_"$9}' 02_results/paml/single.copy.orthologs|\
            sort |\
            uniq -c|\
            awk '$1>10 ' \
            |sed -e 's/^    //g' -e 's/ /\t/g'  > 02_results/scaff.haplo1.haplo2.txt
        if [ ! -s 02_results/scaff.haplo1.haplo2.txt ] ; then 
            echo "WARNING FILE 02_results/scaff.haplo1.haplo2.txt is empty" 
            echo "this may cause problem later "
        fi  

        if [ ! -d Rlogs ] ; then mkdir Rlogs ; fi 
        Rscript 00_scripts/Rscripts/dotplot_paf.R  02_results/minimap_alns/aln."$haplo1"_"$haplo2".paf  2> Rlogs/Rlogs_minimap1 
        Rscript 00_scripts/Rscripts/dotplot_paf.R  02_results/minimap_alns/aln.ancestral_sp_"$haplo1".paf  2> Rlogs/Rlogs_minimap2 
        Rscript 00_scripts/Rscripts/dotplot_paf.R  02_results/minimap_alns/aln.ancestral_sp_"$haplo2".paf  2> Rlogs/Rlogs_minimap3 

        Rscript 00_scripts/Rscripts/synteny_plot.R 02_results/minimap_alns/aln.ancestral_sp_"$haplo1".paf \
             02_results/scaff.anc.haplo1.txt  2> Rlogs/Rlogs_minimap4 
        Rscript 00_scripts/Rscripts/synteny_plot.R 02_results/minimap_alns/aln.ancestral_sp_"$haplo2".paf \
            02_results/scaff.anc.haplo2.txt  2> Rlogs/Rlogs_minimap5 
        Rscript 00_scripts/Rscripts/synteny_plot.R 02_results/minimap_alns/aln."$haplo1"_"$haplo2".paf \
            02_results/scaff.haplo1.haplo2.txt  2> Rlogs/Rlogs_minimap6 
    
    else 
        awk '{gsub("_","\t",$0) ; print $2"_"$3"\t"$5"_"$6}' 02_results/paml/single.copy.orthologs|\
            sort |\
            uniq -c|\
            awk '$1>10 ' \
           |sed -e 's/^    //g' -e 's/ /\t/g'  > 02_results/scaff.haplo1.haplo2.txt
        if [ ! -d Rlogs ] ; then mkdir Rlogs ; fi 
        #then run pafr to generate a whole genome dotplot and eventually dotplot for some target scaffold:
        Rscript 00_scripts/Rscripts/dotplot_paf.R  02_results/minimap_alns/aln."$haplo1"_"$haplo2".paf 2> Rlogs/Rlogs_minimap1 
        Rscript 00_scripts/Rscripts/synteny_plot.R 02_results/minimap_alns/aln."$haplo1"_"$haplo2".paf \
            02_results/scaff.haplo1.haplo2.txt  2> Rlogs/Rlogs_minimap2 
    fi
fi

##------------------- step 5 some data processing ------------------------------#
##optional - to be optimize: 
if [[ $options = "synteny_and_Ds" ]]  || [[ $options = "synteny_only" ]] ; then
    if [ -n "${ancestral_genome}" ]; then
        #ancestral genome exist
        #./00_scripts/extract_singlecopy.sh -h1 "$haplo1" -h2 "$haplo2" -s "$scaffold" -a ancestral_sp
        #preparing scaffold to highlight in dotplot:
        awk '{gsub("_","\t",$0) ; print $2"_"$3"_"$4"\t"$6"_"$7}' 02_results/paml/single.copy.orthologs|\
            sort |\
            uniq -c|\
            awk '$1>4 '  > 02_results/scaff.anc.haplo1.txt
        awk '{gsub("_","\t",$0) ; print $2"_"$3"_"$4"\t"$9"_"$10}' 02_results/paml/single.copy.orthologs|\
            sort |\
            uniq -c\
            |awk '$1>4 ' > 02_results/scaff.anc.haplo2.txt
        awk '{gsub("_","\t",$0) ; print $6"_"$7"\t"$9"_"$10}' 02_results/paml/single.copy.orthologs|\
            sort |\
            uniq -c|\
            awk '$1>4 ' \
            |sed -e 's/^    //g' -e 's/ /\t/g' > 02_results/scaff.haplo1.haplo2.txt

    else
        #ancestral genome not provided  
        awk '{gsub("_","\t",$0) ; print $2"_"$3"\t"$5"_"$6}' 02_results/paml/single.copy.orthologs|\
            sort |\
            uniq -c|\
            awk '$1>4 ' \
           |sed -e 's/^    //g' -e 's/ /\t/g'  > 02_results/scaff.haplo1.haplo2.txt

    fi
fi

if [ -n "${ancestral_genome}" ]; then
    cut  -f3 02_results/paml/single.copy.orthologs > 02_results/paml/sco."$haplo1".txt 
    cut  -f4 02_results/paml/single.copy.orthologs > 02_results/paml/sco."$haplo2".txt 
else
    cut  -f2 02_results/paml/single.copy.orthologs > 02_results/paml/sco."$haplo1".txt 
    cut  -f3 02_results/paml/single.copy.orthologs > 02_results/paml/sco."$haplo2".txt 
fi

if [ "$options" = "Ds_only" ] ; then
    source config/config
    echo "computing only dS and plotting results"
    # first option: single copy orthologs provided by the user:
    if [ -r "$single_copy_file" ] ; then
        echo "single copy ortholog file is : $single_copy_file"
        single_copy_size=$(wc -l "$single_copy_file" |awk '{print $1}' ) 
        echo -e "there is $single_copy_size orthologs in file $single_copy_file" 
        
        #matching column containing the haplotype 1 :
        p1=$(awk -v hap1="$haplo1" '{for(i=1;i<=NF;++i)if($i ~ hap1 )print $i}' "$single_copy_file" )
        size1=$(paste <(echo "$p1") |wc -l )
        
        #matching column containing the haplotype 2 :
        p2=$(awk -v hap2="$haplo2" '{for(i=1;i<=NF;++i)if($i ~ hap2 )print $i}' "$single_copy_file" )
        size2=$(paste <(echo "$p2") |wc -l )
        
        #checking number of orthologs in haplo1 and haplo2, these should be identical
        if [ "$size1" != "$size2" ] ; then 
            echo "error! number of single copy orthologs in $haplo1 and $haplo2 are not identical" ; 
            echo "please check your single copy orthologs file"
        else 
            if [ -n "${ancestral_genome}" ]; then
                p_anc=$(awk -v anc="$ancestral" '{for(i=1;i<=NF;++i)if($i ~ anc )print $i}' "$single_copy_file" )
                size_anc=$(paste <(echo "$p_anc") |wc -l )
                if [ "$size1" != "$size_anc" ] ; then 
                    echo "error! number of single copy orthologs in ancestral species is different from haplo1"
                    exit
                else
                    paste <(echo "$p_anc" ) <(echo "$p1") <(echo "$p2") |\
                        awk '{print "OG\t"$0}' > 02_results/paml/single.copy.orthologs 
                    if [ ! -s "02_results/paml/single.copy.orthologs" ];then
                        echo "error! file 02_results/paml/single.copy.orthologs is empty"
                        echo "please check your single copy orthologs in N0.tsv"
                        echo "please check your input data as well"
                        exit 1
            fi
                    cut  -f3 02_results/paml/single.copy.orthologs > 02_results/paml/sco."$haplo1".txt 
                    cut  -f4 02_results/paml/single.copy.orthologs > 02_results/paml/sco."$haplo2".txt
                fi
            else
                echo "no ancestral species"
                paste <(echo "$p1") <(echo "$p2") |\
                    awk '{print "OG\t"$0}' > 02_results/paml/single.copy.orthologs 
                if [ ! -s "02_results/paml/single.copy.orthologs" ];then
                    echo "error! file 02_results/paml/single.copy.orthologs is empty"
                    echo "please check your single copy orthologs in N0.tsv"
                    echo "please check your input data as well"
                    exit 1
            fi
                cut  -f2 02_results/paml/single.copy.orthologs > 02_results/paml/sco."$haplo1".txt 
                cut  -f3 02_results/paml/single.copy.orthologs > 02_results/paml/sco."$haplo2".txt
            fi         
            #fi
        fi  
    #second option : 
    elif [ -f  "$single_copy_file" ] ; then
        echo "the file $single_copy_file" cannot be read
        echo "check your file"
        exit 1
    #third option: run from a previous genespace run: 
    elif [[ -f "02_results/paml/single.copy.orthologs"  ]] ; then
        #assumption: file "02_results/scaff.haplo1.haplo2.txt and file 02_results/scaff.anc.haploX.txt exist"
        if [ -n "${ancestral_genome}" ]; then
            cut  -f3 02_results/paml/single.copy.orthologs > 02_results/paml/sco."$haplo1".txt 
            cut  -f4 02_results/paml/single.copy.orthologs > 02_results/paml/sco."$haplo2".txt
        else
            cut  -f2 02_results/paml/single.copy.orthologs > 02_results/paml/sco."$haplo1".txt
            cut  -f3 02_results/paml/single.copy.orthologs > 02_results/paml/sco."$haplo2".txt
        fi
    #fourth option: no file provided nor genespace run:
    else 
    echo "error no single copy orthologs file exist"
    echo "check your data" 
    exit 1
    fi
fi
