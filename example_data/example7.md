# Example 7: stickleback data from Peichel et al. 2020

this readme provides command lines to reproduce our analysis based on already published data
from [peichel et al. 2020](https://doi.org/10.1186/s13059-020-02097-x)
In this work the authors clearly identified three evolutionary strata

let's see if can find them in a first pass analysis.'

## 1 download published data

get the genome including the Y chromosome + gff :  

wget https://stickleback.genetics.uga.edu/downloadData/v5.0.1_assembly/stickleback_v5.0.1_assembly.fa.gz  

wget https://stickleback.genetics.uga.edu/downloadData/v5_assembly/stickleback_v5_maker_genes_chrY.gff3.gz  

wget https://stickleback.genetics.uga.edu/downloadData/v5_assembly/stickleback_v5_maker_genes_nath2020.gff3.gz  

```
sed 's/^>/>stickleback_/g' ../stickleback_v5.0.1_assembly.fa > stickleback.fa
zcat stickleback_v5_ensembl_genes.gff3.gz stickleback_v5_maker_genes_chrY.gff3.gz |sed 's/^chr/stickleback_chr/' > stickleback.gff3
```

#run busco on the genome and on the protein prediction:
busco_lineage="actinopterygii_odb10"
busco -c24 -o busco_genome -i stickleback_v4_assembly.fa -l "$busco_lineage" -m genome -f 

 ***** Results: *****

        C:87.9%[S:65.7%,D:22.2%],F:4.4%,M:7.7%,n:3640      
        3202    Complete BUSCOs (C)                        
        2393    Complete and single-copy BUSCOs (S)        
        809     Complete and duplicated BUSCOs (D)         
        159     Fragmented BUSCOs (F)                      
        279     Missing BUSCOs (M)                         
        3640    Total BUSCO groups searched                

Dependencies and versions:
        hmmsearch: 3.1
        python: sys.version_info(major=3, minor=7, micro=12, releaselevel='final', serial=0)
        busco: 5.7.1



```sh
cat stickleback_v5_ensembl_genes.gff3 stickleback_v5_maker_genes_chrY.gff3 > full.gff3

busco -c24 -o busco_genome -i stickleback_v4_assembly.fa -l "$busco_lineage" -m protein -f 

gffread -g stickleback_v5_assembly.fa -y prot.fa full.gff3

    C:88.8%[S:67.0%,D:21.8%],F:4.2%,M:7.0%,n:3640      
        3230    Complete BUSCOs (C)                        
        2438    Complete and single-copy BUSCOs (S)        
        792     Complete and duplicated BUSCOs (D)         
        152     Fragmented BUSCOs (F)                      
        258     Missing BUSCOs (M)                         
        3640    Total BUSCO groups searched                



```
 

we see a lot of duplication so this needs cleaning because there is multiple isoform

```
gffread full.gff3 -T -o full.gtf

#lets remove fully overlapping gene with a very simple approach
awk '$3=="transcript" {print $1"_"$4"_"$5"\t"$1"_"$4"\t"$1"_"$5"\t"$0}'  full.gtf > tmp
awk 'NR == FNR {count[$2]++; next} count[$2]>0 {print $2"\t"$7"\t"$8"\t"$13"\t"$8-$7}' tmp tmp \
    |awk '$5>max[$1]{max[$1]=$5; row[$1]=$4} END{for (i in row) print row[i]}' > longest.to.keep.tmp

grep -Ff longest.to.keep.tmp full.gtf > full.dedup1.gtf

awk '$3=="transcript" {print $1"_"$4"_"$5"\t"$1"_"$4"\t"$1"_"$5"\t"$0}'  full.dedup1.gtf > tmp2

awk 'NR == FNR {count[$3]++; next} count[$3]>0 {print $3"\t"$7"\t"$8"\t"$13"\t"$8-$7}' tmp2 tmp2 \
    |awk '$5>max[$1]{max[$1]=$5; row[$1]=$4} END{for (i in row) print row[i]}' > longest.to.keep2 

grep -Ff longest.to.keep2 full.dedup1.gtf > full.dedup2.gtf

#finally use AGAT to remove some remaining isoform:
agat_sp_keep_longest_isoform.pl -gff full.dedup2.gtf -o longestisoform.gff

gffread longestisoform.gff -T -o test.agat.gtf

gffread -g stickleback_v5_assembly.fa -y agat.prot test.agat.gtf 

busco -c24 -o busco_prot_agat -i agat.prot -l "$busco_lineage" -m protein -f  

sed 's/^chr/stickleback_chr/g' ../test.agat.gtf > tmp.gtf

awk '{gsub("_id \"", "_id \"" $1 "_", $0); print }' tmp.gtf > stickleback.gtf 
```


**Important note:** actually, the Y gff comes from another annotation tool (maker) and seem not
really compatible with the rest of the annotation.
Therefore we will use miniprot to transfer the stickleback X based gene to the Y chromosome.
This approach is similar to the one in Peichel et al. except that miniprot is surely faster than exonerate (and perhaps more accurate?)


## 2 - preparaing data for miniprot:
```
1. extracting prot without the Y:
gffread -g stickleback.fa  -y revalid.prot stickleback.gtf 
grep -v "stickleback_chrY" revalid.prot > stickleback.dedup.noY.peptide.fa
awk_linearize_fasta.sh stickleback.dedup.noY.peptide.fa 


#we extract the X chromosome (chrXIX in the assembly):  
grep -A1 "chrXIX" stickleback.dedup.noY.peptide.fa.lin.fasta > chrXIX.fa

#transfert with miniprot:
miniprot -t8 --gff --aln chrY.fa chrXIX.fa > chrY.miniprot.gff

#reshape:
grep -v "PAF\|##" chrY.miniprot.gff > chrY.gff
gffread chrY.gff -T -o chrY.gtf

grep -v "chrY" stickleback.gtf > no.Y.gtf
grep -v "chrXIX" stickleback.gtf > chrX.gtf

sed 's/^chrY/stickleback_chrY/g' chrY.gtf |awk '{gsub("_id \"", "_id \"" $1 "_", $0); print }' > chrY.renamed.gtf
cat no.Y.gtf chrY.renamed.gtf > stickleback.gtf
```


now we have all the data to perform the analysis. We must still prepare the config file

```
git clone https://github.com/QuentinRougemont/EASYstrata/

cd EASYstrata
```

edit the config file with vim, so it looks like this:

# config file
#_________________
# INPUT FOR ALL STEPS
# ----- Genome information -----
genome1="/home/quentin/04.pipelinevalidate/stickleback/V5/on_miniprot/X_Yonly/stickleback.fa"     #full path to either an assembly containing both haplotypes, or an assembly containing one of two haplotypes. Assemblies can be compressed (.gz) or not.
genome2=""     #full path to the second haplotypes (optional)
haplotype1="stickleback"  #name1 [name of haplotype1 - will be used to rename the genome and the contigs inside the genome]
haplotype2="stickleback_chrY"  #name2 [name of haplotype2 - will be used to rename the genome and the contigs inside the genome] (optional)

# ----- Ancestral reference (optional) -----
ancestral_genome="" #full path to the ancestral genome to be used whenever possible
ancestral_gtf=""   #full path to the associated gtf to be used whehever possible

# ----- Identification of the region of interest -----
scaffold="/home/quentin/04.pipelinevalidate/stickleback/V5/EASYstrata/scaffold" #full path to a 3 columns tab separated file containing information on scaffolds of interest (sex chromosomes or mating-type chromosomes) for plotting.
# from first to last column: genome name, scaffold name, a string ("normal" or "reversed") stating wether the scaffold orientation should be reversed

# ----- Already produced TE and gene annotation -----
# if step I is not performed or annotateTE is set to "NO", then TE can be provided in fasta format (optional):
#TEgenome1="" #full path to the bed of TE for genome1
#TEgenome2="" #full path to the bed of TE for genome2
#TEancestral="" #full path to the bed of TE for the ancestral genome

# if step I is not performed or annotate is set to "NO", then gtf must be provided: 
gtf1="/home/quentin/04.pipelinevalidate/stickleback/V5/on_miniprot/X_Yonly/stickleback.gtf" #full path to gtf for genome1
gtf2="" #full path to gtf for genome2









## now run EASYstrata with option 3 
note : this can run on a laptop in less than an hour, depending on the MCMC size in MCP

```
./master.sh -o3  2>&1 |tee log_stickleback_o3
```


#here are some results and plots 

- the first tool run was Genespace. Here is the riparian plot : 

![Fig1.stickleback.png](https://github.com/QuentinRougemont/EASYstrata/blob/main/.pictures/Fig1.stickleback.png)



- number of single copy orthologs from the N0.tsv file from Orthofinder :


- here s a circos plot of all single copy orthologs: 

![Fig2.stickleback.png](https://github.com/QuentinRougemont/EASYstrata/blob/main/.pictures/Fig2.stickleback.png)


- mean dS: 



- plots of dS values along the whole X chromosome : 


![Fig3.stickleback.png](https://github.com/QuentinRougemont/EASYstrata/blob/main/.pictures/Fig3.stickleback.png)



- some plots extracted from the MCP analysis are provided here:  

![Fig4.stickleback.png](https://github.com/QuentinRougemont/EASYstrata/blob/main/.pictures/Fig4.stickleback.png)



![Fig5.stickleback.png](https://github.com/QuentinRougemont/EASYstrata/blob/main/.pictures/Fig5.stickleback.png)




