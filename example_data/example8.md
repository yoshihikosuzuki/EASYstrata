# running the workflow for multi species/samples  


usefull to reduce the set of single copy orthologs making it more robust

example data are taken from cannabis paper with data downloaded here:

https://resources.michael.salk.edu/root/tools.html?tool=%2Fresources%2Fcannabis_genomes%2Findex.html

after getting the data for all genome and gene models we generated a table of single copy orthologs from 20 genome assemblies

* extract the X/Y peptide data from all peptide file

* run orthofinder

* export the N0.tsv.file

# formatting file

grep -v "," N0.tsv.file |awk 'NF==20' > single_copy_orthologs

next we want to obtain the dS pattern from one randomly chosen assembly pair with separated X/Y chromosomes 

# single copy orthologue file:  
 
the path to this file **MUST** be provided in the config/config file:

see last line:

single_copy_file="/home/quentin/10.collab/hemp/05_EASYstrata_from_SCO/single_copy_orthologs"

full example provided [here](example8.config)

# then the workflow is run as usual
