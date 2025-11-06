# Wickerhamomyces run

# Sequali
# sequali --outdir sequali_reports/ Wickerhamomyces_anomalus.fastq
#   later positions are low quality  

# filtering and cleaning
library(QuasR)

outfiles <- 'above100bp.fastq'

preprocessReads('/Users/mtaillefer00/Documents/BIOL_4315_PROJECT/Wickerhamomyces/Wickerhamomyces_anomalus.fastq', outfiles, minLength = 100)

# Assembly
# flye --nano-hq above100bp.fastq --out-dir flye/ --threads 8 --genome-size 14m

Total length:	13205329
Fragments:	2025
Fragments N50:	9214
Largest frg:	53218
Scaffolds:	0
Mean coverage:	13

#run quast and busco to check assembly quality
docker run --rm \
-v "$PWD":/data \
-w /data \
ezlabgva/busco:v6.0.0_cv1 \
busco -i flye_full/assembly.fasta -o busco_output -l saccharomycetaceae_odb12 -m genome -c 8
#49% COMPLETE BUSCOS

docker run --rm -v "$(pwd)":/data -w /data quay.io/biocontainers/quast:5.3.0--py313pl5321h5ca1c30_2 \
quast.py assembly.fasta \
--fungus \
-r GCF_001661255.1_Wican1_genomic.fna \
-g GCF_001661255.1_Wican1_genomic.gff \
-o quast_output2
# 75% complete genome 

#funannotate
#mask
# docker run --rm -v "$PWD":/data -w /data --platform linux/amd64 nextgenusfs/funannotate funannotate mask -i flye/assembly.fasta -o masked_assembly.fasta  
#Predict
# docker run --rm -v "$PWD":/data -w /data --platform linux/amd64 nextgenusfs/funannotate funannotate predict -i masked_assembly.fasta -o /data/output -s 'Wickerhamomyces anomalus'
#Annotate cd to be in funnanotate
# docker run --rm -v "$PWD":/data -w /data --platform linux/amd64 nextgenusfs/funannotate funannotate annotate -i /data/output -o /data/annotation_output

# parsing my gff file
library(stringr)
library(dplyr)
gffFile <- readLines('/Users/mtaillefer00/Documents/BIOL_4315_PROJECT/Wickerhamomyces/funannotate/output/annotate_results/Wickerhamomyces_anomalus.gff3')
product_lines <- str_subset(gffFile, 'product=')

annotation_df <- tibble(raw = product_lines) %>%
  mutate(
    tID = str_extract(raw, "(?<=ID=)[^;]+"),
    product = str_extract(raw, "(?<=product=)[^;]+"),
    pfams = str_extract(raw, "(?<=Dbxref=)[^;]+")
  ) %>%
  mutate(
    pfams = str_extract_all(pfams, "PF\\d+")
  ) %>%
  select(tID, product, pfams)

annotation_df

nonHypothetical <- annotation_df %>% filter(product != 'hypothetical protein')

pfams <- unique(unlist(annotation_df$pfams))
pfam_file <- 'listofpfams.txt'
writeLines(pfams, pfam_file)

https://www.ebi.ac.uk/interpro/api/entry/pfam/PF24563

# using the python program to get a list of InterPro accessions
./pfam.py pfams.txt iprs.txt

iprs <- readLines('/Users/mtaillefer00/Documents/BIOL_4315_PROJECT/wickerhamomyces/listofiprs.txt')

pfam_to_ipr <- setNames(iprs,pfams)

annotation_df$iprs <- lapply(annotation_df$pfams, function(key_vec) pfam_to_ipr[key_vec])

annotation_df$iprs <- lapply(annotation_df$iprs, unname)

unique_iprs <- unique(iprs)

ipr_file <- 'uniqueIprs.txt'
writeLines(unique_iprs, ipr_file)

# Funannotate attempt 2
#   first go I barely supplied funnanotate with anything: just the assembly
#   this time I will give it more info to work with. 

# using this tutorial for funannotate.
#https://training.galaxyproject.org/training-material/topics/genome-annotation/tutorials/funannotate/tutorial.html#hands-on-2

#eggnogdb mapper -> just used online tool
# -> r package that converts gene ontology numbers to something useful (biomart?)


#interproscan https://interproscan-docs.readthedocs.io/en/v5/HowToUseViaContainer.html

docker run --rm \
-v $PWD/interproscan-5.76-107.0/data:/opt/interproscan/data \
-v $PWD/input:/input \
-v $PWD/temp:/temp \
-v $PWD/output:/output \
--platform linux/amd64 \
interpro/interproscan:5.76-107.0 \
--input /input/Wickerhamomyces_anomalus.proteins.fa \
--output-dir /output \
--tempdir /temp \
--cpu 8

#IT WORKED!!!! took 5 hours

#https://link.springer.com/article/10.1007/s11274-023-03737-7#Sec1
# - paper about xenobiotics enzymes and pathways in fungi 
# - parse annotation for proteins mentioned in this paper
# - look into references