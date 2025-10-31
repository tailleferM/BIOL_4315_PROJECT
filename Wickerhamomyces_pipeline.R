# Wickerhamomyces run

# Sequali
# sequali --outdir sequali_reports/ Wickerhamomyces_anomalus.fastq
#   later positions are low quality  

# filtering and cleaning
library(QuasR)

outfiles <- 'above100bp.fastq'

preprocessReads('/Users/mtaillefer00/Documents/BIOL_4315_PROJECT/realData/Wickerhamomyces_anomalus.fastq', outfiles, minLength = 100)

# Assembly
# flye --nano-hq above100bp.fastq --out-dir flye/ --threads 8 --genome-size 14m

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
gffFile <- readLines('/Users/mtaillefer00/Documents/BIOL_4315_PROJECT/realData/funannotate/output/annotate_results/Wickerhamomyces_anomalus.gff3')
product_lines <- str_subset(gffFile, 'product=')
products <- str_match(product_lines, "product=([^;]+);")[,2]
hypo_count <- sum(str_detect(products, regex("^hypothetical", ignore_case = TRUE)))
filtered_products <- products[!str_detect(products, regex("^hypothetical", ignore_case = TRUE))]
summary_line <- paste(hypo_count, "hypothetical proteins")
final_list <- c(summary_line, filtered_products)


# KEGG