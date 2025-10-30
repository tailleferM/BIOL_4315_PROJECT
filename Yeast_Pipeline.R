# Pipeline for yeast genome assembly and annotation,
# QC
#   sequali
# Assembly
#   Flye
# QC
#   BUSCO
#   QUAST
# Annotation
#   Funannotate

#pipeline test data : https://www.ebi.ac.uk/ena/browser/view/SRR10092046
 # long reads from ONT gridion

#QC 
  # sequali --outdir sequali_reports/ /path to the file
  
# Assembly
  # flye --nano-raw file.name --out-dir /flye_output --threads 8 --genome-size 14m

#assembly stats
Total length:	13858957
Fragments:	166
Fragments N50:	299060
Largest frg:	973058
Scaffolds:	0
Mean coverage:	30

# Annotation
# going off this doc https://funannotate.readthedocs.io/en/latest/tutorials.html

#Funannotate
#make sure you're in the right directory!!!
#step 1 masking (need to do before annotation)
docker run --rm -v "$PWD":/data -w /data --platform linux/amd64 nextgenusfs/funannotate funannotate mask -i trimmed_assembly.fasta -o tr_masked_assembly.fasta  

#step 2 prediction 
docker run --rm -v "$PWD":/data -w /data --platform linux/amd64 nextgenusfs/funannotate funannotate predict -i tr_masked_assembly.fasta -o /data/output -s 'Wickerhamomyces anomalus'

#step 3 annotation
docker run --rm -v "$PWD":/data -w /data --platform linux/amd64 nextgenusfs/funannotate funannotate annotate -i /data/output_attempt_1 -o /data/annotation_output

try giving it rna data -> reference genome 
align the reference rna seqs against my assembly -> use for funannotate





#reading my funannotate gff file output 
library(stringr)
gffFile <- readLines('/Users/mtaillefer00/Documents/BIOL_4315_PROJECT/data/funannotate_stuff/output_attempt_1/annotate_results/Wickerhamomyces_anomalus.gff3')
product_lines <- str_subset(gffFile, 'product=')
products <- str_match(product_lines, "product=([^;]+);")[,2]
hypo_count <- sum(str_detect(products, regex("^hypothetical", ignore_case = TRUE)))
filtered_products <- products[!str_detect(products, regex("^hypothetical", ignore_case = TRUE))]
summary_line <- paste(hypo_count, "hypothetical proteins")
final_list <- c(summary_line, filtered_products)

make list of genes involved in the pathways
 -> look at kegg or cog for pathways interesting for me 
 -> parse my list of products for these genes 
 
 
#busco
docker run --rm \
-v "$PWD":/data \
-w /data \
ezlabgva/busco:v6.0.0_cv1 \
busco -i assembly.fasta -o busco_output -l saccharomycetaceae_odb12 -m genome -c 8

library(blastinR)
make_blast_db(infile = '/Users/mtaillefer00/Documents/BIOL_4315_PROJECT/data/funannotate_stuff/output/predict_results/Wickerhamomyces_anomalus.proteins.fa', dbtype = 'prot')
build blast database of the genes im interested in
input metadata table assoc with genes in database



