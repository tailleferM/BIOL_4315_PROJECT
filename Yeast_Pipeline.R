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

Total length:	13858957
Fragments:	166
Fragments N50:	299060
Largest frg:	973058
Scaffolds:	0
Mean coverage:	30

# Annotation
# going off this doc https://funannotate.readthedocs.io/en/latest/tutorials.html
#masking
docker run --rm -v "$PWD":/data -w /data --platform linux/amd64 nextgenusfs/funannotate funannotate mask -i trimmed_assembly.fasta -o tr_masked_assembly.fasta  

  # funannotate predict -i softMasked.fa -o fun \
    #--species "Pseudogenus specicus" --strain JMP12345 \
    #--busco_seed_species botrytis_cinerea --cpus 12

#docker
docker run --rm \
-v "$PWD":/data \
-w /data \
ezlabgva/busco:v6.0.0_cv1 \
busco -i assembly.fasta -o busco_output -l saccharomycetaceae_odb12 -m genome -c 8


