#Sungouiella akabanensis run


# Sequali
# conda activate sequali-env
#sequali --outdir sequali_reports/ sungouiellaAkabanensis.fastq

# filtering
library(QuasR)

outfiles <- 'Sungouiella/above140bp.fastq'

preprocessReads('/Users/mtaillefer00/Documents/BIOL_4315_PROJECT/Sungouiella/sungouiellaAkabanensis.fastq', outfiles, minLength = 140)

# Assembly
# conda activate fly-env
# flye --nano-hq above140bp.fastq --out-dir flye/ --threads 8 --genome-size 14m
#   Total length:	12064457
#   Fragments:	90
#   Fragments N50:	363453
#   Largest frg:	815303
#   Scaffolds:	0
#   Mean coverage:	18

#run quast and busco to check assembly quality
docker run --rm \
-v "$PWD":/data \
-w /data \
ezlabgva/busco:v6.0.0_cv1 \
busco -i flye/assembly.fasta -o busco_output -l saccharomycetaceae_odb12 -m genome -c 8
# 43% complete  

docker run --rm -v "$(pwd)":/data -w /data quay.io/biocontainers/quast:5.3.0--py313pl5321h5ca1c30_2 \
quast.py /data/flye/assembly.fasta \
--fungus \
-r /data/referenceData/GCA_900106115.1_CBS_141442_assembly_genomic.fna \
-g /data/referenceData/GCA_900106115.1_CBS_141442_assembly_genomic.gff \
-o quast_output
# note that reference data is from Sungouiella intermedia https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_900106115.1/

#funannotate
#mask
# docker run --rm -v "$PWD":/data -w /data --platform linux/amd64 nextgenusfs/funannotate funannotate mask -i flye/assembly.fasta -o masked_assembly.fasta  
#Predict
# docker run --rm -v "$PWD":/data -w /data --platform linux/amd64 nextgenusfs/funannotate funannotate predict -i masked_assembly.fasta -o /data/output -s 'Sungouiella akabanensis'

flye scaffolding
raccoon scaffolding <- de novo scafffolding?? 
medaka <- polish genome

#interpro
docker run --rm \
-v $PWD/interproscan-5.76-107.0/data:/opt/interproscan/data \
-v $PWD/input:/input \
-v $PWD/temp:/temp \
-v $PWD/output:/output \
--platform linux/amd64 \
interpro/interproscan:5.76-107.0 \
--input /input/Sungouiella_akabanensis.proteins.fa \
--output-dir /output \
--tempdir /temp \
--cpu 8

