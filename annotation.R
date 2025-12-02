# Analysis of annotation
# mapping gene ontology numbers to their associated term and ontology

#Load libraries
library(GO.db)
library(AnnotationDbi)
library(dplyr)
library(stringr)
library(readxl)
library(tidyr)
library(ggplot2)
library(reshape2)
library(scales)
library(knitr)
# read eggnogdbmapper table

Wanomalus <- read_excel('data/MM_jpced9_r.emapper.annotations.xlsx')
Sakabanensis <- read_excel('data/out.emapper.annotations.xlsx')
#set column names
colnames(Sakabanensis) <- as.character(Sakabanensis[2,])
#remove empty rows
Sakabanensis <- Sakabanensis[-c(1,2),]

COGs <- readLines('COGs.txt')
# use a regex to extract key and value
keys <- sub("^\\[([A-Z])\\].*$", "\\1", COGs)              # extract the letter inside brackets
values <- sub("^\\[[A-Z]\\]\\s*(.*)$", "\\1", COGs)        # extract everything after the bracket

# create a named vector (acts like a dictionary)
dict <- setNames(values, keys)

Sakabanensis <- Sakabanensis %>%
  mutate(COGs = recode(COG_category, !!!dict))

######### GO TABLE ##########

#Convert to long table
GoTable <- Sakabanensis %>%
  separate_rows(GOs, sep = ",")

#Get the terms and the ontology

GoTable$TERM <- AnnotationDbi::select(GO.db,
                                      keys = GoTable$GOs,
                                      columns = c("TERM"),
                                      keytype = "GOID")$TERM

GoTable$ONTOLOGY <- AnnotationDbi::select(GO.db,
                                          keys = GoTable$GOs,
                                          columns = c("ONTOLOGY"),
                                          keytype = "GOID")$ONTOLOGY

GoTable$DEFINITION <- AnnotationDbi::select(GO.db,
                                            keys = GoTable$GOs,
                                            columns = c("DEFINITION"),
                                            keytype = "GOID")$DEFINITION

####### KEGG TABLE ##########
library(KEGGREST)
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)

#-------------------------------------------------------
# 0. Helper: split comma-separated KEGG ID fields
#-------------------------------------------------------
split_ids <- function(x) {
  if (is.na(x) || x == "" || x == "-") return(character(0))
  unlist(strsplit(x, ","))
}

#-------------------------------------------------------
# 1. Download all KEGG dictionaries in one shot (FAST)
#-------------------------------------------------------

message("Downloading KEGG dictionaries...")

ko_dict <- keggList("ko")
pathway_dict <- keggList("pathway")
module_dict <- keggList("module")
reaction_dict <- keggList("reaction")
rclass_dict <- keggList("rclass")
brite_dict <- keggList("brite")

# Clean dictionaries → each becomes a named vector
clean_dict <- function(x, prefix_pattern = NULL) {
  keys <- names(x)
  vals <- as.character(x)
  if (!is.null(prefix_pattern)) {
    keys <- gsub(prefix_pattern, "", keys)
  }
  setNames(vals, keys)
}

Sakabanensis$KEGG_ko <- gsub("^ko:", "", Sakabanensis$KEGG_ko)

ko_dict      <- clean_dict(ko_dict,      prefix_pattern = "^ko:")
pathway_dict <- clean_dict(pathway_dict, prefix_pattern = "^(path:|ko:|map:)")
module_dict  <- clean_dict(module_dict)
reaction_dict <- clean_dict(reaction_dict)
rclass_dict   <- clean_dict(rclass_dict)
brite_dict    <- clean_dict(brite_dict)

#-------------------------------------------------------
# 2. Vectorized lookup (SUPER FAST)
#-------------------------------------------------------
lookup_many <- function(vec, dict) {
  if (length(vec) == 0) return(NA)
  matches <- dict[vec]
  matches <- matches[!is.na(matches)]
  if (length(matches) == 0) return(NA)
  paste(matches, collapse = "; ")
}

#-------------------------------------------------------
# 3. Apply lookup to your dataframe
#-------------------------------------------------------
message("Annotating your data...")

KeggTableW <- Wanomalus %>%
  mutate(
    KO_name = map_chr(KEGG_ko,      ~ lookup_many(split_ids(.x), ko_dict)),
    Pathway_name = map_chr(KEGG_Pathway, ~ lookup_many(split_ids(.x), pathway_dict)),
    Module_name = map_chr(KEGG_Module,   ~ lookup_many(split_ids(.x), module_dict)),
    Reaction_name = map_chr(KEGG_Reaction, ~ lookup_many(split_ids(.x), reaction_dict)),
    Rclass_name = map_chr(KEGG_rclass, ~ lookup_many(split_ids(.x), rclass_dict)),
    BRITE_name = map_chr(BRITE, ~ lookup_many(split_ids(.x), brite_dict))
  )

message("Done! 🚀 KEGG annotation completed.")

#-------------------------------------------------------
# 4. Preview output
#-------------------------------------------------------
head(annotated)

############### STUFF TO PARSE FOR ##############

# enzymes for bioremediation
# breaking down organic compounds
# according to:     Laccase: A potential biocatalyst for pollutant degradation 
# Peroxidases-ligninolytic peroxidases, laccase, lignin peroxidase, manganese peroxidase horseradish peroxidase (HRP), tyrosinases, esterases, nitrilases, aminohydrolases, lipase, cutinase, and organophosphorus hydrolase, dehalogenases, cytochrome p450 monooxygenase, multicopper oxidase, ascorbate oxidase, 

# xenobiotic metabolism enzymes
# according to: recent advances in fungal xenobiotic metabolism: enzymes and applications
# CYPs, peroxidases, laccases, tyrosinases and unspecific peroxygenase, 

#UPOs unspecific peroxidases : Identification and Expression of New Unspecific Peroxygenases – Recent Advances, Challenges and Opportunities

# metal boremediation stuff: Designing yeast as plant-like hyperaccumulators for heavy metals
#transporters can be divalent metal transporters, permeases, exporters, or have auxiliary metal transport fxn
# metal transporters: ZRT1, ZRT2, CTR1, CTR3, FTR1, FET4, SMF1, SMF2. 
# phosphate transporters can transport arsenate: PHOs 84, 87, 89 
# sulfate transporters can transport chromate: SUL2
# vacuole transporters: CCC1, COT1, ZRC1, SMF3
# Yeast optimizes metal utilization based on metabolic network and enzyme kinetics
# iron-sulfur cluster (ISC) synthesis reactions


# GO terms for metal resistance according to Evolution of cross-tolerance to metals in yeast
#“fungal-type vacuole membrane,” “incipient cellular bud site,” and “vacuolar transporter chaperone complex” among cellular components; 
#“response to stimulus,” “biological regulation,” “cytokinetic process,” and “polyphosphate biosynthetic process” among biological processes; 
#“phosphotransferase activity,” “ubiquitin-protein transferase activity,” and “MAP-kinase scaffold activity” among molecular functions.
# read paper for association between over-representation of thing and tolerance of stuff
# vacuole transporter chaperone genes: VTCs

# GO terms for metal resistance: chemical-genomic profiling identifies gene that protect yeast from aluminium, gallium and indium toxicity
# read this paper for more GOs 

##### PLOTTING COGS n Stuff #########

filtered <- Sakabanensis %>% filter(nchar(COGs) > 4)
metal_related <- c(
  "Inorganic ion transport and metabolism",
  "Posttranslational modification, protein turnover, chaperones",
  'Intracellular trafficking, secretion and vesicular transport',
  'Secondary metabolites biosynthesis, transport and catabolism'
)
filtered$metal_flag <- ifelse(filtered$COGs %in% metal_related, "Metal-related", "Other")
filtered <- filtered %>% filter(COGs != 'Function unknown')


ggplot(filtered, aes(x = COGs)) +
  geom_bar() +
  theme_minimal() +
  xlab("Value") +
  ylab("Count") +
  ggtitle("Frequency of Shared Values")

ggplot(filtered, aes(x = COGs)) +
  geom_bar(fill = "steelblue") +
  coord_flip() +
  theme_minimal() +
  xlab("Value") +
  ylab("Count") +
  ggtitle("COG categories")

ggplot(filtered, aes(y = COGs, fill = metal_flag)) +
  geom_bar() +
  scale_fill_manual(values = c("Metal-related" = "red", "Other" = "grey70")) +
  labs(
    x = "Count",
    y = "COG category",
    fill = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
   # legend.position = "top"
  )

### ============================
### BIOREMEDIATION GENE FILTERING PIPELINE
### ============================

### ---- 1. Keyword lists ----

xeno_keywords <- c(
  "laccase", "peroxidase", "oxidase", "monooxygenase",
  "cytochrome P450", "CYP", "dehydrogenase", "oxidoreductase",
  "glutathione", "GST", "xenobiotic"
)

metal_keywords <- c(
  "ferric", "ferroxidase", "iron", "copper", "zinc",
  "permease", "reductase", "SOD", "superoxide dismutase",
  "metallothionein", "metal transporter", "ATP-binding cassette"
)

### ---- 2. Combine all annotation columns into one searchable string ----

annotation_cols <- c("Description", "Preferred_name", "GOs", 
                     "KEGG_ko", "KEGG_Pathway", "PFAMs", "COG_category")

Sakabanensis$all_annotations <- apply(Sakabanensis[, annotation_cols], 1, paste, collapse = " ")

### ---- 3. FILTERING ----

df <- Sakabanensis

## 3A. Xenobiotic metabolism genes
xeno_hits <- df[grepl(paste(xeno_keywords, collapse="|"),
                      df$all_annotations, ignore.case = TRUE), ]
xeno_hits$bioremediation_category <- "xenobiotic_metabolism"

## 3B. Heavy metal uptake genes
metal_hits <- df[grepl(paste(metal_keywords, collapse="|"),
                       df$all_annotations, ignore.case = TRUE), ]
metal_hits$bioremediation_category <- "metal_uptake"

## 3C. COG category = P (inorganic ion transport/metabolism)
#cog_p <- df[df$COG_category == "P", ]
#cog_p$bioremediation_category <- "COG_P_inorganic_ion_metabolism"


### ---- 4. Combine and deduplicate ----

combined_hits <- do.call(rbind, list(xeno_hits, metal_hits))

# remove identical "query" entries
combined_hits <- combined_hits[!duplicated(combined_hits$query), ]


### ---- 5. Summary barplot ----

library(ggplot2)

ggplot(combined_hits, aes(bioremediation_category)) +
  geom_bar() +
  theme_bw() +
  labs(
    title = "Bioremediation-Related Genes",
    x = "Category",
    y = "Number of Genes"
  )


### ---- 6. OPTIONAL: print summary counts ----

print(table(combined_hits$bioremediation_category))

### ---- DONE ----

alcohol <- data.frame(
  strain = c('W. anomalus', 'S. akabanensis'),
  adh = c(9, 5),
  cyp = c(1,9),
  multicopper = c(6,2),
  
  
)

ggplot(alcohol, aes(y = adh, x = strain)) +
  geom_bar(stat = "identity", fill = "skyblue") + 
  labs(title = "Number of Alcohol Dehydrogenase Copies", x = "Yeast Strain", y = "ADH copies") +
  scale_y_continuous(breaks = pretty_breaks())
       
