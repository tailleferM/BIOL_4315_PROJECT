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

# read eggnogdbmapper table

#Wanomalus <- read_excel('Wickerhamomyces/funannotate/output/eggnog_results/MM_jpced9_r.emapper.annotations.xlsx')
Sakabanensis <- read_excel('Sungouiella/eggnog/out.emapper.annotations.xlsx')
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

#Convert to long table
W_long <- Sakabanensis %>%
  separate_rows(GOs, sep = ",")

#Get the terms and the ontology

W_long$TERM <- AnnotationDbi::select(GO.db,
                                      keys = W_long$GOs,
                                      columns = c("TERM"),
                                      keytype = "GOID")$TERM

W_long$ONTOLOGY <- AnnotationDbi::select(GO.db,
                                          keys = W_long$GOs,
                                          columns = c("ONTOLOGY"),
                                          keytype = "GOID")$ONTOLOGY

W_long$DEFINITION <- AnnotationDbi::select(GO.db,
                                            keys = W_long$GOs,
                                            columns = c("DEFINITION"),
                                            keytype = "GOID")$DEFINITION


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



Sakabanensis <- read_excel('Sungouiella/eggnog/out.emapper.annotations.xlsx')
#set column names
colnames(Sakabanensis) <- as.character(Sakabanensis[2,])
#remove empty rows
Sakabanensis <- Sakabanensis[-c(1,2),]

table(Sakabanensis$COG_category)


# Compute counts
counts <- as.data.frame(table(Sakabanensis$COG_category))

# Write CSV
write.csv(counts, file = 'myNAMEHERECOGs.csv', row.names = FALSE)


