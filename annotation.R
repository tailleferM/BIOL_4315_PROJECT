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



Sakabanensis <- read_excel('data/out.emapper.annotations.xlsx')
#set column names
colnames(Sakabanensis) <- as.character(Sakabanensis[2,])
#remove empty rows
Sakabanensis <- Sakabanensis[-c(1,2),]

table(Sakabanensis$COG_category)


# Compute counts
counts <- as.data.frame(table(Sakabanensis$COG_category))

# Write CSV
write.csv(counts, file = 'DebaryomycesHansenii.csv', row.names = FALSE)

######### READING EVERYONE'S COG CATEGORY COUNTS ################

# Set your folder path
folder <- "data/everyoneCOGs"

# List CSV files
files <- list.files(folder, pattern = "\\.csv$", full.names = TRUE)

# Read all files into named lists
data_list <- lapply(files, function(file) {
  df <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)
  # named vector: variable = value
  setNames(df[[2]], df[[1]])
})

# Row names = file base names
names(data_list) <- tools::file_path_sans_ext(basename(files))

# ---- Determine full set of variables across all files ----
all_vars <- unique(unlist(lapply(data_list, names)))

# ---- Build a clean dataframe with matching columns ----
final_df <- data.frame(matrix(NA, nrow = length(data_list), ncol = length(all_vars)),
                       stringsAsFactors = FALSE)
colnames(final_df) <- all_vars
rownames(final_df) <- names(data_list)

# Fill each row with the variable values for that file
for (i in seq_along(data_list)) {
  vec <- data_list[[i]]
  final_df[i, names(vec)] <- vec
}

#replace NA with 0
final_df[is.na(final_df)] <- 0

#give first column a name

colnames(final_df)[1] <- 'None'

# View result
final_df

# add a sum row
sum_row <- colSums(final_df)

COGsums <- rbind(final_df, Sum = sum_row) 


### splitting multi cog categories. ###

# df is your original dataframe

df <- COGsums
# Identify single-letter columns
single_cols <- grep("^[a-zA-Z]$", names(df), value = TRUE)

# Remove "None" if it's in the list (just in case)
single_cols <- setdiff(single_cols, "None")

# Identify multi-letter columns
multi_cols  <- grep("^[a-zA-Z]{2,}$", names(df), value = TRUE)

# Remove "None" from multi-letter list (in case it appears)
multi_cols <- setdiff(multi_cols, "None")

# Create output dataframe with single-letter columns only
df_out <- df[single_cols]

# For each multi-letter column, distribute values to corresponding letters
for (mc in multi_cols) {
  
  letters_in_col <- strsplit(mc, "")[[1]]
  
  # keep only letters that exist as single-letter columns
  letters_in_col <- letters_in_col[letters_in_col %in% single_cols]
  
  # add values
  for (ltr in letters_in_col) {
    df_out[[ltr]] <- df_out[[ltr]] + df[[mc]]
  }
}

df_out

#renaming categories

names(df_out) <- dict[names(df_out)]

df_out$None <- COGsums$None


new_df <- df_out[, !(names(df_out) %in% "Function unknown")]
df3 <- new_df[, !(names(new_df) %in% "None")]


#### plotting results ####

# Extract the Sum row as a named numeric vector
sum_vec <- df3["Sum", ]

# Convert to a dataframe suitable for ggplot
plot_df <- data.frame(
  variable = names(sum_vec),
  value = as.numeric(sum_vec)
)

fermentation_related <- c(
  "Energy production and conversion",
  "Carbohydrate transport and metabolism",
  "Amino acid transport and metabolism",
  "Coenzyme transport and metabolism"
)

resistance_related <- c(
  "Defense mechanisms",
  "Cell wall/membrane/envelope biogenesis",
  "Inorganic ion transport and metabolism",
  "Secondary metabolites biosynthesis, transport and catabolism"
)

both_related <- c(
  "Lipid transport and metabolism",
  "Signal transduction mechanisms",
  "Posttranslational modification, protein turnover, chaperones"
)

plot_df$category <- case_when(
  plot_df$variable %in% fermentation_related ~ "Fermentation",
  plot_df$variable %in% resistance_related ~ "Resistance",
  plot_df$variable %in% both_related ~ "Both",
  TRUE ~ "Other"
)

plot_df$variable <- factor(
  plot_df$variable,
  levels = plot_df$variable[order(plot_df$value)]
)

ggplot(plot_df, aes(x = variable, y = value, fill = category)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c(
    "Fermentation" = "#1b9e77",
    "Resistance" = "#d95f02",
    "Both" = "#7570b3",
    "Other" = "grey70"
  )) +
  labs(
    title = "Sum of Annotated Genes",
    x = "COG Category",
    y = "Sum",
    fill = "Category"
  )
