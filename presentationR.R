# reading everyone's eggNOG-mapper output

library(dplyr)
library(ggplot2)
library(scales)
library(tidyr)
library(purrr)
library(cowplot)
library(KEGGREST)
library(GO.db)
library(AnnotationDbi)

############ READING FILES / CONSTRUCTING LIST OF DATAFRAMES ###############

# get files
emapper_files <- list.files(path = 'data/', pattern = "\\.xlsx$", full.names = TRUE)

emapper_list <- list()

for (file in emapper_files) {
  
  # Extract the base filename
  fname <- basename(file)
  
  # Remove the emapper suffix
  obj_name <- sub("\\.emapper\\.annotations\\.xlsx$", "", fname)
  
  # Read file
  df <- readxl::read_xlsx(file)
  
  #set column names
  colnames(df) <- as.character(df[2,])
  
  #remove empty rows
  df <- df[-c(1,2),]
  
  # Store in list
  emapper_list[[obj_name]] <- df
}

############ MAPPING KOGS ############

#build COG dictionary

COGs <- readLines('COGs.txt')
# use a regex to extract key and value
keys <- sub("^\\[([A-Z])\\].*$", "\\1", COGs)              # extract the letter inside brackets
values <- sub("^\\[[A-Z]\\]\\s*(.*)$", "\\1", COGs)        # extract everything after the bracket

# create a named vector (acts like a dictionary)
dict <- setNames(values, keys)

# map COGs
for (name in names(emapper_list)) {
  emapper_list[[name]] <- emapper_list[[name]] %>%
    mutate(
      COGs = sapply(COG_category, function(x) {
        letters <- strsplit(x, "")[[1]] #split multi COGs  
        mapped  <- dict[letters]
        paste(mapped, collapse = "; ")
      })
    )
}

########### FILTERING AND PLOTTING ALCOHOL DEHYDROGENASE ############
format_strain <- function(name) {
  genus_initial <- substr(name, 1, 1)
  species <- substring(name, 2)
  species <- tolower(species)
  paste0(genus_initial, ". ", species)
}

geneCounts <- data.frame(
  strain = sapply(names(emapper_list), format_strain),
  ADH_count = sapply(names(emapper_list), function(name) {
    df <- emapper_list[[name]]
    sum(grepl("^ADH[0-9]+$", df$Preferred_name, ignore.case = FALSE))
  })
)

#plot alcohol dehydrogenase counts
ggplot(geneCounts, aes(y = ADH_count, x = strain)) +
  geom_bar(stat = "identity", fill = "skyblue") + 
  labs(title = "Number of Alcohol Dehydrogenase Copies", x = "Yeast Strain", y = "ADH copies") +
  scale_y_continuous(breaks = pretty_breaks())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9))

########## plotting CYP counts ############

geneCounts$P450_count <- sapply(names(emapper_list), function(name) {
    df <- emapper_list[[name]]
    sum(grepl("P-?450", df$Description, ignore.case = TRUE))
  })

#plot cytochrome p450 counts
ggplot(geneCounts, aes(y = P450_count, x = strain)) +
  geom_bar(stat = "identity", fill = "darkred") + 
  labs(title = "Number of Cytochrome P450s", x = "Yeast Strain", y = "P450 count") +
  scale_y_continuous(breaks = pretty_breaks())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9))

####### heat shock protein ########

geneCounts$heatShock_count <- sapply(names(emapper_list), function(name) {
  df <- emapper_list[[name]]
  sum(grepl("Heat Shock Protein", df$Description, ignore.case = TRUE))
})

#plot heat shock protein counts
ggplot(geneCounts, aes(y = heatShock_count, x = strain)) +
  geom_bar(stat = "identity", fill = "skyblue") + 
  labs(title = "Number of Heat Shock Proteins", x = "Yeast Strain", y = "Heat Shock Protein count") +
  scale_y_continuous(breaks = pretty_breaks())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9))


##########plot comparison

# Reshape and plot
gene_counts_long <- pivot_longer(geneCounts, 
                                 cols = c(ADH_count, P450_count),
                                 names_to = "Gene_Type",
                                 values_to = "Count")
                                 
# adh p450

ggplot(gene_counts_long, aes(x = strain, y = Count, fill = Gene_Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Gene Copy Numbers by Yeast Strain", 
       x = "Yeast Strain", 
       y = "Gene Copies",
       fill = "Gene Type") +
  scale_fill_manual(values = c("ADH_count" = "skyblue", "P450_count" = "coral"),
                    labels = c("ADH", "P450")) +
  scale_y_continuous(breaks = pretty_breaks()) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Reshape and plot
gene_counts_long <- pivot_longer(geneCounts, 
                                 cols = c(ADH_count, heatShock_count),
                                 names_to = "Gene_Type",
                                 values_to = "Count")

# adh p450

ggplot(gene_counts_long, aes(x = strain, y = Count, fill = Gene_Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Gene Copy Numbers by Yeast Strain", 
       x = "Yeast Strain", 
       y = "Gene Copies",
       fill = "Gene Type") +
  scale_fill_manual(values = c("ADH_count" = "skyblue", "heatShock_count" = "coral"),
                    labels = c("ADH", "HeatShockProtein")) +
  scale_y_continuous(breaks = pretty_breaks()) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

########## plotting overall COG counts ###########


count_cog_letters <- function(df) {
  
  # Split e.g. "ABC" → c("A","B","C")
  letters <- unlist(strsplit(df$COG_category, ""))
  
  # Count occurrences
  tbl <- table(letters)
  
  # Convert to a 1-row dataframe
  as.data.frame(as.list(tbl))
}

# Run for all strains
cog_letter_df <- lapply(emapper_list, count_cog_letters) %>%
  bind_rows(.id = "strain")

# Replace missing with 0
cog_letter_df[is.na(cog_letter_df)] <- 0

letter_cols <- setdiff(colnames(cog_letter_df), "strain")
colnames(cog_letter_df)[match(letter_cols, colnames(cog_letter_df))] <- dict[letter_cols]
colnames(cog_letter_df)[2] <- 'None'
row.names(cog_letter_df) <- cog_letter_df[,1]
COGcounts_df <- cog_letter_df[,-1]


# Define the order of categories by group
category_order <- c(
  # INFORMATION STORAGE AND PROCESSING - Blues
  "Translation, ribosomal structure and biogenesis",
  "RNA processing and modification",
  "Transcription",
  "Replication, recombination and repair",
  "Chromatin structure and dynamics",
  
  # CELLULAR PROCESSES AND SIGNALING - Greens
  "Cell cycle control, cell division, chromosome partitioning",
  "Nuclear structure",
  "Defense mechanisms",
  "Signal transduction mechanisms",
  "Cell wall/membrane/envelope biogenesis",
  "Cell motility",
  "Cytoskeleton",
  "Extracellular structures",
  "Intracellular trafficking, secretion, and vesicular transport",
  "Posttranslational modification, protein turnover, chaperones",
  
  # METABOLISM - Oranges/Reds
  "Energy production and conversion",
  "Carbohydrate transport and metabolism",
  "Amino acid transport and metabolism",
  "Nucleotide transport and metabolism",
  "Coenzyme transport and metabolism",
  "Lipid transport and metabolism",
  "Inorganic ion transport and metabolism",
  "Secondary metabolites biosynthesis, transport and catabolism",
  
  # POORLY CHARACTERIZED - Purples
  "General function prediction only",
  "Function unknown",
  
  # None
  "None"
)

# Colorblind-friendly palette
cog_colors <- c(
  # None - Grey
  "None" = "#666666",
  
  # INFORMATION STORAGE AND PROCESSING - Blues (deuteranopia/protanopia safe)
  "Translation, ribosomal structure and biogenesis" = "#004488",
  "RNA processing and modification" = "#1965B0",
  "Transcription" = "#4E79A7",
  "Replication, recombination and repair" = "#7AAAD0",
  "Chromatin structure and dynamics" = "#A6CEE3",
  
  # CELLULAR PROCESSES AND SIGNALING - Greens/Teals
  "Cell cycle control, cell division, chromosome partitioning" = "#005A32",
  "Nuclear structure" = "#238B45",
  "Defense mechanisms" = "#41AB5D",
  "Signal transduction mechanisms" = "#78C679",
  "Cell wall/membrane/envelope biogenesis" = "#ADDD8E",
  "Cell motility" = "#006837",
  "Cytoskeleton" = "#31A354",
  "Extracellular structures" = "#74C476",
  "Intracellular trafficking, secretion, and vesicular transport" = "#A1D99B",
  "Posttranslational modification, protein turnover, chaperones" = "#00441B",
  
  # METABOLISM - Oranges/Vermillions (colorblind distinguishable from greens/blues)
  "Energy production and conversion" = "#882255",
  "Carbohydrate transport and metabolism" = "#CC6677",
  "Amino acid transport and metabolism" = "#DD8855",
  "Nucleotide transport and metabolism" = "#DDAA33",
  "Coenzyme transport and metabolism" = "#E8AD45",
  "Lipid transport and metabolism" = "#F0C75E",
  "Inorganic ion transport and metabolism" = "#F4D88D",
  "Secondary metabolites biosynthesis, transport and catabolism" = "#F7E5B2",
  
  # POORLY CHARACTERIZED - Purples/Magentas
  "General function prediction only" = "#AA4499",
  "Function unknown" = "#CC99BB"
)

# Reshape data and set factor order
cog_long <- cog_letter_df %>%
  mutate(strain = format_strain(strain)) %>%  # Format strain names
  pivot_longer(cols = -strain, 
               names_to = "COG_category", 
               values_to = "count") %>%
  mutate(COG_category = factor(COG_category, levels = category_order))

# Create plot WITH legend (for extracting legend)
plot_with_legend <- ggplot(cog_long, aes(x = strain, y = count, fill = COG_category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cog_colors) +
  labs(x = "Strain", 
       y = "Count", 
       fill = "KOG Category",
       title = "KOG Category Distribution by Strain") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9))

# Create plot WITHOUT legend
plot_no_legend <- plot_with_legend + 
  theme(legend.position = "none")

# Extract and save the legend
legend <- get_legend(plot_with_legend)

# Save the plot without legend
ggsave("plots/cog_stacked_barplot.png", plot_no_legend, 
       width = 10, height = 6, dpi = 300)

# Save the legend separately using ggdraw
legend_plot <- ggdraw(legend)

# Display the plot without legend
print(plot_no_legend)

############ MAPPING KEGGS ############

# Helper: split comma-separated KEGG ID fields
split_ids <- function(x) {
  if (is.na(x) || x == "" || x == "-") return(character(0))
  unlist(strsplit(x, ","))
}

# Download all KEGG dictionaries
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

ko_dict      <- clean_dict(ko_dict,      prefix_pattern = "^ko:")
pathway_dict <- clean_dict(pathway_dict, prefix_pattern = "^(path:|ko:|map:)")
module_dict  <- clean_dict(module_dict)
reaction_dict <- clean_dict(reaction_dict)
rclass_dict   <- clean_dict(rclass_dict)
brite_dict    <- clean_dict(brite_dict)

# Vectorized lookup function
lookup_many <- function(vec, dict) {
  if (length(vec) == 0) return(NA)
  matches <- dict[vec]
  matches <- matches[!is.na(matches)]
  if (length(matches) == 0) return(NA)
  paste(matches, collapse = "; ")
}

# PART 4: Apply KEGG annotation to all dataframes

for (name in names(emapper_list)) {
  message(paste0("Processing ", name, "..."))
  
  emapper_list[[name]] <- emapper_list[[name]] %>%
    # Clean KEGG_ko column (remove "ko:" prefix if present)
    mutate(KEGG_ko = gsub("^ko:", "", KEGG_ko)) %>%
    # Add KEGG annotations
    mutate(
      KO_name = map_chr(KEGG_ko,      ~ lookup_many(split_ids(.x), ko_dict)),
      Pathway_name = map_chr(KEGG_Pathway, ~ lookup_many(split_ids(.x), pathway_dict)),
      Module_name = map_chr(KEGG_Module,   ~ lookup_many(split_ids(.x), module_dict)),
      Reaction_name = map_chr(KEGG_Reaction, ~ lookup_many(split_ids(.x), reaction_dict)),
      Rclass_name = map_chr(KEGG_rclass, ~ lookup_many(split_ids(.x), rclass_dict)),
      BRITE_name = map_chr(BRITE, ~ lookup_many(split_ids(.x), brite_dict))
    )
}

############ mapping go terms #############

# Add GO descriptions column to each dataframe in emapper_list
emapper_list <- lapply(emapper_list, function(df) {
  # Extract all unique GO terms from the GOs column
  # Assuming GO terms are separated by commas or semicolons
  all_go_terms <- unique(unlist(strsplit(as.character(df$GOs), "[,;]")))
  all_go_terms <- trimws(all_go_terms)
  all_go_terms <- all_go_terms[all_go_terms != "" & all_go_terms != "-" & !is.na(all_go_terms)]
  
  # Get GO term descriptions
  go_lookup <- AnnotationDbi::select(GO.db,
                                     keys = all_go_terms,
                                     columns = c("TERM"),
                                     keytype = "GOID")
  
  # Create a function to map GO IDs to terms for each row
  df$GO_descriptions <- sapply(df$GOs, function(go_string) {
    if (is.na(go_string) || go_string == "" || go_string == "-") {
      return(NA)
    }
    go_terms <- unlist(strsplit(as.character(go_string), "[,;]"))
    go_terms <- trimws(go_terms)
    
    # Look up each term
    terms <- sapply(go_terms, function(go_id) {
      match_idx <- which(go_lookup$GOID == go_id)
      if (length(match_idx) > 0) {
        return(go_lookup$TERM[match_idx[1]])
      } else {
        return(go_id)
      }
    })
    
    paste(terms, collapse = "; ")
  })
  
  return(df)
})


################## Plotting completeness ############

assemblyMetrics <- read.csv('data/assemblyMetrics.csv', row.names = 1)

strain_mapping <- c(
  "Dhansenii" = "Debaryomyces hansenii",
  "Lfermentati" = "Lachancea fermentati",
  "Lthermotolerans" = "Lachancea thermotolerans",
  "Mguilliermondi" = "Meyerozyma guilliermondii",
  "Nholstii" = "Nakazawaea holstii",
  "Npopuli" = "Nakazawaea populi",
  "Sakabanensis" = "Sungouiella akabanensis",
  "Sbayanus" = "Saccharomyces bayanus",
  "Sparadoxus" = "Saccharomyces paradoxus",
  "Wanomalus" = "Wickerhamomyces anomalus",
  "Zpseurodouxii" = "Zygosaccharomyces pseudorouxii"
)

assemblyMetrics$predictedGenes <- sapply(rownames(assemblyMetrics), function(sp) {
  # Find matching short name
  short_name <- names(strain_mapping)[strain_mapping == sp]
  if (length(short_name) > 0 && short_name %in% names(emapper_list)) {
    return(nrow(emapper_list[[short_name]]))
  } else {
    return(NA)
  }
})

#putting it in the COG count df
COGcounts_df$predictedgenes <- sapply(emapper_list, nrow)

COGcounts_df$annotated <- COGcounts_df$predictedgenes - COGcounts_df$None 
COGcounts_df$annotated100 <- COGcounts_df$annotated/COGcounts_df$predictedgenes
COGcounts_df$functions <- COGcounts_df$annotated - COGcounts_df$`Function unknown` 
COGcounts_df$functions100 <- COGcounts_df$functions / COGcounts_df$predictedgenes

COGsums <- colSums(COGcounts_df)

# Convert to data frame for ggplot
COGsums_df <- data.frame(
  COG = names(COGsums),
  Count = as.numeric(COGsums)
)

COGsums_df <- COGsums_df[-c(27,28,29,30,31,26),]
COGsums_df$percent <- 100*(COGsums_df$Count / sum(COGsums_df$Count))
COGsums_df <- COGsums_df[-c(19),]
COGsums_df <- COGsums_df[-1,]

# Create barplotlibrary(scales)

ggplot(COGsums_df, aes(x = COG, y = percent)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "KOG Category Distribution", 
       x = "KOG Category", 
       y = "Percentage") +
  scale_y_continuous(labels = percent_format(scale = 1)) +  # if Percentage is already 0-100
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

########### heat resistance #############

# Conservative heat resistance GO terms - direct heat response only
heat_resistance_GO <- c(
  # Direct heat response
  "GO:0009408",  # response to heat
  "GO:0034605",  # cellular response to heat
  "GO:0070370",  # cellular heat acclimation
  "GO:0010286",  # heat acclimation
  
  # Heat shock proteins (core)
  "GO:0006986",  # response to unfolded protein
  "GO:0042026",  # protein refolding
  
  # Trehalose (specific to heat protection in yeast)
  "GO:0005992"   # trehalose biosynthetic process
)

# Dataframe with descriptions
heat_GO_terms <- data.frame(
  GO_ID = c("GO:0009408", "GO:0034605", "GO:0070370", "GO:0010286",
            "GO:0006986", "GO:0042026", "GO:0005992"),
  Description = c("response to heat", 
                  "cellular response to heat", 
                  "cellular heat acclimation", 
                  "heat acclimation",
                  "response to unfolded protein", 
                  "protein refolding", 
                  "trehalose biosynthetic process")
)

geneCounts$heatShock_count <- sapply(names(emapper_list), function(name) {
  df <- emapper_list[[name]]
  sum(grepl("Heat Shock Protein", df$Description, ignore.case = TRUE))
})

# Count heat resistance genes per strain
geneCounts$heatTolerance <- sapply(emapper_list, function(df) {
    sum(grepl(paste(heat_resistance_GO, collapse = "|"), df$GOs), na.rm = TRUE)
  })

geneCounts$totalgenes <- COGcounts_df$predictedgenes
geneCounts$heatovertotal <- geneCounts$heatTolerance/geneCounts$totalgenes*100

# Create bar plot
ggplot(geneCounts, aes(x = strain, y = heatovertotal)) +
  geom_bar(stat = "identity", fill = "coral") +
  labs(title = "Genome Fraction of Heat Resistance Genes by Yeast Strain",
       x = "Yeast Strain",
       y = "Heat Resistance Genes fraction") +
  scale_y_continuous(labels = percent_format(scale = 1)) +  # if Percentage is already 0-100
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

######## Cold resistance ##############

# Conservative cold resistance GO terms - direct cold response
# Most conservative - only direct cold response
cold_resistance_GO <- c(
  "GO:0009409",  # response to cold
  "GO:0070417",  # cellular response to cold
  "GO:0009631",  # cold acclimation
  "GO:0050826"   # response to freezing
)

# Dataframe with descriptions
cold_GO_terms <- data.frame(
  GO_ID = c("GO:0009409", "GO:0070417", "GO:0009631", "GO:0050826",
            "GO:0042542", "GO:0005992", "GO:0046526"),
  Description = c("response to cold",
                  "cellular response to cold",
                  "cold acclimation",
                  "response to freezing",
                  "response to hydrogen peroxide",
                  "trehalose biosynthetic process",
                  "D-gluconate catabolic process")
)

# Count unique genes with cold resistance GO terms per strain
cold_gene_counts <- data.frame(
  strain = sapply(names(emapper_list), format_strain),
  cold_gene_count = sapply(emapper_list, function(df) {
    has_cold_GO <- grepl(paste(cold_resistance_GO, collapse = "|"), df$GOs)
    sum(has_cold_GO, na.rm = TRUE)
  })
)

geneCounts$coldResistance <- cold_gene_counts$cold_gene_count

# Plot
ggplot(cold_gene_counts, aes(x = strain, y = cold_gene_count)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  labs(title = "Cold Resistance Genes by Yeast Strain",
       x = "Yeast Strain",
       y = "Number of Cold Resistance Genes") +
  scale_y_continuous(breaks = pretty_breaks()) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# hot cold plot 
# Select only temperature-related genes
temp_counts <- geneCounts[, c("strain", "heatTolerance", "coldResistance")]

temp_counts_long <- pivot_longer(temp_counts,
                                 cols = c(heatTolerance, coldResistance),
                                 names_to = "Gene_Type",
                                 values_to = "Count")

ggplot(temp_counts_long, aes(x = strain, y = Count, fill = Gene_Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Temperature Resistance Genes by Yeast Strain",
       x = "Yeast Strain",
       y = "Gene Count",
       fill = "Temperature Response") +
  scale_fill_manual(values = c("heatTolerance" = "orangered",
                               "coldResistance" = "steelblue"),
                    labels = c("Cold Resistance",'Heat Tolerance')) +
  scale_y_continuous(breaks = pretty_breaks()) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")

########### lactic acid fermentation ##########

geneCounts$laccase <- sapply(names(emapper_list), function(name) {
  df <- emapper_list[[name]]
  sum(grepl("L-lactate dehydrogenase", df$KO_name, ignore.case = TRUE))
})

#plot lac counts
ggplot(geneCounts, aes(y = laccase, x = strain)) +
  geom_bar(stat = "identity", fill = "orchid4") + 
  labs(title = "L-Lactate Dehydrogenase Hits", x = "Yeast Strain", y = "potential LDH count") +
  scale_y_continuous(breaks = pretty_breaks())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9))

######### osmotolerance ##########

geneCounts$osmotolerance <- sapply(names(emapper_list), function(name) {
  df <- emapper_list[[name]]
  sum(grepl("hyperosmotic salinity", df$GO_descriptions, ignore.case = TRUE))
})

#plot hyperosmotic counts
ggplot(geneCounts, aes(y = osmotolerance, x = strain)) +
  geom_bar(stat = "identity", fill = "lightblue") + 
  labs(title = "Hyperosmotic Salinity Response Genes", x = "Yeast Strain", y = "Hyperosmotic response hits")+
  scale_y_continuous(breaks = pretty_breaks())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9))

######### CYP P450 ############

geneCounts$CYP <- sapply(names(emapper_list), function(name) {
  df <- emapper_list[[name]]
  sum(grepl("P450", df$GO_descriptions, ignore.case = TRUE))
})

#plot hyperosmotic counts
ggplot(geneCounts, aes(y = osmotolerance, x = strain)) +
  geom_bar(stat = "identity", fill = "lightblue") + 
  labs(title = "Hyperosmotic Salinity Response Genes", x = "Yeast Strain", y = "Hyperosmotic response hits")+
  scale_y_continuous(breaks = pretty_breaks())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9))


############# Beta-Lyase (for thiol beer) ##############


geneCounts$bLyase <- sapply(names(emapper_list), function(name) {
  df <- emapper_list[[name]]
  sum(grepl("beta-lyase", df$Description, ignore.case = TRUE))
})

#plot beta lyase counts
ggplot(geneCounts, aes(y = bLyase, x = strain)) +
  geom_bar(stat = "identity", fill = "lightblue") + 
  labs(title = "Beta Lyase hits", x = "Yeast Strain", y = "copy count")+
  scale_y_continuous(breaks = pretty_breaks())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9))


geneCounts$STR3 <- sapply(names(emapper_list), function(name) {
  df <- emapper_list[[name]]
  sum(grepl("STR3", df$Preferred_name, ignore.case = TRUE))
})

#plot beta lyase counts
ggplot(geneCounts, aes(y = STR3, x = strain)) +
  geom_bar(stat = "identity", fill = "lightblue") + 
  labs(title = "STR3 hits", x = "Yeast Strain", y = "copy count")+
  scale_y_continuous(breaks = pretty_breaks())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9))


geneCounts$IRC7 <- sapply(names(emapper_list), function(name) {
  df <- emapper_list[[name]]
  sum(grepl("IRC7", df$Preferred_name, ignore.case = TRUE))
})

#plot beta lyase counts
ggplot(geneCounts, aes(y = IRC7, x = strain)) +
  geom_bar(stat = "identity", fill = "lightblue") + 
  labs(title = "IRC7 hits", x = "Yeast Strain", y = "copy count")+
  scale_y_continuous(breaks = pretty_breaks())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(size = 9))


