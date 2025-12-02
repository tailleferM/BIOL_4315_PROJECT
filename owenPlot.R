library(GO.db)
library(AnnotationDbi)
library(dplyr)
library(stringr)
library(tidyr)
library(readxl)
library(ggplot2)
library(scales)  # for scales::percent

akabanensis <- read_xlsx("data/out.emapper.annotations.xlsx", skip = 2)
anomalus  <- read_xlsx("data/MM_jpced9_r.emapper.annotations.xlsx", skip = 2)

process_go <- function(df) {
  df <- as.data.frame(df)
  names(df) <- str_trim(names(df))
  go_col <- names(df)[str_detect(names(df), regex("^GO", ignore_case = TRUE))][1]
  gene_col <- names(df)[str_detect(names(df), regex("^query$", ignore_case = TRUE))][1]
  names(df)[names(df) == go_col] <- "GO_IDs"
  names(df)[names(df) == gene_col] <- "Gene"
  
  df <- df[, c("Gene", "GO_IDs"), drop = FALSE]
  df$GO_IDs <- as.character(df$GO_IDs)
  df <- tidyr::separate_rows(df, GO_IDs, sep = ",|;|\\|")
  df <- df[!is.na(df$GO_IDs) & df$GO_IDs != "", ]
  
  keys <- unique(df$GO_IDs)
  keys <- as.character(keys)
  
  go_info <- AnnotationDbi::select(x = GO.db,
                                   keys = keys,
                                   columns = c("GOID", "TERM", "ONTOLOGY"),
                                   keytype = "GOID")
  
  df <- merge(df, go_info, by.x = "GO_IDs", by.y = "GOID", all.x = TRUE)
  df$ONTOLOGY <- dplyr::case_when(
    df$ONTOLOGY == "MF" ~ "Molecular Function",
    df$ONTOLOGY == "CC" ~ "Cellular Component",
    df$ONTOLOGY == "BP" ~ "Biological Process",
    TRUE ~ NA_character_
  )
  
  return(df)
}

akabanensis_long <- process_go(akabanensis)
anomalus_long  <- process_go(anomalus)

# Combine both species for joint analysis
combined_long <- bind_rows(
  akabanensis_long %>% mutate(Species = "S. akabanensis"),
  anomalus_long  %>% mutate(Species = "W. anomalus")
)

# Keep only relevant ontologies
combined_long <- combined_long %>%
  filter(ONTOLOGY %in% c("Biological Process", "Molecular Function"))

# Define keyword groups
sugar_keywords <- c("glycolysis", "mannose", "carbon metabolism")
osmotic_keywords <- c("osmotic", "turgor", "salt stress", "osmoregulation", "osmolarity", "HOG")

# Categorize GO terms based on TERM text
combined_long <- combined_long %>%
  mutate(
    Pathway = case_when(
      str_detect(str_to_lower(TERM), str_c(sugar_keywords, collapse = "|"))   ~ "sugar",
      str_detect(str_to_lower(TERM), str_c(osmotic_keywords, collapse = "|")) ~ "osmotic",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Pathway))

summary_counts <- combined_long %>%
  group_by(Species, ONTOLOGY, Pathway, TERM) %>%
  summarise(Gene_Count = n(), .groups = "drop") %>%
  group_by(Species, Pathway) %>%
  mutate(Percent = Gene_Count / sum(Gene_Count))

ggplot(summary_counts, aes(x = Pathway, y = TERM, size = Percent, color = Pathway)) +
  geom_point(alpha = 0.7) +
  facet_grid(Pathway ~ Species, scales = "free") +
  scale_size_continuous(range = c(3, 12), labels = scales::percent) +
  # 🔥 color labels now match actual Pathway values ("sugar", "osmotic")
  scale_color_manual(values = c("sugar" = "#0072B2", "osmotic" = "#D55E00")) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  labs(
    title = "Sugar and Osmotic Pathways in GO Terms",
    x = NULL,
    y = "GO Term"
  )
