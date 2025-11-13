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

Wanomalus <- read_excel('data/MM_jpced9_r.emapper.annotations.xlsx')
#set column names
colnames(Wanomalus) <- as.character(Wanomalus[2,])
#remove empty rows
Wanomalus <- Wanomalus[-c(1,2),]

#Convert to long table
W_long <- Wanomalus %>%
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

COGs <- readLines('data/COGcategories')
# use a regex to extract key and value
keys <- sub("^\\[([A-Z])\\].*$", "\\1", COGs)              # extract the letter inside brackets
values <- sub("^\\[[A-Z]\\]\\s*(.*)$", "\\1", COGs)        # extract everything after the bracket

# create a named vector (acts like a dictionary)
dict <- setNames(values, keys)

W_long <- W_long %>%
  mutate(COGs = recode(COG_category, !!!dict))



