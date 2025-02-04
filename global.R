## 2024/09/07 Update: Load all packages after installation
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
library(BiocManager)
if (!requireNamespace("Biostrings", quietly = TRUE)){
  BiocManager::install("Biostrings")
}
library(Biostrings)
if (!requireNamespace("graphite", quietly = TRUE)){
  BiocManager::install("graphite")
}
library(graphite)
if (!requireNamespace("ggm", quietly = TRUE)){
  BiocManager::install("ggm")
}
library(ggm)
if (!requireNamespace("impute", quietly = TRUE)){
  BiocManager::install("impute")
}
library(impute)
if (!requireNamespace("igraph", quietly = TRUE)){
  install.packages("igraph")
}
library(igraph)
if (!requireNamespace("lavaan", quietly = TRUE)){
  install.packages("lavaan")
}
library(lavaan)
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)){
  BiocManager::install("org.Hs.eg.db")
}
library(org.Hs.eg.db)
if (!requireNamespace("openxlsx", quietly = TRUE)){
  install.packages("openxlsx")
}
library(openxlsx)
if (!requireNamespace("samr", quietly = TRUE)){
  install.packages("samr")
}
library(samr)
if (!requireNamespace("semPlot", quietly = TRUE)){
  install.packages("semPlot")
}
library(semPlot)
if (!requireNamespace("semTools", quietly = TRUE)){
  install.packages("semTools")
}
library(semTools)
if (!requireNamespace("shiny", quietly = TRUE)){
  install.packages("shiny")
}
library(shiny)
if (!requireNamespace("SPIA", quietly = TRUE)){
  BiocManager::install("SPIA")
}
library(SPIA)
if (!requireNamespace("stringr", quietly = TRUE)){
  BiocManager::install("stringr")
}
library(stringr)
if (!requireNamespace("SEMgraph", quietly = TRUE)){
  install.packages("SEMgraph", type = "binary")
}
library(SEMgraph)


# library(org.Hs.eg.db)
#setwd("...")

database <- "kegg"
specieslist <- graphite::pathwayDatabases()[graphite::pathwayDatabases()[,2]=='kegg',1] #all species 
speciesname <- "hsapiens"
kegg <- graphite::pathways(speciesname, "kegg")
kegg <- graphite::convertIdentifiers(kegg, "ENTREZID")
