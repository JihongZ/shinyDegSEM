# setwd("gene_app24aug/")
source("BMC_functions0925.R")
source("global.R")

# Data Read-in ------------------------------------------------------------
data <- read.table(file = "MS_entrez_id_alldata.txt", header = TRUE, sep = "\t", row.names = 1)
rownames(data) <- paste("ENTREZID", rownames(data), sep = ":")
seqtmp <- readLines("label.txt")
DEGs_table <- read.csv("DEGs_table.csv")
seq_ <- seqfile <- as.numeric(unlist(strsplit(seqtmp, " ")))
group <- seq_ - 1 # start with 0


## SAM Analysis ----------------------------------------------------------

data_ <- data
data_sam <- list(x = as.matrix(data_), y = seq_, 
                 geneid = as.character(1:nrow(data_)), 
                 genenames = paste("g", as.character(1:nrow(data_)), sep = ""), 
                 logged2 = TRUE)

### SAMR: Significance analysis of microarrays ---------------------------

samr.obj <- samr(data_sam, resp.type = "Two class unpaired", nperms = 400)
delta.table <- samr.compute.delta.table(samr.obj)

delta <- 1 # source: [0, 10, .1]
siggenes.table <- samr.compute.siggenes.table(samr.obj, 
                                              delta, 
                                              data_sam, 
                                              delta.table, 
                                              min.foldchange = 2)
gene_up <- siggenes.table$genes.up # Endpoint: output$deg_analysis_output2
gene_lo <- siggenes.table$genes.lo # Endpoint: output$deg_analysis_output3
pos_sig <- as.numeric(c(gene_up[, 3], gene_lo[, 3])) # Endpoint: output$deg_analysis_output4
pos_sig1 <- as.vector(DEGs_table[, 2])
data_sign <- data[pos_sig, ]
fold_change_sign <- foldchange_log2(data_sign, seq_, TRUE)

# Endpoint: output$deg_analysis_output1
samr.plot(samr.obj, del = delta)          

# Source: input$EnrichOrTable, choices = c("enrich" = "enrich", "tab" = "tab")
enrich_input_type <- "enrich" 


# run SPIA  ---------------------------------------------------------------

prepareSPIA(kegg, "out_kegg2", print.names = FALSE)
r_graph <- runSPIA(de = fold_change_sign, all = rownames(data_), 
                   pathwaySetName = "out_kegg2")
r_graph # I have issue running r_graph

options <- as.vector(r_graph[, 1])
selected_options <- options[1:5]


# Network-SEM Conversation ------------------------------------------------------------

i <- selected_options[1]
pathinfo <- labeledge("kegg", i)
union_path <- merge_pathway("kegg", i)
graph_path <- as(union_path, "graphNEL")
sg_path <- subGraph(intersect(rownames(data), graph_path@nodes), graph_path)
outp_path <- graph_DE(graph_path, sg_path, names(fold_change_sign), "out")

## vcol
vcol1 <- rep("c", length(outp_path[[2]]))
vcol1[find_positions(outp_path[[2]], names(fold_change_sign))] <- "green"
vcol1[-find_positions(outp_path[[2]], names(fold_change_sign))] <- "yellow"

## Data Preparation for SEM
pos_gene_path <- find_positions(rownames(data), outp_path[[2]])
data_pathway <- data[pos_gene_path, ]
data_pathway <- t(data_pathway)
colnames(data_pathway) <- sub("ENTREZID:", "G", colnames(data_pathway))
data_pathway <- cbind(group, data_pathway)
colnames(data_pathway)[1] <- "group"
data_pathway <- as.data.frame(data_pathway)

# SEM Analysis ------------------------------------------------------------

## input$selected_pathway: Select Pathway for SEM
selected_pathway <- selected_options[1]
select_option <- selected_pathway 
i <- select_option


union_path_sem <- merge_pathway("kegg", i)
graph_path_sem <- as(union_path, "graphNEL")
sg_path_sem <- subGraph(intersect(rownames(data), graph_path@nodes), graph_path)
outp_path_sem <- graph_DE(graph_path, sg_path, names(fold_change_sign), "out")


pathinfo_sem <- labeledge("kegg", i)
edges0 <- getedge(outp_path_sem[[1]], pathinfo_sem)
edge0_global <- (edges0)
mod0 <- originalmodel(edges0)
cat(mod0)
pos_gene_path_sem <- find_positions(rownames(data), outp_path_sem[[2]])

data_pathway_sem <- data[pos_gene_path_sem, ]
data_pathway_sem <- t(data_pathway_sem)
colnames(data_pathway_sem) <- sub("ENTREZID:", "G", colnames(data_pathway_sem))
data_pathway_sem <- cbind(group, data_pathway_sem)
colnames(data_pathway_sem)[1] <- "group"
data_pathway_sem <- as.data.frame(data_pathway_sem)
data_pathway_sem


## Initial Modle
data_pathway_sem <- data_pathway_sem_global <- data_pathway_sem
mod0 <- originalmodel(edges0)

fit0 <- sem(mod0,
             data = data_pathway_sem,
             estimator = "ML",
             fixed.x = FALSE, std.ov = F
)
# 对最新的model进行modification indices计算
MI_tbl <- modindices(fit0, minimum.value = 3.84, sort = TRUE)

## edge analysis ---------------------------------------------------------
MITable <- MITablefinal <- NULL

cov_data_pathway_list <-
  lapply(split(data_pathway_sem, 
               data_pathway_sem[["group"]]), 
         function(x) cov(x[, setdiff(colnames(x), "group")]))
all_path <- as.matrix(Matrix::bdiag(cov_data_pathway_list))
colnames(all_path) <- rownames(all_path) <- unlist(lapply(names(cov_data_pathway_list), 
                                                          \(x) paste0(colnames(cov_data_pathway_list[[x]]), "_", x)))
moddiff <- edgemodel(edges0 = edge0_global, MIchoice = MITable)
cat(moddiff)

fitdiff <-
  sem(
    moddiff,
    estimator = "ML",
    sample.cov = all_path,
    sample.nobs = mean(table(data_pathway_sem[["group"]])),
    fixed.x = FALSE
  )
summary(fitdiff)
