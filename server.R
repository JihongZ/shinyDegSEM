source("BMC_functions0925.R")
source("global.R")
options(shiny.maxRequestSize = 100 * 1024^2)

function(input, output, session) {
  data_pathway_sem_global <- reactiveVal(NULL)
  mod0 <- reactiveValues()

  data <- reactive({
    file <- input$data_file
    if (is.null(file)) {
      return(NULL)
    }
    # data <- read.table(file = "MS_entrez_id_alldata.txt", header = TRUE, sep = "\t", row.names = 1)
    data <- read.table(file$datapath, header = TRUE, sep="\t", row.names = 1)
    if (is.null(data)) {
      showNotification("Error reading data file. Please check the file format.", type = "error")
      return(NULL)
    }
    rownames(data) <- paste("ENTREZID", rownames(data), sep = ":")
    return(data)
  })

  input_type <- reactive({
    input_type <- input$fileOrText
  })

  seqtext <- reactive({
    text <- input$case_sequence
    req(text)
    seqtext <- isolate(as.numeric(unlist(strsplit(input$case_sequence, " "))))
  })

  seqfile <- reactive({
    # seqtmp <- readLines("label.txt")
    file1 <- input$sequence_file
    req(file1)
    seqtmp <- readLines(file1$datapath)
    seqfile <- isolate(as.numeric(unlist(strsplit(seqtmp, " "))))
  })

  seq_ <- reactive({
    if (input_type() == "text") {
      seq_ <- seqtext()
    } else {
      seq_ <- seqfile()
    }
  })

  group <- reactive({
    group <- seq_() - 1
  })

  output$input_data1 <- renderText({
    paste("choose species: ", input$species)
  })
  output$input_data2 <- renderText({
    paste("choose data file: ", input$data_file$name)
  })
  output$input_data3 <- renderText({
    paste("choose label: from", input$fileOrText)
  })
  output$input_data4 <- renderText({
    paste(if (input$fileOrText == "text") {
      paste("label: ", input$case_sequence)
    } else {
      paste("label file: ", input$sequence_file$name)
    })
  })
  
  output$input_data5 <- renderText({
    paste("row number: ", nrow(data()))
  })
  output$input_data6 <- renderText({
    paste("column number: ", ncol(data()))
  })
  
  output$previewdata <- renderTable({
    data()[1:20,]
  })
  

  deg_input_type <- reactive({
    deg_input_type <- input$SAMOrTable
  })

  delta <- reactive({
    delta <- input$delta_value
  })

  data_sam <- reactive({
    data_ <- data()
    if (is.null(data_)) {
      return(NULL)
    }
    data_sam <- list(x = as.matrix(data_), y = seq_(), geneid = as.character(1:nrow(data_)), genenames = paste("g", as.character(1:nrow(data_)), sep = ""), logged2 = TRUE)
  })
  
  # DEG Analysis (02/04/2025) ----
  samr.obj <- reactive({
    data_sam_ <- data_sam()
    if (is.null(data_sam_)) return(NULL)
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Running SAM", value = 0)
    
    #try1
    samr.obj <- samr(data_sam_, resp.type = "Two class unpaired", nperms = 400)
    # load("D:\\phd\\project\\gSEM2311\\review\\shinyDegSEM\\shinyDegSEM\\samr.obj")
      
    progress$inc(1, detail = "Finish")
    samr.obj
  })

  delta.table <- reactive({
    samr.obj_ <- samr.obj()
    if (is.null(samr.obj_)) {
      return(NULL)
    }
    # Create a Progress object
    progress <- shiny::Progress$new()
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())
    
    progress$set(message = "Computing delta table", value = 0)
    progress$inc(1/5, detail = "Running")
    
    # try2
    delta.table <- samr.compute.delta.table(samr.obj_)
    # load("D:\\phd\\project\\gSEM2311\\review\\shinyDegSEM\\shinyDegSEM\\delta.table")
    
    progress$inc(4/5, detail = "Finish")
    
    delta.table
  })

  siggenes.table <- reactive({
    samr.obj_ <- samr.obj()
    delta.table_ <- delta.table()
    if (is.null(samr.obj_) || is.null(delta.table_)) {
      return(NULL)
    }
    
    # Create a Progress object
    progress <- shiny::Progress$new()
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())
    
    progress$set(message = "Computing significant genes table", value = 0)
    progress$inc(1/5, detail = "Running")
    
    siggenes.table <- samr.compute.siggenes.table(samr.obj_, delta(), data_sam(), delta.table_, min.foldchange = 2)
    
    progress$inc(4/5, detail = "Finish")
    
    siggenes.table
  })

  gene_up <- reactive({
    siggenes.table_ <- siggenes.table()
    if (is.null(siggenes.table_)) {
      return(NULL)
    }
    gene_up <- siggenes.table_$genes.up
  })

  gene_lo <- reactive({
    siggenes.table_ <- siggenes.table()
    if (is.null(siggenes.table_)) {
      return(NULL)
    }
    gene_lo <- siggenes.table_$genes.lo
  })

  pos_sig <- reactive({
    pos_sig <- as.numeric(c(gene_up()[, 3], gene_lo()[, 3]))
  })

  pos_sig1 <- reactive({
    file2 <- input$DEGs_table
    req(file2)
    pos_sigtmp <- read.csv(file2$datapath)
    pos_sig1 <- as.vector(pos_sigtmp[, 2])
  })

  pos_sig_ <- reactive({
    if (deg_input_type() == "SAM") {
      pos_sig_ <- pos_sig()
    } else {
      pos_sig_ <- pos_sig1()
    }
  })

  data_sign <- reactive({
    data_sign <- data()[pos_sig_(), ]
  })

  fold_change_sign <- reactive({
    data_sign_ <- data_sign()
    if (is.null(data_sign_)) {
      return(NULL)
    }
    fold_change_sign <- foldchange_log2(data_sign_, seq_(), TRUE)
  })

  output$deg_analysis_output1 <- renderPlot({
    samr.obj_ <- samr.obj()
    if (is.null(samr.obj_)) {
      return(NULL)
    }
    samr.plot(samr.obj_, del = delta())
  })

  output$deg_analysis_output2 <- renderTable({
    gene_up()
  })

  output$deg_analysis_output3 <- renderTable({
    gene_lo()
  })

  # output$deg_analysis_output4 <- renderText({
  #   pos_sig_()
  # })

  output$deg_input <- renderText({
    paste("DEGs file uploaded: ", input$DEGs_table$name)
  })

  output$download_DEGs_table <- downloadHandler(
    filename = function() {
      "DEGs_table.csv"
    },
    content = function(file) {
      write.csv(data_sign(), file)
    }
  )

  enrich_input_type <- reactive({
    enrich_input_type <- input$EnrichOrTable
  })

  # try3
  prepareSPIA(kegg, "out_kegg2")
  
  # Enrichment Analysis (02/04/2025)----
  r_graph <- reactive({
    fold_change_sign_ <- fold_change_sign()
    data_ <- data()
    if (is.null(fold_change_sign_) || is.null(data_)) {
      return(NULL)
    }
    
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    
    progress$set(message = "Preparing pathway dataset", value = 0)
    progress$inc(1/5, detail = "Running SPIA analysis")
    
    # try4
    r_graph <- runSPIA(fold_change_sign_, rownames(data_), "out_kegg2")
    # load("D:\\phd\\project\\gSEM2311\\review\\shinyDegSEM\\shinyDegSEM\\r_graph")
    
    r_graph <- r_graph[r_graph$NDE > 2, ]  # Filter out pathways with NDE <= 2
    
    progress$inc(4/5, detail = "Finish! You can run Network Analysis now.")
    
    r_graph
  })

  r_graph1 <- reactive({
    file3 <- input$pathway_table
    req(file3)
    r_graphtmp <- read.csv(file3$datapath)
    r_graph1 <- r_graphtmp[, -1]
  })

  r_graph_ <- reactive({
    if (enrich_input_type() == "enrich") {
      r_graph_ <- r_graph()
    } else {
      r_graph_ <- r_graph1()
    }
    
    return(r_graph_)
  })

  # observeEvent(input$pathway_upd_indicate, {
  observe({
    options <- as.vector(r_graph_()[, 1])
    # selected_options <- options[1:5]
    selected_options <- options[1]
    updateSelectInput(session, "pathway_select", choices = options, selected = selected_options)
    updateSelectInput(session, "selected_pathway", choices = options, selected = selected_options[1])
  })

  output$pathway_ori <- renderTable({
    r_graph_()
  })

  output$enrich_input <- renderText({
    paste("enrichment file uploaded: ", input$pathway_table$name)
  })

  output$download_rgraph <- downloadHandler(
    filename = function() {
      "ig.RData"
    },
    content = function(file) {
      write.csv(r_graph_(), file)
    }
  )

  output$network_analysis <- renderUI({
    select_options <- input$pathway_select

    tabs <- lapply(select_options, function(opt) {
      tabPanel(
        opt,
        verbatimTextOutput(paste0("network_original_", opt)),
        # strong("Click 'Run Network Analysis' on the left panels of Network Analysis, choose  pathways for network analysis."),
        # p(" "),
        # strong("For the next step, complete the left panels of SEM Analysis, then click the button 'Initial Run'."),
        h3("Original Model"),
        h4("Average Degree:"),
        textOutput(paste0("network_original_avdegree_", opt)),
        h4("Density:"),
        textOutput(paste0("network_original_density_", opt)),
        h3("Final Model"),
        h4("Network Graph:"),
        plotOutput(paste0("network_graph_", opt)),
        h4("Average Degree:"),
        textOutput(paste0("network_minpath_avdegree_", opt)),
        h4("Density:"),
        textOutput(paste0("network_minpath_density_", opt)),
        h3("Pathway data fed into SEM analysis"),
        tableOutput(paste0("data_for_SEM_", opt))
      )
    })
    do.call(tabsetPanel, tabs)
  })

  # Network analysis through looping
  process_pathway <- function(i) {
    local({
      pathinfo <- reactive({
        labeledge("kegg", i)
      })

      union_path <- reactive({
        merge_pathway("kegg", i)
      })

      pathinfo <- reactive({
        labeledge("kegg", i)
      })

      union_path <- reactive({
        merge_pathway("kegg", i)
      })
      pathinfo <- reactive({
        pathinfo <- labeledge("kegg", i)
      })

      union_path <- reactive({
        union_path <- merge_pathway("kegg", i)
      })

      graph_path <- reactive({
        graph_path <- as(union_path(), "graphNEL")
      })

      sg_path <- reactive({
        sg_path <- subGraph(intersect(rownames(data()), graph_path()@nodes), graph_path())
      })

      outp_path <- reactive({
        outp_path <- graph_DE(graph_path(), sg_path(), names(fold_change_sign()), "out")
      })

      ig <- reactive({
        ig <- outp_path()[[1]]
        V(ig)$name <- gsub("ENTREZID:", "", V(ig)$name)
        V(ig)$label <- mapIds(org.Hs.eg.db, V(ig)$name, column = 'SYMBOL', keytype = 'ENTREZID')
        return(ig)
      })

      vcol1 <- reactive({
        vcol1 <- rep("c", length(outp_path()[[2]]))
        vcol1[find_positions(outp_path()[[2]], names(fold_change_sign()))] <- "green"
        vcol1[-find_positions(outp_path()[[2]], names(fold_change_sign()))] <- "yellow"
        return(vcol1)
      })



      pos_gene_path <- reactive({
        pos_gene_path <- find_positions(rownames(data()), outp_path()[[2]])
      })

      data_pathway <- reactive({
        data_pathway <- data()[pos_gene_path(), ]
        data_pathway <- t(data_pathway)
        colnames(data_pathway) <- sub("ENTREZID:", "G", colnames(data_pathway))
        data_pathway <- cbind(group(), data_pathway)
        colnames(data_pathway)[1] <- "group"
        data_pathway <- as.data.frame(data_pathway)
        return(data_pathway)
      })



      # 定义 output 对象时，使用 local() 函数创建局部环境
      local({
        output[[paste0("network_original_", i)]] <- renderPrint({
          graph_path()
        })

        output[[paste0("network_graph_", i)]] <- renderPlot({
          # plot(outp_path()[[1]],
          #      vertex.label = outp_path()[[2]],
          #      vertex.color = vcol1(), vertex.size = 5, vertex.label.cex = .7, vertex.label.dist = 1, edge.arrow.size = .4, edge.curved = .1,
          #      main = "Plot"
          # )
          
          plot(ig(),
            vertex.label = V(ig())$label,
            vertex.color = vcol1(), vertex.size = 5, vertex.label.cex = .7, vertex.label.dist = 1, edge.arrow.size = .4, edge.curved = .1,
            main = "Plot"
          )
        })

        output[[paste0("vcol1_test_", i)]] <- renderText({
          vcol1()
        })

        output[[paste0("network_minpath_avdegree_", i)]] <- renderText({
          mean(igraph::degree(outp_path()[[1]]))
        })

        output[[paste0("network_minpath_density_", i)]] <- renderText({
          igraph::graph.density(outp_path()[[1]], loops = TRUE)
        })
        output[[paste0("network_original_", i)]] <- renderPrint({
          graph_path()
        })

        output[[paste0("network_original_avdegree_", i)]] <- renderText({
          mean(degree(igraph.from.graphNEL(graph_path())))
        })

        output[[paste0("network_original_density_", i)]] <- renderText({
          graph.density(igraph.from.graphNEL(graph_path()), loops = TRUE)
        })

        output[[paste0("network_microarray_avdegree_", i)]] <- renderText({
          mean(degree(igraph.from.graphNEL(sg_path())))
        })

        output[[paste0("network_microarray_density_", i)]] <- renderText({
          graph.density(igraph.from.graphNEL(sg_path()), loops = TRUE)
        })



        output[[paste0("data_for_SEM_", i)]] <- renderTable({
          data_pathway()
        })

        output$download_SEM_data <- downloadHandler(
          filename = function() {
            "data_for_SEM.xlsx"
          },
          content = function(file) {
            write.xlsx(data_pathway(), file)
          }
        )
        
        
        # Download Pathway Igraph: ----
        output$download_pathway_igraph <- downloadHandler(
          filename = function() {
            "pathway_igraph.rds"
          },
          content = function(file) {
            saveRDS(ig(), file)
          }
        )
        
      })
    })
  }
  observe({
    select_options <- as.vector(input$pathway_select)
    for (i in select_options) {
      # 调用函数，处理单个 pathway
      process_pathway(i)
    }
  })

  # SEM analysis using the same logic like Network analysis
  data_pathway_sem <- reactiveVal(NULL)
  # 1. construct function like process_pathway() but change to SEM; we call it process_sem
  process_sem <- function(i) {
    local({
      pathinfo_sem <- reactive({
        labeledge("kegg", i)
      })

      union_path_sem <- reactive({
        merge_pathway("kegg", i)
      })

      graph_path_sem <- reactive({
        as(union_path_sem(), "graphNEL")
      })

      sg_path_sem <- reactive({
        subGraph(intersect(rownames(data()), graph_path_sem()@nodes), graph_path_sem())
      })

      outp_path_sem <- reactive({
        graph_DE(graph_path_sem(), sg_path_sem(), names(fold_change_sign()), "out")
      })

      edges0 <- reactive({
        edges0 <- getedge(outp_path_sem()[[1]], pathinfo_sem())
        # edge0_global <<- (edges0)
        # edges0 <- getedge(outp_path_sem()[[1]], pathinfo_sem())
      })

      ############# key steps#########
      #edge0_global <<- reactive({
      #  getedge(outp_path_sem()[[1]], pathinfo_sem())
      #})

      # Create Initial SEM model: mod0 (Removed)---- 
      mod0 <<- reactive({
        originalmodel(edges0())
      })
      
      
      ############# /key steps#########

      pos_gene_path_sem <- reactive({
        pos_gene_path_sem <- find_positions(rownames(data()), outp_path_sem()[[2]])
      })

      data_pathway_sem <- reactive({
        data_pathway_sem <- data()[pos_gene_path_sem(), ]
        data_pathway_sem <- t(data_pathway_sem)
        colnames(data_pathway_sem) <- sub("ENTREZID:", "G", colnames(data_pathway_sem))
        data_pathway_sem <- cbind(group(), data_pathway_sem)
        colnames(data_pathway_sem)[1] <- "group"
        data_pathway_sem <- as.data.frame(data_pathway_sem)
        return(data_pathway_sem)
      })


      data_pathway_sem_global(data_pathway_sem)


      edge0_global_global <<- reactive({
        edges0()
      })


      output$SEM_data <- renderTable({
        data_pathway_sem()
      })

      output$SEM_networkSummary <- renderTable({
        edges0()
      })
    })
  }

  # 2.It has to be observe() select_option
  observeEvent(input$Iinitialrun, {
    select_option <- as.vector(input$selected_pathway)
    process_sem(select_option)
  })


  # ... (其他代码)
  # estimator <- reactive({
  #   input$Iestimator
  # })



  # # Download Pathway Igraph: ----
  # output$download_pathway_igraph <- downloadHandler(
  #   filename = function() {
  #     "pathway_igraph.rds"
  #   },
  #   content = function(file) {
  #     saveRDS(outp_path()[[1]], file)
  #   }
  # )

  # SEM: Output Tabset ----------------------------------------------------
  output$sem_analysis <- renderUI({
    tabs <- list(
      id = "panel1_sem1",
      conditionalPanel(
        condition = !is.null(input$selected_pathway),
        tableOutput("all_model_showcase"),
        h3("Edge List of the Latest Model"),
        tableOutput("edge0_global_display")
      ),
      tabPanel("Step1:Original Model",
        value = "step1",
        # strong("Click the button 'Initial Run' on the left panel"),
        # p(" "),
        verbatimTextOutput("Omod0smry"),
        verbatimTextOutput("Omod0fit")
      ),
      tabPanel("Step2:Modification Indices and Final Model",
        value = "step2",
        tabsetPanel(
          id = "panel1_sem1_2",
          tabPanel("Latest Model",
            value = "step2m",
            # verbatimTextOutput("OmodifiedModfit"),
            # verbatimTextOutput("Omod1fit"),
            verbatimTextOutput("Omod1smry")
          ),
          tabPanel("Final Model (refreshed after selecting your final model)",
            value = "step2f",
            verbatimTextOutput("Omodfinalfit"),
            verbatimTextOutput("Omodfinalsmry"),
            h4("SEM Plot"),
            plotOutput("Plotfinalfit"),
            h4("Parameter Estiamtes"),
            tableOutput("model_parameters")
          )
        )
      ),
      tabPanel("Step3:Model Invariance",
        value = "step3",
        # strong("Choose a model and click the button 'Run Measurement Invariance' on the left panel."),
        # p(" "),
        h4("Model Invariance"),
        verbatimTextOutput("Omodcompr")
      ),
      tabPanel("Step4:Node Analysis",
        value = "step4",
        # strong("Click the button 'Run Node Analysis' on the left panel."),
        # p(" "),
        verbatimTextOutput("Omodnodeanalysis"),
        verbatimTextOutput("Omodnodeanalysisfit"),
        plotOutput("OmodnodeanalysisGraph")
      ),
      tabPanel("Step5:Edge Analysis",
        value = "step5",
        # strong("Click the button 'Run Edge Analysis' on the left panel."),
        # p(" "),
        verbatimTextOutput("Omodedgeanalysis"),
        verbatimTextOutput("Omodedgeanalysisfit"),
        plotOutput("OmodedgeanalysisGraph")
      )
    )
    do.call(tabsetPanel, tabs)
  })


  # Step 1: Initial Run SEM -------------------------------------------------------------
  observeEvent(input$Iinitialrun, {
    updateTabsetPanel(session, inputId = "panel1", selected = "panel1_sem")
    updateTabsetPanel(session, inputId = "panel1_sem1", selected = "step1")

    
    fitList <<- reactiveVal(value = NA) # initialized to empty
    modList <<- reactiveVal(value = NA) # initialized to empty

    # open tab step 1
    data_pathway_sem <- data_pathway_sem_global()
    
    ## 更新model list展示表格 ----
    data_model_showcase <<- reactive({
      tryCatch(
        data.frame(
          Index = paste0("Model ", 0:(length(fitList()) - 1)),
          Fit_Index = unlist(lapply(fitList(), \(x) paste0(
            paste0(
              c("rmsea: ", "rmsea.pvalue: ", "srmr: "),
              round(fitMeasures(x, c("rmsea", "rmsea.pvalue", "srmr")), 3)
            ),
            collapse = "; "
          ))),
          Added_Pathway = c("NA", MIchoiceList())
        ),
        error = function(e) {
          data.frame(
            Index = paste0("Model", 0:(length(fitList()) - 1)),
            Fit_Index = unlist(lapply(fitList()[length(fitList()) - 1], \(x) paste0(
              paste0(
                c("rmsea: ", "rmsea.pvalue: ", "srmr: "),
                round(fitMeasures(x, c("rmsea", "rmsea.pvalue", "srmr")), 3)
              ),
              collapse = "; "
            ))),
            Added_Pathway = c("NA", MIchoiceList())
          )
        }
      )
    })

    tryCatch(
      {
        # Use SEMgraph Package for SEM Analysis: ----
        # run estimation
        fit0 <<- reactive({
          SEMrun(
            graph = lavaan2graph(mod0()),
            data = data_pathway_sem()
          )$fit
        })

        # UnitTest: 模型运行
        output$Omod0smry <- renderPrint({
          fit0_ <- fit0()
          if (is.null(fit0_)) {
            return(NULL)
          }
          summary(fit0_)
        })
        
        output$edge0_global_display <- renderTable({
          fit0_ <- fit0()
          pathway_tbl <- parameterEstimates(fit0_)
          pathway_tbl
        })
        
        output$Omod0fit <- renderPrint({
          fit0_ <- fit0()
          if (is.null(fit0_)) {
            return(NULL)
          }
          fitMeasures(fit0_, c("rmsea", "rmsea.pvalue", "srmr"))
        })

        # 保存模型组和模型参数
        fitList <<- reactiveVal(list(fit0())) # initialized to original model
        modList <<- reactiveVal(list(mod0())) # initialized to original specs
        MIchoiceList <<- reactiveVal(NULL) # initialized to original modification indices table
        MITableList <<- reactiveVal(NULL) # initialized to original modification indices table
        

        ## 迭代修正原始模型 -----
        ### * 选择修正指数较大的路径添加进模型
        # MI0_tbl <- modindices_new(fitList()[[length(fitList())]], minimum.value = 3.84, sort = TRUE)
        alpha<- 0.05; Q <-qchisq(alpha, df=1, lower.tail=F)
        MI0_tbl <- modindices_new(fitList()[[length(fitList())]], minimum.value =Q, sort = TRUE)
        MI0_tbl_1<- subset(MI0_tbl, MI0_tbl $op == "~")
        MI0_tbl_2<- subset(MI0_tbl, MI0_tbl $op == "~~")
        MI0_tbl<- rbind(MI0_tbl_1, MI0_tbl_2)

        ## A table to store Index, FitIndex, Added Pathway ----
        output$all_model_showcase <- renderTable({
          data_model_showcase()
        })

        ## 生成一个下拉菜单，现在可选的modification indices -----
        output$OModifSel <- renderUI({
          MIchoice0 <- as.character(apply(MI0_tbl[, 1:3], 1, \(x) paste0(x, collapse = " ")))
          selectInput(
            "IModifSel",
            "Select your modification indices to improve your model: ",
            choices = MIchoice0
          )
        })
        
        ## 重置model的下拉菜单 ----
        updateSelectInput(inputId = "ImodifiedModSelNum",  choices = paste0("Model ", 0))
      },
      error = function(e) {
        output$Omod0fit <- renderPrint({
          "Fail to run initial model and modification indices! Try re-run network analysis"
        })
      }
    )
  })

  # Step 2: Update Model Based on MI -----
  observeEvent(input$Irerun, {
    # open tab step 2
    updateTabsetPanel(session, inputId = "panel1", selected = "panel1_sem")
    updateTabsetPanel(session, inputId = "panel1_sem1", selected = "step2")
    updateTabsetPanel(session, inputId = "panel1_sem1_2", selected = "step2m")

    ### Revise Model based on selected path ----
    nmod <<- length(modList())
    MIchoice0 <- input$IModifSel
    mod1 <- paste(modList()[[nmod]], MIchoice0, sep = " \n ")
    data_pathway_sem <- data_pathway_sem_global()
    
    # Fit the updated model ----
    tryCatch(
      {
        fit1 <<- SEMrun(lavaan2graph(mod1),
                        data = data_pathway_sem())$fit
      
        # 对最新的model进行modification indices计算
        # MI_tbl <- modindices_new(fit1, minimum.value = 3.84, sort = TRUE)
        
        alpha<- 0.05; Q <-qchisq(alpha, df=1, lower.tail=F)
        MI_tbl <- modindices_new(fit1, minimum.value =Q, sort = TRUE)
        MI_tbl_1<- subset(MI_tbl, MI_tbl $op == "~")
        MI_tbl_2<- subset(MI_tbl, MI_tbl $op == "~~")
        MI_tbl<- rbind(MI_tbl_1, MI_tbl_2)
        
        
        # 将通路(edges)提取出来
        MIchoice_total <- as.character(apply(MI_tbl[, 1:3], 1, \(x) paste0(x, collapse = " ")))
        
        # 将下拉菜单中去掉选择的这个modification indice
        updateSelectInput(inputId = "IModifSel", choices = MIchoice_total)
        
        ## Output latest modified model
        output$Omod1smry <- renderPrint({
          summary(fit1)
        })
        
        output$edge0_global_display <- renderTable({
          pathway_tbl <- parameterEstimates(fit1) 
          pathway_tbl
        })
        
        output$Omod1fit <- renderPrint({
          fitMeasures(fit1, c("rmsea", "rmsea.pvalue", "srmr"))
        })
        fitList(c(fitList(), fit1)) # 将lavaan fit结果更新到最新的
        modList(c(modList(), mod1)) # 将lavaan model文本更新到最新
        MIchoiceList(c(MIchoiceList(), MIchoice0)) # 将参数列表更新到最新
        MITableList(c(MITableList(), list(MI_tbl))) # 将参数列表更新到最新
        nmod <<- length(modList())
        
        # Dropdown Menu：Select the final model
        output$ImodifiedModSel <- renderUI({
          nmod <<- max(1, nmod - 1)
          selectInput("ImodifiedModSelNum", "Select your final model
                      (Note: Only effective/refreshed after selecting your final model; 
                      If node/edge analysis fail, choose previous model with worst model fit): ", 
                      choices = paste0("Model ", 0:nmod))
        })
      }, warning = function(w) {
        fit1 <<- SEMrun(lavaan2graph(modList()[[nmod]]),
                        data = data_pathway_sem())$fit
        showNotification(paste0("Error: Modified Model Fail! Total number of modified models: ", length(fitList())- 1), type = "error")
        output$Omodcount <- renderText({ # total number of modified models
          print(paste0("Total number of modified models: ", nmod - 1))
        })
        output$OmodifiedModfit <- renderPrint({ # fit indices of modified models
          fitMeasures(fit1, c("rmsea", "rmsea.pvalue", "srmr"))
        })
        
        output$Omod1smry <- renderPrint({
          summary(fit1)
        })
        
        output$edge0_global_display <- renderTable({
          pathway_tbl <- parameterEstimates(fit1)
          pathway_tbl
        })
        
        output$Omod1fit <- renderPrint({
          fitMeasures(fit1, c("rmsea", "rmsea.pvalue", "srmr"))
        })
        output$ImodifiedModSel <- renderUI({
          selectInput("ImodifiedModSelNum",
                      "Select your final model 
                      (Note: if node/edge analysis fail, choose previous model with worst model fit): ",
                      choices = paste0("Model ", 0:(nmod-1))
          )
        })
      })
    
  })

  # Selet the final model after MI -----
  observeEvent(input$ImodifiedModSelNum, {
    # open tab step 2: select final model
    updateTabsetPanel(session, inputId = "panel1", selected = "panel1_sem")
    updateTabsetPanel(session, inputId = "panel1_sem1", selected = "step2")
    updateTabsetPanel(session, inputId = "panel1_sem1_2", selected = "step2f")

    # 选择最终模型的编号
    modfinalIndex <- reactive({
      as.numeric(sub(pattern = "Model ", replacement = "", input$ImodifiedModSelNum))
    })
    modfinal <<- reactive({
      modList()[[modfinalIndex() + 1]]
    })
    fitfinal <<- reactive({
      fitList()[[modfinalIndex() + 1]]
    })
    MIchoicefinal <<- reactive({
      if (modfinalIndex() == 0) {
        NULL
      }else{
       MIchoiceList()[1:modfinalIndex()]
      }
    })
    MITablefinal <<- reactive({
      MITableList()[[modfinalIndex()]]
    })  

    # output selected final model
    output$Omodfinalsmry <- renderPrint({
      summary(fitfinal())
    })
    
    output$Omodfinalfit <- renderPrint({
      fitMeasures(fitfinal(), c("rmsea", "rmsea.pvalue", "srmr"))
    })
    
    output$Plotfinalfit <- renderPlot({
     semPlot::semPaths(fitfinal())
    })
    
    output$final_model <- renderPrint({
      cat(modfinal())
    })
    
    # output$Plotfinalfit <- renderPlot({
    #   # semPlot::semPaths(fitfinal())
    #   SEMgraph::gplot(final_model()$graph)
    # })
    
    output$download_FinalSEM_table <- downloadHandler(
      filename = function() {
        "results_for_final_SEM.xlsx"
      },
      content = function(file) {
        write.xlsx(as.data.frame(parameterestimates(fitfinal())), file)
      }
    )
  })
  
  # Step3: Measurement Invariance ----
  observeEvent(input$IRunMI, {
    # open tab step 4: run modification indices
    updateTabsetPanel(session, inputId = "panel1", selected = "panel1_sem")
    updateTabsetPanel(session, inputId = "panel1_sem1", selected = "step3")
    data_pathway_sem <- data_pathway_sem_global()
    tryCatch(
      {
        fit_base <-
          SEMrun(
            lavaan2graph(modfinal()),
            data = data_pathway_sem()
          )$fit
        fit_node <-
          SEMrun(
            lavaan2graph(modfinal()),
            data = data_pathway_sem(),
            group = data_pathway_sem()$group,
            fit = 1
          )$fit
        fit_edge <-
          SEMrun(
            lavaan2graph(modfinal()),
            data = data_pathway_sem(),
            group = data_pathway_sem()$group,
            fit = 2
            )$fit
        modcompr <- semTools::compareFit(fit_base, 
                                         fit_node, 
                                         fit_edge) # output
        # Output of model comparison
        output$model_table <- renderPrint({
          cat("Model Fit Indices of the Base, Group Effects on Node, and Group Effects on Edge Models\n")
          print(modcompr@fit[, c("npar", "chisq", "df", "pvalue", "cfi", "aic", "bic", "rmsea", "srmr")])
        })
        
        output$anova_results <- renderPrint({
          cat("ANOVA (Base vs. Edge Models) \n")
          print(anova(fit_base, fit_edge))  # Assumes nested models
        })
        
        output$node_fit <- renderPrint({
          cat("Chi-square Goodness-of-Fit Test (for Group Effects on Node Model, i.e., fit_node model) \n")
          print(fit_node)  # Replace with actual chi-square test output
        })
      },
      error = function(e) {
        output$Omodcompr <<- renderPrint({
          "Model Invariance Fails. Choose another model to re-try."
        })
      }
    )
  })

  # Step4: SEM perturbated (UP/DOWN) nodes analysis: nodes ~ group ----
  observeEvent(input$IRunNodeAnalysis, {
    # open tab step 4
    updateTabsetPanel(session, inputId = "panel1", selected = "panel1_sem")
    updateTabsetPanel(session, inputId = "panel1_sem1", selected = "step4")
    data_pathway_sem <- data_pathway_sem_global()
    # run model
    # mod1g <- nodemodel(modfinal(), data_pathway_sem())
    modg <- SEMrun(lavaan2graph(modfinal()),
      data = data_pathway_sem(),
      group = data_pathway_sem()$group
    )
    

    # output of node analysis
    output$Omodnodeanalysis <- renderPrint({modg$gest})
    
    output$Omodnodeanalysisfit <- renderPrint({
      fitMeasures(modg$fit, c("rmsea", "rmsea.pvalue", "srmr"))
    })
    
    output$OmodnodeanalysisGraph <- renderPlot({
      SEMgraph::gplot(modg$graph)
    })
    output$download_NodeAnalysis_table <- downloadHandler(
      filename = function() {
        "results_for_node_analysis.xlsx"
      },
      content = function(file) {
        write.xlsx(as.data.frame(parameterestimates(modg$fit)), file)
      }
    )
  })

  # Step5: SEM perturbated (UP/DOWN) edges analysis: edges ~ group ----
  observeEvent(input$IRunEdgeAnalysis, {
    # open tab step 5
    updateTabsetPanel(session, inputId = "panel1", selected = "panel1_sem")
    updateTabsetPanel(session, inputId = "panel1_sem1", selected = "step5")
    
    data_pathway_sem <- data_pathway_sem_global()
    MITable <- MIchoicefinal()
    
    ### 输入数据data_pathway，以及它的group变量名
    cov_data_pathway_list <-
      lapply(split(data_pathway_sem(), 
                   data_pathway_sem()[["group"]]), 
             function(x) cov(x[, setdiff(colnames(x), "group")]))
    
    all_path <- as.matrix(Matrix::bdiag(cov_data_pathway_list))
    colnames(all_path) <- rownames(all_path) <- unlist(lapply(names(cov_data_pathway_list), \(x) paste0(colnames(cov_data_pathway_list[[x]]), "_", x)))
    
    pathway_tbl <- parameterEstimates(fit0())
    pathway_tbl <- pathway_tbl[,c(1,3,2, 4:ncol(pathway_tbl))]
      
    moddiff <- edgemodel(pathway_tbl, MITable)
    
    ## 我这里将edge_pathway改了，改为组1:_1，组0:_0。而非原来的只有组0改了。
    tryCatch({
      modedge <- SEMrun(lavaan2graph(modfinal()),
                        data = data_pathway_sem(),
                        group = data_pathway_sem()$group,
                        fit = 2
      )
      # output of edge analysis
      output$Omodedgeanalysis <- renderPrint({modedge$dest})
      output$Omodedgeanalysisfit <- renderPrint({
        fitMeasures(modedge$fit, c("rmsea", "rmsea.pvalue", "srmr"))
      })
      output$OmodedgeanalysisGraph <- renderPlot({
        SEMgraph::gplot(modedge$graph)
      })
      output$download_EdgeAnalysis_table <- downloadHandler(
        filename = function() {
          "results_for_edge_analysis.xlsx"
        },
        content = function(file) {
          write.xlsx(as.data.frame(parameterestimates(modedge$fit)), file)
        }
      )
    },
    warning = function(w){
      # output of edge analysis
      showNotification(paste0("Error occurs! Select another model", type = "error"))
  
    })
  })
}
