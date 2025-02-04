library(shiny)

fluidPage(
  titlePanel("shinyDegSEM"),
  
  sidebarLayout(
    sidebarPanel(
      h3("Step1: Data Input"),
      selectInput("species", "Choose species", choices = "hsapiens", selected = "hsapiens"),
      fileInput("data_file", "Upload the data file (in .txt format, row represents gene and column represents case)", accept = c(".txt")),
      radioButtons('fileOrText', "Sequence data: Upload file or number sequence?", choices = c("file" = "file", "text" = "text")),
      conditionalPanel(
        condition = "input.fileOrText == 'text'",
        textInput("case_sequence", "Input the sequence (only inputing number 1 and 2, using space to seperate each number)"),
        actionButton("upload_button", "Upload the sequence")
      ),
      conditionalPanel(
        condition = "input.fileOrText == 'file'",
        fileInput("sequence_file", "Upload the sequence file (in .txt format, with only number 1 and 2, using space to seperate each number)", accept = c(".txt"))
      ),
      
      h3("Step2: DEG Analysis"),
      radioButtons('SAMOrTable', "SAM analysis or Uploading the DEGs Table?", choices = c("SAM" = "SAM", "table" = "table")),
      conditionalPanel(
        condition = "input.SAMOrTable == 'SAM'",
        numericInput("delta_value", "delta(in SAM analysis):", value = 1, min = 0, max = 10, step = 0.1),
      ),
      conditionalPanel(
        condition = "input.SAMOrTable == 'table'",
        fileInput("DEGs_table", "Upload the DEGs Table")
      ),
      downloadButton("download_DEGs_table", "DEGs Result"),
      
      h3("Step3: Enrichment Analysis"),
      radioButtons('EnrichOrTable', 
                   "Enrichment Analysis or Selecting the pathway name(s)?", 
                   choices = c("enrich" = "enrich", "tab" = "tab")),
      conditionalPanel(
        condition = "input.EnrichOrTable == 'tab'",
        fileInput("pathway_table", "Upload the table of pathway names")
      ),
      conditionalPanel(
        condition = "input.EnrichOrTable == 'enrich'",
        downloadButton("download_rgraph", "Graph-type Result")
      ),
      
      h3("Step4: Network Analysis"),
      selectInput("pathway_select", 
                  "Select the pathway(s) for network visualization", 
                  choices = NULL, selected = NULL, multiple = TRUE),
      actionButton("pathway_upd_indicate", "Run Network Analysis"),
      downloadButton("download_SEM_data", "Download Pathway Data"),
      downloadButton("download_pathway_igraph", "Download Pathway Igraph"),
      
      h3("Step5: SEM Analysis"),
      selectInput("selected_pathway", 
                  "Select a pathway of interest \n(Note: click the 'Run Initial SEM' button to clean up your previous models.):",
                  choices = NULL, selected = NULL, multiple = F),
      # 其他 SEM 分析选项
      # selectInput("Iestimator", "Select SEM estimator:", 
      #             choices = c("ML", "MLM", "WLS", "WLSM", "ULS", "GLS", "FIML"),
      #             selected = "ML"),
      # 按钮：初始化模型
      actionButton("Iinitialrun", "Run Initial SEM"),
      # 按钮：运行模型修正
      # actionButton("IrunModif", "Run Modification Indices"),
      # 下拉菜单：选择要添加的 MI
      uiOutput("OModifSel"),
      # 按钮：重新运行模型
      actionButton("Irerun", "Add the path and run the model again"),
      # 下拉菜单：选择最终模型
      uiOutput("ImodifiedModSel"),
      ## 按钮：选择最终模型的图
      #actionButton("finalPlotShow", "Show the Final Model"),
      hr(),
      # 按钮：运行测量不变性检验
      actionButton("IRunMI", "Run Measurement Invariance"),
      # 按钮：运行节点分析
      actionButton("IRunNodeAnalysis", "Run Node Analysis"),
      # 按钮：运行边分析
      actionButton("IRunEdgeAnalysis", "Run Edge Analysis"),
      hr(),
      downloadButton("download_FinalSEM_table", "Final SEM Result"),
      downloadButton("download_NodeAnalysis_table", "Node Analysis Result"),
      downloadButton("download_EdgeAnalysis_table", "Edge Analysis Result")
    ),
    
    mainPanel(
      tabsetPanel(
        id = "panel1",
        tabPanel(
          "Data Input",
          strong("Preview of the uploaded data is shown below."),
          # textOutput("input_data"),
          p(" "),
          textOutput("input_data1"), 
          textOutput("input_data2"), 
          textOutput("input_data3"), 
          textOutput("input_data4"),
          p(" "),
          h3("Data preview"),
          textOutput("input_data5"),
          textOutput("input_data6"),
          tableOutput("previewdata")
        ),
        
        tabPanel(
          "DEGs acquirement",
          strong("The analysis isn't run yet."),
          p(" "),
          # strong("For the next step, complete DEG Analysis on the left panel and click the tab Enrichment analysis - Get pathway."),
          # p(" "),
          
          conditionalPanel(
            condition = "input.SAMOrTable == 'SAM'",
            h3("SAM plot"),
            plotOutput("deg_analysis_output1"),
            h3("SAM Table"),
            tableOutput("deg_analysis_output2"),
            h3("Upregulated Genes"),
            tableOutput("deg_analysis_output3"),
            # h3("DEG names are listed here:"),
            # textOutput("deg_analysis_output4")
          ),
          
          conditionalPanel(
            condition = "input.SAMOrTable == 'table'",
            h3("DEG input"),
            textOutput("deg_input")
          )
        ),

        tabPanel(
          "Enrichment analysis - Get pathway",
          strong("The analysis isn't run yet."),
          p(" "),
          # strong("For the next step, click Network Analysis."),
          conditionalPanel(
            condition = "input.EnrichOrTable == 'enrich'",
            h3("Result of Enrichment Analysis"),
            tableOutput("pathway_ori")
          ),
          conditionalPanel(
            condition = "input.EnrichOrTable == 'tab'",
            h3("Enrichment input"),
            textOutput("enrich_input")
          )
        ),

        tabPanel(
          "Network Analysis",
          uiOutput("network_analysis")
        ),
        
        tabPanel(
          "Structural Equation Modeling",
          value = "panel1_sem",
          uiOutput("sem_analysis")
        )  
      )
    )
  )
)
