########################### DILA Pan Cancer shiny app ##########################
## This app takes as input a list of genes of interest in the first column of a
## CSV format file. The user can select Cancer Type, Major Cell Type and Cell
## Sub Type and the app will generate a bar chart displaying the frequences of 
## the cell types for that cancer type as well as a scatter plot displaying
## correlation of cell frequency against the signature score of genes of
## interest

## load libraries
library("Seurat")
library("dplyr")
library("tibble")
library("ggplot2")
library("patchwork")
library("ggpubr")
library("plotly")
library("shiny")
library("shinydashboard")

## TODO
## 1) we need consistant names, e.g. (cell_sub_type has all underscores or none)
## gsub("_", " ", gsub("_sub", "", sort(unique(meta_data$CellType_lev5))))
## SOLVED: just use the current names as they are...
## 2) find some way that all the datasets don't have to be pre-loaded
## 3) dynamic dropdown labels

############################ load datasets ##################################

  ## load toy datasets
  # rds1 = readRDS(file = paste0("../test_data/pancancer_1000a.rds"))
  # rds2 = readRDS(file = paste0("../test_data/pancancer_1000b.rds"))
  
  ## load 2 smallest datasets
  rds1 = readRDS(file = paste0("../test_data/real_datasets/CRC_fixed.rds"))
  rds2 = readRDS(file = paste0("../test_data/real_datasets/HNSCC_fixed.rds"))
  
  ## TODO:
  ## find a way to dynamically load datasets
  # c("BC_BK_fixed.rds", "BC_SYN_fixed.rds",  "CC_fixed.rds",  "CRC_fixed.rds",
  # "GBM_fixed.rds", "HCC_fixed.rds", "HNSCC_fixed.rds", "Melanoma_fixed.rds",
  # "NSCLC_A_fixed.rds", "NSCLC_fixed.rds", "OV_fixed.rds")

############################## User Interface ##################################
ui <- dashboardPage(
        dashboardHeader(title = "Pan Cancer App"),
        dashboardSidebar(
          br(),
          h5("Select Pan Cancer Features ", align = "center"),
          ## choose cancer cell type
          selectInput(inputId = "cancer_type",
                      label =  "Cancer Type:",
                      choices = c("Toy Dataset A" = "rds1",
                                  "Toy Dataset B" = "rds2")),
          # selectInput(inputId = "cancer_type",
          #             label =  "Cancer Type:",
          #             choices = c("BC" = "BC_df",
          #                         "CC" = "CC_df",
          #                         "CRC" = "CRC_df",
          #                         "Endometrial" = "Endometrial_df",
          #                         "GBM" = "GBM_df",
          #                         "HCC" = "HCC_df",
          #                         "HGSOC" = "HGSOC_df",
          #                         "HNSCC" = "HNSCC_df",
          #                         "Melanoma" = "Melanoma_df",
          #                         "NSCLC" = "NSCLC_df",
          #                         "OV" = "OV_df")),
          ## choose major cell type 
          selectInput(inputId = "major_cell_type",
                      label =  "Major Cell Type:",
                      choices = c("B cell" = "B_cell_major",
                                  "Cancer" = "Cancer_major",
                                  "cDC" = "cDC_major",
                                  "DC" = "DC_major",
                                  "Endothelial" = "Endothelial_major",
                                  "Fibroblast" = "Fibroblast_major",
                                  "Hepatocyte" = "Hepatocyte_major",
                                  "Mast Cell" = "Mast_Cell_major",
                                  "Melanocytes" = "Melanocyte_major",
                                  "Myeloid" = "Myeloid_major",
                                  "pDC" = "pDC_major",
                                  "Tcell" = "Tcell_major")),
          ## choose cell sub type 
          selectInput(inputId = "cell_sub_type",
                      label =  "Cell Sub Type:",
                      choices = c("Arterial" = "Arterial_sub",
                                  "CD4+ effector-memory" = "CD4+_effector-memory_sub",
                                  "CD4+ follicular helper" = "CD4+_follicular_helper_sub",
                                  "CD4+ naive" = "CD4+_naive_sub",
                                  "CD4+ regulatory" = "CD4+_regulatory_sub",
                                  "CD4+ T-helper 1" = "CD4+_T-helper_1_sub",
                                  "CD4+ T-helper 17" = "CD4+_T-helper_17_sub",
                                  "CD8- gamma-delta" = "CD8-_gamma-delta_sub",
                                  "CD8+ EMRA" = "CD8+_EMRA_sub",
                                  "CD8+ effector memory" = "CD8+_effector_memory_sub",
                                  "CD8+ gamma-delta" = "CD8+_gamma-delta_sub",
                                  "CD8+ MAIT" = "CD8+_MAIT_sub",
                                  "CD8+ naive" = "CD8+_naive_sub",
                                  "CD8+ resident memory" = "CD8+_resident_memory_sub",
                                  "cDC2" = "cDC2_sub",
                                  "Fibroblast" = "Fibroblast_sub",
                                  "Hepatocyte" = "Hepatocyte_sub",
                                  "igA mature" = "igA_mature_sub",
                                  "igG mature" = "igG_mature_sub",
                                  "Lymphatic" = "Lymphatic_sub",
                                  "Macro CCL18" = "Macro_CCL18_sub",
                                  "Macro CCL2" = "Macro_CCL2_sub",
                                  "Macro CCL2 CCL18" = "Macro_CCL2_CCL18_sub",
                                  "Macro CCR2" = "Macro_CCR2_sub",
                                  "Macro CCR2 CX3CR1" = "Macro_CCR2_CX3CR1_sub",
                                  "Macro CX3CR1" = "Macro_CX3CR1_sub",
                                  "Macro CXCL10" = "macro_CXCL10_sub",
                                  "Macro LYVE1" = "Macro_LYVE1_sub",
                                  "Macro MT1G" = "Macro_MT1G_sub",
                                  "Macro prolif" = "Macro_prolif_sub",
                                  "Macro SLC2A1" = "Macro_SLC2A1_sub",
                                  "MAIT" = "MAIT_sub",
                                  "Mast Cell" = "Mast_Cell_sub",
                                  "Melanocytes" = "Melanocytes_sub",
                                  "Memory igM-" = "Memory_igM-_sub",
                                  "Memory igM+" = "Memory_igM+_sub",
                                  "Mono CD14" = "Mono_CD14_sub",
                                  "Mono CD16" = "Mono_CD16_sub",
                                  "Naive mature" = "Naive_mature_sub",
                                  "Neutrophils" = "Neutrophils_sub",
                                  "NK cytotoxic" = "NK_cytotoxic_sub",
                                  "NK inflammatory" = "NK_inflammatory_sub",
                                  "PCV" = "PCV_sub",
                                  "pDC" = "pDC_sub",
                                  "Proliferative endothelial" = "Proliferative_endothelial_sub",
                                  "Proliferative T-cell" = "Proliferative_T-cell_sub",
                                  "Quiesc_mig_DC" = "Quiesc_mig_DC_sub",
                                  "Reg Memory" = "Reg_Memory_sub",
                                  "Reg_mig_DC" = "Reg_mig_DC_sub",
                                  "Stalk cells" = "Stalk_cells_sub",
                                  "Tip cells" = "Tip_cells_sub",
                                  "Venous" = "Venous")),
          ## pasting list of genes into text box may be added later
          # textInput(inputId = "text_box",
          #           label = "Paste Genes of Interest Here",
          #           width = "400px",
          #           placeholder = "CD8"
          #           ),
          ## upload genes of interest
          tags$hr(),  ## makes line
          fileInput(inputId = "file1",
                    label = "Upload CSV file containing genes of interest",
                    multiple = FALSE,
                    accept = c("text/csv",
                               "text/comma-separated-values,text/plain",
                               ".csv")),
          # tags$hr(),  ## makes line
          h5("CSV file parameters", align = "center"),
          radioButtons(inputId = "sep",
                       label = "Separator",
                       choices = c(Comma = ",",
                                   Semicolon = ";",
                                   Tab = "\t"),
                       selected = ","),
          radioButtons(inputId = "quote",
                       label = "Quote",
                       choices = c(None = "",
                                   "Double Quote" = '"',
                                   "Single Quote" = "'"),
                       selected = '"'),
          checkboxInput(inputId = "header",
                        label = "Header",
                        value = TRUE)),
        ## end dashboardSidebar
        dashboardBody(
          tabsetPanel(type = "tabs",
              tabPanel(title = "Genes of Interest",
                     fluidRow(column(width = 12,
                     # textOutput("cell_sub_type"),
                     h5("This displays the list of uploaded genes of interest"),
                     DT::dataTableOutput("table1"), align = "left"))),
              tabPanel(title = "Metadata",
                     fluidRow(column(width = 12,
                     h5("This displays the metadata of dataset of interest"),
                     DT::dataTableOutput("table2"), align = "left"))),
              tabPanel(title = "Bar Charts",
                       fluidRow(column(width = 12,
                                       h5(paste0("This plot displays the Mean ",
                                       "Expression Score for the uploaded ",
                                       "genes of interest for major cell types ",
                                       "from the selected cancer dataset"),
                                       align = "left"),
                                       plotlyOutput("bar1"),
                                       align = "left")),
                       fluidRow(column(width = 12,
                                       h5(paste0("This plot displays the Mean ",
                                       "Expression Score for the uploaded ",
                                       "genes of interest for the minor cell ",
                                       "sub types from the selected cancer ",
                                       "dataset corresponding to the selected ",
                                       "major cell type"), align = "left"),
                                       plotlyOutput("bar2"),
                                       align = "left"))),
              tabPanel(title = "Scatter Plots",
                       fluidRow(column(width = 12,
                                       h5(paste0("This plot displays a scatter ",
                                       "plot illustrating the signature score ",
                                       "for the uploaded genes of interest ",
                                       "plotted against proportion of the ",
                                       "selected major cell type across the ",
                                       "selected cancer datasets"), align = "left"),
                                       plotlyOutput("scatter1"),
                                       align = "left")),
                       fluidRow(column(width = 12,
                                       h5(paste0("This plot displays a scatter ",
                                       "plot illustrating the signature score ",
                                       "for the uploaded genes of interest ",
                                       "plotted against proportion of the ",
                                       "selected minor cell sub type across the ",
                                       "selected cancer datasets"), align = "left"),
                                       plotlyOutput("scatter2"),
                                       align = "left"))
                       )
          ) ## end tabsetPanel
        ) ## end dashboardBody 
      ) ## end dashboardPage
## end UI
################################## Server ######################################

server <- function(input, output) {
  
  ## switch reactive variable that pipes the selected dataset to sel_df()
  sel_df <- reactive({switch(input$cancer_type,
                             "rds1" = rds1,
                             "rds2" = rds2)})
  
  ############################## READ  Gene of Interest table ##################
   output$table1 <- DT::renderDataTable({
       req(input$file1)
       tryCatch({
          read.csv(input$file1$datapath,
                   header = input$header,
                   sep = input$sep,
                   quote = input$quote) %>%
          dplyr::rename(ID = names(.)[1]) -> df1
                },
      error = function(e) {## return a safeError if a parsing error occurs
        stop(safeError(e))
      })
      return(df1)
  })
  
  ############################ READ CANCER DATASETS ############################
   
  ## define CalculateSignatureScore function
  CalculateSignatureScore = function(object,
                                     features,
                                     assay = "RNA",
                                     slot = "data") {
  
  # Get the data from the Seurat object
  data = GetAssayData(object,
                      slot = slot,
                      assay = assay)[features, ]
  
  # If using log-normalized data then perform average in non-log space
  if (slot == "data")
    data = expm1(data)
  
  # Calculate average per cell
  data = data %>% Matrix::colMeans(data)
  
  # Return named vector with the signature score
  return(data)
  }
  
  output$table2 <- DT::renderDataTable({
       req(input$file1)
       tryCatch({
         ## load selected seurat metadata  
          sel_df()@meta.data %>%
          tibble::rownames_to_column() %>%
          dplyr::rename(ID = names(.)[1]) -> meta_data
         
         ## load genes of interest
         read.csv(input$file1$datapath,
                   header = input$header,
                   sep = input$sep,
                   quote = input$quote) %>%
          dplyr::rename(ID = names(.)[1]) -> df1
          
            ## calculate scores
        as.data.frame(
          CalculateSignatureScore(sel_df(), as.list(df1)[[1]])) %>%
          tibble::rownames_to_column() %>%
          dplyr::rename(ID = names(.)[1],
                        sig_score = names(.)[2]) %>%
          dplyr::left_join(meta_data, by = "ID") %>%
          dplyr::rename(cell_id = ID,
                        sample_id = `orig.ident`) %>%
          dplyr::mutate(sample_id = as.factor(sample_id),
                        CellType_lev5 = as.factor(CellType_lev5),
                        sig_score = round(sig_score, 3)) %>% 
          dplyr::select(cell_id, sig_score, sample_id, Patient, TumorType,
                        CellType, BiopsySite, `CellType_lev2.5`) -> df2
                },
      error = function(e) {## return a safeError if a parsing error occurs
        stop(safeError(e))
      })
      return(df2)
  })
  
  ################################ BAR CHART 1 #################################
  
   output$bar1 <- renderPlotly({
      req(input$file1)
      tryCatch({
         ## load selected seurat metadata
          sel_df()@meta.data %>%
          tibble::rownames_to_column() %>% 
          dplyr::rename(ID = names(.)[1]) -> meta_data
         
         ## load genes of interest
         read.csv(input$file1$datapath,
                   header = input$header,
                   sep = input$sep,
                   quote = input$quote) %>%
          dplyr::rename(ID = names(.)[1]) -> df1
          
            ## calculate scores and plot scores by cell type
        as.data.frame(
          CalculateSignatureScore(sel_df(), as.list(df1)[[1]])) %>%
          tibble::rownames_to_column() %>%
          dplyr::rename(ID = names(.)[1],
                        sig_score = names(.)[2]) %>%
          dplyr::left_join(meta_data, by = "ID") %>%
          dplyr::rename(cell_id = ID,
                        sample_id = `orig.ident`) %>%
          dplyr::mutate(sample_id = as.factor(sample_id),
                        CellType_lev5 = as.factor(CellType_lev5)) %>%
          dplyr::group_by(CellType) %>%
          dplyr::summarise(mean_sig_score = mean(sig_score)) %>%
          dplyr::arrange(desc(mean_sig_score)) %>%
          dplyr::mutate(CellType = factor(CellType, levels = unique(CellType))) %>%
          ggplot(aes(x = mean_sig_score, y = CellType, fill = CellType)) +
          geom_bar(stat = "identity") +
          coord_flip(expand = TRUE) +
          theme_classic() +
          ylab("Cell Types") +
          xlab("Mean Expression Score") +
          theme(legend.position = "none",
                axis.text.x = element_text(angle = 45, hjust = 1)) +
          scale_fill_viridis_d() +
          ggtitle("Mean Expression Score by All Major Cell Types") -> p1
               },
      error = function(e) {## return a safeError if a parsing error occurs
        stop(safeError(e))
      })
      return(p1)
   })

  ################################ BAR CHART 2 #################################
  
   output$bar2 <- renderPlotly({
      req(input$file1)
      tryCatch({
         ## load selected seurat metadata
          sel_df()@meta.data %>%
          tibble::rownames_to_column() %>%
          dplyr::rename(ID = names(.)[1]) -> meta_data

         ## load genes of interest
         read.csv(input$file1$datapath,
                   header = input$header,
                   sep = input$sep,
                   quote = input$quote) %>%
          dplyr::rename(ID = names(.)[1]) -> df1

            ## calculate scores and plot scores by cell type
        as.data.frame(
          CalculateSignatureScore(sel_df(), as.list(df1)[[1]])) %>%
          tibble::rownames_to_column() %>%
          dplyr::rename(ID = names(.)[1],
                        sig_score = names(.)[2]) %>%
          dplyr::left_join(meta_data, by = "ID") %>%
          dplyr::rename(cell_id = ID,
                        sample_id = `orig.ident`) %>%
          dplyr::filter(CellType == gsub(
            "_", " ", gsub("_major", "", input$major_cell_type))) %>%
          dplyr::mutate(sample_id = as.factor(sample_id),
                        CellType_lev6 = as.factor(CellType_lev6)) %>%
          dplyr::group_by(CellType_lev6) %>%
          dplyr::summarise(mean_sig_score = mean(sig_score)) %>%
          dplyr::arrange(desc(mean_sig_score)) %>%
          dplyr::mutate(CellType_lev6 = factor(
            CellType_lev6, levels = unique(CellType_lev6))) %>%
          ggplot(aes(x = mean_sig_score, y = CellType_lev6,
                     fill = CellType_lev6)) +
          geom_bar(stat = "identity") +
          coord_flip(expand = TRUE) +
          theme_classic() +
          ylab("Cell Sub Types") +
          xlab("Mean Expression Score") +
          theme(legend.position = "none",
                axis.text.x = element_text(angle = 45, hjust = 1)) +
          scale_fill_viridis_d() +
          ggtitle(paste0("Mean Expression Score by Selected ",
                         gsub("_", " ", gsub("_major", "", input$major_cell_type
                                             )), " Cell Sub Types")) -> p1
               },
      error = function(e) {## return a safeError if a parsing error occurs
        stop(safeError(e))
      })
      return(p1)
   })

  ################################ SCATTER PLOT 1 ##############################
  
  # input_1 <- reactive(get(input$cell_sub_type))
  
  output$scatter1 <- renderPlotly({
      req(input$file1)
      tryCatch({
         ## load selected seurat metadata
          sel_df()@meta.data %>%
          tibble::rownames_to_column() %>% 
          dplyr::rename(ID = names(.)[1]) -> meta_data
         
         ## load genes of interest
         read.csv(input$file1$datapath,
                   header = input$header,
                   sep = input$sep,
                   quote = input$quote) %>%
          dplyr::rename(ID = names(.)[1]) -> df1
         
          ## calculate scores and plot scores by cell type
        as.data.frame(
          CalculateSignatureScore(sel_df(), as.list(df1)[[1]])) %>%
          tibble::rownames_to_column() %>%
          dplyr::rename(ID = names(.)[1],
                        sig_score = names(.)[2]) %>%
          dplyr::left_join(meta_data, by = "ID") %>%
          dplyr::rename(cell_id = ID,
                        sample_id = `orig.ident`) %>%
          dplyr::mutate(sample_id = as.factor(sample_id),
                        CellType = as.factor(CellType)) -> df2
          
          # Calculate the frequencies
          df2 %>%
            dplyr::group_by(sample_id, TumorType3, CellType, .drop = FALSE) %>%
            dplyr::summarise(n = n()) %>%
            dplyr::mutate(freq = n/sum(n)) %>%
            dplyr::filter(CellType == paste0(
              gsub("_", " ", gsub("_major", "", input$major_cell_type)))
              ) -> res1
          
          # Calculate the average score per sample and merge with other dataframe
          df2 %>%
            dplyr::group_by(sample_id) %>%
            dplyr::summarise_at(vars(sig_score), mean) %>%
            merge(res1, by = "sample_id") %>%
            dplyr::filter(CellType == paste0(
              gsub("_", " ", gsub("_major", "", input$major_cell_type)))) %>%
            ggplot(aes(x = freq, y = sig_score)) +
            geom_point(size = 1, aes(color = TumorType3)) +
            stat_smooth(method = "lm", se = TRUE, fill = "gray",
                        color = "black", formula = y ~ poly(x, 1, raw = TRUE)) +
            ggpubr::stat_cor(method = "spearman") +
            theme_classic() +
            theme(axis.title = element_text(size = 20)) +
            labs(y = "Signature score",
                 x = paste0(gsub("_", " ",
                                 gsub("_major", "",
                                      input$major_cell_type)), " proportion")
                 ) -> p2
               },
      error = function(e) {## return a safeError if a parsing error occurs
        stop(safeError(e))
      })
      return(p2)
   })
  
  
  ################################ SCATTER PLOT 2 ##############################
  
  # input_1 <- reactive(get(input$cell_sub_type))
  
  output$scatter2 <- renderPlotly({
      req(input$file1)
      tryCatch({
         ## load selected seurat metadata
          sel_df()@meta.data %>%
          tibble::rownames_to_column() %>% 
          dplyr::rename(ID = names(.)[1]) -> meta_data
         
         ## load genes of interest
         read.csv(input$file1$datapath,
                   header = input$header,
                   sep = input$sep,
                   quote = input$quote) %>%
          dplyr::rename(ID = names(.)[1]) -> df1
         
          ## calculate scores and plot scores by cell type
        as.data.frame(
          CalculateSignatureScore(sel_df(), as.list(df1)[[1]])) %>%
          tibble::rownames_to_column() %>%
          dplyr::rename(ID = names(.)[1],
                        sig_score = names(.)[2]) %>%
          dplyr::left_join(meta_data, by = "ID") %>%
          dplyr::rename(cell_id = ID,
                        sample_id = `orig.ident`) %>%
          dplyr::mutate(sample_id = as.factor(sample_id),
                        CellType_lev5 = as.factor(CellType_lev5)) -> df2
          
          # Calculate the frequencies
          df2 %>%
            dplyr::group_by(sample_id, TumorType3, CellType_lev5,
                            .drop = FALSE) %>%
            dplyr::summarise(n = n()) %>%
            dplyr::mutate(freq = n/sum(n)) %>%
            dplyr::filter(CellType_lev5 == paste0(
              gsub("_", " ", gsub("_sub", "", input$cell_sub_type)))) -> res1
          
          # Calculate the average score per sample and merge with other dataframe
          df2 %>%
            dplyr::group_by(sample_id) %>%
            dplyr::summarise_at(vars(sig_score), mean) %>%
            merge(res1, by = "sample_id") %>%
            dplyr::filter(CellType_lev5 == paste0(
              gsub("_", " ", gsub("_sub", "", input$cell_sub_type)))) %>%
            ggplot(aes(x = freq, y = sig_score)) +
            geom_point(size = 1, aes(color = TumorType3)) +
            stat_smooth(method = "lm", se = TRUE, fill = "gray",
                        color = "black", formula = y ~ poly(x, 1, raw = TRUE)) +
            ggpubr::stat_cor(method = "spearman") +
            theme_classic() +
            theme(axis.title = element_text(size = 20)) +
            labs(y = "Signature score",
                 x = paste0(gsub("_", " ",
                                 gsub("_sub", "",
                                      input$cell_sub_type)), " proportion")
                 ) -> p2
               },
      error = function(e) {## return a safeError if a parsing error occurs
        stop(safeError(e))
      })
      return(p2)
   })
  
## end server   
}
## Run the application 
shinyApp(ui = ui, server = server)
