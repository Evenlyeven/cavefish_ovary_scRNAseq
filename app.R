#load library
library(shiny)
library(shinythemes)
library(shinyWidgets)
library(scales)
library(ggplot2)
library(Seurat)
library(scCustomize)
library(patchwork)
library(cowplot)
library(dplyr)
library(magrittr)
library(shinyjs)
library(DT)

## ============================== configuration ============================== ##
##seurat obj
seurat_obj <-
  readRDS(file = "./data/seuobj.rds")

##color for UMAP
colors_to_use <- readRDS("./data/colors.rds")
colors_umap <-  c("#979FA2","#93D2DF","#CDB4DB","#94D2BD","#E9D8A6","#EE9B00","#CA6702","#D8AE9B","#DF76BC","#E2696D")
purple_gradient <- rev(c("#231942", "#5E548E", "#9F86C0", "#BE95C4", "#E0B1CB", "#CED4DA"))
red_gradient <- rev(c("#461220", "#8C2F39", "#B23A48", "#FCB9B2", "#FED0BB", "#CED4DA"))
colors_fish_group <- c("#979FA2", "#EE9B00", "#93D2DF")
colors_sample <- c("#486090", "#D6BBCF", "#A8A890", "#7890A8", "#F0D8C0", "#A8C0C0")

##sample genes
smpl_genes_sm <- paste0("ddx4 gsdf")
smpl_genes_bg <- paste0("ddx4 gsdf cyp11a1.1 mpx col1a1a")

##images list
#UMAP page: img1.png, img2.png, img3.png
img1 <- "img1.png"#seurat cluster
img2 <- "img2.png"#sample_type
img3 <- "img3.png"#sample_name
img4 <- "img4.png"#annotation

##add ensembl ID to DE tab
res_j <- readRDS("./data/feat_tb_tidy.rds")

##data structure
###featureplot tab
####seurat_obj$orig.ident: each sample identity is stored at seurat_obj$orig.ident

##helper function
caseinsenFun <- function(input, seurat_obj){
  selected <- unique(unlist(strsplit(input, " ")))
  for (j in seq_along(selected)){
    selected[j] <- ifelse(selected[j] %in% rownames(seurat_obj), selected[j], 
                          grep(pattern = paste0("^", selected[j], "$"), x = rownames(seurat_obj), ignore.case = TRUE, value = TRUE))
  }
  selected <- na.omit(selected)
  return(selected)
}

## =========================================================================== ##



## ================================ define UI ================================ ##
# Define UI -----
ui <- fluidPage(
  theme = shinytheme("spacelab"),
  
  titlePanel("Ovary scRNA-seq Data"),
  
  tags$style(HTML("
                  .fluid-image {
                  display: block;
                  margin: 0 auto;
                  max-width: 95%;
                  height: auto;
                  }")),
  
  tabsetPanel(
    #=====================
    #UMAP tab
    tabPanel(
      "UMAP",
      
      fluidRow(tags$br()),
      
      fluidRow(tags$br()),
      
      column(width = 10,
             align = "center",
             offset = 1,
             column(width = 12,
                    align = "left",
                    tags$p("This app features scRNA-seq data from the ovaries of surface fish, Pachón and Molino populations of ", 
                           tags$i("Astyanax mexicanus"), "."))),
      
      fluidRow(tags$br()),
      
      fluidRow(tags$br()),
      
      tags$hr(),
      
      fluidRow(tags$br()),
      
      fluidRow(tags$br()),
      
      fluidRow(tags$br()),
      
      fluidRow(column(
        width = 6,
        offset = 3,
        align = "center",
        tags$img(
          src = img4,
          class = "fluid-image"
        )
      )),
      
      fluidRow(tags$br()),
      
      fluidRow(tags$br()),
      
      fluidRow(tags$br()),
      
      fluidRow(column(
        width = 6,
        offset = 3,
        align = "center",
        tags$img(
          src = img1,
          class = "fluid-image"
        )
      )),
      
      fluidRow(tags$br()),
      
      fluidRow(tags$br()),
      
      fluidRow(tags$br()),
      
      fluidRow(column(
        width = 6,
        offset = 3,
        align = "center",
        tags$img(
          src = img2,
          class = "fluid-image"
        )
      )),
      
      fluidRow(tags$br()),
      
      fluidRow(tags$br()),
      
      fluidRow(tags$br()),
      
      fluidRow(column(
        width = 8,
        offset = 2,
        align = "center",
        tags$img(
          src = img3,
          class = "fluid-image"
        )
      )),
      
      tags$hr(),
      
      fluidRow(tags$br()),
      
      useShinyjs(),  # Initialize shinyjs
      fluidRow(
        column(
          width = 12,
          align = "center",
          tags$h3("Select Cell Subclusters")
        ),
        
        column(
          width = 12,
          align = "center",
          selectInput(
            "Cell_sub",
            label = "",
            choices = c(
              "Follicular somatic cells"
            ),
            selected = "Follicular somatic cells"
          )
        ),
        
        column(
          width = 12,
          align = "center",
          actionButton("lauchApp", "Go!", icon("arrow-right"))
        )
      ),
      
      fluidRow(tags$br()),
      
      fluidRow(tags$br()),
      
      fluidRow(tags$br()),
      
      # fluidRow(
      #   column(
      #     width = 12,
      #     align = "center",
      #     tags$h3("Select Cell Subclusters")
      #   ),
      #   
      #   column(
      #     width = 12,
      #     align = "center",
      #     selectInput(
      #       "Cell_sub",
      #       label = "",
      #       choices = c(
      #         #"Germ cells",
      #         "Follicular somatic cells"#,
      #         #"Immune cells"
      #       ),
      #       selected = "Follicular somatic cells"
      #     )
      #   ),
      #   
      #   column(
      #     width = 12,
      #     align = "center",
      #     actionButton("lauchApp", "Go!", icon("arrow-right")),
      #     uiOutput("appLink")
      #     
      #   ),
      #   
      #   fluidRow(tags$br()),
      #   
      #   fluidRow(tags$br()),
      #   
      #   fluidRow(tags$br()),
      #   
      #   fluidRow(tags$br())),
      
      column(width = 10,
             align = "center", 
             offset = 1,
             column(width = 12,
                    align = "center",
                    tags$p(
                      'If you encounter issues using this App, please contact dwu@stowers.org'
                    ))),
      
      fluidRow(tags$br()),
      
      fluidRow(tags$br())
    ),
    
    #=====================
    tabPanel(
      "Feature Plots",
      fluid = FALSE,
      
      sidebarLayout(
        fluid = TRUE,
        
        sidebarPanel(
          width = 4,
          
          column(
            width = 12,
            align = "left",
            textInput("featureGenes", "Insert gene name(s):",
                      value = smpl_genes_sm)
          ),
          
          column(
            width = 12,
            align = "center",
            actionButton("runFeatPlot", "Generate Plot(s)",
                         style = "padding:5px; font-size:80%")
          ),
          
          column(width = 12, tags$hr(width = "50%"), align = "center"),
          
          column(
            width = 12,
            align = "center",
            downloadButton("downloadFeaturePlotF", "Download pdf",
                           style = 'padding:5px; font-size:80%')
          ),
          
          column(width = 12, tags$br()),
          
          column(
            width = 12,
            align = "center",
            uiOutput("cellSelectFeat")
          ),
          
          column(width = 12,
                 tags$br()),
          
          column(
            width = 12,
            align = "center",
            uiOutput("cellSelectCol")
          ),
          
          column(width = 12,
                 tags$br()),
          
          column(
            width = 12,
            align = "center",
            numericInput(
              "ptSizeFeature",
              "Input Cell Size:",
              value = 0.50,
              min = 0.25,
              step = 0.25,
              max = 2.00,
              width = "40%"
              
            )
          ),
          
          fluidRow(tags$br()),
          
          fluidRow(tags$br()),
          
          fluidRow(column(
            width = 12,
            align = "center",
            tags$img(
              src = img4,
              class = "fluid-image"
            )
          )),
          
          fluidRow(tags$br()),
          
          fluidRow(column(
            width = 12,
            align = "center",
            tags$img(
              src = img1,
              class = "fluid-image"
            )
          )),
          
          fluidRow(tags$br()),
          
          fluidRow(tags$br())
        ),
        
        mainPanel(
          fluidRow(
            column(width = 8, tags$br()),
            column(
              width = 8,
              tags$b("Mismatches or genes not present"),
              "(if applicable)",
              tags$b(":")
            ),
            column(width = 8, uiOutput("notInFeat")),
            column(width = 8, tags$hr()),
            
            fluidRow(tags$br()),
            column(width = 12, uiOutput("plot.uiFeaturePlotF"))
          )
        )
      )
    ),
    
    #=====================
    tabPanel(
      "Feature Plots by Sample",
      fluid = FALSE,
      
      sidebarLayout(
        fluid = TRUE,
        
        sidebarPanel(
          width = 4,
          
          column(
            width = 12,
            align = "center",
            textInput(
              "featureGenes_split",
              "Insert ONE gene name:",
              value = "ddx4",
              width = "50%"
            )
          ),
          
          column(
            width = 12,
            align = "center",
            actionButton("runFeatPlot_split", "Generate Plot(s)",
                         style = "padding:5px; font-size:80%")
          ),
          
          column(width = 12, tags$hr(width = "50%"), align = "center"),
          
          column(
            width = 12,
            align = "center",
            downloadButton("downloadFeaturePlotF_split", "Download pdf",
                           style = 'padding:5px; font-size:80%')
          ),
          
          column(width = 12, tags$br()),
          
          column(
            width = 12,
            align = "center",
            uiOutput("cellSelectFeat_split")
          ),
          
          column(width = 12,
                 tags$br()),
          
          column(
            width = 12,
            align = "center",
            uiOutput("cellSelectCol_tab3")
          ),
          
          # column(
          #   width = 12,
          #   column(
          #     width = 6,
          #     textInput(
          #       "CellBackCol_split",
          #       "Cell background color:",
          #       value = "lightgrey"
          #     )
          #   ),
          #   column(
          #     width = 6,
          #     textInput("CellForeCol_split",
          #               "Cell foreground color:", value = "tomato")
          #   )
          # ),
          # 
          column(width = 12,
                 tags$br()),
          
          column(
            width = 12,
            align = "center",
            numericInput(
              "ptSizeFeature_split",
              "Input Cell Size:",
              value = 0.50,
              min = 0.25,
              step = 0.25,
              max = 2.00,
              width = "40%"
            )
          ),
          
          fluidRow(tags$br()),
          
          fluidRow(tags$br()),
          
          fluidRow(column(
            width = 12,
            align = "center",
            tags$img(
              src = img4,
              class = "fluid-image"
            )
          )),
          
          fluidRow(tags$br()),
          
          fluidRow(column(
            width = 12,
            align = "center",
            tags$img(
              src = img1,
              class = "fluid-image"
            )
          )),
          
          fluidRow(tags$br()),
          
          fluidRow(tags$br())
        ),
        
        mainPanel(
          fluidRow(
            column(width = 8, tags$br()),
            column(
              width = 8,
              tags$b("Mismatches or genes not present"),
              "(if applicable)",
              tags$b(":")
            ),
            column(width = 8, uiOutput("notInFeat_split")),
            column(width = 8, tags$hr()),
            
            fluidRow(tags$br()),
            column(width = 12, uiOutput("plot.uiFeaturePlotF_split"))
          )
        )
      )
    ),
    
    #=====================
    tabPanel(
      "Feature Plots by Fish Population",
      fluid = FALSE,
      
      sidebarLayout(
        fluid = TRUE,
        
        sidebarPanel(
          width = 4,
          
          column(
            width = 12,
            align = "center",
            textInput(
              "featureGenes_split2",
              "Insert ONE gene name:",
              value = "ddx4",
              width = "50%"
            )
          ),
          
          column(
            width = 12,
            align = "center",
            actionButton("runFeatPlot_split2", "Generate Plot(s)",
                         style = "padding:5px; font-size:80%")
          ),
          
          column(width = 12, tags$hr(width = "50%"), align = "center"),
          
          column(
            width = 12,
            align = "center",
            downloadButton("downloadFeaturePlotF_split2", "Download pdf",
                           style = 'padding:5px; font-size:80%')
          ),
          
          column(width = 12, tags$br()),
          
          column(
            width = 12,
            align = "center",
            uiOutput("cellSelectCol_tab4")
          ),
          
          column(width = 12,
                 tags$br()),
          
          # column(
          #   width = 12,
          #   align = "center",
          #   uiOutput("cellSelectFeat_split2")
          # ),
          # 
          # column(width = 12,
          #        tags$br()),
          
          # column(
          #   width = 12,
          #   column(
          #     width = 6,
          #     textInput(
          #       "CellBackCol_split2",
          #       "Cell background color:",
          #       value = "lightgrey"
          #     )
          #   ),
          #   column(
          #     width = 6,
          #     textInput("CellForeCol_split2",
          #               "Cell foreground color:", value = "tomato")
          #   )
          # ),
          
          column(
            width = 12,
            align = "center",
            numericInput(
              "ptSizeFeature_split2",
              "Input Cell Size:",
              value = 0.50,
              min = 0.25,
              step = 0.25,
              max = 2.00,
              width = "40%"
            )
          ),
          
          column(width = 12,
                 tags$br()),
          
          fluidRow(tags$br()),
          
          fluidRow(tags$br()),
          
          fluidRow(column(
            width = 12,
            align = "center",
            tags$img(
              src = img4,
              class = "fluid-image"
            )
          )),
          
          fluidRow(tags$br()),
          
          fluidRow(column(
            width = 12,
            align = "center",
            tags$img(
              src = img1,
              class = "fluid-image"
            )
          )),
          
          fluidRow(tags$br()),
          
          fluidRow(tags$br())
        ),
        
        mainPanel(
          fluidRow(
            column(width = 8, tags$br()),
            column(
              width = 8,
              tags$b("Mismatches or genes not present"),
              "(if applicable)",
              tags$b(":")
            ),
            column(width = 8, uiOutput("notInFeat_split2")),
            column(width = 8, tags$hr()),
            
            fluidRow(tags$br()),
            column(width = 12, uiOutput("plot.uiFeaturePlotF_split2"))
          )
        )
      )
    ),
    
    #=====================
    tabPanel(
      "Violin Plots",
      sidebarLayout(
        fluid = TRUE,
        
        sidebarPanel(
          fluid = FALSE,
          width = 4,
          
          column(
            width = 12,
            textInput("vlnGenes", "Insert gene name(s):",
                      value = smpl_genes_sm)
          ),
          
          column(
            width = 12,
            align = "center",
            actionButton("runVlnPlot", "Generate Plots",
                         style = "padding:5px; font-size:80%")
          ),
          
          column(width = 12, tags$hr(width = "50%"), align = "center"),
          
          column(
            width = 12,
            align = "center",
            downloadButton("downloadVlnPlot", "Download pdf",
                           style = "padding:5px; font-size:80%")
          ),
          
          column(width = 12, tags$br()),
          
          column(width = 12, tags$br()),
          
          column(
            width = 12,
            align = "center",
            radioGroupButtons(
              "selectGrpVln",
              "Group cells by:",
              choices = list(
                "Cell type" = "annotation",
                "Fish population" = "fish",
                "Sample" = "orig.ident",
                "Clusters" = "seurat_clusters"
              ),
              width = "100%",
              size = "xs"
            )#,
          ),
          
          column(
            width = 12,
            align = "center",
            numericInput(
              "ptSizeVln",
              "Input dot size:",
              value = 0.25,
              min = 0.00,
              step = 0.75,
              max = 2.00,
              width = "40%"
            )
          ),
          
          fluidRow(tags$br()),
          
          fluidRow(tags$br()),
          
          fluidRow(column(
            width = 12,
            align = "center",
            tags$img(
              src = img4,
              class = "fluid-image"
            )
          )),
          
          fluidRow(tags$br()),
          
          fluidRow(column(
            width = 12,
            align = "center",
            tags$img(
              src = img1,
              class = "fluid-image"
            )
          )),
          
          fluidRow(tags$br()),
          
          fluidRow(tags$br())
        ),
        
        mainPanel(fluidRow(
          column(width = 8, tags$br()),
          
          column(width = 8, tags$b("Gene mismatches"), "(if present)", tags$b(":")),
          
          column(width = 8, uiOutput("notInVln")),
          
          column(width = 8, tags$hr()),
          
          column(width = 12, uiOutput("plot.uiVlnPlotF"))
        ))
      )
    ),
    
    #=====================
    tabPanel(
      "Dot Plot",
      sidebarLayout(
        fluid = TRUE,
        
        sidebarPanel(
          fluid = FALSE,
          width = 4,
          
          column(
            width = 12,
            align = "left",
            textInput(
              "dotGenes",
              "Insert gene name(s):",
              value = smpl_genes_bg
            ),
            checkboxInput("dPlotClust_fea",
                          label = "Check box to enable feature clustering", value = FALSE),
            checkboxInput("dPlotClust_ident",
                          label = "Check box to enable identity clustering", value = FALSE)
          ),
          
          column(
            width = 12,
            align = "center",
            actionButton("runDotPlot", "Generate Plots",
                         style = "padding:5px; font-size:80%")
          ),
          
          column(width = 12,
                 tags$hr(width = "50%"),
                 align = "center"),
          
          column(
            width = 12,
            align = "center",
            downloadButton("downloadDotPlot", "Download pdf",
                           style = "padding:5px; font-size:80%")
          ),
          
          column(width = 12, tags$br()),
          
          column(width = 12, tags$br()),
          
          column(
            width = 12,
            align = "center",
            radioGroupButtons(
              "selectGrpDot",
              "Group cells by:",
              choices = list(
                "Cell type" = "annotation",
                "Fish population" = "fish",
                "Sample" = "orig.ident",
                "Clusters" = "seurat_clusters"
              ),
              width = "100%",
              size = "xs"
            )
          ),
          
          column(
            width = 12,
            numericInput(
              "dotScale",
              HTML("<b>Dot size</b><br>(supported only for non-clustered plots):"),
              value = 10,
              min = 4,
              step = 1,
              max = 100,
              width = "80%"
            ),
            align = "center"
          ),
          
          fluidRow(tags$br()),
          fluidRow(tags$br()),
          
          fluidRow(column(
            width = 12,
            align = "center",
            tags$img(
              src = img4,
              class = "fluid-image"
            )
          )),
          
          fluidRow(tags$br()),
          
          fluidRow(column(
            width = 12,
            align = "center",
            tags$img(
              src = img1,
              class = "fluid-image"
            )
          )),
          
          fluidRow(tags$br()),
          
          fluidRow(tags$br())
        ),
        
        mainPanel(
          fluidRow(
            column(width = 8, tags$br()),
            
            column(
              width = 8,
              tags$b("Mismatches or genes not present"),
              "(if applicable)",
              tags$b(":")
            ),
            
            column(width = 8, uiOutput("notInDot")),
            
            column(width = 8, tags$hr()),
            
            column(
              width = 8,
              align = "left",
              column(
                width = 3,
                align = "left",
                numericInput(
                  "manAdjustDotW",
                  label = "Width (inches):",
                  value = 10,
                  step = 1,
                  width = "100%"
                )
              ),
              column(
                width = 3,
                align = "left",
                numericInput(
                  "manAdjustDotH",
                  label = "Height (inches):",
                  value = 6,
                  step = 1,
                  width = "100%"
                )
              )
            ),
            
            column(
              width = 8,
              align = "left",
              column(
                width = 3,
                align = "left",
                numericInput(
                  "identSize",
                  label = "Ident font size:",
                  value = 18,
                  step = 1,
                  width = "100%"
                )
              ),
              column(
                width = 3,
                align = "left",
                numericInput(
                  "FeatSize",
                  label = "Feature font size:",
                  value = 18,
                  step = 1,
                  width = "100%"
                )
              )
            ),
            
            column(
              width = 8,
              align = "left",
              column(
                width = 3,
                align = "left",
                numericInput(
                  "LegLabSize",
                  label = "Legend label font size:",
                  value = 15,
                  step = 1,
                  width = "100%"
                )
              ),
              column(
                width = 3,
                align = "left",
                numericInput(
                  "LegTitSize",
                  label = "Legend title font size:",
                  value = 15,
                  step = 1,
                  width = "100%"
                )
              )
            ),
            
            fluidRow(tags$br()),
            
            column(width = 12,
                   uiOutput("plot.uiDotPlotF")),
            
            fluidRow(
              tags$script(
                "$(document).on('shiny:connected', function(event) {
var myWidth = $(window).width();
Shiny.onInputChange('shiny_width',myWidth)

});"
              ),
              
              tags$script(
                "$(document).on('shiny:connected', function(event) {
var myHeight = $(window).height();
Shiny.onInputChange('shiny_height',myHeight)

});"
              )
            )
          )
        )
      )
    ),
    
    #=====================
    tabPanel(
      "Dot Plot by Fish",
      sidebarLayout(
        fluid = TRUE,
        
        sidebarPanel(
          fluid = FALSE,
          width = 4,
          
          column(
            width = 12,
            align = "left",
            textInput(
              "dotGenesSplit",
              "Insert gene name(s):",
              value = smpl_genes_bg
            )
          ),
          
          column(
            width = 12,
            align = "center",
            actionButton("runDotPlotSplit", "Generate Plots",
                         style = "padding:5px; font-size:80%")
          ),
          
          column(width = 12,
                 tags$hr(width = "50%"),
                 align = "center"),
          
          column(
            width = 12,
            align = "center",
            downloadButton("downloadDotPlotSplit", "Download pdf",
                           style = "padding:5px; font-size:80%")
          ),
          
          column(width = 12, tags$br()),
          
          column(width = 12, tags$br()),
          
          column(
            width = 12,
            align = "center",
            uiOutput("fishSelect")
          ),
          
          column(width = 12, tags$br()),
          
          column(width = 12, tags$br()),
          
          column(
            width = 12,
            align = "center",
            uiOutput("CluSelect")
          ),
          
          column(
            width = 12,
            numericInput(
              "dotScaleSplit",
              "Dot diameter",
              value = 10,
              min = 4,
              step = 1,
              max = 100,
              width = "40%"
            ),
            align = "center"
          ),
          
          fluidRow(tags$br()),
          fluidRow(tags$br()),
          
          fluidRow(column(
            width = 12,
            align = "center",
            tags$img(
              src = img4,
              class = "fluid-image"
            )
          )),
          
          fluidRow(tags$br()),
          
          fluidRow(column(
            width = 12,
            align = "center",
            tags$img(
              src = img1,
              class = "fluid-image"
            )
          )),
          
          fluidRow(tags$br()),
          
          fluidRow(tags$br())
        ),
        
        mainPanel(
          fluidRow(
            column(width = 8, tags$br()),
            
            column(
              width = 8,
              tags$b("Mismatches or genes not present"),
              "(if applicable)",
              tags$b(":")
            ),
            
            column(width = 8, uiOutput("notInDotSplit")),
            
            column(width = 8, tags$hr()),
            
            column(
              width = 8,
              align = "left",
              column(
                width = 3,
                align = "left",
                numericInput(
                  "manAdjustDotWSplit",
                  label = "Width (inches):",
                  value = 10,
                  step = 1,
                  width = "100%"
                )
              ),
              column(
                width = 3,
                align = "left",
                numericInput(
                  "manAdjustDotHSplit",
                  label = "Height (inches):",
                  value = 6,
                  step = 1,
                  width = "100%"
                )
              )
            ),
            
            column(
              width = 8,
              align = "left",
              column(
                width = 3,
                align = "left",
                numericInput(
                  "identSizeSplit",
                  label = "Ident font size:",
                  value = 18,
                  step = 1,
                  width = "100%"
                )
              ),
              column(
                width = 3,
                align = "left",
                numericInput(
                  "FeatSizeSplit",
                  label = "Feature font size:",
                  value = 18,
                  step = 1,
                  width = "100%"
                )
              )
            ),
            
            column(
              width = 8,
              align = "left",
              column(
                width = 3,
                align = "left",
                numericInput(
                  "LegLabSizeSplit",
                  label = "Legend label font size:",
                  value = 15,
                  step = 1,
                  width = "100%"
                )
              ),
              column(
                width = 3,
                align = "left",
                numericInput(
                  "LegTitSizeSplit",
                  label = "Legend title font size:",
                  value = 15,
                  step = 1,
                  width = "100%"
                )
              )
            ),
            
            fluidRow(tags$br()),
            
            column(width = 12,
                   uiOutput("plot.uiDotPlotSplitF")),
            
            fluidRow(
              tags$script(
                "$(document).on('shiny:connected', function(event) {
var myWidth = $(window).width();
Shiny.onInputChange('shiny_width',myWidth)

});"
              ),
              
              tags$script(
                "$(document).on('shiny:connected', function(event) {
var myHeight = $(window).height();
Shiny.onInputChange('shiny_height',myHeight)

});"
              )
            )
          )
        )
      )
    ),
    
    #=====================
    tabPanel("Cluster DE Analysis (fish group)",
             fluid = TRUE,
             sidebarLayout(
               sidebarPanel(
                 column(width = 12, 
                        align = "center",
                        uiOutput("diffOut1_less"),
                        fluidRow(tags$br()),
                        uiOutput("diffOut2_less")),
                 
                 column(width = 12, tags$br()),
                 
                 column(width = 12,
                        align = "center",
                        uiOutput("cellSelectDiff_less")),
                 
                 fluidRow(tags$br()),
                 
                 column(width = 12, 
                        align = "center", 
                        numericInput(
                          "pValCutoff_less",
                          "Set adjusted p-value cutoff:",
                          value = 0.05,
                          min = 0.00,
                          step = 0.001,
                          max = 1.00,
                          width = "210px"
                        )),
                 
                 fluidRow(tags$br()),
                 
                 column(width = 12,
                        align = "center",
                        actionButton("runDiffExp_less", "Run Differential Expression",
                                     style = "padding:5px; font-size:80%")),
                 
                 column(width = 12, 
                        align = "center",
                        tags$hr(width = "50%")),
                 
                 column(width = 12,
                        align = "center",
                        downloadButton("downloadDiffExp_less",
                                       "Download Results",
                                       style = "padding:5px; font-size:80%")),
                 
                 fluidRow(tags$br()),
                 fluidRow(tags$br()),
                 
                 fluidRow(column(
                   width = 12,
                   align = "center",
                   tags$img(
                     src = img4,
                     class = "fluid-image"
                   )
                 )),
                 
                 fluidRow(tags$br()),
                 fluidRow(tags$br())
               ),
               
               mainPanel(
                 
                 fluidRow(
                   column(width = 8, tags$br()),
                   
                   column(width = 12, 
                          align = "left",
                          class = "diffExpMain",
                          column(
                            width = 12,
                            align = "left",
                            fluidRow(tags$br()),
                            tags$b(
                              HTML("Differential expression analysis was performed using the Wilcoxon Rank test. <br> Genes are ranked in descending order based on their average log2 fold change (avg_log2FC)."))
                          ),
                          # column(
                          #   width = 12,
                          #   tags$br()
                          # ),
                          column(width = 12,
                                 tags$hr()),
                          column(
                            width = 12,
                            tags$br()
                          ),
                          uiOutput("diffTable_less")
                   )
                 )
                 
               )
             )),
    
    #=====================
    tabPanel("Cluster DE Analysis (cell type)",
             fluid = TRUE,
             sidebarLayout(
               sidebarPanel(
                 column(width = 12,
                        align = "center",
                        uiOutput("diffOut1_l_tp"),
                        checkboxInput(
                          inputId = "find_all_markers",
                          label = "Enable Group 1 vs All Other Cells",
                          value = FALSE
                        ),
                        fluidRow(tags$br()),
                        uiOutput("diffOut2_l_tp")),

                 column(width = 12, tags$br()),

                 column(width = 12,
                        align = "center",
                        uiOutput("cellSelectDiff_l_tp")),

                 # column(width = 12,
                 #        align = "center",
                 #        tags$hr(width = "50%")),
                 # 
                 # fluidRow(tags$br()),

                 fluidRow(tags$br()),

                 column(width = 12,
                        align = "center",
                        numericInput(
                          "pValCutoff_l_tp",
                          "Input adjusted p-value cutoff:",
                          value = 0.05,
                          min = 0.00,
                          step = 0.001,
                          max = 1.00,
                          width = "210px"
                        )),

                 fluidRow(tags$br()),

                 column(width = 12,
                        align = "center",
                        actionButton("runDiffExp_l_tp", "Run Differential Expression",
                                     style = "padding:5px; font-size:80%")),

                 column(width = 12,
                        align = "center",
                        tags$hr(width = "50%")),

                 column(width = 12,
                        align = "center",
                        downloadButton("downloadDiffExp_l_tp",
                                       "Download Results",
                                       style = "padding:5px; font-size:80%")),

                 fluidRow(tags$br()),
                 fluidRow(tags$br()),

                 fluidRow(column(
                   width = 12,
                   align = "center",
                   tags$img(
                     src = img4,
                     class = "fluid-image"
                   )
                 )),

                 fluidRow(tags$br()),
                 fluidRow(tags$br())
               ),

               mainPanel(

                 fluidRow(
                   column(width = 8, tags$br()),

                   column(width = 12,
                          align = "left",
                          class = "diffExpMain",
                          column(
                            width = 12,
                            align = "left",
                            fluidRow(tags$br()),
                            tags$b(
                              HTML("Differential expression analysis was performed using the Wilcoxon Rank test. <br> Genes are ranked in descending order based on their average log2 fold change (avg_log2FC)."))
                          ),
                          # column(
                          #   width = 12,
                          #   tags$br()
                          # ),
                          column(width = 12,
                                 tags$hr()),
                          column(
                            width = 12,
                            tags$br()
                          ),
                          uiOutput("diffTable_l_tp")
                   )
                 )

               )
             ))
  ))
## =========================================================================== ##

## ============================== define server ============================== ##
server <- function(input, output) {
  #===== UMAP tab =====#
  #UMAP tab - Select Cell Subclusters
  observeEvent(input$lauchApp, {
    selectedApp <- input$Cell_sub
    
    appURLs <- list(
      "Follicular somatic cells" = "https://simrcompbio.shinyapps.io/AstMex_ovary_folsom/"
    )
    
    if (selectedApp %in% names(appURLs)) {
      url <- appURLs[[selectedApp]]
      runjs(sprintf("window.open('%s', '_blank');", url))
    }
  })
  
  # observeEvent(input$lauchApp, {
  #   selectedApp <- input$Cell_sub
  #   
  #   appURLs <- list(
  #     #"Germ cells" = "http://shiny:3431/dw2733/Rhoner_cavefish_germ_cells/",
  #     "Follicular somatic cells" = "http://shiny:3431/dw2733/Rhoner_cavefish_ovarian_follicle/"#,
  #     #"Immune cells" = "http://shiny:3431/dw2733/Rhoner_cavefish_immune_cells/"
  #   )
  #   
  #   if (selectedApp %in% names(appURLs)) {
  #     url <- appURLs[[selectedApp]]
  #     #browseURL(url = unlist(url), browser = getOption("browser"), encodeIfNeeded = TRUE)
  #     #IT said it is not safe to invoke a browser, so we will do below
  #     output$appLink <- renderUI({
  #       tags$a("Click here to launch the app",
  #              href = url,
  #              target = "_blank")
  #     })
  #   }
  # })
  
  #===== FeaturePlot tab =====#
  #FeaturePlot tab - notInFeat
  mismatchFeat <- function() {
    #selected <- unique(unlist(strsplit(input$featureGenes, " ")))
    selected <- tolower(unlist(unique(strsplit(input$featureGenes, " "))))
    
    mismatch <- setdiff(selected, tolower(rownames(seurat_obj)))
    
    return(mismatch)
  }
  
  output$notInFeat <- renderText({
    input$runFeatPlot
    isolate({
      mismatchFeat()
    })
  })
  
  #FeaturePlot tab - cellSelectFeat
  output$cellSelectFeat <- renderUI({
    choices <- c("Surface fish" = "Surface_fish", 
                 "Pachón" = "Pachon", 
                 "Molino" = "Molino")
    
    pickerInput(
      "cellIdentsFeat",
      "Add or remove dataset(s) from plot:",
      choices = choices,             
      multiple = TRUE,
      selected = choices,            
      options = list(`actions-box` = TRUE),
      width = "80%"
    )
  })
  
  #FeaturePlot tab - cellSelectCol
  output$cellSelectCol <- renderUI({
    choices <- c("Purple Gradient" = "Purple",
                 "Red Gradient" = "Red")
    pickerInput(
      "colorScheme",
      "Select color scheme:",
      choices = choices,
      selected = "Purple",
      options = list(`actions-box` = TRUE),
      width = "70%"
    )
  })
  
  #FeaturePlot tab - plot.uiFeaturePlotF
  FeaturePlotF <- reactive({
    # Define color gradients
    # purple_gradient <- rev(c("#231942", "#5E548E", "#9F86C0", "#BE95C4", "#E0B1CB", "#CED4DA"))
    # red_gradient <- rev(c("#461220", "#8C2F39", "#B23A48", "#FCB9B2", "#FED0BB", "#CED4DA"))
    
    # Determine colors based on the selection
    colors <- switch(input$colorScheme,
                     "Purple" = purple_gradient,
                     "Red" = red_gradient)
    
    #to make the gene name case insensitive, we will make both the names in the selected and rownames of seuobj lowercase
    selected <- unique(unlist(strsplit(input$featureGenes, " ")))
    
    for (j in seq_along(selected)){
      selected[j] <- ifelse(selected[j] %in% rownames(seurat_obj), selected[j], 
                            grep(pattern = paste0("^", selected[j], "$"), x = rownames(seurat_obj), ignore.case = TRUE, value = TRUE))
    }
    selected <- na.omit(selected)
    
    cells_to_plt <- names(seurat_obj$fish[seurat_obj$fish %in% input$cellIdentsFeat])
    
    feat <- FeaturePlot(
      object = seurat_obj,
      features = selected,
      raster = FALSE,
      cells = cells_to_plt,
      order = TRUE,
      slot = "data",
      reduction = "umap",
      #cols = colors,
      combine = FALSE,
      pt.size = input$ptSizeFeature
    )
    
    for (k in 1:length(feat)) {
      feat[[k]] <- feat[[k]] +
        scale_colour_gradientn(colors = colors) +
        coord_fixed(ratio = 1) +
        labs(x = "UMAP 1", y = "UMAP 2") +
        theme(
          text = element_text(family = "Helvetica", size = 18),
          axis.text.x = element_blank(),
          legend.position = "right",
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          axis.title = element_text(size = 18),
          panel.border = element_rect(
            colour = "#FFFFFF",
            fill = NA,
            linewidth = 1
          )
        )
    }
    
    pg <- plot_grid(plotlist = feat, ncol = 1) +
      theme(plot.title = element_text(
        face = "bold",
        size = 15,
        hjust = 0
      ))
    
    return(pg)
  })
  
  output$myFeaturePlotF <- renderPlot({
    input$runFeatPlot
    isolate({
      withProgress({
        p <- FeaturePlotF()
        print(p)
      },
      message = "Rendering plot..", min = 0, max = 10, value = 10)
    })
  })
  
  output$plot.uiFeaturePlotF <- renderUI({
    input$runFeatPlot
    isolate({
      # Convert both to lowercase for case-insensitive matching
      features <- tolower(unique(unlist(strsplit(input$featureGenes, " "))))
      feature_names <- tolower(rownames(seurat_obj))
      
      # Check case-insensitively if features are in the seurat_obj row names
      valid_features <- features[features %in% feature_names]
      
      plotOutput("myFeaturePlotF",
                 width = "840px",
                 height = length(valid_features) * 800)  # Set height based on the number of matched features
    })
  })
  
  output$downloadFeaturePlotF <- downloadHandler(
    filename = "Feature_plot.pdf",
    content = function(file) {
      # Convert input features and rownames to lowercase for case-insensitive matching
      feat_N <- tolower(unique(unlist(strsplit(input$featureGenes, " "))))
      feature_names <- tolower(rownames(seurat_obj))
      
      # Match features in a case-insensitive manner
      valid_feat_N <- feat_N[feat_N %in% feature_names]
      
      # Set PDF height based on the number of matched features
      pdf(file,
          width = 14,
          height = 10 * length(valid_feat_N))
      print(FeaturePlotF())
      dev.off()
    }
  )
  
  #===== Split FeaturePlot tab =====#
  #Split FeaturePlot tab - notInFeat_split
  
  mismatchFeat_split <- function() {
    selected <- tolower(input$featureGenes_split)
    
    mismatch <- setdiff(selected, tolower(rownames(seurat_obj)))
    
    return(mismatch)
  }
  
  output$notInFeat_split <- renderText({
    input$runFeatPlot_split
    isolate({
      mismatchFeat_split()
    })
  })
  
  #FeaturePlot tab - cellSelectCol_tab3
  output$cellSelectCol_tab3 <- renderUI({
    choices <- c("Purple Gradient" = "Purple_tab3",
                 "Red Gradient" = "Red_tab3")
    pickerInput(
      "colorScheme_tab3",
      "Select color scheme:",
      choices = choices,
      selected = "Purple_tab3",
      options = list(`actions-box` = TRUE),
      width = "70%"
    )
  })
  
  #Split FeaturePlot tab - plot.uiFeaturePlotF_split
  FeaturePlotF_split <- reactive({
    # Determine colors based on the selection
    colors <- switch(input$colorScheme_tab3,
                     "Purple_tab3" = purple_gradient,
                     "Red_tab3" = red_gradient)
    
    #to make the gene name case insensitive, we will make both the names in the selected and rownames of seuobj lowercase
    selected <- input$featureGenes_split
    
    selected <- ifelse(selected %in% rownames(seurat_obj), selected, 
                            grep(pattern = paste0("^", selected, "$"), x = rownames(seurat_obj), ignore.case = TRUE, value = TRUE))
    
    selected <- na.omit(selected)
    
    # cells_to_plt <-
    #   names(seurat_obj$orig.ident[seurat_obj$orig.ident %in% input$cellIdentsFeat_split])
    
    seurat_obj$orig.ident <- factor(seurat_obj$orig.ident,
                                    levels = c("Surface_fish_1", "Surface_fish_2", "Pachon_1", "Pachon_2", "Molino_1", "Molino_2"),
                                    labels = c("Surface fish 1", "Surface fish 2", "Pachón 1", "Pachón 2", "Molino 1", "Molino 2"))
    
    feat <- FeaturePlot(
      object = seurat_obj,
      features = selected,
      #cells = cells_to_plt,
      order = TRUE,
      slot = "data",
      #cols = colors,
      combine = FALSE,
      split.by = "orig.ident",
      pt.size = input$ptSizeFeature_split,
      ncol = 1
    )
    
    for (k in 1:length(feat)) {
      feat[[k]] <- feat[[k]] +
        scale_colour_gradientn(colors = colors) +
        coord_fixed() +
        labs(x = "", y = "") +
        theme(
          text = element_text(family = "Helvetica", size = 18),
          legend.position = "left",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          axis.title = element_text(size = 20),
          panel.border = element_rect(
            colour = "#FFFFFF",
            fill = NA,
            size = 1
          )
        )
    }
    
    pg <- plot_grid(plotlist = feat, ncol = 1) +
      theme(plot.title = element_text(
        face = "bold",
        size = 15,
        hjust = 0
      ))
    
    return(pg)
  })
  
  output$myFeaturePlotF_split <- renderPlot({
    input$runFeatPlot_split
    isolate({
      withProgress({
        p <- FeaturePlotF_split()
        
        width <- "840px"
        #height <- length(input$cellIdentsFeat_split) * 800
        height <- length(unique(seurat_obj$orig.ident)) * 800
        
        options(repr.plot.width = width, repr.plot.height = height)
        
        print(p)
      },
      message = "Rendering plot..", min = 0, max = 10, value = 10)
    })
  })
  
  output$plot.uiFeaturePlotF_split <-
    renderUI({
      plotOutput(
        "myFeaturePlotF_split",
        width = "840px",
        #height = length(input$cellIdentsFeat_split) * 800
        height = length(unique(seurat_obj$orig.ident)) * 800
      )
    })
  
  #Split FeaturePlot tab - downloadFeaturePlotF
  output$downloadFeaturePlotF_split <- downloadHandler(
    filename = "Feature_plot.pdf",
    content = function(file) {
      # feat_N <-
      #   unique(unlist(strsplit(
      #     as.character(input$cellIdentsFeat_split), " "
      #   )))
      feat_N <- unique(seurat_obj$orig.ident)
      
      pdf(file,
          width = 14,
          height = 10 * length(feat_N))
      print(FeaturePlotF_split())
      dev.off()
    }
  )
  
  #===== Split FeaturePlot by fish tab =====#
  #Split FeaturePlot tab - notInFeat_split2
  mismatchFeat_split2 <- function() {
    selected <- tolower(input$featureGenes_split2)
    
    mismatch <- setdiff(selected, tolower(rownames(seurat_obj)))
    
    return(mismatch)
  }
  
  output$notInFeat_split2 <- renderText({
    input$runFeatPlot_split2
    isolate({
      mismatchFeat_split2()
    })
  })
  
  # #Split FeaturePlot tab - cellSelectFeat_split2
  # output$cellSelectFeat_split2 <- renderUI({
  #   pickerInput(
  #     "cellIdentsFeat_split2",
  #     "Add or remove samples from plot:",
  #     choices = as.character(unique(seurat_obj$fish)),
  #     multiple = TRUE,
  #     selected = as.character(unique(seurat_obj$fish)),
  #     options = list(`actions-box` = TRUE)
  #   )
  # })
  
  #FeaturePlot by fish tab - cellSelectCol_tab4
  output$cellSelectCol_tab4 <- renderUI({
    choices <- c("Purple Gradient" = "Purple_tab4",
                 "Red Gradient" = "Red_tab4")
    pickerInput(
      "colorScheme_tab4",
      "Select color scheme:",
      choices = choices,
      selected = "Purple_tab4",
      options = list(`actions-box` = TRUE),
      width = "70%"
    )
  })
  
  #Split FeaturePlot tab - plot.uiFeaturePlotF_split2
  FeaturePlotF_split2 <- reactive({
    # purple_gradient4 <- rev(c("#231942", "#5E548E", "#9F86C0", "#BE95C4", "#E0B1CB", "#CED4DA"))
    # red_gradient4 <- rev(c("#461220", "#8C2F39", "#B23A48", "#FCB9B2", "#FED0BB", "#CED4DA"))
    
    # Determine colors based on the selection
    colors <- switch(input$colorScheme_tab4,
                     "Purple_tab4" = purple_gradient,
                     "Red_tab4" = red_gradient)
    
    #to make the gene name case insensitive, we will make both the names in the selected and rownames of seuobj lowercase
    selected <- input$featureGenes_split2
    
    selected <- ifelse(selected %in% rownames(seurat_obj), selected, 
                       grep(pattern = paste0("^", selected, "$"), x = rownames(seurat_obj), ignore.case = TRUE, value = TRUE))
    
    selected <- na.omit(selected)
    
    seurat_obj$fish <- factor(seurat_obj$fish,
                              levels = c("Surface_fish", "Pachon", "Molino"),
                              labels = c("Surface fish", "Pachón", "Molino"))
    # cells_to_plt2 <-
    #   names(seurat_obj$fish[seurat_obj$fish %in% input$cellIdentsFeat_split2])
    
    feat2 <- FeaturePlot(
      object = seurat_obj,
      features = selected,
      #cells = cells_to_plt2,
      order = TRUE,
      slot = "data",
      reduction = "umap",
      #cols = colors,
      combine = FALSE,
      split.by = "fish",
      pt.size = input$ptSizeFeature_split2,
      ncol = 1
    )
    
    for (k in 1:length(feat2)) {
      feat2[[k]] <- feat2[[k]] +
        scale_colour_gradientn(colors = colors) +
        coord_fixed() +
        labs(x = "", y = "") +
        theme(
          text = element_text(family = "Helvetica", size = 18),
          axis.text.x = element_blank(),
          legend.position = "left",
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          axis.title = element_text(size = 20),
          panel.border = element_rect(
            colour = "#FFFFFF",
            fill = NA,
            linewidth = 1
          )
        )
    }
    
    pg <- plot_grid(plotlist = feat2, ncol = 1) +
      theme(plot.title = element_text(
        face = "bold",
        size = 15,
        hjust = 0
      ))
    
    return(pg)
  })
  
  output$myFeaturePlotF_split2 <- renderPlot({
    input$runFeatPlot_split2
    isolate({
      withProgress({
        p <- FeaturePlotF_split2()
        
        width <- "840px"
       # height <- length(input$cellIdentsFeat_split2) * 800
        height <- length(unique(seurat_obj$fish)) * 800
        
        options(repr.plot.width = width, repr.plot.height = height)
        
        print(p)
      },
      message = "Rendering plot..", min = 0, max = 10, value = 10)
    })
  })
  
  output$plot.uiFeaturePlotF_split2 <-
    renderUI({
      plotOutput(
        "myFeaturePlotF_split2",
        width = "840px",
        #height = length(input$cellIdentsFeat_split2) * 800
        height = length(unique(seurat_obj$fish)) * 800
      )
    })
  
  #Split FeaturePlot tab - downloadFeaturePlotF
  output$downloadFeaturePlotF_split2 <- downloadHandler(
    filename = "Feature_plot.pdf",
    content = function(file) {
      # feat_N2 <-
      #   unique(unlist(strsplit(
      #     as.character(input$cellIdentsFeat_split2), " "
      #   )))
      feat_N2 <- unique(seurat_obj$fish)
      
      pdf(file,
          width = 14,
          height = 10 * length(feat_N2))
      print(FeaturePlotF_split2())
      dev.off()
    }
  )
  
  #===== VlnPlot tab =====#
  #VlnPlot tab - notInVln
  mismatchVln <- function() {
    selected <- tolower(unique(unlist(strsplit(input$vlnGenes, " "))))
    
    mismatch <- setdiff(selected, tolower(rownames(seurat_obj)))
    
    return(mismatch)
  }
  
  output$notInVln <- renderText({
    input$runVlnPlot
    isolate({
      mismatchVln()
    })
  })
  
  #VlnPlot tab - plot.uiVlnPlotF
  VlnPlotF <- reactive({
    selected <- unique(unlist(strsplit(input$vlnGenes, " ")))
    
    for (j in seq_along(selected)){
      selected[j] <- ifelse(selected[j] %in% rownames(seurat_obj), selected[j], 
                            grep(pattern = paste0("^", selected[j], "$"), x = rownames(seurat_obj), ignore.case = TRUE, value = TRUE))
    }
    selected <- na.omit(selected)
    
    #define colors for each group
    vln_colors <- switch(input$selectGrpVln,
                         "annotation" = colors_umap,
                         "seurat_clusters" = colors_to_use,
                         "fish" = colors_fish_group,
                         "orig.ident" = colors_sample)
    
    seurat_obj$fish <- factor(seurat_obj$fish,
                              levels = c("Surface_fish", "Pachon", "Molino"),
                              labels = c("Surface fish", "Pachón", "Molino"))
    
    seurat_obj$orig.ident <- factor(seurat_obj$orig.ident,
                                    levels = c("Surface_fish_1", "Surface_fish_2", "Pachon_1", "Pachon_2", "Molino_1", "Molino_2"),
                                    labels = c("Surface fish 1", "Surface fish 2", "Pachón 1", "Pachón 2", "Molino 1", "Molino 2"))
    
    g <- VlnPlot(
      object = seurat_obj,
      features = selected,
      pt.size = input$ptSizeVln,
      group.by = input$selectGrpVln,
      cols = vln_colors,
      combine = FALSE
    )
    
    for (k in 1:length(g)) {
      g[[k]] <- g[[k]] + theme(legend.position = "none",
                               text = element_text(family = "Helvetica", size = 18)) +#,
                               #axis.text.x = element_text(family = "Helvetica", angle = 90, hjust = 1, size = 10)) +
        xlab("")
    }
    
    pg <- plot_grid(plotlist = g, ncol = 1)
    
    return(pg)
  })
  
  output$myVlnPlotF <- renderPlot({
    input$runVlnPlot
    isolate({
      withProgress({
        p <- VlnPlotF()
        print(p)
      },
      message = "Rendering plot..",
      min = 0, max = 10, value = 10)
    })
  })
  
  output$plot.uiVlnPlotF <- renderUI({
    input$runVlnPlot
    isolate({
      # Convert both to lowercase for case-insensitive matching
      features <- tolower(unique(unlist(strsplit(input$vlnGenes, " "))))
      feature_names <- tolower(rownames(seurat_obj))
      
      # Check case-insensitively if features are in the seurat_obj row names
      valid_features <- features[features %in% feature_names]
      
      # feat_L <- unique(unlist(strsplit(input$vlnGenes, " ")))
      plotOutput("myVlnPlotF",
                 width = "950px",
                 height = length(valid_features) * 500)
    })
  })
  
  #VlnPlot tab - downloadVlnPlot
  output$downloadVlnPlot <- downloadHandler(
    filename = "Violin_plot.pdf",
    content = function(file) {
      
      # Convert input features and rownames to lowercase for case-insensitive matching
      feat_L <- tolower(unique(unlist(strsplit(input$vlnGenes, " "))))
      feature_names <- tolower(rownames(seurat_obj))
      
      # Match features in a case-insensitive manner
      valid_feat_L <- feat_L[feat_L %in% feature_names]
      
      pdf(file,
          width = 14,
          height = 10 * length(valid_feat_L))
      print(VlnPlotF())
      dev.off()
    }
  )
  
  #===== Dot plot tab =====#
  #Dot plot tab - notInDot
  mismatchFeat_dot <- function() {
    selected <- tolower(unique(unlist(strsplit(input$dotGenes, " "))))
    
    mismatch <- setdiff(selected, tolower(rownames(seurat_obj)))
    
    return(mismatch)
  }
  
  output$notInDot <- renderText({
    input$runDotPlot
    isolate({
      mismatchFeat_dot()
    })
  })
  
  #Dot plot tab - plot.uiDotPlotF
  DotPlotF <- reactive({
    clustering_fea <- input$dPlotClust_fea
    clustering_ident <- input$dPlotClust_ident
    
    seurat_obj$fish <- factor(seurat_obj$fish,
                              levels = c("Surface_fish", "Pachon", "Molino"),
                              labels = c("Surface fish", "Pachón", "Molino"))
    
    seurat_obj$orig.ident <- factor(seurat_obj$orig.ident,
                                    levels = c("Surface_fish_1", "Surface_fish_2", "Pachon_1", "Pachon_2", "Molino_1", "Molino_2"),
                                    labels = c("Surface fish 1", "Surface fish 2", "Pachón 1", "Pachón 2", "Molino 1", "Molino 2"))
    
    if (!clustering_fea & !clustering_ident) {
      
      # no clustering
      # selected <- unique(unlist(strsplit(input$dotGenes, " ")))
      selected <- caseinsenFun(input$dotGenes, seurat_obj)
      
      p <- DotPlot_scCustom(
        seurat_object = seurat_obj,
        features = selected,
        x_lab_rotate = TRUE,
        group.by = input$selectGrpDot,
        flip_axes = TRUE,
        dot.scale = input$dotScale
      ) +
        scale_color_viridis_c(alpha = 0.8, direction = -1) +
        theme(text = element_text(family = "Helvetica", size = 18),
              legend.text = element_text(size = input$LegLabSize),
              legend.title = element_text(size = input$LegTitSize),
              axis.text.x = element_text(size = input$identSize, angle = 70),
              axis.text.y = element_text(size = input$FeatSize))
      
    } else {
      
      if (clustering_fea & !clustering_ident) {
        
        # cluster feature
        #selected <- unique(unlist(strsplit(input$dotGenes, " ")))
        selected <- caseinsenFun(input$dotGenes, seurat_obj)
        
        p <- Clustered_DotPlot(
          seurat_object = seurat_obj,
          features = selected,
          x_lab_rotate = TRUE,
          group.by = input$selectGrpDot,
          flip = TRUE,
          cluster_feature = TRUE,
          cluster_ident = FALSE,
          colors_use_idents = colors_to_use,
          colors_use_exp = viridis_dark_high,
          plot_km_elbow = FALSE,
          row_label_size = input$identSize,
          column_label_size = input$FeatSize,
          legend_label_size = input$LegLabSize,
          legend_title_size = input$LegTitSize
        ) 
        
      } else {
        
        if (!clustering_fea & clustering_ident){
          
          #cluster ident
          # selected <- unique(unlist(strsplit(input$dotGenes, " ")))
          selected <- caseinsenFun(input$dotGenes, seurat_obj)
          
          p <- Clustered_DotPlot(
            seurat_object = seurat_obj,
            features = selected,
            x_lab_rotate = TRUE,
            group.by = input$selectGrpDot,
            flip = TRUE,
            cluster_feature = FALSE,
            cluster_ident = TRUE,
            colors_use_idents = colors_to_use,
            colors_use_exp = viridis_dark_high,
            plot_km_elbow = FALSE,
            row_label_size = input$identSize,
            column_label_size = input$FeatSize,
            legend_label_size = input$LegLabSize,
            legend_title_size = input$LegTitSize
          ) 
        } else{
          
          # cluster both
          # selected <- unique(unlist(strsplit(input$dotGenes, " ")))
          selected <- caseinsenFun(input$dotGenes, seurat_obj)
          
          p <- Clustered_DotPlot(
            seurat_object = seurat_obj,
            features = selected,
            x_lab_rotate = TRUE,
            group.by = input$selectGrpDot,
            flip = TRUE,
            cluster_feature = TRUE,
            cluster_ident = TRUE,
            colors_use_idents = colors_to_use,
            colors_use_exp = viridis_dark_high,
            plot_km_elbow = FALSE,
            row_label_size = input$identSize,
            column_label_size = input$FeatSize,
            legend_label_size = input$LegLabSize,
            legend_title_size = input$LegTitSize
          )
        }
      }
    }
    
    return(p)
  })
  
  output$myDotPlotF <- renderPlot({input$runDotPlot
    isolate({withProgress({p <- DotPlotF(); print(p)},
                          message = "Rendering plot..", min = 0, max = 10, value = 10)
    })
  })
  
  output$plot.uiDotPlotF <- renderUI({input$runDotPlot
    isolate({plotOutput("myDotPlotF",
                        width = paste0(input$manAdjustDotW, "in"),
                        height = paste0(input$manAdjustDotH, "in"))})
  })
  
  output$downloadDotPlot <- downloadHandler(
    filename = "dot_plot.pdf", content = function(file) {
      pdf(file, height = as.numeric(input$manAdjustDotH +3),
          width = as.numeric(input$manAdjustDotW +3))
      print(DotPlotF())
      dev.off()
    }
  )
 
  #===== Dot plot by fish tab =====# 
  #Dot plot by fish tab - notInDotSplit
  mismatchFeat_dotS <- function() {
    selected <- tolower(unique(unlist(strsplit(input$dotGenesSplit, " "))))
    
    mismatch <- setdiff(selected, tolower(rownames(seurat_obj)))
    
    return(mismatch)
  }
  
  output$notInDotSplit <- renderText({
    input$runDotPlotSplit
    isolate({
      mismatchFeat_dotS()
    })
  })
  
  #Dot plot by fish tab - fishSelect
  output$fishSelect <- renderUI({
    choices <- c("Surface fish" = "Surface_fish",
                 "Pachón" = "Pachon",
                 "Molino" = "Molino")
    pickerInput(
      "fishIdents",
      "Add or remove fish group from plot:",
      choices = choices,
      multiple = TRUE,
      selected = choices,
      options = list(`actions-box` = TRUE)
    )
  })
  
  #Dot plot by fish tab - CluSelect
  output$CluSelect <- renderUI({
    pickerInput(
      "cluIdents",
      "Add or remove cluster(s) from plot:",
      choices = as.character(0:(length(unique(seurat_obj$seurat_clusters)) - 1)),
      multiple = TRUE,
      selected = as.character(0:(length(unique(seurat_obj$seurat_clusters)) - 1)),
      options = list(`actions-box` = TRUE)
    )
  })
  
  #Dot plot by fish tab - plot.uiDotPlotSplitF
  DotPlotF_split <- reactive({
    #selected <- unique(unlist(strsplit(input$dotGenesSplit, " ")))
    selected <- caseinsenFun(input$dotGenesSplit, seurat_obj)
    
    fish_to_plt <-
      names(seurat_obj$fish[seurat_obj$fish %in% input$fishIdents])
    
    Idents(seurat_obj) <- seurat_obj$seurat_clusters
    
    seurat_obj$fish <- factor(seurat_obj$fish,
                              levels = c("Surface_fish", "Pachon", "Molino"),
                              labels = c("Surface fish", "Pachón", "Molino"))
    
    p <- DotPlot_scCustom(
      seurat_object = seurat_obj[ , fish_to_plt],
      features = selected,
      x_lab_rotate = TRUE,
      group.by = "seurat_clusters",
      flip_axes = TRUE,
      dot.scale = input$dotScaleSplit,
      idents = input$cluIdents
    ) +
      scale_color_viridis_c(alpha = 0.8, direction = -1) +
      theme(text = element_text(family = "Helvetica", size = 18),
            legend.text = element_text(size = input$LegLabSizeSplit),
            legend.title = element_text(size = input$LegTitSizeSplit),
            axis.text.x = element_text(size = input$identSizeSplit, angle = 70),
            axis.text.y = element_text(size = input$FeatSizeSplit)) +
      ggtitle(paste(input$fishIdents, collapse = ", "))
    
    return(p)
  })
  
  output$myDotPlotF_split <- renderPlot({
    input$runDotPlotSplit
    isolate({
      withProgress({
        p <- DotPlotF_split()
        print(p)
      },
      message = "Rendering plot..", min = 0, max = 10, value = 10)
    })
  })
  
  output$plot.uiDotPlotSplitF <- renderUI({
    input$runDotPlotSplit
    isolate({
      plotOutput("myDotPlotF_split",
                 width = paste0(input$manAdjustDotWSplit, "in"),
                 height = paste0(input$manAdjustDotHSplit, "in"))
    })
  })
  
  output$downloadDotPlotSplit <- downloadHandler(
    filename = "dot_plot_fish.pdf",
    content = function(file) {
      pdf(file, height = as.numeric(input$manAdjustDotHSplit +3),
          width = as.numeric(input$manAdjustDotWSplit +3))
      print(DotPlotF_split())
      dev.off()
    }
  )
  
  #===== cluster DE fish tab =====#
  #cluster DE fish tab - diffTable
  diffExp_less <- reactive({
    #less clusters: fish
    meta <- seurat_obj@meta.data
    
    meta <- meta[meta$annotation %in% input$cellIdentsDiff_less, ]
    
    subset1 <- as.character(input$identText1_less)
    subset2 <- as.character(input$identText2_less)
    
    group1 <- rownames(meta[meta$fish %in% subset1, ])
    group2 <- rownames(meta[meta$fish %in% subset2, ])
    
    diff_results <- FindMarkers(object = seurat_obj,
                                assay = "SCT",
                                slot = "data",
                                ident.1 = group1, 
                                ident.2 = group2)
    
    pval <- as.numeric(input$pValCutoff_less)
    
    diff_results$gene <- rownames(diff_results)
    
    diff_results <- diff_results %>%
      dplyr::filter(p_val_adj <= pval) %>%
      arrange(desc(avg_log2FC)) %>%
      select(gene, pct.1, pct.2, avg_log2FC, p_val_adj, p_val)
    
    #add ensembl ID
    diff_results <- left_join(diff_results, res_j %>%
                                dplyr::select(gene_symbol_uniq, ens_id), by = c("gene" = "gene_symbol_uniq"))
    diff_results %<>%
      rename(ensembl_id = ens_id)
  }
  )
  
  #requires input$identText to execute before diffExp_less()
  diffReact_less <-
    eventReactive(c(input$identText1_less, input$identText2_less, input$cellIdentsDiff_less), diffExp_less())
  
  output$diffTable_less <- renderTable({
    input$runDiffExp_less
    isolate({
      withProgress(
        diffReact_less(),
        message = "Calculating..",
        min = 0,
        max = 10,
        value = 10
      )
    })
  },
  digits = -5,
  spacing = "s",
  bordered = TRUE,
  rownames = FALSE)
  
  #cluster DE fish tab - diffOut1_less and 2
  output$diffOut1_less <- renderUI({
    choices = c("Surface fish" = "Surface_fish",
                "Pachón" = "Pachon",
                "Molino" = "Molino")
    pickerInput(
      "identText1_less",
      tags$b("Group 1 - positive FC"),
      choices = choices,
      multiple = TRUE,
      selected = "Surface_fish",
      options = list(`actions-box` = TRUE),
      width = "80%"
    )
  })
  
  output$diffOut2_less <- renderUI({
    choices = c("Surface fish" = "Surface_fish",
                "Pachón" = "Pachon",
                "Molino" = "Molino")
    pickerInput(
      "identText2_less",
      tags$b("Group 2 - negative FC"),
      choices = choices,
      multiple = TRUE,
      selected = "Pachon",
      options = list(`actions-box` = TRUE),
      width = "80%"
    )
  })
  
  #cluster DE fish tab - cellSelectDiff_less
  output$cellSelectDiff_less <- renderUI({
    pickerInput(
      "cellIdentsDiff_less",
      "Add or remove cell type(s):",
      choices = as.character(unique(seurat_obj$annotation)),
      multiple = TRUE,
      selected =  as.character(unique(seurat_obj$annotation)),
      options = list(`actions-box` = TRUE),
      width = "80%"
    )
  })          
  
  #cluster DE fish tab - downloadDiffExp_less
  output$downloadDiffExp_less <- downloadHandler(
    filename = "DE_results.tsv",
    content = function(file) {
      write.table(
        diffReact_less(),
        file,
        row.names = FALSE,
        col.names = TRUE,
        qmethod = "double",
        sep = "\t"
      )
    }
  )
  
  #===== cluster DE cell type tab =====#
  #cluster DE cell type tab - diffTable
  diffExp_l_tp <- reactive({
    meta <- seurat_obj@meta.data

    meta <- meta[meta$fish %in% input$cellIdentsDiff_l_tp, ]

    subset1 <- as.character(input$identText1_l_tp)
    subset2 <- as.character(input$identText2_l_tp)
    subset2_grp2disabled <- setdiff(as.character(unique(seurat_obj$annotation)), subset1)

    group1 <- rownames(meta[meta$annotation %in% subset1, ])

    if (input$find_all_markers) {
      group2 <- rownames(meta[meta$annotation %in% subset2_grp2disabled, ])
    } else {
      group2 <- rownames(meta[meta$annotation %in% subset2, ])
    }

    diff_results <- FindMarkers(object = seurat_obj,
                                assay = "SCT",
                                slot = "data",
                                ident.1 = group1,
                                ident.2 = group2)

    pval <- as.numeric(input$pValCutoff_l_tp)

    diff_results$gene <- rownames(diff_results)

    diff_results <- diff_results %>%
      dplyr::filter(p_val_adj <= pval) %>%
      arrange(desc(avg_log2FC)) %>%
      select(gene, pct.1, pct.2, avg_log2FC, p_val_adj, p_val)
    
    #add ensembl ID
    diff_results <- left_join(diff_results, res_j %>%
                                dplyr::select(gene_symbol_uniq, ens_id), by = c("gene" = "gene_symbol_uniq"))
    diff_results %<>%
      rename(ensembl_id = ens_id)
  }
  )

  #requires input$identText to execute before diffExp_l_tp()
  diffReact_l_tp <-
    eventReactive(c(input$identText1_l_tp, input$identText2_l_tp, input$cellIdentsDiff_l_tp), diffExp_l_tp())

  output$diffTable_l_tp <- renderTable({
    input$runDiffExp_l_tp
    isolate({
      withProgress(
        diffReact_l_tp(),
        message = "Calculating..",
        min = 0,
        max = 10,
        value = 10
      )
    })
  },
  digits = -5,
  spacing = "s",
  bordered = TRUE,
  rownames = FALSE)

  #cluster DE cell type tab - diffOut1_l_tp and 2
  output$diffOut1_l_tp <- renderUI({
    pickerInput(
      "identText1_l_tp",
      tags$b("Group 1 - positive FC"),
      choices = as.character(unique(seurat_obj$annotation)),
      multiple = TRUE,
      selected = as.character(unique(seurat_obj$annotation))[1],
      options = list(`actions-box` = TRUE),
      width = "80%"
    )
  })

  output$diffOut2_l_tp <- renderUI({
    if (input$find_all_markers) {
      pickerInput(
        "identText2_l_tp",
        tags$b("Group 2 - negative FC (Disabled)"),
        choices = as.character(NA),
        multiple = TRUE,
        selected = as.character(NA),
        options = list(`actions-box` = TRUE),
        width = "80%"
      )
    } else {
      pickerInput(
        "identText2_l_tp",
        tags$b("Group 2 - negative FC"),
        choices = as.character(unique(seurat_obj$annotation)),
        multiple = TRUE,
        selected = as.character(unique(seurat_obj$annotation))[2],
        options = list(`actions-box` = TRUE),
        width = "80%"
      )
    }
  })

  #cluster DE cell type tab - cellSelectDiff_l_tp
  #this is less clusters tab: annotation
  output$cellSelectDiff_l_tp <- renderUI({
    choices = c("Surface fish" = "Surface_fish",
                "Pachón" = "Pachon",
                "Molino" = "Molino")
    pickerInput(
      "cellIdentsDiff_l_tp",
      "Add or remove dataset(s):",
      choices = choices,
      multiple = TRUE,
      selected = choices,
      options = list(`actions-box` = TRUE),
      width = "80%"
    )
  })

  #DE timepoint DE cell type tab - downloadDiffExp_l_tp
  output$downloadDiffExp_l_tp <- downloadHandler(
    filename = "DE_results.tsv",
    content = function(file) {
      write.table(
        diffReact_l_tp(),
        file,
        row.names = FALSE,
        col.names = TRUE,
        qmethod = "double",
        sep = "\t"
      )
    }
  )
}
## =========================================================================== ##



## ================================ run the app ============================== ##
shinyApp(ui = ui, server = server)
