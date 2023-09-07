library(fgsea)
library(ggplot2)
library(metap)
library(ggrepel)
library(shiny)
library(shinybusy)
library(htmltools)
library(htmlwidgets)
library(openxlsx)
library(pheatmap)
library(shinyWidgets)
library(rintrojs)
library(shinythemes)
library(shinydashboard)

load("Data/Alldata_20Sep.RData", envir=.GlobalEnv)
load("Data/DE_GSE179379.RData", envir=.GlobalEnv)
rat_homologs<-read.csv("Data/Human rat homologs.txt")
mouse_homologs<-read.csv("Data/Human mouse homologs.txt")
metadata$Age<-as.numeric(metadata$Age)/365#TODO make this the default and remove division by 365 in server
test_list<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange>0 & DE_GSE179379$padj<0.05)]

ui <- navbarPage( 
  title = "Age reversal tester",
  
  ###### Here : insert shinydashboard dependencies ######
  header = tagList(
    useShinydashboard()
  ),
  
  tabPanel(
    
    #titlePanel("MetaLCM breast"),
    introBox("Test reversal",
             data.step = 1,
             data.intro = "In this tab, you can test the age-reversal impact of an input gene list"),
    actionButton("tour","Start Tour", class = "btn-warning btn-lg"),
    sidebarLayout(
      
      #introjsUI(),
      
      sidebarPanel(width=3,
                   fluidRow(
                     
                     
                     
                     column(6,
                            introBox(
                              
                              data.position = "top",
                              data.step = 2,
                              data.intro = "To start the analysis, select the metadata to filter the samples of human brain aging",
                              h3("Metadata filters"),   
                              
                              sliderInput("age_range", "Age range:",
                                          min = 0, max = max(as.numeric(metadata$Age), na.rm=T),
                                          value = c(0,100)),
                              awesomeCheckboxGroup(inputId = "input1.gender", 
                                                   label = "Gender", 
                                                   choices = c("Male", "Female", "Any"), 
                                                   selected = "Any"),
                              awesomeCheckboxGroup(inputId = "input1.race", 
                                                   label = "Ethnicity", 
                                                   choices = c("African American", "Caucasian", "Any"), 
                                                   selected = "Any"),
                              awesomeCheckboxGroup(inputId = "input1.suicide", 
                                                   label = "Suicide", 
                                                   choices = c("No", "Yes", "Any" ),
                                                   selected = c("No")),
                              awesomeCheckboxGroup(inputId = "input1.diagnosis", 
                                                   label = "Diagnosis", 
                                                   choices = c("Alzheimer"="AD", "Autism"="ASD", 
                                                               "Alcohol Use Disorder"="AUD",
                                                               "Bipolar Disorder"="BD",
                                                               "Cocaine Use Disorder"="CUD",
                                                               "Down Syndrome"="DS",
                                                               "Huntington"="HD",
                                                               "Healthy",
                                                               "Major Depression"="MDD",
                                                               "Obsessive Compulsive Disorder"="OCD",
                                                               "Opioid Use Disorder"="OUD",
                                                               "Schizophrenia"="SCZ",
                                                               "Any"),
                                                   selected = "Any"),
                              awesomeCheckboxGroup(inputId = "input1.brain_region", 
                                                   label = "Brain region", 
                                                   choices = c("Cerebellar Cortex"="CBC", 
                                                               "Dorsolateral Prefrontal Cortex"="DLPFC",
                                                               "Orbitofrontal Cortex"="OFC", 
                                                               "Prefrontal Cortex"="PFC",
                                                               "Superior Frontal Gyrus"="SFG", 
                                                               "Temporal Cortex"="TC",
                                                               "Ventral Frontal Cortex"="VFC",
                                                               "Any" ),
                                                   selected = "Any"),
                              h3("Other parameters"),
                              sliderInput("thr", "Minimum number of samples per dataset",
                                          min = 1, max = max(table(metadata$Dataset)),
                                          value = 10),
                              sliderInput("num_perm", "Number of permutations for GSEA p-value",
                                          min = 100, max = 10000000,
                                          value = 1000)
                              
                            ))
                    
                     
                   ),
                   
                   
                   introBox(
                     data.step = 3,
                     data.intro = "Press \"Run\" to run the analysis and each time you change settings",
                     
                     actionButton("go_heat", "Run",
                                  style="color: #fff; background-color: #3cb043; border-color: #2e6da4")
                   ),
                   
                   hr(), 
                   introBox(
                     data.step = 4,
                     data.intro = "To test the age-reversal potential of a gene list you can either
                 type gene symbols separated by commas
                 or import a csv file containing gene symbols in one single column",
                     
                     h4("Genes' selection"),
                     #to choose a set of genes for filtering
                     radioButtons("organism", h4("Organism"),
                                  choices = list("Rat" = "Rat", "Mouse" = "Mouse",
                                                 "Human" = "Human"),selected = "Rat"),
                     
                     textInput("Genes_net", label = strong("Enter gene symbols"), value = ""),
                     #column(style="padding-top:25px;", 3, actionButton("go_net", "Submit gene list"))),
                     
                     p("OR", align="center"),
                     # Input: Select a file ----
                     fileInput("file1", "Choose CSV File",
                               multiple = FALSE,
                               accept = c(".csv"))
                   ),
                   
                   introBox(
                     data.step = 5,
                     data.intro = "Press \"Submit gene list\" to import gene symbols. Press \"Reset\" to restore the analysis using the example gene list",
                     
                     actionButton("go_net", "Submit gene list", style="color: #fff; background-color: #3cb043; border-color: #2e6da4"),
                     actionButton("reset_net", "Reset", style="color: #fff; background-color: #3cb043; border-color: #2e6da4")
                   )
      ),
      
      
      mainPanel(
        add_busy_spinner(spin = "double-bounce", position="full-page"),
        
        tabsetPanel(type = "tabs",
                    
                    tabPanel(
                    introjsUI(),
                    
                    introBox(
                      data.step = 6,
                      data.intro = "You can save the plots and the table with stats",
                      
                      
                      downloadButton('downloadGSEA', 'Save GSEA plot')),
                    downloadButton('downloadSummary', 'Save summary plot')),
                    downloadButton('downloadTable', 'Save stats table'),
                    
                    fluidRow(
                    
                    #textOutput("message"),
                    
                    valueBox(value=textOutput("p_ratio"), subtitle="Age reversal significance")),
                    
                    fluidRow(
                    box(title="Plot of GSEA", plotOutput("GSEA"))),
                    
                    fluidRow(
                    box(title="Summary plot of GSEA stats", plotOutput("summary")),
                    box(title="Summary table of GSEA stats", dataTableOutput("table"))
                    
                    )
                    ),
                    
                    
                    tabPanel(
                      
                      introBox("Gene-Age correlation",
                               data.step = 7,
                               data.intro = "In this tab, a Venn diagram showing the intersection between differentially expressed genes and genes in your input list will be plotted.
                                      The enrichment is tested with a one-sided Fisher test"
                               
                      ),
                      textOutput("listenrichment"),
                      plotOutput("Venn")
                      
                      
                      
                    )
        )
      )
    )
    
  )
  




server <- function(session, input, output) {
  #Tour 
  observeEvent(input$tour,
               introjs(session, options = list("nextLabel"="Next",
                                               "prevLabel"="Back"
               )))
  
  #######################################
  ##########functions for gene list inputs
  #######################################
  genes<-reactiveValues(g = NULL)
  fileGenes<-reactiveValues(g=NULL)
  user_data<-reactiveValues(selection=NULL)
  out<-reactiveValues(p=NULL, stats=NULL)
  
  ###lista di geni in input
  observeEvent(input$go_net, {
    
    if(length(unlist(strsplit(input$Genes_net, "[, ]+")))>0){
      genes$g <- strsplit(input$Genes_net, "[, ]+")
    } else if (!is.null(fileGenes$g)){
      genes$g<-(c(fileGenes$g[,1]))
    } else {
      genes$g <- NULL
    }
    
  })
  
  observeEvent(input$reset_net, {
    
    genes$g <- NULL
    
  })
  
 
  #TODO test also with csv input
  ##Convert to HGCN
  #look at switch function
  #avoid if-else nested
  #input_genes is a list of character vectors with gene names
  converted_genes <- reactive({
    if (input$organism != "Human") {
      if (input$organism == "Rat") {
        homologs <- rat_homologs
      } else if (input$organism == "Mouse") {
        homologs <- mouse_homologs
      } else {
        stop("USER Unsupported organism")
      }
    }
    
    tryCatch({
      conv_genes <-
        lapply(genes$g, function(x) {
          unique(homologs[which(homologs[, 2] %in% x), 3])
        })
      zero_id <- names(which(unlist(lapply(
        conv_genes, length
      )) == 0))
      if (length(zero_id) > 0) {
        shiny:::reactiveStop(paste("USER No gene IDs could be identified in:", zero_id))
      }
    }, error = function(e) {
      shiny:::reactiveStop("SYSTEM input genes conversion")
    })
    
    return(list(manual_input=conv_genes))
    
  })
  
  # when reading semicolon separated files,
  # having a comma separator causes `read.csv` to error
  observeEvent(input$file1, {
    
    fileGenes$g <- read.csv(input$file1$datapath, stringsAsFactors = F)
    
    
  })
  
  ###################################
  ## Dataset selection based on input
  #####################################
  observeEvent(input$go_heat|input$go_net|input$reset_net, {
    user_data$selection<-list(Sample=NA, Dataset=NA, Age=c(input$age_range[1],input$age_range[2]), Gender=input$input1.gender, Race=input$input1.race, 
                    Suicide=input$input1.suicide, Diagnosis=input$input1.diagnosis, Alcohol=NA, 
                    Smoke=NA, 
                    Platform=NA, Region_simpl=input$input1.brain_region)
    user_data$selection[which(user_data$selection=="Any")]<-NA
    
  })
  
  #Output: vector of samples' indexes
  #togliere il primo tryCatch perché ci sono dei valori di default
  #nella lista usare il nome del fattore invece che i
  #nel for fare direttamente l'intersezione degli indici
  samples_select <- reactive({
    #filter for factors not NA
    tryCatch({
      if (sum(is.na(user_data$selection)) == length(user_data$selection)) {
        #ind_sel_allfact <- c(1:nrow(metadata))
        ind_sel_allfact <- c(1:10)
        return(ind_sel_allfact)
      } else {
        sel_filter <- user_data$selection[-which(is.na(user_data$selection))]
        
        ind_sel <- list()
        i <- 0
        for (fact in names(sel_filter)) {
          i <- i + 1
          if (fact == "Age") {
            #for age, all values comprised in the range should be considered
            ind_sel[[i]] <-
              intersect(which(as.numeric(metadata[, fact]) >= sel_filter[[fact]][1]),
                        which(as.numeric(metadata[, fact]) <= sel_filter[[fact]][2]))
          } else {
            ind_sel[[i]] <- which(metadata[, fact] %in% sel_filter[[fact]])
          }
        }
        names(ind_sel) <- names(sel_filter)
        ind_sel_allfact <- Reduce(intersect, ind_sel)
        return(ind_sel_allfact)
      }
    }, error = function(e) {
      shiny:::reactiveStop(paste("SYSTEM samples selection -", e))
    })
})
  
 
  ##1. selecting the datasets matching all features and age range
  #Input: user selection of metadata conditions and #thr->minimum number of samples per dataset to consider it
  #Output: list of (vector of datasets) and (samples)
  dataset_select <- reactive({
    ind_sel_allfact <- samples_select()
    datasets <-
      names(which(table(metadata$Dataset[ind_sel_allfact]) >= input$thr))#select only datasets with at least thr samples with the selected features
    #if (length(datasets) == 0) {
    #  stop("USER No datasets with the selected features")
    #} else {
      ind_sel_allfact <- intersect(
        ind_sel_allfact, 
        which(metadata$Dataset %in% datasets)
      )#select only samples belonging to the chosen datasets
      samples_sel <- metadata$Sample[ind_sel_allfact] #get sample names
      return(list(datasets = datasets, samples = samples_sel))
 #   }
  })
  
  genes_forGSEA<-reactive({
  if(length(unlist(converted_genes()))>0){
    g<-converted_genes()
  } else {
    g <-
      lapply(test_list, function(x) {
        unique(rat_homologs[which(rat_homologs[, 2] %in% x), 3])
      })
    #g<-list(test_list=g)
  }
    return(g)
  })
  ####################################
  ### Compute GSEA
  ######################################
  
  forgesea<-function(dataset){
    corrs <- corr_age(dataset)
    
    forgsea <- unlist(corrs)
    names(forgsea) <- rownames(corrs)
    forgsea <- forgsea[!is.na(forgsea)]
    return(forgsea)
  }
  
   #TODO do I want to use this? YES
  corr_age <- function(dataset){
    #data <- get(dataset_select()[[1]][1])
    data<-get(dataset)
    
    #TODO if no genes in common return message
    samples <- intersect(dataset_select()[[2]], colnames(data))
    data <- data[, samples]
    
    age <- metadata$Age[match(samples, metadata$Sample)]
    #TODO add correlation type and possibility to choose "pairwise.complete.obs"
    corrs <- cor(t(data), age, use="pairwise.complete.obs")
    return(corrs)
  }
  
  
 #dat per ripetere su tutti i dataset
 GSEA_one_list <- function(dataset){
   #function(gene_list, dataset, forgsea, nPerm = 1000) {
   out_list<-list()
   tryCatch({
     fgsea_res <-
       fgsea(list(my_list = genes_forGSEA()), forgesea(dataset), nPermSimple = input$num_perm)
     fgsea_res$dataset<-dataset
     fgsea_res<-fgsea_res[,c("pathway","dataset", "pval", "padj", "log2err", "ES", "NES", "size", "leadingEdge")]
     
     p1 <- plotEnrichment(genes_forGSEA(), forgesea(dataset))
     df <- data.frame(p1$data, dataset = rep(dataset, nrow(p1$data)))
     
     out_list<-list(plot = df, GSEA = fgsea_res)
   }, error = function(e) {
     warning(paste("SYSTEM one GSEA failed -", e))
   })
   
     return(out_list)
   }
 
 
  GSEA_all <- observeEvent(input$go_heat|input$go_net|input$reset_net, {
    tryCatch({
      GSEA_dat <- list()
      GSEA_stats <- list()
      i <- 0
      failed <- c()
      for (dat in dataset_select()[[1]]) {
        GSEA_tmp <- GSEA_one_list(dat)
        if (length(GSEA_tmp) == 0) {
          failed <- c(failed, dat)
        } else {
          i <- i + 1
          GSEA_dat[[i]] <- GSEA_tmp[[1]]
          GSEA_stats[[i]] <- GSEA_tmp[[2]]
        }
      }
      
      if (i > 0) {
        names(GSEA_dat) <- setdiff(dataset_select()[[1]], failed)
        names(GSEA_stats) <- setdiff(dataset_select()[[1]], failed)
        
        df_tot <- data.frame(GSEA_dat[[1]])
        stats_tot<-data.frame(GSEA_stats[[1]])
        if (length(GSEA_dat) > 1) {
          for (n in 2:length(GSEA_dat)) {
            df_tot <- rbind.data.frame(df_tot, GSEA_dat[[n]])
            stats_tot <- rbind.data.frame(stats_tot, GSEA_stats[[n]])
          }
        }
        
        colnames(df_tot) <-
          c("rank", "enrichment", colnames(df_tot)[-c(1:2)])
        colnames(stats_tot)<-colnames(GSEA_stats[[1]])
        
        out$p <-
          ggplot(df_tot, aes(rank, enrichment)) + geom_line() + geom_hline(yintercept = 0, size =
                                                                                                 0.5) +
          facet_wrap(~ dataset, scale = "free_x") +
          theme(panel.background = element_blank())
        
        out$stats<-stats_tot
        
      } else {
        shiny:::reactiveStop("SYSTEM GSEA failed for all datasets and gene lists")
      }
    }, error = function(e) {
      shiny:::reactiveStop(paste("SYSTEM GSEA failed for all datasets and gene lists -", e))
    })
    
  })
  
  #######################################################
  ### merging of pvalues
  #######################################################
  
  ###plot summary of pvalues for each dataset and list
  pval_summary <- observeEvent(input$go_heat|input$go_net|input$reset_net, {
    #<- function(all_GSEA, input_genes, out_path) {
    tryCatch({
      GSEA_stats <- list()
      padj <- list()
      NES <- list()
      
        if (nrow(out$stats)==0) {
          shiny:::reactiveStop("SYSTEM no enrichment for any GSEA")
        } else {
      
      df <- data.frame(
        NES = unlist(out$stats[,"NES"]),
        padj = unlist(out$stats[,"padj"]),
        gene_list = unlist(out$stats[,"pathway"]),
        dataset = unlist(out$stats[,"dataset"]))
      
      if (sum(is.na(df$NES)) == length(df$NES) |
          sum(is.na(df$padj)) == length(df$padj)) {
        shiny:::reactiveStop("USER/SYSTEM No p-value for any GSEA: try increasing the num_perm parameter")
      } else {
        p <-
          ggplot(df, aes(
            x = NES,
            y = -log10(padj),
            label = dataset
          )) + geom_point() +
          geom_label_repel() +
          geom_vline(xintercept = 0,
                     linetype = "dashed",
                     color = "red") + theme_classic() + facet_wrap(~ gene_list)
        
        out$p_summary<-p
         }
    }
      }, error = function(e) {
        shiny:::reactiveStop(e)
    })
  })
  
  #p<-merged_p(all, input_genes)
  
  meta_one_list_reversal <- reactive({
    #function(GSEA_stats, direction = "reversal") {
    tryCatch({
      datasets <- names(out$stats)
      istwo <- rep(T, length(out$stats[,"dataset"]))
      
        toinvert <-
          ifelse(sign(out$stats[,"NES"]) == (1), T, F)
        p_all <- out$stats[,"pval"]
      
      
      missing <- union(which(is.na(toinvert)), which(is.na(p_all)))
      if (length(missing) == 0) {
        p <-
          sumlog(two2one(out$stats[,"pval"], two = istwo, invert = toinvert))
        return(p[[3]])
      } else if (length(missing) == length(toinvert)) {
        shiny:::reactiveStop(
          "USER/SYSTEM No non-missing p-values for one list of genes: try increasing num_perm"
        )
      } else if (length(missing) == (length(toinvert) - 1)) {
        p <-
          two2one(out$stats[,"pval"][-missing], two = istwo[-missing], invert = toinvert[-missing])
      } else {
        p <-
          sumlog(two2one(out$stats[,"pval"][-missing], two = istwo[-missing], invert = toinvert[-missing]))
        return(p[[3]])
      }
    }, error = function(e) {
      shiny:::reactiveStop(e)
    })
  })
  
  
  meta_one_list_mimick <- reactive({
    #function(GSEA_stats, direction = "reversal") {
    tryCatch({
      datasets <- names(out$stats)
      istwo <- rep(T, length(out$stats[,"dataset"]))
      
        toinvert <-
          ifelse(sign(out$stats[,"NES"]) == (-1), T, F)
        p_all <- out$stats[,"pval"]
      
      
      missing <- union(which(is.na(toinvert)), which(is.na(p_all)))
      if (length(missing) == 0) {
        p <-
          sumlog(two2one(out$stats[,"pval"], two = istwo, invert = toinvert))
        return(p[[3]])
      } else if (length(missing) == length(toinvert)) {
        shiny:::reactiveStop(
          "USER/SYSTEM No non-missing p-values for one list of genes: try increasing num_perm"
        )
      } else if (length(missing) == (length(toinvert) - 1)) {
        p <-
          two2one(out$stats[,"pval"][-missing], two = istwo[-missing], invert = toinvert[-missing])
      } else {
        p <-
          sumlog(two2one(out$stats[,"pval"][-missing], two = istwo[-missing], invert = toinvert[-missing]))
        return(p[[3]])
      }
    }, error = function(e) {
      shiny:::reactiveStop(e)
    })
  })
  
  meta_ratio<-reactive({
    p_ratio<-meta_one_list_reversal()/meta_one_list_mimick()
    return(p_ratio)
  })
  ################################
  ### Output
  #################################
  
  observeEvent(input$go_heat|input$go_net|input$reset_net, {
    output$message<-renderText({class(out$stats)})
    output$table <- renderDataTable({out$stats[,c("dataset", "pval", "padj", "NES")]})
    output$GSEA<-renderPlot({out$p})
    output$summary<-renderPlot({out$p_summary})
    output$p_ratio<-renderText({meta_ratio()})
  })
  
}

# Run app ----
shinyApp(ui, server)