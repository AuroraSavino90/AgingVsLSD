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

load("Data/Alldata_20Sep_Age.RData", envir=.GlobalEnv)
load("Data/DE_GSE179379.RData", envir=.GlobalEnv)
rat_homologs<-read.csv("Data/Human rat homologs.txt")
mouse_homologs<-read.csv("Data/Human mouse homologs.txt")
test_list<-rownames(DE_GSE179379)[which(DE_GSE179379$log2FoldChange>0 & DE_GSE179379$padj<0.05)]



ui <- navbarPage( 
  title = "Age reversal",
  
  tags$style(type="text/css",
             ".shiny-output-error { visibility: hidden; }",
             ".shiny-output-error:before { visibility: hidden; }"
  ),
  
  ###### Here : insert shinydashboard dependencies ######
  header = tagList(
    useShinydashboard()
  ),
  
  tabPanel(
    
    #titlePanel("MetaLCM breast"),
    introBox("Test reversal",
             data.step = 1,
             data.intro = "With this tool, you can test the age-reversal impact of an input gene list"),
    actionButton("tour","Start Tour", class = "btn-warning btn-lg"),
    sidebarLayout(
      
      #introjsUI(),
      
      sidebarPanel(width=3,
                   fluidRow(
                     
                     
                     
                     column(12,
                            introBox(
                              
                              data.position = "top",
                              data.step = 2,
                              data.intro = "To start the analysis, select the metadata to filter the samples of human brain aging",
                              h3("Metadata filters"),   
                              
                              sliderInput("age_range", "Age range:",
                                          min = 0, max = max(as.numeric(metadata$Age), na.rm=T),
                                          value = c(20,106)),
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
                                                   selected = c("Any")),
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
                                                   selected = "Healthy"),
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
                                          value = 10)),
                            introBox(
                              data.step = 3,
                              data.intro = "Here, you can set the number of permutations to compute the GSEA p-value",
                              
                            
                              radioButtons("num_perm", "Number of permutations for GSEA p-value",
                                           choices = c("100"= 100, "1000" = 1000, "10000" = 10000, "100000" = 100000, "1000000" = 1000000 ),
                                          selected = 1000)
                              
                            ))
                    
                     
                   ),
                   
                   hr(), 
                   introBox(
                     data.step = 4,
                     data.intro = "To test the age-reversal potential of a gene list you can either
                 type gene symbols separated by commas
                 or import a csv file containing gene symbols in one single column",
                     
                     h3("Genes' selection"),
                     #to choose a set of genes for filtering
                     fluidRow(
                       column(width=6, radioButtons("organism", h4("Organism"),
                                           choices = list("Rat" = "Rat", "Mouse" = "Mouse",
                                                          "Human" = "Human"),selected = "Rat")),
                       column(width=6, radioButtons("direction", h4("Direction"),
                                           choices = list("Up-regulated" = "Up", "Down-regulated" = "Down"),selected = "Up"))
                       ),
                     
                     
                    
                     textInput("Genes_net", label = strong("Enter gene symbols"), value = ""),
                     #column(style="padding-top:25px;", 3, actionButton("go_net", "Submit gene list"))),
                     
                     p("OR", align="center"),
                     # Input: Select a file ----
                     uiOutput('file1_ui')
                     ),
                   
                   introBox(
                     
                     data.step = 5,
                     data.intro = "Press \"Submit gene list\" to import gene symbols. Press \"Run test list\" to run the analysis using the example gene list.
                     Press one of these buttons each time you change the settings.",
                     
                     actionButton("go_list", "Submit gene list", style="color: #fff; background-color: #3cb043; border-color: #2e6da4"),
                     actionButton("cancel_list", "Delete uploaded list", style="color: #fff; background-color: #3cb043; border-color: #2e6da4"),
                     actionButton("reset_list", "Run test list", style="color: #fff; background-color: #3cb043; border-color: #2e6da4")
                   )),
      
      
      mainPanel(
        add_busy_spinner(spin = "double-bounce", position="full-page"),
        
        tabsetPanel(type = "tabs",
                    
                    tabPanel(
                    introjsUI(),
                    
                    introBox(
                      data.step = 6,
                      data.intro = "Here, you will find the results of the analysis: an overall age-reversal score (between -1 and 1), plots for GSEAs run on the input list and each of the datasets with the selected features,
                      a summary of GSEA p-values and normalized enrichment scores (NES), a table with GSEA statistics for each dataset with the selected features.
                      Once the analysis is completed and plots/tables are produced, you can download the results via the download buttons that will appear.",
      
                    
                    textOutput("error_message"),
                    
                    fluidRow(
                    
                    textOutput("message"),
                    textOutput("message2"),
                    
                    column(12, align="center",
                           valueBox(width = 12, value=textOutput("p_ratio"), subtitle="Age reversal Score")
                    )
                     ),
                    
                    fluidRow(
                    #box(title="Plot of GSEA", )),
                      uiOutput("buttonGSEA"),
                      plotOutput("GSEA",height=420, width=600)),
                      
                    fluidRow(
                    #box(title="Summary plot of GSEA stats", ),
                      uiOutput("buttonSummary"),
                      
                      plotOutput("summary",height=420, width=600),
                    #box(title="Summary table of GSEA stats", )
                    uiOutput("buttonTable"),
                    dataTableOutput("table")
                    
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
  user_data<-reactiveValues(selection=NULL, threshold=NULL, num_perm=NULL)
  out<-reactiveValues(p=NULL, stats=NULL, p_ratio=NULL, meta_one_list_reversal=1, meta_one_list_mimick=1)
  
  ###lista di geni in input
  observeEvent(input$go_list, {
    
    if(length(unlist(strsplit(input$Genes_net, "[, ]+")))>0){
      genes$g <- strsplit(input$Genes_net, "[, ]+")
    } else if (!is.null(fileGenes$g)){
      genes$g<-(c(fileGenes$g[,1]))
    } else {
      genes$g <- NULL
    }
    
  })
  
  observeEvent(input$reset_list, {
    
    genes$g <- test_list
    
  })
  
    converted_genes <- eventReactive(input$go_list|input$reset_list, {
    out$error_message0<-NULL
    tryCatch({
    if (input$organism != "Human") {
      if (input$organism == "Rat") {
        homologs <- rat_homologs
      } else if (input$organism == "Mouse") {
        homologs <- mouse_homologs
      } else {
        stop("USER Unsupported organism")
      }
    
    
    
      conv_genes <-
        lapply(genes$g, function(x) {
          unique(homologs[which(homologs[, 2] %in% x), 3])
        })
      zero_id <- names(which(unlist(lapply(
        conv_genes, length
      )) == 0))
      if (length(zero_id) > 0) {
        out$error_message0<-paste("USER No gene IDs could be identified in:", zero_id)
        shiny:::reactiveStop(paste("USER No gene IDs could be identified in:", zero_id))
      }
    } else {
      conv_genes<-genes$g
    }
      }, error = function(e) {
        conv_genes<-NULL
      out$error_message0<-"SYSTEM input genes conversion"
      shiny:::reactiveStop("SYSTEM input genes conversion")
      
    },
    finally={
      return(list(manual_input=conv_genes))
    })
  })
  
  # when reading semicolon separated files,
  # having a comma separator causes `read.csv` to error
  observeEvent(input$file1, {
    
    fileGenes$g <- read.csv(input$file1$datapath, stringsAsFactors = F)
    
  })
  
  observeEvent(input$cancel_list, {
    
    fileGenes$g <- NULL
    
  })
  
  output$file1_ui <- renderUI({
    input$cancel_list ## Create a dependency with the reset button
    fileInput("file1", "Choose CSV File",
              multiple = FALSE,
              accept = c(".csv"))
  })
  
  ###################################
  ## Dataset selection based on input
  #####################################
  observeEvent(input$go_list|input$reset_list, {
    user_data$selection<-list(Sample=NA, Dataset=NA, Age=c(input$age_range[1],input$age_range[2]), Gender=input$input1.gender, Race=input$input1.race, 
                    Suicide=input$input1.suicide, Diagnosis=input$input1.diagnosis, Alcohol=NA, 
                    Smoke=NA, 
                    Platform=NA, Region_simpl=input$input1.brain_region)
    user_data$selection[which(user_data$selection=="Any")]<-NA
    user_data$threshold<-input$thr
    user_data$num_perm<-as.numeric(input$num_perm)
    
  })
  
  #Output: vector of samples' indexes
  #togliere il primo tryCatch perché ci sono dei valori di default
  #nella lista usare il nome del fattore invece che i
  #nel for fare direttamente l'intersezione degli indici
  samples_select <- reactive({
    out$error_message1<-NULL
    #filter for factors not NA
    tryCatch({
      if (sum(is.na(user_data$selection)) == length(user_data$selection)) {
        ind_sel_allfact <- c(1:nrow(metadata))
        
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
        
      }
    }, error = function(e) {
      shiny:::reactiveStop(paste("SYSTEM samples selection -", e))
      out$error_message1<-paste("SYSTEM samples selection -", e)
      ind_sel_allfact <- NULL
    }, finally={
      return(ind_sel_allfact)
    })
})
  
 
  ##1. selecting the datasets matching all features and age range
  #Input: user selection of metadata conditions and #thr->minimum number of samples per dataset to consider it
  #Output: list of (vector of datasets) and (samples)
  dataset_select <- reactive({
    ind_sel_allfact <- samples_select()
    datasets <-
      names(which(table(metadata$Dataset[ind_sel_allfact]) >= user_data$threshold))#select only datasets with at least thr samples with the selected features
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
    g<-unlist(converted_genes())
  } else {
    g <- NULL
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
   
   out_list<-list()
   tryCatch({
     commgenes<-intersect(unlist(genes_forGSEA()), names(forgesea(dataset)))
     if(length(commgenes)>0){
     fgsea_res <-
       fgsea(list(my_list = genes_forGSEA()), forgesea(dataset), nPermSimple = user_data$num_perm)
     fgsea_res$dataset<-dataset
     fgsea_res<-fgsea_res[,c("pathway","dataset", "pval", "padj", "log2err", "ES", "NES", "size", "leadingEdge")]
     
     p1 <- plotEnrichment(genes_forGSEA(), forgesea(dataset))
     df <- data.frame(p1$data, dataset = rep(dataset, nrow(p1$data)))
     
     out_list<-list(plot = df, GSEA = fgsea_res, ngenes=length(commgenes))
     } else {
       out_list<-list(plot = NULL, GSEA = NULL, ngenes=length(commgenes))
     }
   }, error = function(e) {
     warning(paste("SYSTEM one GSEA failed -", e))
   })
   
     return(out_list)
   }
 
  GSEA_all <- reactive({
    out$error_message2<-NULL
    if(length(converted_genes())>0){
    tryCatch({
      GSEA_dat <- list()
      GSEA_stats <- list()
      i <- 0
      failed <- c()
      for (dat in dataset_select()[[1]]) {
        GSEA_tmp <- GSEA_one_list(dat)
        if (length(GSEA_tmp) == 0) {
          failed <- c(failed, dat)
        } else if (GSEA_tmp[[3]] == 0){
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
        
        df_tot$dataset<-factor(df_tot$dataset, levels=unique(c(stats_tot$dataset[order(stats_tot$NES)])))
        
        p <-
          ggplot(df_tot, aes(rank, enrichment)) + geom_line() + geom_hline(yintercept = 0, size =
                                                                                                 0.5) +
          facet_wrap(~ dataset, scale = "free_x") +
          theme(panel.background = element_blank())
        
        
        return(list(p=p, stats=stats_tot))
        
      } else {
        out$error_message2<-"SYSTEM GSEA failed for all datasets and gene lists"
        return(list(p=NULL, stats=NULL))
        shiny:::reactiveStop("SYSTEM GSEA failed for all datasets and gene lists")
      }
    }, error = function(e) {
      out$error_message2<-paste("SYSTEM GSEA failed for all datasets and gene lists -", e)
      return(list(p=NULL, stats=NULL))
      shiny:::reactiveStop(paste("SYSTEM GSEA failed for all datasets and gene lists -", e))
      
    })
    }
  })
  
  #######################################################
  ### merging of pvalues
  #######################################################
  
  ###plot summary of pvalues for each dataset and list
  pval_summary <- reactive({
    out$error_message3<-NULL
    tryCatch({
      GSEA_stats <- list()
      padj <- list()
      NES <- list()
      
        if (is.null(GSEA_all()[["stats"]]) | nrow(GSEA_all()[["stats"]])==0) {
          out$error_message3<-"SYSTEM no enrichment for any GSEA"
          p_summary<-NULL
          return(p_summary)
          shiny:::reactiveStop("SYSTEM no enrichment for any GSEA")
        } else {
      
      df <- data.frame(
        NES = unlist(GSEA_all()[["stats"]][,"NES"]),
        padj = unlist(GSEA_all()[["stats"]][,"padj"]),
        gene_list = unlist(GSEA_all()[["stats"]][,"pathway"]),
        dataset = unlist(GSEA_all()[["stats"]][,"dataset"]))
      
      if (sum(is.na(df$NES)) == length(df$NES) |
          sum(is.na(df$padj)) == length(df$padj)) {
        out$error_message3<-"USER/SYSTEM No p-value for any GSEA: try increasing the num_perm parameter"
        p_summary<-NULL
        return(p_summary)
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
        
        p_summary<-p
        return(p_summary)
         }
    }
      }, error = function(e) {
        out$error_message3<-e
        stop(e)
        p_summary<-NULL
        return(p_summary)
    })
  })
  
  #p<-merged_p(all, input_genes)
  
  meta_one_list_reversal <- reactive({
    out$error_message4<-NULL
    if(length(unlist(converted_genes()))>0 & length(dataset_select()[[1]])>1){
      if(!is.null(GSEA_all()[["stats"]])){
    #function(GSEA_stats, direction = "reversal") {
    tryCatch({
      
      datasets <- names(GSEA_all()[["stats"]])
      istwo <- rep(T, length(GSEA_all()[["stats"]][,"dataset"]))
      
        toinvert <-
          ifelse(sign(GSEA_all()[["stats"]][,"NES"]) == (1), T, F)
        p_all <- GSEA_all()[["stats"]][,"pval"]
      
      
      missing <- union(which(is.na(toinvert)), which(is.na(p_all)))
      if (length(missing) == 0) {
        p <-
          sumlog(two2one(GSEA_all()[["stats"]][,"pval"], two = istwo, invert = toinvert))
        if(p[[3]]==0){
          p[[3]]<-10^(-320)
        }
        return(-log10(p[[3]]))
      } else if (length(missing) == length(toinvert)) {
        out$error_message4<-"USER/SYSTEM No non-missing p-values for one list of genes: try increasing num_perm"
        shiny:reactiveStop(
          "USER/SYSTEM No non-missing p-values for one list of genes: try increasing num_perm"
        )
        return(NULL)
      } else if (length(missing) == (length(toinvert) - 1)) {
        p <-
          two2one(GSEA_all()[["stats"]][,"pval"][-missing], two = istwo[-missing], invert = toinvert[-missing])
        if(p==0){
          p<-10^(-320)
        }
        return(-log10(p))
      } else {
        p <-
          sumlog(two2one(GSEA_all()[["stats"]][,"pval"][-missing], two = istwo[-missing], invert = toinvert[-missing]))
        if(p[[3]]==0){
          p[[3]]<-10^(-320)
        }
        return(-log10(p[[3]]))
      }
    }, error = function(e) {
      out$error_message4<-e
      stop(e)
    })
    } else if (length(unlist(converted_genes()))>0 & length(dataset_select()[[1]])==1){
      
      tryCatch({
        datasets <- names(GSEA_all()[["stats"]])
        istwo <- rep(T, length(GSEA_all()[["stats"]][,"dataset"]))
        
        toinvert <-
          ifelse(sign(GSEA_all()[["stats"]][,"NES"]) == (1), T, F)
        p_all <- GSEA_all()[["stats"]][,"pval"]
        
        
        missing <- union(which(is.na(toinvert)), which(is.na(p_all)))
        if (length(missing) == 0) {
          p <-two2one(GSEA_all()[["stats"]][,"pval"], two = istwo, invert = toinvert)
          if(p==0){
            p<-10^(-320)
          }
          return(-log10(p))
        } else if (length(missing) == length(toinvert)) {
          out$error_message4<-"USER/SYSTEM No non-missing p-values for one list of genes: try increasing num_perm"
          shiny:reactiveStop(
            "USER/SYSTEM No non-missing p-values for one list of genes: try increasing num_perm"
          )
          return(NULL)
        } 
      }, error = function(e) {
        out$error_message4<-e
        stop(e)
      })
      
    }
    } else {
      return(NULL)
    }
  })
  
  
  meta_one_list_mimick <- reactive({
    out$error_message5<-NULL
    if(!is.null(GSEA_all()[["stats"]])){
    if(length(unlist(converted_genes()))>0 & length(dataset_select()[[1]])>1){
      
    tryCatch({
      datasets <- names(GSEA_all()[["stats"]])
      istwo <- rep(T, length(GSEA_all()[["stats"]][,"dataset"]))
      
        toinvert <-
          ifelse(sign(GSEA_all()[["stats"]][,"NES"]) == (-1), T, F)
        p_all <- GSEA_all()[["stats"]][,"pval"]
      
      
      missing <- union(which(is.na(toinvert)), which(is.na(p_all)))
      if (length(missing) == 0) {
        p <-
          sumlog(two2one(GSEA_all()[["stats"]][,"pval"], two = istwo, invert = toinvert))
        if(p[[3]]==0){
          p[[3]]<-10^(-320)
        }
        return(-log10(p[[3]]))
      } else if (length(missing) == length(toinvert)) {
        out$error_message5<-"USER/SYSTEM No non-missing p-values for one list of genes: try increasing num_perm"
        shiny:::reactiveStop(
          "USER/SYSTEM No non-missing p-values for one list of genes: try increasing num_perm"
        )
        return(NULL)
      } else if (length(missing) == (length(toinvert) - 1)) {
        p <-
          two2one(GSEA_all()[["stats"]][,"pval"][-missing], two = istwo[-missing], invert = toinvert[-missing])
        if(p==0){
          p<-10^(-320)
        }
        return(-log10(p))
        } else {
        p <-
          sumlog(two2one(GSEA_all()[["stats"]][,"pval"][-missing], two = istwo[-missing], invert = toinvert[-missing]))
        if(p[[3]]==0){
          p[[3]]<-10^(-320)
        }
        return(-log10(p[[3]]))
      }
    }, error = function(e) {
      out$error_message5<-e
      stop(e)
    })
    } else if (length(unlist(converted_genes()))>0 & length(dataset_select()[[1]])==1){
      
      tryCatch({
        datasets <- names(GSEA_all()[["stats"]])
        istwo <- rep(T, length(GSEA_all()[["stats"]][,"dataset"]))
        
        toinvert <-
          ifelse(sign(GSEA_all()[["stats"]][,"NES"]) == (-1), T, F)
        p_all <- GSEA_all()[["stats"]][,"pval"]
        
        
        missing <- union(which(is.na(toinvert)), which(is.na(p_all)))
        if (length(missing) == 0) {
          p <-two2one(GSEA_all()[["stats"]][,"pval"], two = istwo, invert = toinvert)
          if(p==0){
            p<-10^(-320)
          }
          return(-log10(p))
        } else if (length(missing) == length(toinvert)) {
          out$error_message5<-"USER/SYSTEM No non-missing p-values for one list of genes: try increasing num_perm"
          shiny:reactiveStop(
            "USER/SYSTEM No non-missing p-values for one list of genes: try increasing num_perm"
          )
          return(NULL)
        } 
      }, error = function(e) {
        out$error_message5<-e
        stop(e)
      })
    }
    } else {
      return(NULL)
    }
  })
  
  meta_ratio<-observeEvent(input$go_list|input$reset_list, {
    out$p_ratio<-NULL
    
    if(!is.null(meta_one_list_reversal()) & !is.null(meta_one_list_mimick())){
    
    if(input$direction=="Up"){
    p_ratio<-meta_one_list_reversal()-meta_one_list_mimick()
    out$p_ratio<- p_ratio/320
    } else if(input$direction=="Down"){
      p_ratio<-meta_one_list_mimick()-meta_one_list_reversal()
      out$p_ratio<- p_ratio/320
    }
    } else {
      return(NULL)
    }
      
  })
  
  ################################
  ### Output
  #################################
  
  observeEvent(input$reset_list, {
    output$message<-NULL
    if(length(dataset_select()[[1]])>=1){
      output$message2<-NULL
      output$error_message<-renderText({paste(out$error_message0,
                                              out$error_message1,
                                              out$error_message2,
                                              out$error_message3,
                                              out$error_message4,
                                              out$error_message5, sep="\t")})
      output$table <- renderDataTable({GSEA_all()[["stats"]][,c("dataset", "padj", "NES")]})
      output$GSEA<-renderPlot({GSEA_all()[["p"]]})
      output$summary<-renderPlot({pval_summary()})
      output$p_ratio<-renderText({out$p_ratio})
    } else if(length(dataset_select()[[1]])==0){
      output$message<-renderText({"No dataset with the selected features"})
      output$message2<-NULL
      output$error_message<-NULL
      output$table <- NULL
      output$GSEA<-NULL
      output$summary<-NULL
      output$p_ratio<-NULL
    }
   
  })
  
  observeEvent(input$go_list, {
    if(length(unlist(converted_genes()))==0){
      output$message<-renderText({"No gene list provided or gene IDs not matching the required input type"})
      output$message2<-renderText({unlist(fileGenes$g)})
      output$error_message<-NULL
     } else {
       output$message<-NULL
       if(length(dataset_select()[[1]])>=1){
         output$message2<-renderText({unlist(fileGenes$g)})
         output$error_message<-renderText({paste(out$error_message0,
                                                 out$error_message1,
                                                 out$error_message2,
                                                 out$error_message3,
                                                 out$error_message4,
                                                 out$error_message5, sep="\t")})
         output$table <- renderDataTable({GSEA_all()[["stats"]][,c("dataset", "padj", "NES")]})
         output$GSEA<-renderPlot({GSEA_all()[["p"]]})
         output$summary<-renderPlot({pval_summary()})
         output$p_ratio<-renderText({out$p_ratio})
       } else if(length(dataset_select()[[1]])==0){
         output$message2<-renderText({"No dataset with the selected features"})
         output$error_message<-NULL
         output$table <- NULL
         output$GSEA<-NULL
         output$summary<-NULL
         output$p_ratio<-NULL
       }
    }
     })
  
  

  output$downloadGSEA <- downloadHandler(
     
     filename = function() {
       paste('GSEA-', Sys.Date(), '.pdf', sep='')
     },
     content = function(file) {
       if(!is.null(GSEA_all()[["p"]])){
       outplotGSEA<-GSEA_all()[["p"]]
       pdf(file)
       grid::grid.draw(outplotGSEA)
       dev.off()
         }
       }
   )



output$downloadSummary <- downloadHandler(
  
  filename = function() {
    paste('GSEA_summary-', Sys.Date(), '.pdf', sep='')
  },
  content = function(file) {
    if(!is.null(pval_summary())){
      outplotsummary<-pval_summary()
      pdf(file)
      grid::grid.draw(outplotsummary)
      dev.off()
    }
  }
)


output$downloadTable <- downloadHandler(
  
  filename = function() {
    paste('GSEA_table-', Sys.Date(), '.csv', sep='')
  },
  content = function(file) {
    if(!is.null(GSEA_all()[["stats"]])){
      outtable<-data.frame(GSEA_all()[["stats"]][,c("dataset", "padj", "NES")])
      write.csv(outtable, file, row.names=T, quote=F)
    }
  }
)

output$buttonGSEA<-renderUI({
  if(!is.null(GSEA_all()[["p"]])){
    downloadButton('downloadGSEA', 'Save GSEA plot')
  }
})

output$buttonSummary<-renderUI({
  if(!is.null(pval_summary())){
    downloadButton('downloadSummary', 'Save summary plot')
  }
})

output$buttonTable<-renderUI({
  if(!is.null(GSEA_all()[["stats"]])){
    downloadButton('downloadTable', 'Save stats table')
  }
})




}
# Run app ----
shinyApp(ui, server)