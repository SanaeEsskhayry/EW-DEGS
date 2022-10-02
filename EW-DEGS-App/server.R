#********************************************************************************
#Download libraries:
#******************************************************************************** 

library(shiny)
library(rtracklayer)
library(shinydashboard)
library(Rqc)
library(QuasR)
library(Rsubread)
library(DT)
library(DESeq2)
source("DEGs.R")

options(shiny.maxRequestSize=1000*1024^2)

# Define server logic.
shinyServer(function(input, output) {
  
  isNA<-reactiveValues(x = NA)
  
  sys.name<-Sys.info()[['user']]
  
  is.finished<-reactiveValues(id=NULL)
  
  show.notification<-reactive({
    # If there's currently a notification, don't add another
    if (!is.null(is.finished$id))
      return()
    # Save the ID for removal later
    is.finished$id <- showNotification(paste("Please wait, this is take a few minutes of processing..."), duration = 0, type="message")
  })
  
  remove.notification<-reactive({
    if (!is.null(is.finished$id))
      removeNotification(is.finished$id)
    is.finished$id = NULL
  })
  
  
  #*****************************************************************************
  ##############################3Downloading Data:#############################
  #*****************************************************************************
  
  import.metadata<-reactive({
    metadata_file<-input$dw_pheno_file$datapath
    metadata <- read.table (metadata_file,sep=",",header=TRUE)
  })
  #-----------------------------------------------------------------------------
  #**1- Get information about the phenotype 
  #*----------------------------------------------------------------------------
  
  output$pheno_description<-DT:: renderDataTable({
    if(is.null(input$dw_pheno_file) ){return (as.matrix("Import the metadata of the study"))}
    DT::datatable(import.metadata(), options = list(scrollX = TRUE)) 
  })
  
  #-----------------------------------------------------------------------------
  #**2-Get metadata about the imported fastq files
  #*----------------------------------------------------------------------------

  output$select_gr_field<-renderUI({
    if(is.null(input$dw_pheno_file) ){return (tags$span(style="color:blue", "Import the metadata of the study !"))}
    selectInput("selectGroup","Select group Field?",
                choices=names(import.metadata()) )
  })
  #-----------------------------------------------------------------------------
  #**3-Downloading data from the folder
  #*----------------------------------------------------------------------------
  
  read.data<-reactive({
    pheno<-import.metadata()
    pheno <- cbind (pheno, group=pheno[[input$selectGroup]])
    group <- factor(pheno$group)
    
    checkpoint("qcResult", path=".", {
      folder<-input$folderName
      files <- list.files(full.names=TRUE, path=folder)
      qcResult <- rqcQA(files,group = group, workers=5)
    },keep="qcResult")
    qcResult
  })
  
   getData<-eventReactive(input$btn_upload,{
     show.notification()
     rqc.result<-read.data()
     perFileInformation(rqc.result)
    }, ignoreNULL = TRUE )
  
   output$file_description <-renderTable( {
    if(input$folderName==""){return ("Enter the directory path to your FASTQ files")}
    if(is.null(input$dw_pheno_file)){return ("Import the metadata of the study")}
    else{
      remove.notification()
      getData()
    } })
   #----------------------------------------------------------------------------
   #**4-perFileTopReads
   #*---------------------------------------------------------------------------
  
   getPerFileTopReads<-eventReactive(input$btn_perFileTopReads,{
    DT::datatable(as.data.frame(perFileTopReads(read.data())), options = list(scrollX = TRUE))
  }, ignoreNULL = TRUE )
   output$perFileTopReads<-DT::renderDataTable({
     if(input$folderName==""){return (as.matrix("Enter the directory path to your FASTQ files"))}
     if(is.null(input$dw_pheno_file)){return (as.matrix("Import the metadata of the study"))}
     else
     getPerFileTopReads()
    })
  
  
  #*******************************************************************************************************************
  ######################################################## Plots of the quality  #####################################
  #*******************************************************************************************************************
  output$select_fastq_file<-renderUI({
    if(input$folderName=="" ||is.null(input$dw_pheno_file) ){return ("Verify if you are imported your data") }
    radioButtons("selected_file","Select a fastq file",choices =c("None", 
                  m <- c(),for(i in 1:length(read.data())) { y <- as.character(perFileInformation(subsetByPair(read.data(), i))$filename)
                  m <- c(m, y)},
                  m )) }) 
  
  output$quality_plot <-renderPlot({

    if(input$folderName==""  || input$selected_file=="None"){
      return (plot(1, type = "n", xlab = "",ylab = "", xlim = c(0, 5), ylim = c(0, 5)))
    }
    selectedPlot = input$select_plot
    selectedFile = input$selected_file
    result = switch(  
      selectedPlot,  
      "Per Read Mean Quality Distribution of Files"= rqcReadQualityBoxPlot(read.data()[selectedFile]),  
      "Average Quality"= rqcReadQualityPlot(read.data()[selectedFile]),  
      "Cycle-specific Average Quality"= rqcCycleAverageQualityPlot(read.data()[selectedFile]),  
      "Read Frequency"= rqcReadFrequencyPlot(read.data()[selectedFile]) ,
      "Read Length Distribution"= rqcReadWidthPlot(read.data()[selectedFile]),
      "Cycle-specific GC Content"= rqcCycleGCPlot(read.data()[selectedFile]),
      "Cycle-specific Quality Distribution"= rqcCycleQualityPlot(read.data()[selectedFile]),
      "Cycle-specific Quality Distribution - Boxplot"=rqcCycleQualityBoxPlot(read.data()[selectedFile]),
      "Cycle-specific Base Call Proportion"=rqcCycleBaseCallsLinePlot(read.data()[selectedFile])
    )
    result
  })
  
  output$quality_data <-DT::renderDataTable ({
    if(input$folderName=="" || input$selected_file=="None"){return (as.matrix("Must have at least one file selected.s"))}
    selectedPlot = input$select_plot
    selectedFile = input$selected_file
    result = switch(  
      selectedPlot,  
      "Per Read Mean Quality Distribution of Files"=  rqcReadQualityBoxCalc(read.data()[selectedFile]),
      "Average Quality"= rqcReadQualityCalc(read.data()[selectedFile]), 
      "Cycle-specific Average Quality"= rqcCycleAverageQualityCalc(read.data()[selectedFile]), 
      "Read Frequency"= rqcReadFrequencyCalc(read.data()[selectedFile]),
      "Read Length Distribution"= rqcReadWidthCalc(read.data()[selectedFile]), 
      "Cycle-specific GC Content"= rqcCycleGCCalc(read.data()[selectedFile]), 
      "Cycle-specific Quality Distribution"= rqcCycleQualityCalc(read.data()[selectedFile]),
      "Cycle-specific Quality Distribution - Boxplot"=rqcCycleQualityBoxCalc(read.data()[selectedFile]), 
      "Cycle-specific Base Call Proportion"= rqcCycleBaseCallsCalc(read.data()[selectedFile])
    )
    result
  })

  #----------------------------------------------------------------------------------------------------------------
  ######################################## Button of Download Quality Report:#####################################
  #----------------------------------------------------------------------------------------------------------------
  generate.report<-eventReactive(input$btn_QR,{
     show.notification()
     sys.name= Sys.info()[['user']]
     outdir = paste("C:/Users/",sys.name,"/Documents/",sep = "")
     rqcReport(read.data(),outdir = outdir ,file = paste( "Quality-report-",Sys.Date()))
  })
  output$qualityReport<-renderText({
    path<-as.character(generate.report())
    remove.notification()
    paste("Check the quality report in this path : \n ",path) 
  })
  
  
  #******************************************************************************************************************
  

  
  #----------------------------------------------------------------------------------------------------------------
  ######################################## Filter and Trimming :#####################################
  #----------------------------------------------------------------------------------------------------------------
  
  #****************************************************************************************************************
  #-----------------------------------------------------------------------------------------------------------------
  #Button of Filter and trimming:                                                                                   |
  #-----------------------------------------------------------------------------------------------------------------
  
  filter_Trim<-eventReactive(input$btn_filter,{
    if(input$folderName==""){return(as.matrix("Import your FASTQ file(s)"))}
    show.notification()
    folder<-input$folderName
    
    fastqFiles <- list.files(full.names=TRUE, path=folder)
    
    file.format=".fastq.gz"
    
    read.files<- stringr::str_remove(fastqFiles,".fastq.gz") 
    
    
    results<-preprocessReads(fastqFiles,paste(read.files,"-Filtered",file.format,sep=""),
                             nBases = as.numeric(input$nBases),
                             truncateStartBases = as.numeric(input$truncateStartBases),
                             truncateEndBases = as.numeric(input$truncateEndBases),
                             Lpattern =as.character(input$Lpattern),
                             Rpattern = as.character(input$Rpattern),
                             complexity = as.numeric(input$complexity),
                             minLength = as.numeric(input$minLength))
    
    DT::datatable(as.data.frame(results))
    })
  
    output$tb_filter_trim<-DT:: renderDataTable({
      remove.notification()
      filter_Trim()
      })
  
  #*********************************************************************************************************************
  
  
  
  
  #*********************************************************************************************************************
  
  
  #---------------------------------------------------------------------------------------------------------------------
  ################################################ Read mapping or alignment  #########################################                                                                                      |
  #---------------------------------------------------------------------------------------------------------------------
    
  
   
    ############################################################################
    run.alignment<-reactive({
      #---------------------------------------------------------------------------
      #Before lunch alignment we have to build the index:                         |
      #---------------------------------------------------------------------------
      
      ref <- input$reference_file$datapath
      buildindex(basename="reference_index",reference=ref)
      
      folder<-input$filtered_data_folder
      files <-list.files(full.names=TRUE, path=folder)
      read.files<- stringr::str_remove(files,".fastq.gz") 
      output_format=".BAM"
      
      #---------------------------------------------------------------------------
      #align reads to the ref:                                                    |
      #---------------------------------------------------------------------------
      
      align.stat <- align(index="reference_index",
                          readfile1=files,
                          input_format = "gzFASTQ",
                          output_format = "BAM",
                          output_file=paste(read.files,"-subread",output_format,sep=""),
                          phredOffset=33 )
      align.stat 
    })

    output$Dwn_finish<-renderTable({
      if(is.null(input$reference_file$datapath)){
      return("Import the reference file with the FASTA format")
       }
      input$reference_file
       })
  
    lunch_alignment<-eventReactive(input$btn_lunch_alignment,{
      show.notification()
      res.alignment <- run.alignment()
      DT::datatable(as.data.frame(res.alignment))
    })
    
    output$alignment_results<-DT:: renderDataTable({
      remove.notification()
      lunch_alignment()
      }) 
  
  #**********************************************************************************************************************
    
    
    
    
    
  
  #**********************************************************************************************************************
  #*Quantification:
  #**********************************************************************************************************************
  
    #----------------------------------------------------------------------
    # Image about the annotation file:                                     |
    #----------------------------------------------------------------------

    output$AnnotationFile<-DT::renderDataTable({
      if(input$dw_GTF==""){
        return(as.matrix("Import an annotation file .GTF"))
      }
      GTF_File <- input$dw_GTF
      #GTF_File <-"C:/Users/DELL/Documents/Master/stage_pfa/NGS/My Work/Data-NGS/genome_assemblies_genome_gtf/ncbi-genomes-2022-09-21/GCF_003203755.1_ASM320375v1_genomic.gtf.gz"
      Annotation<- rtracklayer::import(GTF_File)
      dataGTF<-as.data.frame(Annotation)
      
      DT::datatable(dataGTF, options = list(scrollX = TRUE))
      })
    
    
    #----------------------------------------------------------------------------------------------------------
    # Flatten Features in GTF Annotation Files                                                                |
    #----------------------------------------------------------------------------------------------------------
    
    
  output$flatten_features<-DT::renderDataTable({
      if(input$dw_GTF==""){
        return(as.matrix("Import an annotation file .GTF"))
      }
    file <-input$dw_GTF
    annotatedData<-flattenGTF(GTFfile =  file)
    DT::datatable(as.data.frame(annotatedData),options =list(scrollX = TRUE) )
    
    })
    #--------------------------------------------------------------------------------------------------------------------------------------
    
    
    #----------------------------------------------------------------------------------------------------------
    ############################################### Feature Count  : ##########################################                                                                                       |
    #----------------------------------------------------------------------------------------------------------
    
    runCount<-eventReactive(input$btn_count,{
      show.notification()
      BAM_folder<-input$BAMFiles
     
      BAM_files<-list.files(full.names=TRUE, path=BAM_folder)
      AnnotationFile<-input$dw_GTF
      
      fc_SE <- featureCounts(BAM_files, isGTFAnnotationFile = TRUE, annot.ext=AnnotationFile)
      countMatrix<-as.data.frame(fc_SE$counts)
      colnames(countMatrix)<-t(import.metadata()['Run'])
      return(countMatrix)
      
      })
     
    output$featureCount<- DT::renderDataTable({
      if(input$BAMFiles==""){
        return(as.matrix("Import the BAM file(s)"))
      }
      if(is.null(input$dw_pheno_file)){
        return(as.matrix("Import the metadata of the study"))
      }
       remove.notification()
       DT::datatable(runCount(), options = list(scrollX = TRUE)) 
      })
    
    
    #-----------------------------------------------------------------------------------------------------
    ################################## Download the count matrix : #######################################                     
    #-----------------------------------------------------------------------------------------------------
    output$dwn_UI<-renderUI({
      downloadBttn("dwn_count_matrix",label = "Download the count matrix",icon = icon("download"),size = "sm",color = "warning")
    })
    
    output$dwn_count_matrix<- downloadHandler(
      filename = function() {
        paste("count-Matrix-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        write.csv(runCount(), file)
      }
    )
    #--------------------------------------------------------------------------------------------------------
  
    
    
    
    
    
    
    #**********************************************************************************************************************
    #*************************************Differential expression analysis********************************************    *
    #**********************************************************************************************************************
    
    
    #-----------------------------------------------------------------------------------------------------
    # 1- ############################ LOADING PHENOTYPE DATA: #############################              |
    #-----------------------------------------------------------------------------------------------------
    
    loading.pheno<-reactive({
      phenotype_file<-input$phenotype_data$datapath
      pheno <- read.csv(phenotype_file, row.names=1)
      pheno
    })
    
 
    #-----------------------------------------------------------------------------------------------------
    # 2- ############################ LOADING GENE EXPRESSION DATA: #############################        |
    #-----------------------------------------------------------------------------------------------------
    
    loading.GE.data<-reactive({
      gene_expression_file<-input$gene_expression_data$datapath
      data <- as.matrix(read.csv(gene_expression_file, row.names=1))
      data
    })
    
   
    
    #-----------------------------------------------------------------------------------------------------
    #############################  DESPLAY PHENOTYPE DATA GENE #############################             |
    #-----------------------------------------------------------------------------------------------------
    
    output$desc_phenotype_data<-renderUI({
      if(is.null(input$phenotype_data)){
        return( tags$span(style="color:red", "Select a phenotype data !"))
      }
      DT::dataTableOutput("phenotype_table")
    })
    output$phenotype_table<-DT::renderDataTable({
      DT::datatable(loading.pheno(),options = list(scrollX = TRUE))
    })
    
    #-----------------------------------------------------------------------------------------------------
    #############################  DESPLAY GENE EXPRESSION DATA: #############################           |
    #-----------------------------------------------------------------------------------------------------
    
    output$desc_gene_expression_data<-renderUI({
      if(is.null(input$gene_expression_data)){
        return( tags$span(style="color:red", "Select a count data !"))
      }
      selectedOption<-input$exp_gene_expression_data
      switch (selectedOption,
        "Data table"= DT::dataTableOutput("desc_expression_data"),
        "dimension of data" = tableOutput('dim'),
        "Explore the data distribution using the histogram plot"= plotOutput("data_hist")
         )
     })
    
    ########################################################################################
    # GE DATA : DATA DISTRIBUTION                                                          |
    ########################################################################################
    
    output$data_hist<-renderPlot({
      hist(log2(loading.GE.data()+1), col = "orange", main="The data distribution using the histogram plot")
    })
    
    #######################################################################################
    #GE DATA : DATA DIMENSION                                                             |
    #######################################################################################
    output$dim<-renderTable({
      dimData<-cbind(t(dim(loading.GE.data())))
      colnames(dimData)<-c("Number of rows(genes)","Number of columns(samples)")
      dimData
    })
    
    ########################################################################################
    # GE DATA: loading the gene expression data                                            |
    ########################################################################################
    output$desc_expression_data<-DT::renderDataTable({
      DT::datatable(loading.GE.data(),options = list(scrollX = TRUE))
    })
    
    
    
    
    #-----------------------------------------------------------------------------------------------------------
    ############################# DO the differential EXP analysis using DeSeq2 : ############################  |
    #-----------------------------------------------------------------------------------------------------------
    
    output$Diff_EXP<-renderUI({
      if(is.null(input$phenotype_data)){
        return( selectInput("select_conditions","Specify the conditions that you want to compare according to the phenotypic table",choices = c("select a phenotype data")))
      }
      selectInput("select_conditions","Specify the conditions that you want to compare according to the phenotypic table",choices =c("NULL",names(loading.pheno())) )
    })
    
    output$result_conditions<-renderUI({
      if(is.null(input$phenotype_data)){
        return( )
      }
      
      if(input$select_conditions=="select a phenotype data"){
        return(tags$span(style="color:red", "Select a phenotype data !"))
      }
      if(input$select_conditions=="NULL"){
        return(tags$span(style="color:red", "Select an attribute !"))
      }
      else
      tableOutput("tab_conditions")
    })
    
    output$tab_conditions<-renderTable({
      selected.condition <- input$select_conditions
      stat.condition <- table(loading.pheno()[selected.condition])
      stat.condition
    })
   
    #-----------------------------------------------------------------------------------------------------------
    #********************************** Run DESq2: Differential expression analysis *************************** |
    #-----------------------------------------------------------------------------------------------------------
    

    run.DEGs<-reactive({
      selected.condition<-input$select_conditions
      selected.FC<-input$foldChange
      selected.design <- loading.pheno()[selected.condition]
      res.DESq2<-DEGs(countData=loading.GE.data(),colData=loading.pheno(), design=selected.design, FoldChange=selected.FC, condition= selected.condition)
      res.DESq2
    })
    

    runDESq2<-eventReactive(input$btn_run_DFEXR,{
        show.notification()
        run.DEGs()
    })
    
    output$resultDESq2<- DT::renderDataTable({
      remove.notification()
      DT::datatable(runDESq2(),options = list(scrollX = TRUE))
    })
    
    
    #-----------------------------------------------------------------------------------------------------
    #Download the DEGS Data:                                                                          |
    #-----------------------------------------------------------------------------------------------------
    
    output$btn_dwn_Degs<- downloadHandler(
      filename = function() {
        paste("DEGs-Analysis-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        write.csv(as.matrix(run.DEGs()),file, quote=F,row.names=T)
      }
    )

    #******************************************************************************************************************************************

})

