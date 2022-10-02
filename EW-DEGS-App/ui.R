library(shiny)
library(DT)
library(shinydashboard)
library(shinyWidgets)

#*********************************************************************************************************
###################################Define UI for application. ###########################################
#*********************************************************************************************************

shinyUI(

   dashboardPage( skin="yellow",
     #-----------------------------------------------------------------------            
     #Header:
     #-----------------------------------------------------------------------
     
     dashboardHeader(title="Easy workflow for NGS analysis", titleWidth = 450),
     
     #*****************************************************************************************************
     #dashboardSidebar
     #*****************************************************************************************************
     
     dashboardSidebar(
       width = 250,
       sidebarMenu(
       menuItem("getting DATA" ,tabName = "collected_data", icon = icon("download",lib = "font-awesome")),
       menuItem("Pre-processing" , 
        menuSubItem("Quality check",tabName = "quality", icon= icon("vial-circle-check")),
        menuSubItem("Filtering and trimming", tabName = "filtering", icon= icon("filter"))
        ),
       menuItem("Read mapping or alignment", tabName = "alignment", icon = icon("indent")),
       menuItem("Quantification", tabName = "quantification", icon = icon("table-list")),
       menuItem("Differential expression analysis", tabName = "DESeq2", icon = icon("dna"))
       )
       ),
     #******************************************************************************************************

     
     #******************************************************************************************************
     #*DASHBOARDBODY:
     #******************************************************************************************************
     
   dashboardBody(
     
     tags$head(tags$style(HTML('
      .main-header .logo {
        font-family: "Georgia", Times, "Times New Roman", serif;
        font-weight: bold;
        font-size: 23px;
      }'))),
     tags$head(tags$style(".progress-bar{background-color:#222d32;}")),
     tags$head(
       tags$style(
         HTML(".shiny-notification {
              height: 100px;
              width: 800px;
              position:fixed;
              top: calc(50% - 50px);;
              left: calc(50% - 400px);;} "))),
     
     tabItems(
       #*****************************************************************************************************
       ####################################Collecting data from the Computer :##############################
       #*****************************************************************************************************
       tabItem(tabName = "collected_data",
               #------------------------------------------------------------------------------------------------------------------
               ####################################################Box ONE#######################################################
               #------------------------------------------------------------------------------------------------------------------
                   box( width= FALSE, status = "warning", solidHeader = TRUE, title="Import the SRA metadata of the study",
                        tags$h4(tags$b(tags$tspan(style="color:gray","The SRA metadata"))),
                        tags$i(class = "fa-solid fa-thumbtack"),
                        tags$span(style="font-size:17px","The SRA metadata describes the technical aspects of sequencing experiments: the sequencing libraries, preparation techniques and data files.It is therefore imperative that submitters provide clear and informative Title and Description for each EXPERIMENT. according to the metadata fields your data will be split into groups of samples."), 
                        br(),hr(),
                        sidebarLayout(
                          sidebarPanel(fileInput("dw_pheno_file","Import phenotype data or metadata of the study"),uiOutput("select_gr_field")),
                          mainPanel(DT::dataTableOutput("pheno_description"))
                        ), #END SIDEBARLAYOUT.
                        br()),
               
               #------------------------------------------------------------------------------------------------------------------
               ####################################################Box TWO#######################################################
               #------------------------------------------------------------------------------------------------------------------
                   box(width=FALSE, title="Import the SRR data (SRA RUN files(FASTQ files)) ",status = "warning", solidHeader = TRUE,
                       textInputIcon("folderName", label="Enter the directory path that contains fastq input files",placeholder = "Enter the directory path",icon = icon("folder-open")),
                       (tags$span(style="color:blue", "Make sure that the directory path contains the specific FASTQ file only!")),br(),
                       tags$span(style="color:blue", "Attention! the file path must contain '/' , Example  :  C:/Users/username/path/to/file(s).fastq.gz"),
                       br(),br(),
                       actionBttn("btn_upload",label = "Upload fastq files", icon = icon("download"),color = "warning",size = "sm",no_outline = TRUE),br(),br(),
                       tableOutput("file_description"),
                       ),
               
               #------------------------------------------------------------------------------------------------------------------
               ####################################################Box THREE#######################################################
               #------------------------------------------------------------------------------------------------------------------
                   box(width=FALSE, title="per File Top Reads",status = "warning", solidHeader = TRUE,
                   tags$span(style="font-size:17px"," PerFileTopReads Most represented sequencing reads and their counts."),br(),br(), 
                   actionBttn("btn_perFileTopReads",label = "perFileTopReads", icon = icon("play"),color = "warning", size = "sm", no_outline = TRUE),  
                   hr(),
                   DT::dataTableOutput("perFileTopReads"))
               ),
      
       #*********************************************************************************************************************************
       
       
       #******************************************************************************************************
       ###################################### Quality tab: ##################################################
       #******************************************************************************************************
       tabItem(tabName = "quality",
             box( width = FALSE, title="Quality check on sequencing reads", status = "warning", solidHeader = TRUE, height = "100%",
             tags$h4(tags$b(tags$tspan(style="color:gray","Why doing quality check ?"))),
             tags$i(class = "fa-solid fa-thumbtack"),
             tags$span(style="font-size:17px","The sequencing technologies usually produce basecalls with varying quality.In addition, there could be sample-specific issues in your sequencing run, such as adapter contamination. It is standard procedure to check the quality of the reads and identify problems before doing further analysis. Checking the quality and making some decisions for the downstream analysis can influence the outcome of your project."), 
             hr(),br(),
             sidebarLayout(
                   # Sidebar with a slider input
                   sidebarPanel(
                     #----------------------------------------------------------
                     #Select a FASTQ File                                       |
                     #----------------------------------------------------------
                     uiOutput("select_fastq_file"),
                     hr(),
                     #----------------------------------------------------------
                     #Setect a Plot                                             |
                     #----------------------------------------------------------
                     
                     selectInput("select_plot", "Available Graphics",choices = (
                       c("Per Read Mean Quality Distribution of Files", 
                         "Average Quality",
                         "Cycle-specific Average Quality",
                         "Read Frequency",
                         "Read Length Distribution",
                         "Cycle-specific GC Content",
                         "Cycle-specific Quality Distribution",
                         "Cycle-specific Quality Distribution - Boxplot",
                         "Cycle-specific Base Call Proportion"))),
                     hr(),
                     tags$span(style="font-size:17px","Generates an HTML report file."), 
                     br(),br(),
                     actionBttn("btn_QR","Generate",icon=icon("play"), size = "sm",color = "warning"),
                     br(),br(),
                     textOutput("qualityReport"),
                     br()
                   ),
                   #--------------------------------------------------------------
                   # Show a plot of the generated distribution                   |
                   #--------------------------------------------------------------
                   mainPanel(
                     fluidPage( tabsetPanel(
                       tabPanel("Data graphic", plotOutput("quality_plot")),
                       tabPanel("Data", DT::dataTableOutput("quality_data")),
                       ))
                     
                     ))) ),
       #************************************************************************************************************************************************
       
       
       
       
       
       #****************************************************************************************************
       ############################################## Filtering tab: #######################################
       #****************************************************************************************************
       tabItem(tabName = "filtering",
               box(width= FALSE, status = "warning", solidHeader = TRUE, title=" Filtering and trimming reads",
                   tags$i(class = "fa-solid fa-thumbtack"),
                   tags$span(style="font-size:16px","Based on the results of the quality check, you may want to trim or filter the reads. The quality check might have shown the number of reads that have low quality scores. These reads will probably not align very well because of the potential mistakes in base calling, or they may align to wrong places in the genome.Therefore, you may want to remove these reads from your fastq file."), 
                   hr(),
                   tags$h4(tags$b(tags$tspan(style="color:gray","Truncate sequences, remove parts matching to adapters and filter out low quality or low complexity sequences from (compressed) ’fastq’ files."))),
                   br(),
                   splitLayout( 
                     box( width= FALSE, solidHeader = FALSE,
                     textInput("truncateStartBases","The number of bases to be truncated (removed) from the beginning of each sequence.",placeholder = "truncate Start Bases"),
                     textInput("truncateEndBases","The number of bases to be truncated (removed) from the end of each sequence.",placeholder = "truncate End Bases"),
                     textInput("Lpattern","The left (5’-end) adapter sequence.",placeholder = "Left pattern"),
                     textInput("Rpattern","The right (3’-end) adapter sequence.",placeholder = "Right pattern")),
                     box( width= FALSE, solidHeader = FALSE,
                      textInput("minLength","The minimal allowed sequence length.",placeholder = "min Length"),
                     textInput("nBases"," The maximal number of Ns allowed per sequence.",placeholder = "N Bases"),
                     textInput("complexity","The minimal sequence complexity",placeholder = "complexity"))),
                   helpText("Check your filtered data in the same directory that contains non-filtered fastq input file(s) "),
                   actionBttn("btn_filter",label = "Filter and Trim reads", icon = icon("filter"),color = "warning",size = "sm", no_outline = TRUE),
                   br(),br(),
                   #-------------------------------------------------Enter the directory path that contains fastq input files--------------------
                   #Show the results of the filtering and trimming:                      |
                   #---------------------------------------------------------------------
                   DT::dataTableOutput("tb_filter_trim"),                                
                   #---------------------------------------------------------------------
                   
                    )
       ),
       #----------------------------------------------------------------------------------------------------------------------------
       #**************************************************************************************************************************************
       
       
       
       
       
       
       
       #*****************************************************************************************************************************
       #Read mapping or alignment
       #*****************************************************************************************************************************
       
       tabItem(tabName = "alignment",
               box( width= FALSE, title = "Mapping your reads to a reference genome", status= "warning", solidHeader = TRUE,
               fileInput("reference_file","Download a reference file:fasta file"),
               tags$b(span(style="color:gray", "Bref description of the reference file")),br(),
               tableOutput("Dwn_finish"),
               hr(),
               textInputIcon("filtered_data_folder","Enter the directory path that contains your filtered data from the last step of trimming and filtering",placeholder = "Enter the directory path",icon = icon("folder-open")),
               tags$span(style="color:blue", "Make sure that the directory path contains the specific FASTQ file(s) only!"),br(),
               tags$span(style="color:blue", "Attention! the file path must contain '/' , Example  :  C:/Users/username/path/to/Filtered-file(s).fastq.gz"),br(),br(),
               actionBttn('btn_lunch_alignment',"Launch alignment",icon=icon("indent"), size="sm", color="warning"),
                    hr(),
                    #Show the result of alignment:
                    DT::dataTableOutput("alignment_results")
               )
       ),
       #--------------------------------------------------------------------------------------------------------------------------------
       
       
       
       #********************************************************************************************************
       #*-------------------------------------------------------------------------------------
       ######################################## Quantification: ##############################
       #--------------------------------------------------------------------------------------
       #********************************************************************************************************
       tabItem(tabName = "quantification",
               #----------------------------------------------------------------------------------------------------------------
               #Box One                                                                                                         |
               #----------------------------------------------------------------------------------------------------------------
               box( width= FALSE, title = "Quatification : A program for assigning sequence reads to genomic features", status= "warning", solidHeader = TRUE,
                verticalLayout(
                box( width = FALSE,
                      textInputIcon("dw_GTF", label="Enter the absolute path to the .gtf ( gtf.gz ) annotation file",placeholder = "Enter the path of the gtf or gtf.gz  file",icon = icon("folder-open")),
                      tags$span(style="color:red", "Attention! the file path must contain '/' , Example  :  C:/Users/username/path/to/file.gtf") 
                ),
                box( width = FALSE,
                      tags$b(span(style="color:gray", "Description of the annotation file uploaded")),br(),br(),
                      DT::dataTableOutput("AnnotationFile"),
                      hr(), 
                ))),
               #-----------------------------------------------------------------------------------------------------------------
               
               #-------------------------------------------------------------------------------------------------------
               #Box two                                                                                                          |
               #-----------------------------------------------------------------------------------------------------------------
               box( width= FALSE, title = "Flatten Features in GTF  Annotation File", status= "warning", solidHeader = TRUE,
               tags$h3(tags$b(tags$tspan("Looking for 'exon' features... (grouped by 'gene_id')"))),
               DT::dataTableOutput("flatten_features"),
               hr()
               ),
               #-----------------------------------------------------------------------------------------------------------------
               
               
               #-----------------------------------------------------------------------------------------------------------------
               #Box Three                                                                                                         |
               #-----------------------------------------------------------------------------------------------------------------
               box( width= FALSE, title = "Feature Count", status= "warning", solidHeader = TRUE,
               textInputIcon("BAMFiles", label="The function takes as input a set of SAM or BAM files containing read mapping results.",placeholder = "Enter the path of the folder",icon = icon("folder-open")),
               tags$span(style="color:red", "Attention! make sure that the folder contains .BAM files only"),
               br(),
               br(),
               actionBttn("btn_count",label = "Run Count", icon = icon("play"),color = "warning",size = "sm", no_outline = TRUE),
               br(),
               br(),
               tags$b(span(style="color:gray", "The result of the featureCounts")),br(),br(),
               DT::dataTableOutput("featureCount"),
               hr(),
               uiOutput("dwn_UI")
              
               )
    
       ),
       
       #***************************************************************************************************************************
       
       
       
       #***************************************************************************************************************************
       #*-----------------------------------------------------------------------------------------------------------------------
       #*Differential expression genes:                                                                                         |
       #*-----------------------------------------------------------------------------------------------------------------------
       #***************************************************************************************************************************
       
       tabItem(tabName="DESeq2",
               
               #------------------------------------------------------------------------------------------------------------------
               ############################################# BOX INFORMATION : #################################################  |
               #------------------------------------------------------------------------------------------------------------------
               
               box( width= FALSE, title = "Differential expression analysis", status= "warning", solidHeader = TRUE,
                    tags$span(style="font-size:17px","A basic task in the analysis of count data from RNA-seq is the detection of differentially expressed genes.An important analysis question is the quantification and statistical inference of systematic changes between conditions, as compared to within-condition variability. The package DESeq2 provides methods to test for differential expression by use of negative binomial generalized linear models; the estimates of dispersion and logarithmic fold changes incorporate data-driven prior distributions. "),
                    tags$h3(tags$b(tags$tspan("Praparing data"))),
                    tags$h4(tags$b(tags$tspan(style="color:gray","Count data/gene expression data"))),
                    tags$span(style="font-size:17px","The count data are presented as a table which reports, for each sample, the number of sequence fragments that have been assigned to each gene."),
                    tags$h4(tags$b(tags$tspan(style="color:gray","Phenotype data fo analysis : "))),
                    tags$span(style="font-size:17px","A table of sample information"),
                    hr(),
                    tags$i(class = "fa-solid fa-warning"),
                    tags$b(style="color:red","It is very important to check manually that the columns of the count matrix correspond to the rows of the sample information table."),
                    hr()),
               
               splitLayout( 
               #------------------------------------------------------------------------------------------------------------------
               ##############################GENE EXPRESSION DATA : ######################################                       |
               #------------------------------------------------------------------------------------------------------------------
               
               box( width= FALSE, title = "Loading the gene expression data", status= "warning", solidHeader = TRUE,
               hr(),
               fileInput("gene_expression_data","Loading the gene expression data"),
               radioButtons( "exp_gene_expression_data", "Explore gene expression data",choices =c("Data table","dimension of data","Explore the data distribution using the histogram plot" )),
               tags$h4(tags$b(tags$tspan("Description of gene expression data"))),
               hr(),
               uiOutput("desc_gene_expression_data"),
               hr(),
               br(),
               ),
               
               #------------------------------------------------------------------------------------------------------------------
               ############################################### PHENOTYPE DATA : #########################################        |
               #------------------------------------------------------------------------------------------------------------------
               box( width= FALSE, title = "Loading the phenotype data", status= "warning", solidHeader = TRUE,
                   hr(),
                   fileInput("phenotype_data","Loading the phenotype data"),
                   tags$h4(tags$b(tags$tspan("Description of the phenotype data"))),
                   hr(),
                   uiOutput("desc_phenotype_data"),
                   hr(),
                   br(),
               ), 
               ), ############### END Split Layout ##############
               #-------------------------------------------------------------------------------------------------------------------         
      
               #------------------------------------------------------------------------------------------------------------------
               ################################DO the differential EXP analysis using DeSeq2######################################                                          |
               #------------------------------------------------------------------------------------------------------------------
               
               box( width= FALSE, title = "DO the differential EXP analysis using DeSeq2", status= "warning", solidHeader = TRUE,
                    
                    sidebarLayout(
                      sidebarPanel(
                        uiOutput("Diff_EXP"),
                        br(),
                        textInput("foldChange","Chose the value of the Fold Change",placeholder = "Fold change"),
                        tags$details("Chose the statstical significant differentaily expressed genes (DEGs) based on the p adjusted value less than 0.05 and biological significance  based on the fold change more than 'Fold Change'"),
                        br(),
                        actionBttn("btn_run_DFEXR","run analysis",icon=icon("play"),size = "sm",color = "warning"),
                        downloadBttn("btn_dwn_Degs","export the Degs",icon=icon("download"),size = "sm",color = "warning")
                      ),
                      mainPanel(
                        tags$h4(tags$b(tags$tspan("RUN"))),
                        hr(),
                        uiOutput("result_conditions"),
                        br(),
                        
                        DT::dataTableOutput("resultDESq2")
                      )),
                  
                    
               )
       )
       
       
       #***************************************************************************************************************************
       
       
       
      
       
     )
     
     )))


