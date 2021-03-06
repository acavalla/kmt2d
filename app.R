#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
###################################################################################################################

library(shiny)
library(shinyjs)
library(DT)
library(fmsb)
library(gridExtra)
library(stringi)
library(dplyr)
library(ggplot2)
library(ggrepel)

##load csvs of data and amalgamate
file.names <- c("\\\\ye2svr02\\Erlwood_NTEG\\CavallaA\\Microglia\\bigsevenmicrog.txt",
                "\\\\ye2svr02\\Erlwood_NTEG\\Psychiatric Genetics\\Projects\\Phago\\bigseven.txt",
                "\\\\ye2svr02\\Erlwood_NTEG\\CavallaA\\Microglia\\othersources.txt",
                "\\\\ye2svr02\\Erlwood_NTEG\\Psychiatric Genetics\\Projects\\Phago\\extras.txt")
for (i in 1:length(file.names)){
  x <- read.csv(file.names[[i]], sep="\t", row.names = NULL, stringsAsFactors = FALSE)
  if(i == 1){
    allgenes <- x
  } else {
    colnames(x) <- colnames(allgenes)
    allgenes <- rbind(allgenes, x)
  }
  rm(x)
}

##round to two sig figs

allgenes[,7:8] <- allgenes[,7:8] %>% 
  mutate_if(is.numeric, signif, digits = 3)
colnames(allgenes) <- gsub("\\_", " ", colnames(allgenes))
allgenes$"Paper title" <- paste0("<a href='",allgenes$Link,"' target='_blank'>",allgenes$"Paper title","</a>")
allgenes$LogFC <- gsub("#VALUE!", "Inf", allgenes$LogFC)
allgenes <- allgenes[order(as.character(allgenes$Gene), allgenes$Year),]
allgenes$Gene <- trimws(as.character(allgenes$Gene))
allgenes <- unique(allgenes)
rownames(allgenes) <- 1:nrow(allgenes)
allgenes <- rbind(allgenes[329,], allgenes[c(1:328, 330:nrow(allgenes)),])
rownames(allgenes) <- 1:nrow(allgenes)

#Create datatable for other sources
extras <- allgenes[!grepl("Kunkle|Ellis|Levey|Lancour|Marioni|Galatro|Jansen", allgenes$Author),1:4]

#Kelley radar plots data
noma <- read.csv("\\\\ye2svr02\\Erlwood_NTEG\\CavallaA\\Microglia\\All_kelley.txt", sep="\t")

#Load cano data
cano <- read.csv("\\\\ye2svr02\\Erlwood_NTEG\\CavallaA\\Microglia\\galatro_diffex.txt", sep="\t", row.names = NULL, stringsAsFactors = FALSE)
cano <- cano[,c(2,12,13)]

# logifySlider javascript function
JS.logify <-
  "
// function to logify a sliderInput
function logifySlider (sliderId, sci = false) {
if (sci) {
// scientific style
$('#'+sliderId).data('ionRangeSlider').update({
'prettify': function (num) { return ('10<sup>'+num+'</sup>'); }
})
} else {
// regular number style
$('#'+sliderId).data('ionRangeSlider').update({
'prettify': function (num) { return (Math.pow(10, num)); }
})
}
}"

# call logifySlider for each relevant sliderInput
JS.onload <-
  "
// execute upon document loading
$(document).ready(function() {
// wait a few ms to allow other scripts to execute
setTimeout(function() {
// include call for each slider
logifySlider('p', sci = false)
logifySlider('log_slider2', sci = true)
}, 5)})
"

# Define UI for application that draws a histogram
ui <- fluidPage(
  useShinyjs(),
  
  # Application title
  tabsetPanel(
    tabPanel("Dropdown list",
             titlePanel("Alzheimer Disease Genes"),
             
             # Show a plot of the generated distribution
             
             verticalLayout(selectizeInput("gene",
                                           "Select gene of interest",
                                           choices = allgenes$Gene,
                                           options = list(maxOptions = 20000)),
                            em("Type to search, or click to select from the list")),
             h1("microglia"),
             splitLayout(cellWidths = c("17%", "4%", "66%", "15%"),
                         HTML("study"),
                         HTML(""),
                         HTML("Results"),
                         HTML("Significant?"),
                         cellArgs = list(style='white-space: normal;')),
             splitLayout(cellWidths = c("17%", "4%", "71%", "10%"),
                         HTML(paste0("<a href='",allgenes[grepl("Galatro", allgenes$Author), 5][1],"' target='_blank'>Galatro</a>")),
                         uiOutput("booGal"),
                         textOutput("resGal"),
                         uiOutput("sigGal"),
                         cellArgs = list(style='white-space: normal;')),
             splitLayout(cellWidths = c("17%", "4%", "71%", "10%"),
                         HTML(paste0("<a href='",allgenes[grepl("Ellis", allgenes$Author), 5][1],"' target='_blank'>Olah 2018a</a>")),
                         uiOutput("booOlah"),
                         textOutput("resOlah"),
                         uiOutput("sigOlah"),
                         cellArgs = list(style='white-space: normal;')),
             
             h1("genetics"),
             splitLayout(cellWidths = c("17%", "3%", "20%", "10%", "17%", "3%", "20%", "10%"),
                         HTML("study"),
                         HTML(""),
                         HTML("Results"),
                         HTML("Significant?"),
                         HTML("study"),
                         HTML(""),
                         HTML("Results"),
                         HTML("Significant?"),
                         cellArgs = list(style='white-space: normal;')),
             splitLayout(cellWidths = c("17%", "3%", "25%", "5%", "17%", "3%", "25%", "5%"),
                         HTML(paste0("<a href='",allgenes[grepl("Levey", allgenes$Author), 5][1],"' target='_blank'>IGAP MAGMA</a>")),
                         uiOutput("booIgapMagma"),
                         textOutput("resIgapMagma"),
                         uiOutput("sigIgapMagma"),
                         HTML(paste0("<a href='",allgenes[grepl("Marioni", allgenes$Author), 5][1],"' target='_blank'>Marioni famhis</a>")),
                         uiOutput("booMfamhis"),
                         textOutput("resMfamhis"),
                         uiOutput("sigMfamhis"),
                         cellArgs = list(style='white-space: normal;')),
             splitLayout(cellWidths = c("17%", "3%", "25%", "5%", "17%", "3%", "25%", "5%"),
                         HTML(paste0("<a href='",allgenes[grepl("Lancour", allgenes$Author), 5][1],"' target='_blank'>Lancour</a>")),
                         uiOutput("booLancour"),
                         HTML("N/A"),
                         HTML("N/A"),
                         HTML(paste0("<a href='",allgenes[grepl("Kunkle", allgenes$Author), 5][1],"' target='_blank'>Kunkle GWAS</a>")),
                         uiOutput("booK"),
                         textOutput("resK"),
                         uiOutput("sigK"),
                         cellArgs = list(style='white-space: normal;')),
             splitLayout(cellWidths = c("17%", "3%", "25%", "5%", "17%", "3%", "25%", "5%"),
                         HTML(paste0("<a href='",allgenes[grepl("Jansen", allgenes$Author), 5][1],"' target='_blank'>Jansen</a>")),
                         uiOutput("booJans"),
                         textOutput("resJans"),
                         uiOutput("sigJans"),
                         HTML(""),
                         HTML(""),
                         HTML(""),
                         HTML(""),
                         cellArgs = list(style='white-space: normal;')),
             
             tags$head(tags$script(HTML(JS.logify))),
             tags$head(tags$script(HTML(JS.onload))),
             br(),
             splitLayout(sliderInput("p", "Significant P value",
                                     min = -10, max = 1, value = -4, step = 1),
                         sliderInput("FC", "Significant log fold change",
                                     min = 0, max = 6, value = 1),
                         cellArgs = list(style='white-space: normal;')),
             DTOutput('selectedgene'),
             splitLayout(cellWidths = c("75%", "15%", "10%"),
                         plotOutput("radarPlot"),
                         htmlOutput("mytext"),
                         uiOutput("booKel"),
                         cellArgs = list(style='white-space: normal;')),
             actionLink("button", "What is a radar plot?"),
             hidden(
               div(id='text_div',
                   verbatimTextOutput("text")
               )
             ),
             #conditionalPanel(condition = "cano[cano$Gene.Symbol==input.gene,]) > 0",
             plotOutput("cano")),
    tabPanel("Batch input",
             textInput("caption", "Input gene list", ""),
             fileInput("file1", "Choose file",
                       multiple = FALSE,
                       accept = c("text/csv",
                                  "text/comma-separated-values,text/plain",
                                  ".csv")),
             tags$head(tags$style("#totalnrow{color: red;
                                  font-size: 20px;
                                  font-style: italic;
                                  }"
                         )
             ),
             textOutput("totalnrow"),
             DTOutput("tab2table")),
    tabPanel("Documentation",
             "If you would like to add a new gene list, it must contain the following columns:",
             br(),
             "1. Gene",
             br(),
             "2. Paper title",
             br(),
             "3. Author",
             br(),
             "4. Year",
             br(),
             "5. Link",
             br(),
             "Lancour is a meta-analysis.",
             br(),
             "Input files for the batch input tab should be in .csv format, with no header row."
    )
             )
    )


# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  #### --- IN TAB 1 --- ####
  
  ### Is selected gene present in any of the seven specified papers?
  ### ["&#x2715;" is the cross; "&#x2713;" the tick]
  output$booGal <- renderText({
    ifelse(nrow(allgenes[grepl("Galatro", allgenes$Author) & allgenes$Gene == input$gene,])==0, "&#x2715;", "&#x2713;")
  })
  output$booOlah <- renderText({
    ifelse(nrow(allgenes[grepl("Ellis", allgenes$Author) & allgenes$Gene == input$gene,])==0, "&#x2715;", "&#x2713;")
  })
  output$booIgapMagma <- renderText({
    ifelse(nrow(allgenes[grepl("Levey", allgenes$Author) & allgenes$Gene == input$gene,])==0, "&#x2715;", "&#x2713;")
  })
  output$booMfamhis <- renderText({
    ifelse(nrow(allgenes[grepl("Marioni", allgenes$Author) & allgenes$Gene == input$gene,])==0, "&#x2715;", "&#x2713;")
  })
  output$booLancour <- renderText({
    ifelse(nrow(allgenes[grepl("Lancour", allgenes$Author) & allgenes$Gene == input$gene,])==0, "&#x2715;", "&#x2713;")
  })
  output$booK <- renderText({
    ifelse(nrow(allgenes[grepl("Kunkle", allgenes$Author) & allgenes$Gene == input$gene,])==0, "&#x2715;", "&#x2713;")
  })
  output$booJans <- renderText({
    ifelse(nrow(allgenes[grepl("Jansen", allgenes$Author) & allgenes$Gene == input$gene,])==0, "&#x2715;", "&#x2713;")
  })
  
  ### What are the results for selected gene, if present?
  output$resGal <- renderText({
    req(nrow(allgenes[grepl("Galatro", allgenes$Author) & allgenes$Gene == input$gene,]) > 0)
    paste0("RPKM = ",
           allgenes[grepl("Galatro", allgenes$Author) & allgenes$Gene == input$gene, 6],
           "; LogFC = ",
           allgenes[grepl("Galatro", allgenes$Author) & allgenes$Gene == input$gene, 7],
           "; P = ",
           allgenes[grepl("Galatro", allgenes$Author) & allgenes$Gene == input$gene, 8])
  })
  
  output$resOlah <- renderText({
    req(nrow(allgenes[grepl("Ellis", allgenes$Author) & allgenes$Gene == input$gene,]) > 0)
    paste0("RPKM = ",
           allgenes[grepl("Ellis", allgenes$Author) & allgenes$Gene == input$gene, 6],
           "; LogFC = ",
           allgenes[grepl("Ellis", allgenes$Author) & allgenes$Gene == input$gene, 7],
           "; P = ",
           allgenes[grepl("Ellis", allgenes$Author) & allgenes$Gene == input$gene, 8])
  })
  
  output$resIgapMagma <- renderText({
    req(nrow(allgenes[grepl("Levey", allgenes$Author) & allgenes$Gene == input$gene,]) > 0)
    paste0("P = ",
           allgenes[grepl("Levey", allgenes$Author) & allgenes$Gene == input$gene, 8])
  })
  output$resMfamhis <- renderText({
    req(nrow(allgenes[grepl("Marioni", allgenes$Author) & allgenes$Gene == input$gene,]) > 0)
    paste0("P = ",
           allgenes[grepl("Marioni", allgenes$Author) & allgenes$Gene == input$gene, 8])
  })
  output$resLancour <- renderText({
    req(nrow(allgenes[grepl("Lancour", allgenes$Author) & allgenes$Gene == input$gene,]) > 0)
    paste0("P = ",
           allgenes[grepl("Lancour", allgenes$Author) & allgenes$Gene == input$gene, 8])
  })
  output$resK <- renderText({
    req(nrow(allgenes[grepl("Kunkle", allgenes$Author) & allgenes$Gene == input$gene,]) > 0)
    paste0("P = ",
           allgenes[grepl("Kunkle", allgenes$Author) & allgenes$Gene == input$gene, 8])
  })
  output$resJans <- renderText({
    req(nrow(allgenes[grepl("Jansen", allgenes$Author) & allgenes$Gene == input$gene,]) > 0)
    paste0("P = ",
           allgenes[grepl("Jansen", allgenes$Author) & allgenes$Gene == input$gene, 8])
  })
  
  ### Are selected gene's results significant?
  output$sigGal <- renderText({
    ifelse(nrow(allgenes[grepl("Galatro", allgenes$Author) & allgenes$Gene == input$gene & allgenes$P <= 10^input$p & allgenes$LogFC>=input$FC,])==0, "&#x2715;", "&#x2713;")
  })
  output$sigOlah <- renderText({
    ifelse(nrow(allgenes[grepl("Ellis", allgenes$Author) & allgenes$Gene == input$gene & allgenes$P <= 10^input$p & allgenes$LogFC>=input$FC,])==0, "&#x2715;", "&#x2713;")
  })
  output$sigIgapMagma <- renderText({
    ifelse(nrow(allgenes[grepl("Levey", allgenes$Author) & allgenes$Gene == input$gene & allgenes$P <= 10^input$p & allgenes$LogFC>=input$FC,])==0, "&#x2715;", "&#x2713;")
  })
  output$sigMfamhis <- renderText({
    ifelse(nrow(allgenes[grepl("Marioni", allgenes$Author) & allgenes$Gene == input$gene & allgenes$P <= 10^input$p & allgenes$LogFC>=input$FC,])==0, "&#x2715;", "&#x2713;")
  })
  output$sigK <- renderText({
    ifelse(nrow(allgenes[grepl("Kunkle", allgenes$Author) & allgenes$Gene == input$gene & allgenes$P <= 10^input$p & allgenes$LogFC>=input$FC,])==0, "&#x2715;", "&#x2713;")
  })
  output$sigJans <- renderText({
    ifelse(nrow(allgenes[grepl("Jansen", allgenes$Author) & allgenes$Gene == input$gene & allgenes$P <= 10^input$p & allgenes$LogFC>=input$FC,])==0, "&#x2715;", "&#x2713;")
  })
  
  ### Radar plot
  #output$rplotboo <- nrow(noma[noma$Gene==input$gene,])
  
  output$radarPlot <- renderPlot({
    req(nrow(noma[noma$Gene==input$gene,]) > 0)
    a <- noma[noma$Gene==input$gene,]
    a <- a[,2:5]
    # fmsb requires max and min of each category
    a <- rbind(rep(200, ncol(a)) , rep(0, ncol(a)) , a)
    labels <- gsub("_.*","", colnames(a))
    
    sizes <- as.data.frame(sort(a[3,]))
    sizes[1:3] <- 1
    sizes[4] <- 1.5
    sizes <- c(sizes$Astrocyte_Fidelity, sizes$Microglia_Fidelity, sizes$Oligodendrocyte_Fidelity, sizes$Neuron_Fidelity)
    colors_border = ifelse(sizes[1] != 1, rgb(0.141, 0.690, 0.760, alpha = 0.7),
                           ifelse(sizes[2] != 1, rgb(0.925, 0.572, 0.294, alpha = 0.7),
                                  ifelse(sizes[3]!=1, rgb(0.168, 0.592, 0.403, alpha = 0.7),
                                         rgb(0.780, 0.321, 0.878, alpha = 0.7))))
    colors_in = ifelse(sizes[1] != 1, rgb(0.360, 0.768, 0.819, alpha = 0.6),
                       ifelse(sizes[2] != 1, rgb(0.960, 0.768, 0.619, alpha = 0.6),
                              ifelse(sizes[3]!=1, rgb(0.650, 0.929, 0.807, alpha = 0.6),
                                     rgb(0.949, 0.827, 0.972, alpha = 0.6))))
    
    radarchart(a,
               axistype = 4 ,
               #custom polygon
               pcol = colors_border ,
               pfcol = colors_in ,
               plwd = 4 ,
               plty = 1,
               #custom the grid
               cglty = 3,
               cglcol = "light grey",
               caxislabels = c(0, "","","", 200),
               axislabcol = "black",
               cglwd = 0.8,
               vlabels = labels,
               vlcex = sizes)
    
  })
  
  ##Kelley data 
  output$mytext <- renderUI({
    mylist <- c("", "", "", "", "", "", "", "", "", "", "", "Is the gene in the top 50 of the Kelley data?")
    HTML(paste(mylist, sep = "", collapse = '<br/>'))
  })
  
  output$booKel <- renderUI({
    x <- ifelse(nrow(allgenes[grepl("Kelley", allgenes$Author) & allgenes$Gene == input$gene,])==0, "&#x2715;", "&#x2713;")
    mylist <- c("", "", "", "", "", "", "", "", "", "", "", x)
    HTML(paste(mylist, sep = "", collapse = '<br/>'))
  })
  
  ### Explanatory radar plot button
  observeEvent(input$button, {
    toggle('text_div')
    output$text <- renderText({"The radar plots are created by the data reported by Kelley et al (2018), which is presented at http://oldhamlab.ctec.ucsf.edu/. The axes represent a 'fidelity score', another name for a Z-score, which is defined as the number of standard deviations a data point is from the mean, which is calculated over all brain tissue. A z-score of a gene of 150 in microglia would therefore represent an upregulation of 150 standard deviations over the mean for this gene in microglia in combined brain tissue."})
    
    if (input$button %% 2 == 1) {
      txt <- "Hide"
    } else {
      txt <- "What is a radar plot?"
    }
    updateActionButton(session, "button", label = txt)
  })
  
  ### Volcano plot
  #output$canplotboo <- nrow(cano[cano$Gene.Symbol==input$gene,])
  
  output$cano <- renderPlot({
    req(nrow(cano[cano$Gene.Symbol==input$gene,]) > 0)
    cano<-cano%>%mutate(threshold = ifelse(abs(glia.brain_logFC)>=input$FC & P.Value<=10^input$p, "A",
                                           ifelse(abs(glia.brain_logFC)>=input$FC, "B",
                                                  ifelse(P.Value<=10^input$p, "C", "D"))))
    
    cano<-cano%>%mutate(callit = ifelse(cano$Gene.Symbol==input$gene, as.character(Gene.Symbol), ""))
    cano<-cano%>%mutate(size = ifelse(cano$Gene.Symbol==input$gene, 3, 1))
    
    cano[cano$Gene.Symbol==input$gene, 4] <- ifelse(abs(cano[cano$Gene.Symbol==input$gene, 2])>=1 & cano[cano$Gene.Symbol==input$gene, 3]<=10^input$p, "B",
                                                    ifelse(abs(cano[cano$Gene.Symbol==input$gene, 2])>=input$FC, "C",
                                                           ifelse(cano[cano$Gene.Symbol==input$gene, 3]<=10^input$p, "D", "A")))
    
    ggplot(cano, aes(x=glia.brain_logFC, y=-log10(P.Value), label = callit)) +
      geom_point(aes(color = threshold), alpha=0.4, size=cano$size) +
      scale_colour_manual(values = c("A"= "green", "B"="red",  "C"="blue", "D"= "black")) +
      xlab("log fold change") + ylab("-log10 p-value") +
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                         plot.title = element_text(hjust = 1),
                         legend.position = "none") +
      ggtitle("Galatro et al, 2017") +
      geom_text_repel() +
      geom_vline(aes(xintercept = input$FC), col="black", lty=3, lwd=1.0) +
      geom_vline(aes(xintercept = -(input$FC)), col="black", lty=3, lwd=1.0) +
      geom_hline(aes(yintercept = -(input$p)), col="black", lty=3, lwd=1.0) +
      scale_x_continuous(breaks=seq(-4,7,1))
  })
  
  ### If it's not present in the specified papers, why is it in the database? 
  ### (Sometimes-present larger table of where evidence is found)
  
  output$selectedgene <- DT::renderDataTable({
    req(nrow(extras[extras$Gene==input$gene,]) > 0)
    selectedgene <- extras[extras$Gene == input$gene,]
    DT::datatable(data = selectedgene, rownames= FALSE,
                  escape = FALSE,
                  selection = 'none',
                  options = list(columnDefs = list(list(
                    targets = 2,
                    render = JS(
                      "function(data, type, row, meta) {",
                      "return type === 'display' && data.length > 20 ?",
                      "'<span title=\"' + data + '\">' + data.substr(0, 20) + '...</span>' : data;",
                      "}")
                  )), dom='t',bPaginate=FALSE)
    )
  })
  
  
  #### --- IN TAB 2 --- ####
  
  ##counting the rows in the input file or text
  file_row <- reactive({
    inFile <- input$file1
    if (is.null(inFile)){
      tmp <- as.data.frame(strsplit(input$caption, " "))
      colnames(tmp) <- "Gene"
      y <- merge(allgenes, tmp)
      return (list(nrow(tmp), length(unique(y$Gene))))
    }else{
      tmp <- read.csv(inFile$datapath, header = F)
      colnames(tmp) <- "Gene"
      y <- merge(allgenes, tmp)
      return (list(nrow(tmp), length(unique(y$Gene))))
    }
  })
  
  ##textoutput of number of genes in list found in papers
  output$totalnrow <- renderText({
    if(length(file_row) < 1) {
      return(NULL)
    } else {
    inFile <- input$file1
    if (is.null(inFile)){
      return(paste0(file_row()[[2]], " (of your list of ", file_row()[[1]], ") genes are related to microglia or Alzheimer disease"))
    }else{
      paste0(file_row()[[2]], " (of your list of ", file_row()[[1]], ") genes are related to microglia or Alzheimer disease")
    }
    }
  })
  
  
  ###upload file and compare to database
  output$tab2table <- DT::renderDataTable({
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame.
    inFile <- input$file1
    
    if (is.null(inFile)) {
      x <- as.data.frame(strsplit(input$caption, " "))
      
  } else {
    
    x <- read.csv(inFile$datapath, header=F)
  }
    colnames(x) <- "Gene"
    y <- merge(allgenes, x, 
               all.y = TRUE
    )
    y <- y[,1:4]
  
    y$`Paper title`[is.na(y$`Paper title`)] <- "No known link"
    y[is.na(y)] <- ""
    
    DT::datatable(data = y, rownames= FALSE,
                  escape = FALSE,
                  selection = 'none',
                  options = list(columnDefs = list(list(
                    targets = 2,
                    render = JS(
                      "function(data, type, row, meta) {",
                      "return type === 'display' && data.length > 20 ?",
                      "'<span title=\"' + data + '\">' + data.substr(0, 20) + '...</span>' : data;",
                      "}")
                  )), dom='t',bPaginate=FALSE)
    )
  
    })
  
  # output$textcontent <- DT::renderDataTable({
  #   inText <- as.data.frame(input$caption)
  #   colnames(x) <- "Gene"
  #   y <- merge(allgenes, x, 
  #              all.y = TRUE
  #   )
  #   y <- y[,1:4]
  #   
  #   y$`Paper title`[is.na(y$`Paper title`)] <- "No known link"
  #   y[is.na(y)] <- ""
  #   
  #   DT::datatable(data = y, rownames= FALSE,
  #                 escape = FALSE,
  #                 selection = 'none',
  #                 options = list(columnDefs = list(list(
  #                   targets = 2,
  #                   render = JS(
  #                     "function(data, type, row, meta) {",
  #                     "return type === 'display' && data.length > 20 ?",
  #                     "'<span title=\"' + data + '\">' + data.substr(0, 20) + '...</span>' : data;",
  #                     "}")
  #                 )), dom='t',bPaginate=FALSE)
  #   )
  # })
  # 
  # output$tab2table <- reactive({
  #   inFile <- input$file1
  #   if (is.null(inFile)){
  #     return(output$textcontent)
  #   }else{
  #     return (output$filecontent)
  #   }
  # })

  
  }




# Run the application 
shinyApp(ui = ui, server = server)

