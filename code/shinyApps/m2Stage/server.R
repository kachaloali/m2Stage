# Install package if it does not exist yet
library(ggplot2)
if (!require("UpSetR")){
  install.packages("UpSetR")
  library(UpSetR)
}else{
  library(UpSetR)
}
if (!require(futile.logger)){
  install.packages("futile.logger")
  library(futile.logger)
}else{
  library(futile.logger)
}

options(shiny.sanitize.errors = TRUE)
flog.threshold(ERROR)

# This app will have the following main components:
# => User need to load the data that is the output of analyzers scripts. 
# the contents of this data is a list of ec numbers predicted by a given program
# of prediction such as: PRIAM, E2P2, IPRSCAN, Blasta2Go, etc.
# => This app load the data and copy it in a temprary directory named tmp.
# => The data is plotted on the dashboardBody of the app showing a circle that
# represent the numbers of the ec numbers contained in the input file.
# As and when you add a file, it is plotted in the dashboardBody and shows the 
# corresponding venn Diagram.

# Define server logic required
shinyServer(
  function(input, output) {
    mainDir <- "."
    tmpDir <- "tmp"
    # delete temp dir if it exists already
    if (dir.exists(file.path(mainDir, tmpDir))){
      unlink(file.path(mainDir, tmpDir), recursive=TRUE)
    }
    nbOfPrograms = 5
    insertedNums = c()
    insertedFileNames = c()
    insertedProgs = c()
    #%%%%%%%%%%%%%%%%%%%%%%%%%
    insertedNumsTab = c()
    insertedFileNamesTab = c()
    insertedLevels = c()
    #%%%%%%%%%%%%%%%%%%%%%%%%%
    insertedNumsHeat = c()
    insertedFileNamesHeat = c()
    insertedMesures = c()
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # This function is called when the user loads a file on the application and 
    # specify the program associated (PRIAM, E2P2, IPRSCAN, Blasta2Go, etc) to this file.
    # its main role is to store the file in a temporary directory called tmp
    # and display the venn diagram corresponding to this file
    save_file <- reactive({
      idFile = paste("filesUpload", input$file, sep = "")
      inFile = input[[paste("filesUpload", input$file, sep = "")]]
      if (is.null(inFile)){
        return(NULL)
      }
      
      if (input$file > 0 && !is.null(inFile)){
        dir.create(file.path(mainDir, tmpDir), showWarnings = FALSE)
        inFile_without_ext = tools::file_path_sans_ext(basename(inFile$name))
        num = input$file
        switch (input[[paste("programSelect", input$file, sep="")]],
            prog1 = {
              prog <- 'sp'
            },
            prog2 = {
              prog <- 'priam'
            },
            prog3 = {
              prog <- 'e2p2'
            },
            prog4 = {
              prog <- 'iprscan'
            },
            prog5 = {
              prog <- 'kaas'
            },
            prog6 = {
              prog <- 'Blast2go'
            },
            prog7 = {
              prog <- 'koala'
            },
            prog8 = {
              prog <- 'Detect'
            },
            {
              prog = 'None'
            }
        )
        
        if(prog != 'None'){
          filename <- file.path(mainDir, paste0(tmpDir, "/", inFile_without_ext, "_", prog, ".", num,".txt"))
          file.copy(inFile$datapath, filename)
          if(length(insertedNums) > 0){
            prevNum = insertedNums[length(insertedNums)]
          }else{
            prevNum = 0
          }
          insertedNums <<- c(insertedNums, num)
          insertedFileNames <<- c(insertedFileNames, inFile_without_ext)
          pos = which(insertedNums == num)
          if (prevNum != num) {
            insertedProgs <<- c(insertedProgs, prog)
          }else{
            replace(insertedProgs, pos, prog)
          }
        }
       
      if (length(insertedProgs) > 0){
          # making venn diagram
          library(VennDiagram)
          l1 <- list()
          for (i in 1:length(insertedProgs)){
            IF = insertedFileNames[i]
            IP = insertedProgs[i]
            IN = insertedNums[i]
            file <- file.path(mainDir, paste0(tmpDir, "/", IF, "_", IP, ".", IN,".txt"))
            df <- read.table(file)
            l1[paste0(toupper(insertedProgs[i]))] <- df
          }
          
          # The venn.diagram() function in the package VennDiagram only works with sizes up to 5.
          if (length(insertedProgs) > 5){
            l2 = l1[1:5]
            
            vp <- venn.diagram(l2, 
                 fill = 1:length(l2), alpha = 0.3, filename = NULL,
                 output = FALSE , compression = "lzw", category = names(l2),
                 cat.fontface = "bold", cat.fontfamily = "time", lwd = "2", 
                 main = "Diagramme de Venn de prédiction des numéros ECs", main.fontface = "bold"
            )
          }else{
            vp <- venn.diagram(l1, 
                 fill = 1:length(l1), alpha = 0.3, filename = NULL,
                 output = FALSE , compression = "lzw", category = names(l1),
                 cat.fontface = "bold", cat.fontfamily = "time", lwd = "2", 
                 main = "Diagramme de Venn de prédiction des numéros ECs", main.fontface = "bold"
            )
          }
         
          # Ploting table of ec numbers
          # list.len = sapply(l1, function(x) nrow(unique(data.frame(x))))
          # max.len = max(list.len)
          # l3 = list()
          # for (i in 1:length(insertedProgs)){
          #   data <- data.frame(unlist(l1[paste0(toupper(insertedProgs[i]))]))
          #   data.nrow = nrow(unique(data))
          #   l4 = list()
          #   for(j in 1: data.nrow){
          #     l4 <- c(l4, paste(unique(data)[[1]][j]))
          #   }
          #   l4 <- sort(unlist(l4))
          #   if(max.len - data.nrow > 0){
          #     l4 <- c(l4, rep(NA, max.len - data.nrow))
          #   }
          #   l3[paste0(toupper(insertedProgs[i]))] <- data.frame(unlist(l4))
          # }
          
          output$outPlot <- renderPlot({ grid.draw(vp) })
          output$summary <- renderPrint({ intersect_all(l1) })
          #output$ectab <- renderTable(data.frame(l3))
          unlink(file.path(mainDir, paste0("*.log")))
          #
          # making upset ploting
          library(UpSetR)
          vect <- c()
          l2 <- list()
          for (i in 1:length(insertedProgs)){
            IF = insertedFileNames[i]
            IP = insertedProgs[i]
            IN = insertedNums[i]
            file <- file.path(mainDir, paste0(tmpDir, "/", IF, "_", IP, ".", IN,".txt"))
            df <- read.table(file)
            vect <- c(vect, as.vector(df[,1]))
          }
          vect <- vect[!duplicated(vect)] # remove duplicated ec numbers
          for (i in 1:length(insertedProgs)){
            elts = l1[paste0(toupper(insertedProgs[i]))]
            elts <- as.character(unlist(elts))
            v <- c()
            for (j in 1:length(vect)){
              if (is.element(vect[j], elts) == TRUE){
                v <- c(v, 1)
              }else{
                v <- c(v, 0)
              }
            }
            l2[paste0(toupper(insertedProgs[i]))] <- data.frame(v)
          }
          df <- as.data.frame(l2)
          if (length(insertedProgs) > 1){
            output$upsetPlot <- renderPlot({
              upset(df, sets = colnames(df), 
                sets.bar.color = "#56B4E9", 
                order.by = "freq", 
                empty.intersections = "on",
                text.scale = 2
              )
            },
            width = "auto", height = "auto", res = 72)
          }else{
            output$upsetPlot <- renderPlot({ paste0("The number of set must be greater than 1!") })
          }
        }else{
          output$outPlot <- renderPlot({ paste0("")})
          output$summary <- renderPrint({ paste("Please, you have to load the data file") })
          output$ectab <- renderPrint({ paste("Please, you have to load the data file") })
        }
      }
    })
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # This function is called when the user remove a file on the application 
    # its main role is to remove the file in a temporary directory called tmp
    # and deleted the venn diagram corresponding to this file
    remove_file <- reactive({
      idFile = paste("filesUpload", input$remove, sep = "")
      inFile = input[[paste("filesUpload", input$remove, sep="")]]
      if (is.null(inFile)){
        return(NULL)
      }
      if (input$file > 0){
        inFile_without_ext = tools::file_path_sans_ext(basename(inFile$name))
        num = input$remove
        switch (input[[paste("programSelect", input$remove, sep="")]],
            prog1 = {
              prog <- 'sp'
            },
            prog2 = {
              prog <- 'priam'
            },
            prog3 = {
              prog <- 'e2p2'
            },
            prog4 = {
              prog <- 'iprscan'
            },
            prog5 = {
              prog <- 'kaas'
            },
            prog6 = {
              prog <- 'Blast2go'
            },
            prog7 = {
              prog <- 'koala'
            },
            prog8 = {
              prog <- 'Detect'
            },
            {
              prog <- 'None'
            }
        )
        if (prog != 'None'){
          unlink(file.path(tmpDir, paste0("*.", num, ".txt")))
          pos = which(insertedNums == num)
          insertedNums <<- insertedNums[-which(insertedNums == num)]
          insertedFileNames <<- insertedFileNames[-pos]
          insertedProgs <<- insertedProgs[-which(insertedProgs == prog)]
        }
        
        if (length(insertedProgs) > 0){
          # making venn diagram
          library(VennDiagram)
          l1 <- list()
          for (i in 1:length(insertedProgs)){
            IF = insertedFileNames[i]
            IP = insertedProgs[i]
            IN = insertedNums[i]
            file <- file.path(mainDir, paste0(tmpDir, "/", IF, "_", IP, ".", IN,".txt"))
            df <- read.table(file)
            l1[paste0(toupper(insertedProgs[i]))] <- df
          }
          
          # The venn.diagram() function in the package VennDiagram only works with sizes up to 5.
          if (length(insertedProgs) > 5){
            l2 = l1[1:5]
            vp <- venn.diagram(l2, 
                 fill = 1:5, alpha = 0.3, filename = NULL,
                 output = FALSE , compression = "lzw", category = names(l2),
                 cat.fontface = "bold", cat.fontfamily = "time", lwd = "2", 
                 main = "Diagramme de Venn de prédiction des numéros ECs", main.fontface = "bold"
            )
          }else{
            vp <- venn.diagram(l1, 
                 fill = 1:length(l1), alpha = 0.3, filename = NULL,
                 output = FALSE , compression = "lzw", category = names(l1),
                 cat.fontface = "bold", cat.fontfamily = "time", lwd = "2", 
                 main = "Diagramme de Venn de prédiction des numéros ECs", main.fontface = "bold"
            )
          }
          
          #
          # Ploting table of ec numbers
          # list.len = sapply(l1, function(x) nrow(unique(data.frame(x))))
          # max.len = max(list.len)
          # l3 = list()
          # for (i in 1:length(insertedProgs)){
          #   data <- data.frame(unlist(l1[paste0(toupper(insertedProgs[i]))]))
          #   data.nrow = nrow(unique(data))
          #   l4 = list()
          #   for(j in 1: data.nrow){
          #     l4 <- c(l4, paste(unique(data)[[1]][j]))
          #   }
          #   l4 <- sort(unlist(l4))
          #   if(max.len - data.nrow > 0){
          #     l4 <- c(l4, rep(NA, max.len - data.nrow))
          #   }
          #   l3[paste0(toupper(insertedProgs[i]))] <- data.frame(unlist(l4))
          # }
          output$outPlot <- renderPlot({ grid.draw(vp) })
          output$summary <- renderPrint({ intersect_all(l1) })
          #output$ectab <- renderTable(data.frame(l3))

          # making upset ploting
          library(UpSetR)
          vect <- c()
          l2 <- list()
          for (i in 1:length(insertedProgs)){
            IF = insertedFileNames[i]
            IP = insertedProgs[i]
            IN = insertedNums[i]
            file <- file.path(mainDir, paste0(tmpDir, "/", IF, "_", IP, ".", IN,".txt"))
            df <- read.table(file)
            vect <- c(vect, as.vector(df[,1]))
          }
          vect <- vect[!duplicated(vect)] # remove duplicated ec numbers
          for (i in 1:length(insertedProgs)){
            elts = l1[paste0(toupper(insertedProgs[i]))]
            elts <- as.character(unlist(elts))
            v <- c()
            for (j in 1:length(vect)){
              if (is.element(vect[j], elts) == TRUE){
                v <- c(v, 1)
              }else{
                v <- c(v, 0)
              }
            }
            l2[paste0(toupper(insertedProgs[i]))] <- data.frame(v)
          }
          df <- as.data.frame(l2)
          if (length(insertedProgs) > 1){
            output$upsetPlot <- renderPlot({
              upset(df, sets = colnames(df), 
                sets.bar.color = "#56B4E9", 
                order.by = "freq", 
                empty.intersections = "on",
                text.scale = 2
                )
            },
            width = "auto", height = "auto", res = 72)
          }else{
            output$upsetPlot <- renderPlot({ paste0("The number of set must be greater than 1!") })
          }
        }else{
          output$outPlot <- renderPlot({ paste0("")})
          output$summary <- renderPrint({ paste("Please, you have to load the data file") })
          #output$ectab <- renderPrint({ paste("Please, you have to load the data file") })
        }
      }
    })
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # This function is called when the user loads a file on the application and 
    # specify the program associated (PRIAM, E2P2, IPRSCAN, Blasta2Go, etc) to this file.
    # its main role is to store the file in a temporary directory called tmp
    # and display the venn diagram corresponding to this file
    save_fileTab <- reactive({
      idFileTab = paste("filesUploadTab", input$fileTab, sep = "")
      inFileTab = input[[paste("filesUploadTab", input$fileTab, sep = "")]]
      if (is.null(inFileTab)){
        return(NULL)
      }
      
      if (input$fileTab > 0 && !is.null(inFileTab)){
        dir.create(file.path(mainDir, tmpDir), showWarnings = FALSE)
        inFileTab_without_ext = tools::file_path_sans_ext(basename(inFileTab$name))
        numTab = input$fileTab
        switch (input[[paste("levSelectTab", input$fileTab, sep="")]],
            lev1 = {
              level <- '1'
            },
            lev2 = {
              level <- '2'
            },
            lev3 = {
              level <- '3'
            },
            lev4 = {
              level <- '4'
            },
            {
              level <- '0'
            }
        )
        
        if(level != '0'){
          filenameTab <- file.path(mainDir, paste0(tmpDir, "/", "ec_to_pick", ".", level, ".tab"))
          file.copy(inFileTab$datapath, filenameTab)
          if(length(insertedNumsTab) > 0){
            prevNumTab = insertedNumsTab[length(insertedNumsTab)]
          }else{
            prevNumTab = 0
          }
          insertedNumsTab <<- c(insertedNumsTab, numTab)
          insertedFileNamesTab <<- c(insertedFileNamesTab, inFileTab_without_ext)
          posTab = which(insertedNumsTab == numTab)
          if (prevNumTab != numTab) {
            insertedLevels <<- c(insertedLevels, level)
          }else{
            replace(insertedLevels, posTab, level)
          }
        }
        
        if (length(insertedLevels) > 0){
          df = read.table(filenameTab, header = TRUE)
          output$ectab <- renderTable(df, striped =TRUE, bordered = TRUE, hover=TRUE, width = '100%')
        }else{
          output$ectab <- renderPrint({ paste("Please, you have to load the data file") })
        }
      }
    })
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # This function is called when the user loads a file on the application and 
    # specify the program associated (PRIAM, E2P2, IPRSCAN, Blasta2Go, etc) to this file.
    # its main role is to store the file in a temporary directory called tmp
    # and display the venn diagram corresponding to this file
    remove_fileTab <- reactive({
      filenameTab <- file.path(mainDir, paste0(tmpDir, "/", "ec_to_pick.", input$removeTab, ".tab"))
      if (file.exists(filenameTab)){
        unlink(file.path(tmpDir, paste0("*.", input$removeTab, ".tab")))
      }
    })
    
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # This function is called when the user loads a file on the application and 
    # specify the program associated (PRIAM, E2P2, IPRSCAN, Blasta2Go, etc) to this file.
    # its main role is to store the file in a temporary directory called tmp
    # and display the venn diagram corresponding to this file
    save_fileHeat <- reactive({
      idFileHeat = paste("filesUploadHeat", input$fileHeat, sep = "")
      inFileHeat = input[[paste("filesUploadHeat", input$fileHeat, sep = "")]]
      
      if (is.null(inFileHeat)){
        return(NULL)
      }
      if (input$fileHeat > 0 && !is.null(inFileHeat)){
        dir.create(file.path(mainDir, tmpDir), showWarnings = FALSE)
        inFileHeat_without_ext = tools::file_path_sans_ext(basename(inFileHeat$name))
        numHeat = input$fileHeat
        switch (input[[paste("mesSelectHeat", input$fileHeat, sep="")]],
                sens = {
                  type <- 'sens'
                },
                prec = {
                  type <- 'prec'
                },
                f1score = {
                  type <- 'f1score'
                },
                {
                  type <- 'none'
                }
        )
        
        if(type != 'none'){
          filenameHeat <- file.path(mainDir, paste0(tmpDir, "/", "out.heatmap", ".", type, ".tab"))
          file.copy(inFileHeat$datapath, filenameHeat)
          
          if(length(insertedNumsHeat) > 0){
            prevNumHeat = insertedNumsHeat[length(insertedNumsHeat)]
          }else{
            prevNumHeat = 0
          }
          
          insertedNumsHeat <<- c(insertedNumsHeat, numHeat)
          insertedFileNamesHeat <<- c(insertedFileNamesHeat, inFileHeat_without_ext)
          posHeat = which(insertedNumsHeat == numHeat)
          if (prevNumHeat != numHeat) {
            insertedMesures <<- c(insertedMesures, type)
          }else{
            replace(insertedMesures, posHeat, type)
          }
        }
        
        if (length(insertedFileNamesHeat) > 0){
          df = read.table(filenameHeat, header = TRUE, sep = "\t", dec = ",")
          column1 = df$EC
          names1 = c()
          for (i in 1:nrow(df)){
            names1 <- c(names1, c(trimws(column1[i])))
          }
          row.names(df) = names1
          for (j in 2:7){
            df[,j] = as.numeric(df[,j])
          }
          df <- df[-1]
          df = as.matrix(df)
          x = colnames(df)
          y = rownames(df)
          library(RColorBrewer)
          library(ComplexHeatmap)
          library(iheatmapr)
          library(d3heatmap)
          output$heatmapPlot <- renderPlot({ 
            main_heatmap(df, name = type, colors =colorRampPalette(brewer.pal(8, "Reds"))(10)) %>%
              add_col_labels(ticktext = x, font = list(size = 10)) %>%
              add_row_labels(size = 0.3,font = list(size = 6)) %>%
              add_row_clustering() %>%
              add_col_clustering()
            # d3heatmap(df, symm=FALSE, scale = "none", col = colorRampPalette(brewer.pal(8, "Reds"))(10))
            })
        }else{
          output$heatmapPlot <- renderPrint({ paste("Please, you have to load the data file") })
        }
      }
    })
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # This function is called when the user loads a file on the application and 
    # specify the program associated (PRIAM, E2P2, IPRSCAN, Blasta2Go, etc) to this file.
    # its main role is to store the file in a temporary directory called tmp
    # and display the venn diagram corresponding to this file
    remove_fileHeat <- reactive({
      filenameHeat <- file.path(mainDir, paste0(tmpDir, "/", "out.heatmap.", input$removeHeat, ".tab"))
      if (file.exists(filenameHeat)){
        unlink(file.path(tmpDir, paste0("*.", input$removeHeat, ".tab")))
      }
    })
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    intersect_all <- function(l){
      writeLines("The common ec numbers of all programs are:")
      if (length(l) > 1){
        partitions = get.venn.partitions(l)
        dic = vector(mode = "list", length = 0)
        for (i in 1:nrow(partitions)){ 
          lab <- partitions[i,][,length(partitions)-2]
          values = partitions[i,][,length(partitions)-1]
          dic[lab] = values
          #writeLines(lab)
          #writeLines(paste0(unlist(values, use.names=FALSE)))
          #writeLines('\n')
        }
        dic
      }else{
        dic = vector(mode = "list", length = 0)
        dic[insertedProgs[1]] = l
        dic
        #writeLines(paste0(unlist(l, use.names=FALSE)))
      }
    }
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # It is an event that is called continuously. Whenever the user clicks on the buton "add file"
    # a file loading field is added. each field column is identified by a unique ID, this will allow 
    # to delete exactly the desired field
    cpt = 0
    observeEvent(input$addSelect, {
      cpt <<- cpt + 1
      insertUI(
        selector = '#placeholder',
        ## wrap element in a div with id for ease of removal
        ui = span(
          fluidRow(
            column(6, 
             fileInput(
              paste("filesUpload", cpt, sep = ""), 
              "Select file", 
              accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv'), 
              width = "100%"),
              onchange = paste0("Shiny.onInputChange('file', ", cpt, ")")
            ),
            column(4, 
             selectInput(
               paste("programSelect", cpt, sep = ""), 
               "Programs",
               list(
                 "None"         = "prog0",
                 "Kaas"         = "prog5",
                 "E2P2"         = "prog3",
                 "Priam"        = "prog2",
                 "Koala"        = "prog7",
                 "Detect"       = "prog8",
                 "IprScan"      = "prog4",
                 "Blast2go"     = "prog6",
                 "SwissProt"    = "prog1"
                 ),
               width = "95%"
             )
            ),
            column(1, 
             actionButton(
               paste("removeFactor", cpt, sep = ""), 
               label = "",
               width = "60%",
               icon = icon("times", class = NULL, lib = "font-awesome"),
               onclick = paste0("Shiny.onInputChange('remove', ", cpt, ")")
             )
            )
          ), 
          id = paste0('insert', cpt)
        )
      )
    })
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # It is an event that is called continuously. Whenever the user clicks on the buton "add file"
    # a file loading field is added. each field column is identified by a unique ID, this will allow 
    # to delete exactly the desired field
    cptTab = 0
    observeEvent(input$addSelectTab, {
      cptTab <<- cptTab + 1
      insertUI(
        selector = '#placeholderTab',
        ## wrap element in a div with id for ease of removal
        ui = span(
          column(1, 
           fileInput(
             paste("filesUploadTab", cptTab, sep = ""), 
             "Select file", 
             accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv'), 
             width = "100%"),
             onchange = paste0("Shiny.onInputChange('fileTab', ", cptTab, ")"
           )
          ),
          column(1, 
           selectInput(
             paste("levSelectTab", cptTab, sep = ""), 
             "Level",
             list(
               "0"         = "lev0",
               "1"         = "lev1",
               "2"         = "lev2",
               "3"         = "lev3",
               "4"         = "lev4"
             ),
             width = "100%"
           )
          ),
          column(1, 
           actionButton(
             paste("removeFactorTab", cptTab, sep = ""), 
             label = "",
             width = "50%",
             icon = icon("times", class = NULL, lib = "font-awesome"),
             onclick = paste0("Shiny.onInputChange('removeTab', ", cptTab, ")")
           )
          ), 
          id = paste0('insertTab', cptTab)
        )
      )
    })
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # It is an event that is called continuously. Whenever the user clicks on the buton "add file"
    # a file loading field is added. each field column is identified by a unique ID, this will allow 
    # to delete exactly the desired field
    cptHeat = 0
    observeEvent(input$addSelectHeat, {
      cptHeat <<- cptHeat + 1
      insertUI(
        selector = '#placeholderHeat',
        ## wrap element in a div with id for ease of removal
        ui = span(
          column(1, 
                 fileInput(
                   paste("filesUploadHeat", cptHeat, sep = ""), 
                   "Select file", 
                   accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv'), 
                   width = "100%"),
                 onchange = paste0("Shiny.onInputChange('fileHeat', ", cptHeat, ")"
                 )
          ),
          column(1, 
                 selectInput(
                   paste("mesSelectHeat", cptHeat, sep = ""), 
                   "Mesures",
                   list(
                     "None"         = "none",
                     "Sens"         = "sens",
                     "Prec"         = "prec",
                     "F1-SCore"     = "f1score"
                   ),
                   width = "100%"
                 )
          ),
          column(1, 
                 actionButton(
                   paste("removeFactorHeat", cptHeat, sep = ""), 
                   label = "",
                   width = "50%",
                   icon = icon("times", class = NULL, lib = "font-awesome"),
                   onclick = paste0("Shiny.onInputChange('removeHeat', ", cptHeat, ")")
                 )
          ), 
          id = paste0('insertHeat', cptHeat)
        )
      )
    })
    
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # for keeping track of elements inserted and not yet removed
    inserted <- c()
    # Each time a file is loaded, the corresponding field is stored
    observeEvent(input[[paste("filesUpload", input$file, sep = "")]], {
      idFile = paste("filesUpload", input$file, sep = "")
      inFile = input[[paste("filesUpload", input$file, sep = "")]]
      if (is.null(inFile)){
        return(NULL)
      }
      inserted <<- c(inserted, cpt)
    })
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Event that listens to the loading of the file and the selection of 
    # the associated program. Whenever this event occurs, the save_file() 
    # function is called to see above
    observeEvent(input[[paste("programSelect", input$file, sep = "")]], {
      idFile = paste("filesUpload", input$file, sep = "")
      inFile = input[[paste("filesUpload", input$file, sep = "")]]
      if (is.null(inFile)){
        return(NULL)
      }else{
        save_file()
      }
    })
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Event that listens to the removing of the loading file field. Whenever 
    # this event occurs, the remove_file() function is called to see above
    observeEvent(input$remove, {
      # removing of all the files created with the select id input$remove
      remove_file()
      j <- input$remove
      removeUI(
        # select the appropriate div id and remove it
        selector = paste0('#insert', j)
      )
      inserted <<- inserted[-length(inserted)]
    })
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # for keeping track of elements inserted and not yet removed
    insertedTab <- c()
    # Each time a file is loaded, the corresponding field is stored
    observeEvent(input[[paste("filesUploadTab", input$fileTab, sep = "")]], {
      idFileTab = paste("filesUploadTab", input$fileTab, sep = "")
      inFileTab = input[[paste("filesUploadTab", input$fileTab, sep = "")]]
      if (is.null(inFileTab)){
        return(NULL)
      }
      insertedTab <<- c(insertedTab, cptTab)
    })
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Event that listens to the loading of the file and the selection of 
    # the associated level. Whenever this event occurs, the save_fileTab() 
    # function is called to see above
    observeEvent(input[[paste("levSelectTab", input$fileTab, sep = "")]], {
      idFileTab = paste("filesUploadTab", input$fileTab, sep = "")
      inFileTab = input[[paste("filesUploadTab", input$fileTab, sep = "")]]
      if (is.null(inFileTab)){
        return(NULL)
      }else{
        save_fileTab()
      }
    })
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    observeEvent(input$removeTab, {
      # removing of all the files created with the select id input$removeTab
      remove_fileTab()
      k <- input$removeTab
      removeUI(
        # select the appropriate div id and remove it
        selector = paste0('#insertTab', k)
      )
      insertedTab <<- insertedTab[-length(insertedTab)]
    })
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Event that listens to the entering of an EC number. Whenever this event occurs, 
    # all the informations about the entered EC number are searched from the file which
    # contain this EC number and printed in a DataTable. The files are to be download with 
    # the button «Select a file» that is in the tabPanel of «ECs Tab» space.
    observeEvent(input$inputECnumber, {
      ec = unlist(strsplit(input$inputECnumber, ".", fixed = TRUE))
      if (length(ec) > 0 & length(ec) < 5){
        filenameTab <- file.path(mainDir, paste0(tmpDir, "/", "ec_to_pick.", length(ec), ".tab"))
        normal_ec = paste(ec, collapse = ".")
        if (file.exists(filenameTab)){
          df1 = read.table(filenameTab, header = TRUE, dec = ",", numerals = c("no.loss"))
          df2 <- data.frame(EC = character(),
                            PROG = character(),
                            TP = character(),
                            FP = character(),
                            FN = character(),
                            SENS = character(),
                            PREC = character(),
                            F1SCORE = character(),
                            N = character(),
                            ZSCORE = character(), 
                            stringsAsFactors=FALSE) 
          df2 <- setNames(df2, names(df1))
          for (l in 1:nrow(df1)){
            ec_tab = as.character(unlist(df1[l,][1]))
            if (ec_tab == normal_ec){
              df2 <- rbind(df2, df1[l,], stringsAsFactors = FALSE)
            }
          }
          if (nrow(df2) == 0){
            df2 <- data.frame("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")
            names(df2) <- c("EC", "PROG","TP", "FP", "FN", "SENS", "PREC", "F1SCORE", "N", "ZSCORE")
            output$ectab <- renderTable(df2, striped =TRUE, bordered=TRUE, hover=TRUE, width='100%')
          }else{
            output$ectab <- renderTable(df2, striped =TRUE, bordered=TRUE, hover=TRUE, width='100%')
          }
        }else{
          output$ectab <- renderPrint({ 
            paste("Please, you have to load the file that contain ECs numbers on level ",length(ec)) 
          })
        }
      }else if(length(insertedNumsTab) > 0){
        if(length(ec) > 0){
          output$ectab <- renderPrint({ 
            paste("Please, you have to load the file that contain ECs numbers on level ",length(ec)) 
          })
        }else{
          df2 <- data.frame("NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA")
          names(df2) <- c("EC", "PROG","TP", "FP", "FN", "SENS", "PREC", "F1SCORE", "N", "ZSCORE")
          output$ectab <- renderTable(df2, striped =TRUE, bordered=TRUE, hover=TRUE, width='100%')
        }
      }
    })
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # for keeping track of elements inserted and not yet removed
    insertedHeat <- c()
    # Each time a file is loaded, the corresponding field is stored
    observeEvent(input[[paste("filesUploadHeat", input$fileHeat, sep = "")]], {
      idFileHeat = paste("filesUploadHeat", input$fileHeat, sep = "")
      inFileHeat = input[[paste("filesUploadHeat", input$fileHeat, sep = "")]]
      if (is.null(inFileHeat)){
        return(NULL)
      }
      insertedHeat <<- c(insertedHeat, cptHeat)
    })
    
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Event that listens to the loading of the file and the selection of 
    # the associated level. Whenever this event occurs, the save_fileTab() 
    # function is called to see above
    observeEvent(input[[paste("mesSelectHeat", input$fileHeat, sep = "")]], {
      idFileHeat = paste("filesUploadHeat", input$fileHeat, sep = "")
      inFileHeat = input[[paste("filesUploadHeat", input$fileHeat, sep = "")]]
      if (is.null(inFileHeat)){
        return(NULL)
      }else{
        save_fileHeat()
      }
    })
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    observeEvent(input$removeHeat, {
      # removing of all the files created with the select id input$removeTab
      #remove_fileHeat()
      k <- input$removeHeat
      removeUI(
        # select the appropriate div id and remove it
        selector = paste0('#insertHeat', k)
      )
      insertedHeat <<- insertedHeat[-length(insertedHeat)]
    })
  }
)