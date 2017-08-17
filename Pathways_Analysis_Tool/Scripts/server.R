options(shiny.maxRequestSize=90*1024^2)

shinyServer(function(input,output){
  GE_FILE <- NULL
  CNV_FILE <- NULL
  GE_Data <- NULL
  CNV_Data <- NULL
  
  output$SampsInfo <- renderDataTable({
    input$goButton
    inFile1 <- input$file1
    if (is.null(inFile1)){
          return(NULL)
    }else{
      x <- read.delim(inFile1$datapath, header=T, sep="\t",check.names = F,stringsAsFactors = F)
      x2 <- as.data.frame(table(x$Sample_Type))
      rownames(x2) <- x2[,1]
      x2[,1] <- NULL
      x2 <- as.data.frame(t(x2))
      withProgress(message = 'TEST', value = 0, {
        n <- 10
        for (i in 1:n) {
          incProgress(1/n, detail = paste("Doing part", i))
          Sys.sleep(0.1)
        }
      })
      total <- 20
      # create progress bar
      # pb <- txtProgressBar(min = 0, max = total, style = 3)
      # for(i in 1:total){
      #   Sys.sleep(0.1)
      #   # update progress bar
      #   setTxtProgressBar(pb, i)
      # }
      # close(pb)
      
      # pb <- tkProgressBar(title = "progress bar", min = 0,
      #                     max = total, width = 300)
      # 
      # for(i in 1:total){
      #   Sys.sleep(0.1)
      #   setTkProgressBar(pb, i, label=paste( round(i/total*100, 0),
      #                                        "% done"))
      # }
      # close(pb)
      
      return(x2)
    }
  })
  
  output$Process_table <- renderDataTable(as.data.frame(rbind(c("Data Input","pending"),c("InFlo Preprocessing","pending"),c("InFlo Running","Pending"))))
  
  
  output$ExpInfo <- renderDataTable({
    input$goButton
    inFile2 <- input$file2
    GE_FILE <<- inFile2$datapath
    print(inFile2$name)
    if (is.null(inFile2)){
      return(NULL)
    }else{
      Data2 <- Data_Read(inFile2$datapath)
      GE_Data <<- Data2
      return(GE_Data)
      withProgress(message = 'TEST', value = 0, {
        n <- 10
        
        for (i in 1:n) {
          incProgress(1/n, detail = paste("Doing part", i))
          Sys.sleep(0.1)
        }
      })
    }
  },outputArgs = c("GE_FILE","GE_Data"))
  
  output$CNVInfo <- renderDataTable({
    input$goButton
    inFile3 <- input$file3
    CNV_FILE <<- inFile3$datapath
    print(inFile3$name)
    if (is.null(inFile3)){
      return(NULL)
    }else{
      Data2 <- Data_Read(inFile3$datapath)
      CNV_Data <<- Data2
      return(CNV_Data)
      withProgress(message = 'TEST', value = 0, {
        n <- 10
        for (i in 1:n) {
          incProgress(1/n, detail = paste("Doing part", i))
          Sys.sleep(0.1)
        }
      })
    }
  },outputArgs = c("CNV_FILE","CNV_Data"))
  
  output$PathInfo <- renderDataTable({
    input$goButton
    inFile4 <- input$file4
    if (is.null(inFile4)){
      return(NULL)
    }else{
      x <- read.delim(inFile4$datapath, header=T, sep="\t",check.names = F,stringsAsFactors = F)
      #x2 <- as.data.frame(table(x$Sample_Type))
      #rownames(x2) <- x2[,1]
      #x2[,1] <- NULL
      #x2 <- as.data.frame(t(x2))
      withProgress(message = 'TEST', value = 0, {
        n <- 10
        
        for (i in 1:n) {
          incProgress(1/n, detail = paste("Doing part", i))
          Sys.sleep(0.1)
        }
      })
      return(x)
    }
  })
  # Initial_chk <- File_chk()
  # GE_FILE <- input$file2$datapath
  # CNV_FILE <- input$file3$datapath
  # GE_Data <- dataTableOutput('ExpInfo')
  # CNV_Data <- dataTableOutput('CNVInfo')
  # print(CNV_FILE)
  # print(GE_FILE)
  #Dir_create(Initial_chk)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  #will need to change conditionals in most functions so they do not rely on inputting the 
  # files from a computer to run 
  
  
  totNodes <- NULL
  totEdges <- NULL
  nodes <- NULL
  edges <- NULL
  inflo <- NULL
  newNodes <- NULL
  palette <- NULL
  rampbreaks <- NULL
  prevN1 <- NULL
  prevE1 <- NULL
  nextN1 <- NULL
  nextE1 <- NULL
  samps <- NULL
  mean <- NULL
  sd <- NULL
  survival <- NULL
  surviveG1 <- NULL
  surviveG2 <- NULL
  level <- NULL



  
  #reads in files from the computer, orders totNodes so id is the first column and totEdges so to and from are
  # first two columns.  creates all the columns that will be manipulated later and assigns placeholder values.
  # creates palette based on maximum and minimum activity in the table.
  observe({
    if(is.null(input$file5) || is.null(input$file6) || is.null(input$file7) || is.null(input$file8)){
      #
    }else{
      nodes1 <- read.csv(file=input$file5$datapath, header=T, sep="\t")
      totNodes <- nodes1[,c(1,2,3,15:ncol(nodes1))]
      colnames(totNodes)[1] <- "Parent"
      colnames(totNodes)[3] <- "id"
      totNodes <- totNodes[, c(3,1,2, 4:(ncol(totNodes)))]
      min <- min(totNodes[, c(4:ncol(totNodes))])

      max <- max(totNodes[, c(4:ncol(totNodes))])
      col1 <- colorRampPalette(colors = c("darkgreen", "white"), space="Lab")(nrow(totNodes/2))
      col2 <- colorRampPalette(colors = c("white", "darkred"), space="Lab")(nrow(totNodes/2))
      palette <<- c(col1, col2)
      rb1 <- seq(min, 0, length.out=nrow(totNodes/2)+1)
      rb2 <- seq(0, max, length.out=nrow(totNodes/2)+1)[-1]
      rampbreaks <<- c(rb1, rb2)
      
      totNodes$means <- rowMeans(totNodes[,4:ncol(totNodes)])
      totNodes$SDs <- apply(totNodes[,4:(ncol(totNodes-1))], 1, sd)
      totNodes$avg1 <- 0
      totNodes$avg2 <- 0
      totNodes$Difference <- 0
      totNodes$p <- 1
      mean <<- mean(totNodes[['means']])
      sd <<- sd(totNodes$means)
      level <<- 1.5
      parTypes <- unique(totNodes$Parent_Type)
      shapes <- c("dot", "square", "diamond", "triangle", "triangleDown", "star")

      totNodes <- totNodes[order(totNodes$Parent_Type),]
      
      #assign aesthetics
      totNodes$shape <- shapes[totNodes$Parent_Type]
      totNodes$title <- totNodes$Parent
      totNodes$color.background <- palette[as.numeric(cut((totNodes$means), breaks=rampbreaks))]

      totNodes$color.highlight.background <- "yellow"
      totNodes$color.highlight.border <- "black"
      totNodes$color.border <- "black"
      totNodes$physics <- T
      totNodes$hidden <- F
      totNodes$label <- totNodes$Parent
      totNodes <<- totNodes
      nodes <<- totNodes[(abs(totNodes$SDs)>1.5), ]
      edges1 <<- read.csv(file=input$file6$datapath, header=T, sep="\t")

      colors <- c("darkviolet", "darkgreen", "darkblue", "gold", "lightseagreen",  "deeppink", "indianred", "black")
      totEdges <- edges1[, c(5,2,1,4,7)]
      totEdges <- totEdges[(!duplicated(totEdges)), ]
      totEdges$id <- 1:nrow(totEdges)
      colnames(totEdges)[1] <- "to"
      colnames(totEdges)[2] <- "from"
      totEdges <- totEdges[order(totEdges$Interaction),]
      totEdges$hidden <- F
      totEdges$arrows <- "to"
      totEdges$color <- (colors[totEdges$Interaction])
      totEdges <<- totEdges
      edges <- totEdges[(totEdges$to %in% nodes$id), ]
      edges <- edges[(edges$from %in% nodes$id), ]
      edges <<- edges


      
      inflo <<- read.csv(file=input$file7$datapath, header=T, sep="\t") 
      survival <<- read.csv(file=input$file8$datapath, header=T, sep="\t")
    }
  })
  
  #allows user to choose which category to group nodes by, only works for current survival format
  output$Gr <- renderUI({
    if(is.null(input$file5) || is.null(input$file6) || is.null(input$file7) || is.null(input$file8)){
      #
    }else{
      names <- colnames(survival)[c(3,4,5,8,10,11)]
      selectInput("group", "Group By", names, selected=names[[2]])
    }
  })
  
  #gives options to select for group1 of samples (group2 would be samples with nonselected attributes)
  output$Gr1 <- renderUI({
    if(is.null(input$file5) || is.null(input$file6) || is.null(input$file7) || is.null(input$file8) || is.null(input$group)){
      #
    }else{
      options <- unique(survival[[input$group]])
      
      options <- sort(options[!is.na(options)])

      checkboxGroupInput("g1", "Group 1", choices=options, selected=options[[1]])
    }
  })
  
  
  #allows nodes to be selected by Parent column in select bar
  output$selectNames <- renderUI({
    if(is.null(input$file5) || is.null(input$file6) || is.null(input$file7) || is.null(input$file8)){
      #
    }else{
      names <- totNodes$Parent
      selectInput("selectNodes", "Add Nodes to Network", names, multiple=T)
    }
  })



  
  #allows sample names to be selected from InFlo file
  output$selectSamp <- renderUI({
    if(is.null(input$file5) || is.null(input$file6) || is.null(input$file7) || is.null(input$file8)){
      #
    }else{
      tochoice <- character()
      for(i in 15:ncol(inflo)){
        tochoice <- append(tochoice, as.character(colnames(inflo)[[i]]))
      }
      selectInput("selectSample", "Samples", tochoice, multiple = T)
    }
  })
  
  #when most deviant is selected, shows SD box
  output$spread1 <- renderUI({
    if(input$type == "Most Deviant"){
      numericInput(inputId="spread", label="Show Nodes with SD Greater Than", value=1.5, step=.25, min=0, max=15)
    }
    else{
      return(NULL)
    }
  })
  
  #when most disparate is selected, shows p-value box
  output$dis1 <- renderUI({
    if(input$type == "Most Disparate"){
      div(style="display: inline-block;vertical-align:center; width: 250px;", numericInput(inputId="P", label="Nodes with P-Value Less Than", value=.05, step=.01, min=0))
    }else{
      return(NULL)
    }
  })
  
  #when most disparate is selected, shows difference box
  output$dis2 <- renderUI({
    if(input$type == "Most Disparate"){
      div(style="display: inline-block;vertical-align:center; width: 260px;", numericInput(inputId="disparate", label="Nodes with Difference Greater Than", value=1, step=.5, min=0))
    }else{
      return(NULL)
    }
  })
  
  #makes the original network - only changed through proxy so the whole network isn't redrawn with every change
  output$network <- renderVisNetwork({
    
    if(is.null(input$file5) || is.null(input$file6) || is.null(input$file7) || is.null(input$file8)){
      #
    }else{
      #allows legend to work for various inputs (up to 6 types of nodes and 8 types of edges)
      intTypes <- unique(totEdges$Interaction)
      colors <- c("darkviolet", "darkgreen", "darkblue", "gold", "lightseagreen",  "deeppink", "indianred", "black")
      parTypes <- unique(totNodes$Parent_Type)
      shapes <- c("dot", "square", "diamond", "triangle", "triangleDown", "star")

      lnodes <- data.frame(label = c(as.character(parTypes), "High Activity", "Low Activity"), 
                 shape=c(shapes[1:length(parTypes)], "dot", "dot"), 
                 color.background=c(rep(c("slategrey"), times=length(parTypes)), "darkred", "darkgreen"), 
                 color.border=c(rep(c("black"), times=length(parTypes)), "darkred", "darkgreen"))
      ledges <- data.frame(label = intTypes, color = colors[1:length(intTypes)])

      nodes <<- nodes
      edges <<- edges

      print("creating network")
      
      ##### Switch for legends ###########
      # network <<- visNetwork(nodes, edges, width = "100%",height = "100%") %>% visIgraphLayout() %>% visLegend(addNodes=lnodes, addEdges=ledges, width=.2, 
      #           useGroup=F, stepY=70, stepX=120, ncol=2) %>% visOptions(highlightNearest=F) %>% visInteraction(
      #           multiselect=T, navigationButtons=T, hideEdgesOnDrag=T)  %>% visPhysics(stabilization=F) %>% visLayout(improvedLayout=T)
      # 
      network <<- visNetwork(nodes, edges, width = "100%",height = "100%") %>% visIgraphLayout() %>% visOptions(highlightNearest=F) %>% visInteraction(
                                                                                                                 multiselect=T, navigationButtons=T, hideEdgesOnDrag=T)  %>% visPhysics(stabilization=F) %>% visLayout(improvedLayout=T)

    }
  })
  
  
 

  #create a timer for 250 ms, used for remembering selected nodes
  autoInvalidate <- reactiveTimer(250)
  
  #Set selected nodes as newNodes every 250 ms
  observe({
    autoInvalidate()
    if(is.null(input$file5) || is.null(input$file6) || is.null(input$file7) || is.null(input$file8)){
      #
    }else{
      visNetworkProxy('network') %>% visGetSelectedNodes()
      newNodes <<- input$network_selectedNodes
    }
  })
  
  

  #create new network of only selected nodes
  observeEvent(input$refreshNet, {
    if(!is.null(newNodes) && !is.null(nodes)){
      print("creating new network")
      prevN1 <<- nodes
      prevE1 <<- edges
      nextN1 <<- NULL
      nextE1 <<- NULL
      

      #hides nodes not in new network by checking if each node is in newNodes (is selected)
      nodes$hidden <- apply(nodes, 1, function(i){
        return(!(as.numeric(i['id']) %in% newNodes))
      })
      

      #hides edges not in new network by checking if both to and from ids for each edge are in newNodes 
      edges$hidden <- apply(edges, 1, function(i){
        idTo <- as.numeric(i['to'])
        idFrom <- as.numeric(i['from'])
        return( !((idTo %in% newNodes) && (idFrom %in% newNodes)))
      })
      nodes <<- nodes
      edges <<- edges
      
      print("ending new network")
      visNetworkProxy('network') %>% visUpdateNodes(nodes) %>% visUpdateEdges(edges)
    }
  })
  
  #allows removal of selected nodes and their connected edges
  observeEvent(input$remove, {
    if(!is.null(newNodes)){
      prevN1 <<- nodes
      prevE1 <<- edges
      nextN1 <<- NULL
      nextE1 <<- NULL
      
      rmEdges <- edges[(edges[['to']] %in% newNodes || edges[['from']] %in% newNodes), ]
      nodes <- nodes[(nodes[['id']] %in% newNodes) == F, ]
      edges <- edges[(edges[['id']] %in% rmEdges) == F, ]
      nodes <<- nodes
      edges <<- edges
      visNetworkProxy('network') %>% visRemoveNodes(newNodes) %>% visUpdateEdges(rmEdges)
    }
  })

  
  
  #selects nodes chosen from selectbar and puts them in network
  observe({
    if(is.null(input$file5) || is.null(input$file6) || is.null(input$file7) || is.null(input$file8)){
      #
    }else{
      prevN1 <<- nodes
      prevE1 <<- edges
      nextN1 <<- NULL
      nextE1 <<- NULL
      nodesAdd <- input$selectNodes
      toAdd <- totNodes[totNodes$Parent %in% nodesAdd,]
      nodes <- rbind(nodes, toAdd)
      nodes <- nodes[!duplicated(nodes),]
      nodes$color.background <- palette[as.numeric(cut(as.numeric(nodes$means),breaks=rampbreaks))]
      edges <- totEdges[(totEdges$to %in% nodes$id), ]
      edges <- edges[(edges$from %in% nodes$id), ]
      

      hiding = logical()
      hideTo <- T
      hideFrom <- T

      for(i in 1:nrow(edges)){
 
        if(as.numeric(edges[i,][['to']]) %in% nodes$id && !nodes[which(nodes$id==edges[i,][['to']]), ][['hidden']]){
          hideTo <- F
        }
        if(as.numeric(edges[i,][['from']]) %in% nodes$id && !nodes[which(nodes$id==edges[i,][['from']]), ][['hidden']]){
          hideFrom <- F
        }
        
        #add F only if both hideTo and hideFrom are F, otherwise add T
        hiding <- append(hiding, !(!hideTo && !hideFrom))
        
        hideTo <- T
        hideFrom <- T
      }
      
      edges$hidden <- hiding


      edges <<- edges
      nodes <<- nodes
    
      
      visNetworkProxy('network') %>% visUpdateNodes(nodes) %>% visUpdateEdges(edges)
    }
  })

  
  #change coloring depending on sample selection, needs updated to take most disparate into account
  observe({
    if(is.null(input$file5) || is.null(input$file6) || is.null(input$file7) || is.null(input$file8) || input$type=='Most Disparate'){
      #
    }else{
      print("creating new coloring")
      prevN1 <<- nodes
      prevE1 <<- edges
      nextN1 <<- NULL
      nextE1 <<- NULL
      nodes1 <- nodes
      if(input$allSamps || is.null(input$selectSample)){
        samps <<- colnames(totNodes)[4:(ncol(totNodes)-15)]
      }else{
        samps <<- input$selectSample
      }
      
      
      totNodes$means <<- rowMeans(subset(totNodes, select=samps))

      mean <<- mean(totNodes[['means']])
      sd <<- sd(totNodes$means)
      if(!is.null(input$spread) && !is.na(input$spread)  && input$spread>=0){
        level <<- input$spread
      }
      if(length(samps)>1){
        totNodes$SDs <- apply(subset(totNodes, select=samps), 1, sd)
        nodes1 <- totNodes[(abs(totNodes$SDs)>level), ]
        ids <- totNodes[(abs(totNodes$SDs)<=level), ][['id']]
      }else{
        nodes1 <- totNodes[(abs(totNodes$means)>4*sd+mean), ]
        ids <- totNodes[(abs(totNodes$means)<=4*sd+mean), ][['id']]
      }
      
      edges <- totEdges[(totEdges$to %in% nodes1$id), ]
      edges <- edges[(edges$from %in% nodes1$id), ]
      edgeId <- edges[((totEdges$to %in% nodes1$id)==F && (totEdges$from %in% nodes1$id)==F), ][['id']]

      isolate({
        phys <- input$physics
      })
      
      if(nrow(nodes1)>0){
        nodes1$color.background <- palette[as.numeric(cut(as.numeric(nodes1$means),breaks=rampbreaks))]
        nodes1$physics <- phys == 'On'
      }
      edges <<- edges
      nodes <<- nodes1

      
      print("ending new coloring")
      visNetworkProxy('network')  %>% visRemoveNodes(ids) %>% visUpdateNodes(nodes) %>% visRemoveEdges(edgeId) %>% visUpdateEdges(edges)
    }
  })
  
  

  
  #when button is pressed, shows entire original network that fits selected parameters
  observeEvent(input$revert, {
    if(is.null(nodes)){
      #
    }else{
      print("original")
      prevN1 <<- nodes
      prevE1 <<- edges

      nextN1 <<- NULL
      nextE1 <<- NULL
      
      isolate({
        phys <- input$physics
        type <- input$type

        if(!is.null(input$spread) && !is.na(input$spread) && input$spread<=15 && input$spread>=0){
          level <<- input$spread
        }
        group <- input$group
        g1 <- input$g1
      })
      
      if(is.null(samps)){
        samps = colnames(inflo)[[15]]
      }
      totNodes$means <<- rowMeans(subset(totNodes, select=samps))
      totNodes$SDs <<- apply(subset(totNodes, select=samps), 1, sd)
      
      mean <<- mean(totNodes[['means']])
      sd <<- sd(totNodes$means)
      if(type == 'Most Deviant'){

        if(length(samps)>1){
          totNodes$SDs <- apply(subset(totNodes, select=samps), 1, sd)
          nodes1 <- totNodes[(abs(totNodes$SDs)>level), ]
          ids <- totNodes[(abs(totNodes$SDs)<=level), ][['id']]
        }else{
          nodes1 <- totNodes[(abs(totNodes$means)>4*sd+mean), ]
          ids <- totNodes[(abs(totNodes$means)<=4*sd+mean), ][['id']]
        }
        edges <- totEdges[(totEdges$to %in% nodes1$id), ]
        edges <- edges[(edges$from %in% nodes1$id), ]
        edgeId <- edges[((totEdges$to %in% nodes1$id)==F && (totEdges$from %in% nodes1$id)==F), ][['id']]
      }else{
        survival <- survival[!is.na(survival[[input$group]]), ]
        group1Samps <- as.character(survival$Sample_title[(survival[[group]] %in% g1)])
        
        group2Samps <- as.character(survival$Sample_title[!(survival[[group]] %in% g1)])

        if(length(group1Samps)>1 && length(group2Samps)>1){
          totNodes$avg1 <- rowMeans(totNodes[,group1Samps])
          totNodes$avg2 <- rowMeans(totNodes[,group2Samps])
         
        }else if(length(group1Samps)==1 && length(group2Samps)>1){
          totNodes$avg1 <- (totNodes[,group1Samps])
          totNodes$avg2 <- rowMeans(totNodes[,group2Samps])
        }else if(length(group2Samps)==1 && length(group1Samps)>1){
          totNodes$avg1 <- rowMeans(totNodes[,group1Samps])
          totNodes$avg2 <- (totNodes[,group2Samps])
        }else{
          totNodes$avg1 <- totNodes$means
          totNodes$avg2 <- 0
        }
        
        totNodes$Difference <- (totNodes$avg1-totNodes$avg2)

        #catches error thrown when trying to do a t test between two groups with the same values
        totNodes$p <- apply(totNodes[,4:(ncol(totNodes)-15)], 1, function(i){
          obj <- try(t.test(i[group1Samps], i[group2Samps]), silent=T)
          if(is(obj, "try-error")){
            return(1)
          }else{
            return(obj$p.value)
          }
        })

        totNodes <- totNodes[order(-(abs(totNodes$Difference))),]
        

        toAddNodes <- totNodes[totNodes$p < input$P, ]
        toAddNodes <- toAddNodes[toAddNodes$Difference > input$disparate, ]
        ids <- nodes$id[(nodes$id %in% toAddNodes$id)==F]
        toAddEdges <- totEdges[(totEdges$to %in% toAddNodes$id), ]
        toAddEdges <- toAddEdges[(toAddEdges$from %in% toAddNodes$id), ]
        edgeId <- edges$id[(edges$id %in% toAddEdges$id)==F]
        nodes1 <- toAddNodes
        edges <- toAddEdges
      }

      totNodes <<- totNodes
      nodes1$color.background <- palette[as.numeric(cut(as.numeric(nodes1$means),breaks=rampbreaks))]
      nodes1$physics <- phys=='On'
      edges <<- edges
      nodes <<- nodes1

      visNetworkProxy('network')  %>% visRemoveNodes(ids) %>% visUpdateNodes(nodes) %>% visRemoveEdges(edgeId) %>% visUpdateEdges(edges)
    }
  }) 
  
  #show previous network when '<' is pressed
  observeEvent(input$back, {
    if(!is.null(prevN1)){
      print("back")
      ids <- nodes[(nodes$id %in% prevN1$id) == F, ][['id']]
      edgeId <- edges[(edges$id %in% prevE1$id)==F, ][['id']]
      nextN1 <<- nodes
      nextE1 <<- edges
      nodes <<- prevN1
      prevN1 <<- NULL
      edges <<- prevE1
      prevE1 <<- NULL
      nodes$color.background <<- palette[as.numeric(cut(as.numeric(nodes$means),breaks=rampbreaks))]

      isolate({
        phys <- input$physics
      })
      nodes$physics <- phys=='On'
      
      visNetworkProxy('network') %>% visUpdateNodes(nodes) %>% visUpdateEdges(edges) %>% visRemoveNodes(ids) %>% visRemoveEdges(edgeId)
    }
  })
  
  #show next network when '>' is pressed
  observeEvent(input$forward, {
    if(!is.null(nextN1)){
      print("forward")

      ids <- nodes[(nodes$id %in% nextN1$id) == F, ][['id']]
      edgeId <- edges[(edges$id %in% nextE1$id)==F, ][['id']]
      
      prevN1 <<- nodes
      prevE1 <<- edges
      nodes <<- nextN1
      nextN1 <<- NULL
      edges <<- nextE1
      nextE1 <<-NULL
      
      isolate({
        phys <- input$physics
      })
      nodes$physics <- phys=='On'
      
      nodes$color.background <<- palette[as.numeric(cut(as.numeric(nodes$means),breaks=rampbreaks))]
      
      visNetworkProxy('network') %>% visUpdateNodes(nodes) %>% visUpdateEdges(edges) %>% visRemoveNodes(ids) %>% visRemoveEdges(edgeId)
    }
  })
  
  #based on button selected, shows either the most disparate nodes between groups or nodes with SD higher than selected
  observe({
    if(!is.null(input$g1) && !is.null(input$group) && input$type=="Most Disparate"){

      survival <- survival[!is.na(survival[[input$group]]), ]
      group1Samps <- as.character(survival$Sample_title[(survival[[input$group]] %in% input$g1)])

      group2Samps <- as.character(survival$Sample_title[!(survival[[input$group]] %in% input$g1)])
      
      

      if(length(group1Samps)>1 && length(group2Samps)>1){
        totNodes$avg1 <- rowMeans(totNodes[,group1Samps])
        totNodes$avg2 <- rowMeans(totNodes[,group2Samps])
      }else if(length(group1Samps)==1 && length(group2Samps)>1){
        totNodes$avg1 <- (totNodes[,group1Samps])
        totNodes$avg2 <- rowMeans(totNodes[,group2Samps])
      }else if(length(group2Samps)==1 && length(group1Samps)>1){
        totNodes$avg1 <- rowMeans(totNodes[,group1Samps])
        totNodes$avg2 <- (totNodes[,group2Samps])
      }else{
        totNodes$avg1 <- totNodes$means
        totNodes$avg2 <- 0
      }
      
      totNodes$p <- apply(totNodes[,4:(ncol(totNodes)-15)], 1, function(i){
        obj <- try(t.test(i[group1Samps], i[group2Samps]), silent=T)
        if(is(obj, "try-error")){
          return(1)
        }else{
          return(obj$p.value)
        }
      })
      
      totNodes$Difference <- (totNodes$avg1-totNodes$avg2)
      totNodes <- totNodes[order((abs(totNodes$p))),]
      totNodes <<- totNodes

      toAddNodes <- totNodes[totNodes$p < input$P, ]
      toAddNodes <- toAddNodes[toAddNodes$Difference > input$disparate, ]
      toAddNodes <- na.omit(toAddNodes)
      rmNodes <- nodes$id[(nodes$id %in% toAddNodes$id)==F]
      toAddEdges <- totEdges[(totEdges$to %in% toAddNodes$id), ]
      toAddEdges <- toAddEdges[(toAddEdges$from %in% toAddNodes$id), ]
      rmEdges <- edges$id[(edges$id %in% toAddEdges$id)==F]
      isolate({
        phys <- input$physics
      })

      if(nrow(toAddNodes)>0){

        toAddNodes$physics <- phys=='On'
      }
      nodes <<- toAddNodes
      edges <<- toAddEdges

      visNetworkProxy('network') %>% visUpdateNodes(nodes) %>% visUpdateEdges(edges) %>% visRemoveNodes(rmNodes) %>% visRemoveEdges(rmEdges)
    }else if(!is.null(nodes)){
      if(is.null(samps)){
        samps = colnames(inflo)[[15]]
      }
      totNodes$means <<- rowMeans(subset(totNodes, select=samps))
      totNodes$SDs <<- apply(subset(totNodes, select=samps), 1, sd)
      
      mean <<- mean(totNodes[['means']])
      sd <<- sd(totNodes$means)
      if(!is.null(input$spread) && !is.na(input$spread) && input$spread<=15 && input$spread>=0){
        level <<- input$spread
      }
      if(length(samps)>1){
        totNodes$SDs <- apply(subset(totNodes, select=samps), 1, sd)
        nodes1 <- totNodes[(abs(totNodes$SDs)>level), ]
        ids <- totNodes[(abs(totNodes$SDs)<=level), ][['id']]
      }else{
        nodes1 <- totNodes[(abs(totNodes$means)>4*sd+mean), ]
        ids <- totNodes[(abs(totNodes$means)<=4*sd+mean), ][['id']]
      }
      edges <- totEdges[(totEdges$to %in% nodes1$id), ]
      edges <- edges[(edges$from %in% nodes1$id), ]
      edgeId <- edges[((totEdges$to %in% nodes1$id)==F && (totEdges$from %in% nodes1$id)==F), ][['id']]

      isolate({
        phys <- input$physics
      })
      
      nodes1$color.background <- palette[as.numeric(cut(as.numeric(nodes1$means),breaks=rampbreaks))]
      nodes1$physics <- phys=='On'
      edges <<- edges
      nodes <<- nodes1

      visNetworkProxy('network')  %>% visRemoveNodes(ids) %>% visUpdateNodes(nodes) %>% visRemoveEdges(edgeId) %>% visUpdateEdges(edges)
    }
  })
  
  
  #shows first neighbors of selected nodes
  observeEvent(input$neighbors, {
    if(!is.null(newNodes)){

      prevN1 <<- nodes
      prevE1 <<- edges
      nextN1 <<- NULL
      nextE1 <<- NULL
      
      
      totNodes$hidden <- F
      totEdges$hidden <- F
      toAddEdge <- totEdges[totEdges$to %in% newNodes, ]
      toAddEdge <- rbind(toAddEdge, totEdges[totEdges$from %in% newNodes, ])
      
      toAddNodes <- totNodes[totNodes$id %in% toAddEdge$to, ]
      toAddNodes <- rbind(toAddNodes, totNodes[totNodes$id %in% toAddEdge$from, ])
      nodes <- rbind(nodes, toAddNodes)
      nodes <- nodes[order(nodes$hidden), ]

      nodes1 <- nodes[!duplicated(nodes$id), ]
      
      from <- (totEdges[totEdges$from %in% nodes1$id, ])
      to <- (totEdges[totEdges$to %in% nodes1$id, ])
      keep <- merge(to, from)
      toAddEdge <- rbind(toAddEdge, keep)
      edges <- rbind(edges, toAddEdge)

      
      edges$hidden <- apply(edges, 1, function(i){
          idTo <- as.numeric(i['to'])
          idFrom <- as.numeric(i['from'])
          return( !(!(nodes1[nodes1$id==idTo, ][['hidden']]) && !(nodes1[nodes1$id==idFrom, ][['hidden']])))
        })

      edges <- edges[order(edges$hidden), ]
      edges1 <- edges[!duplicated(edges$id), ]


      isolate({
        phys <- input$physics
      })
      nodes1$physics <- phys=='On'
      nodes1$color.background <- palette[as.numeric(cut(as.numeric(nodes1$means),breaks=rampbreaks))]
      nodes = nodes1
      edges = edges1
      nodes <<- nodes
      edges <<- edges

      
      visNetworkProxy('network') %>% visUpdateNodes(nodes) %>% visUpdateEdges(edges)
    }
  })

  
  
  
  #changes layout between hierarchical and web when button selected
  #crashes the app with large networks, so only attempted if less that 500 nodes
  observeEvent(input$layout, {
    if(is.null(input$file5) || is.null(input$file6) || is.null(input$file7) || is.null(input$file8)){
      #
    }else{
      isolate({
        phys <- input$physics
      })

      if(input$layout == 'Hierarchical' && nrow(nodes)<500){
        visNetworkProxy('network') %>% visHierarchicalLayout()
      }else{
        nodes$physics <- phys=='On'
        visNetworkProxy('network') %>% visHierarchicalLayout(enabled = F) %>% visUpdateNodes(nodes)
      }
    }
  })
  
  

  #assigns physics based on option selected
  observe({
    if(is.null(input$file5) || is.null(input$file6) || is.null(input$file7) || is.null(input$file8)){
      #
    }else{
      nodes$physics <<- input$physics == 'On'
      visNetworkProxy('network') %>% visUpdateNodes(nodes)
    }
  })
  
  
  
  #shows activity information for each node for selected samples, ordered depending on selection
  output$InFlo <- renderDataTable({
    input$refreshInflo
    if(is.null(input$file5) || is.null(input$file6) || is.null(input$file7) || is.null(input$file8)){
      #
    }else{

      if(input$type == 'Most Deviant'){
        if(is.null(samps)){
          cols <- subset(totNodes, select=c('Parent', 'Parent_Type', colnames(totNodes)[4:(ncol(totNodes)-9)]))
          cols <- cols[order(-(cols$SDs)),]
          return(cols)
        }else if(length(samps)==1){
          cols <- subset(totNodes, select=c('Parent', 'Parent_Type', samps, 'means'))
          cols <- cols[order(-abs(cols$means)),]
          return(cols)
        }else{
          cols <- subset(totNodes, select=c('Parent', 'Parent_Type', samps, 'means', 'SDs'))
          cols <- cols[order(-(cols$SDs)),]
          return(cols) 
        }
      }else{
          if(is.null(samps)){
            cols <- subset(totNodes, select=c('Parent', 'Parent_Type', colnames(totNodes)[c(4:(ncol(totNodes)-9))]))
            cols <- cols[order(-(cols$Difference)),]
            return(cols)
          }else if(length(samps)==1){
            cols <- subset(totNodes, select=c('Parent', 'Parent_Type', samps, 'means', 'Difference', 'p'))
            cols <- cols[order(-(cols$Difference)),]
            return(cols)
          }else{
            cols <- subset(totNodes, select=c('Parent', 'Parent_Type', samps, 'means', 'SDs', 'Difference', 'p'))
            cols <- cols[order(-(cols$Difference)),]
            return(cols) 
          }
        }
      }

  },
  options = list(scrollY = '300px', pageLength=30, scrollX = TRUE))
  
  
  
  #Do survival analysis for all samples and plot survival curve
  output$allSurvivalCurve <- renderPlot({
    input$plot
    if(is.null(input$file5) || is.null(input$file6) || is.null(input$file7) || is.null(input$file8) || nrow(nodes)==0){
      #
    }else{
      print("survive")
      toUse <- nodes[nodes$hidden == F, ]
      means <- colMeans(toUse[,4:(ncol(toUse)-15)])
      highAll <- character()
      lowAll <- character()
      for(i in 1:length(means)){
        if(as.numeric(means[[i]])>0){
          highAll <- append(highAll, colnames(toUse)[[i+3]])
          
        }else{
          lowAll <- append(lowAll, colnames(toUse)[[i+3]])
          
        }
      }
      survival <- survival[is.na(survival$SurvivalDays)==F, ]
      survival$LevelAll <- "Low Activity"
      for(i in 1:nrow(survival)){
        if(survival[i,][['Sample_title']] %in% highAll){
          survival[i,][['LevelAll']] = "High Activity"
        }
      }
      survival <- survival[order(survival$LevelAll),]

      isolate({
      
       if(!is.null(input$group)){
          survival <- survival[!is.na(survival[[input$group]]), ]
          group1Samps <- as.character(survival$Sample_title[(survival[[input$group]] %in% input$g1)])
  
          group2Samps <- as.character(survival$Sample_title[!(survival[[input$group]] %in% input$g1)])
          
          if(length(group1Samps)>0){
            surviveG1 <<- survival[survival[[1]] %in% group1Samps, ]
          }else{
            surviveG1 <<- NULL
          }
          
          if(length(group2Samps)>0){
            surviveG2 <<- survival[survival[[1]] %in% group2Samps, ]
          }else{
            surviveG2 <<- NULL
          }
       }
      })
      
      plot(survfit(Surv(time=survival$SurvivalDays, survival$SurvivalStatus_0_Alive_1_Dead)~survival$LevelAll), conf.int=F, mark.time=T, col=c(2,4), 
           ylab="Survival Probability", xlab="Days", xlim=c(0, max(survival$SurvivalDays+200, na.rm=T)))
      legend('bottomleft', legend=unique(survival$LevelAll), col=c(2,4), lty = 1)
      title("All Sample Survival")
      axis(1, at=seq(0,max(survival$SurvivalDays+200), by=500))
      grid()

    }
  })

  
  
  #Plot group 1 survival curve
  output$G1survivalCurve <- renderPlot({
    input$plot
    if(is.null(input$file5) || is.null(input$file6) || is.null(input$file7) || is.null(input$file8) || is.null(surviveG1)){
      #
    }else{
      max <- max(survival$SurvivalDays, na.rm=T)
      plot(survfit(Surv(time=surviveG1$SurvivalDays, surviveG1$SurvivalStatus_0_Alive_1_Dead)~surviveG1$LevelAll), conf.int=F, mark.time=T, col=c(2,4), 
           ylab="Survival Probability", xlab="Days", xlim=c(0, max+200))
      legend('bottomleft', legend=unique(surviveG1$LevelAll), col=c(2,4), lty = 1)
      title("Group 1 Survival")
      axis(1, at=seq(0,as.numeric(max)+200, by=500))
      grid()
    }
  })
  
  #Plot group 2 survival curve
  output$G2survivalCurve <- renderPlot({
    input$plot
    if(is.null(input$file5) || is.null(input$file6) || is.null(input$file7) || is.null(input$file8) || is.null(surviveG2)){
      #
    }else{
      max <- max(survival$SurvivalDays, na.rm=T)
      plot(survfit(Surv(time=surviveG2$SurvivalDays, surviveG2$SurvivalStatus_0_Alive_1_Dead)~surviveG2$LevelAll), conf.int=F, mark.time=T, col=c(2,4), 
           ylab="Survival Probability", xlab="Days", xlim=c(0, max+200))
      legend('bottomleft', legend=unique(surviveG2$LevelAll), col=c(2,4), lty = 1)
      title("Group 2 Survival")
      axis(1, at=seq(0,as.numeric(max)+200, by=500))
      grid()
    }
  })
  
  #plot bar graph showing relative number of samples with high and low activity for each group
  output$bar <- renderPlot({
    input$plot
    if(is.null(input$file5) || is.null(input$file6) || is.null(input$file7) || is.null(input$file8) || is.null(surviveG1) || is.null(surviveG2)){
      #
    }else{
      high1 <- nrow(surviveG1[surviveG1$LevelAll == "High Activity", ])/nrow(surviveG1)
      low1 <- nrow(surviveG1[surviveG1$LevelAll == "Low Activity", ])/nrow(surviveG1)
      high2 <- nrow(surviveG2[surviveG2$LevelAll == "High Activity", ])/nrow(surviveG2)
      low2 <- nrow(surviveG2[surviveG2$LevelAll == "Low Activity", ])/nrow(surviveG2)
      cols <- c("red", "blue")
      labels <- c("Group1", "Group2")
      leg <- c("High", "Low")
      vals <- matrix(c(high1, high2, low1, low2), nrow=2, ncol=2, byrow=T)
      barplot(vals, names.arg=labels, ylab="Proportion", col=cols)
      legend("topleft", leg, cex=1.3, fill=cols)
    }
  })
  
})