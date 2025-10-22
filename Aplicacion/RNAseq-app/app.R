library(shiny)
library(ggplot2)
library(ggpubr)
library(corrplot)
library(RColorBrewer)
library(cluster)
library(factoextra)
library(ggdendro)
library(scatterplot3d)
library(DESeq2)
library(edgeR)
library(pheatmap)

library(DT)
library(shinydashboard)
library(shinyWidgets)


ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      html, body {
        height: 100%;
        margin: 0;
        padding: 0;
      }

      .container-fluid {
        display: flex;
        flex-direction: column;
        height: 100vh;
      }

      .custom-header {
        flex: 0 0 25vh;
        background-image: url('fondo.png');
        background-size: cover;
        background-position: center;
        color: white;
        padding: 20px;
        position: relative;
      }

      .custom-header h1 {
        font-size: 3.5em;
        margin-top: 30px;
        font-weight: bold;  /* corregido 'weigth' */
        text-shadow: 1px 1px 2px #000;
      }

      .logo {
        position: absolute;
        top: 10px;
        right: 10px;
        height: 100px; /* antes 60px */
      }


      .main-panel {
        flex: 1 1 auto;
        display: flex;
        overflow: hidden;
      }

      .main-panel .col-sm-3 {
        background-color: #f7f7f7;
        border-right: 1px solid #ddd;
        overflow-y: auto;
        padding-top: 20px;
        height: 100%;
      }

      .main-panel .col-sm-9 {
        overflow-y: auto;
        padding-top: 20px;
        height: 100%;
      }

      .footer {
        flex: 0 0 auto;
        background: linear-gradient(to right, #f0f4f8, #dbe9f4);
        color: #333;
        text-align: center;
        padding: 15px;
        border-top: 1px solid #ddd;
        font-size: 0.9em;
      }
    
    "))
  ),
  
  # Estructura principal
  fluidRow(class = "custom-header",
           column(12,
                  h1("Análisis multivariante de datos RNA-seq"),
                  tags$img(src = "logo.png", class = "logo")
           )
  ),
  
  fluidRow(class = "main-panel",
           column(3,
                  
                  helpText("Paso 1: Seleccione su archivo (.csv/.tsv) con la matriz de expresión."),
                  
                  fileInput("file", "Sube tu archivo CSV/TSV",
                            accept = c(".csv", ".tsv")),
                  
                  helpText("Paso 2: Indique el separador de columnas de su fichero."),
                  
                  radioButtons("sep", "Separador",
                               choices = c("Coma" = ",",
                                           "Punto y coma" = ";",
                                           "Tab" = "\t"),
                               selected = character(0)),
                  
                  
                  helpText("Paso 3: Seleccione cómo se escriben los decimales en su fichero."),
                  
                  radioButtons("dec", "Símbolo decimal",
                               choices = c("Punto" = ".", "Coma" = ","),
                               selected = character(0)),
      
                  helpText("Paso 4: Pulse este botón para cargar los datos."),
                  
                  actionButton("cargar", "Cargar datos"),
                  br(), br(),
                  actionButton("reset", "Reset", icon = icon("redo"))
                  
           ),
           column(9,
                  tabsetPanel(
                    tabPanel("PCA",
                             plotOutput("pcaPlot"),
                             plotOutput("varPlot")
                    ),
                    tabPanel("Análisis Cluster",
                             tabsetPanel(
                               tabPanel("Jerárquico", plotOutput("hclustPlot")),
                               tabPanel("No Jerárquico",
                                        selectInput("opt_method", "Método para estimar número de clusters",
                                                    choices = c("Elbow", "Silhouette")),
                                        plotOutput("optClusterPlot"),
                                        sliderInput("k_input", "Número de clusters", min = 2, max = 10, value = 3),
                                        actionButton("run_kmeans", "Ejecutar K-means"),
                                        plotOutput("kmeansPlot")
                               )
                             )
                    ),
                    tabPanel("Expresión diferencial", tableOutput("tablaDEG"))
                  )
           )
  ),
  
  fluidRow(class = "footer",
           column(12, HTML('
           <p style="text-align: center;">
             <a href="plantilla-tfg.pdf#69" target="_blank" style="margin-right: 30px;">Manual de Usuario</a>
             <a href="https://drive.google.com/uc?export=download&id=1VQpCryXb_OB54bNhIEOC5kvkDHCkW5nJ" style="margin-right: 30px;">Descargar código fuente</a>
             <a href="contacto.html" target="">Contacto</a>
           </p>
         '))
  )
  
  
)


server <- function(input, output, session) {
  datos <- reactiveVal(NULL)
  condiciones <- c("LM", "MT", "PN", "LM", "MT", "PN", "LM", "MT", "PN", "NM", "NM", "NM")
  group <- factor(condiciones)
  
  observeEvent(input$cargar, {
    if (is.null(input$file)) {
      sendSweetAlert(session, "¡Falta el archivo!", "Por favor, seleccione un fichero CSV o TSV antes de cargar.", "error")
      return()
    }
    if (is.null(input$sep)) {
      sendSweetAlert(session, "Seleccione separador", "Debe indicar qué separador de columnas usa su fichero.", "warning")
      return()
    }
    if (is.null(input$dec)) {
      sendSweetAlert(session, "Seleccione símbolo decimal", "Indique si su fichero usa punto o coma para los decimales.", "warning")
      return()
    }
    
    # Leer las primeras líneas para comprobar el separador
    lineas <- readLines(input$file$datapath, n = 5)
    sep <- switch(input$sep, "," = ",", ";" = ";", "\t" = "\t")
    
    # Contar número de columnas en primera línea según separador elegido
    n_col_elegido <- length(strsplit(lineas[1], sep)[[1]])
    
    # Si solo hay una columna, probablemente está mal el separador
    if (n_col_elegido == 1) {
      sendSweetAlert(
        session,
        title = "Separador sospechoso",
        text = paste0("Parece que el separador que ha indicado ('", sep, "') no está presente en el archivo. ¿Seguro que es el correcto?"),
        type = "warning"
      )
      return()
    }
    
    # Si todo está OK, cargar datos
    dec <- switch(input$dec, "." = ".", "," = ",")
    df <- read.csv(input$file$datapath, sep = sep, dec = dec, row.names = 1)
    df <- df[rowSums(df) != 0, ]
    datos(df)
    
    sendSweetAlert(
      session,
      title = "Datos cargados correctamente",
      text = "Ya puede ver los resultados del análisis multivariante en las pestañas de la derecha.",
      type = "success"
    )
  })
  
  
  
  # PCA
  output$pcaPlot <- renderPlot({
    req(datos())
    dgelist <- DGEList(counts = datos(), group = group)
    dgelist <- calcNormFactors(dgelist)
    dgelist <- estimateDisp(dgelist)
    expr_genes <- cpm(dgelist, normalized = TRUE, log = TRUE)
    
    pca <- prcomp(t(expr_genes))
    pca_df <- as.data.frame(pca$x)
    pca_df$cond <- condiciones
    
    ggplot(pca_df, aes(x = PC1, y = PC2, color = cond)) +
      geom_point(size = 4) +
      theme_minimal() +
      ggtitle("PCA - muestras") +
      scale_color_brewer(palette = "Dark2")
  })
  
  output$varPlot <- renderPlot({
    req(datos())
    dgelist <- DGEList(counts = datos(), group = group)
    dgelist <- calcNormFactors(dgelist)
    expr_genes <- cpm(dgelist, normalized = TRUE, log = TRUE)
    pca <- prcomp(t(expr_genes))
    var_exp <- pca$sdev^2 / sum(pca$sdev^2)
    barplot(var_exp, main = "Varianza explicada", col = "skyblue",
            names.arg = paste0("PC", 1:length(var_exp)))
  })
  
  # Clustering jerárquico
  output$hclustPlot <- renderPlot({
    req(datos())
    esc <- scale(datos())
    data_clust <- t(esc)
    dend <- hclust(dist(data_clust), method = "ward.D2")
    plot(dend, 
         main = "Dendrograma - Clustering jerárquico",
         xlab = "Muestras", ylab = "Distancia", sub="")
  })
  
  # Selección de clusters óptimos
  output$optClusterPlot <- renderPlot({
    req(datos())
    data_clust <- t(scale(datos()))
    method <- ifelse(input$opt_method == "Elbow", "wss", "silhouette")
    g <- fviz_nbclust(data_clust, kmeans, method = method)
    g + labs(title = ifelse(method == "wss", "Método del Codo", "Método de la Silueta"),
             x = "Número de clusters", y = ifelse(method == "wss", "Suma de cuadrados intra-cluster", "Anchura de silueta media"))
  })
  
  
  # KMeans clustering
  kmeans_result <- eventReactive(input$run_kmeans, {
    req(datos())
    data_clust <- t(scale(datos()))
    kmeans(data_clust, centers = input$k_input, nstart = 25)
  })
  
  output$kmeansPlot <- renderPlot({
    req(kmeans_result())
    data_clust <- t(scale(datos()))
    fviz_cluster(kmeans_result(), data = data_clust)
  })
  
  # Análisis de expresión diferencial
  output$tablaDEG <- renderTable({
    req(datos())
    
    expr_data <- datos()
    top_genes_idx <- order(rowSums(expr_data), decreasing = TRUE)[1:500]
    expr_top <- expr_data[top_genes_idx, ]
    
    sampleTable <- data.frame(
      row.names = colnames(expr_top),
      condition = factor(condiciones)
    )
    
    dds <- DESeqDataSetFromMatrix(countData = expr_top,
                                  colData = sampleTable,
                                  design = ~ condition)
    
    dds <- DESeq(dds)
    resultados <- results(dds)
    
    resultados$gene <- rownames(resultados)
    
    resultados <- resultados[!is.na(resultados$padj), ]
    resultados <- resultados[order(resultados$padj), ]
    
    head(as.data.frame(resultados), 10)
  })
  
  
  
  # Reset
  observeEvent(input$reset, {
    datos(NULL)
    updateRadioButtons(session, "sep", selected = ",")
    updateRadioButtons(session, "dec", selected = ".")
    updateSliderInput(session, "k_input", value = 3)
  })
}

shinyApp(ui = ui, server = server)

