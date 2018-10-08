

  datos = read.csv("DATA/PalmCurrentDataset.csv")
# Define UI for app that draws a histogram ----
ui <- fluidPage(
  # App title ----
  titlePanel("Accumulation Curves per Palm"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Slider for the number of bins ----
      selectInput(inputId = "bins",
                  label = "Palm Species", 
                  choices = sort(droplevels(unique(datos$PALM)))
                  )
        
    ),
    
    
    # Main panel for displaying outputs ----
    mainPanel(
      # Output: Histogram ----
      plotOutput(outputId = "distPlot"),
      h5("Plot above represents the individual accumulation curves constructed by randomizing 
         the number of frugivores (y axis) in function of the number of unique studies a palm species has been found 
         red line shows the expected assymptote calculated with Chao1, confidence intervals are
         represented with the dashed red lines
         Sampling Completeness (i.e. SC) is calculated as the number of frugivore species observed / expected"),
      h4("Interaction data from Zona and Henderson have been omitted")
    )
    )
)


server <- function(input, output) { 
  
  
  makeFrugPlot = function(dataset, x){ 
    SROm = droplevels(dataset[dataset$PALM == x,])
    SROm1 = table(SROm$FRUGIVORE,SROm$referenceKey)
    Acum = vegan::specpool(SROm1)[c("Species", "chao", "chao.se")]
    plot(vegan::specaccum(SROm1), 
         xlab = "No Studies",
         ylab = "Frugivores", 
         ylim = c(0, Acum$chao + Acum$chao.se + 2),
         main = paste(x,"from:", unique(SROm$biogeographicRegion)))
    abline(h=Acum$chao, col = "red")
    abline(h=c(Acum$chao-Acum$chao.se, Acum$chao+Acum$chao.se), lty = 2, col = "red") 
    legend("topleft" , paste("SC = ", round(Acum$Species/Acum$chao, 3) * 100, "%"),
           bty = "n")}
  
    datos2 = reactive({read.csv("DATA/PalmCurrentDataset.csv")})
  
  output$distPlot <- renderPlot({
    makeFrugPlot(datos2(), input$bins)
})

}

shinyApp(ui, server)






  
  
  
  