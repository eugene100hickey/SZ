library(shiny)
library(ggvis)
library(readr)

df <- read_csv("ABA_df.csv")
df <- df[,-1]

gene_names <- function(data){
    ifelse(data$Schiz == 50, 
           paste0("<div style = 'color:red'>", "Gene: ", data$Hgnc_symbol),
           paste0("<div style = 'color:black'>", "Gene: ", data$Hgnc_symbol)
           )
           }

server <- shinyServer(function(input, output) {
    plotData <- reactive({
        df <- df[,c(input$xVariable,
                    input$yVariable, 
                    "Hgnc_symbol",
                    "Colour",
                    "Schiz")]
        names(df) <- c("x","y", "Hgnc_symbol","Colour","Schiz")
        df
    })
    reactive({ plotData() %>%  
            ggvis(x=~x,y=~y, 
                  opacity := input$opacity, 
                  fill := ~Colour,
                  size := ~Schiz,
                  key := ~Hgnc_symbol) %>%
            layer_points() %>%
            add_axis("x", title = input$xVariable) %>%
            add_axis("y", title = input$yVariable) %>% 
            add_tooltip(gene_names, "hover")
        
    }) %>%  bind_shiny("ggvisPlot")
})


ui <- shinyUI(fluidPage(
    titlePanel("ABA Gene Expression for SZ"),
    sidebarLayout(
        sidebarPanel(
            selectInput("xVariable", "X Variable:",
                        names(df)[6:419], selected = names(df)[6]),
            selectInput("yVariable", "Y Variable:",
                        names(df)[6:419], selected = names(df)[200]),
            sliderInput(inputId = "opacity",
                        label = "Set Opacity",
                        min = 0,
                        max = 1,
                        value = 0.3)
        ),
        
        mainPanel(
            ggvisOutput("ggvisPlot")
        )
    )
))


# Create a Shiny app object
shinyApp(ui = ui, server = server)