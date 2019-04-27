library(shiny)
library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)

my_edges <- read_csv("ABA_edges.csv")
g <- graph_from_edgelist(as.matrix(my_edges), directed = T)
areas <- V(g)$name
areas <- areas[degree(g, mode = "out") > 2]

trim_and_plot <- function(graph = g, node = "Br_Brain", depth = 3){
    v_trimmed <- V(g)$name[distances(g, node, mode = "out") < depth]
    g1 <- as_tbl_graph(induced_subgraph(g, v_trimmed))
    text_colour <- distances(g1, node, mode = "out")
    text_colour <- as.numeric(text_colour)[order(as.numeric(text_colour))]
    ggraph(g1, 'dendrogram') + 
        geom_edge_link(colour = "lightblue") + 
        geom_node_text(aes(label = name, colour = text_colour), cex = 5, angle = -10) + 
        xlim(-2, 10) +
        ylim(-1, 3) +
        
        theme(axis.line = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              legend.position = "none",
              panel.background = element_blank(),
              panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              plot.background = element_blank())
}

# Define UI for application that draws a histogram
ui <- shinyUI(fluidPage(

    # Application title
    titlePanel("Allen Brain Atlas Area Tree"),

    # Sidebar with a drop down menu
    sidebarLayout(
        sidebarPanel(
            selectInput(inputId = "area1",
                    label = "Brain Area:",
                    choices = areas,
                    selected = areas[1]),
            radioButtons(inputId = "depth1", 
                     label = "How Many Levels Down?",
                     choices = c("2", "3", "4", "5"), 
                     selected = "3", inline = FALSE)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput(outputId = "brain_tree")
        )
    )
))

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$brain_tree <- renderPlot({
        trim_and_plot(depth = input$depth1, node = input$area1, g = g)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
