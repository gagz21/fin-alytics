require(shiny)
require(mvtnorm)
require(psych)
require(ggplot2)
require(quantreg)
require(quadprog)
#
risk.sim <- function(input) {
  # Begin enterprise risk simulation
  # test: input <- c(.5, .25, -.5, 1000)
  set.seed(1016) # Freezes the random seed to reproduce results exactly
  n.risks <- 3 # Number of risk factors
  m <- n.risks
  n.sim <- 1000 # pull slider settings into the sigma correlation matrix
  sigma <- matrix(c(1, input[1], input[2], input[1], 1, input[3], input[2], input[3], 1), nrow = m)
  z <- rmvt(n.sim, delta = rep(0, nrow(sigma)),sigma = sigma, df = 6, type = "shifted")
  u <- pt(z, df = 6)
  x1 <- qgamma(u[,1],shape=2,scale=1)
  x2 <- qbeta(u[,2],2,2)
  x3 <- qt(u[,3],df=6)
  factors.df <- cbind(x1/10,x2,x3/10)
  colnames(factors.df) <- c("Revenue", "Variable Cost", "Fixed Cost")
  revenue <- 1000*(1+factors.df[,1])
  variable.cost <- revenue * factors.df[,2]
  fixed.cost <- revenue * factors.df[,3]
  total.cost  <- variable.cost + fixed.cost
  operating.margin  <- revenue - variable.cost - fixed.cost
  analysis.t  <- data.frame(Revenue = revenue, Expense = total.cost, Margin = operating.margin)
  threshold <- input[4]
  analysis.t <- subset(analysis.t, subset = (operating.margin > threshold))
  return(analysis.t)
  }

# End enterprise risk simulation

ui <- fluidPage(
  titlePanel("Enterprise Risk Analytics"),
  sidebarLayout(
    sidebarPanel(
      sliderInput(inputId = "cor.1", 
              label = "Set the Revenue - Variable Cost Correlation", 
              value = 0.0, min = -0.9, max = 0.9),
      sliderInput(inputId = "cor.2", 
              label = "Set the Revenue - Variable Cost Correlation", 
              value = 0.0, min = -0.9, max = 0.9),
      sliderInput(inputId = "cor.3", 
              label = "Set the Variable - Fixed Cost Correlation", 
              value = 0.0, min = -0.9, max = 0.9),
      sliderInput(inputId = "threshold",
              label = "Choose loss threshold",
              value = 0.0, min = 0, max = 2500)
    ),
    mainPanel(
      h3("Metrics Relationships"),
      plotOutput("pairs.1")
    )
  )
)

server <- function(input, output) {
  output$pairs.1 <- renderPlot({
    analysis.t <- risk.sim(c(input$cor.1, input$cor.2, input$cor.3, input$threshold))
    pairs.panels(analysis.t)
  })
}

shinyApp(ui = ui, server = server)
