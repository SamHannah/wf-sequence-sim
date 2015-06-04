shinyUI(fluidPage(
  titlePanel("Word frequency sequence effect simulator--Quartets"),
  
  sidebarLayout(position = "left",
    sidebarPanel(
      fluidRow(
          column(6,
             helpText("This app runs a simulation of a word frequency sequence experiment (with LF/HF pairs forming structured quartets) using Minerva-AL and the criterion calibration rule for
              shifting criterion. User supplies a parameter (L) controlling the encoding fidelity of study words, a parameter (A)
              controlling the accessibility of pre-study words, a criterion shift factor (SF) controlling the lability of the criterion, and a 
              decrement parameter (dp) controlling the rate at which this shift factor decreases across trials is set by the decrement parameter (dp)."),
            
           br(),
              helpText("As the 'L' value increases, the proportion of study words encoded into memory increases.
              As the 'A' value increases, words encountered prior to study become more accessible at test. ")
          ),
                  
          column(6,
             helpText("If SF = 0, then criterion is fixed. If dp = 0, then the lability of the 
              criterion is constant; if dp > 1, then SF quickly comes to alternate between positive and negative
              values. A negative SF value in theory implements response effects")
          ),
                    
          fluidRow(
            column(6,
                   br(),
                   sliderInput("L", "Learning (L)", min = 0.0, max=1.0, step=0.05, value=0.75)
            ), 
            column(6,
                   sliderInput("A", "Pre-study accessibility (A)", min = 0.0, max=1.0, step=0.05, value=0.5)
            ), 
            column(6,
                   sliderInput("SF", "Shift factor (SF)", min = 0.0, max=20.0, step=0.25, value=10.0)
            ),    
            column(6,
                   selectInput("dp", "Decrement parameter (dp)",choices=c(0.0, 0.01, 0.025, 0.05, 0.1, 0.5, 1.25, 2.0), selected = 0.01)
            )
        ),
        fluidRow(
          column(8,
                 selectInput("N_subjects", "Number of iterations per simulation",choices=c(5, 25, 100, 200), selected = 5)
            ),
          column(4,style="margin-top: 25px",
                 actionButton("runSim","Run simulation")
            )
          )
      )
    ),
    # plot the output of simulations given SF and dp values specified by the user
    # plot the criterion for one simulation where the criterion that starts at .3 or less, and the criterion 
    # from a second simulation where the criterion starts at .8 or higher
    mainPanel(
      h3("Mirror-sequence simulation"),
      plotOutput("all", height="900px")
    )
  )
))
