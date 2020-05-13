#
# This is a Shiny web app implementing the Drake equation with context and the latest numbers.
# - Bayesianify: add Gaussian distributions across each parameter to get probability distribution over number of civs
# - add simulation with frequency distribution of our distance from neighboring stars


# http://www.atlasoftheuniverse.com/50lys.html
# 133 stars similar to the Sun within 50 light-years; probably many Earth-like planets around these stars. 
# There are roughly 1400 star systems within this volume of space containing 2000 stars.

# 33 stars within 12.5 light years
# 250,000 stars within 250 light years


# stars within X light-years: 
# estimated 0.120 stars/cubic parsec number, and using a volume for a distance 
# 100 light-years = 100/3.26 = 30.7 parsecs
# Number = density * volume = 0.120 stars/cubic parsec  * 4/3 pi (30.7 parsecs)^3 = 14,600 stars
# http://teacherlink.ed.usu.edu/tlnasa/reference/imaginedvd/files/imagine/docs/ask_astro/answers/980123d.html

library(shiny)
library(truncnorm)

drake_eq <- function(R, fp, ne, fl, fi, fc, L) {
    expected_num_civs = R*fp*ne*fl*fi*fc*L
    return(expected_num_civs)
}

# Define UI for application
ui <- fluidPage(
    # Application title
    titlePanel("Drake Equation: How many communicating intelligent civilizations might there be in our galaxy?"),
    # Sidebar with a slider input for number of bins 
    withMathJax(),
    tags$div(HTML("<script type='text/x-mathjax-config' >
            MathJax.Hub.Config({
            tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}
            });
            </script >")),
    sidebarLayout(
        sidebarPanel(
            #mathjax.tags$div(HTML("")),
            # should include number of stars currently in the galaxy...
            sliderInput("R", "$R =$ number of stars born in the Milky Way each year", min=0, max=25, value=10),
            helpText("Currently there are 100 billion stars in the Milky Way, with an estimated star formation rate of 1.5-3 stars per year."),
            sliderInput("fp", "$f_p =$ Proportion of stars that have planetary systems", min=0, max=1, value=1, step=.01),
            helpText("Recent estimates put $f_p$ close to 1, i.e. all stars have a planetary system."),
            sliderInput("ne", "$n_e =$ Average number of planets per solar system that are suitable for life", min=.05, max=5, value=.4, step=.05),
            helpText("Estimated 40 billion Earth-sized planets in the habitable zones of sun-like/red dwarf stars, meaning $f_p \\cdot n_e = 0.4$."),
            sliderInput("fl", "$f_l =$ Proportion of planets suitable for life on which life actually evolves.", min=0, max=1, value=.5),
            helpText("?? Perhaps very likely: Life on Earth began around the same time as favorable conditions arose. On the other hand, N=1 and anthropic bias."),
            sliderInput("fi", "$f_i =$ Proportion of life-bearing planets where intelligence evolves", min=0, max=1, value=.5, step=.01),
            helpText("?? Perhaps very likely: Arguably many forms of intelligent life on Earth."),
            sliderInput("fc", "$f_c =$ Proportion of intelligent species that produce interstellar communications.", min=0, max=1, value=.1),
            helpText("?? Unknown."),
            sliderInput("L", "$L =$ Average lifetime of a communicating civilization in years.", min=100, max=10000, value=1000, step=100),
            helpText("?? Historically, human empires tend to last 200 years, but as society and technology advances we hope we are reaching greater stability.")
        ),

        #N = The number of civilizations in the Milky Way Galaxy whose electromagnetic emissions are detectable.
        #R* = The rate of formation of stars suitable for the development of intelligent life.
        #fp = The fraction of those stars with planetary systems.
        #ne = The number of planets, per solar system, with an environment suitable for life.
        #fl = The fraction of suitable planets on which life actually appears.
        #fi = The fraction of life bearing planets on which intelligent life emerges.
        #fc = The fraction of civilizations that develop a technology that releases detectable signs of their existence into space.
        #L = The length of time such civilizations release detectable signals into space.
        
        # Show a plot of the generated distribution
        mainPanel(
           withMathJax(),
           textOutput("printDrakeEqn"),
           textOutput("expectedComCivs"),
           plotOutput("distPlot"),
           textOutput("probMoreThanUs")
        )
    )
)

# Define server logic 
server <- function(input, output) {
    get_drake = reactive({
        drake_eq(input$R, input$fp, input$ne, input$fl, input$fi, input$fc, input$L)
    })
    
    get_drake_distro = reactive({
        minp = 1e-7 # none of the probabilities / Ns can be 0, since Earth/humans exist..
        x = 10000 # number of worlds to simulate
        Rd = rtruncnorm(n=x, a=minp, b=Inf, mean=input$R, sd=1)
        fpd = rtruncnorm(n=x, a=minp, b=1, mean=input$fp, sd=.1) # fraction of stars with planetary systems
        ned = rtruncnorm(n=x, a=minp, b=Inf, mean=input$ne, sd=.5) # number of planets suitable for life per solar system
        fld = rtruncnorm(n=x, a=minp, b=1, mean=input$fl, sd=.1) # fraction of suitable planets on which life actually appears.
        fid = rtruncnorm(n=x, a=minp, b=1, mean=input$fi, sd=.1) # fraction of life-bearing planets on which intelligent life emerges
        fcd = rtruncnorm(n=x, a=minp, b=1, mean=input$fc, sd=.1) # fraction of civilizations that develop and use a detectable communication technology
        Ld = rtruncnorm(n=x, a=50, b=Inf, mean=input$L, sd=50) # years such civilizations release detectable signals into space.
        return(Rd * fpd * ned * fld * fid * fcd * Ld)
    })
    
    output$distPlot <- renderPlot({
        civdist = get_drake_distro()
        civdist[which(civdist<0)] = 0 # can't have negative civilizations
        hist(civdist, col = 'darkgray', border = 'white', xlab="Number of Civilizations", 
             main="Distribution of Broadcasting Civilizations in 10,000 Simulations", 
             freq=T) #  , xlim=c(0,1000)
        abline(v=median(civdist), lty="dotted")
        #text(x=(median(civdist)+20), y=200, paste("Median:\n",round(median(civdist))))
        abline(v=mean(civdist), lty="dashed", col="red")
        #text(x=(mean(civdist)+40), y=400, paste("Mean:\n",round(mean(civdist))), col="red")
        legend("topright",c(paste("Median =",round(median(civdist))),
                            paste("Mean =",round(mean(civdist)))),
               lty = c("dotted", "dashed"), col=c("black", "red"))
        print(mean(civdist))
        #ggplot() + geom_histogram(data = civdist, aes(x = )) + geom_hline(x=mean(civdist))
    })
    
    output$probMoreThanUs <- renderText({
        civdist = get_drake_distro()
        paste("In this simulation, the probability that there is more than one communicating civilization in our galaxy is ", 
              round(sum(civdist>1)/length(civdist),2), ".",sep='')
    })
    
    output$printDrakeEqn <- renderText({ 
        paste("The Drake equation gives the expected number of communicating intelligent civilizations in the Milky Way Galaxy: $N = R*f_p*n_e*f_l*f_i*f_c*L = $") 
    })
    
    output$expectedComCivs <- renderText({ 
        N = get_drake()
        paste(round(N, 2)) 
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
