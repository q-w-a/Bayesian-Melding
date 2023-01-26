


#########################
# NEW VERSION
#########################
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(latex2exp)
library(future)
library(truncdist)
library(shinyscreenshot)
library(shinyjs)
library(patchwork)


###############################################################
# BETA PARAMETERS FROM DESIRED MEAN AND VARIANCE
###############################################################
get_beta_params <- function(mu, sd) {
    var = sd^2
    alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
    beta <- alpha * (1 / mu - 1)
    return(params = list(alpha = alpha,
                         beta = beta))
}



###########
# TESTING
###########
# library(truncdist)
# mean <- .5
# sd <- .02
#
# params <- get_beta_params(mean,sd)
#
# samp <- rtrunc(1e4, spec = "beta", a = .48, b = .9,
#         shape1= params$alpha,
#         shape2 = params$beta)
#
# mean(samp)
# sd(samp)



#######################################################################
# BETA DENSITY WITH DESIRED MEAN AND VARIANCE (WITH TRUNCATION OPTION)
#######################################################################
# beta_density <- function(x, mean, sd, bounds=NULL) {
#
#     shape_params <-  get_beta_params(
#         mu = mean,
#         sd = sd)
#
#     if(!is.null)
#     dbeta(x,
#           shape1 = shape_params$alpha,
#           shape2 = shape_params$beta)
#
#     dtrunc(1e4, spec = "beta", a = .48, b = .9,
#            #         shape1= params$alpha,
#            #         shape2 = params$beta)
#
# }





###############################################################
# GAMMA PARAMETERS FROM DESIRED MEAN AND VARIANCE
###############################################################
get_gamma_params <- function(mu, sd) {
    var = (mu/sd)^2
    shape = (mu/sd)^2
    scale = sd^2/mu
    return(params = list(shape = shape,
                         scale = scale))
}



################################
# TEST GAMMA PARAM FUNCTION
################################
# mean <- 2
# std <- 1.5
#
# params <- get_gamma_params(mu = mean, sd = std)
#
#
# samp <- rgamma(1e4,
#                shape = params$shape,
#                scale = params$scale)
#
# mean(samp)
# sd(samp)




###############################################################
# BETA DENSITY WITH DESIRED MEAN AND VARIANCE
###############################################################
beta_density <- function(x, mean, sd, bounds=NA) {
    shape_params <-  get_beta_params(
        mu = mean,
        sd = sd)

    if(!length(bounds) == 1){
        # message("here")
        dtrunc(x,
               spec = "beta",
               a = bounds[1],
               b = bounds[2],
              shape1 = shape_params$alpha,
              shape2 = shape_params$beta) %>%
            return()
    }else{
        dbeta(x,
          shape1 = shape_params$alpha,
          shape2 = shape_params$beta)  %>%
            return()
        }
}




###############################################################
# SAMPLE FROM BETA DENSITY WITH DESIRED MEAN AND VARIANCE
###############################################################

sample_beta_density <- function(n, mean, sd, bounds = NA) {

    shape_params <-  get_beta_params(
        mu = mean,
        sd = sd)

    rbeta(n,
          shape1 = shape_params$alpha,
          shape2 = shape_params$beta)

    if(!length(bounds) == 1){
        # message("here")
        rtrunc(n,
               spec = "beta",
               a = bounds[1],
               b = bounds[2],
               shape1 = shape_params$alpha,
               shape2 = shape_params$beta) %>%
            return()
    }else{
        rbeta(n,
              shape1 = shape_params$alpha,
              shape2 = shape_params$beta)  %>%
            return()
    }
}


#
# x <- seq(0,1,length =1000)
#
# tibble(d = beta_density(x,
#                         mean = .9, sd = .04, bounds = c(.8,.95))) %>%
#     ggplot(aes(x=x,y=d)) +
#     geom_line()


###############################################################
# GAMMA DENSITY WITH DESIRED MEAN AND VARIANCE
###############################################################
gamma_density <- function(x, mean, sd, bounds=NA) {

    shape_params <-  get_gamma_params(
        mu = mean,
        sd = sd)

    if(!length(bounds) == 1){
        #message("here")
        dtrunc(x,
               spec = "gamma",
               a = bounds[1],
               b = bounds[2],
               shape = shape_params$shape,
               scale = shape_params$scale) %>%
            return()
    }else{
        dgamma(x,
               shape = shape_params$shape,
               scale = shape_params$scale) %>%
            return()
    }
}
#
# x <- seq(0,1.5,length =1000)
#
# tibble(d = gamma_density(x,
#                         mean = 1.2, sd = .04, bounds = c(.8,1.3))) %>%
#     ggplot(aes(x=x,y=d)) +
#     geom_line()



sample_gamma_density <- function(n, mean, sd, bounds = NA) {

    shape_params <-  get_gamma_params(
        mu = mean,
        sd = sd)

    if(!length(bounds) == 1){
        #message("here")
        rtrunc(n,
               spec = "gamma",
               a = bounds[1],
               b = bounds[2],
               shape = shape_params$shape,
               scale = shape_params$scale) %>%
            return()
    }else{
        rgamma(n,
               shape = shape_params$shape,
               scale = shape_params$scale) %>%
            return()
    }
}


# sample_gamma_density(10, mean = .9, sd= .04)
#
# sample_gamma_density <- function(n, mean, sd) {
#
#     shape_params <-  get_gamma_params(
#         mu = mean,
#         sd = sd)
#
#     rgamma(n,
#            shape = shape_params$shape,
#            scale = shape_params$scale)
#
#
# }





###############################################################
# INDUCED PRIOR ON ASYMPTOMATIC RATE  P(S_0|test+,untested)
###############################################################
# q_1^*(\theta)

# input sampled values of theta and compute M(\theta)
est_P_A_testpos = function(P_S_untested, alpha, beta){
    beta * (1 - P_S_untested) / (( beta * (1 - P_S_untested)) + (alpha * P_S_untested))
}


get_melded <- function(alpha_mean = 0.9,
                       alpha_sd = 0.04,
                       alpha_bounds = NA,
                       beta_mean = .15,
                       beta_sd =.09,
                       beta_bounds = NA,
                       s_untested_mean = .025,
                       s_untested_sd = .0225,
                       s_untested_bounds = NA,
                       p_s0_pos_mean = .4,
                       p_s0_pos_sd = .1225,
                       p_s0_pos_bounds = NA,
                       nsamp = 1e3) {

  given_args <- as.list(environment())
  cat("Arguments to get_melded:\n")
  print(given_args)


    theta <- tibble(alpha = sample_gamma_density(nsamp,
                                                mean = alpha_mean,
                                                sd = alpha_sd,
                                                bounds = alpha_bounds),
                    beta= sample_beta_density(nsamp,
                                              mean = beta_mean,
                                              sd = beta_sd,
                                              bounds = beta_bounds),
                    P_S_untested = sample_beta_density(nsamp,
                                                       mean = s_untested_mean,
                                                       sd = s_untested_sd,
                                                       bounds = s_untested_bounds)) %>%
        mutate(phi_induced = est_P_A_testpos(P_S_untested = P_S_untested,
                                             alpha = alpha,
                                             beta=beta))

    # theta contains values sampled from alpha, beta, P_S_untested, and M(theta) = phi_induced
    # induced phi
    phi <- theta$phi_induced

    # approximate induced distribution via a density approximation
    phi_induced_density <- density(x = phi, n = nsamp, adjust = 2, kernel = "gaussian")


    incProgress(.2, detail = paste("Computing induced density of phi..."))

    # future::plan(multisession, workers = 3)
    # tictoc::tic()
    # indexes <- furrr::future_map(phi, ~{
    #     which(phi_induced_density$x > .x)[1] }) %>%
    #     unlist()
    # tictoc::toc()

    # FASTER IMPLEMENTATION
    indexes <- findInterval(phi, phi_induced_density$x)


    phi_sampled_density <- phi_induced_density$y[indexes]

    dp_s0_pos <- function(x) {

      beta_density(x,
                   mean=p_s0_pos_mean,
                   sd = p_s0_pos_sd,
                   bounds=p_s0_pos_bounds)
    }

  #  message("CLASS----", class(dp_s0_pos))

    incProgress(.6, detail = paste("Calculating weights..."))


    # weights <- purrr::map2_dbl(
    #     phi_sampled_density,
    #     phi,
    #     function(phi_sampled_density_i, phi_i) {
    #         # pooling weight
    #         alpha = .5
    #         (phi_sampled_density_i/ dp_s0_pos(phi_i))^(1-alpha)
    #     }
    # )

    weights <- (phi_sampled_density/ dp_s0_pos(phi))^(.5)


    post_samp_ind <-sample.int(n=nsamp,
                               size=nsamp,
                               prob=1/weights,
                               replace=TRUE)


    pi_samp <- cbind(theta[post_samp_ind,],
                     P_A_testpos =  phi[post_samp_ind]) %>%
        select(-phi_induced)

    pi_samp_long <- pi_samp %>%
        pivot_longer(cols=everything()) %>%
        mutate(type = "After Melding")


    melded <- theta %>%
        mutate(P_A_testpos = sample_beta_density(nsamp,
                                                 mean = p_s0_pos_mean,
                                                 sd = p_s0_pos_sd,
                                                 bounds = p_s0_pos_bounds)) %>%
        pivot_longer(cols=everything()) %>%
        mutate(type = ifelse(
            name == "phi_induced",
            "Induced", "Before Melding")) %>%
        mutate(name = ifelse(name == "phi_induced",
                             "P_A_testpos",
                             name)) %>%
        bind_rows(pi_samp_long) %>%
        mutate(name = case_when(
            name == "alpha" ~"$\\alpha$",
            name == "beta" ~"$\\beta$",
            name == "P_A_testpos" ~ "$P(S_0|test+,untested)$",
            name == "P_S_untested" ~ "$P(S_1|untested)$")
        ) %>%
        mutate(name = factor(name,
                             levels = c(
                                 "$\\alpha$",
                                 "$\\beta$",
                                 "$P(S_1|untested)$",
                                 "$P(S_0|test+,untested)$")))

    incProgress(.2, detail = paste("Generating plot..."))


    return(melded)

}

# test <- get_melded()



# Define UI for application that draws a histogram
ui <- fluidPage(

    tags$head(
        # Note the wrapping of the string in HTML()
        tags$style(HTML("
        #beta_plot, #alpha_plot, #s_untested_plot, #asymp_plot { overflow:hidden;  }
        #submit_melding, #reset_input {color: #fff;
                        background-color: #337ab7;
                        border-color: #2e6da4;
                        font-weight:bold;
                        size:2vw;}
        #submit_melding:hover, #reset_input:hover {
                        background-color:#78CFDE;
                        border-color:gold;
                        border-width:2px;}
        #screenshot {
                        float:right;
                        color: #fff;
                        background-color: #4A706D;
                        border-color: #2e6da4;
                        font-weight:bold;
                        size:2vw;}
        #screenshot:hover {
                        background-color:#659D99;
                        border-color:gold;
                        border-width: 2px;}
        input[type=checkbox] { transform: scale(1.8);}
        #loadmessage {
           //  position: fixed;
             float:center;
             top: 50%;
             left: 0px;
             width: 100%;
             padding: 10px 0px 10px 0px;
             text-align: center;
             font-weight: bold;
             font-size: 180%;
             color: #000000;
             background-color: #659D99;
             z-index: 105;
             border-radius: 15px;
           }",
        ))),


    # Application title
    titlePanel("Bayesian Melding"),
    withMathJax(),
    useShinyjs(),


    # Sidebar with a slider inputs for bias parameters
    sidebarLayout(

        sidebarPanel(
          # conditionalPanel(condition="$('html').hasClass('shiny-busy')",
          #                  tags$div("Loading...",id="loadmessage") ),

          selectInput("nsamp",
                      label =  "$$\\textbf{Number of Samples for } \\\\ \\textbf{Bayesian Melding:}$$",
                      choices = list(1e3,
                                     1e4,
                                     1e5,
                                     1e6),
                      selected = 1e5),

            # sliderInput("nsamp",
            #             "$$\\textbf{Number of Samples for } \\\\ \\textbf{Bayesian Melding:}$$",
            #             min = 1e3,
            #             max = 1e5,
            #             value = 1e3),
            h4("$$\\underline{\\;\\textbf{Distributions for }\\theta\\;\\;}$$"),
            h5("$$\\textbf{Distribution of } \\alpha\\\\(\\textit{Gamma Distribution})$$"),

            sliderInput("alpha_mean",
                        "$$\\text{ Mean for } \\alpha$$",
                        min = 0,
                        max = 5,
                        value = 0.9,
                        step = 0.1),
            sliderInput("alpha_sd",
                        "$$\\text{ Standard deviation for } \\alpha$$",
                        min = 0,
                        max = 2,
                        value = 0.04,
                        step = 0.01),
        #    p("$$\\text{ Truncation Bounds for } \\alpha\\\\\\textit{(Default is No Truncation)}$$"),
            p("$$\\text{ Truncation Bounds for } \\alpha$$"),

            checkboxInput(inputId = "trunc_alpha", label = "Truncate", value = TRUE),
            sliderInput(inputId = "alpha_bounds",
                        label = "",
                        #"$$\\text{ Truncation Bounds for } \\alpha\\\\\\textit{(Default is No Truncation)}$$",
                        min = 0,
                        max = 4,
                        value=c(.8, 1),
                        step = .1),
            br(),

            h5("$$\\textbf{Distribution of } \\beta\\\\(\\textit{Beta Distribution})$$"),
            sliderInput(inputId = "beta_mean",
                        "$$\\text{ Mean for } \\beta$$",
                        min = 0,
                        max = 1,
                        value = 0.15),
            sliderInput("beta_sd",
                        "$$\\text{ Standard deviation for } \\beta$$",
                        min = 0,
                        max = 0.5,
                        value = 0.09),
 #           p("$$\\text{ Truncation Bounds for } \\beta\\\\\\textit{(Default is No Truncation)}$$"),
            p("$$\\text{ Truncation Bounds for } \\beta$$"),

            checkboxInput(inputId = "trunc_beta", label = "Truncate", value = TRUE),
            sliderInput(inputId = "beta_bounds",
                        label = "",
                        min =0,
                        max = 1,
                        value=c( 0.002, 0.4),
                        step = 0.001),
            br(),

            # h5("$$\\textbf{Distribution of } P(S_1|untested)\\\\(\\textit{Beta Distribution})$$"),
            h5("$$\\textbf{Distribution of } P(S_1|untested)\\\\(\\textit{Beta Distribution})$$"),
            sliderInput("s_untested_mean",
                        "$$\\text{ Mean for } P(S_1|untested)$$",
                        min = 0,
                        max = 1,
                        value = 0.025,
                        step = .001),
            sliderInput("s_untested_sd",
                        "$$\\text{ Standard deviation for }\\\\ P(S_1|untested)$$",
                        min = 0,
                        max = 0.5,
                        value =  0.0225,
                        step = .0001),
          #  p("$$\\text{ Truncation Bounds for } P(S_1|untested)\\\\\\textit{(Default is No Truncation)}$$"),
            p("$$\\text{ Truncation Bounds for } P(S_1|untested)$$"),

            checkboxInput(inputId = "trunc_s_untested", label = "Truncate", value = TRUE),
            sliderInput(inputId = "s_untested_bounds",
                        label = "",
                        min = 0,
                        max = 1,
                        value=c(0, .15)),

            h4("$$\\underline{\\;\\textbf{Distribution for }\\phi\\;\\;}\\\\(\\textit{Beta Distribution})$$"),
            sliderInput("p_s0_pos_mean",
                        "$$\\text{ Mean for }\\\\ P(S_0|untested, test +)$$",
                        min = 0,
                        max = 1,
                        value =  0.4),
            sliderInput("p_s0_pos_sd",
                        "$$\\text{ Standard deviation for }\\\\ P(S_0|untested, test +)$$",
                        min = 0,
                        max = 0.5,
                        value =  0.1225,
                        step = .0001),
            # p(paste0("$$\\text{ Truncation Bounds for }\\\\ P(S_0|untested, test +)",
            # "\\\\\\textit{(Default is No Truncation)}$$")),
            p("$$\\text{ Truncation Bounds for }\\\\ P(S_0|untested, test +)$$"),
            checkboxInput(inputId = "trunc_p_s0_pos",
                          label = "Truncate",
                          value = TRUE),
           # checkboxInput(inputId = "trunc_p_s0_pos", label = "Truncate", value = FALSE),
            sliderInput(inputId = "p_s0_pos_bounds",
                        label = "",
                        min = 0,
                        max = 1,
                        value=c(.25, .7),
                        step = .01),
            actionButton(inputId ="reset_input",
                         label ="Reset inputs"),
            width = 3
        ),

        # Show a plot of the generated distribution
        mainPanel(
         # sidebarPanel(
            conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                             tags$div("Loading...",id="loadmessage") ),
            br(),

            # run melding upon click of button
            actionButton(
                inputId = "submit_melding",
                label = "Run Melding"
            ),
            actionButton(inputId = "screenshot",
                         label ="Take a screenshot"),
            br(), br(),
            p(HTML(paste0("<b></i>Instructions:</b></i> Change input parameters and press submit melding to generate ",
            "the melded distributions. Melding is only performed once <i>Submit Melding</i> is pressed."))),
            h4(paste0("$$\\textbf{Distributions for }\\mathbf{\\boldsymbol{\\theta}",
            "= \\{ \\boldsymbol{\\alpha, \\beta,} P(S_1|untested)\\}}$$")),
            fluidRow(
                splitLayout(cellWidths = c("33%", "33%", "33%"),
                            plotOutput("alpha_plot",  height = "200px"),
                            plotOutput("beta_plot",   height = "200px"),
                            plotOutput("s_untested_plot",   height = "200px"))),
            fluidRow(
                splitLayout(
                    cellWidths = c("50%"),
                    plotOutput("induced_plot", height = "300px"),
                    plotOutput("asymp_plot", height = "300px"))),
            br(),
            h4("$$\\textbf{Distributions After Bayesian Melding:}$$"),
            fluidRow(
                plotOutput("melded_plot", height = "700px")),
            width = 9

        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {

    observeEvent(input$screenshot, {
        time <- gsub(":", ".", Sys.time())
        time <- gsub(" ", "_", time)
        file_name=paste0("shiny_", time)
        screenshot(filename=file_name)
    })

    observeEvent( eventExpr = input$reset_input,
                  priority = 1,
                  handlerExpr = {
      #  input$reset_input

        updateNumericInput(session, "alpha_mean", value = 0.9)
        updateNumericInput(session, "alpha_sd", value = 0.04)
        updateNumericInput(session, "alpha_bounds", value=c(.8, 1))
        updateCheckboxInput(session, "trunc_alpha", value = TRUE)


        updateNumericInput(session, "beta_mean", value = 0.15)
        updateNumericInput(session, "beta_sd", value = 0.09)
        updateNumericInput(session, "beta_bounds", value = c(0.002,.4))
        updateCheckboxInput(session, "trunc_beta", value = TRUE)


        updateNumericInput(session, "s_untested_mean", value = .025)
        updateNumericInput(session, "s_untested_sd", value = 0.0225)
        updateNumericInput(session, "s_untested_bounds", value = c(0,.15))
        updateCheckboxInput(session, "trunc_s_untested", value = TRUE)

        updateNumericInput(session, "p_s0_pos_mean", value = 0.4)
        updateNumericInput(session, "p_s0_pos_sd", value = 0.1225)
        updateNumericInput(session, "p_s0_pos_bounds", value = c(0.25,.7))
        updateCheckboxInput(session, "trunc_p_s0_pos", value = TRUE)

     #   input$submit_melding <<- TRUE
       # hide("melded_plot")

    })


    output$alpha_plot <- renderPlot({

        # print(input$alpha_mean)
        # print(input$alpha_sd)


        inp_bounds <- if(input$trunc_alpha) input$alpha_bounds else NA


        # message(input$alpha_bounds)
        # message(input$trunc_alpha)
        # message("alpha bounds", inp_bounds)
        # message(is.na(inp_bounds))

        tibble(Value = sample_gamma_density(1e5,
                                            mean = input$alpha_mean,
                                            sd = input$alpha_sd,
                                            bounds = inp_bounds)) %>%
            ggplot(aes(x = Value)) +
            geom_density(fill = "#4F4B68",
                         color = "#4F4B68",
                         size = .2) +
            theme_bw() +
            theme(plot.title = element_text(hjust = .5, size = 16),
                  plot.subtitle = element_text(hjust = .5, size = 11),
                  axis.title = element_text(size = 16),
                  axis.text.x = element_text(size = 12)) +
            labs(title = latex2exp::TeX("Prior for $\\alpha$"),
                 subtitle = ifelse(length(inp_bounds) ==2,
                                   paste0("Mean = ",
                                   input$alpha_mean,
                                   ", Standard Deviation = ",
                                   input$alpha_sd,
                                   "\nTruncation Bounds: ",
                                   "(", inp_bounds[1], ", ", inp_bounds[2], ")"),
                                   paste0("Mean = ",
                                          input$alpha_mean,
                                          ", Standard Deviation = ",
                                          input$alpha_sd)),

                 y = "Density",
                 x = "Value")  +
            # add point just to ensure axis starts at 0
            geom_point(aes(x=0, y = 0), size=0) +
            scale_x_continuous(n.breaks = 6,
                               limits = c(0,
                                          input$alpha_mean + 8*input$alpha_sd))

        # PREVIOUS VERSION
        # ggplot() +
        #     stat_function(fun =gamma_density,
        #                   geom="area",
        #                   args = list("mean" = input$alpha_mean,
        #                               "sd" = input$alpha_sd)
        #                   # ,
        #                   # xlim= c(0,3),
        #                   # n = 300
        #                   ) +
        #     theme_bw() +
        #     theme(plot.title = element_text(hjust = .5, size = 16),
        #           plot.subtitle = element_text(hjust = .5, size = 11),
        #           axis.title = element_text(size = 16)) +
        #     labs(title = latex2exp::TeX("Prior for $\\alpha$"),
        #          subtitle = paste0("Mean = ",
        #                            input$alpha_mean,
        #                            ", Standard Deviation = ",
        #                            input$alpha_sd ),
        #          y = "Density",
        #          x = "Value")
    })

    output$beta_plot <- renderPlot({

        inp_bounds <- if(input$trunc_beta) input$beta_bounds else NA


        ggplot() +
            stat_function(fun = beta_density,
                          geom="area",
                          args = list("mean" = input$beta_mean,
                                      "sd" = input$beta_sd,
                                      "bounds" = inp_bounds),
                          fill = "#4F4B68") +
            theme_bw() +
            theme(plot.title = element_text(hjust = .5, size = 16),
                  plot.subtitle = element_text(hjust = .5, size = 11),
                  axis.title = element_text(size = 16),
                  axis.text.x = element_text(size = 12)) +
            labs(title = latex2exp::TeX("Prior for $\\beta$"),
                 # subtitle = paste0("Mean = ",
                 #                   input$beta_mean,
                 #                   ", Standard Deviation = ",
                 #                   input$beta_sd ),
                 subtitle = ifelse(length(inp_bounds) ==2,
                                   paste0("Mean = ",
                                          input$beta_mean,
                                          ", Standard Deviation = ",
                                          input$beta_sd,
                                          "\nTruncation Bounds: ",
                                          "(", inp_bounds[1], ", ", inp_bounds[2], ")"),
                                   paste0("Mean = ",
                                          input$beta_mean,
                                          ", Standard Deviation = ",
                                          input$beta_sd)),
                 y = "Density",
                 x = "Value") +
            xlim(0,1)
    })

    output$s_untested_plot <- renderPlot({
        inp_bounds <- if(input$trunc_s_untested) input$s_untested_bounds else NA

        ggplot() +
            stat_function(fun = beta_density,
                          geom="area",
                          args = list("mean" = input$s_untested_mean,
                                      "sd" = input$s_untested_sd,
                                      "bounds" = inp_bounds),
                          fill = "#4F4B68") +
            theme_bw() +
            theme(plot.title = element_text(hjust = .5, size = 16),
                  plot.subtitle = element_text(hjust = .5, size = 11),
                  axis.title = element_text(size = 16),
                  axis.text.x = element_text(size = 12)) +
            labs(title = latex2exp::TeX("Prior for $\\P(S_1|untested)$"),
                 # subtitle = paste0("Mean = ",
                 #                   input$s_untested_mean,
                 #                   ", Standard Deviation = ",
                 #                   input$s_untested_sd ),
                subtitle = ifelse(length(inp_bounds) ==2,
                        paste0("Mean = ",
                               input$s_untested_mean,
                               ", Standard Deviation = ",
                               input$s_untested_sd,
                               "\nTruncation Bounds: ",
                               "(", inp_bounds[1], ", ", inp_bounds[2], ")"),
                        paste0("Mean = ",
                               input$s_untested_mean,
                               ", Standard Deviation = ",
                               input$s_untested_sd)),
                 y = "Density",
                 x = "Value") +
            xlim(0,1)
    })


    output$induced_plot <- renderPlot({

        nsamp <- 1e5
        theta <- tibble(alpha = sample_gamma_density(nsamp,
                                                    mean = input$alpha_mean,
                                                    sd = input$alpha_sd,
                                                    bounds = input$alpha_bounds),
                        beta = sample_beta_density(nsamp,
                                                   mean = input$beta_mean,
                                                   sd = input$beta_sd,
                                                   bounds = input$beta_bounds),
                        s_untested = sample_beta_density(nsamp,
                                                         mean = input$s_untested_mean,
                                                         sd = input$s_untested_sd,
                                                         bounds = input$s_untested_bounds),
                        induced = est_P_A_testpos(P_S_untested = s_untested,
                                                  alpha = alpha,
                                                  beta=beta))
        theta %>%
            ggplot(aes(x = induced)) +
            geom_density(fill ="#DEB578", color = "#DEB578", alpha = .7) +
            theme_bw() +
            theme(plot.title = element_text(hjust = .5, size = 16),
                  plot.subtitle = element_text(size = 14, hjust = .5),
                  axis.title = element_text(size = 16),
                  axis.text.x = element_text(size = 12)) +
            labs( title = latex2exp::TeX("Induced Prior for $P(S_0|untested,test+)$"),
                  subtitle = latex2exp::TeX("Density of $M(\\theta)$"),
                  y = "Density",
                  x = "Value")
    })


    output$asymp_plot <- renderPlot({

     # message("checked?",input$trunc_p_s0_pos)

      inp_bounds <- if(input$trunc_p_s0_pos) input$p_s0_pos_bounds else NA


        ggplot() +
            stat_function(fun = beta_density,
                          geom="area",
                          args = list("mean" = input$p_s0_pos_mean,
                                      "sd" = input$p_s0_pos_sd,
                                      "bounds" = inp_bounds),
                          fill = "#166C89",
                          alpha = .8) +
            theme_bw() +
            theme(plot.title = element_text(hjust = .5, size = 16),
                  axis.title = element_text(size = 16),
                  plot.subtitle = element_text(size = 11, hjust = .5),
                  axis.text.x = element_text(size = 12)) +
            labs(
                title = latex2exp::TeX("Prior for $\\P(S_0|untested, test +)$"),
                # subtitle = paste0("Based on Meta-analyses on Asymptomatic Rate",
                #                   "\n", "Mean = ",  input$p_s0_pos_mean,
                #                   ", Standard Deviation = ", input$p_s0_pos_sd),
               subtitle = ifelse(length(inp_bounds) ==2,
                       paste0("Based on Meta-analyses on Asymptomatic Rate\n",
                              "Mean = ",
                              input$p_s0_pos_mean,
                              ", Standard Deviation = ",
                              input$p_s0_pos_sd,
                              "\nTruncation Bounds: ",
                              "(", inp_bounds[1], ", ", inp_bounds[2], ")"),
                       paste0("Mean = ",
                              input$p_s0_pos_mean,
                              ", Standard Deviation = ",
                              input$p_s0_pos_sd)),
                y = "Density",
                x = "Value")
    })


    observeEvent(
    #  (input$submit_melding | input$reset_input) ,
      #eventExpr = (input$submit_melding !=0 | input$reset_input != 0),
      eventExpr = input$submit_melding,
      priority = -1,
      handlerExpr = {
      message("reset input: ", input$reset_input)
      message("submit melding: ", input$submit_melding)

        if(input$submit_melding == 0 && input$reset_input==0){
          return()
        }


      # showModal(modalDialog("Doing a function", footer=NULL))
     # show("melded_plot")

      output$melded_plot <- renderPlot({


      #  show(
        # input$submit_melding
        # req(input$submit_melding)

        # if(output$melding_plot == "hide"){
        #   return()
        # } else {

        isolate( {

          alpha_bounds <- if(input$trunc_alpha) input$alpha_bounds else NA
          beta_bounds <- if(input$trunc_beta) input$beta_bounds else NA
          p_s0_pos_bounds <- if(input$trunc_p_s0_pos) input$p_s0_pos_bounds else NA
          s_untested_bounds <- if(input$trunc_s_untested) input$s_untested_bounds else NA

          cat("----| Bounds: |----\n",
              "alpha:", alpha_bounds, "\n",
              "beta:", beta_bounds, "\n",
              "P(S1|untested)", s_untested_bounds, "\n",
              "P(S0|untested,+)",p_s0_pos_bounds, "\n",
              "----------------")

            cat(paste0("-----\n",
                       "Sample size:", input$nsamp, "\n",
                       "Alpha mean: ", input$alpha_mean, "\n",
                       "Alpha sd: ", input$alpha_sd, "\n",
                       "Beta mean: ", input$beta_mean, "\n",
                       "Beta sd: ", input$beta_sd, "\n",
                       "P(S1|untested) mean: ", input$s_untested_mean, "\n",
                       "P(S1|untested) sd: ", input$s_untested_sd, "\n",
                       "P(S0|untested,+) mean: ", input$p_s0_pos_mean, "\n",
                       "P(S0|untested,+) sd: ", input$p_s0_pos_sd, "\n",
                       "-----\n"))


            withProgress(message = "Performing melding computation",
                         value = 0,
                         {

                             melded <- get_melded(alpha_mean =input$alpha_mean,
                                                  alpha_sd = input$alpha_sd,
                                                  alpha_bounds = alpha_bounds,
                                                  beta_mean = input$beta_mean,
                                                  beta_sd = input$beta_sd,
                                                  beta_bounds = beta_bounds,
                                                  s_untested_mean = input$s_untested_mean,
                                                  s_untested_sd = input$s_untested_sd,
                                                  s_untested_bounds = s_untested_bounds,
                                                  p_s0_pos_mean = input$p_s0_pos_mean,
                                                  p_s0_pos_sd = input$p_s0_pos_sd,
                                                  p_s0_pos_bounds = p_s0_pos_bounds,
                                                  nsamp = as.numeric(input$nsamp))


                            p1 <- melded %>%
                              filter(name != "$P(S_0|test+,untested)$") %>%
                                 ggplot(aes(x = value, fill = type)) +
                                 geom_density(alpha = .5, show.legend=FALSE) +
                                 facet_wrap(~name,
                                            labeller = as_labeller(
                                                TeX,
                                                default = label_parsed),
                                            ncol = 3,
                                            scales = "free_y") +
                                 theme_bw() +
                                 theme(
                                       # axis.text.y = element_blank(),
                                       # axis.ticks.y = element_blank(),
                                       axis.title = element_text(size = 18),
                                       axis.text.x = element_text(size = 11),
                                       plot.title =element_text(size = 18,
                                                                margin =margin(0,0, .5,0, 'cm')),
                                       axis.text.y = element_text(size = 11),
                                       strip.text = element_text(size = 16),
                                       legend.text = element_text(size = 16)) +
                                 labs(title = paste0("Number of Samples: ", input$nsamp),
                                      fill = "",
                                      y = "Density") +
                                 scale_fill_manual(values = c("#5670BF",
                                                              "#418F6A",
                                                              "#B28542")) +
                                 guides(fill = guide_legend(keyheight = 2,
                                                            keywidth = 2))  +
                                scale_x_continuous(n.breaks = 5,
                                                   limits = c(0,input$alpha_mean + 2*input$alpha_sd))

                          p2 <- melded %>%
                            filter(name == "$P(S_0|test+,untested)$") %>%
                            ggplot(aes(x = value, fill = type)) +
                            geom_density(alpha = .5) +
                            facet_wrap(~name,
                                       labeller = as_labeller(
                                         TeX,
                                         default = label_parsed),
                                       ncol = 3,
                                       scales = "fixed") +
                            theme_bw() +
                            theme(
                              # axis.text.y = element_blank(),
                              # axis.ticks.y = element_blank(),
                              axis.title = element_text(size = 18),
                              axis.text.x = element_text(size = 11),
                              axis.text.y = element_text(size = 11),
                              plot.title =element_text(size = 18,
                                                       margin =margin(0,0, .5,0, 'cm')),
                              strip.text = element_text(size = 16),
                              legend.text = element_text(size = 19)) +
                            labs(
                              # title = paste0("Number of Samples: ", input$nsamp),
                                 fill = "",
                                 y = "Density") +
                            scale_fill_manual(values = c("#5670BF",
                                                         "#418F6A",
                                                         "#B28542")) +
                            guides(fill = guide_legend(keyheight = 3,
                                                       keywidth = 3))  +
                            scale_x_continuous(n.breaks = 6, limits = c(0, 1))
                          p1 / p2 +  plot_layout(nrow =2,widths = c(4,1))
                         }  # withProgress bracket
            ) # close withProgress
          }) # close isolate
       # } # else
    }) # close renderPlot

     # }

    })

#    observeEvent(input$reset_input, {hide("melded_plot")})
 #   observeEvent(input$submit_melded, {message('here');show("melded_plot")})

} # end server function



#  ) # end observe event
# } # end server function
#
#
# melded %>%
#     ggplot(aes(x = value, fill = type)) +
#     geom_density(alpha = .5) +
#     facet_wrap(~name,
#                labeller = as_labeller(
#                    TeX,
#                    default = label_parsed),
#                ncol = 3) +
#     theme_bw() +
#     theme(axis.text.y = element_blank(),
#           axis.ticks.y = element_blank(),
#           axis.title = element_text(size = 18),
#           axis.text.x = element_text(size = 12),
#           plot.title =element_text(size = 25,
#                                    margin =margin(.5,.5,.5,.5)),
#           strip.text = element_text(size = 16),
#           legend.text = element_text(size = 16)) +
#     labs(title = paste0("Number of Samples: ", input$nsamp),
#         fill = "",
#         y = "Density") +
#     scale_fill_manual(values = c("#5670BF",
#                                  "#418F6A",
#                                  "#B28542")) +
#     guides(fill = guide_legend(keyheight = 2,
#                                keywidth = 2))
# # + coord_cartesian(xlim=c(0,2), ylim = c(0,2),
# #                     clip = "off")

#
# tibble(Value = sample_gamma_density(1e5,
#                                     mean = input$alpha_mean,
#                                     sd = input$alpha_sd)) %>%
#     ggplot(aes(x = Value)) +
#     geom_density(fill = "black",
#                  size = .8) +
#     # stat_function(fun = gamma_density,
#     #               geom="area",
#     #               args = list("mean" = .9,
#     #                           "sd" = .04),
#     #               xlim = c(0,5)) +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = .5, size = 16),
#           plot.subtitle = element_text(hjust = .5, size = 11),
#           axis.title = element_text(size = 16)) +
#     labs(title = latex2exp::TeX("Prior for $\\alpha$"),
#          subtitle = paste0("Mean = ",
#                            input$alpha_mean,
#                            ", Standard Deviation = ",
#                            input$alpha_sd ),
#          y = "Density",
#          x = "Value") +
#     # add point just to ensure axis starts at 0
#     geom_point(aes(x=0, y = 0), size=0) +
#     scale_x_continuous(n.breaks = 6, limits = c(0, input$alpha_mean + 8*input$alpha_sd))

# Run the application
shinyApp(ui = ui, server = server)









#############
# RESOURCES
#############
# https://stackoverflow.com/questions/70002107/how-to-build-a-shiny-app-that-shows-the-output-only-when-the-user-clicks-a-click
# strange things happen with rounding
# https://stackoverflow.com/questions/59923693/avoid-sliderinput-rounding
