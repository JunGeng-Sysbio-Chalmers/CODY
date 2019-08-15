library(shiny)
library(shinyWidgets)
Cohort_catergory <- c("","Infant", "Adult")

all_paras<-c("cohort_catergory","diet_regime","lumen_variable","mucus_variable","feces_variable","variable_3d")

ui<-shinyUI(navbarPage("CODY",
                    
                       tabPanel("Home",
                                div(img(src='Method_Fig.jpg',height=720, width=884), style="text-align: center;"),
                                br(),
                                wellPanel(
                                  div("COmupting microbial DYnamics (CODY)", style="color: #333333;font-family:sans-serif; font-weight: bold;
                                             font-size: 33px; text-align: left;"),
                                  br(),
                                  div("CODY is a computational framework that can visulize spatio-temporal scale microbial dynamic
                                                     variability across in vivo site-specific colon regions and in vitro feces. This is achieved by 
                                             three multi-scale frameworks of modeling biomimetically the large intestine physiology (A), modeling gut microbial
                                             ecosystem interactions (B), as well as species-level adaptiveness (C). CODY is currently available demonstrated
                                             for two independent longitudinal cohort studies.",style="color: #6e6e6e;font-family: helvetica; font-weight: lighter;
                                             font-size: 18px; text-align: left;margin-left:auto;margin-right:auto"),style = "overflow-y:scroll;
                                             background: #e4eff2; border-color: white;")
                                ),
                       tabPanel("ReadMe",
                                div(img(src='ReadMe_Fig.png',height=672, width=930), style="text-align: center;")
                                                       ),
                       tabPanel("Quantitative Abundance Profiles",
                                fluidRow(
                                  column(4,
                                         h4(strong("Spatiotemporal Visualization")),
                                         tags$img(src='CODY.png',height=275,width=275, align="center")
                                         ),
                                         
                                  column(8,
                                         fluidRow(
                                           
                                           column(6,selectInput('cohort_catergory', "Cohort Catergory", Cohort_catergory),selected=FALSE),

                                           
                                           column(6, actionButton("do", "Run Calculation",
                                                                  # style = "color: #1f4852; 
                                                                  style = "color: #fffdf5;
                                                                  # background-color: #4e7d8a; 
                                                                  background-color: #8da7ae; 
                                                                  height: 43px;
                                                                  width: 120px;
                                                                  font-size: 14px;
                                                                  font-family: helvetica ;font-weight: 450;
                                                                  text-align:left;"))
                                         ),
                                         h5(strong("Diet Input [g/day]")),

                                         fluidRow(
                                           column(3,uiOutput("HMO_BSL_table")),
                                           column(3,uiOutput("Mucin_BSL_table")),
                                           column(3,uiOutput("Fiber_BSL_table")),
                                           column(3,uiOutput("RS_BSL_table"))
                                         ), 
                                         fluidRow(
                                           column(3,uiOutput("HMO_Intv_table")),
                                           column(3,uiOutput("Mucin_Intv_table")),
                                           column(3,uiOutput("Fiber_Intv_table")),
                                           column(3,uiOutput("RS_Intv_table"))
                                         )
                                         
                                  )),
                                
                                  # setBackgroundColor("ghostwhite")),
                                fluidRow(       
                                  column(4,
                                         br(),
                                         br(),
                                         br(),
                                         h4(strong('Three-dimension visualization')),
                                         div(strong('Please select a variable:',style="color: #1d1d1d;font-family: helvetica; font-weight: lighter;
                                             font-size: 16px; text-align: left;")),
                                         br(),

                                         fluidRow(
                                           column(7,h5(strong('3D-Visualization'))),
                                           column(5, selectInput("variable_3d", NA,
                                                                 c("Blg","Bfr","Fpr","Acetate","Butyrate","Propionate")))
                                           ),
                                         
                                         uiOutput("plot6")),
                                         
                                         
                                  column(8,
                                         div("The calculation would proceed for 10 ~ 20 minutes, Please wait...",
                                              style="color: #6e6e6e;font-family: helvetica; font-weight: lighter;
                                             font-size: 19px; text-align: left;margin-left:auto;margin-right:auto"),
                                         h4(strong("Longitudinal Prediction of Absolute Abundance")),
                                         
                                         div(strong('Please select Colon site and Species:',style="color: #1d1d1d;font-family: helvetica; font-weight: lighter;
                                             font-size: 16px; text-align: left;")),
                                         fluidRow(
                                           column(6,h5(strong('Colon Site'))),
                                           column(6,selectInput('site', NA, c("Lumen_I","Lumen_II","Lumen_III","Lumen_IV",
                                                                              "Mucosa_I","Mucosa_II","Mucosa_III","Mucosa_IV",
                                                                              "feces","blood")))
                                         ),
                                         fluidRow(
                                           column(6,h5(strong('Species'))),
                                           column(3, selectInput("plot_variable", NA,
                                                                 c("microbial","metabolite")), selected="microbial")
                                         ),
                                         
                                         plotOutput("lineplot",height = "700px", width = "450px"),
                                  
                                         h4(strong('Longitudinal in vitro feces relative abundance predictions')),
                                         br(),
                                         plotOutput('plot4',height = "500px")
                                                    # height = "500px")))
                                         
                                  
                                  )
                                )),
                                
                                  
              
  
                  
                  tabPanel("Longitudinal Site-specific Variability Profiles",
                             column(8,
                                    h3(strong("In vivo variability across colon sites")),
                                    br(),
                                    plotOutput("Plot_site")),
                             column(4,
                                    h3(strong("Plot Control")),
                                    br(),
                                  
                                  fluidRow(
                                  column(5, selectInput("site_variable", "Site variable",
                                                          c("microbial","metabolite"))),
                                  column(7, selectInput("plot_site", "In vivo Site-specific",
                                                        c("lumen","mucosa","feces","blood")), selected="lumen")
                                  
                                  )
                                  )
                  ),
                   
                   tabPanel("More Information",   # Information about data collection.
                            br(),
                            br(),
                            "Please visit",
                            br(),
                            br(),
                            br(),
                            br(),
                            "Any questions or comments can be sent to",
                            br(),
                            "Chalmers-Sysbio-Jun Geng: " ,
                            a("gejun@chalmers.se", href="mailto:gejun@chalmers.se"),
                            br()),
                
                  tags$style(type = 'text/css', 
                             '.navbar { color: white;background-color: lightblue;  font-family: helvetica; font-size: 16px}',
                             '.navbar-default .navbar-brand{color: white;}',
                             '.tab-panel{ background-color: lightblue; color: white;font-family: helvetica; font-size: 16px}',
                             '.nav navbar-nav li.active:hover a, .nav navbar-nav li.active a {
                              color: white;
                              background-color: lightyellow;
                              border-color: yellow;
                              font-family: BentonSans Book; 
                              font-size: 50px;
                             }',
                             '.navbar-nav li a:hover, .navbar-nav > .active > a {
                             color: #fff !important;
                             background-color:#2b8cc4 !important;
                             background-image: none !important;
                             }'

                  )
                 
                  
)
)

