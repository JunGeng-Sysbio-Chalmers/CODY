library(shiny)
library(scales)
library(grid)
# devtools::install_github("schmidtchristoph/reach/reach", force = TRUE)
library(shiny)
library(readxl)
library(plyr)
library(tibble)
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(stringr)
library(ggpubr)
library(tidyverse)
library(ggsci)
library(ggplotify)
library(ggthemes)
library(ggrepel)
library(ggplot2)

library(grid)
library(gridExtra)
library(cowplot)
library(MASS)
library(ggpmisc)
library(reach)



Bac8<-read.csv("Bac8_species.txt",sep="\t", header=F, stringsAsFactors=F)
Bac8<-c("Bth", "Bfr", "Blg", "Bbv", "Bad", "Ehal", "Fpr", "Rint")

parameterNames<-c("cohort_catergory","diet_regime","lumen_variable","mucus_variable","feces_variable","variable_3d",
                  "HMO_BSL","Mucin_BSL","Fiber_BSL","RS_BSL","HMO_Intv","Mucin_Intv","Fiber_Intv","RS_Intv")



# if(file.exists("./www/Metabolites_3d.jpg"))file.remove("./www/Metabolites_3d.jpg")
# if(file.exists("./Metabolites_3d.jpg"))file.remove("./Metabolites_3d.jpg")

shinyServer(function(input, output,session) {
  v <- reactiveValues(doCalc = FALSE) ## control calculation time
  # diet_regime<-reactiveVal()
  d3_plot<-reactiveValues(plot = FALSE)
  


    
output$HMO_BSL_table<-renderUI({
  if(input$cohort_catergory=="")return()
  switch(input$cohort_catergory,
         "Infant"=return(sliderInput('HMO_BSL', 'HMO.BSL',min = 13, max = 15, step=0.1, value = 14, ticks = FALSE,round = FALSE)),
         "Adult"=return(sliderInput('HMO_BSL1', 'HMO.BSL',min = 0, max = 0, step=0, value = 0, ticks = FALSE,round = FALSE))
         )
})
  
output$Mucin_BSL_table<-renderUI({
  if(input$cohort_catergory=="")return()
  switch(input$cohort_catergory,
         "Infant"=return(sliderInput('Mucin_BSL', 'Mucin.BSL',min = 0.9, max = 1.1, step=0.05, value = 1, ticks = FALSE,round = FALSE)),
         "Adult"=return(sliderInput('Mucin_BSL1', 'Mucin.BSL',min = 1.8, max = 2.2, step=0.1, value = 2, ticks = FALSE,round = FALSE))
  )
})
output$Fiber_BSL_table<-renderUI({
  if(input$cohort_catergory=="")return()
  switch(input$cohort_catergory,
         "Infant"=return(sliderInput('Fiber_BSL', 'Fiber.BSL',min = 0, max = 0, step=0, value = 0, ticks = FALSE,round = FALSE)),
         "Adult"=return(sliderInput('Fiber_BSL1', 'Fiber.BSL',min = 30, max = 40, step=0.1, value = 33.2, ticks = FALSE,round = FALSE))
  )
})
output$RS_BSL_table<-renderUI({
  if(input$cohort_catergory=="")return()
  switch(input$cohort_catergory,
         "Infant"=return(sliderInput('RS_BSL', 'RS.BSL',min = 0, max = 0, step=0, value = 0, ticks = FALSE,round = FALSE)),
         "Adult"=return(sliderInput('RS_BSL1', 'RS.BSL',min = 15, max = 25, step=0.1, value = 23.8, ticks = FALSE,round = FALSE))
  )
})
output$HMO_Intv_table<-renderUI({
  if(input$cohort_catergory=="")return()
  switch(input$cohort_catergory,
         "Infant"=return(sliderInput('HMO_Intv', 'HMO.Intv',min = 0, max = 0, step=0, value = 0, ticks = FALSE,round = FALSE)),
         "Adult"=return(sliderInput('HMO_Intv1', 'HMO.Intv',min = 0, max = 0, step=0, value = 0, ticks = FALSE,round = FALSE))
  )
})
output$Mucin_Intv_table<-renderUI({
  if(input$cohort_catergory=="")return()
  switch(input$cohort_catergory,
         "Infant"=return(sliderInput('Mucin_Intv', 'Mucin.Intv',min = 1.45, max = 1.55, step=0.05, value = 1.5, ticks = FALSE,round = FALSE)),
         "Adult"=return(sliderInput('Mucin_Intv1', 'Mucin.Intv',min = 1.9, max = 2.1, step=0.05, value = 2, ticks = FALSE,round = FALSE))
  )
})
output$Fiber_Intv_table<-renderUI({
  if(input$cohort_catergory=="")return()
  switch(input$cohort_catergory,
         "Infant"=return(sliderInput('Fiber_Intv', 'Fiber.Intv',min = 10, max = 30, step=0.5, value = 12.5, ticks = FALSE,round = FALSE)),
         "Adult"=return(sliderInput('Fiber_Intv1', 'Fiber.Intv',min = 12, max = 50, step=0.1, value = 19.7, ticks = FALSE,round = FALSE))
  )
})
output$RS_Intv_table<-renderUI({
  if(input$cohort_catergory=="")return()
  switch(input$cohort_catergory,
         "Infant"=return(sliderInput('RS_Intv', 'RS.Intv',min = 5, max = 15, step=0.5, value = 10, ticks = FALSE,round = FALSE)),
         "Adult"=return(sliderInput('RS_Intv1', 'RS.Intv',min = 1, max = 20, step=0.1, value = 4, ticks = FALSE,round = FALSE))
  )
})

source("Draw_infant_Fig3b.R")
source("Draw_adult_Fig3e.R")
source("Draw_Fig3a_1.R")
source("Draw_Fig3d.R")
source("Tanks_cmp.R")
source("Tanks_cmp_mets.R")
source("Feces_microbial.R")
# source("Feces_metabolites.R")
source("Feces_metabolites_infant.R")
source("Plasma_metabolites_barplot.R")
source("speciesPlot_12M_Jun.R")
source("metsPlot_12M_Jun.R")
source("lineplot_infant.R")
source("lineplot_adult.R")
source("Draw_dotplot_infant.R")
source("Feces_metabolites_adult.R")

diet_regime<-reactive({
  # req(input$cohort_catergory)
  if(input$cohort_catergory=="")return(NULL)
  if(input$cohort_catergory=="Infant"){
    check_diet<-c(input$HMO_BSL,input$Mucin_BSL,input$Fiber_BSL,input$RS_BSL,input$HMO_Intv,input$Mucin_Intv,input$Fiber_Intv,input$RS_Intv)
    print(check_diet)
    if(length(check_diet)<8)return(NULL)
    # return(rbind(c(hmo_bsl,mucin_bsl,fiber_bsl,rs_bsl),c(hmo_intv,mucin_intv,fiber_intv,rs_intv)))
    # return(as.matrix(check_diet,nrow=2))
    return(rbind(check_diet[1:4],check_diet[5:8]))
  }
  if(input$cohort_catergory=="Adult"){
    check_diet<-c(input$HMO_BSL1,input$Mucin_BSL1,input$Fiber_BSL1,input$RS_BSL1,input$HMO_Intv1,input$Mucin_Intv1,input$Fiber_Intv1,input$RS_Intv1)
    print(check_diet)
    if(length(check_diet)<8)return(NULL)
    # return(rbind(c(hmo_bsl,mucin_bsl,fiber_bsl,rs_bsl),c(hmo_intv,mucin_intv,fiber_intv,rs_intv)))
    # return(as.matrix(check_diet,nrow=2))
    return(rbind(check_diet[1:4],check_diet[5:8]))
  }
  
})


plot_3d_index<-reactive({
  if(input$cohort_catergory==""||input$variable_3d==""||is.null(diet_regime()))return(NULL)
  if(input$cohort_catergory=="Infant"){
    switch(input$variable_3d,
           "Blg" = return(3),
           "Bfr" = return(2),
           "Fpr" = return(7),
           "Acetate" = return(12),
           "Butyrate" = return(17),
           "Propionate" = return(13)
    )
  }
  
  if(input$cohort_catergory=="Adult"){
    switch(input$variable_3d,
           "Blg" = return(2),
           "Bfr" = return(1),
           "Fpr" = return(6),
           "Acetate" = return(11),
           "Butyrate" = return(16),
           "Propionate" = return(12)
    )
  }
})


observeEvent(input$do,{
  # req(diet_regime(),input$cohort_catergory)
  v$doCalc<-input$do
  })
  
observeEvent(input$cohort_catergory,{
  v$doCalc <- FALSE
})  

observe({
  # req(input$cohort_catergory,diet_regime())
  req(v$doCalc)
  req(diet_regime())
  # print(diet_regime())
  # print(input$cohort_catergory)
  
  withProgress(message = "Calculation in progress",
  detail = "This may take 10 ~ 20 minutes...", {
               diet_input<-diet_regime()
               cohort_info<-input$cohort_catergory
               results1<-runMatlabFct("restul_total=run_Simulation(diet_input,cohort_info)")
               
             })
})


  


  
  exp_data<-reactive({
    if(input$cohort_catergory==""||is.null(diet_regime()))return(NULL)
    if(input$cohort_catergory=="Infant"){
      return(read.table("68_BF_sample_8spc_RA.csv",sep="\t", header=T, stringsAsFactors=F ))
    }
    if(input$cohort_catergory=="Adult"){
      read.table("Mean_RA_species.txt",sep="\t", header=T, stringsAsFactors=F )
     }
  })
  
  prediction_data<-reactive({
    if(input$cohort_catergory==""||is.null(diet_regime())||v$doCalc=="")return(NULL)
    if(input$cohort_catergory=="Infant"){
      return(read.table(file="Preidiction_RA_exp_infant.csv", sep=";", header=T, stringsAsFactors=F ))
    }
    if(input$cohort_catergory=="Adult"){
      return(read.table(file="prediction_format.csv", sep=",", header=T, stringsAsFactors=F))
    }
  })
  
    
  col_8<-c(brewer.pal(n=8,"Dark2")[1:5],c("#BC80BD", "#E6AB02", "#80B1D3"))
  
  Colon_Site<-c("I","II","III","IV")
  

  output$plot4<-renderPlot({
    req(v$doCalc)
    
    if(input$cohort_catergory==""||is.null(diet_regime()))return(NULL)
    if(is.null(prediction_data()))return(NULL)
    if(input$cohort_catergory=="Infant"){
      exp_input<-exp_data()
      pre_input_infant<-prediction_data()
      plot4<-Draw_infant_Fig3b(exp_input,pre_input_infant)
      return(plot4)
      
    }    ##
    if(input$cohort_catergory=="Adult"){
      exp_input<-exp_data()
      pre_input_adult<-prediction_data()
      plot4<-Draw_adult_Fig3e(exp_input,pre_input_adult)
      return(plot4)
      
    }
  })
  
  # # Validation 2, pie plot for infant and dot plot for adult
  # output$plot5<-renderPlot({
  #   req(v$doCalc)
  #   
  #   if(is.null(prediction_data()))return(NULL)
  #   if(input$cohort_catergory=="Infant"){
  #     prediction_exp_file<-"Preidiction_RA_exp_mdf.csv"
  #     pre_input<-prediction_data()
  #    
  #     plot5<-Draw_dotplot_infant(pre_input)         # can save as png and then insert figures later, if this does not work well
  #     
  #     return(plot5)
  #     
  #   }    ##
  #   if(input$cohort_catergory=="Adult"){
  #     exp_input<-"Bac7_sample_ra.txt"
  #     pre_input<-prediction_data()
  #     print(pre_input)
  #     plot5<-Draw_Fig3d(exp_input,pre_input)
  #     return(plot5)
  #     
  #   }
  # })
  
  
 
  data_pre_lumen1<-reactive({
    req(v$doCalc)
    if(input$cohort_catergory==""||is.null(diet_regime()))return(NULL)
    switch(input$cohort_catergory,
           "Infant" = return(read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Lumen_1.csv"))),
           "Adult" = return(read.csv(paste0("Community_10tank_adult_concate_profiles_","Lumen_1.csv")))
    )
      })

  data_pre_lumen2<-reactive({
    req(v$doCalc)
    
    if(input$cohort_catergory==""||is.null(diet_regime()))return(NULL)
    switch(input$cohort_catergory,
           "Infant" = return(read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Lumen_2.csv"))),
           "Adult" = return(read.csv(paste0("Community_10tank_adult_concate_profiles_","Lumen_2.csv"))))
  })
  data_pre_lumen3<-reactive({
    req(v$doCalc)
    
    if(input$cohort_catergory==""||is.null(diet_regime()))return(NULL)
    switch(input$cohort_catergory,
           "Infant" = return(read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Lumen_3.csv"))),
           "Adult" = return(read.csv(paste0("Community_10tank_adult_concate_profiles_","Lumen_3.csv"))))
  })
  data_pre_lumen4<-reactive({
    req(v$doCalc)
    
    if(input$cohort_catergory==""||is.null(diet_regime()))return(NULL)
    switch(input$cohort_catergory,
           "Infant" = return(read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Lumen_4.csv"))),
           "Adult" = return(read.csv(paste0("Community_10tank_adult_concate_profiles_","Lumen_4.csv"))))
  })
  data_pre_Mucosa1<-reactive({
    req(v$doCalc)
    
    if(input$cohort_catergory==""||is.null(diet_regime()))return(NULL)
    switch(input$cohort_catergory,
           "Infant" = return(read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Mucosa_1.csv"))),
           "Adult" = return(read.csv(paste0("Community_10tank_adult_concate_profiles_","Mucosa_1.csv"))))
  })
  data_pre_Mucosa2<-reactive({
    req(v$doCalc)
    
    if(input$cohort_catergory==""||is.null(diet_regime()))return(NULL)
    switch(input$cohort_catergory,
           "Infant" = return(read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Mucosa_2.csv"))),
           "Adult" = return(read.csv(paste0("Community_10tank_adult_concate_profiles_","Mucosa_2.csv"))))
  })
  data_pre_Mucosa3<-reactive({
    req(v$doCalc)
    
    if(input$cohort_catergory==""||is.null(diet_regime()))return(NULL)
    switch(input$cohort_catergory,
           "Infant" = return(read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Mucosa_3.csv"))),
           "Adult" = return(read.csv(paste0("Community_10tank_adult_concate_profiles_","Mucosa_3.csv"))))
  })
  data_pre_Mucosa4<-reactive({
    req(v$doCalc)
    
    if(input$cohort_catergory==""||is.null(diet_regime()))return(NULL)
    switch(input$cohort_catergory,
           "Infant" = return(read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Mucosa_4.csv"))),
           "Adult" = return(read.csv(paste0("Community_10tank_adult_concate_profiles_","Mucosa_4.csv"))))
  })
  
  
  data_pre_blood<-reactive({
    req(v$doCalc)
    
    if(input$cohort_catergory==""||is.null(diet_regime()))return(NULL)
    switch(input$cohort_catergory,
           "Infant" = return(read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Blood.csv"))),
           "Adult" = return(read.csv(paste0("Community_10tank_adult_concate_profiles_","Blood.csv"))))
  })
  data_pre_feces<-reactive({
    req(v$doCalc)
    
    if(input$cohort_catergory==""||is.null(diet_regime()))return(NULL)
    switch(input$cohort_catergory,
           "Infant" = return(read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Feces.csv"))),
           "Adult" = return(read.csv(paste0("Community_10tank_adult_concate_profiles_","Feces.csv"))))
  })


  # Microbial RA in the feces
  # deal with absolute abundance profiles of lumen-spatial data, mucosa-spatial data
  Time_str_BF<-294.9   ## endpoint of feces
  Time_str_SF<-1032
  Time_str_BSL<-588  ## endpoint of feces, for adult
  
  Time_str_Intv<-1032
  Bac8
  nl<-reactive({
    req(v$doCalc)
    
    if(input$cohort_catergory==""||is.null(diet_regime()))return(NULL)
    ncol(data_pre_lumen1())
  })

  BF_index<-reactive({
    req(v$doCalc)
    
    if(input$cohort_catergory==""||is.null(diet_regime()))return(NULL)
      data_pre_lumen1_input<-data_pre_lumen1()
      return(match(sort(data_pre_lumen1_input$Time[data_pre_lumen1_input$Time-Time_str_BSL>=-1e-20])[1],data_pre_lumen1_input$Time)-1)
  })
  SF_index<-reactive({
    req(v$doCalc)
    
    if(input$cohort_catergory==""||is.null(diet_regime()))return(NULL)
      data_pre_lumen1_input<-data_pre_lumen1()
      return(nrow(data_pre_lumen1_input))
  })
  
  
  data_lumen_microbial_all_BF<-reactive({
    req(v$doCalc)
    
    if(input$cohort_catergory==""||is.null(diet_regime()))return(NULL)
    if(input$cohort_catergory=="Infant"){
      data_1<-read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Lumen_1.csv"))
      data_2<-read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Lumen_2.csv"))
      data_3<-read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Lumen_3.csv"))
      data_4<-read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Lumen_4.csv"))
      BFIndex<-match(sort(data_1$Time[data_1$Time-Time_str_BF>=-1e-20])[1],data_1$Time)-1
      BF_ncol<-ncol(data_1)
      # data_lumen_microbial_all_BF<-rbind(data_1[BFIndex,2:9],data_2[BFIndex,2:9],data_3[BFIndex,2:9],data_4[BFIndex,2:9])
      return(rbind(data_1[BFIndex,2:9],data_2[BFIndex,2:9],data_3[BFIndex,2:9],data_4[BFIndex,2:9]))
      }
    if(input$cohort_catergory=="Adult"){
      # data_1<-read_excel(path = "Community_10tanks_output_Feces_0to12m_68sample_20181031.xlsx", sheet = "Lumen_1")
      data_1<-read.csv(paste0("Community_10tank_adult_concate_profiles_","Lumen_1.csv"))
      data_2<-read.csv(paste0("Community_10tank_adult_concate_profiles_","Lumen_2.csv"))
      data_3<-read.csv(paste0("Community_10tank_adult_concate_profiles_","Lumen_3.csv"))
      data_4<-read.csv(paste0("Community_10tank_adult_concate_profiles_","Lumen_4.csv"))
      BFIndex<-match(sort(data_1$Time[data_1$Time-Time_str_BSL>=-1e-20])[1],data_1$Time)-1
      BF_ncol<-ncol(data_1)
      return(rbind(data_1[BFIndex,2:8],data_2[BFIndex,2:8],data_3[BFIndex,2:8],data_4[BFIndex,2:8]))
    }
  })
  
  data_lumen_metabolite_all_BF<-reactive({
    req(v$doCalc)
    
    if(input$cohort_catergory==""||is.null(diet_regime()))return(NULL)
    if(input$cohort_catergory=="Infant"){
      data_1<-read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Lumen_1.csv"))
      data_2<-read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Lumen_2.csv"))
      data_3<-read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Lumen_3.csv"))
      data_4<-read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Lumen_4.csv"))
      BFIndex<-match(sort(data_1$Time[data_1$Time-Time_str_BF>=-1e-20])[1],data_1$Time)-1
      BF_ncol<-ncol(data_1)
      # data_lumen_metabolite_all_BF<-rbind(data_1[BFIndex,10:BF_ncol],data_2[BFIndex,10:BF_ncol],data_3[BFIndex,10:BF_ncol],data_4[BFIndex,10:BF_ncol])
      return(rbind(data_1[BFIndex,10:BF_ncol],data_2[BFIndex,10:BF_ncol],data_3[BFIndex,10:BF_ncol],data_4[BFIndex,10:BF_ncol]))
    }
    if(input$cohort_catergory=="Adult"){
      data_1<-read.csv(paste0("Community_10tank_adult_concate_profiles_","Lumen_1.csv"))
      data_2<-read.csv(paste0("Community_10tank_adult_concate_profiles_","Lumen_2.csv"))
      data_3<-read.csv(paste0("Community_10tank_adult_concate_profiles_","Lumen_3.csv"))
      data_4<-read.csv(paste0("Community_10tank_adult_concate_profiles_","Lumen_4.csv"))
      BFIndex<-match(sort(data_1$Time[data_1$Time-Time_str_BSL>=-1e-20])[1],data_1$Time)-1
      BF_ncol<-ncol(data_1)
      return(rbind(data_1[BFIndex,9:BF_ncol],data_2[BFIndex,9:BF_ncol],data_3[BFIndex,9:BF_ncol],data_4[BFIndex,9:BF_ncol]))
    }
  })
  data_mucosa_microbial_all_BF<-reactive({
    req(v$doCalc)
    
    if(input$cohort_catergory==""||is.null(diet_regime()))return(NULL)
    if(input$cohort_catergory=="Infant"){
      data_1<-read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Mucosa_1.csv"))
      data_2<-read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Mucosa_2.csv"))
      data_3<-read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Mucosa_3.csv"))
      data_4<-read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Mucosa_4.csv"))
      BFIndex<-match(sort(data_1$Time[data_1$Time-Time_str_BF>=-1e-20])[1],data_1$Time)-1
      BF_ncol<-ncol(data_1)
      return(rbind(data_1[BFIndex,2:9],data_2[BFIndex,2:9],data_3[BFIndex,2:9],data_4[BFIndex,2:9]))
    }
    if(input$cohort_catergory=="Adult"){
      data_1<-read.csv(paste0("Community_10tank_adult_concate_profiles_","Mucosa_1.csv"))
      data_2<-read.csv(paste0("Community_10tank_adult_concate_profiles_","Mucosa_2.csv"))
      data_3<-read.csv(paste0("Community_10tank_adult_concate_profiles_","Mucosa_3.csv"))
      data_4<-read.csv(paste0("Community_10tank_adult_concate_profiles_","Mucosa_4.csv"))
      BFIndex<-match(sort(data_1$Time[data_1$Time-Time_str_BSL>=-1e-20])[1],data_1$Time)-1
      BF_ncol<-ncol(data_1)
      return(rbind(data_1[BFIndex,2:8],data_2[BFIndex,2:8],data_3[BFIndex,2:8],data_4[BFIndex,2:8]))
    }
  })
  
  data_mucosa_metabolite_all_BF<-reactive({
    req(v$doCalc)
    
    if(input$cohort_catergory==""||is.null(diet_regime()))return(NULL)
    if(input$cohort_catergory=="Infant"){
      data_1<-read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Mucosa_1.csv"))
      data_2<-read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Mucosa_2.csv"))
      data_3<-read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Mucosa_3.csv"))
      data_4<-read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Mucosa_4.csv"))
      BFIndex<-match(sort(data_1$Time[data_1$Time-Time_str_BF>=-1e-20])[1],data_1$Time)-1
      BF_ncol<-ncol(data_1)
      return(rbind(data_1[BFIndex,10:BF_ncol],data_2[BFIndex,10:BF_ncol],data_3[BFIndex,10:BF_ncol],data_4[BFIndex,10:BF_ncol]))
    }
    if(input$cohort_catergory=="Adult"){
      data_1<-read.csv(paste0("Community_10tank_adult_concate_profiles_","Mucosa_1.csv"))
      data_2<-read.csv(paste0("Community_10tank_adult_concate_profiles_","Mucosa_2.csv"))
      data_3<-read.csv(paste0("Community_10tank_adult_concate_profiles_","Mucosa_3.csv"))
      data_4<-read.csv(paste0("Community_10tank_adult_concate_profiles_","Mucosa_4.csv"))
      BFIndex<-match(sort(data_1$Time[data_1$Time-Time_str_BSL>=-1e-20])[1],data_1$Time)-1
      BF_ncol<-ncol(data_1)
      return(rbind(data_1[BFIndex,9:BF_ncol],data_2[BFIndex,9:BF_ncol],data_3[BFIndex,9:BF_ncol],data_4[BFIndex,9:BF_ncol]))
    }
  })

  
  
  data_lumen_microbial_all_SF<-reactive({
    req(v$doCalc)
    
    if(input$cohort_catergory==""||is.null(diet_regime()))return(NULL)
    if(input$cohort_catergory=="Infant"){
      data_1<-read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Lumen_1.csv"))
      data_2<-read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Lumen_2.csv"))
      data_3<-read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Lumen_3.csv"))
      data_4<-read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Lumen_4.csv"))
      SFIndex<-nrow(data_1)
      SF_ncol<-ncol(data_1)
      # data_lumen_microbial_all_SF<-rbind(data_1[SFIndex,2:9],data_2[SFIndex,2:9],data_3[SFIndex,2:9],data_4[SFIndex,2:9])
      return(rbind(data_1[SFIndex,2:9],data_2[SFIndex,2:9],data_3[SFIndex,2:9],data_4[SFIndex,2:9]))
    }
    if(input$cohort_catergory=="Adult"){
      data_1<-read.csv(paste0("Community_10tank_adult_concate_profiles_","Lumen_1.csv"))
      data_2<-read.csv(paste0("Community_10tank_adult_concate_profiles_","Lumen_2.csv"))
      data_3<-read.csv(paste0("Community_10tank_adult_concate_profiles_","Lumen_3.csv"))
      data_4<-read.csv(paste0("Community_10tank_adult_concate_profiles_","Lumen_4.csv"))
      SFIndex<-nrow(data_1)
      SF_ncol<-ncol(data_1)
      return(rbind(data_1[SFIndex,2:8],data_2[SFIndex,2:8],data_3[SFIndex,2:8],data_4[SFIndex,2:8]))
    }
  })
  
  
  data_lumen_metabolite_all_SF<-reactive({
    req(v$doCalc)
    
    if(input$cohort_catergory==""||is.null(diet_regime()))return(NULL)
    if(input$cohort_catergory=="Infant"){
      data_1<-read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Lumen_1.csv"))
      data_2<-read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Lumen_2.csv"))
      data_3<-read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Lumen_3.csv"))
      data_4<-read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Lumen_4.csv"))
      SFIndex<-nrow(data_1)
      SF_ncol<-ncol(data_1)
      # data_lumen_metabolite_all_SF<-rbind(data_1[SFIndex,10:SF_ncol],data_2[SFIndex,10:SF_ncol],data_3[SFIndex,10:SF_ncol],data_4[SFIndex,10:SF_ncol])
      return(rbind(data_1[SFIndex,10:SF_ncol],data_2[SFIndex,10:SF_ncol],data_3[SFIndex,10:SF_ncol],data_4[SFIndex,10:SF_ncol]))
    }
    if(input$cohort_catergory=="Adult"){
      data_1<-read.csv(paste0("Community_10tank_adult_concate_profiles_","Lumen_1.csv"))
      data_2<-read.csv(paste0("Community_10tank_adult_concate_profiles_","Lumen_2.csv"))
      data_3<-read.csv(paste0("Community_10tank_adult_concate_profiles_","Lumen_3.csv"))
      data_4<-read.csv(paste0("Community_10tank_adult_concate_profiles_","Lumen_4.csv"))
      SFIndex<-nrow(data_1)
      SF_ncol<-ncol(data_1)
      return(rbind(data_1[SFIndex,9:SF_ncol],data_2[SFIndex,9:SF_ncol],data_3[SFIndex,9:SF_ncol],data_4[SFIndex,9:SF_ncol]))
    }
  })
  
  
  data_mucosa_microbial_all_SF<-reactive({
    req(v$doCalc)
    
    if(input$cohort_catergory==""||is.null(diet_regime()))return(NULL)
    if(input$cohort_catergory=="Infant"){
      data_1<-read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Mucosa_1.csv"))
      data_2<-read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Mucosa_2.csv"))
      data_3<-read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Mucosa_3.csv"))
      data_4<-read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Mucosa_4.csv"))
      SFIndex<-nrow(data_1)
      SF_ncol<-ncol(data_1)
      return(rbind(data_1[SFIndex,2:9],data_2[SFIndex,2:9],data_3[SFIndex,2:9],data_4[SFIndex,2:9]))
    }
    if(input$cohort_catergory=="Adult"){
      data_1<-read.csv(paste0("Community_10tank_adult_concate_profiles_","Mucosa_1.csv"))
      data_2<-read.csv(paste0("Community_10tank_adult_concate_profiles_","Mucosa_2.csv"))
      data_3<-read.csv(paste0("Community_10tank_adult_concate_profiles_","Mucosa_3.csv"))
      data_4<-read.csv(paste0("Community_10tank_adult_concate_profiles_","Mucosa_4.csv"))
      SFIndex<-nrow(data_1)
      SF_ncol<-ncol(data_1)
      return(rbind(data_1[SFIndex,2:8],data_2[SFIndex,2:8],data_3[SFIndex,2:8],data_4[SFIndex,2:8]))
    }
  })
  
  

  data_mucosa_metabolite_all_SF<-reactive({
    req(v$doCalc)
    
    if(input$cohort_catergory==""||is.null(diet_regime()))return(NULL)
    if(input$cohort_catergory=="Infant"){
      data_1<-read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Mucosa_1.csv"))
      data_2<-read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Mucosa_2.csv"))
      data_3<-read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Mucosa_3.csv"))
      data_4<-read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Mucosa_4.csv"))
      SFIndex<-nrow(data_1)
      SF_ncol<-ncol(data_1)
      return(rbind(data_1[SFIndex,10:SF_ncol],data_2[SFIndex,10:SF_ncol],data_3[SFIndex,10:SF_ncol],data_4[SFIndex,10:SF_ncol]))
    }
    if(input$cohort_catergory=="Adult"){
      data_1<-read.csv(paste0("Community_10tank_adult_concate_profiles_","Mucosa_1.csv"))
      data_2<-read.csv(paste0("Community_10tank_adult_concate_profiles_","Mucosa_2.csv"))
      data_3<-read.csv(paste0("Community_10tank_adult_concate_profiles_","Mucosa_3.csv"))
      data_4<-read.csv(paste0("Community_10tank_adult_concate_profiles_","Mucosa_4.csv"))
      SFIndex<-nrow(data_1)
      SF_ncol<-ncol(data_1)
      return(rbind(data_1[SFIndex,9:SF_ncol],data_2[SFIndex,9:SF_ncol],data_3[SFIndex,9:SF_ncol],data_4[SFIndex,9:SF_ncol]))
    }
  })
  
  
  data_feces_microbial_all<-reactive({
    req(v$doCalc)
    
    if(input$cohort_catergory==""||is.null(diet_regime()))return(NULL)
    if(input$cohort_catergory=="Infant"){
      data_1<-read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Feces.csv"))
      BFIndex<-match(sort(data_1$Time[data_1$Time-Time_str_BF>=-1e-20])[1],data_1$Time)-1
      SFIndex<-nrow(data_1)
      SF_ncol<-ncol(data_1)
      # data_feces_microbial_all<-rbind(data_1[BFIndex,2:9]*19.5,data_1[SFIndex,2:9]*6.5)
      return(rbind(data_1[BFIndex,2:9]*19.5,data_1[SFIndex,2:9]*6.5))
    }
    if(input$cohort_catergory=="Adult"){
      data_1<-read.csv(paste0("Community_10tank_adult_concate_profiles_","Feces.csv"))
      
      BFIndex<-match(sort(data_1$Time[data_1$Time-Time_str_BSL>=-1e-20])[1],data_1$Time)-1
      
      SFIndex<-nrow(data_1)
      SF_ncol<-ncol(data_1)
      # pre_Data<-rbind(data_1[BFIndex,2:8]*3.3,data_1[SFIndex,2:8]*3.3)
      return(rbind(data_1[BFIndex,2:8]*3.3,data_1[SFIndex,2:8]*3.3))
    }
  })
  

  
  Time_vector<-reactive({
    req(v$doCalc)
    
    if(input$cohort_catergory==""||is.null(diet_regime()))return(NULL)
    if(input$cohort_catergory=="Infant")
      return(c("Breastfeeding","Solidfood"))
    if(input$cohort_catergory=="Adult")
      return(c("Baseline","Intervention"))
  })
  
  
  data_feces_metabolite_all<-reactive({
    req(v$doCalc)
    
    if(input$cohort_catergory==""||is.null(diet_regime()))return(NULL)
    if(input$cohort_catergory=="Infant"){
      data_1<-read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Feces.csv"))
      
      BFIndex<-match(sort(data_1$Time[data_1$Time-Time_str_BF>=-1e-20])[1],data_1$Time)-1
      SFIndex<-nrow(data_1)
      SF_ncol<-ncol(data_1)
      # data_feces_metabolite_all<-rbind(data_1[BFIndex,10:SF_ncol]*19.5,data_1[SFIndex,10:SF_ncol]*6.5)
      return(rbind(data_1[BFIndex,10:SF_ncol]*19.5,data_1[SFIndex,10:SF_ncol]*6.5))
    }
    if(input$cohort_catergory=="Adult"){
      data_1<-read.csv(paste0("Community_10tank_adult_concate_profiles_","Feces.csv"))
      
      BFIndex<-match(sort(data_1$Time[data_1$Time-Time_str_BSL>=-1e-20])[1],data_1$Time)-1
      
      SFIndex<-nrow(data_1)
      SF_ncol<-ncol(data_1)
      return(rbind(data_1[BFIndex,9:SF_ncol]*3.3,data_1[SFIndex,9:SF_ncol]*3.3))
    }
  })
  
  
  
  data_blood_metabolite_all<-reactive({
    req(v$doCalc)
    
    if(input$cohort_catergory==""||is.null(diet_regime()))return(NULL)
    if(input$cohort_catergory=="Infant"){
      data_1<-read.csv(paste0("Community_10tanks_infant_BF_SF_new_","Blood.csv"))
      
      BFIndex<-match(sort(data_1$Time[data_1$Time-Time_str_BF>=-1e-20])[1],data_1$Time)-1
      SFIndex<-nrow(data_1)
      SF_ncol<-ncol(data_1)
      return(rbind(data_1[BFIndex,2:SF_ncol],data_1[SFIndex,2:SF_ncol]))
    }
    if(input$cohort_catergory=="Adult"){
      data_1<-read.csv(paste0("Community_10tank_adult_concate_profiles_","Blood.csv"))
      
      BFIndex<-match(sort(data_1$Time[data_1$Time-Time_str_BSL>=-1e-20])[1],data_1$Time)-1
      
      SFIndex<-nrow(data_1)
      SF_ncol<-ncol(data_1)
      return(rbind(data_1[BFIndex,2:SF_ncol],data_1[SFIndex,2:SF_ncol]))
    }
  })
  
  
  
  output$lineplot<-renderPlot({
    req(v$doCalc)
    
    req(input$site)
    if(input$site_variable==""||input$plot_site==""||input$cohort_catergory==""||is.null(diet_regime()))return(NULL)
    if(input$cohort_catergory=="Infant"&&input$site=="Lumen_I"&&input$plot_variable=="microbial"){
      plot1<-lineplot_infant(data_pre_lumen1(),"microbial",1:9)
      return(plot1)
    }
    if(input$cohort_catergory=="Infant"&&input$site=="Lumen_I"&&input$plot_variable=="metabolite"){
      plot1<-lineplot_infant(data_pre_lumen1(),"metabolite",c(1, 10:19))
      # plot1<-lineplot_infant(data_1,"metabolite",c(1, 10:19))
      return(plot1)
    }
    if(input$cohort_catergory=="Infant"&&input$site=="Lumen_II"&&input$plot_variable=="microbial"){
      plot1<-lineplot_infant(data_pre_lumen2(),"microbial",1:9)
      return(plot1)
    }
      if(input$cohort_catergory=="Infant"&&input$site=="Lumen_II"&&input$plot_variable=="metabolite")
      {
        plot1<-lineplot_infant(data_pre_lumen2(),"metabolite",c(1, 10:19))
        return(plot1)
      }
 
      if(input$cohort_catergory=="Infant"&&input$site=="Lumen_III"&&input$plot_variable=="microbial")
      {
        plot1<-lineplot_infant(data_pre_lumen3(),"microbial",1:9)
        return(plot1)
      }
   
      if(input$cohort_catergory=="Infant"&&input$site=="Lumen_III"&&input$plot_variable=="metabolite")
      {
        plot1<-lineplot_infant(data_pre_lumen3(),"metabolite",c(1, 10:19))
        return(plot1)
      }
      if(input$cohort_catergory=="Infant"&&input$site=="Lumen_IV"&&input$plot_variable=="microbial")
      {
        plot1<-lineplot_infant(data_pre_lumen4(),"microbial",1:9)
        return(plot1)
      }
      if(input$cohort_catergory=="Infant"&&input$site=="Lumen_IV"&&input$plot_variable=="metabolite")
      {
        plot1<-lineplot_infant(data_pre_lumen4(),"metabolite",c(1, 10:19))
        return(plot1)
      }
      if(input$cohort_catergory=="Infant"&&input$site=="Mucosa_I"&&input$plot_variable=="microbial")
      {
        plot1<-lineplot_infant(data_pre_Mucosa1(),"microbial",1:9)
        return(plot1)
      }
      if(input$cohort_catergory=="Infant"&&input$site=="Mucosa_I"&&input$plot_variable=="metabolite")
      {
        plot1<-lineplot_infant(data_pre_Mucosa1(),"metabolite",c(1, 10:19))
        return(plot1)
      }
      if(input$cohort_catergory=="Infant"&&input$site=="Mucosa_II"&&input$plot_variable=="microbial")
      {
        plot1<-lineplot_infant(data_pre_Mucosa2(),"microbial",1:9)
        return(plot1)
      }
      if(input$cohort_catergory=="Infant"&&input$site=="Mucosa_II"&&input$plot_variable=="metabolite")
      {
        plot1<-lineplot_infant(data_pre_Mucosa2(),"metabolite",c(1, 10:19))
        return(plot1)
      }
      if(input$cohort_catergory=="Infant"&&input$site=="Mucosa_III"&&input$plot_variable=="microbial")
      {
        plot1<-lineplot_infant(data_pre_Mucosa3(),"microbial",1:9)
        return(plot1)
      }
      if(input$cohort_catergory=="Infant"&&input$site=="Mucosa_III"&&input$plot_variable=="metabolite")
      {
        plot1<-lineplot_infant(data_pre_Mucosa3(),"metabolite",c(1, 10:19))
        return(plot1)
      }
      if(input$cohort_catergory=="Infant"&&input$site=="Mucosa_IV"&&input$plot_variable=="microbial")
      {
        plot1<-lineplot_infant(data_pre_Mucosa4(),"microbial",1:9)
        return(plot1)
      }
      if(input$cohort_catergory=="Infant"&&input$site=="Mucosa_IV"&&input$plot_variable=="metabolite")
      {
        plot1<-lineplot_infant(data_pre_Mucosa4(),"metabolite",c(1, 10:19))
        return(plot1)
      }
      if(input$cohort_catergory=="Infant"&&input$site=="feces"&&input$plot_variable=="microbial")
      {
        plot1<-lineplot_infant(data_pre_feces(),"microbial",1:9)
        return(plot1)
      }
      if(input$cohort_catergory=="Infant"&&input$site=="feces"&&input$plot_variable=="metabolite")
      {
        plot1<-lineplot_infant(data_pre_feces(),"metabolite",c(1, 10:19))
        return(plot1)
      }
    if(input$cohort_catergory=="Infant"&&input$site=="blood"&&input$plot_variable=="microbial") return("No Microbial in Blood")
      if(input$cohort_catergory=="Infant"&&input$site=="blood"&&input$plot_variable=="metabolite")
      {
        plot1<-lineplot_infant(data_pre_blood(),"metabolite",1:11)
        return(plot1)
      }
    
      if(input$cohort_catergory=="Adult"&&input$site=="Lumen_I"&&input$plot_variable=="microbial")
      {
        plot1<-lineplot_adult(data_pre_lumen1(),"microbial",1:8)
        return(plot1)
      }
      if(input$cohort_catergory=="Adult"&&input$site=="Lumen_I"&&input$plot_variable=="metabolite")
      {
        plot1<-lineplot_adult(data_pre_lumen1(),"metabolite",c(1,9:18))
        return(plot1)
      }
      if(input$cohort_catergory=="Adult"&&input$site=="Lumen_II"&&input$plot_variable=="microbial")
      {
        plot1<-lineplot_adult(data_pre_lumen2(),"microbial",1:8)
        return(plot1)
      }
      if(input$cohort_catergory=="Adult"&&input$site=="Lumen_II"&&input$plot_variable=="metabolite")
      {
        plot1<-lineplot_adult(data_pre_lumen2(),"metabolite",c(1,9:18))
        return(plot1)
      }
      if(input$cohort_catergory=="Adult"&&input$site=="Lumen_III"&&input$plot_variable=="microbial")
      {
        plot1<-lineplot_adult(data_pre_lumen3(),"microbial",1:8)
        return(plot1)
      }
      if(input$cohort_catergory=="Adult"&&input$site=="Lumen_III"&&input$plot_variable=="metabolite")
      {
        plot1<-lineplot_adult(data_pre_lumen3(),"metabolite",c(1,9:18))
        return(plot1)
      }
      if(input$cohort_catergory=="Adult"&&input$site=="Lumen_IV"&&input$plot_variable=="microbial")
      {
        plot1<-lineplot_adult(data_pre_lumen4(),"microbial",1:8)
        return(plot1)
      }
      if(input$cohort_catergory=="Adult"&&input$site=="Lumen_IV"&&input$plot_variable=="metabolite")
      {
        plot1<-lineplot_adult(data_pre_lumen4(),"metabolite",c(1,9:18))
        return(plot1)
      }
      if(input$cohort_catergory=="Adult"&&input$site=="Mucosa_I"&&input$plot_variable=="microbial")
      {
        plot1<-lineplot_adult(data_pre_Mucosa1(),"microbial",1:8)
        return(plot1)
      }
      if(input$cohort_catergory=="Adult"&&input$site=="Mucosa_I"&&input$plot_variable=="metabolite")
      {
        plot1<-lineplot_adult(data_pre_Mucosa1(),"metabolite",c(1,9:18))
        return(plot1)
      }
      if(input$cohort_catergory=="Adult"&&input$site=="Mucosa_II"&&input$plot_variable=="microbial")
      {
        plot1<-lineplot_adult(data_pre_Mucosa2(),"microbial",1:8)
        return(plot1)
      }
      if(input$cohort_catergory=="Adult"&&input$site=="Mucosa_II"&&input$plot_variable=="metabolite")
      {
        plot1<-lineplot_adult(data_pre_Mucosa2(),"metabolite",c(1,9:18))
        return(plot1)
      }
      if(input$cohort_catergory=="Adult"&&input$site=="Mucosa_III"&&input$plot_variable=="microbial")
      {
        plot1<-lineplot_adult(data_pre_Mucosa3(),"microbial",1:8)
        return(plot1)
      }
      if(input$cohort_catergory=="Adult"&&input$site=="Mucosa_III"&&input$plot_variable=="metabolite")
      {
        plot1<-lineplot_adult(data_pre_Mucosa3(),"metabolite",c(1,9:18))
        return(plot1)
      }
      if(input$cohort_catergory=="Adult"&&input$site=="Mucosa_IV"&&input$plot_variable=="microbial")
      {
        plot1<-lineplot_adult(data_pre_Mucosa4(),"microbial",1:8)
        return(plot1)
      }
      if(input$cohort_catergory=="Adult"&&input$site=="Mucosa_IV"&&input$plot_variable=="metabolite")
      {
        plot1<-lineplot_adult(data_pre_Mucosa4(),"metabolite",c(1,9:18))
        return(plot1)
      }
      if(input$cohort_catergory=="Adult"&&input$site=="feces"&&input$plot_variable=="microbial")
      {
        plot1<-lineplot_adult(data_pre_feces(),"microbial",1:8)
        return(plot1)
      }
      if(input$cohort_catergory=="Adult"&&input$site=="feces"&&input$plot_variable=="metabolite")
      {

        plot1<-lineplot_adult(data_pre_feces(),"metabolite",c(1,9:18))
        return(plot1)
      }
    if(input$cohort_catergory=="Adult"&&input$site=="blood"&&input$plot_variable=="microbial") return("No Microbial in Blood")
    
      if(input$cohort_catergory=="Adult"&&input$site=="blood"&&input$plot_variable=="metabolite")
      {
        plot1<-lineplot_adult(data_pre_blood(),"metabolite",1:11)
        return(plot1)
      }

  })
  
  
  
  output$Plot_site<-renderPlot({
    if(is.null(data_lumen_microbial_all_BF())||is.null(data_lumen_metabolite_all_BF()))return(NULL)
    if(input$site_variable==""||input$plot_site=="")return(NULL)

    if(input$plot_site=="lumen" && input$site_variable=="microbial"){
      data_lumen_microbial_all_BF_input<-data_lumen_microbial_all_BF()
      data_lumen_microbial_all_SF_input<-data_lumen_microbial_all_SF()
      plot_BSL<-Tanks_cmp(data_lumen_microbial_all_BF_input,Bac8[1:5],"Baseline")
      plot_Intv<-Tanks_cmp(data_lumen_microbial_all_SF_input,colnames(data_lumen_microbial_all_SF_input),"Intervention")
      multi_plot <- ggarrange(plot_BSL+
                              # +rremove("ylab")+rremove("x.text")+rremove("x.ticks")+rremove("legend"), 
                              # +rremove("ylab")+
                                rremove("legend"), 
                              
                              plot_Intv,
                              # +rremove("ylab"), 
                              nrow=2, ncol=1,
                              # widths=c(5,3.1),
                              # heights=c(4,5.5),
                              # heights=c(10,18),
                              heights=c(6,10.5),
                              
                              align="v")
      return(multi_plot)
    }   ## supplementary fig.
    
    if(input$plot_site=="lumen" && input$site_variable=="metabolite"){
      if(is.null(data_lumen_metabolite_all_BF()))return(NULL)
      
      data_lumen_metabolite_all_BF_input<-data_lumen_metabolite_all_BF()
      data_lumen_metabolite_all_SF_input<-data_lumen_metabolite_all_SF()
      plot_BSL<-Tanks_cmp_mets(data_lumen_metabolite_all_BF_input,Bac8[1:5],"Baseline")
      plot_Intv<-Tanks_cmp_mets(data_lumen_metabolite_all_SF_input,colnames(data_lumen_metabolite_all_SF_input),"Intervention")
      
    multi_plot <- ggarrange(plot_BSL+
                            # +rremove("x.text")+rremove("x.ticks"), 
                            # +rremove("ylab")+
                              rremove("legend"), 
                            
                            plot_Intv
                            # +rremove("ylab")+rremove("x.text")+rremove("x.ticks"), 
                            # +rremove("ylab"),
                            ,
                            nrow=2, ncol=1,
                            # widths=c(5,3.1),
                            heights=c(6,7.5),
                            align="v")
    return(multi_plot)}
    
    if(input$plot_site=="mucosa" && input$site_variable=="microbial"){
      plot_BSL<-Tanks_cmp(data_mucosa_microbial_all_BF(),Bac8[1:5],"Baseline")
      plot_Intv<-Tanks_cmp(data_mucosa_microbial_all_SF(),colnames(data_mucosa_microbial_all_SF()),"Intervention")
      
      multi_plot <- ggarrange(plot_BSL
                              # +rremove("x.text")+rremove("x.ticks"), 
                              # +rremove("ylab")
                              +rremove("legend"), 
                              
                              plot_Intv
                              # +rremove("ylab"),
                              ,
                              # +rremove("x.text")+rremove("x.ticks"), 
                              nrow=2, ncol=1,
                              # widths=c(5,3.1),
                              # heights=c(4,5.5),
                              # heights=c(10,18),
                              heights=c(6,10.5),
                              
                              
                              align="v")
      return(multi_plot)
    }
    
    if(input$plot_site=="mucosa" && input$site_variable=="metabolite"){
      plot_BSL<-Tanks_cmp_mets(data_mucosa_metabolite_all_BF(),Bac8[1:5],"Baseline")
      plot_Intv<-Tanks_cmp_mets(data_mucosa_metabolite_all_SF(),colnames(data_mucosa_metabolite_all_SF()),"Intervention")
      
      multi_plot <- ggarrange(plot_BSL
                              # +rremove("x.text")+rremove("x.ticks"), 
                              # +rremove("ylab")
                              +rremove("legend"), 
                              
                              plot_Intv
                              # +rremove("ylab")+rremove("x.text")+rremove("x.ticks"), 
                              # +rremove("ylab"),
                              ,
                              nrow=2, ncol=1,
                              # widths=c(5,3.1),
                              heights=c(6,7.5),
                              align="v")
      return(multi_plot)
    }
    
    if(input$plot_site=="feces" && input$site_variable=="microbial"){
      if(is.null(data_feces_microbial_all()))return(NULL)
      
      plot2<-Feces_microbial(data_feces_microbial_all(),Time_vector(),0)
      return(plot2)
    }
    
    if(input$cohort_catergory=="Infant"&&input$plot_site=="feces" && input$site_variable=="metabolite"){
      if(is.null(data_feces_metabolite_all()))return(NULL)
      plot2<-Feces_metabolites_infant(data_feces_metabolite_all(),Time_vector(),1)
      return(plot2)
    }
    
    if(input$cohort_catergory=="Adult"&&input$plot_site=="feces" && input$site_variable=="metabolite"){
      if(is.null(data_feces_metabolite_all()))return(NULL)
     
        # plot2<-Feces_metabolites_adult(pre_Data,Time_vector(),1)
      
      plot2<-Feces_metabolites_adult(data_feces_metabolite_all(),Time_vector(),1)
      return(plot2)
    }
    
    if(input$plot_site=="blood"){
      if(is.null(data_blood_metabolite_all()))return(NULL)
      
      data_blood<-data_blood_metabolite_all()
      # data_blood$Regime<-Time_vector()
      plot3<-Plasma_metabolites_barplot(data_blood,Time_vector(),0)
      return(plot3)
    }
    
  })
  

    
  


  # 3d plot for infant, select several bacteria and metabolites, directly call matlab function to plot

  # index<-reactive({
  #   if(input$cohort_catergory=="Adult")return()
  #   return(switch(input$variable_3d,
  #          "Blg" = 4,
  #          "Bfr" = 3,
  #          "Acetate" = 13,
  #          "Butyrate" = 18,
  #          "Propionate" = 14)
  #   )
  # })
 
  # variable_3d_plot<-reactive({
  #   input$variable_3d
  # })
  # plot_3d<-reactive(do.call(run_3d_plot,list(input$cohort_catergory,index())))
  # navA <- reactive(do.call(simulate_nav, getParams("a")))
  
      
  # output$plot6<-renderImage({
  #   # type<-cohort_catergory()
  #   # type2<-diet_regime()
  #   if(is.null(input$cohort_catergory)||is.null(input$variable_3d)||is.null(diet_regime())) return()
  #   # plot_3d()
  #   run_matlab_code(run_3d_plot_code())
  #   
  #   outfile<-"Metabolites_3d.jpg"
  #   
  #   list(src = outfile,
  #        contentType = 'image/jpg',
  #        width = 400,
  #        height = 300,
  #        alt = "Spatiotemporal Variability")
  # })
  
  # Sys.time()
  
  outfile_3d<-reactive({
    req(input$variable_3d)
    my_timestamp<-gsub("[\\: \\-]","_",Sys.time())
    
    if(input$cohort_catergory==""||is.null(diet_regime()))return(NULL)
    if(input$cohort_catergory=="Infant"){
      return(paste0("Metabolite_3d_",my_timestamp,".jpg"))
    }
    if(input$cohort_catergory=="Adult"){
      return(paste0("Metabolite_3d_adult_",my_timestamp,".jpg"))
    }
  })
  
  # list.files(pattern = glob2rx(pattern = "Metabolite_*.jpg")) 
  
  
  # observe({
  #   req(input$cohort_catergory,diet_regime(), input$variable_3d,v$doCalc)
  #   index_input<-plot_3d_index()
  #   cohort_info<-input$cohort_catergory
  #   results2<-runMatlabFct("restul_total=run_Simulatin_3d(cohort_info,index_input)")
  # })
  
  
  observeEvent(input$variable_3d,{
    # req(diet_regime(),input$cohort_catergory)
    d3_plot$plot<-input$variable_3d
  })
  observeEvent(input$cohort_catergory,{
    d3_plot$plot<-FALSE
  })

  observe({
    # req(input$variable_3d,v$doCalc)
    req(d3_plot$plot)
    index_input<-plot_3d_index()
    cohort_info<-input$cohort_catergory
    outfile<-outfile_3d()
    
    met_files<-list.files(pattern = glob2rx("Metabolite_*.jpg"))
    unlink(met_files)
    
    met_files<-list.files(path="www/", pattern = glob2rx("Metabolite_*.jpg"))
    unlink(paste0("www/",met_files))
    
    results2<-runMatlabFct("restul_total=run_Simulatin_3d(cohort_info,index_input,outfile)")
  })
  
  
  # outfile_3d<-reactive({
  #   req(input$variable_3d)
  #   if(input$cohort_catergory==""||is.null(diet_regime()))return(NULL)
  #   if(input$cohort_catergory=="Infant"){
  #     return(paste0("Metabolite_3d_",plot_3d_index(),".jpg"))
  #   }
  #   if(input$cohort_catergory=="Adult"){
  #     return(paste0("Metabolite_3d_adult_",plot_3d_index(),".jpg"))
  #   }
  # })
  # 
  output$plot6<-renderUI({
    # req(v$doCalc)
    req(d3_plot$plot)
    
    if(input$variable_3d=="")return(NULL)
    # req(input$variable_3d)
    # outfile<-"Metabolites_3d.jpg"
    outfile<-outfile_3d()
    
    h6("Spatiotemporal Variability")
    br()
    tags$img(src = outfile,contentType = 'image/jpg',
             width = 300,
             height = 300
             )
    # list(src = outfile,
    #      contentType = 'image/jpg',
    #      width = 400,
    #      height = 300,
    #      alt = "Spatiotemporal Variability")
  })
  
  
  
})

