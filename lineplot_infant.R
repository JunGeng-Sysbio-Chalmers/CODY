lineplot_infant<-function(data,plot_variable,index){
  # data<-data_1
  # index<-1:9
  # # index<-c(1, 10:19)
  # plot_variable<-"microbial"
  
  sample_index <- seq(1, nrow(data), by=10)
    # data<-data[sample_index,]
    data<-data[sample_index,index]
    
  
  if(plot_variable=="microbial"){
    # colnames(data) <- str_replace(colnames(data),"BIOMASS_","")
    # colnames(data)[1] <- "t"
    # colnames(data) <- c("t", "Btheta", "Bfragilis", "Blongum", "Bbreve", "Bado", "Ehalli", "Fpransnitzii", "Rintestinalls")
    colnames(data) <- c("t", "Bth", "Bfr", "Blg", "Bbv", "Bad", "Ehal", "Fpr", "Rint")
    
    bac_data <- gather(data, -t, key="species", value="concentration")
    # colnames(bac_data) <- c("t", "species", "concentration")
    
    max_value <- ifelse(ceiling(max(bac_data$concentration)) %% 2 == 0, ceiling(max(bac_data$concentration)), ceiling(max(bac_data$concentration))+1)
    
    bac_data$species <- as.factor(bac_data$species)
    # bac_data$species <- factor(bac_data$species,levels = c("Bth", "Bfr", "Blg",  "Bbv", "Bad", "Ehal", "Fpr", "Rint"),ordered = TRUE)
    
        # levels(bac_data$species)
    # levels(bac_data$species) <- c("B. breve", "B. fragilis", "B. longum", "B. thetaiotaomicron", "E. halli", "F. pransnitzii", "R. intestinalls")
    # levels(bac_data$species) <- c("Bad", "Bbv", "Bfr", "Blg", "Bth", "Ehal", "Fpr", "Rint")
    
    # species <- as.character(unique(bac_data[, c("species")]))
    species <- unique(bac_data[, "species"])
    
    species_col_set <- data.frame(species=species)
    species_col_set$t <- species_col_set$concentration <- 1
    
    
    bac_plot <-  ggplot(bac_data, aes(x=t, y=concentration)) +
      geom_line(lwd=0.5) +
      facet_grid(species~.) +
      geom_rect(data = species_col_set, aes(fill = species),xmin = -Inf,xmax = Inf,
                ymin = -Inf,ymax = Inf,alpha = 0.2) +
      xlab("Time (h)") +
      ylab("Concentration (g/L)") +
      # scale_x_continuous(expand=c(0,0.01)) +
      scale_x_continuous(expand=c(0,0), limits=c(0, 1033.52), breaks=c(0, 300,  1033.52 ), labels=c('NB','4M','12M')) +
      scale_y_continuous(expand=c(0.05,0.1), limits=c(0,max_value), breaks=c(0, max_value/2, max_value), labels=c(0, max_value/2, max_value)) +
      #scale_y_continuous(expand=c(0.05,0.1)) +
      # ggtitle(figLabel) +
      theme_pander(base_family = "Helvetica", base_size = 10) +
      theme(legend.position = "none",
            plot.title = element_text(color="grey5", size=13, face="plain", family="Helvetica"),
            axis.title.x = element_text(color="black", size=13, face="plain",family="Helvetica"),
            axis.title.y = element_text(color="black", size=13, face="plain",family="Helvetica"),
            plot.margin = unit( c(0.5,0.5,0.5,0.5) , units = "lines" ),
            panel.spacing = unit(0.5, "lines"),
            panel.border = element_rect(colour = "black", size=0.5),
            panel.grid = element_line(size = 0, linetype = 5, color="#EFEFEF"),
            axis.line = element_line(size = 0.3),
            strip.background = element_rect(fill="#F2F3F4"),
            axis.text.x = element_text(size=9, family="Helvetica",angle=0,hjust=0.65),
            strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm"), family="Helvetica", size=6)
            
      )
    bac_plot
    # ggsave('1.pdf', bac_plot, height=12, width=10, units = "cm")
    
    return(bac_plot)
  }
  else {
    colnames(data) <- c("t", "macs", "hexose", "succinate", "acetate", "propionate", "lactate", "formate", "ethanol", "butyrate", "h2")
    mets_data <- gather(data, -t, key="metabolites", value="concentration")
    # colnames(mets_data) <- c("t", "metabolites", "concentration")
    
    max_value <- ifelse(ceiling(max(mets_data$concentration)) %% 2 == 0, ceiling(max(mets_data$concentration)), ceiling(max(mets_data$concentration))+1)
    
    mets_data$metabolites <- as.factor(mets_data$metabolites)
    # levels(mets_data$metabolites)
    # "acetate"    "butyrate"   "ethanol"    "fiber"      "formate"    "h2"         "hexose"     "lactate"    "propionate" "succinate"
    # levels(mets_data$metabolites) <- c("Acetate", "Butyrate", "Ethanol", "Formate", "H2", "Hexose", "Lactate", "MACs", "Propionate", "Succinate")
    
    
    # metabolites <- as.character(unique(mets_data[, c("metabolites")]))
    metabolites <- unique(mets_data[, c("metabolites")])
    
    metabolites_col_set <- data.frame(metabolites=metabolites)
    metabolites_col_set$t <- metabolites_col_set$concentration <- 1
    
    
    mets_plot <-  ggplot(mets_data, aes(x=t, y=concentration)) +
      geom_line(lwd=0.25) +
      facet_wrap(~metabolites, ncol=2, scale="free_y") +
      geom_rect(data = metabolites_col_set, aes(fill = metabolites),xmin = -Inf,xmax = Inf,
                ymin = -Inf,ymax = Inf,alpha = 0.2) +
      xlab("Time (h)") +
      ylab("Concentration (mM)") +
      scale_x_continuous(expand=c(0,0), limits=c(0, 1033.52), breaks=c(0, 302,  1033.52 ), labels=c('NB','M4','M12')) +
      scale_y_continuous(expand=c(0.05,0.1)) +
      # ggtitle(figLabel) +
      theme_pander(base_family = "Helvetica", base_size=8) +
      theme(legend.position = "none",
            plot.title = element_text(color="grey5", size=13, face="plain", family="Helvetica"),
            axis.title.x = element_text(color="black", size=13, face="plain",family="Helvetica"),
            axis.title.y = element_text(color="black", size=13, face="plain",family="Helvetica"),
            plot.margin = unit( c(0.5,1,0.5,0.5) , units = "lines" ),
            panel.spacing = unit(0.5, "lines"),
            panel.border = element_rect(colour = "black", size=0.5),
            panel.grid = element_line(size = 0, linetype = 5, color="#EFEFEF"),
            axis.line = element_line(size = 0.3),
            axis.text.x = element_text(size=9, family="Helvetica",hjust=0.65),
            strip.background = element_rect(fill="#F2F3F4"),
            strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm"), family="Helvetica", size=7)
      ) 
    # ggsave('2.pdf', bac_plot, height=12, width=10, units = "cm")
    
    return(mets_plot)
    
  }
  
# return(p)

  
  
}
# ####################
# # loading functions for plotting
# ####################
# source("speciesPlot_12M_Jun.R")
# source("metsPlot_12M_Jun.R")
# 
# ####################
# # loading data
# ####################
# 
# load("./Rdata_12M/12M_lumen_T1.Rdata")
# load("./Rdata_12M/12M_lumen_T2.Rdata")
# load("./Rdata_12M/12M_lumen_T3.Rdata")
# load("./Rdata_12M/12M_lumen_T4.Rdata")
# 
# load("./Rdata_12M/12M_mucosa_T1.Rdata")
# load("./Rdata_12M/12M_mucosa_T2.Rdata")
# load("./Rdata_12M/12M_mucosa_T3.Rdata")
# load("./Rdata_12M/12M_mucosa_T4.Rdata")
# 
# load("./Rdata_12M/12M_blood.Rdata")
# load("./Rdata_12M/12M_feces.Rdata")
# load("./Rdata_12M/12M_rectum.Rdata")
# 
# #############
# # preprocessing data
# #############
# 
# # only select 10% of samples from the dataset, and only keep 17 columns
# 
# ##############################
# # lumen
# #
# ##########
# # lumen Tank 1 
# lumen_T1_s <- lumen_T1[sample_index, 1:19]
# lumen_T1_bac <- lumen_T1_s[,1:9]
# colnames(lumen_T1_bac) <- c("t", "Btheta", "Bfragilis", "Blongum", "Bbreve", "Bado", "Ehalli", "Fpransnitzii", "Rintestinalls")
# lumen_T1_mets <- lumen_T1_s[, c(1, 10:19)]
# colnames(lumen_T1_mets) <- c("t", "macs", "hexose", "succinate", "acetate", "propionate", "lactate", "formate", "ethanol", "butyrate", "h2")
# 
# speciesPlot_12M_Jun(lumen_T1_bac, "./fig_12M_Jun/12M_lumen_T1_species", "Ascending Lumen")
# metsPlot_12M_Jun(lumen_T1_mets, "./fig_12M_Jun/12M_lumen_T1_mets", "Ascending Lumen")
# 
# ##########
# # lumen Tank 2
# 
# 
# ##########
# # lumen Tank 3
# lumen_T3_s <- lumen_T3[sample_index, 1:19]
# lumen_T3_bac <- lumen_T3_s[,1:9]
# colnames(lumen_T3_bac) <- c("t", "Btheta", "Bfragilis", "Blongum", "Bbreve", "Bado", "Ehalli", "Fpransnitzii", "Rintestinalls")
# lumen_T3_mets <- lumen_T3_s[, c(1, 10:19)]
# colnames(lumen_T3_mets) <- c("t", "macs", "hexose", "succinate", "acetate", "propionate", "lactate", "formate", "ethanol", "butyrate", "h2")
# 
# speciesPlot_12M_Jun(lumen_T3_bac, "./fig_12M_Jun/12M_lumen_T3_species", "Descending Lumen")
# metsPlot_12M_Jun(lumen_T3_mets, "./fig_12M_Jun/12M_lumen_T3_mets", "Descending Lumen")
# 
# ##########
# # lumen Tank 4
# lumen_T4_s <- lumen_T4[sample_index, 1:19]
# lumen_T4_bac <- lumen_T4_s[,1:9]
# colnames(lumen_T4_bac) <- c("t", "Btheta", "Bfragilis", "Blongum", "Bbreve", "Bado", "Ehalli", "Fpransnitzii", "Rintestinalls")
# lumen_T4_mets <- lumen_T4_s[, c(1, 10:19)]
# colnames(lumen_T4_mets) <- c("t", "macs", "hexose", "succinate", "acetate", "propionate", "lactate", "formate", "ethanol", "butyrate", "h2")
# 
# speciesPlot_12M_Jun(lumen_T4_bac, "./fig_12M_Jun/12M_lumen_T4_species", "Sigmoid Lumen")
# metsPlot_12M_Jun(lumen_T4_mets, "./fig_12M_Jun/12M_lumen_T4_mets", "Sigmoid Lumen")
# 
# 
# ##############################
# # mucosa
# #
# ##########
# # mucosa Tank 1 
# mucosa_T1_s <- mucosa_T1[sample_index, 1:19]
# mucosa_T1_bac <- mucosa_T1_s[,1:9]
# colnames(mucosa_T1_bac) <- c("t", "Btheta", "Bfragilis", "Blongum", "Bbreve", "Bado", "Ehalli", "Fpransnitzii", "Rintestinalls")
# mucosa_T1_mets <- mucosa_T1_s[, c(1, 10:19)]
# colnames(mucosa_T1_mets) <- c("t", "macs", "hexose", "succinate", "acetate", "propionate", "lactate", "formate", "ethanol", "butyrate", "h2")
# 
# speciesPlot_12M_Jun(mucosa_T1_bac, "./fig_12M_Jun/12M_mucosa_T1_species", "Ascending Mucus")
# metsPlot_12M_Jun(mucosa_T1_mets, "./fig_12M_Jun/12M_mucosa_T1_mets", "Ascending Mucus")
# 
# ##########
# # mucosa Tank 2
# mucosa_T2_s <- mucosa_T2[sample_index, 1:19]
# mucosa_T2_bac <- mucosa_T2_s[,1:9]
# colnames(mucosa_T2_bac) <- c("t", "Btheta", "Bfragilis", "Blongum", "Bbreve", "Bado", "Ehalli", "Fpransnitzii", "Rintestinalls")
# mucosa_T2_mets <- mucosa_T2_s[, c(1, 10:19)]
# colnames(mucosa_T2_mets) <- c("t", "macs", "hexose", "succinate", "acetate", "propionate", "lactate", "formate", "ethanol", "butyrate", "h2")
# 
# speciesPlot_12M_Jun(mucosa_T2_bac, "./fig_12M_Jun/12M_mucosa_T2_species", "Transverse Mucus")
# metsPlot_12M_Jun(mucosa_T2_mets, "./fig_12M_Jun/12M_mucosa_T2_mets", "Transverse Mucus")
# 
# ##########
# # mucosa Tank 3
# mucosa_T3_s <- mucosa_T3[sample_index, 1:19]
# mucosa_T3_bac <- mucosa_T3_s[,1:9]
# colnames(mucosa_T3_bac) <- c("t", "Btheta", "Bfragilis", "Blongum", "Bbreve", "Bado", "Ehalli", "Fpransnitzii", "Rintestinalls")
# mucosa_T3_mets <- mucosa_T3_s[, c(1, 10:19)]
# colnames(mucosa_T3_mets) <- c("t", "macs", "hexose", "succinate", "acetate", "propionate", "lactate", "formate", "ethanol", "butyrate", "h2")
# 
# speciesPlot_12M_Jun(mucosa_T3_bac, "./fig_12M_Jun/12M_mucosa_T3_species", "Descending Mucus")
# metsPlot_12M_Jun(mucosa_T3_mets, "./fig_12M_Jun/12M_mucosa_T3_mets", "Descending Mucus")
# 
# ##########
# # mucosa Tank 4
# mucosa_T4_s <- mucosa_T4[sample_index, 1:19]
# mucosa_T4_bac <- mucosa_T4_s[,1:9]
# colnames(mucosa_T4_bac) <- c("t", "Btheta", "Bfragilis", "Blongum", "Bbreve", "Bado", "Ehalli", "Fpransnitzii", "Rintestinalls")
# mucosa_T4_mets <- mucosa_T4_s[, c(1, 10:19)]
# colnames(mucosa_T4_mets) <- c("t", "macs", "hexose", "succinate", "acetate", "propionate", "lactate", "formate", "ethanol", "butyrate", "h2")
# 
# speciesPlot_12M_Jun(mucosa_T4_bac, "./fig_12M_Jun/12M_mucosa_T4_species", "Sigmoid Mucus")
# metsPlot_12M_Jun(mucosa_T4_mets, "./fig_12M_Jun/12M_mucosa_T4_mets", "Sigmoid Mucus")
# 
# ##############################
# # blood
# #
# blood_s <- blood[sample_index, 1:11]
# blood_mets <- blood_s
# colnames(blood_mets) <- c("t", "macs", "hexose", "succinate", "acetate", "propionate", "lactate", "formate", "ethanol", "butyrate", "h2")
# 
# metsPlot_12M_Jun(blood_mets, "./fig_12M_Jun/12M_blood_mets", "Blood")
# 
# ##############################
# # rectum
# #
# rectum_s <- rectum[sample_index, 1:19]
# rectum_bac <- rectum_s[,1:9]
# colnames(rectum_bac) <- c("t", "Btheta", "Bfragilis", "Blongum", "Bbreve", "Bado", "Ehalli", "Fpransnitzii", "Rintestinalls")
# rectum_mets <- rectum_s[, c(1, 10:19)]
# colnames(rectum_mets) <- c("t", "macs", "hexose", "succinate", "acetate", "propionate", "lactate", "formate", "ethanol", "butyrate", "h2")
# 
# speciesPlot_12M_Jun(rectum_bac, "./fig_12M_Jun/12M_rectum_species", "Rectum")
# metsPlot_12M_Jun(rectum_mets, "./fig_12M_Jun/12M_rectum_mets", "Rectum")
# 
# 
# ##############################
# # feces
# #
# feces_s <- feces[sample_index, 1:19]
# feces_bac <- feces_s[,1:9]
# colnames(feces_bac) <- c("t", "Btheta", "Bfragilis", "Blongum", "Bbreve", "Bado", "Ehalli", "Fpransnitzii", "Rintestinalls")
# feces_mets <- feces_s[, c(1, 10:19)]
# colnames(feces_mets) <- c("t", "macs", "hexose", "succinate", "acetate", "propionate", "lactate", "formate", "ethanol", "butyrate", "h2")
# 
# speciesPlot_12M_Jun(feces_bac, "./fig_12M_Jun/12M_feces_species", "Feces")
# metsPlot_12M_Jun(feces_mets, "./fig_12M_Jun/12M_feces_mets", "Feces")

