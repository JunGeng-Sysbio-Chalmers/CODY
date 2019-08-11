
lineplot_adult<-function(data,plot_variable,index){
  sample_index <- seq(1, nrow(data), by=10)
  data<-data[sample_index,index]

if(plot_variable=="microbial"){
	    colnames(data) <- c("t",  "Bfr", "Blg", "Bbv", "Bad", "Ehal", "Fpr", "Rint")
	    bac_data <- gather(data, -t, key="species", value="concentration")
  # colnames(bac_data) <- c("t", "species", "concentration")
  
  max_value <- ifelse(ceiling(max(bac_data$concentration)) %% 2 == 0, ceiling(max(bac_data$concentration)), ceiling(max(bac_data$concentration))+1)
  
  # bac_data$species <- as.factor(bac_data$species)
  
  # levels(bac_data$species)
  # levels(bac_data$species) <- c("B. breve", "B. fragilis", "B. longum", "B. thetaiotaomicron", "E. halli", "F. pransnitzii", "R. intestinalls")
  #levels(bac_data$species) <- c("Bad", "Bbv", "Bfr", "Blg", "Ehal", "Fpr", "Rint")
  bac_data$species <- factor(bac_data$species, levels = c("Bfr", "Blg", "Bbv", "Bad", "Ehal", "Fpr", "Rint") )
  
  
  # species <- as.character(unique(bac_data[, c("species")]))
  species <- unique(bac_data[, "species"])
  
  species_col_set <- data.frame(species=species)
  species_col_set$t <- species_col_set$concentration <- 1
  
  
  
  bac_plot <-  ggplot(bac_data, aes(x=t, y=concentration)) +
    geom_line(lwd=0.3) +
    facet_grid(species~.) +
    geom_rect(data = species_col_set, aes(fill = species),xmin = -Inf,xmax = Inf,
              ymin = -Inf,ymax = Inf,alpha = 0.2) +
    xlab("Time (h)") +
    ylab("Concentration (g/L)") +
    scale_x_continuous(expand=c(0,0.01),breaks = c(15,612,1200), labels=c("0" = "Initial", "612" = "Baseline", "1200" = "Intervention")) +
    scale_y_continuous(expand=c(0.05,0.1), limits=c(0,max_value), breaks=c(0, max_value/2, max_value), labels=c(0, max_value/2, max_value)) +
    #scale_y_continuous(expand=c(0.05,0.1)) +
    # ggtitle(figLabel) +
    theme_pander(base_family = "Helvetica", base_size = 10) +
    theme(legend.position = "none",
          plot.margin = unit( c(0.5,0.5,0.5,0.5) , units = "lines" ),
          panel.spacing = unit(0.8, "lines"),
          panel.border = element_rect(colour = "black", size=0.45),
          panel.grid = element_line(size = 0, linetype = 5, color="#EFEFEF"),
          axis.line = element_line(size = 0.3),
          strip.background = element_rect(fill="#F2F3F4"),
          strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")),
          # text=element_text(size=16,  family="Helvetica"),
          # text = element_text(size=14),
          # axis.text.x = element_text(),
          
          plot.title = element_text(color="grey5", size=14, face="plain"),
          # axis.text.x = element_text(size=14, family="Helvetica",angle=0, hjust=1),
          axis.text.x = element_text(size=9, family="Helvetica",angle=0,hjust=0.65),
          axis.title.x = element_text( size=14, family="Helvetica"),
          axis.title.y = element_text( size=14, family="Helvetica")
          
          
    )
  return(bac_plot)
}
else if(plot_variable=="metabolite"){
	colnames(data) <- c("t", "macs", "hexose",  "succinate", "acetate", "propionate", "lactate", "formate", "ethanol", "butyrate", "h2")
	mets_data <- gather(data, -t, key="metabolites", value="concentration")
  # colnames(mets_data) <- c("t", "metabolites", "concentration")
  
  max_value <- ifelse(ceiling(max(mets_data$concentration)) %% 2 == 0, ceiling(max(mets_data$concentration)), ceiling(max(mets_data$concentration))+1)
  
  mets_data$metabolites <- factor(mets_data$metabolites, levels=c("macs", "hexose",  "succinate", "acetate", "propionate", "lactate", "formate", "ethanol", "butyrate", "h2"))
  # levels(mets_data$metabolites)
  # "acetate"    "butyrate"   "ethanol"    "fiber"      "formate"    "h2"         "hexose"     "lactate"    "propionate" "succinate"
  #levels(mets_data$metabolites) <- c("Acetate", "Butyrate", "Ethanol", "Formate", "H2", "Hexose", "Lactate", "MACs", "Propionate", "Succinate")
  
  
  metabolites <- unique(mets_data[, "metabolites"])
  metabolites_col_set <- data.frame(metabolites=metabolites)
  metabolites_col_set$t <- metabolites_col_set$concentration <- 1
  
  
  mets_plot <-  ggplot(mets_data, aes(x=t, y=concentration)) +
    geom_line(lwd=0.3) +
    facet_wrap(~metabolites, ncol=2, scale="free_y") +
    geom_rect(data = metabolites_col_set, aes(fill = metabolites),xmin = -Inf,xmax = Inf,
              ymin = -Inf,ymax = Inf,alpha = 0.2) +
    xlab("Time (h)") +
    ylab("Concentration (mmol/L)") +
    scale_x_continuous(expand=c(0,0),breaks = c(5,612,1200), labels=c("0" = "Initial", "612" = "Baseline", "1200" = "Intervention")) +
    scale_y_continuous(expand=c(0.05,0.1)) +
    # ggtitle(figLabel) +
    theme_pander(base_family = "Helvetica", base_size=8) +
    theme(legend.position = "none",
          plot.margin = unit( c(0.5,1,0.5,0.5) , units = "lines" ),
          panel.spacing = unit(1, "lines"),
          panel.border = element_rect(colour = "black", size=0.45),
          panel.grid = element_line(size = 0, linetype = 5, color="#EFEFEF"),
          axis.line = element_line(size = 0.3),
          strip.background = element_rect(fill="#F2F3F4"),
          strip.text.x = element_text(margin = margin(0.05,0,0.05,0, "cm")),
          # text=element_text(size=14,  family="Helvetica"),
          plot.title = element_text(color="grey5", size=14, face="plain"),
          # axis.text.x = element_text(angle=45, hjust=1),
          # axis.text.x = element_text(size=13, family="Helvetica",angle=15, hjust=1),
          axis.text.x = element_text(size=9, family="Helvetica",hjust=0.65),
          axis.title.x = element_text( size=14, family="Helvetica"),
          axis.title.y = element_text( size=14, family="Helvetica")
    ) 
  return(mets_plot)
	}
# return(p)
}
#source("speciesPlot_new.R")
#source("metsPlot_new.R")

####################
# loading data
####################

# load("./Rdata/lumen_T1.Rdata")
# load("./Rdata/lumen_T2.Rdata")
# load("./Rdata/lumen_T3.Rdata")
# load("./Rdata/lumen_T4.Rdata")

# load("./Rdata/mucosa_T1.Rdata")
# load("./Rdata/mucosa_T2.Rdata")
# load("./Rdata/mucosa_T3.Rdata")
# load("./Rdata/mucosa_T4.Rdata")

# load("./Rdata/blood.Rdata")
# load("./Rdata/feces.Rdata")
# load("./Rdata/rectum.Rdata")

# #############
# # preprocessing data
# #############

# # only select 10% of samples from the dataset, and only keep 17 columns
# sample_index <- seq(1, nrow(lumen_T1), by=20)
# # 613

# ##############################
# # lumen
# #
# ##########
# # lumen Tank 1 
# lumen_T1_s <- lumen_T1[sample_index, 1:18]
# head(lumen_T1_s)
# lumen_T1_bac <- lumen_T1_s[,1:8]
# head(lumen_T1_bac)
# colnames(lumen_T1_bac) <- c("t",  "Bfr", "Blg", "Bbv", "Bad", "Ehal", "Fpr", "Rint")

# lumen_T1_mets <- lumen_T1_s[, c(1, 9:18)]
# head(lumen_T1_mets)
# colnames(lumen_T1_mets) <- c("t", "macs", "hexose",  "succinate", "acetate", "propionate", "lactate", "formate", "ethanol", "butyrate", "h2")

# speciesPlot(lumen_T1_bac, "./fig_20190529/lumen_T1_species", "Ascending Lumen")
# metsPlot(lumen_T1_mets, "./fig_20190529/lumen_T1_mets", "Ascending Lumen")

# ##########
# # lumen Tank 2
# lumen_T2_s <- lumen_T2[sample_index, 1:18]
# lumen_T2_bac <- lumen_T2_s[,1:8]
# colnames(lumen_T2_bac) <- c("t",  "Bfr", "Blg", "Bbv", "Bad", "Ehal", "Fpr", "Rint")
# lumen_T2_mets <- lumen_T2_s[, c(1, 9:18)]
# colnames(lumen_T2_mets) <- c("t",  "macs", "hexose",  "succinate", "acetate", "propionate", "lactate", "formate", "ethanol", "butyrate", "h2")

# speciesPlot(lumen_T2_bac, "./fig_20190529/lumen_T2_species", "Transverse Lumen")
# metsPlot(lumen_T2_mets, "./fig_20190529/lumen_T2_mets", "Transverse Lumen")

# ##########
# # lumen Tank 3
# lumen_T3_s <- lumen_T3[sample_index, 1:18]
# lumen_T3_bac <- lumen_T3_s[,1:8]
# colnames(lumen_T3_bac) <- c("t",  "Bfr", "Blg", "Bbv", "Bad", "Ehal", "Fpr", "Rint")
# lumen_T3_mets <- lumen_T3_s[, c(1, 9:18)]
# colnames(lumen_T3_mets) <- c("t",  "macs", "hexose",  "succinate", "acetate", "propionate", "lactate", "formate", "ethanol", "butyrate", "h2")

# speciesPlot(lumen_T3_bac, "./fig_20190529/lumen_T3_species", "Descending Lumen")
# metsPlot(lumen_T3_mets, "./fig_20190529/lumen_T3_mets", "Descending Lumen")

# ##########
# # lumen Tank 4
# lumen_T4_s <- lumen_T4[sample_index, 1:18]
# lumen_T4_bac <- lumen_T4_s[,1:8]
# colnames(lumen_T4_bac) <- c("t",  "Bfr", "Blg", "Bbv", "Bad", "Ehal", "Fpr", "Rint")
# lumen_T4_mets <- lumen_T4_s[, c(1, 9:18)]
# colnames(lumen_T4_mets) <- c("t",  "macs", "hexose",  "succinate", "acetate", "propionate", "lactate", "formate", "ethanol", "butyrate", "h2")

# speciesPlot(lumen_T4_bac, "./fig_20190529/lumen_T4_species", "Sigmoid Lumen")
# metsPlot(lumen_T4_mets, "./fig_20190529/lumen_T4_mets", "Sigmoid Lumen")


# ##############################
# # mucosa
# #
# ##########
# # mucosa Tank 1 
# mucosa_T1_s <- mucosa_T1[sample_index, 1:18]
# mucosa_T1_bac <- mucosa_T1_s[,1:8]
# colnames(mucosa_T1_bac) <- c("t",  "Bfr", "Blg", "Bbv", "Bad", "Ehal", "Fpr", "Rint")
# mucosa_T1_mets <- mucosa_T1_s[, c(1, 9:18)]
# colnames(mucosa_T1_mets) <- c("t",  "macs", "hexose",  "succinate", "acetate", "propionate", "lactate", "formate", "ethanol", "butyrate", "h2")

# speciesPlot(mucosa_T1_bac, "./fig_20190529/mucosa_T1_species", "Ascending Mucus")
# metsPlot(mucosa_T1_mets, "./fig_20190529/mucosa_T1_mets", "Ascending Mucus")

# ##########
# # mucosa Tank 2
# mucosa_T2_s <- mucosa_T2[sample_index, 1:18]
# mucosa_T2_bac <- mucosa_T2_s[,1:8]
# colnames(mucosa_T2_bac) <- c("t",  "Bfr", "Blg", "Bbv", "Bad", "Ehal", "Fpr", "Rint")
# mucosa_T2_mets <- mucosa_T2_s[, c(1, 9:18)]
# colnames(mucosa_T2_mets) <- c("t",  "macs", "hexose",  "succinate", "acetate", "propionate", "lactate", "formate", "ethanol", "butyrate", "h2")

# speciesPlot(mucosa_T2_bac, "./fig_20190529/mucosa_T2_species", "Transverse Mucus")
# metsPlot(mucosa_T2_mets, "./fig_20190529/mucosa_T2_mets", "Transverse Mucus")

# ##########
# # mucosa Tank 3
# mucosa_T3_s <- mucosa_T3[sample_index, 1:18]
# mucosa_T3_bac <- mucosa_T3_s[,1:8]
# colnames(mucosa_T3_bac) <- c("t",  "Bfr", "Blg", "Bbv", "Bad", "Ehal", "Fpr", "Rint")
# mucosa_T3_mets <- mucosa_T3_s[, c(1, 9:18)]
# colnames(mucosa_T3_mets) <- c("t",  "macs", "hexose",  "succinate", "acetate", "propionate", "lactate", "formate", "ethanol", "butyrate", "h2")

# speciesPlot(mucosa_T3_bac, "./fig_20190529/mucosa_T3_species", "Descending Mucus")
# metsPlot(mucosa_T3_mets, "./fig_20190529/mucosa_T3_mets", "Descending Mucus")

# ##########
# # mucosa Tank 4
# mucosa_T4_s <- mucosa_T4[sample_index, 1:18]
# mucosa_T4_bac <- mucosa_T4_s[,1:8]
# colnames(mucosa_T4_bac) <- c("t",  "Bfr", "Blg", "Bbv", "Bad", "Ehal", "Fpr", "Rint")
# mucosa_T4_mets <- mucosa_T4_s[, c(1, 9:18)]
# colnames(mucosa_T4_mets) <- c("t",  "macs", "hexose",  "succinate", "acetate", "propionate", "lactate", "formate", "ethanol", "butyrate", "h2")

# speciesPlot(mucosa_T4_bac, "./fig_20190529/mucosa_T4_species", "Sigmoid Mucus")
# metsPlot(mucosa_T4_mets, "./fig_20190529/mucosa_T4_mets", "Sigmoid Mucus")

# ##############################
# # blood
# #
# blood_s <- blood[sample_index, 1:11]
# blood_mets <- blood_s
# colnames(blood_mets) <- c("t",  "macs", "hexose",  "succinate", "acetate", "propionate", "lactate", "formate", "ethanol", "butyrate", "h2")

# metsPlot(blood_mets, "./fig_20190529/blood_mets", "Blood")

# ##############################
# # rectum
# #
# rectum_s <- rectum[sample_index, 1:18]
# rectum_bac <- rectum_s[,1:8]
# colnames(rectum_bac) <- c("t",  "Bfr", "Blg", "Bbv", "Bad", "Ehal", "Fpr", "Rint")
# rectum_mets <- rectum_s[, c(1, 9:18)]
# colnames(rectum_mets) <- c("t",  "macs", "hexose",  "succinate", "acetate", "propionate", "lactate", "formate", "ethanol", "butyrate", "h2")

# speciesPlot(rectum_bac, "./fig_20190529/rectum_species", "Rectum")
# metsPlot(rectum_mets, "./fig_20190529/rectum_mets", "Rectum")


# ##############################
# # feces
# #
# feces_s <- feces[sample_index, 1:18]
# feces_bac <- feces_s[,1:8]
# colnames(feces_bac) <- c("t",  "Bfr", "Blg", "Bbv", "Bad", "Ehal", "Fpr", "Rint")
# feces_mets <- feces_s[, c(1, 9:18)]
# colnames(feces_mets) <- c("t",  "macs", "hexose",  "succinate", "acetate", "propionate", "lactate", "formate", "ethanol", "butyrate", "h2")

# speciesPlot(feces_bac, "./fig_20190529/feces_species", "Feces")
# metsPlot(feces_mets, "./fig_20190529/feces_mets", "Feces")

