speciesPlot_12M_Jun <- function(data, outLabel, figLabel) {
  
  # set output file
  outpdf <- paste(outLabel, ".pdf", sep="")
  outjpg <- paste(outLabel, ".jpg", sep="")
  outtif <- paste(outLabel, ".tif", sep="")
  outpng <- paste(outLabel, ".png", sep="")
  
  bac_data <- gather(data, -t, key="species", value="concentration")
  # colnames(bac_data) <- c("t", "species", "concentration")
  
  max_value <- ifelse(ceiling(max(bac_data$concentration)) %% 2 == 0, ceiling(max(bac_data$concentration)), ceiling(max(bac_data$concentration))+1)
  
  bac_data$species <- as.factor(bac_data$species)
  # levels(bac_data$species)
  # levels(bac_data$species) <- c("B. breve", "B. fragilis", "B. longum", "B. thetaiotaomicron", "E. halli", "F. pransnitzii", "R. intestinalls")
  levels(bac_data$species) <- c("Bad", "Bbv", "Bfr", "Blg", "Bth", "Ehal", "Fpr", "Rint")
  
  
  species <- as.character(unique(bac_data[, c("species")]))
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
    ggtitle(figLabel) +
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
  ggsave(outpdf, bac_plot, height=12, width=10, units = "cm")
  # ggsave(outjpg, bac_plot, height=12, width=10, units = "cm", dpi=300)
  # ggsave(outtif, bac_plot, height=12, width=10, units = "cm", dpi=300, device="tiff")
  # ggsave(outpng, bac_plot, height=12, width=10, units = "cm", dpi=300, device="png")
}
