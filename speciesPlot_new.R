speciesPlot <- function(data, outLabel, figLabel) {
  
  # set output file
  outpdf <- paste(outLabel, ".pdf", sep="")
  outjpg <- paste(outLabel, ".jpg", sep="")
  outtif <- paste(outLabel, ".tif", sep="")
  outpng <- paste(outLabel, ".png", sep="")
  
  bac_data <- gather(data, -t, key="species", value="concentration")
  # colnames(bac_data) <- c("t", "species", "concentration")
  
  max_value <- ifelse(ceiling(max(bac_data$concentration)) %% 2 == 0, ceiling(max(bac_data$concentration)), ceiling(max(bac_data$concentration))+1)
  
  bac_data$species <- factor(bac_data$species, levels = c("Bfr", "Blg", "Bbv", "Bad", "Ehal", "Fpr", "Rint"), order=TRUE )
  # levels(bac_data$species)
  # levels(bac_data$species) <- c("B. breve", "B. fragilis", "B. longum", "B. thetaiotaomicron", "E. halli", "F. pransnitzii", "R. intestinalls")
  #levels(bac_data$species) <- c("Bad", "Bbv", "Bfr", "Blg", "Ehal", "Fpr", "Rint")
  
  
  species <- as.character(unique(bac_data[, c("species")]))
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
    ggtitle(figLabel) +
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
  bac_plot
  ggsave(outpdf, bac_plot, height=12, width=10, units = "cm")
  # ggsave(outjpg, bac_plot, height=12, width=10, units = "cm", dpi=300)
  # ggsave(outtif, bac_plot, height=12, width=10, units = "cm", dpi=300, device="tiff")
  ggsave(outpng, bac_plot, height=12, width=10, units = "cm", dpi=300, device="png")
}