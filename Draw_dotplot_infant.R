Draw_dotplot_infant <- function(pre_data) {

    color8<-c(brewer.pal(n=8,"Dark2")[1:5],c("#BC80BD", "#E6AB02", "#80B1D3"))
    
    data<-pre_data

experiment_data<-data[c(1,4,7),]
experiment_data[1:3,2:ncol(experiment_data)]<-experiment_data[1:3,2:ncol(experiment_data)]/100
# rownames(experiment_data)[2:4]<-c('Newborn', '4 Month', '12 Month')
experiment_data$Time<-c('Newborn', '4 Months', '12 Months')

prediction_data<-data[c(1,3,6),]
prediction_data[1:3,2:ncol(experiment_data)]<-prediction_data[1:3,2:ncol(experiment_data)]/100
# rownames(prediction_data)[2:4]<-c('Newborn', '4 Month', '12 Month')
prediction_data$Time<-c('Newborn', '4 Months', '12 Months')

# prediction_data<-gather(prediction_data, Type, RA, Newborn:'12 Month')
prediction_data<-melt(prediction_data, id='Time')
colnames(prediction_data)[1:3]<-c('Type','Species','prediction')

experiment_data<-melt(experiment_data, id='Time')
colnames(experiment_data)[1:3]<-c('Type','Species','experiment')


all_data <- merge(experiment_data, prediction_data, by = c("Type", "Species"))



rsqred<-function (x,y) {
  y_hat <- y
  y <- x
  y_bar = mean(y)
  RSS = sum((y-y_hat)^2)
  TSS = sum((y-y_bar)^2)
  R2 = 1 - RSS/TSS
  return(R2)
}




# change order of species
species_name <- factor(all_data$Species, labels = c("Bfr", "Bth", "Bad", "Bbv", "Blg", "Ehal", "Fpr", "Rint") )
species_name_order <- factor(species_name, levels = c("Bfr", "Bth", "Blg", "Bbv", "Bad", "Ehal", "Fpr", "Rint"), ordered = TRUE)

all_data$Species <- species_name_order


colnames(all_data)[3:4]<-c("RelativeAbundance","prediction")
M4_data <- all_data %>% filter(Type == "4 Months")
M12_data <- all_data %>% filter(Type == "12 Months")

# M4_data$RelativeAbundance<-M4_data$RelativeAbundance/100
# baseline_data$Prediction<-baseline_data$Prediction/100
# 
# intervention_data <- all_data %>% filter(type == "Average Intervention")
# intervention_data$Metagenomics<-intervention_data$Metagenomics/100
# intervention_data$Prediction<-intervention_data$Prediction/100

M4_lm_yx <- lm(RelativeAbundance ~ 0 + offset(1*prediction), data=M4_data)
M4_lm_yx_summary <-  summary(M4_lm_yx)
M4_r2<-rsqred(M4_data$RelativeAbundance,M4_data$prediction)   ## 


M4_lm <- lm(RelativeAbundance ~ prediction, data=M4_data)
M4_lm_summary <-  summary(M4_lm)
M4_lm_r2 <- M4_lm_summary$adj.r.squared

M4_cor <- cor.test(M4_data$RelativeAbundance, M4_data$prediction)
M4_cor
# Pearson's product-moment correlation
# 
# data:  M4_data$RelativeAbundance and M4_data$prediction
# t = 38.897, df = 6, p-value = 1.929e-08
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.9886420 0.9996572
# sample estimates:
#      cor 
# 0.998023 
  
M4_RSS <- c(crossprod(M4_lm$residuals))
M4_MSE <- M4_RSS / length(M4_lm$residuals)
M4_RMSE <- sqrt(M4_MSE)
M4_RMSE
# 0.003789035

## Add an alpha value to a colour
add.alpha <- function(col, alpha=1) {
  if (missing(col)) {
    stop("Please provide a vector of colours.")
  }
  apply(
    sapply(col, col2rgb) / 255, 2,
    function(x)
      rgb(x[1], x[2], x[3], alpha = alpha)
  )
}

scaleFUN <- function(x) sprintf("%.2f", x)

col_used<-color8

M4_plot <- ggplot(M4_data, aes(x=prediction, y=RelativeAbundance)) +
  geom_abline(intercept = 0, color="black", size=0.5) +
  geom_point(aes(color=Species), size=2.2) +
  scale_color_manual(values=col_used) +
  # geom_text(aes(x=0.021, y=0.097, label = paste("R^2 == ", round(M4_lm_r2, 2))), parse=T, color="black",size=3,family="Helvetica") +
  # geom_text(aes(x=0.021, y=0.097, label = paste("R^2 == ", round(M4_r2, 2))), parse=T, color="black",size=3,family="Helvetica") +
  
    #coord_fixed() +
  # xlim(0, 10) +
  ylab("Exp. Relative Abundance")+
  xlab("Pre. Relative Abundance")+
  
  scale_x_continuous(limits=c(0,0.1),breaks=c(0,0.05,0.1),labels=scaleFUN)+
  # ylim(0, 10) +
  scale_y_continuous(limits=c(0,0.1),breaks=c(0,0.05,0.1),labels=scaleFUN)+
  #geom_smooth(method='lm') +
  theme(
    legend.title = element_blank(), legend.position = "none",
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    text = element_text(size = 7, family = "Helvetica"),
    axis.line = element_line(color = "black", size = 0.12),
    axis.text = element_text(size = 8, colour="black", family = "Helvetica"),
    axis.title = element_text(size = 8,  family = "Helvetica", colour = "black"),
    axis.text.x = element_text(size = 8, colour="black",family = "Helvetica", vjust=1),
    axis.ticks.length=unit(0.18, "line"),
    axis.ticks = element_line(colour = "black", size = 0.2),
    axis.title.y = element_text(size = 10, colour="black",family = "Helvetica", vjust=0.018),
    axis.title.x = element_text(size = 10, colour="black",family = "Helvetica"),
    # axis.ticks = element_blank()
    plot.margin = unit(c(0, 1.5, 0, 0), "pt")
  )

M4_plot


M12_lm_yx <- lm(RelativeAbundance ~ 0 + offset(1*prediction), data=M12_data)
M12_lm_yx_summary <-  summary(M12_lm_yx)
M12_r2<-rsqred(M12_data$RelativeAbundance,M12_data$prediction)   ## 
# M12_r2<-0.85

M12_lm <- lm(RelativeAbundance ~ prediction, data=M12_data)
M12_lm_summary <-  summary(M12_lm)
M12_lm_r2 <- M12_lm_summary$adj.r.squared

M12_cor <- cor.test(M12_data$RelativeAbundance, M12_data$prediction)
M12_cor
# Pearson's product-moment correlation
# 
# data:  M12_data$RelativeAbundance and M12_data$prediction
# t = 7.6891, df = 6, p-value = 0.0002534
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.7552245 0.9916637
# sample estimates:
#       cor 
# 0.9528198 

M12_RSS <- c(crossprod(M12_lm$residuals))
M12_MSE <- M12_RSS / length(M12_lm$residuals)
M12_RMSE <- sqrt(M12_MSE)
M12_RMSE
# 0.002453427

M12_plot <- ggplot(M12_data, aes(x=prediction, y=RelativeAbundance)) +
  geom_abline(intercept = 0, color="black", size=0.5) +
  geom_point(aes(color=Species), size=2.2) +
  scale_color_manual(values=col_used) +
  # geom_text(aes(x=0.021, y=0.097, label = paste("R^2 == ", round(M12_lm_r2, 2))), parse=T, colour="black",size=3,family="Helvetica") +
  # geom_text(aes(x=0.021, y=0.097, label = paste("R^2 == ", round(M12_r2, 2))), parse=T, colour="black",size=3,family="Helvetica") +
  
    #coord_fixed() +
  ylab("Exp. Relative Abundance")+
  xlab("Pre. Relative Abundance")+
  
  scale_x_continuous(limits=c(0,0.1),breaks=c(0,0.05,0.1),labels=scaleFUN)+
  # ylim(0, 10) +
  scale_y_continuous(limits=c(0,0.1),breaks=c(0,0.05,0.1),labels=scaleFUN)+
  #geom_smooth(method='lm') +
  theme(
    legend.title = element_blank(), legend.position = "none",
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
    text = element_text(size = 7, colour="black",family = "Helvetica"),
    axis.line = element_line(colour = "black", size = 0.12),
    axis.text = element_text(size = 8, colour = "black",family = "Helvetica"),
    axis.title = element_text(size = 8,  family = "Helvetica", colour = "black"),
    axis.text.x = element_text(size = 8, colour="black", family = "Helvetica",vjust=1),
    axis.title.y = element_text(size = 10, colour="black",family = "Helvetica", vjust=0.028),
    axis.title.x = element_text(size = 10, colour="black",family = "Helvetica"),
    
    axis.ticks.length=unit(0.18, "line"),
    axis.ticks = element_line(colour = "black", size = 0.2),
    plot.margin = unit(c(0, 0, 0, 3), "pt")
  )

M12_plot

all_plot <- ggarrange(M4_plot, M12_plot, ncol=2, align="v")
return(all_plot)

# ggsave("./fig/M4_M12_dotplot_formula_mixed_Jun_0517-1.pdf", all_plot, height = 4.6, width = 10.3, units = "cm")
# # ggsave("./fig/M4_M12_dotplot_formula_mixed_Jun_0517.jpg", all_plot, height = 4.6, width = 10.3, units = "cm", dpi = 300)
# # ggsave("./fig/M4_M12_dotplot_formula_mixed_Jun_0517.tif", all_plot, height = 4.6, width = 11.3, units = "cm", dpi = 300, device = "tiff")
# ggsave("./fig/M4_M12_dotplot_formula_mixed_Jun_0517-1.png", all_plot, height = 4.6, width = 10.3, units = "cm", dpi = 300, device = "png")
}

