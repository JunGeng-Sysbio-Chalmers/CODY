Draw_Fig3d <- function(experiment_file,prd_data) {

rsqred<-function (x,y) {
  y_hat <- y
  y <- x
  y_bar = mean(y)
  RSS = sum((y-y_hat)^2)
  TSS = sum((y-y_bar)^2)
  R2 = 1 - RSS/TSS
  return(R2)
}
##
# exp_file <- "./data/experiment.txt"
# prd_file <- "./data/prediction.txt"
experiment_file <- "Bac7_sample_ra.txt"
all_data <- read_tsv(experiment_file)

# only keep data at timepoint 0 & 14
all_data <- all_data %>% filter(Time == 0 | Time == 14)

# select the main 4 strain by strain names
main_data <- all_data %>% filter(strain %in% c("Bfr", "Blg", "Bbv", "Bad", "Ehal", "Fpr", "Rint"))

# other3_data <- all_data %>% filter(! strain %in% c("Bfr", "Blg", "Bbv"))

all_summary <- aggregate(RA~ Time+strain, mean, data=main_data)


color_7spc<-c(brewer.pal(n=8,"Dark2")[2:5],"#BC80BD", "#E6AB02", "#80B1D3")

all_summary$strain <- factor(all_summary$strain, levels = c("Bfr", "Blg", "Bbv", "Bad", "Ehal", "Fpr", "Rint"), order=TRUE)
exp_data<-all_summary
colnames(exp_data)<-c("type","Species", "Metagenomics")
# exp_data <- read.table(file = experiment_file, sep = "\t", header = T, stringsAsFactors = F)
# exp_data$type <- rownames(exp_data)
# prd_data <- read.table(file=prediction_file, sep="\t", header=T, stringsAsFactors=F)

# prd_data <- read.table(file=prediction_file, sep="\t", header=T, stringsAsFactors=F)

# prd_data <- read.table(file = prd_file, sep = "\t", header = T, stringsAsFactors = F)
# prd_data$type <- rownames(prd_data)
# prd_data$type <- prd_data$Time

# head(exp_data)
# 
# exp_data_list <- exp_data %>% gather("Species", "Metagenomics", -type)
# head(exp_data_list)
colnames(prd_data)<-c("Time","Bfr","Blg","Bbv","Bad","Ehal","Fpr","Rint")
prd_data$Time<-c("D0","D14")
prd_data_list <- prd_data %>% gather("Species", "Prediction",-Time)
prd_data_list$Prediction<-prd_data_list$Prediction
colnames(prd_data_list)[1]<-"type"

all_data <- full_join(exp_data, prd_data_list, by=c("type", "Species"))
head(all_data)

# change order of species
species_name <- factor(all_data$Species, labels = c("Bad", "Bbv", "Bfr", "Blg", "Ehal", "Fpr", "Rint") )
species_name_order <- factor(species_name, levels = c("Bfr", "Blg", "Bbv", "Bad", "Ehal", "Fpr", "Rint"), ordered = TRUE)
all_data$Species <- species_name_order

baseline_data <- all_data %>% filter(type == "D0")
# baseline_data$Metagenomics<-baseline_data$Metagenomics/100
# baseline_data$Prediction<-baseline_data$Prediction/100

intervention_data <- all_data %>% filter(type == "D14")
# intervention_data$Metagenomics<-intervention_data$Metagenomics/100
# intervention_data$Prediction<-intervention_data$Prediction/100
rsqred(intervention_data$Metagenomics,intervention_data$Prediction)   ## 


# baseline_lm_yx <- lm(Metagenomics ~ 0 + offset(1*Prediction), data=baseline_data)
# baseline_lm_yx_summary <-  summary(baseline_lm_yx)
# # glance(baseline_lm_yx)
# 
# baseline_lm <- lm(Metagenomics ~ Prediction, data=baseline_data)
# glance(baseline_lm)
baseline_lm_r2<-rsqred(baseline_data$Metagenomics,baseline_data$Prediction)   ## 

# baseline_lm_summary <-  summary(baseline_lm)
# baseline_lm_r2 <- baseline_lm_summary$adj.r.squared
# 
# baseline_cor <- cor.test(baseline_data$Metagenomics, baseline_data$Prediction)
# baseline_cor
# 
# baseline_RSS <- c(crossprod(baseline_lm$residuals))
# baseline_MSE <- baseline_RSS / length(baseline_lm$residuals)
# baseline_RMSE <- sqrt(baseline_MSE)
# baseline_RMSE
# 0.005141222

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


col_8 <- brewer.pal(8, "Set1")
col_8 <- add.alpha(col_8, alpha = 0.9)
col_used <- as.character(col_8[c(1:5, 7, 8)])
# col_used<-c("#b73276", "#228fa0", "#325780", "#3a8f62", "#b48939","#8f5360", "#754577", "#a7624f")
color_org<-c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#BC80BD", "#E6AB02", "#80B1D3")
col_used <- as.character(color_org[c(2:8)])
species_name_order <- factor(species_name, levels = c("Bfr", "Blg", "Bbv", "Bad", "Ehal", "Fpr", "Rint"), ordered = TRUE)

# baseline_lm_r2<-0.93

baseline_plot <- ggplot(baseline_data, aes(x=Prediction, y=Metagenomics)) +
  # geom_abline(intercept = 0, color="#4292c6", size=0.55) +
  geom_abline(intercept = 0, color="black", size=0.5) +
  geom_point(aes(color=Species), size=2.2) +
  scale_color_manual(values=col_used) +
  # geom_text(aes(x=0.021, y=0.097, label = paste("R^2 == ", round(baseline_lm_r2, 2))), parse=T, color="black",size=3,family="Helvetica") +
  # geom_text(aes(x=0.000021, y=0.27, label = paste("R^2 == ", round(baseline_lm_r2, 2))), parse=T, colour="black",size=3,family="Helvetica") +
  #coord_fixed() +
  # xlim(0, 10) +
  ylab("Experiment")+
  scale_x_continuous(limits=c(0,0.1),breaks=c(0,0.05,0.1),labels=scaleFUN)+
  # scale_x_log10(limits=c(0.000001,0.5),breaks=c(0.00001,0.01,0.5),labels=scaleFUN)+
  
  # ylim(0, 10) +
  scale_y_continuous(limits=c(0,0.1),breaks=c(0,0.05,0.1),labels=scaleFUN)+
  # scale_y_log10(limits=c(0.000001,0.5),breaks=c(0.00001,0.01,0.5),labels=scaleFUN)+
  
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

baseline_plot


# intervention_lm_yx <- lm(Metagenomics ~ 0 + offset(1*Prediction), data=intervention_data)
# intervention_lm_yx_summary <-  summary(intervention_lm_yx)
# 
# intervention_lm <- lm(Metagenomics ~ Prediction, data=intervention_data)
# intervention_lm_summary <-  summary(intervention_lm)
# intervention_lm_r2 <- intervention_lm_summary$adj.r.squared
# 
# intervention_cor <- cor.test(intervention_data$Metagenomics, intervention_data$Prediction)
# intervention_cor

# Pearson's product-moment correlation
# 
# data:  intervention_data$Metagenomics and intervention_data$Prediction
# t = 4.8324, df = 5, p-value = 0.004746
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.4880277 0.9864385
# sample estimates:
#       cor 
# 0.9075494 

# intervention_RSS <- c(crossprod(intervention_lm$residuals))
# intervention_MSE <- intervention_RSS / length(intervention_lm$residuals)
# intervention_RMSE <- sqrt(intervention_MSE)
# intervention_RMSE
# 0.007552208
# intervention_lm_r2<-0.78
intervention_lm_r2<-rsqred(intervention_data$Metagenomics,intervention_data$Prediction)   ## 


intervention_plot <- ggplot(intervention_data, aes(x=Prediction, y=Metagenomics)) +
  # geom_abline(intercept = 0, color="#4292c6", size=0.55) +
  geom_abline(intercept = 0, color="black", size=0.5) +
  geom_point(aes(color=Species), size=2.2) +
  scale_color_manual(values=col_used) +
  # geom_text(aes(x=0.021, y=0.097, label = paste("R^2 == ", round(intervention_lm_r2, 2))), parse=T, colour="black",size=3,family="Helvetica") +
  # geom_text(aes(x=0.000021, y=0.27, label = paste("R^2 == ", round(intervention_lm_r2, 2))), parse=T, colour="black",size=3,family="Helvetica") +
  
    #coord_fixed() +
  # xlim(0, 10) +
  # ylim(0, 10) +
  ylab("Experiment")+

  scale_x_continuous(limits=c(0,0.1),breaks=c(0,0.05,0.1),labels=scaleFUN)+
  # scale_x_log10(limits=c(0.000001,0.5),breaks=c(0.000001,0.01,0.5),labels=scaleFUN)+
  
  # ylim(0, 10) +
  scale_y_continuous(limits=c(0,0.1),breaks=c(0,0.05,0.1),labels=scaleFUN)+
  # scale_y_log10(limits=c(0.000001,0.5),breaks=c(0.000001,0.01,0.5),labels=scaleFUN)+
  
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
    axis.title.y = element_text(size = 10, colour="black",family = "Helvetica", vjust=0.018),
    axis.title.x = element_text(size = 10, colour="black",family = "Helvetica"),
    
    axis.ticks.length=unit(0.18, "line"),
    axis.ticks = element_line(colour = "black", size = 0.2),
    plot.margin = unit(c(0, 0, 0, 3), "pt")
    
    # axis.ticks = element_blank()
  )

intervention_plot

all_plot <- ggarrange(baseline_plot, intervention_plot, ncol=2, 
                      widths = c(1.1, 1.1), align="v")

# ggsave("./fig/baseline_intervention_plot_v4_0710.png", all_plot, height = 4.5, width = 9.5, units = "cm", dpi = 300, device = "png")

return(all_plot)

}
