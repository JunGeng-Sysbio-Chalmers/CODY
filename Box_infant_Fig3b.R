Box_infant_Fig3b <- function(exp_file,prediction_file) {
  
####################
# load data
####################
B_file <- "B_formula_mixed_8_ALL.txt"
B_data <- read.table(file=B_file, sep="\t", header=T, stringsAsFactors=F)

M4_file <- "M4_formula_mixed_8_ALL.txt"
M4_data <- read.table(file=M4_file, sep="\t", header=T, stringsAsFactors=F)

M12_file <- "M12_formula_mixed_8_ALL.txt"
M12_data <- read.table(file=M12_file, sep="\t", header=T, stringsAsFactors=F)

# merge 3 data set into one
all_data <- rbind(B_data, M4_data, M12_data)
all_data$Type<-factor(all_data$Type, levels=c("Newborn", "4 Months", "12 Months"))

all_data$RelativeAbundance <- all_data$RelativeAbundance + 1e-10

# read species info
species_file <- "Bac8_species.txt"
species_data <- read.table(file=species_file, sep="\t", header=F, stringsAsFactors=F)
species_name <- species_data[,1]
num_species <- nrow(species_data)

# select the main 4 species by species names
main5_data <- all_data %>% filter(Species %in% c("Bacteroides fragilis", "Bacteroides thetaiotaomicron", "Bifidobacterium longum", "Bifidobacterium breve", "Bifidobacterium adolescentis"))

other3_data <- all_data %>% filter(! Species %in% c("Bacteroides fragilis", "Bacteroides thetaiotaomicron", "Bifidobacterium longum", "Bifidobacterium breve", "Bifidobacterium adolescentis"))

# Eubacterium
Eha <- all_data %>% filter(Species %in% c("Eubacterium hallii"))
summary(Eha)
summary(Eha[which(Eha$Type=="Newborn"),])
summary(Eha[which(Eha$Type=="4 Months"),])
summary(Eha[which(Eha$Type=="12 Months"),])

####################
# plot data by groups
####################

all_summary <- aggregate(RelativeAbundance~ Type+Species, mean, data=all_data)
main5_summary <- aggregate(RelativeAbundance~ Type+Species, mean, data=main5_data)
other3_summary <- aggregate(RelativeAbundance~ Type+Species, mean, data=other3_data)

main5_colset <- data.frame(Species=unique(main5_data[,c('Species')]), A=1, B=1)
other3_colset <- data.frame(Species=unique(main5_data[,c('Species')]), A=1, B=1)

# change the long name
# Helper function for string wrapping. 
# Default 20 character target width.
# swr = function(string, nwrap=10) {
#   paste(strwrap(string, width=nwrap), collapse="\n")
# }
# swr = Vectorize(swr)
# 
# main5_data$Species = swr(main5_data$Species)
# main5_data$Species <- factor(main5_data$Species, labels = c("B. fragilis", "B. thetaiotaomicron", "B. breve", "B. longum"))
# main5_summary$Species <- factor(main5_summary$Species, labels = c("B. fragilis", "B. thetaiotaomicron", "B. breve", "B. longum"))
main5_data$Species <- factor(main5_data$Species, labels = c("Bfr", "Bth", "Bad", "Bbv", "Blg"))
main5_data$Species <- factor(main5_data$Species, levels = c("Bth", "Bfr", "Blg", "Bbv", "Bad"))
main5_summary$Species <- factor(main5_summary$Species, labels = c("Bfr", "Bth", "Bad", "Bbv", "Blg"))
main5_summary$Species <- factor(main5_summary$Species, levels = c("Bth", "Bfr", "Blg", "Bbv", "Bad"))

# main5_data$Species <- factor(main5_data$Species, labels = c("Bfr", "Bth", "Bad", "Bbv", "Blg"))
# main5_data$Species <- factor(main5_data$Species, levels = c("Bfr", "Blg", "Bbv", "Bad","Bth"))
# main5_summary$Species <- factor(main5_summary$Species, labels = c("Bfr", "Bth", "Bad", "Bbv", "Blg"))
# main5_summary$Species <- factor(main5_summary$Species, levels = c("Bfr", "Blg", "Bbv", "Bad","Bth"))


color_1225<-c("#e87163", "#fc9a8f",  "#eda244","#ceece1","#8DD3C7",  "#80B1D3", "#88cec2", "#007eb6")
color_12255<-c("#e2764b",  "#007db6", "#5faa85", "#a2c4d3",  "#dd8c63")
color_org<-c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#BC80BD", "#E6AB02", "#80B1D3")
# col_used <- as.character(color_org[c(2:8)])
# maxy<-max(main5_plot$RelativeAbundance)*1.75
# miny<-(maxy/50)*-1
# mody<-(maxy-maxy%%3)/3

main5_plot <- ggplot(main5_data, aes(y=RelativeAbundance, x=Type, shape=Type)) + 
  xlab("") + 
  ylab("Experiment") +
  stat_boxplot(aes(Type, RelativeAbundance,colour=Species),
               geom='errorbar', linetype=1, width = 0.3, size=0.6)+
  geom_jitter(size=1.25, position=position_jitter(0.1), alpha=0.45, aes(colour=Species), stroke=0.5, fill="white"
              # fill="#7d7d7d"
              # colour="#fcfcfc"
  ) +
  
  
  geom_boxplot(aes(color=Species),alpha=1, size=0.77, width=0.8, outlier.size = -1, fatten=NULL
               # linetype="dotted"
               ) +
  geom_errorbar(data = main5_summary,
                aes(ymin = RelativeAbundance, ymax = RelativeAbundance, colour=Species),
                width = 0.7,
                # col="#fb8072",
                # col="white",
                # color="#535353",
                # col="#f4f4f4",
                # alpha=0.8,
                size=1,linetype=1) +
  geom_line(data = main5_summary,
            aes(y = RelativeAbundance, x = Type, group=Species),
            # colour="#fb8072"
            colour = "#737373",
            # colour = "#114363",
            size=0.9,linetype=1
  ) +
  scale_shape_manual(values=c(21,24,22)) +
  facet_grid(.~Species) +
  scale_x_discrete(labels=c("Newborn" = "0", "4 Months" = "4", "12 Months" = "12")) +
  # scale_fill_manual(values=c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00")) +
  # scale_color_manual(values=c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00")) +
  # scale_fill_manual(values=c("#b73276", "#228fa0", "#325780", "#3a8f62", "#b48939")) +
  # scale_color_manual(values=c("#b73276", "#228fa0", "#325780", "#3a8f62", "#b48939")) +
  # scale_fill_manual(values=brewer.pal(n=8,"Dark2")[1:5]) +
  # scale_color_manual(values=brewer.pal(n=8,"Dark2")[1:5]) +
  scale_fill_manual(values=color_org[1:5]) +
  scale_color_manual(values=color_org[1:5]) +
  
  
  # scale_y_continuous(expand=c(0.03,0), limits=c(10^(-8), 1)) +
  # scale_y_continuous(expand=c(0.03,0), breaks=c(0, 0,.25, 0.5), limits=c(10^(-8), 0.7)) +
  scale_y_continuous(expand=c(0.03,0), breaks=c(0, 0.2,0.4,0.6), limits=c(10^(-8), 0.65)) +
  
    theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "#fcfcfc"),
        text = element_text(size=18, family="Helvetica", colour="black"), 
        panel.spacing = unit(0.15, "lines"),
        axis.text=element_text(size=18, family="Helvetica", colour="black"),
        axis.line.y=element_line(size=0.5),
        axis.line.x=element_line(size=0.5),
        legend.position="none", 
        legend.box = "horizontal",
        legend.background = element_blank(),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        strip.text.x = element_text(size=18),
        strip.background = element_rect(fill="#f3f3f3"),
        axis.title=element_text(size=18),
        axis.ticks.length=unit(0.275, "line"),
        axis.ticks = element_line(colour = "black", size = 0.45),
        # plot.margin = unit(c(0, 0, 0, 0.15), "pt")
        plot.margin = unit(c(0, 0.15, 0.1, 0.1), "pt")
        
          ) 
main5_plot

# other3_data$Species <- factor(other3_data$Species, labels = c("E. hallii", "F. prausnitzii", "R. intestinalis"))
# other3_summary$Species <- factor(other3_summary$Species, labels = c("E. hallii", "F. prausnitzii", "R. intestinalis"))
other3_data$Species <- factor(other3_data$Species, labels = c("Ehal", "Fpr", "Rint"))
other3_summary$Species <- factor(other3_summary$Species, labels = c("Ehal", "Fpr", "Rint"))

#other3_plot <- ggplot(other3_data, aes(y=RelativeAbundance, x=Type, fill=Species, shape=Type)) + 
  # other3_plot <- ggplot(other3_data, aes(y=RelativeAbundance, x=Type, shape=Type)) + 
  # xlab("") + 
  # ylab("Relative abundance") +
  # #geom_jitter(size=2, position=position_jitter(0.2), alpha=0.7, colour="grey2") +
  # geom_boxplot(aes(color=Species), size=0.6, width=0.5, outlier.size = 0.4, fatten=NULL) +
  # #geom_boxplot(color="grey2", size=0.5, width=0.4, outlier.size = 0.5, fatten=NULL) +
  # geom_line(data = other3_summary,
  #           aes(y = RelativeAbundance, x = Type, group=Species),
  #           colour="#fb8072"
  #           ) +
other3_plot <- ggplot(other3_data, aes(y=RelativeAbundance, x=Type, shape=Type)) + 
  xlab("") + 
  ylab("Experiment") +
  stat_boxplot(aes(Type, RelativeAbundance,colour=Species),
               geom='errorbar', linetype=1, width = 0.3, size=0.6)+
  geom_jitter(size=1.25, position=position_jitter(0.1), alpha=0.45, aes(colour=Species), stroke=0.5, fill="white"
              # fill="#7d7d7d"
              # colour="#fcfcfc"
  ) +
  
  
  geom_boxplot(aes(color=Species),alpha=1, size=0.77, width=0.8, outlier.size = -1, fatten=NULL
               # linetype="dotted"
  ) +
  geom_errorbar(data = other3_summary,
                aes(ymin = RelativeAbundance, ymax = RelativeAbundance, colour=Species),
                width = 0.7,
                # col="#fb8072",
                # col="white",
                # color="#535353",
                # col="#f4f4f4",
                # alpha=0.8,
                size=1,linetype=1) +
  geom_line(data = other3_summary,
            aes(y = RelativeAbundance, x = Type, group=Species),
            # colour="#fb8072"
            colour = "#737373",
            # colour = "#114363",
            size=0.9,linetype=1
  ) +
  scale_shape_manual(values=c(21,24,22)) +
  facet_grid(.~Species) +
  scale_x_discrete(labels=c("Newborn" = "0", "4 Months" = "4", "12 Months" = "12")) +
  # scale_fill_manual(values=c("#ffff33", "#a65628", "#f781bf")) +
  # scale_color_manual(values=c("#ffff33", "#a65628", "#f781bf")) +
  # scale_fill_manual(values=c("#ffcc42", "#a65628", "#f781bf")) +
  # scale_color_manual(values=c("#ffcc42", "#a65628", "#f781bf")) +
    # scale_fill_manual(values=c("#8f5360", "#754577", "#a7624f")) +
    # scale_color_manual(values=c("#8f5360", "#754577", "#a7624f")) +
  # scale_fill_manual(values=c("#BC80BD", "#E6AB02", "#80B1D3")) +
  # scale_color_manual(values=c("#BC80BD", "#E6AB02", "#80B1D3")) +
  scale_fill_manual(values=color_org[6:8]) +
  scale_color_manual(values=color_org[6:8]) +
  # scale_fill_manual(values=brewer.pal(n=8,"Dark2")[c(6,8,7)]) +
  # scale_color_manual(values=brewer.pal(n=8,"Dark2")[c(6,8,7)]) +

  #scale_y_continuous(expand=c(0.03,0), limits=c(10^(-8), 0.25), sec.axis = sec_axis(~.*1, name="Relative abundance")) +
  # scale_y_continuous(expand=c(0.03,0), limits=c(10^(-8), 0.25), sec.axis = sec_axis(~.*1)) +
    # scale_y_continuous(expand=c(0.03,0), limits=c(10^(-8), 0.15), position = "right") +
  scale_y_continuous(expand=c(0.03,0),  breaks=c(0, 0.05, 0.1), limits=c(10^(-8), 0.12), position = "right") +
  
  ylab("")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "#fcfcfc"),
        text = element_text(size=18,  colour="black", family = "Helvetica"), 
        panel.spacing = unit(0.15, "lines"),
        axis.text=element_text(size=18,  colour="black", family = "Helvetica"),
        axis.line.y=element_line(size=0.5),
        axis.line.x=element_line(size=0.5),
        legend.position="none", 
        legend.box = "horizontal",
        legend.background = element_blank(),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        strip.text.x = element_text(size=18),
        strip.background = element_rect(fill="#f3f3f3"),
        axis.title=element_text(size=20),
        # plot.margin = unit(c(0, 0, 0, 0.15), "pt")
        axis.ticks.length=unit(0.275, "line"),
        axis.ticks = element_line(colour = "black", size = 0.45),
        plot.margin = unit(c(0, 0, 0.1, 2), "pt")
        
          ) 
other3_plot


ggarrange(main5_plot+rremove("x.text"), other3_plot+rremove("ylab")+rremove("x.text"),
          nrow=1, ncol=2,
          widths=c(5.1,3))

# prediction_file <- "./data_prediction/Prediction.txt"  ## refer to "/Users/gejun/Desktop/Projects/Dynamic_Infant_Microbiome/Data_Visulization/20190422/Figure3/b/box/data_prediction/prediction.txt"

prediction_data <- read.table(file=prediction_file, sep=";", header=T, stringsAsFactors=F )
# prediction_data$Species<-rownames(prediction_data)
experiment_data<-prediction_data[c(1,4,7),]
experiment_data[1:3,2:ncol(experiment_data)]<-experiment_data[1:3,2:ncol(experiment_data)]/100
# rownames(experiment_data)[2:4]<-c('Newborn', '4 Month', '12 Month')
experiment_data$Time<-c('Newborn', '4 Months', '12 Months')

prediction_data<-prediction_data[c(1,3,6),]
prediction_data[1:3,2:ncol(experiment_data)]<-prediction_data[1:3,2:ncol(experiment_data)]/100
# rownames(prediction_data)[2:4]<-c('Newborn', '4 Month', '12 Month')
prediction_data$Time<-c('Newborn', '4 Months', '12 Months')

# prediction_data<-gather(prediction_data, Type, RA, Newborn:'12 Month')
prediction_data<-melt(prediction_data, id='Time')
colnames(prediction_data)[1:3]<-c('Type','Species','prediction')

experiment_data<-melt(experiment_data, id='Time')
colnames(experiment_data)[1:3]<-c('Type','Species','experiment')

# main5_prediction <- prediction_data %>% filter(Species %in% c("Bacteroides fragilis", "Bacteroides thetaiotaomicron", "Bifidobacterium longum", "Bifidobacterium breve", "Bifidobacterium adolescentis"))
# other3_prediction <- prediction_data %>% filter(! Species %in% c("Bacteroides fragilis", "Bacteroides thetaiotaomicron", "Bifidobacterium longum", "Bifidobacterium breve", "Bifidobacterium adolescentis"))

main5_prediction <- prediction_data %>% filter(Species %in% c("Bacteroides_fragilis", "Bacteroides_theta", "Bifido_longum", "Bifido_breve", "Bifido_adolescentis"))
other3_prediction <- prediction_data %>% filter(! Species %in% c("Bacteroides_fragilis", "Bacteroides_theta", "Bifido_longum", "Bifido_breve", "Bifido_adolescentis"))


main5_prediction$Type <- factor(main5_prediction$Type, levels=c("Newborn", "4 Months", "12 Months"))
other3_prediction$Type <- factor(other3_prediction$Type, levels=c("Newborn", "4 Months", "12 Months"))

# main5_prediction$Species <- factor(main5_prediction$Species, labels = c("Bfr", "Bth", "Bad", "Bbv", "Blg"))
main5_prediction$Species <- factor(main5_prediction$Species, labels = c("Bth", "Bfr", "Blg", "Bbv", "Bad"))
main5_prediction$Species <- factor(main5_prediction$Species, levels = c("Bth", "Bfr", "Blg", "Bbv", "Bad"))
other3_prediction$Species <- factor(other3_prediction$Species, labels = c("Ehal", "Fpr", "Rint"))
other3_prediction$Species <- factor(other3_prediction$Species, labels = c("Ehal", "Fpr", "Rint"))

main5_predict_plot <- ggplot(main5_prediction, aes(y=prediction, x=Type, colour=Species, fill=Species, shape=Type)) + 
  xlab("") + 
  ylab("Prediction") +
  # geom_bar(position="dodge", stat="identity", alpha=0.7, width=0.4, colour="black") +
  geom_bar(position="dodge", stat="identity", alpha=1, width=0.61) +
  
    facet_grid(.~Species) +
  # geom_line(data = main5_prediction,
  #           aes(y = prediction, x = Type, group=Species),
  #           # colour="#fb8072"
  #           colour = "#535353",
  #           # colour = "#114363",
  #           size=1.5,linetype=1
  # ) +
  scale_x_discrete(labels=c("Newborn" = "0", "4 Months" = "4", "12 Months" = "12")) +
  # scale_fill_manual(values=c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00")) +
  # # scale_color_manual(values=c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00")) +
  # scale_fill_manual(values=brewer.pal(n=8,"Dark2")[1:5]) +
  # scale_color_manual(values=brewer.pal(n=8,"Dark2")[1:5]) +
  scale_fill_manual(values=color_org[1:5]) +
  scale_color_manual(values=color_org[1:5]) +
  # scale_fill_manual(values=c("#b73276", "#228fa0", "#325780", "#3a8f62", "#b48939")) +
  # scale_color_manual(values=c("#b73276", "#228fa0", "#325780", "#3a8f62", "#b48939")) +
  scale_y_continuous(expand=c(0.003,0), breaks = c(0,0.1,0.2,0.3), limits=c(0, 0.32)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "#fcfcfc"),
        text = element_text(size=18, family="Helvetica", colour="black"), 
        panel.spacing = unit(0.15, "lines"),
        axis.text=element_text(size=16, colour="black", family = "Helvetica"),
        axis.line.y=element_line(size=0.5),
        axis.line.x=element_line(size=0.5),
        legend.position="none", 
        legend.box = "horizontal",
        legend.background = element_blank(),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        strip.text.x = element_text(size=18),
        strip.background = element_rect(fill="#f3f3f3"),
        # plot.margin = unit(c(-1, 0, 0, 0.1), "pt"),
        plot.margin = unit(c(-17, 0.15, 0, 0.1), "pt"),
        axis.text.x = element_text(size = 18,vjust=-0.5),
  
        axis.ticks.length=unit(0.275, "line"),
        axis.ticks = element_line(colour = "black", size = 0.45)
  ) 
main5_predict_plot

other3_predict_plot <- ggplot(other3_prediction, aes(y=prediction, x=Type, colour=Species,fill=Species, shape=Type)) + 
  xlab("") + 
  ylab("Prediction") +
  geom_bar(position="dodge", stat="identity", alpha=1, width=0.61) + # size=1, for the border width 
  # geom_line(data = other3_prediction,
  #           aes(y = prediction, x = Type, group=Species),
  #           # colour="#fb8072"
  #           colour = "#535353",
  #           # colour = "#114363",
  #           size=1.5,linetype=1
  # ) +
  facet_grid(.~Species) +
  scale_x_discrete(labels=c("Newborn" = "0", "4 Months" = "4", "12 Months" = "12")) +
  # scale_fill_manual(values=c("#ffff33", "#a65628", "#f781bf")) +
  # scale_fill_manual(values=c("#ffcc42", "#a65628", "#f781bf")) +
  # scale_color_manual(values=c("#ffcc42", "#a65628", "#f781bf")) +
  # scale_fill_manual(values=c("#BC80BD", "#E6AB02", "#80B1D3")) +
  # scale_color_manual(values=c("#BC80BD", "#E6AB02", "#80B1D3")) +
  scale_fill_manual(values=color_org[6:8]) +
  scale_color_manual(values=color_org[6:8]) +
  # scale_y_continuous(expand=c(0.003,0), limits=c(0, 0.05), sec.axis = sec_axis(~.*1)) +
  scale_y_continuous(expand=c(0.003,0), breaks=c(0,0.05),limits=c(0, 0.05), position="right") +

  ylab("")+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "#fcfcfc"),
        text = element_text(size=18, colour="black", family="Helvetica"), 
        panel.spacing = unit(0.15, "lines"),
        axis.text=element_text(size=18,colour="black", family = "Helvetica"),
        axis.line.y=element_line(size=0.5),
        axis.line.x=element_line(size=0.5),
        legend.position="none", 
        legend.box = "horizontal",
        legend.background = element_blank(),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        strip.text.x = element_text(size=18),
        strip.background = element_rect(fill="#f3f3f3"),
        axis.title=element_text(size=20),
        axis.ticks.length=unit(0.275, "line"),
        axis.ticks = element_line(colour = "black", size = 0.45),
        axis.text.x = element_text(size = 18,vjust=-0.5),
        
        # plot.margin = unit(c(-1, 0, 0, 0.1), "pt")
        plot.margin = unit(c(-17, 0, 0, 2), "pt") ## top, right, bottom, left
        
          ) 

other3_predict_plot


multi_plot <- ggarrange(main5_plot+rremove("x.text")+rremove("x.ticks"), 
          other3_plot+rremove("ylab")+rremove("x.text")+rremove("x.ticks"), 
          main5_predict_plot, 
          other3_predict_plot+rremove("ylab"),
          nrow=2, ncol=2,
          widths=c(5,3.1),
          align="v"
          )

return(multi_plot)


# ggexport(multi_plot, filename="multi_plot_formula_mixed_box_1226.jpg", width=1480, height=790, res=120)
# ggsave(plot=multi_plot, filename="multi_plot_formula_mixed_box_1226.pdf", width=29.6, height=15.8, units="cm")
# write_tsv(all_summary, "Bac8_mean_formula_mixed.txt")
# write_tsv(experiment_data,"experimentdata_20181031.txt")
# write_tsv(prediction_data,"prediction_20181031.txt")

}