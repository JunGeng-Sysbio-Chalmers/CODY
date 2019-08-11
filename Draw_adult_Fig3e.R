Draw_adult_Fig3e <- function(exp_file,prediction_data) {
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

col_8 <- brewer.pal(8, "Set1")
col_8 <- add.alpha(col_8, alpha = 0.9)
col_used_4 <- as.character(col_8[c(1:3)])
col_used_3 <- as.character(col_8[c(4, 5, 7, 8)])


data_file <- "Bac7_sample_ra.txt"
all_data <- read_tsv(data_file)

# only keep data at timepoint 0 & 14
all_data <- all_data %>% filter(Time == 0 | Time == 14)

# select the main 4 strain by strain names
main5_data <- all_data %>% filter(strain %in% c("Bfr", "Blg", "Bbv"))

other3_data <- all_data %>% filter(! strain %in% c("Bfr", "Blg", "Bbv"))

# Eubacterium
# Eha <- all_data %>% filter(strain %in% c("Eubacterium hallii"))
# summary(Eha)
# summary(Eha[which(Eha$Type=="Newborn"),])
# summary(Eha[which(Eha$Type=="4 Months"),])
# summary(Eha[which(Eha$Type=="12 Months"),])

####################
# plot data by groups
####################

# main5_boxplot <- ggplot(aes(y=RelativeAbundance, x=Type, fill=strain, col=strain), data=min4_data)

all_summary <- aggregate(RA~ Time+strain, mean, data=all_data)
main5_summary <- aggregate(RA~ Time+strain, mean, data=main5_data)
other3_summary <- aggregate(RA~ Time+strain, mean, data=other3_data)

main5_colset <- data.frame(strain=unique(main5_data[,c('strain')]), A=1, B=1)
other3_colset <- data.frame(strain=unique(main5_data[,c('strain')]), A=1, B=1)

color_7spc<-c(brewer.pal(n=8,"Dark2")[2:5],"#BC80BD", "#E6AB02", "#80B1D3")

# change the long name
# Helper function for string wrapping. 
# Default 20 character target width.
# swr = function(string, nwrap=10) {
#   paste(strwrap(string, width=nwrap), collapse="\n")
# }
# swr = Vectorize(swr)
# 
# main5_data$strain = swr(main5_data$strain)
# main5_data$strain <- factor(main5_data$strain, labels = c("B. fragilis", "B. thetaiotaomicron", "B. breve", "B. longum"))
# main5_summary$strain <- factor(main5_summary$strain, labels = c("B. fragilis", "B. thetaiotaomicron", "B. breve", "B. longum"))
main5_data$strain <- factor(main5_data$strain, levels = c("Bfr", "Blg", "Bbv"), order=TRUE)
#main5_data$strain <- factor(main5_data$strain, levels = c("Bfr", "Blg", "Bbv", "Bad"))

#main5_summary$strain <- factor(main5_summary$strain, labels = c("Bfr",  "Bad", "Bbv", "Blg"))
main5_summary$strain <- factor(main5_summary$strain, levels = c( "Bfr", "Blg", "Bbv"), order=TRUE)

main5_data$Time <- factor(main5_data$Time)
main5_summary$Time <- factor(main5_summary$Time)

scaleFUN <- function(x) sprintf("%.2f", x)
scaleFUN1 <- function(x) sprintf("%.1f", x)

#main5_plot <- ggplot(main5_data, aes(y=RelativeAbundance, x=Type, fill=strain, shape=Type)) + 
main5_plot <- ggplot(main5_data, aes(y=RA, x=Time, shape=Time)) + 
  xlab("") + 
  ylab("Experiment") +
  stat_boxplot(aes(Time, RA,colour=strain),
               geom='errorbar', linetype=1, width = 0.25, size=0.4)+
  geom_jitter(size=1, position=position_jitter(0.1), alpha=0.45, aes(colour=strain), stroke=0.45, fill="white"
              # fill="#7d7d7d"
              # colour="#fcfcfc"
  ) +
  
  geom_boxplot(aes(color=strain),alpha=1, size=0.5, width=0.6, outlier.size = -1, 
               # weight=0.1,
               fatten=NULL
               # linetype="dotted"
  ) +
  geom_errorbar(data = main5_summary,
                aes(ymin = RA, ymax = RA, colour=strain),
                width = 0.65,
                # col="#fb8072",
                # col="white",
                # color="#535353",
                # col="#f4f4f4",
                # alpha=0.8,
                size=1,linetype=1) +
  scale_shape_manual(values=c(21,24,22)) +
  geom_line(data = main5_summary,
            aes(y = RA, x = Time, group=strain),
            # colour="#fb8072"
            colour = "#737373",
            # colour = "#114363",
            size=0.7,linetype=1
  ) +
  facet_grid(.~strain) +
  #scale_x_discrete(labels=c("Newborn" = "0", "4 Months" = "4", "12 Months" = "12")) +
  #scale_fill_manual(values=c( "#377eb8", "#4daf4a", "#984ea3", "#ff7f00")) +
  #scale_color_manual(values=c( "#377eb8", "#4daf4a", "#984ea3", "#ff7f00")) +
  # scale_fill_manual(values=col_used_4) +
  # scale_color_manual(values=col_used_4) +             
  scale_x_discrete(labels=c("0" = "Day0", "14" = "Day14")) +
  # scale_fill_manual(values=c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00")) +
  # scale_color_manual(values=c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00")) +
  # scale_fill_manual(values=c("#b73276", "#228fa0", "#325780", "#3a8f62", "#b48939")) +
  # scale_color_manual(values=c("#b73276", "#228fa0", "#325780", "#3a8f62", "#b48939")) +
  scale_fill_manual(values=color_7spc[1:3]) +
  scale_color_manual(values=color_7spc[1:3]) +
  
  scale_y_continuous(expand=c(0.003,0), breaks = c(0,0.01,0.02), limits=c(0, 0.02),labels=scaleFUN) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "#fcfcfc"),
        text = element_text(size=12, family="Helvetica"), 
        panel.spacing = unit(0.15, "lines"),
        axis.text=element_text(size=12, colour="black", family = "Helvetica"),
        axis.line.y=element_line(size=0.33),
        axis.line.x=element_line(size=0.33),
        legend.position="none", 
        legend.box = "horizontal",
        legend.background = element_blank(),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        strip.text.x = element_text(size=13),
        strip.background = element_rect(fill="#f3f3f3"),
        axis.title=element_text(size=14),
        axis.ticks.length=unit(0.25, "line"),
        axis.ticks = element_line(colour = "black", size = 0.3),
        plot.margin = unit(c(0, 0, 0, 1.5), "pt")
  ) 

main5_plot

# other3_data$strain <- factor(other3_data$strain, labels = c("E. hallii", "F. prausnitzii", "R. intestinalis"))
# other3_summary$strain <- factor(other3_summary$strain, labels = c("E. hallii", "F. prausnitzii", "R. intestinalis"))
other3_data$strain <- factor(other3_data$strain, levels = c("Bad", "Ehal", "Fpr", "Rint"), order=TRUE)
other3_summary$strain <- factor(other3_summary$strain, levels = c("Bad", "Ehal", "Fpr", "Rint"), order=TRUE)

other3_data$Time <- factor(other3_data$Time)
other3_summary$Time <- factor(other3_summary$Time)


  other3_plot <- ggplot(other3_data, aes(y=RA, x=Time, shape=Time)) + 
  xlab("") + 
  ylab("Experiment RA") +
    stat_boxplot(aes(Time, RA,colour=strain),
                 geom='errorbar', linetype=1, width = 0.27, size=0.5)+
    geom_jitter(size=1.2, position=position_jitter(0.1), alpha=0.45, aes(colour=strain), stroke=0.5, fill="white"
                # fill="#7d7d7d"
                # colour="#fcfcfc"
    ) +
    
    
    geom_boxplot(aes(color=strain),alpha=1, size=0.5, width=0.6, outlier.size = -1, fatten=NULL
                 # linetype="dotted"
    ) +
    geom_errorbar(data = other3_summary,
                  aes(ymin = RA, ymax = RA, colour=strain),
                  width = 0.65,
                  # col="#fb8072",
                  # col="white",
                  # color="#535353",
                  # col="#f4f4f4",
                  # alpha=0.8,
                  size=1,linetype=1) +
    geom_line(data = other3_summary,
              aes(y = RA, x = Time, group=strain),
              # colour="#fb8072"
              colour = "#737373",
              # colour = "#114363",
              size=0.8,linetype=1
    ) +
    scale_shape_manual(values=c(21,24,22)) +
    facet_grid(.~strain) +
    scale_x_discrete(labels=c("0" = "BSL", "14" = "ITV", "12 Months" = "12")) +
    
    scale_fill_manual(values=color_7spc[4:7]) +
    scale_color_manual(values=color_7spc[4:7]) +
    
    scale_y_continuous(expand=c(0.03,0),  breaks=c(0, 0.05, 0.1), limits=c(10^(-8), 0.1), position = "right",labels=scaleFUN1) +
    
    ylab("")+
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "#fcfcfc"),
          text = element_text(size=12, family="Helvetica"), 
          panel.spacing = unit(0.15, "lines"),
          axis.text=element_text(size=12, colour="black", family = "Helvetica"),
          axis.line.y=element_line(size=0.33),
          axis.line.x=element_line(size=0.33),
          legend.position="none", 
          legend.box = "horizontal",
          legend.background = element_blank(),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          strip.text.x = element_text(size=13),
          strip.background = element_rect(fill="#f3f3f3"),
          axis.title=element_text(size=14),
          axis.ticks.length=unit(0.25, "line"),
          axis.ticks = element_line(colour = "black", size = 0.3),
          plot.margin = unit(c(0, 0, 0, 1.5), "pt")
    ) 
other3_plot

ggarrange(main5_plot+rremove("x.text"), other3_plot+rremove("ylab")+rremove("x.text"),
          nrow=1, ncol=2,
          widths=c(3,4))

# prediction_data <- read.table(file=prediction_file, sep="\t", header=T, stringsAsFactors=F)
prediction_data<-read.table(file="prediction_format.csv", sep=",", header=T, stringsAsFactors=F )
prediction_data[,2:ncol(prediction_data)]<-prediction_data[,2:ncol(prediction_data)]
# reformat data
colnames(prediction_data)<-c("Time","Bfr","Blg","Bbv","Bad","Ehal","Fpr","Rint")
prediction_data$Time<-c("D0","D14")
prediction_data <- prediction_data %>% gather("strain", "prediction", -Time)

main5_prediction <- prediction_data %>% filter(strain %in% c("Bfr", "Blg", "Bbv"))
other3_prediction <- prediction_data %>% filter(! strain %in% c("Bfr", "Blg", "Bbv"))
max_value5<-max(main5_prediction[,3])
if(max_value5*1.75<0.1){
  max_lab5<-round(max_value5*1.75,2)
}else if(max_value5*1.75>0.1){
  max_lab5<-round(max_value5*1.75,1)
}
max_value3<-max(other3_prediction[,3])
if(max_value3*1.75<0.1){
  max_lab3<-round(max_value3*1.75,2)
}else if(max_value3*1.75>0.1){
  max_lab3<-round(max_value3*1.75,1)
}
main5_prediction$Time <- factor(main5_prediction$Time)
other3_prediction$Time <- factor(other3_prediction$Time)

main5_prediction$strain <- factor(main5_prediction$strain, levels = c("Bfr", "Blg", "Bbv"), order=TRUE)


main5_predict_plot <- ggplot(main5_prediction, aes(y=prediction, x=Time, fill=strain, shape=Time)) + 
  xlab("") + 
  ylab("Prediction") +
  geom_bar(position="dodge", stat="identity", alpha=1, width=0.56) +
  
  facet_grid(.~strain) +
  # geom_line(data = main5_prediction,
  #           aes(y = prediction, x = Type, group=Species),
  #           # colour="#fb8072"
  #           colour = "#535353",
  #           # colour = "#114363",
  #           size=1.5,linetype=1
  # ) +
  # scale_x_discrete(labels=c("0" = "BSL", "14" = "ITV", "12 Months" = "12")) +
  # scale_x_discrete(labels=c("0" = "Day0", "14" = "Day14")) +
  
  scale_fill_manual(values=color_7spc[1:3]) +
  scale_color_manual(values=color_7spc[1:3]) +
  labs(fill ="Left axis")+ guides(colour=FALSE)+
  
  # scale_fill_manual(values=c("#b73276", "#228fa0", "#325780", "#3a8f62", "#b48939")) +
  # scale_color_manual(values=c("#b73276", "#228fa0", "#325780", "#3a8f62", "#b48939")) +
  # scale_y_continuous(expand=c(0.003,0), breaks = c(0,0.01,0.02), limits=c(0, 0.02),labels=scaleFUN) +
  scale_y_continuous(expand=c(0.003,0), breaks = c(0,max_lab5/2,max_lab5), limits=c(0, max_lab5)) +
  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "#fcfcfc"),
        text = element_text(size=15, family="Helvetica"), 
        panel.spacing = unit(0.15, "lines"),
        axis.text=element_text(size=12, colour="black", family = "Helvetica"),
        axis.line.y=element_line(size=0.33),
        axis.line.x=element_line(size=0.33),
        # legend.position="none", 
        legend.position="bottom", 
        legend.title = element_text(size=15,colour="black", family = "Helvetica"),
        
        legend.box = "horizontal",
        legend.background = element_blank(),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        strip.text.x = element_text(size=13),
        strip.background = element_rect(fill="#f3f3f3"),
        axis.title=element_text(size=14),
        axis.ticks.length=unit(0.25, "line"),
        axis.text.x = element_text(size = 12,vjust=0.5,angle=45),
        
        axis.ticks = element_line(colour = "black", size = 0.3),
        plot.margin = unit(c(-1, 1.5, 0, 1.5), "pt")
  ) 
main5_predict_plot

other3_predict_plot <- ggplot(other3_prediction, aes(y=prediction, x=Time, fill=strain, shape=Time)) + 
  xlab("") + 
  ylab("Prediction RA") +
  geom_bar(position="dodge", stat="identity", alpha=1, width=0.56) + # size=1, for the border width 
  # geom_line(data = other3_prediction,
  #           aes(y = prediction, x = Type, group=Species),
  #           # colour="#fb8072"
  #           colour = "#535353",
  #           # colour = "#114363",
  #           size=1.5,linetype=1
  # ) +
  facet_grid(.~strain) +
  # scale_x_discrete(labels=c("0" = "Day0", "14" = "Day14")) +
  
  scale_fill_manual(values=color_7spc[4:7]) +
  scale_color_manual(values=color_7spc[4:7]) +
  # scale_y_continuous(expand=c(0.003,0), breaks=c(0,0.05,0.1),limits=c(0, 0.1), position="right",labels=scaleFUN1) +
  # scale_y_continuous(expand=c(0.003,0), breaks=c(0,0.05,0.1),limits=c(0, 0.1), position="right") +
  # scale_y_continuous(expand=c(0.003,0), limits=c(0, 0.1), position="right") +
  scale_y_continuous(expand=c(0.003,0), breaks = c(0,max_lab3/2,max_lab3), limits=c(0, max_lab3), position="right") +
  
  
  ylab("")+
  labs(fill ="Right axis")+ guides(colour=FALSE)+
  
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "#fcfcfc"),
        text = element_text(size=15, family="Helvetica", colour="black"), 
        panel.spacing = unit(0.15, "lines"),
        axis.text=element_text(size=12, colour="black", family = "Helvetica"),
        axis.line.y=element_line(size=0.33),
        axis.line.x=element_line(size=0.33),
        legend.position="bottom", 
        legend.title = element_text(size=15,colour="black", family = "Helvetica"),
        
        legend.box = "horizontal",
        legend.background = element_blank(),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        strip.text.x = element_text(size=13),
        strip.background = element_rect(fill="#f3f3f3"),
        axis.title=element_text(size=16),
        axis.ticks.length=unit(0.25, "line"),
        axis.ticks = element_line(colour = "black", size = 0.3),
        axis.text.x = element_text(size = 12,vjust=0.5,angle=45),
        
        plot.margin = unit(c(-1, 0, 0, 0), "pt")
  ) 

other3_predict_plot


# multi_plot <- ggarrange(main5_plot+rremove("x.text")+rremove("x.ticks"), 
#           other3_plot+rremove("ylab")+rremove("x.text")+rremove("x.ticks"), 
#           main5_predict_plot, 
#           other3_predict_plot+rremove("ylab"),
#           nrow=2, ncol=2,
#           widths=c(3.5,4),
#           align="v"
#           )

multi_plot <- ggarrange(
                        main5_predict_plot, 
                        other3_predict_plot+rremove("ylab"),
                        nrow=2, ncol=2,
                        widths=c(3.5,4)
                        # align="v"
)

return(multi_plot)

}

