Draw_Fig3a <- function(prediction_exp_file) {
  
  color_4M<-brewer.pal(n=8,"Dark2")[1:5]
  color_12M<-c(brewer.pal(n=8,"Dark2")[1:5],c("#BC80BD", "#E6AB02", "#80B1D3"))
  # exp_file <- "experimentdata_20181031.txt"
  # prd_file<-"Preidiction_RA_exp_mdf.csv"
  data<- read.table(file=prediction_exp_file, sep=";", header=T, stringsAsFactors=F )

  # exp_data <- read.table(file = exp_file, sep = "\t", header = T, stringsAsFactors = F)
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
  # all_data$Species <- factor(all_data$Species, labels = c(
  #   "Bfr", "Bth", "Bad", "Bbv",
  #   "Blg", "Ehal", "Fpr", "Rint"
  # ))
  all_data$Species <- factor(all_data$Species, labels = c(
    "Bth", "Bfr", "Blg",  "Bbv","Bad",
    "Ehal", "Fpr", "Rint"
  ))
  all_data$Species <- factor(all_data$Species, levels = c(
    "Bth", "Bfr", "Blg", "Bbv",
    "Bad", "Ehal", "Fpr", "Rint"
  ))
  
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
  color_4M<-add.alpha(color_4M, alpha=0.75)
  color_12M<-add.alpha(color_12M, alpha=0.75)
  

  # draw 4M data
 
   M4_data <- all_data %>% filter(Type == "4 Months") %>% filter(prediction != 0)
  # M4_others <- data.frame(Type="4 Months", Species="Others", experiment=1-sum(M4_data[,3]), prediction=1-sum(M4_data[,4]))
  # M4_data_new <- rbind(M4_data, M4_others)
  M4_data_new <- with(M4_data, M4_data[order(Species), ])
  
  M4_sum_experiment <- round(sum(M4_data$experiment) * 100, 2)
  M4_sum_experiment_label <- paste(M4_sum_experiment, "%", sep = " ")
  M4_sum_prediction <- round(sum(M4_data$prediction) * 100, 2)
  M4_sum_prediction_label <- paste(M4_sum_prediction, "%", sep = " ")
  
  
  ## below is to convert plots to ggplot and return
donut<-function(){
  # opar <- par(no.readonly = T)
  # # pdf("4M_two_donut_all_v2_1226_75alpha.pdf", width = 11 / 2.54, height = 5 / 2.54, pointsize = 7)
  # 
  # par(mfrow = c(2, 2), mar = c(0, 0, 0, 0))
  
  pie(
    M4_data_new[, 3],
    labels = M4_data_new$Species,
    # col = c(col_8[1:5]),
    col = color_4M,alpha=0.7,
    
    border="white",
    clockwise = T,
    init.angle = 0,cex=1.5
  )
  par(new = TRUE)
  pie(1, radius = 0.4, col = "white", border = "white", lwd = 2, labels = "")
  # text(0, 0, labels = M4_sum_experiment_label, cex = 1, font = 2)
  text(0, 0, labels = "4M", cex = 2, family ="Helvetica")
  # dev.off()


  pie(
    M4_data_new[, 4],
    labels = M4_data_new$Species,
    # col = c(col_8[1:5]),
    col = color_4M,
    border="white",
    clockwise = T,
    init.angle = 0,cex=1.5
  )
  par(new = TRUE)
  pie(1, radius = 0.4, col = "white", border = "white", lwd = 2, labels = "")
  # text(0, 0, labels = M4_sum_prediction_label, cex = 1, font = 2)
  # text(0, 0, labels = "BF_pre", cex = 1, font =list(family="Arial", face=1))
  text(0, 0, labels = "4M", cex = 2, family ="Helvetica")
  

  
  M12_data <- all_data %>% filter(Type == "12 Months") %>% filter(prediction != 0)
  # M12_others <- data.frame(Type="12 Months", Species="Others", experiment=1-sum(M12_data[,3]), prediction=1-sum(M12_data[,4]))
  # M12_data_new <- rbind(M12_data, M12_others)
  M12_data_new <- with(M12_data, M12_data[order(Species), ])
  
  M12_sum_experiment <- round(sum(M12_data$experiment) * 100, 2)
  M12_sum_experiment_label <- paste(M12_sum_experiment, "%", sep = " ")
  M12_sum_prediction <- round(sum(M12_data$prediction) * 100, 2)
  M12_sum_prediction_label <- paste(M12_sum_prediction, "%", sep = " ")
  

  # opar <- par(no.readonly = T)
  # pdf("12M_two_donut_all_v2_1226_75.pdf", width = 11 / 2.54, height = 5 / 2.54, pointsize = 7)
  # par(mfrow = c(1, 2), mar = c(0, 0, 0, 0))
  pie(
    M12_data_new[, 3],
    labels = M12_data_new$Species,
    # col = c(col_8, "gray80"),
    col = color_12M,
    border="white",
    clockwise = T,cex=1.5
  )
  par(new = TRUE)
  pie(1, radius = 0.4, col = "white", border = "white", lwd = 2, labels = "")
  # text(0, 0, labels = M12_sum_experiment_label, cex = 1, font = 2)
  text(0, 0, labels = "12M", cex = 2, family="Helvetica")

  # dev.off()

   pie(
    M12_data_new[, 4],
    labels = M12_data_new$Species,
    # col = c(col_8, "gray80"),
    col = color_12M,
    border="white",
    clockwise = T,cex=1.5
  )
  par(new = TRUE)
  pie(1, radius = 0.4, col = "white", border = "white", lwd = 2, labels = "")
  # text(0, 0, labels = M12_sum_experiment_label, cex = 1, font = 2)
  text(0, 0, labels = "12M", cex = 2, family="Helvetica")
}
print(donut())
a<-list(donut)
a
p1<-as.grob(function() donut_4m_1())
p2<-as.grob(function() donut_4m_2())
p3<-as.grob(function() donut_12m_1())
p4<-as.grob(function() donut_12m_2())
# plot_grid(p1, p2, p3, p4, ncol=2, labels=c("Exp.M4","Pre.M4","Exp.M12","Pre.M12"))
multi_plot <- ggarrange(p1,p2,p3,p4,
                        nrow=2, ncol=2,
                        # widths=c(5,3.1),
                        heights=c(2.5,2.5),
                        align="v")
multi_plot
return(multi_plot)
}
