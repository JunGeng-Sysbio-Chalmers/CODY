Tanks_cmp <- function(data,header,ytext) {
# data<-data_lumen_microbial_all_SF
# bacteria_12m_1102 <- read_excel("~/Desktop/Projects/Dynamic_Infant_Microbiome/Data_Visulization/20181027/Figure4/d/bacteria_12m_1102.xlsx")
# bacteria<-bacteria_12m_1102[1:4,]
# colnames(bacteria)[1]<-"Position"
# mets_12m<-read_excel("~/Desktop/Projects/Dynamic_Infant_Microbiome/Data_Visulization/20181027/Figure4/d/mets_12m_1102.xlsx")
# mets<-mets_12m[1:4,c(1,6,10,5)]   ## prop, but, ac
# colnames(mets)[1]<-"Position"
  Colon_Site<-c("I","II","III","IV")
  col_8<-c(brewer.pal(n=8,"Dark2")[1:5],c("#BC80BD", "#E6AB02", "#80B1D3"))
  
data$Position<-Colon_Site
# ## 4month
# bacteria_4m_1102 <- read_excel("~/Desktop/Projects/Dynamic_Infant_Microbiome/Data_Visulization/20181027/Figure4/d/bacteria_4m_1102.xlsx")
# bacteria<-bacteria_4m_1102[1:4,]
# colnames(bacteria)[1]<-"Position"
# mets_4m<-read_excel("~/Desktop/Projects/Dynamic_Infant_Microbiome/Data_Visulization/20181027/Figure4/d/mets_4m_1102.xlsx")
# mets<-mets_4m[1:4,c(1,6,5)]
# colnames(mets)[1]<-"Position"

data$Position <- factor(data$Position, levels=unique(data$Position))
colnames(data)[1:ncol(data)-1]<-str_replace(colnames(data)[1:ncol(data)-1],'BIOMASS_','')
colnames(data)[1:ncol(data)-1]<-str_replace(colnames(data)[1:ncol(data)-1],'Bacteroides ','B.')
colnames(data)[1:ncol(data)-1]<-str_replace(colnames(data)[1:ncol(data)-1],'Bifido ','B.')
colnames(data)[1:ncol(data)-1]<-str_replace(colnames(data)[1:ncol(data)-1],'Eubacterium ','E.')
colnames(data)[1:ncol(data)-1]<-str_replace(colnames(data)[1:ncol(data)-1],'Faecalibacterium ','F.')
colnames(data)[1:ncol(data)-1]<-str_replace(colnames(data)[1:ncol(data)-1],'Roseburia ','R.')

# data<-data[c(1,7,5,6,4,2,9,8,3)]
# colnames(data)<-c("Position","Ehallii" ,"Bbr", "Bad", "Blg", "Bth","Rint", "Fpr","Bfr"   )
# bacteria<-bacteria[c(1,6,2,5,3,4)]
# colnames(bacteria)<-c("Position","Bad" ,"Bth", "Bbr", "Bfr", "Blg"  )
# mets$Position <- factor(mets$Position, levels=unique(mets$Position))
# mets<-melt(mets)
maxy<-max(apply(data[,1:ncol(data)-1],2,max))*1.75
miny<-(maxy/50)*-1
mody<-(maxy-maxy%%3)/3
data<-melt(data)




# col_8<-col_8[c(5,1,4,2,3)]

# dat <- data.frame(
#   position_colon = factor(c("Acending","Transverse","Descending","Sigmoid"), 
#                           levels=c("Acending","Transverse","Descending","Sigmoid")),
#   total_bill = c(14.89, 17.23)
# )
# levels(position_colon)=c('ascending','transverse','descending',"sigmoid")

# qplot(name,score,data = a, geom = 'bar')
colnames(data)[2]<-"Species"
p <-ggplot(data, aes(x=Position,y=value,fill=Species))
# p<-p + theme_minimal()
p<-p+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.1,colour="#5e5e5e"),
           # axis.text.x = element_text(angle = 13, hjust = 1))
           axis.text.x = element_text(angle = 0, hjust = 0.2))

p<-p+theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  # panel.border = element_rect(
  #   colour = "#5e5e5e",
  #   fill = NA,
  #   size = 0.3)
  text = element_text(size = 24, family = "Helvetica"),
  axis.text = element_text(size = 20, family = "Helvetica"),
  # panel.border = element_blank(),
  panel.border = element_rect(colour = "transparent",size=3,fill = "transparent"),
  
  axis.line = element_line(colour = "#272727", size=1),
  axis.text.x = element_text(angle = 0, hjust = 0.2)
)

p<-p +geom_bar(stat = "identity", position=position_dodge())+ggtitle(ytext)
p<-p + scale_y_continuous(expand = c(0, 0), limits = c(miny, maxy), breaks = c(0, mody,2*mody,3*mody ))

p <- p + labs(
  # title = "Microbial profiles downstream the colon",
  x="Colon Site",
  # y="Microbial level",
  y="g/L",
    size = 19, color = "black",family = "Helvetica")
# p<-p + theme(plot.title = element_text(size = 11, color = "black", 
#                                     hjust = 0.5, face = "bold", 
#                                     angle = 30))

p<-p+scale_fill_manual(values = col_8)+guides(fill=guide_legend(ncol=2,byrow=FALSE))


# if(legend_label==1)
  p<- p + theme(legend.position="bottom",legend.text=element_text(size=22),legend.title=element_blank(),
                plot.title = element_text(hjust = 0.5,size=22,colour = "#dddddd"))
# p<- p + theme(legend.position="none") 

return(p)
}