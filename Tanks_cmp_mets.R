Tanks_cmp_mets <- function(data,header,ytext) {
# data<-data_lumen_metabolite_all_SF
# bacteria_12m_1102 <- read_excel("~/Desktop/Projects/Dynamic_Infant_Microbiome/Data_Visulization/20181027/Figure4/d/bacteria_12m_1102.xlsx")
# bacteria<-bacteria_12m_1102[1:4,]
# colnames(bacteria)[1]<-"Position"
# mets_12m<-read_excel("~/Desktop/Projects/Dynamic_Infant_Microbiome/Data_Visulization/20181027/Figure4/d/mets_12m_1102.xlsx")
# mets<-mets_12m[1:4,c(1,6,10,5)]   ## prop, but, ac
# colnames(mets)[1]<-"Position"\
  
  # data<-rbind(data_1[SFIndex,10:SF_ncol],data_2[SFIndex,10:SF_ncol],data_3[SFIndex,10:SF_ncol],data_4[SFIndex,10:SF_ncol])
  
  
Colon_Site<-c("I","II","III","IV")

data$Position<-Colon_Site
# ## 4month
# bacteria_4m_1102 <- read_excel("~/Desktop/Projects/Dynamic_Infant_Microbiome/Data_Visulization/20181027/Figure4/d/bacteria_4m_1102.xlsx")
# bacteria<-bacteria_4m_1102[1:4,]
# colnames(bacteria)[1]<-"Position"
# mets_4m<-read_excel("~/Desktop/Projects/Dynamic_Infant_Microbiome/Data_Visulization/20181027/Figure4/d/mets_4m_1102.xlsx")
# mets<-mets_4m[1:4,c(1,6,5)]
# colnames(mets)[1]<-"Position"


data$Position <- factor(data$Position, levels=unique(data$Position))
colnames(data)[1:ncol(data)]<-str_replace(colnames(data)[1:ncol(data)],'_','.')
colnames(data)[1:ncol(data)]<-str_replace(colnames(data)[1:ncol(data)],'HMO','MACs')
data<-data[,-which(colnames(data) %in% c("FORMATE","H2"))]
# bacteria<-bacteria[c(1,7,5,6,4,2,9,8,3)]
# colnames(bacteria)<-c("Position","Ehallii" ,"Bbr", "Bad", "Blg", "Bth","Rint", "Fpr","Bfr"   )
# bacteria<-bacteria[c(1,6,2,5,3,4)]
# colnames(bacteria)<-c("Position","Bad" ,"Bth", "Bbr", "Bfr", "Blg"  )
# data$Position <- factor(data$Position, levels=unique(data$Position))
maxy<-max(apply(data[,1:ncol(data)-1],2,max))*1.5
miny<-(maxy/50)*-1
mody<-(maxy-maxy%%3)/3
mets<-melt(data)

# col_8<-col_8[c(5,1,4,2,3)]
color_mets<-c(brewer.pal(n=8,"Dark2"),brewer.pal(n=10,"Set3"))

# dat <- data.frame(
#   position_colon = factor(c("Acending","Transverse","Descending","Sigmoid"), 
#                           levels=c("Acending","Transverse","Descending","Sigmoid")),
#   total_bill = c(14.89, 17.23)
# )
# levels(position_colon)=c('ascending','transverse','descending',"sigmoid")

# qplot(name,score,data = a, geom = 'bar')
colnames(mets)[2]<-"Metabolites"
mets$Metabolites<-factor(mets$Metabolites,levels = unique(mets$Metabolites))
p <-ggplot(mets, aes(x=Position,y=value,fill=Metabolites))
p<-p + theme_minimal()+theme(panel.border = element_rect(colour = "transparent",size=2,fill = "transparent"))
p<-p +geom_bar(stat = "identity", position=position_dodge())+ggtitle(ytext)
p<-p + scale_y_continuous(expand = c(0, 0), limits = c(miny, maxy))
p<-p+theme(plot.background = NULL)
p<-p+theme(panel.grid=element_blank(),axis.line=element_line(size=1,colour="black"))

p <- p + labs(
  # title = "Metabolite profiles downstream the colon",
              x="Colon Site",
              # y="In vivo Absolute Metabolite level[mmol/L]",
              y="[mM]",
              
              size = 13, color = "black",family = "Helvetica")
p<-p + theme(plot.title = element_text(size = 22, colour = "#6e6e6e", 
                                       hjust = 0.5, face = "bold", 
                                       angle = 0),
             text = element_text(size = 20, family = "Helvetica"),
             axis.text = element_text(size = 20, family = "Helvetica"),
             # panel.border = element_blank(),
             axis.text.x = element_text(angle = 0, hjust = 0.2))

p<-p+scale_fill_manual(values=color_mets)
# p<-p+scale_fill_brewer(palette="Accent")
# p<-p+scale_fill_manual(values = c("lightsalmon2", "darkcyan","lightpink4"))
# p<- p + theme(legend.position="top") 
# if(legend_label==1)
  p<- p + theme(legend.position="bottom",legend.text=element_text(size=20),legend.title=element_blank()) 
  # p<- p + theme(legend.position="none") 

p
return(p)
}