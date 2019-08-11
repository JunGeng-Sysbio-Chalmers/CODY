Plasma_metabolites_barplot <- function(data,Regime,legend_label) {
  # data<-data_blood_metabolite_all
  # data<-rbind(data_1[BFIndex,2:SF_ncol],data_1[SFIndex,2:SF_ncol])
  # Regime<-c("Breastfeeding","Solidfood")

data$Regime<-Regime
  # data$Regime<-str_replace(data$Regime,"Breastfeeding","BF")
  # data$Regime<-str_replace(data$Regime,"Solidfood","SF")
  scfa_list<-c("ACETATE","PROPIONATE","BUTYRATE","SUCCINATE","LACTATE","FORMATE","REGIME")
  colnames(data)<-toupper(colnames(data))
  data_list<-data[,scfa_list]
  colnames(data_list)[ncol(data_list)]<-"Regime"
  colnames(data)[2:ncol(data)]<-str_replace(colnames(data)[2:ncol(data)],'HMO','MACs')
  
  data_list<-gather(data_list,Metabolites,Prediction,ACETATE:FORMATE,factor_key=TRUE)
  
  
  # data_list<-data_list%>%filter(!Prediction<=0)
  # data_list$Metabolites<-paste0(data_list$Metabolites,sep="_",data_list$Regime)
  
  data_list$Regime<-factor(data_list$Regime,levels = unique(data_list$Regime))
  # data_list<-data_list[ordered(data_list$Regime),]

  # data_list<-inner_join(data_list,scfa_infant_exp,by="Metabolites") 
  

  maxy<-max(data_list$Prediction)*1.75
  miny<-(maxy/50)*-1
  mody<-(maxy-maxy%%3)/3
  
  p <-ggplot(data_list, aes(x=Metabolites,y=Prediction,fill=Regime,width=.5))
  p<-p + theme_minimal()
  p<-p +geom_bar(stat = "identity", position=position_dodge())
  p<-p + scale_y_continuous(expand = c(0, 0), limits = c(miny, maxy))+
      # scale_x_discrete(labels=c("ACETATE" = "AC", "PROPIONATE"="PROP","BUTYRATE"="BUTY","SUCCINATE"="SUCC","LACTATE"="LAC","FORMATE"="FOR")
       theme(plot.background = NULL)
  p <- p +
    labs(x="",y="Absolute level [mM]",size = 20, color = "black",family = "Helvetica")+ggtitle("Plasma")
  p
  p<-p + theme(
    #           plot.title = element_text(size = 11, color = "black", 
    #                                        hjust = 0.5, face = "bold", 
    #                                        angle = 0),
    text = element_text(size = 23, colour="black",family = "Helvetica"),
    axis.text = element_text(size = 23, colour="black", family = "Helvetica"),
    # panel.border = element_blank(),
    panel.border = element_rect(colour = NA,fill = "transparent",size=2.5),
    axis.line=element_line(size=1,colour="black"),
    panel.grid=element_blank(),
    axis.text.x = element_text(size = 19, colour="black", family = "Helvetica", vjust = 0.5, hjust = 0.25,angle=330),
    axis.title.y = element_text(size = 19, colour="black", family = "Helvetica"),
    
    axis.ticks.length=unit(0.4, "line"),
    axis.ticks = element_line(colour = "black", size = 0.7)
    # axis.title.x = element_text(vjust=-10)
  )
  p
  # p<-p+theme(panel.grid=element_blank(),panel.border=element_blank())
  
  # p<-p+scale_fill_brewer(palette="Dark2")
  # p<-p+scale_fill_brewer(palette="Accent")
  # p<-p+scale_fill_manual(values = c("#954A45","#986DB2","#AF5F3C","#C18A26"))
  p<-p+scale_fill_manual(values = c("#954A45","#6F3381"))
  
  #
  # p<-p+scale_fill_manual(values = c("#f384b4", "#32c5ca"))
  # p<- p + theme(legend.position="top",legend.text=element_text(size=26)) 
  p
  # if(legend_label==1){
    p<- p + theme(legend.position="bottom",legend.text=element_text(size=23),legend.title=element_blank(),
                                    # plot.margin = margin(2,.8,2,.8, "cm"),
                  plot.title = element_text(hjust = 0.5,size=21,colour = "black",family = "Helvetica")) 
  
  # }
  # else{p<- p + theme(legend.position="none")
  # }
    p
  return(p)
  # ggsave("SCFA_profile_feces_0419.pdf",device="pdf",width=11,height=5.1)
  
}











