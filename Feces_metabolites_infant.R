Feces_metabolites_infant <- function(data,Regime,legend_label) {
# data<-data_feces_metabolite_all
# Regime<-c("Breastfeeding","Solidfood")
  # data<-data_1
data$Regime<-Regime
color_mets<-c(brewer.pal(n=8,"Dark2"),brewer.pal(n=10,"Set3"))

# data<-data_feces_metabolite_all
# colnames(data)[ncol(data)]<-"Regime"
# 
# data$Regime <- factor(data$Regime, levels=unique(data$Regime))
colnames(data)[1:ncol(data)-1]<-str_replace(colnames(data)[1:ncol(data)-1],'HMO','MACs')
# 
# colnames(data)[2]<-"Metabolites"
# p <-ggplot(data, aes(x=Regime,y=value,fill=Metabolites))
# # p<-p + theme_minimal()
# p<-p+theme(panel.grid=element_blank(),panel.border=element_blank(),axis.line=element_line(size=0.1,colour="#5e5e5e"),
#            # axis.text.x = element_text(angle = 13, hjust = 1))
#            axis.text.x = element_text(angle = 0, hjust = 0.2))
# 
# p<-p+theme(
#   panel.grid.major = element_blank(),
#   panel.grid.minor = element_blank(),
#   # panel.border = element_rect(
#   #   colour = "#5e5e5e",
#   #   fill = NA,
#   #   size = 0.3)
#   text = element_text(size = 24, family = "Helvetica"),
#   axis.text = element_text(size = 20, family = "Helvetica"),
#   panel.border = element_blank(),
#   axis.line = element_line(colour = "#272727", size=1)
# )
# 
# p<-p +geom_bar(stat = "identity", position=position_dodge())
# p<-p + scale_y_continuous(expand = c(0, 0), limits = c(miny, maxy), breaks = c(0, mody,2*mody,3*mody ))
# 
# p <- p + labs(
#   # title = "Microbial profiles downstream the colon",
#   x="Diet Regime",y="Absolute Metabolites level",
#   size = 19, color = "black",family = "Helvetica")
# # p<-p + theme(plot.title = element_text(size = 11, color = "black", 
# #                                     hjust = 0.5, face = "bold", 
# #                                     angle = 30))
# 
# p<-p+scale_fill_manual(values = col_8)
# 
# if(legend_label==1){p<- p + theme(legend.position="top") 
# }
# else{p<- p + theme(legend.position="none") 
# }
# 
# return(p)



rsqred<-function (x,y) {
  y_hat <- y
  y <- x
  y_bar = mean(y)
  RSS = sum((y-y_hat)^2)
  TSS = sum((y-y_bar)^2)
  R2 = 1 - RSS/TSS
  return(R2)
}

# data<-data_feces_metabolite_all

data$Regime<-str_replace(data$Regime,"Breastfeeding","BF")
data$Regime<-str_replace(data$Regime,"Solidfood","SF")
scfa_list<-c("ACETATE","PROPIONATE","BUTYRATE","Regime")
data_list<-data[,scfa_list]
data_list<-gather(data_list,Metabolites,Prediction,ACETATE:BUTYRATE,factor_key=TRUE)
           

data_list<-data_list%>%filter(!Prediction<=0)
data_list$Metabolites<-paste0(data_list$Metabolites,sep="_",data_list$Regime)

data_list$Metabolites<-factor(data_list$Metabolites,levels = c("ACETATE_BF","PROPIONATE_BF","ACETATE_SF",
                                                                 "BUTYRATE_SF","PROPIONATE_SF"))
data_list<-data_list[ordered(data_list$Metabolites),]
scfa_infant_exp<-data.frame(Metabolites=c("ACETATE_BF","PROPIONATE_BF","ACETATE_SF","BUTYRATE_SF","PROPIONATE_SF"),
                            Experiment=c(53.9,2.9,60,16,19))
scfa_infant_exp$Metabolites<-as.factor(scfa_infant_exp$Metabolites)

data_list<-inner_join(data_list,scfa_infant_exp,by="Metabolites") 


R2_feces_mets_infant<-rsqred(data_list$Experiment,data_list$Prediction)   ## 

data_list<-gather(data_list,data,value,Prediction:Experiment,factor_key=TRUE)
data_list$Metabolites<-factor(data_list$Metabolites)
# metabolite$Metabolites<-factor(metabolite$Metabolites,levels = c("Acetate (4m)","Propionate (4m)","Acetate (12m)",
#                                "Butyrate (12m)","Propionate (12m)"))
data_list$Metabolites<-factor(data_list$Metabolites,levels = c("ACETATE_BF","PROPIONATE_BF","ACETATE_SF","BUTYRATE_SF","PROPIONATE_SF"))
maxy<-max(data_list$value)*1.75
miny<-(maxy/50)*-1
mody<-(maxy-maxy%%3)/3

p <-ggplot(data_list, aes(x=Metabolites,y=value,fill=data,width=.5))
p<-p + theme_minimal()
p<-p +geom_bar(stat = "identity", position=position_dodge())
p<-p + scale_y_continuous(expand = c(0, 0), limits = c(0, 90))
p<-p+theme(plot.background = NULL)
p <- p +
  labs(x="",y="SCFA (mmol/kg feces)",size = 28, color = "black",family = "Helvetica")
p
p<-p + theme(
  #           plot.title = element_text(size = 11, color = "black", 
  #                                        hjust = 0.5, face = "bold", 
  #                                        angle = 0),
  text = element_text(size = 24, colour="black",family = "Helvetica"),
  axis.text = element_text(size = 24, colour="black", family = "Helvetica"),
  # panel.border = element_blank(),
  panel.border = element_rect(colour = "transparent",size=1.5,fill = "transparent"),
  axis.line=element_line(size=0.5,colour="black"),
  panel.grid=element_blank(),
  # axis.text.x = element_text(size = 24, colour="black", family = "Helvetica",vjust=-1.5, hjust = 0.5),
  # axis.title.y = element_text(vjust=5),
  axis.text.x = element_text(size = 19, colour="black", family = "Helvetica", vjust = 0.5, hjust = 0.25,angle=330),
  axis.title.y = element_text(size = 19, colour="black", family = "Helvetica"),
  
  axis.ticks.length=unit(0.42, "line"),
  axis.ticks = element_line(colour = "black", size = 0.7)
  # axis.title.x = element_text(vjust=-10)
)
p<-p+theme(panel.grid=element_blank())

# p<-p+scale_fill_manual(values=color_mets)
p<-p+scale_fill_manual(values = c("#954A45","#6F3381"))

# p<-p+scale_fill_brewer(palette="Accent")
# p<-p+scale_fill_manual(values = c("lightsalmon2", "#5d8e92","lightpink4"))
# p<- p + theme(legend.position="top",legend.text=element_text(size=26)) 

# if(legend_label==1){
  p<- p + theme(legend.position="bottom",legend.text=element_text(size=22,color = "black", family = "Helvetica"),
                legend.title=element_blank(),
                plot.margin = margin(2,.8,2,.8, "cm"))

# }
# else{p<- p + theme(legend.position="none")
# }

return(p)
# ggsave("SCFA_profile_feces_0419.pdf",device="pdf",width=11,height=5.1)

}