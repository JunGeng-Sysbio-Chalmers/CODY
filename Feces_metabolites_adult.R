Feces_metabolites_adult <- function(data,Regime,legend_label) {
  # data<-pre_Data
  # Regime<-c("Baseline","Intervention")
  # data<-rbind(data_1[BFIndex,9:SF_ncol]*3,data_1[SFIndex,9:SF_ncol]*3)
  
  data$Regime<-Regime
  
colnames(data)[ncol(data)]<-"Regime"
colnames(data)[1:ncol(data)-1]<-str_replace(colnames(data)[1:ncol(data)-1],'HMO','MACs')

# data$Regime <- factor(data$Regime, levels=unique(data$Regime))
scfa_list<-c("ACETATE","PROPIONATE","BUTYRATE","SUCCINATE","LACTATE","REGIME")
colnames(data)<-toupper(colnames(data))
data_list<-data[,scfa_list]
colnames(data_list)[ncol(data_list)]<-"Regime"
data_list<-gather(data_list,Metabolites,Prediction,ACETATE:LACTATE,factor_key=TRUE)

color_mets<-c(brewer.pal(n=8,"Dark2"),brewer.pal(n=10,"Set3"))


# data_list<-data_list%>%filter(!Prediction<=0)
# data_list$Metabolites<-paste0(data_list$Metabolites,sep="_",data_list$Regime)

data_list$Regime<-factor(data_list$Regime,levels = unique(data_list$Regime))
# data_list<-data_list[ordered(data_list$Regime),]

# data_list<-inner_join(data_list,scfa_infant_exp,by="Metabolites") 


maxy<-max(data_list$Prediction)*1.75
miny<-(maxy/50)*-1
mody<-(maxy-maxy%%3)/3
# data<-melt(data)




# col_8<-col_8[c(5,1,4,2,3)]

# dat <- data.frame(
#   position_colon = factor(c("Acending","Transverse","Descending","Sigmoid"), 
#                           levels=c("Acending","Transverse","Descending","Sigmoid")),
#   total_bill = c(14.89, 17.23)
# )
# levels(position_colon)=c('ascending','transverse','descending',"sigmoid")

# qplot(name,score,data = a, geom = 'bar')
colnames(data_list)[2]<-"Metabolites"
p <-ggplot(data_list, aes(x=Regime,y=Prediction,fill=Metabolites))
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
  panel.border = element_blank(),
  axis.line = element_line(colour = "#272727", size=1)
)

p<-p +geom_bar(stat = "identity", position=position_dodge())
p<-p + scale_y_continuous(expand = c(0, 0), limits = c(miny, maxy), breaks = c(0, mody,2*mody,3*mody ))

p <- p + labs(
  # title = "Microbial profiles downstream the colon",
  x="Diet Regime",y="Metabolites level[mM]",
  size = 19, color = "black",family = "Helvetica")
# p<-p + theme(plot.title = element_text(size = 11, color = "black", 
#                                     hjust = 0.5, face = "bold", 
#                                     angle = 30))

p<-p+scale_fill_manual(values = color_mets)

# if(legend_label==1){
  p<- p + theme(legend.position="bottom", legend.title = element_blank(),
                legend.text = element_text(size = 20, color = "black",family = "Helvetica")) +
    guides(fill=guide_legend(nrow=2,byrow=TRUE))
# }
# else{p<- p + theme(legend.position="none") 
# }
p
return(p)

}

rsqred<-function (x,y) {
  y_hat <- y
  y <- x
  y_bar = mean(y)
  RSS = sum((y-y_hat)^2)
  TSS = sum((y-y_bar)^2)
  R2 = 1 - RSS/TSS
  return(R2)
}
# ################################
# # step 2: load related data
# ################################
# ## 12month
# 
# metabolite<-read.delim(file = "~/Desktop/Projects/Dynamic_Infant_Microbiome/Data_Visulization/20181217/Figure3/d/barplot/Fece_metabolites", sep = "\t", header = T, stringsAsFactors = F)
# metabolite$Time<-str_replace(metabolite$Time,"m4","(4m)")
# metabolite$Time<-str_replace(metabolite$Time,"m12","(12m)")
# metabolite$Metabolites<-paste0(metabolite$Metabolites,sep=" ",metabolite$Time)
# # metabolite$Time<-factor(metabolite$Time)
# rsqred(metabolite$Experiment,metabolite$Prediction)   ## 
# 
# metabolite<-gather(metabolite,data,value,Prediction:Experiment,factor_key=TRUE)
# metabolite$Metabolites<-factor(metabolite$Metabolites)
# # metabolite$Metabolites<-factor(metabolite$Metabolites,levels = c("Acetate (4m)","Propionate (4m)","Acetate (12m)",
# #                                "Butyrate (12m)","Propionate (12m)"))
# metabolite$Metabolites<-factor(metabolite$Metabolites,levels = c("Act (4m)","Prop (4m)","Act (12m)",
#                                                                  "Buty (12m)","Prop (12m)"))
# 
# 
# 
# p <-ggplot(metabolite, aes(x=Metabolites,y=value,fill=data,width=.5))
# p<-p + theme_minimal()
# p<-p +geom_bar(stat = "identity", position=position_dodge())
# p<-p + scale_y_continuous(expand = c(0, 0), limits = c(0, 90))
# p<-p+theme(plot.background = NULL)
# p <- p +
#   labs(x="",y="SCFA (mmol/kg feces)",size = 28, color = "black",family = "Helvetica")
# p
# p<-p + theme(
#   #           plot.title = element_text(size = 11, color = "black", 
#   #                                        hjust = 0.5, face = "bold", 
#   #                                        angle = 0),
#   text = element_text(size = 24, colour="black",family = "Helvetica"),
#   axis.text = element_text(size = 24, colour="black", family = "Helvetica"),
#   # panel.border = element_blank(),
#   panel.border = element_rect(colour = "black",size=1.25,fill = "transparent"),
#   axis.line=element_line(size=0.5,colour="black"),
#   panel.grid=element_blank(),
#   axis.text.x = element_text(size = 24, colour="black", family = "Helvetica",vjust=-1.5),
#   axis.title.y = element_text(vjust=5),
#   
#   axis.ticks.length=unit(0.42, "line"),
#   axis.ticks = element_line(colour = "black", size = 0.7)
#   # axis.title.x = element_text(vjust=-10)
# )
# # p<-p+theme(panel.grid=element_blank(),panel.border=element_blank())
# 
# # p<-p+scale_fill_brewer(palette="Dark2")
# # p<-p+scale_fill_brewer(palette="Accent")
# p<-p+scale_fill_manual(values = c("lightsalmon2", "#5d8e92","lightpink4"))
# p<- p + theme(legend.position="top",legend.text=element_text(size=26)) 
# p
# 
# 
# }