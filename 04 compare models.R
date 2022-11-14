#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Lattice data CPUE standardization NOAA YFT  #
# Code to compare the three INLA models       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Francisco Izquierdo #
#~~~~~~~~~~~~~~~~~~~~~~
# 25/10/2022 #        
#~~~~~~~~~~~~~

## press Ctrl + shift + O to see the script outline

# model iid --------------------------------------------------------------------

rm(list=ls()) 

## create dir
out_04<-paste(getwd(),"/output/04 compare models", sep="")
dir.create(out_04)

load("./output/01 iid model/01 iid model.RData")
mod1<-model
idx_iid_1<-data_idx_1
idx_iid_4<-data_idx_4

# model besag ------------------------------------------------------------------

load("./output/02 besag model/02 besag model.RData")
mod2<-model
idx_besag_1<-data_idx_1
idx_besag_4<-data_idx_4

# model besag group ------------------------------------------------------------

load("./output/03 besag model group/03 besag model group.RData")
mod3<-model
idx_besag_st_1<-data_idx_1
idx_besag_st_4<-data_idx_4

# compare fit ------------------------------------------------------------------

COMP<-data.frame(
   model=c("model iid","model besag","model besag st"),
   dic=c(mod1$dic$dic,mod2$dic$dic,mod3$dic$dic),
   waic=c(mod1$waic$waic,mod2$waic$waic,mod3$waic$waic),
   lcpo=c(-mean(log(mod1$cpo$cpo),na.rm=T),-mean(log(mod2$cpo$cpo),na.rm=T),-mean(log(mod3$cpo$cpo),na.rm=T)),
   failure=c(sum((mod1$cpo$failure>0)*1,na.rm=T),sum((mod2$cpo$failure>0)*1,na.rm=T),sum((mod3$cpo$failure>0)*1,na.rm=T)),
   time=c(mod1$cpu.used[4],mod2$cpu.used[4],mod3$cpu.used[4]))

COMP

write.csv(COMP, paste0(out_04,"/04 compare all models.csv"))

# compare idx ------------------------------------------------------------------

# 1 area -----------------------------------------------------------------------

library(ggplot2)

## create model var
idx_iid_1$Model<-rep("iid")
idx_besag_1$Model<-rep("besag")
idx_besag_st_1$Model<-rep("besag st")

## raw observed data
library(plyr)
idx_obs<-ddply(data, .(pseudoyear), summarize, "meanX50"=mean(mcpue))
idx_obs$meanX5<-rep(NA)
idx_obs$meanX95<-rep(NA)
idx_obs$Model<-rep("observed")

## join data frames
idx_data_1<-rbind(idx_iid_1, idx_besag_1, idx_besag_st_1)

ggplot(data = idx_data_1, aes(x = pseudoyear, y = meanX50, col=Model))+
   #geom_ribbon(aes(x = pseudoyear, ymin = meanX5, ymax = meanX95), 
   #            alpha = 0.3,fill="grey", size=0.3)+
   geom_line(size=0.9)+
   scale_color_manual(values = c("orange","dodgerblue3","black"))+
   #geom_line(data=idx_obs, aes(x=as.numeric(pseudoyear), y=meanX50, group=1))+
   ylab("Predicted median CPUE")+
   xlab("Pseudoyear")+
   theme_light() + scale_x_continuous(breaks= seq(from=1, to=length(unique(data$psyearID)), by=12),
                                      labels=seq(from=81, to=256, by=12),
                                      limits=c(1, length(unique(data$pseudoyear))))
## save
ggsave(paste0(out_04,"/04 all model predicted CPUE indices_1area.png"), dpi=300, 
       width=6, height=4)

write.csv(idx_data_1, paste0(out_04,"/04 CPUE indices models 1 area.csv"))

## observed data
ggplot(data=idx_obs, aes(x=as.numeric(pseudoyear), y=meanX50, group=1))+
   geom_line()+theme_light()
ggsave(paste0(out_04,"/04 observed CPUE_1area obs.png"), dpi=300, 
       width=6, height=4)

# 4 areas ----------------------------------------------------------------------

## create model var
idx_iid_4$Model<-rep("iid")
idx_besag_4$Model<-rep("besag")
idx_besag_st_4$Model<-rep("besag st")

## raw observed data
idx_obs<-ddply(data, .(pseudoyear,regionID), summarize, "meanX50"=mean(mcpue))
idx_obs$meanX5<-rep(NA)
idx_obs$meanX95<-rep(NA)
idx_obs$Model<-rep("observed")

## join data frames
idx_data_4<-rbind(idx_iid_4, idx_besag_4, idx_besag_st_4)

ggplot(data = idx_data_4, aes(x = pseudoyear, y = meanX50, col=Model))+
   #geom_ribbon(aes(x = pseudoyear, ymin = meanX5, ymax = meanX95), 
   #            alpha = 0.3,fill="grey", size=0.3)+
   geom_line(size=0.9)+
   scale_color_manual(values = c("orange","dodgerblue3","black","grey"))+
   #geom_line(data=idx_obs, aes(x=as.numeric(pseudoyear), y=meanX50, group=1))+
   ylab("Predicted median CPUE")+
   xlab("Pseudoyear")+
   facet_wrap(~regionID)+
   theme_light() + scale_x_continuous(breaks= seq(from=1, to=length(unique(data$psyearID)), by=20),
                                      labels=seq(from=81, to=256, by=20),
                                      limits=c(1, length(unique(data$pseudoyear))))
## save
ggsave(paste0(out_04,"/04 all model predicted CPUE indices_4areas.png"), dpi=300, 
       width=6, height=4)

write.csv(idx_data_4, paste0(out_04,"/04 CPUE indices models 4 areas.csv"))

## observed data
ggplot(data=idx_obs, aes(x=as.numeric(pseudoyear), y=meanX50, group=1))+
   geom_line()+theme_light() + facet_wrap(~regionID)

ggsave(paste0(out_04,"/04 observed CPUE_4areas obs.png"), dpi=300, 
       width=6, height=4)

# prepare SS input year --------------------------------------------------------

## 1 area, besag st
## year seas index obs se
library(dplyr)
CPUE_bst1<-idx_data_1%>%select(year, seas, meanX50, meanSd, Model)%>%filter(Model=="besag st")
CPUE_bst1$index<-rep(8) # cpue fleet
levels(CPUE_bst1$seas)<-c(1,4,7,10) ## season levels
CPUE_bst1$se<-(CPUE_bst1$meanSd/sqrt(100)) # standard error
#CPUE_bst1$log_se<-(log(CPUE_bst1$se)) # log standard error
CPUE_bst1<-CPUE_bst1[,c(1,2,6,3,7)] # order
colnames(CPUE_bst1)<-c("year","seas","index","obs","se") # rename cols
CPUE_bst1$year<-as.numeric(as.character(CPUE_bst1$year))
CPUE_bst1$seas<-as.numeric(as.character(CPUE_bst1$seas))
CPUE_bst1<-round(CPUE_bst1,4)
write.table(CPUE_bst1, paste0(out_04,"/04 CPUE_bst1_year.txt"), row.names = F) # save txt

## 4 area, besag st
## year seas index obs se
library(dplyr)
CPUE_bst4<-idx_data_4%>%select(year, seas,regionID, meanX50, meanSd, Model)%>%filter(Model=="besag st")
levels(CPUE_bst4$regionID) # rename regions 1 2 3 4 with fleet 17 18 19 20
levels(CPUE_bst4$regionID)<-c(17,18,19,20) ## region ID
levels(CPUE_bst4$seas)<-c(1,4,7,10) ## season levels
CPUE_bst4$se<-(CPUE_bst4$meanSd/sqrt(100)) # standard error
#CPUE_bst4$log_se<-(log(CPUE_bst4$se)) # log standard error
CPUE_bst4<-CPUE_bst4[,c(1,2,3,4,7)] # order
colnames(CPUE_bst4)<-c("year","seas","index","obs","se") # rename cols
CPUE_bst4$year<-as.numeric(as.character(CPUE_bst4$year))
CPUE_bst4$seas<-as.numeric(as.character(CPUE_bst4$seas))
CPUE_bst4$index<-as.numeric(as.character(CPUE_bst4$index))
CPUE_bst4<-round(CPUE_bst4,4)
write.table(CPUE_bst4, paste0(out_04,"/04 CPUE_bst4_year.txt"), row.names = F) # save 


## 1 area, iid
## year seas index obs se
library(dplyr)
CPUE_bst1<-idx_data_1%>%select(year, seas, meanX50, meanSd, Model)%>%filter(Model=="iid")
CPUE_bst1$index<-rep(8) # cpue fleet
levels(CPUE_bst1$seas)<-c(1,4,7,10) ## season levels
CPUE_bst1$se<-(CPUE_bst1$meanSd/sqrt(100)) # standard error
#CPUE_bst1$log_se<-(log(CPUE_bst1$se)) # log standard error
CPUE_bst1<-CPUE_bst1[,c(1,2,6,3,7)] # order
colnames(CPUE_bst1)<-c("year","seas","index","obs","se") # rename cols
CPUE_bst1$year<-as.numeric(as.character(CPUE_bst1$year))
CPUE_bst1$seas<-as.numeric(as.character(CPUE_bst1$seas))
CPUE_bst1<-round(CPUE_bst1,4)
write.table(CPUE_bst1, paste0(out_04,"/04 CPUE_iid1_year.txt"), row.names = F) # save txt

## 4 area, iid
## year seas index obs se
library(dplyr)
CPUE_bst4<-idx_data_4%>%select(year, seas,regionID, meanX50, meanSd, Model)%>%filter(Model=="iid")
levels(CPUE_bst4$regionID) # rename regions 1 2 3 4 with fleet 17 18 19 20
levels(CPUE_bst4$regionID)<-c(17,18,19,20) ## region ID
levels(CPUE_bst4$seas)<-c(1,4,7,10) ## season levels
CPUE_bst4$se<-(CPUE_bst4$meanSd/sqrt(100)) # standard error
#CPUE_bst4$log_se<-(log(CPUE_bst4$se)) # log standard error
CPUE_bst4<-CPUE_bst4[,c(1,2,3,4,7)] # order
colnames(CPUE_bst4)<-c("year","seas","index","obs","se") # rename cols
CPUE_bst4$year<-as.numeric(as.character(CPUE_bst4$year))
CPUE_bst4$seas<-as.numeric(as.character(CPUE_bst4$seas))
CPUE_bst4$index<-as.numeric(as.character(CPUE_bst4$index))
CPUE_bst4<-round(CPUE_bst4,4)
write.table(CPUE_bst4, paste0(out_04,"/04 CPUE_iid4_year.txt"), row.names = F) # save 


# prepare SS input psyear ------------------------------------------------------

# psyear hay que sumarle 1000 para el modelo de Giancarlo
# season siempre es 7

## 1 area, besag st
## year seas index obs se
library(dplyr)
CPUE_bst1<-idx_data_1%>%select(pseudoyear, seas, meanX50, meanSd, Model)%>%filter(Model=="besag st")
CPUE_bst1$index<-rep(8) # cpue fleet
CPUE_bst1$seas<-rep(7) # seas 7
CPUE_bst1$se<-(CPUE_bst1$meanSd/sqrt(100)) # standard error
CPUE_bst1$pseudoyear<-(1080+CPUE_bst1$pseudoyear)
CPUE_bst1<-CPUE_bst1[,c(1,2,6,3,7)] # order
colnames(CPUE_bst1)<-c("year","seas","index","obs","se") # rename cols
CPUE_bst1$year<-as.numeric(as.character(CPUE_bst1$year))
CPUE_bst1$seas<-as.numeric(as.character(CPUE_bst1$seas))
CPUE_bst1<-round(CPUE_bst1,4)
write.table(CPUE_bst1, paste0(out_04,"/04 CPUE_bst1_psyear.txt"), row.names = F) # save txt

## 4 area, besag st
## year seas index obs se
library(dplyr)
CPUE_bst4<-idx_data_4%>%select(pseudoyear, seas,regionID, meanX50, meanSd, Model)%>%filter(Model=="besag st")
CPUE_bst4$index<-rep(8) # cpue fleet
CPUE_bst4$seas<-rep(7) # seas 7
CPUE_bst4$se<-(CPUE_bst4$meanSd/sqrt(100)) # standard error
CPUE_bst4$pseudoyear<-(1080+CPUE_bst4$pseudoyear)
CPUE_bst4<-CPUE_bst4[,c(1,2,3,4,7)] # order
colnames(CPUE_bst4)<-c("year","seas","index","obs","se") # rename cols
CPUE_bst4$year<-as.numeric(as.character(CPUE_bst4$year))
CPUE_bst4$seas<-as.numeric(as.character(CPUE_bst4$seas))
CPUE_bst4$index<-as.numeric(as.character(CPUE_bst4$index))
CPUE_bst4<-round(CPUE_bst4,4)
write.table(CPUE_bst4, paste0(out_04,"/04 CPUE_bst4_psyear.txt"), row.names = F) # save 


## 1 area, iid
## year seas index obs se
library(dplyr)
CPUE_bst1<-idx_data_1%>%select(pseudoyear, seas, meanX50, meanSd, Model)%>%filter(Model=="iid")
CPUE_bst1$index<-rep(8) # cpue fleet
CPUE_bst1$seas<-rep(7) # seas 7
CPUE_bst1$se<-(CPUE_bst1$meanSd/sqrt(100)) # standard error
CPUE_bst1$pseudoyear<-(1080+CPUE_bst1$pseudoyear)
CPUE_bst1<-CPUE_bst1[,c(1,2,6,3,7)] # order
colnames(CPUE_bst1)<-c("year","seas","index","obs","se") # rename cols
CPUE_bst1$year<-as.numeric(as.character(CPUE_bst1$year))
CPUE_bst1$seas<-as.numeric(as.character(CPUE_bst1$seas))
CPUE_bst1<-round(CPUE_bst1,4)
write.table(CPUE_bst1, paste0(out_04,"/04 CPUE_iid1_psyear.txt"), row.names = F) # save txt

## 4 area, iid
## year seas index obs se
library(dplyr)
CPUE_bst4<-idx_data_4%>%select(pseudoyear, seas,regionID, meanX50, meanSd, Model)%>%filter(Model=="iid")
CPUE_bst4$index<-rep(8) # cpue fleet
CPUE_bst4$seas<-rep(7) # seas 7
CPUE_bst4$se<-(CPUE_bst4$meanSd/sqrt(100)) # standard error
CPUE_bst4$pseudoyear<-(1080+CPUE_bst4$pseudoyear)
CPUE_bst4<-CPUE_bst4[,c(1,2,3,4,7)] # order
colnames(CPUE_bst4)<-c("year","seas","index","obs","se") # rename cols
CPUE_bst4$year<-as.numeric(as.character(CPUE_bst4$year))
CPUE_bst4$seas<-as.numeric(as.character(CPUE_bst4$seas))
CPUE_bst4$index<-as.numeric(as.character(CPUE_bst4$index))
CPUE_bst4<-round(CPUE_bst4,4)
write.table(CPUE_bst4, paste0(out_04,"/04 CPUE_iid4_psyear.txt"), row.names = F) # save 
