#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Lattice data CPUE standardization NOAA YFT  #
# Code for INLA spatial besag model           #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Francisco Izquierdo #
#~~~~~~~~~~~~~~~~~~~~~~
# 25/10/2022 #        
#~~~~~~~~~~~~~

## press Ctrl + shift + O to see the script outline

## note that you must run the script "00 prepare data.R" first

## besag model:
## response variable = CPUE, family=gamma
## this model has a spatial structure and separated temporal trend
## spatial effect (polyID) is modelled as besag effect
## temporal effect (psyearID) is modelled as a RW2 effect

## assumptions:
## it assumes the same spatial distribution for all years
## it assumes that all time steps have the same spatial distribution
## the different polys will progress exactly the same than the neighbours along time

## prediction is = to fitted values
## but we still need posterior.samples to account for uncertainty

# start here -------------------------------------------------------------------

rm(list=ls()) 
library(INLA)
load(file="./data/datasets.RData")

## create dir
out_02<-paste(getwd(),"/output/02 besag model", sep="")
dir.create(out_02)

# poly2nb ----------------------------------------------------------------------

## for besag model we need to calculate neighbors relationship 
library("spdep")
grid.nb <- poly2nb(poly_grid) # just a single point (i.e., queen adjacency) by contiguity
xy<-unique(data[,2:3])
png(paste0(out_02,"/02 besag model queen adjacency nb.png"), width = 1000, height = 800, res = 150)
plot(st_geometry(data$geometry),main="queen adjacency", axes=TRUE)
plot(grid.nb, xy, add=T, col="blue")
dev.off()

# graph ------------------------------------------------------------------------

## create the graph 
nb2INLA(paste0(out_02,"/grid_graph.graph"), grid.nb)
grid_graph <- inla.read.graph(paste0(out_02,"/grid_graph.graph"))

## plot
png(paste0(out_02,"/02 besag model graph mat.png"), width = 1000, height = 800, res = 150)
image(inla.graph2matrix(grid_graph),xlab="",ylab="")
dev.off()

png(paste0(out_02,"/02 besag model graph fig.png"), width = 1000, height = 800, res = 150)
plot(grid_graph)
dev.off()

# model ------------------------------------------------------------------------

formula <- mcpue ~ 1 + 
   f(polyID, model = "besag", graph = grid_graph, scale.model = TRUE) +
   f(psyearID, model="rw2")

model <- inla(formula,
              family          = "gamma",
              data            = data,
              control.compute = list(dic = TRUE, config = TRUE, waic = TRUE, cpo = TRUE),
              control.predictor = list(compute=TRUE, cdf=c(log(1))),
              control.inla = list(strategy = 'adaptive'), 
              verbose=TRUE, num.threads = 1)

summary(model)

## save plots for all model effects
pdf(paste0(out_02,"/02 besag model 1A_psy output.pdf"), 
    width = 8, height = 7, 
    bg = "white",         
    colormodel = "cmyk",    
    paper = "A4")          

plot(model)
dev.off() 

# plot effects -----------------------------------------------------------------

poly_grid_out<-poly_grid
## posterior distribution of the besag effects
poly_grid_out$besagMean <- round(model$summary.random$polyID[["mean"]], 4)
poly_grid_out$besagSd <- round(model$summary.random$polyID[["sd"]],5)

## posterior fitted values
data_out<-data
data_out$CPUE_mean <- model$summary.fitted.values$mean # mean
data_out$CPUE_sd <- model$summary.fitted.values$sd #s
data_out$CPUE_median <- model$summary.fitted.values$`0.5quant` # median
data_out$CPUE_q025 <- model$summary.fitted.values$`0.025quant` # quantile
data_out$CPUE_q975 <- model$summary.fitted.values$`0.975quant` # quantile

library(sf)
library(ggspatial)
library(rnaturalearth)
library(maps)
library(ggplot2)

## world map 
world_map_data <- ne_countries(scale = "medium", returnclass = "sf")

## besag mean ------------------------------------------------------------------

library(ggplot2)
# map SPmean
ggplot()+
   geom_sf(data=poly_grid_out, aes(geometry=geometry, fill=besagMean))+
   scale_fill_viridis_c(option="H") + 
   geom_sf(data = world_map_data,alpha=0.8) +
   theme_light() +
   coord_sf(expand=FALSE, xlim = c(10, 115), ylim = c(-45,40))+
   labs(fill="Mean besag")

## save
ggsave(paste0(out_02,"/02 besag model 1A_psy mean besag effect.png"), dpi=300, width=6, height=6)

## besag sd --------------------------------------------------------------------

# map SPmean
ggplot()+
   geom_sf(data=poly_grid_out, aes(geometry=geometry, fill=besagSd))+
   scale_fill_viridis_c(option="H") + 
   geom_sf(data = world_map_data,alpha=0.8) +
   theme_light() +
   coord_sf(expand=FALSE, xlim = c(10, 115), ylim = c(-45,40))+
   labs(fill="Sd besag")

## save
ggsave(paste0(out_02,"/02 besag model 1A_psy sd besag effect.png"), dpi=300, width=6, height=6)

## besag shape -------------------------------------------------------------------

library(ggplot2)
suabinm <- model$summary.random$polyID$mean
suabin2 <- model$summary.random$polyID$`0.025quant`
suabin9 <-model$summary.random$polyID$`0.975quant`
suabinID<-model$summary.random$polyID$ID
suabin<-data.frame(suabinm, suabin2,suabin9,suabinID)

ggplot(data = suabin, aes(x = factor(suabinID), y = suabinm, group=1))+
   geom_line(aes(x = suabinID, y = suabinm), color="#29AF7FFF", size=0.9)+ 
   geom_ribbon(aes(x = suabinID, ymin = (suabin2), ymax = (suabin9)), 
               alpha = 0.25, fill="gray70", linetype=1)+
   ggtitle(" ")+
   xlab("Besag ID")+
   ylab("Besag effect ")+
   theme_light() +
   theme(plot.title = element_text(hjust=0.5))

## save
ggsave(paste0(out_02,"/02 1A_psy shape besag effect.png"), dpi=300, width=6, height=6)


## RW2 year --------------------------------------------------------------------

library(ggplot2)
suabinm <- model$summary.random$psyearID$mean
suabin2 <- model$summary.random$psyearID$`0.025quant`
suabin9 <-model$summary.random$psyearID$`0.975quant`
suabinID<-model$summary.random$psyearID$ID
suabin<-data.frame(suabinm, suabin2,suabin9,suabinID)

ggplot(data = suabin, aes(x = suabinID, y = suabinm))+
   
   geom_line(aes(x = suabinID, y = suabinm), color="#33638DFF", size=0.9)+
   geom_line(aes(suabinID, (suabin2)), color = "grey50", size = 0.1, linetype="dashed") + 
   geom_line(aes(suabinID, (suabin9)), color = "grey50", size = 0.1, linetype="dashed") +
   ggtitle(" ")+
   xlab("Pseudoyear ID")+
   ylab("Pseudoyear effect ")+
   theme_light()# + ylim(-0.6,0.6) + theme(plot.title = element_text(hjust=0.5))

## save
ggsave(paste0(out_02,"/02 besag model 1A_psy psyear RW2 effect.png"), dpi=300, width=6, height=6)

## CPUE mean -------------------------------------------------------------------

# map SPmean
ggplot()+
   geom_sf(data=data_out, aes(geometry=geometry, fill=CPUE_mean))+
   scale_fill_viridis_c(option="H") + 
   geom_sf(data = world_map_data,alpha=0.8) +
   theme_light() +
   coord_sf(expand=FALSE, xlim = c(10, 115), ylim = c(-45,40))+
   labs(fill="Mean CPUE")

## save
ggsave(paste0(out_02,"/02 besag model 1A_y fitted CPUE mean.png"), dpi=300, width=6, height=6)

## CPUE sd ---------------------------------------------------------------------

# map SPmean
ggplot()+
   geom_sf(data=data_out, aes(geometry=geometry, fill=CPUE_sd))+
   scale_fill_viridis_c(option="H") + 
   geom_sf(data = world_map_data,alpha=0.8) +
   theme_light() +
   coord_sf(expand=FALSE, xlim = c(10, 115), ylim = c(-45,40))+
   labs(fill="Sd CPUE")

## save
ggsave(paste0(out_02,"/02 besag model 1A_y fitted CPUE sd.png"), dpi=300, width=6, height=6)

# ips --------------------------------------------------------------------------

n=100
ips = inla.posterior.sample(n=n,model)

# prediction -------------------------------------------------------------------

## we take samples of the posterior distribution for all model components
## we combine the components to build the linear predictor and des-transform

t_fun=function(x){exp(x)} ## exponent for gamma

library(dplyr)

psam <- sapply(ips, function(x){
   ## we define function inside of sapply ()
   
   # x<-ips[[1]] # ips object
   
   ## intercept
   intercept <- x$latent %>% rownames(.) %>% stringr::str_detect("^\\(Intercept\\)") %>% x$latent[.,]
   
   ## rw2 year effect
   library(dplyr)
   rw<-x$latent %>% rownames(.) %>% stringr::str_detect("^psyearID") %>% x$latent[.,]
   
   ## besag poly effect
   besag<-x$latent %>% rownames(.) %>% stringr::str_detect("^polyID") %>% x$latent[.,]
   
   ## df predictor
   predictor<-matrix(0, nrow=length(unique(data$polyID)), 
                     ncol=length(unique(data$psyearID)) )
   
   ##  loop through each year + iid for each cell
   for (s in 1:length(unique(data$polyID))){
      
      for (t in 1:length(unique(data$psyearID))){ 
         
         ## build linear predictor
         predictor[s,t] = intercept + rw[t] + besag[s]   
         
      }   
   }
   
   ## loop to pass from matrix to vector
   pred=0
   for (i in 1:length(unique(data$psyearID))){
      pred <- c(pred, predictor[,i])
      
   }
   ## remove first 0 row
   pred=pred[-1]
   ## apply des-transform function
   t_fun(pred)
   
})

## in each column we have each of the 100 samples
## in each row we loop each value of s and t
psam

## calculate median and quantiles of the 100 samples
q.sam_al_a <- apply(psam, 1, quantile,
                    c(.025, 0.05, 0.5, 0.95, .975), na.rm =TRUE)

## transpose matrix and pass to data.frame
df_plot <- data.frame(t(q.sam_al_a)) # quantile (cols) and space and time in rows
df_plot$sd<-apply(psam, 1, sd)
df_plot$pseudoyear<-sort(rep(unique(data$pseudoyear), (length(unique(data$polyID)))))
df_plot$polyID<-(rep(unique(data$polyID),length(unique(data$psyearID))))

## add geometry polys to predicted dataset
data_pred<-left_join(df_plot, poly_grid, by="polyID")

## create levels for year, season and combined yearseas
lyear<-rep(seq(from=1972, to=2015), each=4); length(lyear) # 176
lseas<-rep(1:4, times=44); length(lseas) # 176
lyearseas<-paste(lyear, lseas, sep=".");  length(lyearseas) 

## replicate factor columns in data frame
data_pred$year<-as.factor(data_pred$pseudoyear); levels(data_pred$year)
data_pred$seas<-as.factor(data_pred$pseudoyear); levels(data_pred$seas)
data_pred$yearseas<-as.factor(data_pred$pseudoyear); levels(data_pred$yearseas)

## rename factor levels
levels(data_pred$year) <- lyear
levels(data_pred$seas)<-lseas
levels(data_pred$yearseas)<-lyearseas

## add region factor
df_reg<-unique(data[,c(1,11)]) # poly ID (160) by region
colnames(df_reg)<-c("polyID","regionID")
data_pred<-left_join(data_pred, df_reg, by="polyID")

# plot pred --------------------------------------------------------------------

## CPUE idx --------------------------------------------------------------------

## 1 area ----------------------------------------------------------------------
library(plyr)
data_idx_1<-ddply(data_pred, .(pseudoyear, yearseas, year, seas), summarize,  "meanX50"=mean(X50.),
                  "meanX5"=mean(X5.), "meanX95"=mean(X95.), "meanSd"=mean(sd))
data_idx_1$pseudoyear<-as.numeric(data_idx_1$pseudoyear)

ggplot(data = data_idx_1, aes(x = pseudoyear, y = meanX50, group=1))+
   geom_line(aes(x=pseudoyear, y=meanX50), size=1)+
   geom_ribbon(aes(x = pseudoyear, ymin = meanX5, ymax = meanX95), 
               alpha = 0.3,fill="grey")+
   ylab("Predicted median CPUE")+
   xlab("Pseudoyear")+
   theme_light() + scale_x_continuous(breaks= seq(from=1, to=length(unique(data$psyearID)), by=12),
                                      labels=seq(from=81, to=256, by=12),
                                      limits=c(1, length(unique(data$pseudoyear))))

## save
ggsave(paste0(out_02,"/02 besag model predicted CPUE index_1 area.png"), dpi=300, 
       width=6, height=6)

## 4 area ----------------------------------------------------------------------

library(plyr)
data_idx_4<-ddply(data_pred, .(pseudoyear, yearseas, year, seas, regionID), summarize,  "meanX50"=mean(X50.),
                  "meanX5"=mean(X5.), "meanX95"=mean(X95.), "meanSd"=mean(sd))
data_idx_4$pseudoyear<-as.numeric(data_idx_4$pseudoyear)

ggplot(data = data_idx_4, aes(x = pseudoyear, y = meanX50, group=1))+
   geom_line(aes(x=pseudoyear, y=meanX50), size=1)+
   geom_ribbon(aes(x = pseudoyear, ymin = meanX5, ymax = meanX95), 
               alpha = 0.3,fill="grey")+
   ylab("Predicted median CPUE")+
   xlab("Pseudoyear")+
   facet_wrap(~regionID)+
   theme_light() + scale_x_continuous(breaks= seq(from=1, to=length(unique(data$psyearID)), by=20),
                                      labels=seq(from=81, to=256, by=20),
                                      limits=c(1, length(unique(data$pseudoyear))))

## save
ggsave(paste0(out_02,"/02 besag model predicted CPUE index_ 4 areas.png"), dpi=300, 
       width=6, height=6)

## CPUE st ---------------------------------------------------------------------

library(sf)
library(ggspatial)
library(rnaturalearth)
library(maps)
library(ggplot2)

## world map 
world_map_data <- ne_countries(scale = "medium", returnclass = "sf")

## define scale limits and breaks
lmin=min(data_pred$X50.)
lmax=max(data_pred$X50.)
bks=c(0,50,100,150,200,250,300,350,400)

## plot all years predicted CPUE
ggplot()+
   geom_sf(data=data_pred,aes(geometry=geometry, fill=X50.), alpha=0.8) +
   scale_fill_viridis_c(option="H", limits = c(lmin, lmax),
                        breaks=bks,
                        oob = scales::squish) + 
   geom_sf(data = world_map_data,alpha=0.8) +
   theme_light() +
   coord_sf(expand=FALSE, xlim = c(10, 115), ylim = c(-45,40))+
   xlab("longitude") + ylab("latitude") +   labs(fill=' Pred. CPUE')+

## save
ggsave(paste0(out_02, "/02 besag model predicted CPUE median (years collapsed).png")
       , dpi=300,  width=6, height=6)



## now plot predicted CPUE map for each year

## create dir
out_022<-paste(getwd(),"/output/02 besag model/predicted CPUE st", sep="")
dir.create(out_022)

## loop plot
for (j in 81:256) {
   
   data_pred_y<-subset(data_pred, data_pred$pseudoyear==j)
   ggplot()+
      geom_sf(data=data_pred_y,aes(geometry=geometry, fill=X50.), alpha=0.8) +
      scale_fill_viridis_c(option="H", limits = c(lmin, lmax),
                           breaks=bks,
                           oob = scales::squish) + 
      geom_sf(data = world_map_data,alpha=0.8) +
      theme_light() +
      coord_sf(expand=FALSE, xlim = c(10, 115), ylim = c(-45,40))+
      xlab("longitude") + ylab("latitude") +   labs(fill=' Pred. CPUE')+
      ggtitle(paste0(j))
   
   ggsave(paste0(out_022,"/predicted CPUE psyear",j,".png"))
   
}

## magick ----------------------------------------------------------------------

## use magick package to create a gift with all year maps
library(magick)

## list file names and read in
imgs <- list.files(out_022, full.names = TRUE)

## join the images together
img_joined <- image_join(image_read(imgs[1:176])) # 1 to 176 pseudoyears

## animate at 2 frames per second
img_animated <- image_animate(img_joined, fps = 2)

## view animated image
#img_animated

## save to disk
image_write(image = img_animated,
            path = paste0(out_02,"/predicted CPUE by pseudoyear.gif"))

# save/load --------------------------------------------------------------------

save.image(file = paste0(out_02,"/02 besag model.RData"))
#load("./output/02 besag model/02 besag model.RData")

# IMAGE GIFT https://www.nagraj.net/notes/gifs-in-r/

# references -------------------------------------------------------------------

# https://becarioprecario.bitbucket.io/inla-gitbook/ch-temporal.html#sec:spacetime
# https://becarioprecario.bitbucket.io/inla-gitbook/ch-spatial.html
# https://cran.r-project.org/web/packages/spdep/vignettes/nb_sf.html
# https://github.com/jmartinez-minaya/inla-course/blob/master/r/4_besag.R
