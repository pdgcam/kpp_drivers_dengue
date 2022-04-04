library(readxl)
library(PBSmapping)
library(raster)
library(ggplot2)
library(INLA)
library(geoR)
library(FNN)
library(lme4)

rm(list=ls())

kfcshh = read.csv("~/LandUseKFCSApr2021.csv", header = TRUE, sep = ",", as.is=TRUE, stringsAsFactors = FALSE)
load("~/mosquito_data.Rdata")

kfcshh = read.csv("~/LandUseKFCSApr2021.csv", header = TRUE, sep = ",", as.is=TRUE, stringsAsFactors = FALSE)
load("~/0_Data/mosquito_data.Rdata")


#Removing newborns
kfcshh = kfcshh[which(kfcshh$ageyrs>=1),]
length(unique(kfcshh$houseNo))

#Removing >30 yo
kfcshh = kfcshh[which(kfcshh$ageyrs<=30),]
length(unique(kfcshh$houseNo))

#Adding meter coordinates
coords_home = kfcshh[,c("gpslong","lat")]
colnames(coords_home)<-c("X","Y")
attr(coords_home,"projection")<-"LL"
coords_home<-convUL(coords_home,km=F)
kfcshh[,c("x_proj","y_proj")] = coords_home

#Run models


#Constant FOI-----------------------------------------
formula = seroposenr ~ 1 + offset(log(ageyrs))
output_constant_foi = inla(formula, data = kfcshh, family = "binomial", control.family = list(link="cloglog"), control.predictor = list(compute=T), verbose = T)
#save(kfcshh, output_constant_foi, file = "~/output_constant_foi.RData")
est_constant_foi = output_constant_foi$summary.fixed
#Fit
age_vect = seq(0,100,by=0.1)
fit_age = 1-exp(-exp(est_constant_foi$mean) * age_vect)
fit_age.loci = 1-exp(-exp(est_constant_foi$`0.025quant`) * age_vect)
fit_age.hici = 1-exp(-exp(est_constant_foi$`0.975quant`) * age_vect)
prop_by_age = tapply(kfcshh$seroposenr, round(kfcshh$ageyrs), mean)


#Adjust for households------------------------
load("~/3_INLA/output_hhadjusted.RData")
output_hhadjusted$.args$formula
est_hhadjusted = output_hhadjusted$summary.fixed


#Proportion captured by household variability - getting residuals
load("~/3_INLA/output_hhadjusted.RData")
hh_vect = sapply(kfcshh$houseNo, function(x) which(output_hhadjusted$summary.random$houseNo$ID == x))
foi_vect = exp(output_hhadjusted$summary.fixed$mean + output_hhadjusted$summary.random$houseNo$mean[hh_vect])
prob_vect = 1 - exp(- foi_vect * kfcshh$ageyrs)
hist(prob_vect, 100)
residuals_hhadjusted = kfcshh$seroposenr - prob_vect
variance_hhadjusted = mean(residuals_hhadjusted^2)
plot(residuals_hhadjusted)
variance_hhadjusted = mean(residuals_hhadjusted^2)

load("~/3_INLA/output_constant_foi.RData")
output_constant_foi$.args$formula
foi_vect = exp(output_constant_foi$summary.fixed$mean)
prob_vect = 1 - exp(- foi_vect * kfcshh$ageyrs) 
hist(prob_vect, 100)
residuals_constantfoi = kfcshh$seroposenr - prob_vect
plot(residuals_constantfoi)
variance_constantfoi = mean(residuals_constantfoi^2)
(1 - (variance_hhadjusted / variance_constantfoi))*100 #Percentage of variance explained by household adjustment


#Including older poeple
load("~/3_INLA/output_hhadjusted_witholder.RData")
hh_vect = sapply(kfcshh$houseNo, function(x) which(output_hhadjusted_witholder$summary.random$houseNo$ID == x))
foi_vect = exp(output_hhadjusted_witholder$summary.fixed$mean + output_hhadjusted_witholder$summary.random$houseNo$mean[hh_vect])
prob_vect = 1 - exp(- foi_vect * kfcshh$ageyrs)
hist(prob_vect, 100)
residuals_hhadjusted = kfcshh$seroposenr - prob_vect
variance_hhadjusted = mean(residuals_hhadjusted^2)
plot(residuals_hhadjusted)
variance_hhadjusted = mean(residuals_hhadjusted^2)

kfcshh = read.csv("~/LandUseKFCSApr2021.csv", header = TRUE, sep = ",", as.is=TRUE, stringsAsFactors = FALSE)
kfcshh = kfcshh[which(kfcshh$ageyrs>=1),]

formula = seroposenr ~ 1 + offset(log(ageyrs))
output_constant_foi = inla(formula, data = kfcshh, family = "binomial", control.family = list(link="cloglog"), control.predictor = list(compute=T), verbose = T)
foi_vect = exp(output_constant_foi$summary.fixed$mean)
prob_vect = 1 - exp(- foi_vect * kfcshh$ageyrs)
hist(prob_vect, 100)
residuals_constantfoi = kfcshh$seroposenr - prob_vect
plot(residuals_constantfoi)
variance_constantfoi = mean(residuals_constantfoi^2)
(1 - (variance_hhadjusted / variance_constantfoi))*100 #Percentage of variance explained by household adjustment



#Proportion captured by spatial field- getting residuals------------------------------------------
load("~/3_INLA/output_constant_foi.RData")
output_constant_foi$.args$formula
nrow(output_constant_foi$.args$data)
foi_vect = exp(output_constant_foi$summary.fixed$mean)
prob_vect = 1 - exp(- foi_vect * kfcshh$ageyrs)
hist(prob_vect, 100)
residuals_constantfoi = kfcshh$seroposenr - prob_vect
plot(residuals_constantfoi)
variance_constantfoi = mean(residuals_constantfoi^2)

load("~/3_INLA/64_univariate_nohh_baseline_model_1.RData")
output$.args$formula
mesh_res = 4000
locs = as.matrix(kfcshh[,c("x_proj","y_proj")])
mesh = inla.mesh.2d(loc=kfcshh[,c("x_proj","y_proj")], max.edge= c(mesh_res,mesh_res*5),
                    cutoff=mesh_res)
A.est = inla.spde.make.A(mesh=mesh,
                         loc=as.matrix(kfcshh[,c("x_proj","y_proj")]))
spde = spde = inla.spde2.matern(mesh=mesh,alpha=2)
output$.args$formula
length(output$.args$data$seroposenr)
spatial.field.contrib = output$summary.random$`spatial.field`[,c("ID","mean")]
foi_vect = exp(output$summary.fixed$mean + as.vector(A.est%*%spatial.field.contrib$mean))
prob_vect = 1 - exp(- foi_vect * kfcshh$ageyrs)
hist(prob_vect, 100)
residuals_spacedep = kfcshh$seroposenr - prob_vect
variance_spacedep = mean(residuals_spacedep^2)
(1 - (variance_spacedep/ variance_constantfoi))*100 #Percentage of variance explained by spatial field when adjusting by household

#Including older people
load("~/3_INLA/output_spacedep_foi_nohh_witholder.RData")
formula = seroposenr ~ 1 + offset(log(ageyrs))
output_constant_foi = inla(formula, data = kfcshh, family = "binomial", control.family = list(link="cloglog"), control.predictor = list(compute=T), verbose = T)
foi_vect = exp(output_constant_foi$summary.fixed$mean)
prob_vect = 1 - exp(- foi_vect * kfcshh$ageyrs)
hist(prob_vect, 100)
residuals_constantfoi = kfcshh$seroposenr - prob_vect
plot(residuals_constantfoi)
variance_constantfoi = mean(residuals_constantfoi^2)


output = output_spacedep_foi
mesh_res = 4000
locs = as.matrix(kfcshh[,c("x_proj","y_proj")])
mesh = inla.mesh.2d(loc=kfcshh[,c("x_proj","y_proj")], max.edge= c(mesh_res,mesh_res*5),
                    cutoff=mesh_res)
A.est = inla.spde.make.A(mesh=mesh,
                         loc=as.matrix(kfcshh[,c("x_proj","y_proj")]))
spde = spde = inla.spde2.matern(mesh=mesh,alpha=2)
output$.args$formula
length(output$.args$data$seroposenr)
spatial.field.contrib = output$summary.random$`spatial.field`[,c("ID","mean")]
foi_vect = exp(output$summary.fixed$mean + as.vector(A.est%*%spatial.field.contrib$mean))
prob_vect = 1 - exp(- foi_vect * kfcshh$ageyrs)
hist(prob_vect, 100)
residuals_spacedep = kfcshh$seroposenr - prob_vect
variance_spacedep = mean(residuals_spacedep^2)
(1 - (variance_spacedep/ variance_constantfoi))*100 #Percentage of variance explained by spatial field when adjusting by household


output_constant_foi$neffp
output.constant.summ = inla.collect.results(inla=output_constant_foi)
variance_hhadjusted = a
variance_spacedep = inla.emarginal(function(x) x,output.field$marginals.variance.nominal[[1]])
(1 - (variance_spacedep/ variance_constantfoi))*100 #Percentage of variance explained by spatial field when adjusting by household



#Variogram-----------------------------------------
plot(kfcshh[,c("x_proj","y_proj")])
coords = kfcshh[,c("x_proj","y_proj")]
data = kfcshh$seroposenr
variogram = variog(coords = coords, data = data,  breaks = seq(0,7e4,by=5000))
plot(variogram)

#Get minimum distance between two different households
unique_hh = kfcshh[sapply(unique(kfcshh$houseNo), function(x) which(kfcshh$houseNo == x)[1]),]
plot(unique_hh[,c("x_proj","y_proj")])
nn_hh = get.knn(unique_hh[,c("x_proj","y_proj")], 1)
sort(nn_hh$nn.dist) #Some households are  co-located => Need to differentiate ppl from the same hh and ppl from different 
nn_hh$nn.index[nn_hh$nn.dist==0]
which(nn_hh$nn.dist==0)
#Jitter by household
jitter_x = runif(length(unique(kfcshh$houseNo)[which(nn_hh$nn.dist==0)]),-1,1)
jitter_y = runif(length(unique(kfcshh$houseNo)[which(nn_hh$nn.dist==0)]),-1,1)

for(x in 1:length(jitter_x)){
  kfcshh$x_proj[kfcshh$houseNo == unique(kfcshh$houseNo)[which(nn_hh$nn.dist==0)][x]] = kfcshh$x_proj[kfcshh$houseNo == unique(kfcshh$houseNo)[which(nn_hh$nn.dist==0)][x]] + jitter_x[x]
  kfcshh$y_proj[kfcshh$houseNo == unique(kfcshh$houseNo)[which(nn_hh$nn.dist==0)][x]] = kfcshh$y_proj[kfcshh$houseNo == unique(kfcshh$houseNo)[which(nn_hh$nn.dist==0)][x]] + jitter_y[x]
}
#Check again
unique_hh = kfcshh[sapply(unique(kfcshh$houseNo), function(x) which(kfcshh$houseNo == x)[1]),]
plot(unique_hh[,c("x_proj","y_proj")])
nn_hh = get.knn(unique_hh[,c("x_proj","y_proj")], 1)
sort(nn_hh$nn.dist) #Some households are  co-located => Need to differentiate ppl from the same hh and ppl from different 
nn_hh$nn.index[nn_hh$nn.dist==0]
which(nn_hh$nn.dist==0)

#2nd variogram
coords = kfcshh[,c("x_proj","y_proj")]
data = kfcshh$seroposenr
variogram = variog(coords = coords, data = data,  breaks = seq(0,7e4,by=5000))
plot(variogram)

#Get range
semivar = variogram$v
dist = variogram$u
fit <- nls(formula = semivar ~ A + (S-A)*exp(-dist/delta), start = list(S = semivar[1], A = semivar[length(semivar)], delta=1e4))
S = coef(fit)[1]
A = coef(fit)[2]
delta = coef(fit)[3]
dist_fit = seq(0,10e4,by=100)
fitted = A + (S-A)*exp(-dist_fit/delta)

ggplot() +
  geom_point(aes(x=dist, y=semivar), color = "black", size=3) +
  geom_line(aes(x=dist_fit, y=fitted), color = "black", size=1) +
  geom_hline(yintercept = A) +
  geom_vline(xintercept = signif(3*delta, digits = 3), col=2) +
  #geom_abline(intercept = S, slope = (A-S)/delta) +
  theme_minimal() +
  scale_x_continuous(name="Distance (m)", limits = c(0,7e4), breaks=c(0,3*delta), labels = c(0,signif(3*delta, digits = 3))) +
  scale_y_continuous(name="Variance") +
  theme(axis.title.x = element_text(size=25), axis.text.x = element_text(angle=0, size=15)) +
  theme(axis.title.y = element_text(size=25), axis.text.y = element_text(angle=0, size=15)) +
  theme(panel.grid.minor = element_blank())

dev.copy(png, '~/2_Plots/variogram.png', width=7.5*75, height=6.81*75)
dev.off()


#Table univariate----------------------
list = list.files("~/3_INLA/univariate")

covar_vect = c("popdens","percity","toilets","hhsize","ncontainers","ncontainersplastic","ncontainersbottles","housetype","garbagemanagement","gender","occupation",
               "hheducation","constructionmaterial","rooftype","watertype","waterressources","doorscreens","devnot","immature_mosquitoes","immature_strict","aegypti_captured","agecat","plastic_burn")


list_effects = list()
for(i in 1:length(list)){
  print(i)
  load(paste0("~/3_INLA/univariate/",list[i]))
  covar = covar_vect[covar_id]
  if(covar %in% c("popdens","percity","toilets","ncontainersbottles","hhsize","garbagemanagement","ncontainers","ncontainersplastic","gender",
                  "constructionmaterial","rooftype","watertype","waterressources","doorscreens","devnot","immature_mosquitoes","immature_strict","aegypti_captured")){
    list_effects[[covar]][[model_id]] = cbind(data.frame(ID="fixed"),output$summary.fixed[2,c("mean","0.025quant","0.975quant")])
  }else{
    if(model_id!=3){
      list_effects[[covar]][[model_id]] = output$summary.random[[3]][,c("ID","mean","0.025quant","0.975quant")]
    }else{
      list_effects[[covar]][[model_id]] = output$summary.random[[3]][,c("ID","mean","0.025quant","0.975quant")]
    }
  }
}

list_effects$agecat[[1]] = list_effects$agecat[[2]] = list_effects$agecat[[3]]
nlevels = sapply(list_effects, function(x) nrow(x[[1]]))
df_list_effects = rbind(data.frame(covariate = rep(names(list_effects), nlevels), model = 1),data.frame(covariate = rep(names(list_effects), nlevels), model = 2),data.frame(covariate = rep(names(list_effects), nlevels), model = 3))
df_list_effects$ID = c(unlist(lapply(list_effects, function(x) x[[1]][,"ID"])),unlist(lapply(list_effects, function(x) x[[2]][,"ID"])),unlist(lapply(list_effects, function(x) x[[3]][,"ID"])))
df_list_effects$mean = c(unlist(lapply(list_effects, function(x) x[[1]][,"mean"])),unlist(lapply(list_effects, function(x) x[[2]][,"mean"])),unlist(lapply(list_effects, function(x) x[[3]][,"mean"])))
df_list_effects$cilo = c(unlist(lapply(list_effects, function(x) x[[1]][,"0.025quant"])),unlist(lapply(list_effects, function(x) x[[2]][,"0.025quant"])),unlist(lapply(list_effects, function(x) x[[3]][,"0.025quant"])))
df_list_effects$cihi = c(unlist(lapply(list_effects, function(x) x[[1]][,"0.975quant"])),unlist(lapply(list_effects, function(x) x[[2]][,"0.975quant"])),unlist(lapply(list_effects, function(x) x[[3]][,"0.975quant"])))

write.csv(df_list_effects, file="~/3_INLA/univariate/univariate_table.csv")

df_list_selec = df_list_effects[df_list_effects$ID=="fixed",]

y_min = df_list_selec$cilo
y_max = df_list_selec$cihi
y_mean = df_list_selec$mean
x_ax = rep(1:(nrow(df_list_selec)/3),3) + rep(c(-0.1,0,0.1),each=nrow(df_list_selec)/3)

ggplot() +
  geom_pointrange(aes(x=x_ax, y = y_mean, ymin =y_min, ymax = y_max), size=1, col=rep(1:3, each=nrow(df_list_selec)/3)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  scale_x_continuous(name=NULL, breaks=1:(nrow(df_list_selec)/3), labels = df_list_selec$covariate[1:(nrow(df_list_selec)/3)]) +
  scale_y_continuous(name=NULL) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle=45, size=20, vjust=1, hjust=1), axis.title.x = element_text(size=20)) +
  theme(axis.title.y = element_text(size=20), axis.text.y = element_text(size=20))  +
  theme(plot.title = element_text(size=15, hjust=0.5)) +
  theme(panel.grid.minor = element_blank())

dev.copy(png, paste0('~/2_Plots/covariates/univariate_fixed.png'), width=7.5*75, height=6.81*75)
dev.off()

df_list_selec = df_list_effects[df_list_effects$ID!="fixed",]
for(i in 1:length(unique(df_list_selec$covariate))){
  df_list_selecc = df_list_selec[which(df_list_selec$covariate == unique(df_list_selec$covariate)[i]),]
  
  y_min = df_list_selecc$cilo
  y_max = df_list_selecc$cihi
  y_mean = df_list_selecc$mean
  x_ax = rep(1:(nrow(df_list_selecc)/3),3) + rep(c(-0.1,0,0.1),each=nrow(df_list_selecc)/3)
  
  p = ggplot() +
    geom_pointrange(aes(x=x_ax, y = y_mean, ymin =y_min, ymax = y_max), size=1, col=rep(1:3, each=nrow(df_list_selecc)/3)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    theme_minimal() +
    scale_x_continuous(name=NULL, breaks=1:(nrow(df_list_selecc)/3), labels = df_list_selecc$ID[1:(nrow(df_list_selecc)/3)]) +
    scale_y_continuous(name=NULL) +
    theme(legend.position = "none") +
    ggtitle(df_list_selecc$covariate[1]) +
    theme(axis.text.x = element_text(angle=45, size=20, vjust=1, hjust=1), axis.title.x = element_text(size=20)) +
    theme(axis.title.y = element_text(size=20), axis.text.y = element_text(size=20))  +
    theme(plot.title = element_text(size=30, hjust=0)) +
    theme(panel.grid.minor = element_blank())
  print(p)
  
  dev.copy(png, paste0('~/2_Plots/covariates/univariate_',df_list_selecc$covariate[1],'.png'), width=7.5*75, height=6.81*75)
  dev.off()
}




#Table univariate no hh----------------------
list = list.files("~/3_INLA/univariate_nohh")

covar_vect = c("popdens","percity","toilets","hhsize","ncontainers","ncontainersplastic","ncontainersbottles","housetype","garbagemanagement","gender","occupation",
               "hheducation","constructionmaterial","rooftype","watertype","waterressources","doorscreens","devnot","immature_mosquitoes","immature_strict","aegypti_captured","agecat")


list_effects = list()
for(i in 1:length(list)){
  print(i)
  load(paste0("~/3_INLA/univariate_nohh/",list[i]))
  covar = covar_vect[covar_id]
  if(covar %in% c("popdens","percity","toilets","ncontainersbottles","hhsize","garbagemanagement","ncontainers","ncontainersplastic","gender","constructionmaterial","rooftype","watertype","waterressources","doorscreens","devnot","immature_mosquitoes","immature_strict","aegypti_captured")){
    list_effects[[covar]][[model_id]] = cbind(data.frame(ID="fixed"),output$summary.fixed[2,c("mean","0.025quant","0.975quant")])
  }else{
    if(model_id!=3){
      list_effects[[covar]][[model_id]] = output$summary.random[[2]][,c("ID","mean","0.025quant","0.975quant")]
    }else{
      list_effects[[covar]][[model_id]] = output$summary.random[[2]][,c("ID","mean","0.025quant","0.975quant")]
    }
  }
}

list_effects$agecat[[1]] = list_effects$agecat[[2]] = list_effects$agecat[[3]]
nlevels = sapply(list_effects, function(x) nrow(x[[1]]))
df_list_effects = rbind(data.frame(covariate = rep(names(list_effects), nlevels), model = 1),data.frame(covariate = rep(names(list_effects), nlevels), model = 2),data.frame(covariate = rep(names(list_effects), nlevels), model = 3))
df_list_effects$ID = c(unlist(lapply(list_effects, function(x) x[[1]][,"ID"])),unlist(lapply(list_effects, function(x) x[[2]][,"ID"])),unlist(lapply(list_effects, function(x) x[[3]][,"ID"])))
df_list_effects$mean = c(unlist(lapply(list_effects, function(x) x[[1]][,"mean"])),unlist(lapply(list_effects, function(x) x[[2]][,"mean"])),unlist(lapply(list_effects, function(x) x[[3]][,"mean"])))
df_list_effects$cilo = c(unlist(lapply(list_effects, function(x) x[[1]][,"0.025quant"])),unlist(lapply(list_effects, function(x) x[[2]][,"0.025quant"])),unlist(lapply(list_effects, function(x) x[[3]][,"0.025quant"])))
df_list_effects$cihi = c(unlist(lapply(list_effects, function(x) x[[1]][,"0.975quant"])),unlist(lapply(list_effects, function(x) x[[2]][,"0.975quant"])),unlist(lapply(list_effects, function(x) x[[3]][,"0.975quant"])))

write.csv(df_list_effects, file="~/3_INLA/univariate_nohh/univariate_nohh_table.csv")

df_list_selec = df_list_effects[df_list_effects$ID=="fixed",]

y_min = df_list_selec$cilo
y_max = df_list_selec$cihi
y_mean = df_list_selec$mean
x_ax = rep(1:(nrow(df_list_selec)/3),3) + rep(c(-0.1,0,0.1),each=nrow(df_list_selec)/3)

ggplot() +
  geom_pointrange(aes(x=x_ax, y = y_mean, ymin =y_min, ymax = y_max), size=1, col=rep(1:3, each=nrow(df_list_selec)/3)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  scale_x_continuous(name=NULL, breaks=1:(nrow(df_list_selec)/3), labels = df_list_selec$covariate[1:(nrow(df_list_selec)/3)]) +
  scale_y_continuous(name=NULL) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle=45, size=20, vjust=1, hjust=1), axis.title.x = element_text(size=20)) +
  theme(axis.title.y = element_text(size=20), axis.text.y = element_text(size=20))  +
  theme(plot.title = element_text(size=15, hjust=0.5)) +
  theme(panel.grid.minor = element_blank())

dev.copy(png, paste0('~/2_Plots/covariates/univariate_nohh_fixed.png'), width=7.5*75, height=6.81*75)
dev.off()

df_list_selec = df_list_effects[df_list_effects$ID!="fixed",]
for(i in 1:length(unique(df_list_selec$covariate))){
  df_list_selecc = df_list_selec[which(df_list_selec$covariate == unique(df_list_selec$covariate)[i]),]
  
  y_min = df_list_selecc$cilo
  y_max = df_list_selecc$cihi
  y_mean = df_list_selecc$mean
  x_ax = rep(1:(nrow(df_list_selecc)/3),3) + rep(c(-0.1,0,0.1),each=nrow(df_list_selecc)/3)
  
  p = ggplot() +
    geom_pointrange(aes(x=x_ax, y = y_mean, ymin =y_min, ymax = y_max), size=1, col=rep(1:3, each=nrow(df_list_selecc)/3)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    theme_minimal() +
    scale_x_continuous(name=NULL, breaks=1:(nrow(df_list_selecc)/3), labels = df_list_selecc$ID[1:(nrow(df_list_selecc)/3)]) +
    scale_y_continuous(name=NULL) +
    theme(legend.position = "none") +
    ggtitle(df_list_selecc$covariate[1]) +
    theme(axis.text.x = element_text(angle=45, size=20, vjust=1, hjust=1), axis.title.x = element_text(size=20)) +
    theme(axis.title.y = element_text(size=20), axis.text.y = element_text(size=20))  +
    theme(plot.title = element_text(size=30, hjust=0)) +
    theme(panel.grid.minor = element_blank())
  print(p)
  
  dev.copy(png, paste0('~/2_Plots/covariates/univariate_nohh_',df_list_selecc$covariate[1],'.png'), width=7.5*75, height=6.81*75)
  dev.off()
}




#Table univariate sensitivity analysis lowering threshold----------------------
list = list.files("~/3_INLA/univariate_senshai")

covar_vect = c("popdens","percity","toilets","hhsize","ncontainers","ncontainersplastic","ncontainersbottles","housetype","garbagemanagement","gender","occupation",
               "hheducation","constructionmaterial","rooftype","watertype","waterressources","doorscreens","devnot","immature_mosquitoes","immature_strict","aegypti_captured","agecat","plastic_burn")


list_effects = list()
for(i in 1:length(list)){
  print(i)
  load(paste0("~/3_INLA/univariate_senshai/",list[i]))
  covar = covar_vect[covar_id]
  if(covar %in% c("popdens","percity","toilets","ncontainersbottles","hhsize","garbagemanagement","ncontainers","ncontainersplastic","gender",
                  "constructionmaterial","rooftype","watertype","waterressources","doorscreens","devnot","immature_mosquitoes","immature_strict","aegypti_captured")){
    list_effects[[covar]][[model_id]] = cbind(data.frame(ID="fixed"),output$summary.fixed[2,c("mean","0.025quant","0.975quant")])
  }else{
    if(model_id!=3){
      list_effects[[covar]][[model_id]] = output$summary.random[[3]][,c("ID","mean","0.025quant","0.975quant")]
    }else{
      list_effects[[covar]][[model_id]] = output$summary.random[[3]][,c("ID","mean","0.025quant","0.975quant")]
    }
  }
}

list_effects$agecat[[1]] = list_effects$agecat[[2]] = list_effects$agecat[[3]]
nlevels = sapply(list_effects, function(x) nrow(x[[1]]))
df_list_effects = rbind(data.frame(covariate = rep(names(list_effects), nlevels), model = 1),data.frame(covariate = rep(names(list_effects), nlevels), model = 2),data.frame(covariate = rep(names(list_effects), nlevels), model = 3))
df_list_effects$ID = c(unlist(lapply(list_effects, function(x) x[[1]][,"ID"])),unlist(lapply(list_effects, function(x) x[[2]][,"ID"])),unlist(lapply(list_effects, function(x) x[[3]][,"ID"])))
df_list_effects$mean = c(unlist(lapply(list_effects, function(x) x[[1]][,"mean"])),unlist(lapply(list_effects, function(x) x[[2]][,"mean"])),unlist(lapply(list_effects, function(x) x[[3]][,"mean"])))
df_list_effects$cilo = c(unlist(lapply(list_effects, function(x) x[[1]][,"0.025quant"])),unlist(lapply(list_effects, function(x) x[[2]][,"0.025quant"])),unlist(lapply(list_effects, function(x) x[[3]][,"0.025quant"])))
df_list_effects$cihi = c(unlist(lapply(list_effects, function(x) x[[1]][,"0.975quant"])),unlist(lapply(list_effects, function(x) x[[2]][,"0.975quant"])),unlist(lapply(list_effects, function(x) x[[3]][,"0.975quant"])))

write.csv(df_list_effects, file="~/3_INLA/univariate_senshai/univariate_senshai_table.csv")

df_list_selec = df_list_effects[df_list_effects$ID=="fixed",]

y_min = df_list_selec$cilo
y_max = df_list_selec$cihi
y_mean = df_list_selec$mean
x_ax = rep(1:(nrow(df_list_selec)/3),3) + rep(c(-0.1,0,0.1),each=nrow(df_list_selec)/3)

ggplot() +
  geom_pointrange(aes(x=x_ax, y = y_mean, ymin =y_min, ymax = y_max), size=1, col=rep(1:3, each=nrow(df_list_selec)/3)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  scale_x_continuous(name=NULL, breaks=1:(nrow(df_list_selec)/3), labels = df_list_selec$covariate[1:(nrow(df_list_selec)/3)]) +
  scale_y_continuous(name=NULL) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle=45, size=20, vjust=1, hjust=1), axis.title.x = element_text(size=20)) +
  theme(axis.title.y = element_text(size=20), axis.text.y = element_text(size=20))  +
  theme(plot.title = element_text(size=15, hjust=0.5)) +
  theme(panel.grid.minor = element_blank())

dev.copy(png, paste0('~/2_Plots/covariates/univariate_senshai_fixed.png'), width=7.5*75, height=6.81*75)
dev.off()

df_list_selec = df_list_effects[df_list_effects$ID!="fixed",]
for(i in 1:length(unique(df_list_selec$covariate))){
  df_list_selecc = df_list_selec[which(df_list_selec$covariate == unique(df_list_selec$covariate)[i]),]
  
  y_min = df_list_selecc$cilo
  y_max = df_list_selecc$cihi
  y_mean = df_list_selecc$mean
  x_ax = rep(1:(nrow(df_list_selecc)/3),3) + rep(c(-0.1,0,0.1),each=nrow(df_list_selecc)/3)
  
  p = ggplot() +
    geom_pointrange(aes(x=x_ax, y = y_mean, ymin =y_min, ymax = y_max), size=1, col=rep(1:3, each=nrow(df_list_selecc)/3)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    theme_minimal() +
    scale_x_continuous(name=NULL, breaks=1:(nrow(df_list_selecc)/3), labels = df_list_selecc$ID[1:(nrow(df_list_selecc)/3)]) +
    scale_y_continuous(name=NULL) +
    theme(legend.position = "none") +
    ggtitle(df_list_selecc$covariate[1]) +
    theme(axis.text.x = element_text(angle=45, size=20, vjust=1, hjust=1), axis.title.x = element_text(size=20)) +
    theme(axis.title.y = element_text(size=20), axis.text.y = element_text(size=20))  +
    theme(plot.title = element_text(size=30, hjust=0)) +
    theme(panel.grid.minor = element_blank())
  print(p)
  
  dev.copy(png, paste0('~/2_Plots/covariates/univariate_senshai',df_list_selecc$covariate[1],'.png'), width=7.5*75, height=6.81*75)
  dev.off()
}




#Table multivariate----------------------
rm(list = ls())

load("~/3_INLA/multivariate/multiv_1_model_1.RData")
output = output_gridcode

df_list_effects = data.frame(covariate = rownames(output$summary.fixed[-1,c("mean","0.025quant","0.975quant")]),
                             model=1, ID = "fixed")
df_list_effects[,c("mean","0.025quant","0.975quant")] = output$summary.fixed[-1,c("mean","0.025quant","0.975quant")]
nlevels = sapply(output$summary.random[c(-1,-2)], function(x) nrow(x))
df_temp = data.frame(covariate = rep(names(output$summary.random[c(-1,-2)]), nlevels),
                     model=1, ID = unlist(sapply(output$summary.random[c(-1,-2)], function(x) x$ID)))
df_temp[,c("mean","0.025quant","0.975quant")] = rbind(output$summary.random[-1]$houseType.y[,c("mean","0.025quant","0.975quant")],output$summary.random[-1]$occupation[,c("mean","0.025quant","0.975quant")],output$summary.random[-1]$education_household[,c("mean","0.025quant","0.975quant")])
df_list_effects = rbind(df_list_effects, df_temp)

df_list_multiv = df_list_effects

load("~/3_INLA/multivariate/multiv_1_model_2.RData")
output = output_gridcode 

df_list_effects = data.frame(covariate = rownames(output$summary.fixed[-1,c("mean","0.025quant","0.975quant")]),
                             model=2, ID = "fixed")
df_list_effects[,c("mean","0.025quant","0.975quant")] = output$summary.fixed[-1,c("mean","0.025quant","0.975quant")]
nlevels = sapply(output$summary.random[c(-1,-2)], function(x) nrow(x))
df_temp = data.frame(covariate = rep(names(output$summary.random[c(-1,-2)]), nlevels),
                     model=2, ID = unlist(sapply(output$summary.random[c(-1,-2)], function(x) x$ID)))
df_temp[,c("mean","0.025quant","0.975quant")] = rbind(output$summary.random[-1]$houseType.y[,c("mean","0.025quant","0.975quant")],output$summary.random[-1]$occupation[,c("mean","0.025quant","0.975quant")],output$summary.random[-1]$education_household[,c("mean","0.025quant","0.975quant")])
df_list_effects = rbind(df_list_effects, df_temp)

df_list_multiv = rbind(df_list_multiv, df_list_effects)

load("~/3_INLA/multivariate/multiv_1_model_3.RData")
output = output_gridcode 

df_list_effects = data.frame(covariate = rownames(output$summary.fixed[-1,c("mean","0.025quant","0.975quant")]),
                             model=3, ID = "fixed")
df_list_effects[,c("mean","0.025quant","0.975quant")] = output$summary.fixed[-1,c("mean","0.025quant","0.975quant")]
nlevels = sapply(output$summary.random[c(-1,-2)], function(x) nrow(x))
df_temp = data.frame(covariate = rep(names(output$summary.random[c(-1,-2)]), nlevels),
                     model=3, ID = unlist(sapply(output$summary.random[c(-1,-2)], function(x) x$ID)))
df_temp[,c("mean","0.025quant","0.975quant")] = rbind(output$summary.random[-1]$houseType.y[,c("mean","0.025quant","0.975quant")],output$summary.random[-1]$occupation[,c("mean","0.025quant","0.975quant")],
                                                      output$summary.random[-1]$education_household[,c("mean","0.025quant","0.975quant")], output$summary.random[-1]$age_cat[,c("mean","0.025quant","0.975quant")])
df_list_effects = rbind(df_list_effects, df_temp)


df_list_multiv = rbind(df_list_multiv, df_list_effects)

add_fake_agecat = rbind(df_list_multiv[df_list_multiv$covariate=="age_cat",],df_list_multiv[df_list_multiv$covariate=="age_cat",])
add_fake_agecat$model = rep(c(1,2),each=6)
df_list_multiv = rbind(df_list_multiv, add_fake_agecat)
write.csv(df_list_multiv, file="~/3_INLA/multivariate/multivariate_table.csv")

df_list_selec = df_list_multiv[df_list_multiv$ID=="fixed",]

y_min = df_list_selec$`0.025quant`
y_max = df_list_selec$`0.975quant`
y_mean = df_list_selec$mean
x_ax = rep(1:(nrow(df_list_selec)/3),3) + rep(c(-0.1,0,0.1),each=nrow(df_list_selec)/3)
lab_x = c("popdens","toilets","ncontainers","garbagemanagement","sex","constructionmaterial","rooftype","watertype","waterressources","doorscreens")


ggplot() +
  geom_pointrange(aes(x=x_ax, y = y_mean, ymin =y_min, ymax = y_max), size=1, col=rep(1:3, each=nrow(df_list_selec)/3)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  scale_x_continuous(name=NULL, breaks=1:(nrow(df_list_selec)/3), labels = lab_x) +
  scale_y_continuous(name=NULL) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle=45, size=20, vjust=1, hjust=1), axis.title.x = element_text(size=20)) +
  theme(axis.title.y = element_text(size=20), axis.text.y = element_text(size=20))  +
  theme(plot.title = element_text(size=15, hjust=0.5)) +
  theme(panel.grid.minor = element_blank())

dev.copy(png, paste0('~/2_Plots/covariates/multivariate_fixed.png'), width=7.5*75, height=6.81*75)
dev.off()

df_list_selec = df_list_multiv[df_list_multiv$ID!="fixed",]
for(i in 1:length(unique(df_list_selec$covariate))){
  print(i)
  df_list_selecc = df_list_selec[which(df_list_selec$covariate == unique(df_list_selec$covariate)[i]),]
  
  y_min = df_list_selecc$`0.025quant`
  y_max = df_list_selecc$`0.975quant`
  y_mean = df_list_selecc$mean
  x_ax = rep(1:(nrow(df_list_selecc)/3),3) + rep(c(-0.1,0,0.1),each=nrow(df_list_selecc)/3)
  
  p = ggplot() +
    geom_pointrange(aes(x=x_ax, y = y_mean, ymin =y_min, ymax = y_max), size=1, col=rep(1:3, each=nrow(df_list_selecc)/3)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    theme_minimal() +
    scale_x_continuous(name=NULL, breaks=1:(nrow(df_list_selecc)/3), labels = df_list_selecc$ID[1:(nrow(df_list_selecc)/3)]) +
    scale_y_continuous(name=NULL) +
    theme(legend.position = "none") +
    ggtitle(df_list_selecc$covariate[1]) +
    theme(axis.text.x = element_text(angle=45, size=20, vjust=1, hjust=1), axis.title.x = element_text(size=20)) +
    theme(axis.title.y = element_text(size=20), axis.text.y = element_text(size=20))  +
    theme(plot.title = element_text(size=30, hjust=0)) +
    theme(panel.grid.minor = element_blank())
  print(p)
  
  dev.copy(png, paste0('~/2_Plots/covariates/multivariate_',df_list_selecc$covariate[1],'.png'), width=7.5*75, height=6.81*75)
  dev.off()
}

#Table multivariate sensitivity analysis lowering threshold----------------------
rm(list = ls())

load("~/3_INLA/multivariate_senshai/multiv_1_model_1_senshai.RData")
output = output_gridcode

df_list_effects = data.frame(covariate = rownames(output$summary.fixed[-1,c("mean","0.025quant","0.975quant")]),
                             model=1, ID = "fixed")
df_list_effects[,c("mean","0.025quant","0.975quant")] = output$summary.fixed[-1,c("mean","0.025quant","0.975quant")]
nlevels = sapply(output$summary.random[c(-1,-2)], function(x) nrow(x))
df_temp = data.frame(covariate = rep(names(output$summary.random[c(-1,-2)]), nlevels),
                     model=1, ID = unlist(sapply(output$summary.random[c(-1,-2)], function(x) x$ID)))
df_temp[,c("mean","0.025quant","0.975quant")] = rbind(output$summary.random[-1]$houseType.y[,c("mean","0.025quant","0.975quant")],output$summary.random[-1]$occupation[,c("mean","0.025quant","0.975quant")],output$summary.random[-1]$education_household[,c("mean","0.025quant","0.975quant")])
df_list_effects = rbind(df_list_effects, df_temp)

df_list_multiv = df_list_effects

load("~/3_INLA/multivariate_senshai/multiv_1_model_2_senshai.RData")
output = output_gridcode
df_list_effects = data.frame(covariate = rownames(output$summary.fixed[-1,c("mean","0.025quant","0.975quant")]),
                             model=2, ID = "fixed")
df_list_effects[,c("mean","0.025quant","0.975quant")] = output$summary.fixed[-1,c("mean","0.025quant","0.975quant")]
nlevels = sapply(output$summary.random[c(-1,-2)], function(x) nrow(x))
df_temp = data.frame(covariate = rep(names(output$summary.random[c(-1,-2)]), nlevels),
                     model=2, ID = unlist(sapply(output$summary.random[c(-1,-2)], function(x) x$ID)))
df_temp[,c("mean","0.025quant","0.975quant")] = rbind(output$summary.random[-1]$houseType.y[,c("mean","0.025quant","0.975quant")],output$summary.random[-1]$occupation[,c("mean","0.025quant","0.975quant")],output$summary.random[-1]$education_household[,c("mean","0.025quant","0.975quant")])
df_list_effects = rbind(df_list_effects, df_temp)

df_list_multiv = rbind(df_list_multiv, df_list_effects)

load("~/3_INLA/multivariate_senshai/multiv_1_model_3_senshai.RData")
output = output_gridcode
df_list_effects = data.frame(covariate = rownames(output$summary.fixed[-1,c("mean","0.025quant","0.975quant")]),
                             model=3, ID = "fixed")
df_list_effects[,c("mean","0.025quant","0.975quant")] = output$summary.fixed[-1,c("mean","0.025quant","0.975quant")]
nlevels = sapply(output$summary.random[c(-1,-2)], function(x) nrow(x))
df_temp = data.frame(covariate = rep(names(output$summary.random[c(-1,-2)]), nlevels),
                     model=3, ID = unlist(sapply(output$summary.random[c(-1,-2)], function(x) x$ID)))
df_temp[,c("mean","0.025quant","0.975quant")] = rbind(output$summary.random[-1]$houseType.y[,c("mean","0.025quant","0.975quant")],output$summary.random[-1]$occupation[,c("mean","0.025quant","0.975quant")],
                                                      output$summary.random[-1]$education_household[,c("mean","0.025quant","0.975quant")], output$summary.random[-1]$age_cat[,c("mean","0.025quant","0.975quant")])
df_list_effects = rbind(df_list_effects, df_temp)


df_list_multiv = rbind(df_list_multiv, df_list_effects)

add_fake_agecat = rbind(df_list_multiv[df_list_multiv$covariate=="age_cat",],df_list_multiv[df_list_multiv$covariate=="age_cat",])
add_fake_agecat$model = rep(c(1,2),each=6)
df_list_multiv = rbind(df_list_multiv, add_fake_agecat)
write.csv(df_list_multiv, file="~/3_INLA/multivariate_senshai/multivariate_senshai_table.csv")

df_list_selec = df_list_multiv[df_list_multiv$ID=="fixed",]

y_min = df_list_selec$`0.025quant`
y_max = df_list_selec$`0.975quant`
y_mean = df_list_selec$mean
x_ax = rep(1:(nrow(df_list_selec)/3),3) + rep(c(-0.1,0,0.1),each=nrow(df_list_selec)/3)
lab_x = c("popdens","toilets","hhsize","ncontainers","garbagemanagement","sex","constructionmaterial","rooftype","watertype","waterressources","doorscreens","devnot","immature_mosquitoes","aegypti_captured")



#Table proportions-------------------------------------------------------------------
kfcshh = read.csv("~/LandUseKFCSApr2021.csv", header = TRUE, sep = ",", as.is=TRUE, stringsAsFactors = FALSE)
#Removing newborns
kfcshh = kfcshh[which(kfcshh$ageyrs>=1),]
table(kfcshh$seroposenr, kfcshh$gender)
load("~/0_Data/mosquito_data.Rdata")
#Remove the .x columns
duplicated =  (1:ncol(kfcshh)) %in% grep(".x", names(kfcshh))
kfcshh = kfcshh[,!duplicated]

#Format covariates
kfcshh$toilets_outside = kfcshh$numToiletsOutside.y>0
kfcshh$hh_size = cut(kfcshh$numSubjects, breaks = c(0,3,6,13))
#agecat
agecat_vect = seq(5,80,by=5)
kfcshh$agecat = cut(kfcshh$ageyrs, breaks=1+c(-2,agecat_vect,200), right=F)
#Nwatercontainers
cov_watercontainers = kfcshh[,grep("waterContainer", names(kfcshh))]
containers_total = rowSums(cov_watercontainers[,-1], na.rm=T)
containers_plastic = cov_watercontainers$waterContainerSoftDrink.y + cov_watercontainers$waterContainerGallons.y
kfcshh$containers_total_uncut = containers_total
kfcshh$containers_total = cut(containers_total, breaks = c(seq(0,50,by=10),250))
#Nwatercontainers - plastic
cov_watercontainers = kfcshh[,grep("waterContainer", names(kfcshh))]
containers_plastic = cov_watercontainers$waterContainerSoftDrink.y + cov_watercontainers$waterContainerGallons.y
kfcshh$containers_plastic = cut(containers_plastic, breaks = c(-1,0,10,20,200))
kfcshh$containers_bottles = kfcshh$waterContainerSoftDrink.y>0
#Garbage management
garbage = as.factor(kfcshh$garbageManagement.y)
levels(garbage) = c("Burn/Bury/Dump","Burn/Bury/Dump","Car","Burn/Bury/Dump")
kfcshh$garbage = garbage!="Car"
kfcshh$gender_bin = kfcshh$gender != "F"
#Occupation
cov_occupation = kfcshh[,grep("occupation", names(kfcshh))]
for(col in 1:ncol(cov_occupation)){cov_occupation[,col] = as.numeric(cov_occupation[,col])}
for(col in 1:ncol(cov_occupation)){cov_occupation[,col][cov_occupation[,col]!=1] = 0}
is_farmer = cov_occupation$occupationFarmer+cov_occupation$occupationAnimalFarmer>=1
is_student = cov_occupation$occupationStudent >= 1
is_unemployed = cov_occupation$occupationUnemployed >= 1
is_employee = is_farmer + is_student + is_unemployed ==0
occupation = rep(NA, nrow(kfcshh))
occupation[is_farmer] = "Farmer"
occupation[is_unemployed] = "Unemployed"
occupation[is_student] = "Student"
occupation[is_employee] = "Employee"
kfcshh$occupation = as.factor(occupation)
#hh_education
educ_per_hh = lapply(unique(kfcshh$houseNo), function(hh) kfcshh$education[kfcshh$houseNo==hh])
higher_educ_per_hh = sapply(educ_per_hh, function(x){
  if("Master" %in% x) return("Master")
  else if("Bachelor" %in% x) return("Bachelor")
  else if("High" %in% x) return("High")
  else if("Secondary" %in% x) return("Secondary")
  else if("Primary" %in% x) return("Primary")
  else if("NoSchool" %in% x) return("No School")
  else return("Unknown")
})
kfcshh$education_household = as.factor(sapply(kfcshh$houseNo, function(x) higher_educ_per_hh[which(unique(kfcshh$houseNo)==x)]))
levels(kfcshh$education_household) = c("bachelor/master","high","bachelor/master","noschool/primary","noschool/primary","secondary")
#construction material
kfcshh$not_concrete = kfcshh$constructionMaterialsConcrete.y == 0
#roof type
kfcshh$zinc_roof = kfcshh$roofTypeZinc.y>0
#water type
kfcshh$water_type_pipe = kfcshh$waterTypePipe.y>0
#Is there a year where no immature mosquitoes were found in the house?
immature_presence_summary = sapply(unique(immature_presence$houseNo), function(x) prod(immature_presence$larvae_pupae_presence[which(immature_presence$houseNo==x)])>0)
kfcshh$immature_mosquitoes = sapply(kfcshh$houseNo, function(x){
  if(x %in% unique(immature_presence$houseNo)) return(immature_presence_summary[which(unique(immature_presence$houseNo)==x)])
  else return(NA)})
#No immature mosquitoes were ever found in this house
immature_presence_summary = sapply(unique(immature_presence$houseNo), function(x) sum(immature_presence$larvae_pupae_presence[which(immature_presence$houseNo==x)])>0)
kfcshh$immature_mosquitoes_strict = sapply(kfcshh$houseNo, function(x){
  if(x %in% unique(immature_presence$houseNo)) return(immature_presence_summary[which(unique(immature_presence$houseNo)==x)])
  else return(NA)})

aegypti_captured = survey_tab_adult_16_to_20[survey_tab_adult_16_to_20$`Genus species`=="aegypti",]
aegypti_captured$`# specimen` = as.numeric(aegypti_captured$`# specimen`)
#add houses where no aegypti were captured
summary_houses = summary_houses[!is.na(summary_houses$houseNo),]
empty_houses = lapply(2016:2020, function(year){
  summary_houses$houseNo[which(summary_houses$year==year)][which(!summary_houses$houseNo[which(summary_houses$year==year)] %in% aegypti_captured$`House code`[aegypti_captured$Year==year])]
})
length_by_year = sapply(empty_houses, length)
empty_df = data.frame(year=rep(2016:2020,length_by_year),
                      houseNo=unlist(empty_houses),
                      Genus = "Aedes",
                      species = "aegypti",
                      Gender = "both",
                      captured = 0)
names(empty_df) = names(aegypti_captured)
aegypti_captured = rbind(aegypti_captured, empty_df)

#Normalise by total mosquitoes captured during a given year
aegypti_captured$norm_count = aegypti_captured$`# specimen`/tapply(as.numeric(aegypti_captured$`# specimen`), aegypti_captured$Year, mean,na.rm=T)[sapply(aegypti_captured$Year, function(x) which(2016:2020==x))]
#Take average across years per household
aegypti_captured_summary = sapply(unique(aegypti_captured$`House code`), function(x) mean(aegypti_captured$norm_count[which(aegypti_captured$`House code`==x)]))
aegypti_captured_notnorm_summary = sapply(unique(aegypti_captured$`House code`), function(x) mean(aegypti_captured$`# specimen`[which(aegypti_captured$`House code`==x)]))

kfcshh$aegypti_captured_norm = sapply(kfcshh$houseNo, function(x){
  if(x %in% unique(aegypti_captured$`House code`)) return(aegypti_captured_summary[which(unique(aegypti_captured$`House code`)==x)])
  else return(NA)})
kfcshh$aegypti_captured = sapply(kfcshh$houseNo, function(x){
  if(x %in% unique(aegypti_captured$`House code`)) return(aegypti_captured_notnorm_summary[which(unique(aegypti_captured$`House code`)==x)])
  else return(NA)})
kfcshh$aegypti_captured_cat = cut(round(kfcshh$aegypti_captured), breaks=c(-1,0,2,5,40))
table(kfcshh$aegypti_captured_cat)

names(kfcshh)
tab_prop = sapply(names(kfcshh), function(covar) table(kfcshh[,covar], kfcshh$seroposenr ))
of_interest = c("gender", "agecat", "occupation", "hh_size", "containers_total","containers_plastic","containers_bottles",
                "houseType.y", "garbageManagement.y", "education_household", "not_concrete","zinc_roof", "waterResourcesNearby.y",
                "water_type_pipe", "doorScreens.y", "immature_mosquitoes" , "immature_mosquitoes_strict", "aegypti_captured_cat", "toilets_outside",
                "gridcodecut", "PerCityCat", "devnot")
df_prop = data.frame(covariate=NA, id=NA,"0"=NA,"1"=NA)
for(covar in of_interest){
  print(covar)
  df_prop = rbind(df_prop, data.frame(covariate=covar,id=rownames(tab_prop[covar][[1]]), "0"=tab_prop[covar][[1]][,1], "1" = tab_prop[covar][[1]][,2]))
}                
df_prop$total = df_prop$X0 + df_prop$X1                
df_prop$prop_neg = df_prop$X0/df_prop$total
df_prop$prop_pos = df_prop$X1/df_prop$total
write.csv(df_prop, file="~/0_Data/tab_proportions.csv")





#Other results-------------------
kfcshh = read.csv("~/LandUseKFCSApr2021.csv", header = TRUE, sep = ",", as.is=TRUE, stringsAsFactors = FALSE)
#Removing newborns
kfcshh = kfcshh[which(kfcshh$ageyrs>=1),]

#Adding meter coordinates
coords_home = kfcshh[,c("gpslong","lat")]
colnames(coords_home)<-c("X","Y")
attr(coords_home,"projection")<-"LL"
coords_home<-convUL(coords_home,km=F)
kfcshh[,c("x_proj","y_proj")] = coords_home

#Age and sex dist
mean(kfcshh$ageyrs)
table(kfcshh$gender)/nrow(kfcshh)
tapply(kfcshh$ageyrs, kfcshh$gender, mean)
#Prop infected per year
load("~/0_Data/infections_df.RData")
df_pop$primary_infections = rowSums(matrix(unlist(lapply(infections_df, function(x) x$n1)), ncol=nrow(thai_structure), byrow=F))
df_pop$secondary_infections = rowSums(matrix(unlist(lapply(infections_df, function(x) x$n2)), ncol=nrow(thai_structure), byrow=F))
df_pop$tertiary_infections = rowSums(matrix(unlist(lapply(infections_df, function(x) x$n3)), ncol=nrow(thai_structure), byrow=F))
df_pop$quaternary_infections = rowSums(matrix(unlist(lapply(infections_df, function(x) x$n4)), ncol=nrow(thai_structure), byrow=F))
df_pop = df_pop[!is.na(df_pop$pop_density),]
df_pop$total_infections = df_pop$primary_infections + df_pop$secondary_infections + df_pop$tertiary_infections + df_pop$quaternary_infections
sum(df_pop$total_infections)/sum(df_pop$pop_density) * 100
sum(df_pop$total_infections)
sum(df_pop$primary_infections) + sum(df_pop$secondary_infections)

#Reproductive base number
lambda = c(0.116,0.126,0.136, 0.115,0.124,0.133)
age_vect = 0:99
f = rep(thai_structure$prop/5, each=5)
x = sapply(1:6, function(i) exp(-lambda[i]*age_vect))
w1 = sapply(1:6, function(i) 4*exp(-0.75*lambda[i]*age_vect)*(1-exp(-0.25*lambda[i]*age_vect)))
w2 = sapply(1:6, function(i) 6*exp(-0.5*lambda[i]*age_vect)*(1-exp(-0.25*lambda[i]*age_vect))^2)
w3 = sapply(1:6, function(i) 4*exp(-0.25*lambda[i]*age_vect)*(1-exp(-0.25*lambda[i]*age_vect))^3)
R0 = sapply(1:6, function(i) 1 / sum(f*(x[,i]+0.75*w1[,i]+0.5*w2[,i]+0.25*w3[,i])))
R0

#Prop detected by hospitals
hosp_cases = read.csv("~/0_Data/tac.csv", header = TRUE, sep = ",", as.is=TRUE, stringsAsFactors = FALSE)
hosp_cases_kpp = hosp_cases[which(hosp_cases$area=="KamphaengPhet"),]
tapply(hosp_cases_kpp$count, hosp_cases_kpp$year, sum)
mean(tapply(hosp_cases_kpp$count, hosp_cases_kpp$year, sum)[-length(tapply(hosp_cases_kpp$count, hosp_cases_kpp$year, sum))]) #exclude 2021
mean(tapply(hosp_cases_kpp$count, hosp_cases_kpp$year, sum)[-length(tapply(hosp_cases_kpp$count, hosp_cases_kpp$year, sum))])/sum(df_pop$primary_infections + df_pop$secondary_infections) * 100 #exclude 2021
mean(tapply(hosp_cases_kpp$count, hosp_cases_kpp$year, sum)[-length(tapply(hosp_cases_kpp$count, hosp_cases_kpp$year, sum))])/sum(df_pop$total_infections) * 100 #exclude 2021

#Most common professions
sum(as.numeric((kfcshh$occupationAnimalFarmer[kfcshh$ageyrs>18]>0|kfcshh$occupationFarmer[kfcshh$ageyrs>18]>0)))/sum(kfcshh$ageyrs>18)*100
sum(as.numeric((kfcshh$occupationGeneralEmployee[kfcshh$ageyrs>18]>0|kfcshh$occupationCompanyEmployee[kfcshh$ageyrs>18]>0)))/sum(kfcshh$ageyrs>18)*100
sum(as.numeric(kfcshh$occupationUnemployed[kfcshh$ageyrs>18]>0))
#Zinc roof
hh_df = kfcshh[!duplicated(kfcshh$houseNo),]
sum(hh_df$roofTypeZinc.y>0)/nrow(hh_df)*100
sum(hh_df$numToiletsOutside.y>0)/nrow(hh_df)*100
table(hh_df$garbageManagement.y)/nrow(hh_df)*100
#Water containers
mean(hh_df$containers_total_uncut, na.rm=T)
#Devnot
mean(hh_df$devnot)
#Pop density per household
pop_hh = get.knnx(df_pop[,c("x","y")], hh_df[,c("x_proj", "y_proj")], 1)
pop_hh = df_pop$pop_density[pop_hh$nn.index]
mean(pop_hh)
range(pop_hh)
#Prop population living within 1km of a household
pop_around_hh = get.knnx(df_pop[,c("x","y")], hh_df[,c("x_proj", "y_proj")], nrow(df_pop))
pop_around_hh$nn.dist[1,]
dist_vect = seq(500,30000,by=500)
prop_by_dist = sapply(dist_vect, function(radius) sum(df_pop$pop_density[unique(pop_around_hh$nn.index[pop_around_hh$nn.dist<radius])])/sum(df_pop$pop_density)*100)
plot(dist_vect, prop_by_dist)
prop_by_dist[which(dist_vect==1000)]

#Prop seropositive and p value by gender
mean(kfcshh$seroposenr)
sum(kfcshh$seroposenr)
nrow(kfcshh)
table(kfcshh$seroposenr, kfcshh$gender)
1236/(163+1236)
798/(179+798)
prop.test(x = c(1236, 798), n = c((163+1236), (179+798)))
#Prop seropositive by age
age_cat = cut(kfcshh$ageyrs, breaks=1+c(-2,seq(5,30,by=5),105), right=F)#
table(age_cat)
prop_by_age = tapply(kfcshh$seroposenr, age_cat, mean)
length_age = tapply(kfcshh$seroposenr, age_cat, length)
ci_asym = Hmisc::binconf(x=c(prop_by_age*length_age), n=c(length_age))
ci_asym
#constant foi
load("~/3_INLA/output_constant_foi.RData")
exp(output_constant_foi$summary.fixed)
#range spatial field
range(df_pop$lambda)
#Proportion of cases happening within muang
kpp_city =  data.frame(x=555784.61884666, y=1822411.20896305)
dist_from_city = get.knnx(kpp_city, df_pop[,c("x","y")], 1)
vect_dist = seq(1,100,by=1)*1000
case_by_dist = sapply(vect_dist, function(dist) sum(df_pop$total_infections[dist_from_city$nn.dist<dist]))
prop_by_dist = case_by_dist/sum(df_pop$total_infections)*100
plot(prop_by_dist)
prop_by_dist[which(vect_dist==5000)]
#Doorscreens
mean(hh_df$doorScreens.y)
#Proportion hh with mosquito study
length(unique(summary_houses$houseNo[summary_houses$mosquito_stage=="Adult"]))
length(unique(summary_houses$houseNo[summary_houses$mosquito_stage=="Adult"]))/length(unique(kfcshh$houseNo))
length(unique(survey_tab_adult_16_to_20$`House code`))
length(unique(survey_tab_adult_16_to_20$`House code`))/length(unique(summary_houses$houseNo[summary_houses$mosquito_stage=="Adult"]))

#Prop hh with burnt plastic
sum(hh_df$garbageManagement.y=="Burn"&hh_df$waterContainerSoftDrink.y>0)/nrow(hh_df)
17.9
#Mosquito data
load("~/0_Data/mosquito_data.Rdata")
table(survey_tab_adult_16_to_20$Year)
#add houses where no aegypti were captured
aegypti_captured = survey_tab_adult_16_to_20[survey_tab_adult_16_to_20$`Genus species`=="aegypti",]
aegypti_captured$`# specimen` = as.numeric(aegypti_captured$`# specimen`)
tapply(as.numeric(aegypti_captured$`# specimen`),aegypti_captured$Gender,sum)/sum(as.numeric(aegypti_captured$`# specimen`)) #57% aegypti captured were female
table(aegypti_captured$Gender, aegypti_captured$`# specimen`)
summary_houses = summary_houses[!is.na(summary_houses$houseNo),]
empty_houses = lapply(2016:2020, function(year){
  summary_houses$houseNo[which(summary_houses$year==year)][which(!summary_houses$houseNo[which(summary_houses$year==year)] %in% aegypti_captured$`House code`[aegypti_captured$Year==year])]
})

length_by_year = sapply(empty_houses, length)
names(aegypti_captured)
empty_df = data.frame(year=rep(2016:2020,length_by_year),
                      houseNo=unlist(empty_houses),
                      Genus = "Aedes",
                      species = "aegypti",
                      Gender = "both",
                      captured = 0)
names(empty_df) = names(aegypti_captured)
aegypti_captured = rbind(aegypti_captured, empty_df)
mean(aegypti_captured$`# specimen`, na.rm=T)
range(aegypti_captured$`# specimen`, na.rm=T)
median(aegypti_captured$`# specimen`, na.rm=T)

