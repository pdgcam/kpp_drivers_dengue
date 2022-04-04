library(readxl)
library(PBSmapping)
library(raster)
library(ggplot2)
library(INLA)

kfcshh = read.csv("LandUseKFCSApr2021.csv", header = TRUE, sep = ",", as.is=TRUE, stringsAsFactors = FALSE)

#Removing newborns
kfcshh = kfcshh[which(kfcshh$ageyrs>=1),]

#Adding meter coordinates
coords_home = kfcshh[,c("gpslong","lat")]
colnames(coords_home)<-c("X","Y")
attr(coords_home,"projection")<-"LL"
coords_home<-convUL(coords_home,km=F)
kfcshh[,c("x_proj","y_proj")] = coords_home

formula = seroposenr ~ 1 + offset(log(ageyrs)) + f(houseNo, model="iid")
output_hhadjusted_witholder = inla(formula, data = kfcshh, family = "binomial", control.family = list(link="cloglog"), control.predictor = list(compute=T), verbose = T)
save(kfcshh, output_hhadjusted_witholder, file = "output_hhadjusted_witholder.RData")

formula = seroposenr ~ -1 + Intercept + offset(log(ageyrs)) + f(spatial.field, model=spde)
mesh_res = 4000
locs = as.matrix(kfcshh[,c("x_proj","y_proj")])
mesh = inla.mesh.2d(loc=kfcshh[,c("x_proj","y_proj")], max.edge= c(mesh_res,mesh_res*5),
                    cutoff=mesh_res)
spde = inla.spde2.matern(mesh=mesh,alpha=2)
s.index <- inla.spde.make.index(name="spatial.field", n.spde=spde$n.spde) 
A.est = inla.spde.make.A(mesh=mesh,
                         loc=as.matrix(kfcshh[,c("x_proj","y_proj")]))
stack.est = inla.stack(data=list(seroposenr=kfcshh$seroposenr),
                       A=list(A.est,
                              1),
                       effects=list(c(s.index, list(Intercept=1)),
                                    list(kfcshh[,which(names(kfcshh)!= "seroposenr")])),
                       tag='stdata')

output_spacedep_foi = inla(formula,
                           data =inla.stack.data(stack.est),
                           family = "binomial",
                           control.family = list(link="cloglog"),
                           control.predictor=list(A=inla.stack.A(stack.est), compute = F),
                           verbose=T)

save(output_spacedep_foi, kfcshh, A.est, file="output_spacedep_foi_nohh_witholder.RData")

#Removing >30 yo
kfcshh = kfcshh[which(kfcshh$ageyrs<=30),]

formula = seroposenr ~ 1 + offset(log(ageyrs)) + f(houseNo, model="iid")
output_hhadjusted = inla(formula, data = kfcshh, family = "binomial", control.family = list(link="cloglog"), control.predictor = list(compute=T), verbose = T)
save(kfcshh, output_hhadjusted, file = "output_hhadjusted.RData")

formula = seroposenr ~ -1 + Intercept + offset(log(ageyrs)) + f(spatial.field, model=spde) + f(houseNo, model="iid")
mesh_res = 4000
locs = as.matrix(kfcshh[,c("x_proj","y_proj")])
mesh = inla.mesh.2d(loc=kfcshh[,c("x_proj","y_proj")], max.edge= c(mesh_res,mesh_res*5),
                    cutoff=mesh_res)
spde = inla.spde2.matern(mesh=mesh,alpha=2)
s.index <- inla.spde.make.index(name="spatial.field", n.spde=spde$n.spde) 
A.est = inla.spde.make.A(mesh=mesh,
                         loc=as.matrix(kfcshh[,c("x_proj","y_proj")]))
stack.est = inla.stack(data=list(seroposenr=kfcshh$seroposenr),
                       A=list(A.est,
                              1),
                       effects=list(c(s.index, list(Intercept=1)),
                                    list(kfcshh[,which(names(kfcshh)!= "seroposenr")])),
                       tag='stdata')

output_spacedep_foi = inla(formula,
                           data =inla.stack.data(stack.est),
                           family = "binomial",
                           control.family = list(link="cloglog"),
                           control.predictor=list(A=inla.stack.A(stack.est), compute = F),
                           verbose=T)

save(output_spacedep_foi, kfcshh, A.est, file="output_spacedep_foi.RData")