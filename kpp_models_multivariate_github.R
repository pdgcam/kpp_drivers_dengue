library(readxl)
library(INLA)
library(PBSmapping)
library(ggplot2)

load("mosquito_data.Rdata")
kfcshh = read.csv("LandUseKFCSApr2021.csv", header = TRUE, sep = ",", as.is=TRUE, stringsAsFactors = FALSE)
kfcshh$gridcode = kfcshh$gridcode/1000*100 #Make it from per 1hab/100m^2 to per 1000hab/km^2


#Removing newborns
kfcshh = kfcshh[which(kfcshh$ageyrs>=1),]
#Adding meter coordinates
coords_home = kfcshh[,c("gpslong","lat")]
colnames(coords_home)<-c("X","Y")
attr(coords_home,"projection")<-"LL"
coords_home<-convUL(coords_home,km=F)
kfcshh[,c("x_proj","y_proj")] = coords_home
#Remove the .x columns
duplicated =  (1:ncol(kfcshh)) %in% grep(".x", names(kfcshh))
kfcshh = kfcshh[,!duplicated]

#Format covariates
kfcshh$toilets_outside = kfcshh$numToiletsOutside.y>0
kfcshh$hh_size = cut(kfcshh$numSubjects, breaks = c(0,3,6,13))
#Nwatercontainers
cov_watercontainers = kfcshh[,grep("waterContainer", names(kfcshh))]
containers_total = rowSums(cov_watercontainers[,-1], na.rm=T)
containers_plastic = cov_watercontainers$waterContainerSoftDrink.y + cov_watercontainers$waterContainerGallons.y
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
kfcshh$aegypti_captured_norm = sapply(kfcshh$houseNo, function(x){
  if(x %in% unique(aegypti_captured$`House code`)) return(aegypti_captured_summary[which(unique(aegypti_captured$`House code`)==x)])
  else return(NA)})
#Plastic_Litter
kfcshh$plastic_litter = kfcshh$garbageManagement.y=="Burn" & kfcshh$waterContainerSoftDrink.y>0

#Run
#Model 2
formula = seroposenr ~ -1 + Intercept + f(spatial.field, model=spde) + offset(log(ageyrs)) + f(houseNo, model = 'iid') +
  #(as.numeric(gridcode)) + 
  #(as.numeric(perCity)) + 
  (as.numeric(toilets_outside)) +
  #(as.numeric(hh_size)) +
  (as.numeric(containers_total)) +
  #(as.numeric(containers_plastic)) +
  #(as.numeric(containers_bottles)) +
  f(houseType.y, model='iid', constr = F, extraconstr = list(A=matrix(c(1e4,rep(1,length(unique(kfcshh$houseType.y))-1)), nrow=1), e=rep(0,1))) +
  (as.numeric(garbage)) +
  (as.numeric(gender_bin)) +
  f(occupation, model='iid', constr = F, extraconstr = list(A=matrix(c(1e4,rep(1,length(unique(kfcshh$occupation))-1)), nrow=1), e=rep(0,1))) +
  f(education_household, model='iid', constr = F, extraconstr = list(A=matrix(c(1e4,rep(1,length(unique(kfcshh$education_household))-1)), nrow=1), e=rep(0,1))) +
  (as.numeric(not_concrete)) +
  (as.numeric(zinc_roof)) +
  (as.numeric(water_type_pipe)) +
  (as.numeric(waterResourcesNearby.y)) +
  (as.numeric(doorScreens.y)) +
  (as.numeric(devnot))
#(as.numeric(immature_mosquitoes)) +
#(as.numeric(immature_mosquitoes_strict)) +
#(as.numeric(aegypti_captured_norm)) +
#(as.numeric(plastic_litter))

formula_gridcode = seroposenr ~ -1 + Intercept + f(spatial.field, model=spde) + offset(log(ageyrs)) + f(houseNo, model = 'iid') +
  (as.numeric(gridcode)) + 
  #(as.numeric(perCity)) + 
  (as.numeric(toilets_outside)) +
  #(as.numeric(hh_size)) +
  (as.numeric(containers_total)) +
  #(as.numeric(containers_plastic)) +
  #(as.numeric(containers_bottles)) +
  f(houseType.y, model='iid', constr = F, extraconstr = list(A=matrix(c(1e4,rep(1,length(unique(kfcshh$houseType.y))-1)), nrow=1), e=rep(0,1))) +
  (as.numeric(garbage)) +
  (as.numeric(gender_bin)) +
  f(occupation, model='iid', constr = F, extraconstr = list(A=matrix(c(1e4,rep(1,length(unique(kfcshh$occupation))-1)), nrow=1), e=rep(0,1))) +
  f(education_household, model='iid', constr = F, extraconstr = list(A=matrix(c(1e4,rep(1,length(unique(kfcshh$education_household))-1)), nrow=1), e=rep(0,1))) +
  (as.numeric(not_concrete)) +
  (as.numeric(zinc_roof)) +
  (as.numeric(water_type_pipe)) +
  (as.numeric(waterResourcesNearby.y)) +
  (as.numeric(doorScreens.y))
#(as.numeric(devnot))
#(as.numeric(immature_mosquitoes)) +
#(as.numeric(immature_mosquitoes_strict)) +
#(as.numeric(aegypti_captured_norm)) +
#(as.numeric(plastic_litter))


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


output = inla(formula,
              data =inla.stack.data(stack.est),
              family = "binomial",
              control.family = list(link="cloglog"),
              control.predictor=list(A=inla.stack.A(stack.est), compute = F),
              verbose=T)
output_gridcode = inla(formula_gridcode,
                       data =inla.stack.data(stack.est),
                       family = "binomial",
                       control.family = list(link="cloglog"),
                       control.predictor=list(A=inla.stack.A(stack.est), compute = F),
                       verbose=T)
save(output, output_gridcode, file = "multiv_1_model_2.RData")

#Model 1
#Removing >30 yo
kfcshh = kfcshh[which(kfcshh$ageyrs<=30),]
A.est = inla.spde.make.A(mesh=mesh,
                         loc=as.matrix(kfcshh[,c("x_proj","y_proj")]))
stack.est = inla.stack(data=list(seroposenr=kfcshh$seroposenr),
                       A=list(A.est,
                              1),
                       effects=list(c(s.index, list(Intercept=1)),
                                    list(kfcshh[,which(names(kfcshh)!= "seroposenr")])),
                       tag='stdata')
output = inla(formula,
              data =inla.stack.data(stack.est),
              family = "binomial",
              control.family = list(link="cloglog"),
              control.predictor=list(A=inla.stack.A(stack.est), compute = F),
              verbose=T)
output_gridcode = inla(formula_gridcode,
                       data =inla.stack.data(stack.est),
                       family = "binomial",
                       control.family = list(link="cloglog"),
                       control.predictor=list(A=inla.stack.A(stack.est), compute = F),
                       verbose=T)
save(output, output_gridcode, file = "multiv_1_model_1.RData")

#Model 3
agecat_vect = seq(5,25,by=5)
agecat = cut(kfcshh$ageyrs,  breaks=1+c(-2,agecat_vect,200), right=F)
kfcshh$age_cat = agecat

formula = seroposenr ~ -1 + Intercept + f(spatial.field, model=spde) + f(houseNo, model = 'iid') +
  #(as.numeric(gridcode)) + 
  #(as.numeric(perCity)) + 
  (as.numeric(toilets_outside)) +
  #(as.numeric(hh_size)) +
  (as.numeric(containers_total)) +
  #(as.numeric(containers_plastic)) +
  #(as.numeric(containers_bottles)) +
  f(houseType.y, model='iid') +
  (as.numeric(garbage)) +
  (as.numeric(gender_bin)) +
  f(occupation, model='iid') +
  f(education_household, model='iid') +
  (as.numeric(not_concrete)) +
  (as.numeric(zinc_roof)) +
  (as.numeric(water_type_pipe)) +
  (as.numeric(waterResourcesNearby.y)) +
  (as.numeric(doorScreens.y)) +
  (as.numeric(devnot)) +
  #(as.numeric(immature_mosquitoes)) +
  #(as.numeric(immature_mosquitoes_strict)) +
  #(as.numeric(aegypti_captured_norm)) +
  f(age_cat, model="rw1", constr = F, extraconstr = list(A=matrix(c(rep(1,length(agecat_vect)),1e4), nrow=1), e=rep(0,1)))
formula_gridcode = seroposenr ~ -1 + Intercept + f(spatial.field, model=spde) + f(houseNo, model = 'iid') +
  (as.numeric(gridcode)) + 
  #(as.numeric(perCity)) + 
  (as.numeric(toilets_outside)) +
  #(as.numeric(hh_size)) +
  (as.numeric(containers_total)) +
  #(as.numeric(containers_plastic)) +
  #(as.numeric(containers_bottles)) +
  f(houseType.y, model='iid') +
  (as.numeric(garbage)) +
  (as.numeric(gender_bin)) +
  f(occupation, model='iid') +
  f(education_household, model='iid') +
  (as.numeric(not_concrete)) +
  (as.numeric(zinc_roof)) +
  (as.numeric(water_type_pipe)) +
  (as.numeric(waterResourcesNearby.y)) +
  (as.numeric(doorScreens.y)) +
  #(as.numeric(devnot)) +
  #(as.numeric(immature_mosquitoes)) +
  #(as.numeric(immature_mosquitoes_strict)) +
  #(as.numeric(aegypti_captured_norm)) +
  f(age_cat, model="rw1", constr = F, extraconstr = list(A=matrix(c(rep(1,length(agecat_vect)),1e4), nrow=1), e=rep(0,1)))

A.est = inla.spde.make.A(mesh=mesh,
                         loc=as.matrix(kfcshh[,c("x_proj","y_proj")]))
stack.est = inla.stack(data=list(seroposenr=kfcshh$seroposenr),
                       A=list(A.est,
                              1),
                       effects=list(c(s.index, list(Intercept=1)),
                                    list(kfcshh[,which(names(kfcshh)!= "seroposenr")])),
                       tag='stdata')
output = inla(formula,
              data =inla.stack.data(stack.est),
              family = "binomial",
              control.family = list(link="logit"),
              control.predictor=list(A=inla.stack.A(stack.est), compute = F),
              verbose=T)
output_gridcode = inla(formula_gridcode,
                       data =inla.stack.data(stack.est),
                       family = "binomial",
                       control.family = list(link="logit"),
                       control.predictor=list(A=inla.stack.A(stack.est), compute = F),
                       verbose=T)
save(output,output_gridcode, file = "multiv_1_model_3.RData")
