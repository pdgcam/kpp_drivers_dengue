library(readxl)
library(INLA)
library(PBSmapping)
library(ggplot2)

for(task_id in 1:70){ #This loop is not optimal - if someone were to run this, I'd recommend pararelising the tasks no a cluster instead
  load("mosquito_data.Rdata")
  load("mosquito_data_revised.Rdata")
  
  kfcshh = read.csv("LandUseKFCSApr2021.csv", header = TRUE, sep = ",", as.is=TRUE, stringsAsFactors = FALSE)
  kfcshh$container_index = container_index_per_ind
  kfcshh$aegypti_captured = aegypti_captured_by_ind
  
  covar_vect = c("popdens","percity","toilets","hhsize","ncontainers","ncomtainersplastic","ncontainersbottles","housetype","garbagemanagement","gender","occupation",
                 "hheducation","constructionmaterial","rooftype","watertype","waterressources","doorscreens","devnot","container_index","immature_strict","aegypti_captured","baseline","plastic_burn")
  
  covar_id = rep(1:length(covar_vect),each=3)[task_id]
  model_id = rep(1:3,length(covar_vect))[task_id] #Model 1 : baseline / Model 2 : including >30 yo / Model 3 : removing age offset
  
  
  #Removing newborns
  kfcshh = kfcshh[which(kfcshh$ageyrs>=1),]
  
  if(model_id !=2){
    #Removing >30 yo
    kfcshh = kfcshh[which(kfcshh$ageyrs<=30),]
  }
  
  #Adding meter coordinates
  coords_home = kfcshh[,c("gpslong","lat")]
  colnames(coords_home)<-c("X","Y")
  attr(coords_home,"projection")<-"LL"
  coords_home<-convUL(coords_home,km=F)
  kfcshh[,c("x_proj","y_proj")] = coords_home
  
  #Remove the .x columns
  duplicated =  (1:ncol(kfcshh)) %in% grep("\\.x", names(kfcshh))
  kfcshh = kfcshh[,!duplicated]
  
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
  
  #Run univariates
  if(covar_id == 1){#Pop density
    kfcshh$gridcode = kfcshh$gridcode/1000*100 #Make it from per 1hab/100m^2 to per 1000hab/km^2
    add_cov = "(as.numeric(gridcode))"} 
  if(covar_id == 2){#perCity
    add_cov = "(as.numeric(PerCity))"}
  if(covar_id == 3){#Toilets
    kfcshh$toilets_outside = kfcshh$numToiletsOutside.y>0
    add_cov = "(as.numeric(toilets_outside))"}
  if(covar_id == 4){#Hh_size
    #table(cut(kfcshh$numSubjects, breaks = c(0,3,6,13)), kfcshh$numSubjects)
    # kfcshh$hh_size = cut(kfcshh$numSubjects, breaks = c(0,3,6,13))
    # add_cov = "f(hh_size, model='rw1')"
    add_cov = "(as.numeric(numSubjects))"
  }
  if(covar_id == 5){
    #Nwatercontainers
    cov_watercontainers = kfcshh[,grep("waterContainer", names(kfcshh))]
    containers_total = rowSums(cov_watercontainers[,-1], na.rm=T)
    containers_plastic = cov_watercontainers$waterContainerSoftDrink.y + cov_watercontainers$waterContainerGallons.y
    #kfcshh$containers_total = cut(containers_total, breaks = c(seq(0,50,by=10),250))
    #table(kfcshh$containers_total, containers_total)
    #add_cov = "f(containers_total, model='rw1')"
    
    kfcshh$containers_total = containers_total
    add_cov = "(as.numeric(containers_total))"
  }
  
  if(covar_id == 6){
    #Nwatercontainers - plastic
    cov_watercontainers = kfcshh[,grep("waterContainer", names(kfcshh))]
    containers_plastic = cov_watercontainers$waterContainerSoftDrink.y + cov_watercontainers$waterContainerGallons.y
    # kfcshh$containers_plastic = cut(containers_plastic, breaks = c(-1,0,10,20,200))
    # #table(kfcshh$containers_plastic, containers_plastic)
    # add_cov = "f(containers_plastic, model='rw1')"
    kfcshh$containers_plastic = containers_plastic
    add_cov = "(as.numeric(containers_plastic))"
  }
  if(covar_id == 7){#Nwatercontainers - plastic bottles
    kfcshh$containers_bottles = kfcshh$waterContainerSoftDrink.y>0
    add_cov = "(as.numeric(containers_bottles))"
  }
  if(covar_id == 8){#Housetype
    add_cov = "f(houseType.y, model='iid', constr = F, extraconstr = list(A=matrix(c(1e4,rep(1,length(unique(kfcshh$houseType.y))-1)), nrow=1), e=rep(0,1)))"
  }
  if(covar_id == 9){
    #Garbage management
    garbage = as.factor(kfcshh$garbageManagement.y)
    levels(garbage) = c("Burn/Bury/Dump","Burn/Bury/Dump","Car","Burn/Bury/Dump")
    kfcshh$garbage = garbage!="Car"
    add_cov = "(as.numeric(garbage))"
  }
  if(covar_id == 10){
    #Gender
    kfcshh$gender_bin = kfcshh$gender != "F"
    add_cov = "(as.numeric(gender_bin))"
  }
  if(covar_id == 11){#Occupation
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
    add_cov = "f(occupation, model='iid', constr = F, extraconstr = list(A=matrix(c(1e4,rep(1,length(unique(kfcshh$occupation))-1)), nrow=1), e=rep(0,1)))"
  }
  if(covar_id == 12){#hh_education
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
    add_cov = "f(education_household, model='iid', constr = F, extraconstr = list(A=matrix(c(1e4,rep(1,length(unique(kfcshh$education_household))-1)), nrow=1), e=rep(0,1)))"
  }
  if(covar_id == 13){#construction material
    kfcshh$not_concrete = kfcshh$constructionMaterialsConcrete.y == 0
    add_cov = "(as.numeric(not_concrete))"
  }
  if(covar_id == 14){
    #roof type
    kfcshh$zinc_roof = kfcshh$roofTypeZinc.y>0
    add_cov = "(as.numeric(zinc_roof))"
  }
  if(covar_id == 15){#water type
    kfcshh$water_type_pipe = kfcshh$waterTypePipe.y>0
    add_cov = "(as.numeric(water_type_pipe))"
  }
  if(covar_id == 16){#water ressources
    add_cov = "(as.numeric(waterResourcesNearby.y))"
  }
  if(covar_id == 17){#door screens
    add_cov = "(as.numeric(doorScreens.y))"
    
  }
  if(covar_id == 18){#devnot
    add_cov = "(as.numeric(devnot))"}
  if(covar_id == 19){#immature mosqu
    kfcshh = kfcshh[!is.na(kfcshh$container_index),]
    add_cov = "(as.numeric(container_index))"}
  if(covar_id == 20){#immature mosqu srict
    add_cov = "(as.numeric(immature_mosquitoes_strict))"}
  if(covar_id == 21){#nb aegypti captured
    kfcshh = kfcshh[!is.na(kfcshh$aegypti_captured),]
    add_cov = "(as.numeric(aegypti_captured))"}
  if(covar_id==23){
    kfcshh$plastic_litter = as.factor((kfcshh$garbageManagement.y=="Burn" & kfcshh$waterContainerSoftDrink.y>0 & !kfcshh$devnot) +
                                        (kfcshh$garbageManagement.y=="Burn" & kfcshh$waterContainerSoftDrink.y>0))
    levels(kfcshh$plastic_litter) = c("No burnt plastic", "Burnt plastic in city", "Rural burnt plastic")
    add_cov = "f(plastic_litter, model='iid', constr = F, extraconstr = list(A=matrix(c(1e4,rep(1,length(unique(kfcshh$plastic_litter))-1)), nrow=1), e=rep(0,1)))"
  }
  
  #Create spatial field
  mesh_res = 4000
  locs = as.matrix(kfcshh[,c("x_proj","y_proj")])
  mesh = inla.mesh.2d(loc=kfcshh[,c("x_proj","y_proj")], max.edge= c(mesh_res,mesh_res*5),
                      cutoff=mesh_res)
  spde = inla.spde2.matern(mesh=mesh,alpha=2)
  s.index <- inla.spde.make.index(name="spatial.field", n.spde=spde$n.spde) 
  A.est = inla.spde.make.A(mesh=mesh,
                           loc=as.matrix(kfcshh[,c("x_proj","y_proj")]))
  
  #Model 1 & 2
  if(model_id!=3){
    if(covar_id == 22){formula = seroposenr ~ -1 + Intercept + offset(log(ageyrs)) + f(houseNo, model = 'iid') + f(spatial.field, model=spde)} #Baseline model
    else{formula = as.formula(paste0("seroposenr ~ -1 + Intercept + offset(log(ageyrs)) + f(houseNo, model = 'iid') + f(spatial.field, model=spde) +", add_cov))}
    
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
    
  }
  
  
  #Model 3
  if(model_id==3){
    
    if(covar_id == 22){formula = seroposenr ~ -1 +Intercept + f(spatial.field, model=spde) + f(houseNo, model = 'iid')} #Baseline model
    else{formula = as.formula(paste0("seroposenr ~ -1 +Intercept + f(spatial.field, model=spde) + f(houseNo, model = 'iid') + ", add_cov))}
    
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
  }
  
  save(output, task_id, covar_id, model_id, file = paste0(task_id,"_univariate_",covar_vect[covar_id],"_model_",model_id,".RData"))
}