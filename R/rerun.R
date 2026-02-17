library(JABBA)
library(ggplot2)

library(dplyr)
library(knitr)
library(kableExtra)
library(tidyr)
library(purrr)
library(reshape2)

library(coda)

library(foreach)
library(doParallel)

load("C:/active/flrpapers/jabba_eval/data/fits/hist/JABBA.RData")
load("C:/active/flrpapers/rebuild/supMat/data/results/lls.RData")
load("C:/active/flrpapers/RebuildAge/data/eql.RData")
load("C:/active/flrpapers/rebuild/supMat/data/results/pfs.RData")
load("C:/active/flrpapers/hist/data/diags/diagSmry.RData")

## Best shape parameter
best=ddply(merge(pfs,subset(lls,bestfit==6)),.(.id), with, 
           data.frame(shape=eb[which.max(catch)]/eb[which.max(eb)]))


ggplot(best)+geom_histogram(aes(shape))


##### MCMC #####################################################################
rerunH=subset(h,!Heidel_passed)[,1:2]
rerunG=subset(diagSmry,!Geweke_passed)[,1:2]

for (i in seq(dim(rerunH)[1]))
  hist[[rerunH[i,"Stock"]]][[rerunH[i,"Scenario"]]][["fit"]]=
  fit_jabba(hist[[rerunH[i,"Stock"]]][[rerunH[i,"Scenario"]]][["input"]],
            quickmcmc=F,verbose=F)

for (i in seq(dim(rerunH)[1]))
  hist[[rerunH[i,"Stock"]]][[rerunH[i,"Scenario"]]][["fit"]]$pars$Heidel.p

for (i in 15:dim(rerunG)[1]){
  print(i)
  hist[[rerunG[i,"Stock"]]][[rerunG[i,"Scenario"]]][["fit"]]=
    fit_jabba(hist[[rerunG[i,"Stock"]]][[rerunG[i,"Scenario"]]][["input"]],
              quickmcmc=F,verbose=F)}

save(hist, file = "C:/active/flrpapers/hist/data/diags/histMCMC.RData")

## AR failures
rerunAR=subset(diagSmry,!Autocorr_passed&Scenario=="EBiomass 1903")[,1:2]
rerunAR=merge(rerunAR,best,by.x="Stock",by.y=".id")

histARFail=list()
for(i in rerunAR$Stock) 
  histARFail[[i]]=hist[[i]][["EBiomass 1903"]][["fit"]]
names(histARFail)=rerunAR$Stock

histAR=list()
for(i in rerunAR$Stock) 
  histAR[[i]]=hist[[i]][["EBiomass 1903"]]
names(histAR)=rerunAR$Stock

for (i in rerunAR$Stock){
  print(i)
  histAR[[i]][["input"]][["settings"]][["BmsyK"]]=
    subset(rerunAR,Stock==i)$shape}


# Set up cluster
nc=parallel::detectCores() - 4
cl=parallel::makeCluster(nc)
doParallel::registerDoParallel(cl)
parallel::clusterExport(cl, c("histAR"))

# Run fits in parallel, returning a list of results
histAR=foreach(id = rerunAR$Stock,
              .packages = c("JABBA")) %dopar% {
              
               fit_jabba(histAR[[id]][["input"]],quickmcmc=F,verbose=F)}
names(histAR)=rerunAR$Stock

parallel::stopCluster(cl)
registerDoSEQ()

diags     = runAll(histAR)
diagsFail = runAll(histARFail)

calcSmry(diags)
calcSmry(diagsFail)

### OBS Err
rerunAR1=subset(diagSmry,!Autocorr_passed&Scenario=="EBiomass 1903")[,1:2]

for (i in 1:dim(rerunAR1)[1]){
  print(rerunAR1[i,])
  
  jbInp=hist[[rerunAR1[i,"Stock"]]][[rerunAR1[i,"Scenario"]]][["input"]]
 
  suppressMessages(jbinput_new <- build_jabba(
    catch      = jbInp$data$catch,
    cpue       = jbInp$data$cpue,
    se         = jbInp$data$se,
    assessment = jbInp$settings$assessment,
    scenario   = jbInp$settings$scenario,
    model.type = jbInp$settings$model.type,
    fixed.obsE = 0.3,
    projection = jbInp$settings$projection,
    TACs       = jbInp$settings$TACs))
  
  
  hist[[rerunAR1[i,"Stock"]]][[rerunAR1[i,"Scenario"]]][["fit"]]=
    fit_jabba(jbInp,quickmcmc=F,verbose=F)}

save(hist, file = "C:/active/flrpapers/hist/data/diags/histFixObs.RData")

#### PE 
load("C:/active/flrpapers/jabba_eval/data/fits/hist/JABBA.RData")
for (i in 1:dim(rerunAR1)[1]){
  print(rerunAR1[i,])
  
  jbInp=hist[[rerunAR1[i,"Stock"]]][[rerunAR1[i,"Scenario"]]][["input"]]
  
  jbInp$fixed.procE =0.3                 
  jbInp$sigma.proc  =FALSE  
  
  hist[[rerunAR1[i,"Stock"]]][[rerunAR1[i,"Scenario"]]][["fit"]]=
    fit_jabba(jbInp,quickmcmc=F,verbose=F)}

save(hist, file = "C:/active/flrpapers/hist/data/diags/histFixProc.RData")
