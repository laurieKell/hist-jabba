# Load required packages with error checking
pkgs=c("doParallel",  # parallel processing
       "parallel",    # detectCores()
       "stringr",
       "latex2exp",
       "plyr",        # data manipulation
       "dplyr",       # data manipulation
       "purrr",       # data manipulation
       "reshape",     # reshaping data
       "JABBA",       # stock assessment
       "FLCore",      # FLR
       "FLBRP",       # biological reference points
       "FLRebuild",
       "icesdata")    

for(pkg in pkgs) 
  if(!require(pkg, character.only = TRUE)) 
    stop(paste("Package", pkg, "is required but not installed"))

data("spp")
data("info")
data("icesdata")
data("ctc1903")

ts   =tseries(icesdata)
excld=c("aru.27.5b6a","bss.27.4bc7ad-h","meg.27.7b-k8abd","mon.27.8c9a","ple.27.21-23","spr.27.3a4")

ctc=rbind(cbind(What="ICES",transmute(ts,     .id=.id,year=year,data=catch)),
          cbind(What="1903",transmute(ctc1903,.id=.id,year=year,data=catch)))

chk=merge(ddply(ts,.(.id), with, data.frame(ICES=min(year))),
          ddply(ctc1903,.(.id), with, data.frame(Hist=min(year)),by="year"),by=".id")
ids=subset(chk,Hist<=ICES)$.id

ctc=subset(ctc,.id%in%ids)

icesdata=icesdata[names(icesdata)%in%ids]

priors=ldply(icesdata, function(x) 
  tryIt(c(eqsim(x)[c("bmsy","b0"),drop=TRUE],"fmsy"=benchmark(x)["fmsy"])))
priors=ddply(priors, .(.id), with, 
  tryIt(FLRebuild::pellatParams(FLPar(msy=bmsy*(1-exp(-fmsy)),bmsy=bmsy,k=b0))[drop=TRUE]))
priors=merge(priors,ldply(icesdata, function(x) data.frame("psi"=median(ssb(x)[,1:2,drop=TRUE]))))
priors=transform(priors,psi=psi/k,shape=bmsy/k)
priorICES=priors[,c(".id","r","p","k","msy","bmsy","fmsy","virgin","psi","shape")]

dat=transform(priorICES,.id=factor(.id,levels=unlist(priorICES[order(priorICES$psi),".id"])))
ggplot(dat)+
  geom_point(aes(.id,psi))+
  geom_hline(aes(yintercept=1))+
  theme(axis.text.x=element_text(angle=45))

# JABBA assessment function
runHist<-function(id, prior, icesdata=icesdata, ctc1903=ctc1903, quick=TRUE) {    
  if(missing(id) || missing(prior)) 
    stop("id and prior are required parameters")

  # Set up priors and standard deviations
  pr=tryCatch({
    FLPar(subset(prior, .id==id)[,-1])
  }, error = function(e) {
    stop("Failed to create FLPar object: ", e$message)})
  
  # Initialize prior SDs with documented values
  # Default 30% CV for all parameters
  pr.sd             =0.3*pr/pr  
  pr.sd["psi"]      =0.01
  pr.sd["r"]        =0.05    
  pr.sd["shape"]    =0.3   
  pr.sd["k"]        =5.0 

  # Process with reporting
  message(paste("Processing stock:", id))
  
  ts   =tseries(icesdata[id]) 

  catch=subset(ts,.id==id)[, c("year", "catch")]
  small=mean(catch$catch,na.rm=TRUE)*1e-6
  catch[is.na(catch)]=small
  
  c1903=try(subset(ctc1903,.id==id)[, c("year", "catch")])
  c1903$catch[is.na(c1903$catch)|c1903$catch<=0]=small
  maxYear=max(catch$year)
  c1903  =subset(c1903,year<=maxYear)

  eb    =transmute(subset(ts,.id==id), year=year, index=eb)
  eb[is.na(eb)]=NA
  
  eb1903=merge(eb, c1903,by="year",all.y=TRUE)[,seq(dim(eb )[2])]
  eb1903[is.na(eb1903)]=NA
  
  rtn=list()

  message(paste(id,": ICES"))                
  rtn[["ICES"]] =try(FLRebuild:::jabba(catch,pr,pr.sd=pr.sd,index=eb,assessment=id,
                                       q_bounds=c(0.5,2),quick=quick))
  
  pr[   "psi"]=0.9
  pr.sd["psi"]=0.5
  
  message(paste(id,": 1903"))                
  rtn[["1903"]] =try(FLRebuild:::jabba(c1903,pr,pr.sd=pr.sd,index=eb1903,assessment=id,
                                       q_bounds=c(0.5,2),quick=quick))
  return(rtn)
  
  return(rtn)}

runParallel<-function(ids, priors, icesdata, ctc1903, quick=TRUE) {
  message(paste("Starting"))
  
  rtn=foreach(id=ids, 
              .combine=list, 
              .maxcombine=length(ids), 
              .packages=pkgs[-(seq(4))],
              .export=c("runHist"),
              .errorhandling='pass') %dopar% {runHist(id, priors, icesdata, ctc1903, quick=quick)}
  
  names(rtn)=ids
  
  return(rtn)}

# Parallel processing 
library(foreach)
library(doParallel)

ncore=detectCores()-4  
registerDoParallel(ncore)

options(foreachOutputs=TRUE)

histICES=runParallel(ids, priorICES, icesdata, ctc1903, quick=FALSE)

stopImplicitCluster()

jbplot_summary(llply(histICES[ids[34]], function(x) x[["1903"]]$fit))

x=histICES[[6]]
jbplot_summary(list("1903"=x[["1903"]]$fit,
                    "ICES"=x[["ICES"]]$fit))

jbICES=jabbaExtractLists(histICES)
## Priors used for initial depletion, r, ...
prior =jbICES[["priors"]]
## Posteriors, estimates of main parameters and derived quantities
post  =jbICES[["posteriors"]]
## Time series of B/BMSY & F/FMSY
trjc  =jbICES[["trajectory"]]
## Current Status
kb    =ddply(trjc,.(.id), subset, year==max(year))[,c(".id","Scenario","stock","harvest")]

postMed=ddply(post, .(.id,Scenario), with, 
              data.frame(stock=median(stock),
                         psi  =median(psi),
                         r    =median(r),
                         m    =median(m),
                         BmsyK=median(BmsyK)))
trjcMed=ddply(trjc, .(.id,Scenario,year), with,
              data.frame(stock=median(stock,na.rm=TRUE),
                         harvest=median(harvest,na.rm=TRUE),
                         psi  =median(BB0,na.rm=TRUE)))
kbMed=ddply(kb,    .(.id,Scenario), with, 
             data.frame(stock  =median(stock),
                        harvest=median(harvest)))

### B/Bmsy in 1st year of ICES assessment data
minyear=map_dfr(icesdata, ~data.frame(year = dims(.)$minyear), .id=".id")
psi    =rename(merge(trjc, minyear)[,c(".id","Scenario","year","BB0")],c(BB0="psi"))
psiMed =ddply(psi,.(.id,Scenario), with, data.frame(psi=median(psi)))

fn<-function(df) df[sample(seq_len(nrow(df)), size=min(100, nrow(df))), ]
trjc =ddply(trjc,  .(.id,Scenario,year), fn)
post =ddply(post,  .(.id,Scenario),      fn)
kb   =ddply(kb,    .(.id,Scenario),      fn)
psi  =ddply(psi,   .(.id,Scenario),      fn)

save(priorICES,   file="../data/fits/priorICES.RData")
save(prior,       file="../data/fits/prior.RData")
save(trjc,trjcMed,file="../data/fits/trj.RData")
save(post,postMed,file="../data/fits/post.RData")
save(kb,  kbMed,  file="../data/fits/kb.RData")
save(psi, psiMed, file="../data/fits/psi.RData")

save(histICES,file="../data/results/histICES.RData")

bbmsy=ldply(histICES, function(x) ldply(x, function(x) 
  cbind(year=an(dimnames(histICES[[2]][[1]]$fit$timeseries)[[1]]),x[["fit"]]$timeseries[,,"BBmsy"]),
        .id="Scenario"),.id="stock")
ffmsy=ldply(histICES, function(x) ldply(x, function(x) 
  cbind(year=an(dimnames(histICES[[2]][[1]]$fit$timeseries)[[1]]),x[["fit"]]$timeseries[,,"FFmsy"]),
  .id="Scenario"),.id="stock")




