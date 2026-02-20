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


pt<-function(r,p,k,b=k*seq(0,1,length.out=101)) 
  data.frame(biomass=b,
             production=(r / p) * b * (1 - (b / k)^p))

ggplot(pt(0.16,1.38,100000))+
  geom_path(aes(biomass,production))

dat=ddply(priorICES,.(.id), with, pt(r,p,k))
ggplot(dat)+
  geom_path(aes(biomass,production))+
  facet_wrap(~.id,scale="free")+
  theme_minimal(base_size=12)+
  theme(
    legend.position="bottom",
    legend.title=element_text(face="bold"),
    panel.grid.major.x=element_line(color="gray90", linewidth=0.2),
    panel.grid.minor=element_blank(),
    strip.text=element_text(face="bold", size=10, color="gray30"),
    plot.title=element_text(face="bold", size=14, margin=margin(b=10)),
    plot.subtitle=element_text(size=10, color="gray40", margin=margin(b=15)),
    axis.title.x=element_text(size=12, margin=margin(t=10)),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    panel.spacing=unit(0, "lines"),
    plot.margin=margin(0,0,0,0))

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
  pr.sd["psi"]      =0.5
  pr.sd["r"]        =0.5    
  pr.sd["shape"]    =0.5   
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
                                       q_bounds   =c(0.5,2),quick=quick,
                                       sigma.est  =FALSE,
                                       fixed.obsE =0.15, 
                                       sigma.proc =TRUE,
                                       fixed.procE=0.3,
                                       igamma     =c(0.001,0.001)))
  
  pr[   "psi"]=0.9
  pr.sd["psi"]=0.5
  
  message(paste(id,": 1903"))                
  rtn[["1903"]] =try(jabba(c1903,pr,pr.sd=pr.sd,index=eb1903,assessment=id,
                                       q_bounds=c(0.5,2),quick=quick,
                                       sigma.est  =FALSE,
                                       fixed.obsE =0.15, 
                                       sigma.proc =TRUE,
                                       fixed.procE=0.3,
                                       igamma     =c(0.001,0.001)))
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

save(histICES,file="data/results/histICES.RData")

jbplot_summary(llply(histICES[ids[34]], function(x) x[["1903"]]$fit))

x=histICES[[6]]
jbplot_summary(list("1903"=x[["1903"]]$fit,
                    "ICES"=x[["ICES"]]$fit))

jbICES=FLRebuild:::jabbaExtractLists(histICES)
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

save(priorICES,   file="data/fits/priorICES.RData")
save(prior,       file="data/fits/prior.RData")
save(trjc,trjcMed,file="data/fits/trj.RData")
save(post,postMed,file="data/fits/post.RData")
save(kb,  kbMed,  file="data/fits/kb.RData")
save(psi, psiMed, file="data/fits/psi.RData")

bbmsy=ldply(histICES, function(x) ldply(x, function(x) 
  cbind(year=an(dimnames(x$fit$timeseries)[[1]]),x[["fit"]]$timeseries[,,"BBmsy"]),
        .id="Scenario"),.id="stock")
ffmsy=ldply(histICES, function(x) ldply(x, function(x) 
  cbind(year=an(dimnames(x$fit$timeseries)[[1]]),x[["fit"]]$timeseries[,,"FFmsy"]),
  .id="Scenario"),.id="stock")
save(bbmsy,ffmsy, file="data/fits/msy.RData")



