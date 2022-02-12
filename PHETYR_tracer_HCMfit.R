graphics.off()
rm(list=ls())
install.packages("readxl")
install.packages("minpack.lm")
library(readxl)
library(minpack.lm)

source("functions/getoutlier.R")
source("functions/initial_value_constrained.R")
source("functions/expfit_MS.R")

##################################################
##### importing the data 
##################################################
#PHE6
PHE6 = read_excel("data/15-OCERA -SUM-Pulse-TTR.xlsx",
                  sheet = "PHE6",
                  col_names=FALSE
)
subject_day= paste(PHE6[2,-1], PHE6[1,-1], sep="")

colnames(PHE6)= c("Time", subject_day)
PHE6.m= data.matrix(PHE6[-c(1,2), ])
#rownames(PHE6.m)=NULL
#rm(subject_day)
PHE6data= as.data.frame(PHE6.m)
time= PHE6data[,"Time"]
rm(subject_day, PHE6.m)

#TYR4
TYR4 = read_excel("data/15-OCERA -SUM-Pulse-TTR.xlsx",
                  sheet = "TYR4",
                  col_names=FALSE
)
subject_day= paste(TYR4[2,-1], TYR4[1,-1], sep="")

colnames(TYR4)= c("Time", subject_day)
TYR4.m= data.matrix(TYR4[-c(1,2), ])
#rownames(TYR4.m)=NULL
#rm(subject_day)
TYR4data= as.data.frame(TYR4.m)
rm(subject_day, TYR4.m)

#TYR6
TYR6 = read_excel("data/15-OCERA -SUM-Pulse-TTR.xlsx",
                  sheet = "TYR6",
                  col_names=FALSE
)
subject_day= paste(TYR6[2,-1], TYR6[1,-1], sep="")

colnames(TYR6)= c("Time", subject_day)
TYR6.m= data.matrix(TYR6[-c(1,2), ])
TYR6data= as.data.frame(TYR6.m)
rm(subject_day, TYR6.m)

##################################################
# A preliminary fit on PHE6 and TYR4 
##################################################
# This part detects and removes the outliers, and 
#  estimates the initial parameter values

f_2exp= function(pars2, t){
  expr = expression(pars2["P1"]* exp(-pars2["p1"]*t)+
                      pars2["P2"]*exp(-pars2["p2"]*t)#+
                    #pars2["P"]
  )
  eval(expr)
}

typ= c(NA, rep("pulse", ncol(PHE6data)-1))
phe6_initial= initial_value_const(data_df = PHE6data, typ = typ)
phe6_fit= expfit_MS(d = PHE6data, start_par = phe6_initial, 
                    rob_method = "lorentz", typ = typ, 
                    constraint = T, shared_pars = NULL)

typ= c(NA, rep("pulse", ncol(TYR4data)-1))
tyr4_initial= initial_value_const(data_df = TYR4data, typ = typ)
tyr4_fit= expfit_MS(d = TYR4data, start_par = tyr4_initial, 
                    rob_method = "lorentz", typ = typ, 
                    constraint = T, shared_pars = NULL)

##################################################
##### Model's data and indicators 
##################################################
# This part creates the data for HCM model along with indicators.
subjectid= c(rep(colnames(PHE6data)[-1], each= nrow(PHE6data)-1),    
             rep(colnames(TYR4data)[-1], each= nrow(TYR4data)-1),
             rep(colnames(TYR6data)[-1], each= nrow(TYR6data)))  

sampletruetime= c(rep(PHE6data[-1,"Time"], ncol(PHE6data)-1),
                  rep(TYR4data[-1,"Time"], ncol(TYR4data)-1),
                  rep(TYR6data[,"Time"], ncol(TYR6data)-1))
ttrphe6= unlist(PHE6data[-1,-1])
ttrtyr4= unlist(TYR4data[-1,-1])
ttrtyr6= unlist(TYR6data[,-1])

is.phe6= c(
  rep(1, length(colnames(PHE6data)[-1])*(nrow(PHE6data)-1)),
  rep(0, length(colnames(TYR4data)[-1])*(nrow(TYR4data)-1)),
  rep(0, length(colnames(TYR6data)[-1])*(nrow(TYR6data)))
)

is.tyr4= c(
  rep(0, length(colnames(PHE6data)[-1])*(nrow(PHE6data)-1)),
  rep(1, length(colnames(TYR4data)[-1])*(nrow(TYR4data)-1)),
  rep(0, length(colnames(TYR6data)[-1])*(nrow(TYR6data)))
)

is.tyr6= c(
  rep(0, length(colnames(PHE6data)[-1])*(nrow(PHE6data)-1)),
  rep(0, length(colnames(TYR4data)[-1])*(nrow(TYR4data)-1)),
  rep(1, length(colnames(TYR6data)[-1])*(nrow(TYR6data)))
)


ddf=data.frame( subjectid= subjectid, 
                sampletruetime= sampletruetime, 
                isotope=c(ttrphe6, ttrtyr4, ttrtyr6), 
                is.phe6= is.phe6,
                is.tyr4= is.tyr4,
                is.tyr6= is.tyr6)

dfac= as.factor(subjectid)
sublev= levels(dfac)
sublev_mat= as.matrix(sublev, ncol= 1)   # create one column containing  unique(subjectids)

# store data in a list where each element is a dataframe for one subject 
dlist= apply (sublev_mat, MARGIN = 1, function(lev)
{
  indx= which(ddf[,"subjectid"] == lev)
  ddf[indx, ]
})
names(dlist)= sublev

##################################################
# Defining model's functions (equations 3 in the manuscript) 
##################################################
phe6expr= expression( pars["A"]*exp(-pars["lambda1"]*t) + pars["B"]*exp(-pars["lambda2"]*t))
tyr4expr= expression( pars["C"]*exp(-pars["lambda3"]*t) + pars["D"]*exp(-pars["lambda4"]*t))
tyr6expr= expression ( pars["k31"]*(
  (pars["A"]*pars["C"]/(pars["lambda3"]-pars["lambda1"]) + 
     pars["A"]*pars["D"]/(pars["lambda4"]-pars["lambda1"])) * exp(-pars["lambda1"]*t) +
    (pars["B"]*pars["C"]/(pars["lambda3"]-pars["lambda2"]) + 
       pars["B"]*pars["D"]/(pars["lambda4"]-pars["lambda2"])) * exp(-pars["lambda2"]*t) - 
    (pars["B"]*pars["C"]/(pars["lambda3"]-pars["lambda2"]) + 
       pars["A"]*pars["C"]/(pars["lambda3"]-pars["lambda1"])) * exp(-pars["lambda3"]*t) - 
    (pars["B"]*pars["D"]/(pars["lambda4"]-pars["lambda2"]) + 
       pars["A"]*pars["D"]/(pars["lambda4"]-pars["lambda1"])) * exp(-pars["lambda4"]*t)
)
)

##################################################
# Defining model's function and initial parameters' values
  # for the simultaneous fit on the equations 3 of the manuscript
##################################################
f_phetyr= function (dd, pars, t) {
  dd$is.phe6 * eval(phe6expr) + dd$is.tyr4 * eval(tyr4expr) + dd$is.tyr6 * eval(tyr6expr)
} 

model_initial= vector(mode = "list", length = length(dlist))
names(model_initial)= names(dlist)

for (i in 1:length(model_initial)) {
  #sbj= names(model_initial)[i]
  A= phe6_fit$all.f_2[[i+1]]$par[["P1"]]
  lambda1= phe6_fit$all.f_2[[i+1]]$par[["p1"]]
  B= phe6_fit$all.f_2[[i+1]]$par[["P2"]]
  lambda2= phe6_fit$all.f_2[[i+1]]$par[["p2"]]
  
  C= tyr4_fit$all.f_2[[i+1]]$par[["P1"]]
  lambda3= tyr4_fit$all.f_2[[i+1]]$par[["p1"]]
  D= tyr4_fit$all.f_2[[i+1]]$par[["P2"]]
  lambda4= tyr4_fit$all.f_2[[i+1]]$par[["p2"]]
  
  k31= 0.005  # equivalent to the k_phe-->tyr
              # approximated from a non-compartmental approach 
              # could also be a random initial value between 0 and 1 
  
  model_initial[[i]]= c(A= A, lambda1= lambda1, B= B, lambda2= lambda2,
                        C= C, lambda3= lambda3, D= D, lambda4= lambda4, 
                        k31= k31)
} 

##############################
# fitting the model's function on the model's data per subject
##############################
for (ii in 1:length(dlist)) {
  dd= dlist[[ii]]
  print(dd[1,1])
  
  indx= which(!is.na(dd[,"isotope"]))  
  dnna= dd[indx,]
  
  # defining the model's residuals for a non-linear least square regression
  res_ord= function(true_value, pars, t) {
    (true_value - f_phetyr(dnna, pars, t)) / diff(range(true_value))
  }
  model_ord= nls.lm(par = model_initial[[ii]], 
                    lower = rep(1e-10, 9), 
                    upper = rep(1e10, 9),
                    fn = res_ord, 
                    true_value= dnna[,"isotope"],
                    t= dnna[,"sampletruetime"],
                    control = nls.lm.control(nprint=1, 
                                             maxiter = 500, 
                                             ftol = 1e-15, 
                                             ptol = 1e-15, 
                                             gtol = 1e-15)
  )
  
  indx_phe6= dd[,"is.phe6"]==1
  indx_tyr4= dd[,"is.tyr4"]==1
  indx_tyr6= dd[,"is.tyr6"]==1
  
 # plotting fits per subject 
  par(mfrow=c(3,1))
  plot(dd[indx_phe6,"sampletruetime"], dd[indx_phe6,"isotope"],
       xlab= "Time(min)", ylab = "phe6")
  lines(dd[indx_phe6,"sampletruetime"], 
        f_phetyr(dd[indx_phe6, ], pars = model_ord$par, t = dd[indx_phe6,"sampletruetime"]), 
        col="red")
  title(dd[1,1])
  
  plot(dd[indx_tyr4,"sampletruetime"], dd[indx_tyr4,"isotope"], 
       xlab= "Time(min)", ylab = "tyr4")
  lines(dd[indx_tyr4,"sampletruetime"], 
        f_phetyr(dd[indx_tyr4, ], pars = model_ord$par, t = dd[indx_tyr4,"sampletruetime"]), 
        col="red")
  
  plot(dd[indx_tyr6,"sampletruetime"], dd[indx_tyr6,"isotope"], 
       xlab= "Time(min)", ylab = "tyr6")
  lines(dd[indx_tyr6,"sampletruetime"], 
        f_phetyr(dd[indx_tyr6, ], pars = model_ord$par, t = dd[indx_tyr6,"sampletruetime"]), 
        col="red")
}
