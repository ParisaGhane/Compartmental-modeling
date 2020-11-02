#first column in data_df must be Time
## func should be a function with arguments par= initial parameters and t= time --> f(par, t)
###Type can be "metabolite" or "pulse" (default is pulse)
#### rob_method can be "lad" or "lorentz" (default is lad)

# function uses nonlinear least square with Levenberge-Marquardt algorithm (function nls.lm)

#Parisa Ghane (pghane@tamu.edu)

###################################################################################################

expfit_MS = function (d, start_par, rob_method= "lad", typ, constraint, shared_pars= NULL, lower=NULL, upper=NULL)  
{
  # install.packages("minpack")
  ## Initializaiton of outputs #####
  data_original= d
  
  indx_l = vector(mode = "list", length = (ncol(data_original)))
  new_data_l = NULL
  new_data_l= vector(mode = "list", length = (ncol(data_original)))
  all.model_fitted= vector(mode = "list", length = (ncol(data_original)))
  fit_est= vector(mode = "list", length = (ncol(data_original)))
  AUC_l= AUC_l_2exp= vector(mode = "list", length = (ncol(data_original)))
  all.f_1= vector(mode = "list", length = (ncol(data_original)))
  all.f_2= vector(mode = "list", length = (ncol(data_original)))
  all.f_3= vector(mode = "list", length = (ncol(data_original)))

  pars1= start_par[["pars1"]]
  pars2= start_par[["pars2"]]
  pars3= start_par[["pars3"]]
  # if (is.null(start)) next
  ##define functions ########
  f_1exp= function(pars1, t){
    expr = expression( (pars1["P1"]* exp(-pars1["p1"]*t)
                        # pars1["P"]
    ))
    eval(expr)
  }
  
  f_2exp= function(pars2, t){
    expr = expression( ( pars2["P1"]* exp(-pars2["p1"]*t)+
                        pars2["P2"]*exp(-pars2["p2"]*t)
                        # pars2["P"]
    ))
    eval(expr)
  }
  
  f_3exp= function(pars3, t){
    expr = expression((pars3["P1"]* exp(-pars3["p1"]*t)+
                        pars3["P2"]*exp(-pars3["p2"]*t)+
                        pars3["P3"]*exp(-pars3["p3"]*t)
                        # pars3["P"]
    ))
    eval(expr)
  }
  
  ## functions fitting procedure#####
  fit= function (data, func, start, shared_const= NULL) {
    f=func
    ##### boundary constraints #####
    if (is.null(lower) | is.null(upper)) {
      if(constraint) {
        if(is.null(shared_const)) {
          if (type == "pulse") {
            # lower = c(rep(1e-3, (length(start))-1), 1e-8)
            # upper= c(rep(1, (length(start))))
            lower = c(rep(1e-6, (length(start))-1), 1e-6) #1e-10
            upper= c(rep(10, (length(start)))) #10
                                      # upper[c(F,T)]= upper[c(F,T)]*0.02
          } else if (type == "metabolite") {
            lower= c(-2e4, rep(1e-6, (length(start)-2)), 1e-8)
            upper= c(-1e-6, rep(1e4, (length(start)-2)), 1e2)
          }
        }else if (!is.null(shared_const)) {
          #aa= which(rownames(shared_const) == names(type)[i])
          shared_const_aa= shared_const[names(type), ]
          
          if (length(shared_const_aa) %in% 4:5) {
            if (type == "pulse") {
              lower= c(1e-4, unname(shared_const_aa["p1"])-1e-7,
                       rep(1e-4, (length(start)-2)))
              upper= c(1e3, unname(shared_const_aa["p1"])+1e-7, 
                       rep(c(1e3, unname(shared_const_aa["p1"])), (length(start)-3)/2),
                       1e3)
            }else if (type == "metabolite") {
              lower= c(-1e4, unname(shared_const_aa["p1"])-1e-7,
                       1e-6, unname(shared_const_aa["p2"])-1e-7,
                       0)
              upper= c(-1e-7, unname(shared_const_aa["p1"])+1e-7, 
                       1e3, unname(shared_const_aa["p2"])+1e-7, 
                       1e3)
            }
          }
          
            if (length(shared_const_aa) %in% 6:7) {
              if(type == "pulse") {
                lower= c(1e-5, unname(shared_const_aa["p1"])-1e-7,
                         1e-5, unname(shared_const_aa["p2"])-1e-7,
                         1e-5, unname(shared_const_aa["p3"])-1e-7,
                         -1e-7)
                upper= c(1e3, unname(shared_const_aa["p1"])+1e-7, 
                         1e3, unname(shared_const_aa["p2"])+1e-7, 
                         1e3, unname(shared_const_aa["p3"])+1e-7,
                         1e3)
            } else if (type == "metabolite") {
              lower= c(-1e4, unname(shared_const_aa["p1"])-1e-7,
                       1e-6, unname(shared_const_aa["p2"])-1e-7,
                       1e-6, unname(shared_const_aa["p3"])-1e-7,
                       0)
              upper= c(-1e-7, unname(shared_const_aa["p1"])+1e-7, 
                       1e4, unname(shared_const_aa["p2"])+1e-7, 
                       1e4, unname(shared_const_aa["p3"])+1e-7,
                       1e4)
              
            }
            }
      } 
        } else {
        lower= NULL
        upper= NULL
        }
      }
    
    ##### Ordinary fit ########## 
    res_ord= function(true_value, par, t) {
      # nn= round(length(true_value)/3)
      # if(diff(range(true_value)) < 1) {
      #   ww= 1/diff(range(true_value))
      # } else { ww= diff(range(true_value))}
      # ww[1:(nn-1)]= 0.6 * ww[1:(nn-1)]
      # ww[(2*nn):(length(true_value))]= 0.3*ww[(2*nn):(length(true_value))]
      # r= (true_value - f(par, t)) * ww  #/ diff(range(true_value))) #* (c(0.5, rep(1, length(true_value)-1))) 
      r= (true_value - f(par, t))   / diff(range(true_value))
      r
      }
    
    model_ord= nls.lm(par = start, 
                      lower = lower, 
                      upper = upper,
                      fn = res_ord, 
                      true_value= data[,"y"],
                      t= data[ , "t"],
                      control = nls.lm.control(nprint=-1, 
                                               maxiter = 500, 
                                               ftol = 1e-25, 
                                               ptol = 1e-25, 
                                               gtol = 1e-25)
    )
    
    ##### LAD (or LORENTZ) robust fit with New initial parameters == fitted parameters from ordinary fit ################     
    new_par= coef(model_ord)

    if (rob_method=="lad") {
      res_rob= function(true_value, par, t) {
        (sqrt(abs(true_value - f(par, t))))
      } 
    } else if (rob_method=="lorentz") {
      P.68= 68.27
      K= length(new_par)
      N= length(data[,"y"])
      RSDR= P.68*(N)/(N-K)
      res_rob= function(true_value, par, t) {
        D= true_value - f(par, t)
        tmp= log (1+ (D/RSDR)^2)
        sqrt(tmp)
      }
    }
    model_rob= nls.lm(par = new_par, 
                      fn = res_rob, 
                      lower = lower, 
                      upper = upper,
                      true_value= data[,"y"],
                      t= data[ , "t"],
                      control = nls.lm.control(nprint=-1, 
                                               maxiter = 500, 
                                               ftol = 1e-25, 
                                               ptol = 1e-25, 
                                               gtol = 1e-25)
    )

        #### get outlier from robust approach ########
        res= resid(model_rob)
        std_rob= sqrt(var(res))
        m= mean(res)
        thr1= m - 1.96* std_rob
        thr2= m + 1.96*std_rob

        if ((thr2-thr1) > 0) {
          indx= which(res<thr1 | res > thr2)

          if (length(indx) == 0) 
          {
            new_data=data
            std_new= std_rob
            indx_l[[i]]= integer(0)
          }else {
            new_data= data.frame(y= data[,"y"][-c(indx)], t= data[,"t"][-c(indx)])
            std_new= std_rob[-indx]
            indx_l[[i]]= c(indx)
          }
          indx_l[[i]]= c(indx)
          rm(indx)
        } else {
          new_data= data
          indx_l[[i]]= NULL
        }
    ########################
#     new_data= data
    #########################
    
    
    
    
    new_data_l[[i]]= new_data
    names(new_data_l [[i]])= names(type)
    ##### Fit again with ordinary approach and new_data (data without outliers) and   ###########
    ### initial parameters == parameters produced in first ord fit and
    ### for pulse: weights= 1/y --> wres= res/y_hat
    ### for metabolite: weight = 1 --> wres= res
    
    # w= max(new_data[ ,"y"]) - min(new_data[ ,"y"])
    wres_ord_new= function (true_value, par, t) {
      (true_value - f(par, t)) / sd(true_value, na.rm = T)
      # (true_value - log(f(par, t), base = 10)) / diff(range(true_value))
    }
    
    # wres_ord_new= if (type == "metabolite") {
    #   function (true_value, par, t) {
    #     true_value - f(par, t)
    #   }
    # } else {
    #   function(true_value, par, t) {
    #     (true_value - f(par, t))
    #   }
    # } 
    
    
    model_nout= nls.lm(par = new_par,
                       lower = lower, 
                       upper = upper,
                       fn = wres_ord_new, 
                       true_value= (new_data[,"y"]),
                       t= new_data[ , "t"],
                       control = nls.lm.control(nprint=-1, 
                                                maxiter = 500, 
                                                ftol = 1e-50, 
                                                ptol = 1e-50, 
                                                gtol = 1e-50)
    )
    
    #cnout= coef(model_nout)
    #rss= resid (model_nout)
    model_nout
    
  }
  # install.packages("minpack.lm")
  library(minpack.lm)
  
  #f= func
  
  #### preprocessing ########
  if (is.list(data_original)) {
    data_df= matrix(unlist(data_original), 
                    ncol = ncol(data_original))
    colnames(data_df)= colnames(data_original)
    rownames(data_df)= rownames(data_original)
  } else {
    data_df= matrix(as.numeric(data_original), 
                    nrow = nrow(data_original), 
                    ncol= ncol(data_original))
    rownames(data_df)= rownames(data_original)
  }
 ######### 
  for (i in 2:ncol(data_df))                                  
  {
    #### preprocessing per amino acid ########
    print(i)
    type= typ[i]
    
    if (is.na(type)) next
    if (type != "metabolite" & type != "pulse") next      # skip columns if type is not metabolite or pulse 
    
    print(type)
    if (type == "pulse") {
      data= data.frame(y= as.numeric(data_df[,i]), #was data.frame(y= as.numeric(data_df[,i][-1])
                       t= as.numeric(as.character(data_df[, 1]))) #was as.numeric(as.character(data_df[, 1][-1])))
      rownames(data)= rownames(data_df) #was rownames(data_df)[-1]
      
    } else if (type == "metabolite") {
      data= data.frame(y= as.numeric(data_df[,i]),
                       t= as.numeric(as.character(data_df[, 1]))+1e-7)
      rownames(data)= rownames(data_df)
      
    }
    #data[,"t"] = data[,"t"]-min(data[,"t"]) + 1e-7
    

    
    if (length(which(is.na(data[,"y"]))) == 0) {                     # remove NA values
      data= data
    }else {
      ind= which(is.na(data[,"y"]))
      data_name= rownames(data)
      na_name= data_name[ind]
      data= data.frame(y= as.numeric(as.character(data_df[(setdiff(data_name, na_name)) , i])), 
                       t= as.numeric(as.character(data_df[(setdiff(data_name, na_name)) , 1]))
      )
    } 
    
    new_data_l[[i]]= data
    #### fit conditioned on l and pick a model according to F_test####
    l= length(data[,"y"])
    if (l %in% 0:2) next
    if (l %/% 3 >=3) {  # is true if l>9  # 2_exp & 3_exp
      #shared_const= shared_pars$shared_pars_3exp
      f_2= fit(data, func= f_2exp, start= pars2[[i]],
               shared_const= shared_pars$shared_pars_2exp)
      rss2= f_2$deviance
      p2=length(pars2[[i]])
      all.f_2[[i]]= f_2
      
      if(l %/% 3 >3) { 
      f_3= fit(data, func= f_3exp, start= pars3[[i]], 
               shared_const= shared_pars$shared_pars_3exp)
      rss3= f_3$deviance
      p3= length(pars3[[i]])
      all.f_3[[i]]= f_3
      
      # rm(shared_const)
      # shared_const= shared_pars$shared_pars_2exp
      
      f_test= ((rss2-rss3)/(p3-p2)) / (rss3/(l-p3))
      critical_val= qf(0.95, p3-p2, l-p3)     #actual value of F-distribution with 
      
      if (f_test >= critical_val) {
        f_fitted= f_3
        fint=function(t) {f_3exp(f_3$par, t)}
        
        fitted_est= f_3exp(pars3 = f_fitted$par, t= data[,"t"])
      } else {
        f_fitted= f_2
        fint=function(t) {f_2exp(f_2$par, t)}
        
        fitted_est= f_2exp(pars2 = f_fitted$par, t= data[,"t"])
        }
      }
      f_fitted= f_2
      fint=function(t) {f_2exp(f_2$par, t)}
      
      fitted_est= f_2exp(pars2 = f_fitted$par, t= data[,"t"])
      }
    
    
    if ( l%/%3 < 3 & l%/% 2 >= 2 & l != 4 ) { # is true if l= 5:8  # 2_exp & 1_exp
      # rm(shared_const)
      # shared_const= shared_pars$shared_pars_2exp
      f_2= fit(data, func= f_2exp, start= pars2[[i]], 
               shared_const= shared_pars$shared_pars_2exp)
      rss2= f_2$deviance
      p2=length(pars2[[i]])
      all.f_2[[i]]= f_2
      
      dl= log(data[,"y"]+ 1e-7)
      dd= data.frame(dl=dl, t= data[,"t"])
      lin_f_1= lm(dl~t, dd)
      pars1= c(P1= unname(exp(lin_f_1$coefficients["(Intercept)"])), 
               p1= -unname(lin_f_1$coefficients["t"]), 
               P= -1e-7)
      # rm(shared_const)
      # shared_const= NULL
      f_1= fit(data, func = f_1exp, start= pars1, shared_const = NULL)
      all.f_1[[i]]= f_1
      
      rss1= f_1$deviance
      p1= length(pars1)
      num= ((rss1-rss2)/(p2-p1))
      if (l-p2 == 0) {
        den = (rss2/(l-p2+1))
        critical_val= qf(0.95, p2-p1, l-p2+1)     #actual value of F-distribution with df1= p2-p1 & df2= 1
      } else {
          den= (rss2/(l-p2))
          critical_val= qf(0.95, p2-p1, l-p2)     #actual value of F-distribution with
          }
      f_test= num/den
       
      
      if (f_test >= critical_val) {
        f_fitted= f_2
        fint=function(t) {f_2exp(f_2$par, t)}
        fitted_est= f_2exp(pars2 = f_fitted$par, t= data[,"t"])
      } else {
        f_fitted= f_1
        fint=function(t) {f_1exp(f_1$par, t)}
        fitted_est= f_1exp(pars1 = f_fitted$par, t= data[,"t"])
      }
    } 
    
    if (l %in% 3:4) { # no fit for l=1 or l=2 # 1_exp
      dl= log(data[,"y"]+ 1e-7)
      dd= data.frame(dl=dl, t= data[,"t"])
      lin_f_1= lm(dl~t, dd)
      pars1= c(P1= unname(exp(lin_f_1$coefficients["(Intercept)"])), 
               p1= -unname(lin_f_1$coefficients["t"]), 
               P= -1e-7)
      # rm(shared_const)
      # shared_const= NULL
      f_1= fit(data, func = f_1exp, start= pars1, shared_const= NULL)
      f_fitted= f_1
      fint=function(t) {f_1exp(f_1$par, t)}
      fitted_est= f_1exp(pars1 = f_fitted$par, t= data[,"t"])
    }
    
    all.model_fitted[[i]]= f_fitted
    fit_est[[i]]= fitted_est
    #names(all.model_fitted[[i]])= names(type)
    
    #### Area Under the Curve ####
    
    cnout= f_fitted$par
    #fint=function(t) {f_fitted(cnout, t)}
                    # AUC1= integrate(fint, max(1e-7, min(data[,"t"])), max(data[,"t"]))$value # if data starts from 0, then integrate from 1e-7
                    aucf= NULL
                    # tmp= floor(length(cnout)/2)   # could use integer division ( length(cnout)%/%2 )
                    # for (j in 2:tmp) 
                    # {
                    #   tmp2= 2*(j-1)+1
                    #   A=cnout[tmp2]
                    #   k=cnout[tmp2+1]
                    #   aucf[j-1]= (A/k)*exp(-k*max(data[,"t"]))
                    # }
                    # AUC2= sum(aucf, na.rm = T)
                    # AUC= AUC1 + AUC2
                    # rm(AUC1, aucf, AUC2, tmp, tmp2)
    
    for(jj in 1: (length(cnout)/2)) {
      A= cnout[2*jj-1]
      k= cnout[2*jj]
      aucf[jj]= A/k
      rm(A, k)
    }
    AUC_fitted= sum(aucf, na.rm = T)
    
    auc2= NULL
    cnout_2exp= f_2$par
    for(jj in 1: (length(cnout_2exp)/2)) {
      A= cnout_2exp[2*jj-1]
      k= cnout_2exp[2*jj]
      auc2[jj]= A/k
      rm(A, k)
    }
    AUC_2exp= sum(auc2, na.rm = T)
    #} 
    # else {
    #   #auc= NULL
    #   ff= function(t) {f(cnout, t)}
    #   AUC= integrate(f = ff, lower = 0, upper = 1e+03)
    # }
    AUC_l[[i]]= AUC_fitted
    AUC_l_2exp[[i]]= AUC_2exp
    #names(AUC_l[[i]])= names(type)

    #### assign weights==1/res^2 and do ordinary fit
    # cnout= coef(model_nout)
    # res_nout= resid(model_nout)
    # y_hat= f(cnout, t)
    # wres= res_nout/y_hat
    
    #rm(f_fitted, AUC, new_data)
  }
  
  names(AUC_l)= names(AUC_l_2exp)= colnames(data_original)
  names(all.model_fitted)= colnames(data_original)
  names(new_data_l)= colnames(data_original)
  names(fit_est)= colnames(data_original)
  names(all.f_1)= colnames(data_original)
  names(all.f_2)= colnames(data_original)
  names(all.f_3)= colnames(data_original)
  
  return(list(new_data= new_data_l, 
              fitted_model= all.model_fitted, 
              fitted_estimate= fit_est,
              AUC_fitted= AUC_l,
              AUC_2exp= AUC_l_2exp,
              outlier_indx= indx_l,
              all.f_1= all.f_1, all.f_2= all.f_2, all.f_3= all.f_3))
  
  # return(list(model_ord= all.model_ord, coeff_ord= all.coeff_ord,
  #             model_rob= all.model_rob, coeff_rob= all.coeff_rob,
  #             model_nout= all.model_nout, coeff_nout= all.coeff_nout,
  #             outlier_indx= indx_l,
  #             data_narm= data_l, data_orig= data_orig, data_nout= new_data_l,
  #             AUC= AUC_l)
  # )
  
}
