

#1st RUN allfitall_fit_exp.R

# arguments should be lists with same length (elements sorted respective to each other!!)
fit.goodness= function(fits, pars, fit_names) 
{
  all_mean_res= NULL
  all_std_res= NULL
  wrss= NULL
  aic= NULL
  
  for (ii in 1:length(fits)) 
  {
    print(paste("fit ", ii))
    fit= fits[[ii]]
    model= fit$model_nout
    #w= fit$weight
    data_orig= fit$data_orig
    outlier_indx= fit$outlier_indx
    par= pars[[ii]]

    wres= matrix(rep(0,(length(model)-1)*length(time)),     # matrix of wres: cols==subjects, rows== time, 
                  ncol = length(model)-1)
    
    
    for (j in 2:length(model))
    {
      print(paste("modelcol", j))
      dev= model[[j]]$deviance
      wrss[[ii]]= dev/sqrt(var(resid(model[[j]])))
      aic[[ii]]= wrss[[ii]] + 2* length(par)
     
      
       #print(j)
      r= (resid(model[[j]]))
      res= r/sqrt(var(resid(model[[j]])))
      indx= outlier_indx[[j]]           #set outliers to NA
      if (length(indx) != 0) {
        wres[indx, j-1] = NA
        orig_na_indx= which(is.na(data_orig[[j]]))
        if (length(orig_na_indx) == 0) {
          #w_new= w[-indx]
          wres[-indx, j-1] = res
        } else {
          wres[orig_na_indx, j-1] = NA
          #w_new= w[-c(indx, orig_na_indx)]
          wres[-c(indx, orig_na_indx), j-1] = res
        }

      } else {
        orig_na_indx= which(is.na(data_orig[[j]]))
        if (length(orig_na_indx) == 0) {
          wres[ ,j-1] = res
        } else {
          wres[orig_na_indx, j-1] = NA
          wres[-orig_na_indx, j-1] = res
        }

      }
    }

    #meanwres_nor= -1 + 2*(mean_res - min(mean_res)) / diff(range(mean_res))    #normalized to [-1,1]

    #AIC criterion
    # wrss_norm= (wrss-mean(wrss))/ sqrt(var(wrss))
    # aic= wrss + 2* length(par)
    # aic_mean= mean(aic)
    # aic_std= sqrt(var(aic))


     mean_res_tmp= apply(wres, MARGIN = 1, mean, na.rm=T)
     mean_res= data.frame(res= mean_res_tmp, fnc= fit_names[ii])
    all_mean_res= rbind(all_mean_res, mean_res)
    std_res_tmp= apply(wres, MARGIN = 1, function(x){sqrt(var(x, na.rm = T))})     # mean & std per timepoint for all patients
    std_res= data.frame(std= std_res_tmp, fnc= paste( fit_names[ii]))
    all_std_res= rbind(all_std_res, std_res)

    #mean_res[[ii]]= apply(wres, MARGIN = 2, mean, na.rm=T)
    #std_res[[ii]]= apply(wres, MARGIN = 2, function(x){sqrt(var(x, na.rm = T))})     # mean & std per timepoint for all patients
    #wrss[[ii]]= apply(wres, MARGIN = 1, function(x) {sum(x^2, na.rm = TRUE)})             # wrss per person
    #aic[[ii]]= wrss[[ii]] + 2* length(par)

    #rm(model, mean_res_tmp, mean_res, std_res_tmp, std_res)

  }

  ggp_d= cbind(all_mean_res, all_std_res)
  #meanwres_nor= -1 + 2*(mean_res - min(mean_res)) / diff(range(mean_res))    #normalized to [-1,1]
  
  #AIC criterion
  #wrss_norm= (wrss-mean(wrss))/ sqrt(var(wrss))
 
  #aic_mean= mean(aic)
  #aic_std= sqrt(var(aic))
  
  new_time= rep(time, length(fits))

  # plot(time, mean_res[[1]],
  #      ylim=range(c(mean_res[[1]]-std_res[[1]], mean_res[[1]]+std_res[[1]])),
  #      xlab= "Time (min)",
  #      ylab= "weighted Residuals",
  #      main= "average of residuals from all subjects",
  #      pch=19)
  # lines(time, mean_res)
  # arrows(time, mean_res[[1]]-std_res[[1]], time, mean_res[[1]]+std_res[[1]], 
  #        length=0.05, angle=90, code=3, col= colors()[+13])
  # abline(h =0 , col= "gray")
  # 
  # for (i in 1: length(fits))
  #   
  # points()
  

  # plot(time, mean_df, 
  #      ylim=range(c(mean_df-std_df, mean_df+std_df)),
  #      xlab= "Time (min)",
  #      ylab= "PHE6",
  #      main= "average fit with error bars",
  #      pch=19)
  # arrows(time, mean_df-std_df, time, mean_df+std_df, length=0.05, angle=90, code=3)
  
  list(wrss=wrss, aic= aic, ggp_d= ggp_d, 
       new_time= new_time)
  
}

