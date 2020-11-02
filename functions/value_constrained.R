#initial values from data points using quadratic programming with constraints

#Parisa Ghane (pghane@tamu.edu) 
#==================================================================================================
initial_value_const = function (data_df, typ, constant=FALSE) {
  # install.packages("quadprog")
  library(quadprog)
  #library(quadprogXT)
  pars2= vector(mode = "list", length = ncol(data_df)-1)
  pars3= vector(mode = "list", length = ncol(data_df)-1)
  initial_est_3exp= vector(mode = "list", length = ncol(data_df)-1)
  initial_est_2exp= vector(mode = "list", length = ncol(data_df)-1)
  
  for (i in 2:ncol(data_df))
  {
    print(i)
    type= typ[i]
    if (is.na(type)) next
    if (type != "metabolite" & type != "pulse") next      # skip columns if type is not metabolite or pulse 
    
    if (type == "pulse") {
      data= data.frame(y= as.numeric(data_df[,i][-1]),
                       t= as.numeric(data_df[, 1][-1]))
    } else if (type == "metabolite") {
      data= data.frame(y= as.numeric(data_df[,i]),
                       t= as.numeric(data_df[, 1]))
    }
    indx0= which(is.na(data[,"t"]))
    if (length(indx0)>0) {
      data= data.frame(y= data[,"y"][-indx0],
                       t= data[,"t"][-indx0])
    }
    rm(indx0)

    data[,"t"] = data[,"t"]- min(data[,"t"])

    
    if (length(which(is.na(data[,"y"]))) == 0) {                     # remove NA values
      data= data
      #std= std_df
    }else {
      ind= which(is.na(data[,"y"]))
      data_name= rownames(data)
      na_name= data_name[ind]
      data= data.frame(y= data[,"y"] [as.numeric(setdiff(data_name, na_name))],
                             t= data[,"t"] [as.numeric(setdiff(data_name, na_name))])
      # data= data.frame(y= as.numeric(data_df[as.numeric(setdiff(data_name, na_name)) , i]), 
      #                  t= as.numeric(data_df[as.numeric(setdiff(data_name, na_name)) , 1]))
    }
    
    if (all(is.na(data[,"y"]))) next        #ignore cols with all NA values

    t= data[,"t"]
    d= (data[,"y"])
    l= length(d)
    
    if (l%/%3 >=3) {
      if (type == "metabolite") {
        ymax_indx= which(d == max(d))
        if (l == 9) {
          indx_1=1:(min(4, ymax_indx))
          indx_2= (min(4, ymax_indx)):6
          indx_3= 7:9
        } else {
          l_3rd= l%/%3
          indx_1= 1:(min(ymax_indx, l_3rd+1))
          indx_3= (l-l_3rd+1):l
          smpl_2_indx= !is.element(el = 1:l, set = union(indx_1,indx_3))
          indx_2= (1:l)[smpl_2_indx]
          rm(smpl_2_indx)
        }
      
    } else if (type == "pulse") {
        #d= DD
        l_3rd= l%/%3
        indx_1= 1: l_3rd
        indx_3= (l-l_3rd):l
        indx_2= l_3rd: (l-l_3rd)
        #smpl_2_indx= !is.element(el = 1:l, set = union(indx_1,indx_3))
        #indx_2= (1:l)[smpl_2_indx]
        #rm(smpl_2_indx)
    }
  
  #### 3rd exponential term ######
  #  take log(last 3 points), fit a line on them, 
  #  C3=exp(intercept) and c3=coeffitient of t, subtract C3*exp(c3*t) from all points
  t_3= t[indx_3]
  d_3= d[indx_3]
  if (any(d_3 == 0)) {d_3[d_3 == 0]= 1e-7}
  if (any(d_3 < 0)) {
    d_3= d_3+ (abs(min(d_3))) + 1e-7
  }
  
  d_3l= log(d_3)
  d_3l_ave= mean(d_3l)
  t_3_ave= mean(t_3)
  num= (d_3l_ave - d_3l)
  den= (t_3 - t_3_ave)
  den[den==0]= 1
  k_rng= range(num/den)
  lnA_rng= range(d_3l_ave - k_rng * t_3_ave)
  
  x= t(rbind(rep(1, length(t_3)), t_3))
  Rinv= solve(chol(t(x) %*% x))     # inverse of uppertriangular factor of the choleski decomposition
  Amat= t(cbind(c(1, -1, 0, 0),
                c(0, 0, 1, -1)))
  #Amat= -diag(ncol(x))
  # Amat[1]=0  # no constrain on log(coef) because coef= exp(fit$sol[1]) is always positive.  
  # bvec= c(-max(lnA_rng),
  #         min(k_rng))
  
  # bvec= c(min(lnA_rng),
  #             -max(lnA_rng),
  #             -max(k_rng), 
  #             min(k_rng))
  bvec= c(min(lnA_rng),
          -max(lnA_rng),
           -max(k_rng), 
          min(k_rng))
  
  dvec= t(d_3l) %*% x
  fit3= solve.QP(Dmat= Rinv, factorized = T, dvec, Amat, bvec)
  coef3= exp(fit3$solution[1])
  exp3= fit3$solution[2]
  
  d_fitted3= coef3* exp(exp3 * t)   #exp3 is already constrained to be negative
  d_updated3= d - d_fitted3
  
  # rm(k_rng, lnA_rng, x, Rinv, Amat, bvec, dvec)
  ## 2nd Exponential ##########
  t_2= t[c(indx_2)]
  d_2= d_updated3[c(indx_2)]
  
  # if any point in d_2 became negative --> remove negative points or shift all point to positives
  if (any(d_2 == 0)) {d_2[d_2 == 0]= 1e-7}
  if(any(d_2 < 0)) {
    d_2= d_2 + (abs(min(d_2))) + 1e-7
  } 
  
  d_2l= log(d_2)
  
  d_2l_ave= mean(d_2l)
  t_2_ave= mean(t_2)
  num= d_2l_ave - d_2l
  den= t_2 - t_2_ave
  den[den==0]= 1
  k_rng= range(num/den)
  lnA_rng= range(d_2l_ave - k_rng * t_2_ave)
  
  x= t(rbind(rep(1, length(t_2)), t_2))
  Rinv=solve(chol(t(x) %*% x))    # inverse of uppertriangular factor of the choleski decomposition
  Amat= t(cbind(c(1, -1, 0, 0),
                c(0, 0, 1, -1)))
  # Amat= -diag(ncol(x))
  # Amat[1]=0
  # bvec= c(-max(lnA_rng),
  #         min(k_rng))
  bvec= c(min(lnA_rng),
          -max(lnA_rng),
          -max(k_rng), 
          min(k_rng))
  
  dvec= t(d_2l) %*% x
  fit2= solve.QP(Dmat= Rinv, factorized = T, dvec, Amat, bvec)
  coef2= exp(fit2$solution[1])
  exp2= fit2$solution[2]
  
  d_fitted2= coef2* exp(exp2 * t)
  d_updated2= d_updated3 - d_fitted2
  
  # rm(k_rng, lnA_rng, x, Rinv, Amat, bvec, dvec)
  
  ###first exponential ######  
  if (type == "pulse") {
    d_1= d_updated2[c(indx_1)]
    t_1= t[c(indx_1)]
    if (any(d_1 == 0)) {d_1[d_1 == 0]= 1e-7}
    if(any(d_1 < 0)) {
      d_1= d_1 + (abs(min(d_1))) + 1e-7
    } 
    d_1l= log(d_1)
    
    if (length(d_1l) == 2) {
      fit1= list(solution= NA)
      fit1$solution[2]= diff(d_1l)/diff(t_1) 
      fit1$solution[1]= d_1l[1]-(fit1$solution[2]*t_1[1])
    } else if((length(d_1l) > 2)) {
      d_1l_ave= mean(d_1l)
      t_1_ave= mean(t_1)
      num= d_1l_ave - d_1l
      den= t_1 - t_1_ave
      den[den==0]= 1
      k_rng= range(num/den)
      lnA_rng= range(d_1l_ave - k_rng * t_1_ave)
      
      x= t(rbind(rep(1, length(t_1)), t_1))
      Rinv= solve(chol(t(x) %*% x))     # inverse of uppertriangular factor of the choleski decomposition
      Amat= t(cbind(c(1, -1, 0, 0),
                    c(0, 0, 1, -1)))
      # Amat= -diag(ncol(x))
      # Amat[1]=0
      # bvec= c(-max(lnA_rng),
      #         min(k_rng))
      bvec= c(min(lnA_rng),
              -max(lnA_rng),
              -max(k_rng), 
              min(k_rng))
      dvec= t(d_1l) %*% x
      fit1= solve.QP(Dmat= Rinv, factorized = T, dvec, Amat, bvec)
    } #else stop("inconsistent with type metabolite, first point is not smaller than other points")

      coef1= exp(fit1$solution[1])
      exp1= (fit1$solution[2])
    
    # rm(k_rng, lnA_rng, x, Rinv, Amat, bvec, dvec)
    
    d_fitted1= coef1*exp(exp1 * t)
  } else if (type == "metabolite") {
    d_1= d_updated2[c(indx_1)]
    t_1= t[c(indx_1)]
    #if (any(d_1 == 0)) {d_1[d_1 == 0]= 1e-7}
    if(any(d_1 <= 0)) {
      d_1= d_1 + (abs(min(d_1))) + 1e-7
    } 
    d_1_neg= -d_1
    rm(d_1)
    d_1= d_1_neg + abs(min(d_1_neg)) + 1e-7
    d_1l= log(d_1)
    
    if (length(d_1l) == 2) {
      fit1= list(solution= NA)
      fit1$solution[2]= diff(d_1l)/diff(t_1) 
      fit1$solution[1]= d_1l[1]-(fit1$solution[2]*t_1[1])
    } else if((length(d_1l) > 2)) {
      d_1l_ave= mean(d_1l)
      t_1_ave= mean(t_1)
      num= d_1l_ave - d_1l
      den= t_1 - t_1_ave
      den[den==0]= 1
      k_rng= range(num/den)
      lnA_rng= range(d_1l_ave - k_rng * t_1_ave)
      
      x= t(rbind(rep(1, length(t_1)), t_1))
      Rinv= solve(chol(t(x) %*% x))     # inverse of uppertriangular factor of the choleski decomposition
      Amat= t(cbind(c(1, -1, 0, 0),
                    c(0, 0, 1, -1)))
      # Amat= -diag(ncol(x))
      # Amat[1]=0
      # bvec= c(-max(lnA_rng),
      #         min(k_rng))
      bvec= c(min(lnA_rng),
              -max(lnA_rng),
              -max(k_rng), 
              max(0, min(k_rng)))
      dvec= t(d_1l) %*% x
      fit1= solve.QP(Dmat= Rinv, factorized = T, dvec, Amat, bvec)
    } #else stop("inconsistent with type metabolite, first point is not smaller than other points")
    
    #coef1= - exp((fit1$solution[1])/log(abs(min(d_1_neg))))
    exp1= (fit1$solution[2])
    if(exp1 == 0) {exp1= -1e-7}
    tmp0= abs(min(d_1_neg))/exp(exp1*t_1)
    tmp1= median(tmp0)
    coef1= -exp(fit1$solution[1]) - tmp1  # - abs(min(d_1_neg))/2)
    
    #rm(d_1_neg)
    #rm(k_rng, lnA_rng, x, Rinv, Amat, bvec, dvec)
    
    d_fitted1= coef1*exp(exp1 * t)
      }
    
    #### estimates & plot######
    initial_est_3exp [[i]]= d_fitted1 + d_fitted2 + d_fitted3
    
  if(constant) {
    pars3[[i]]= c(P1= unname(coef1), p1= -unname(exp1),  # -a1 because we have A1exp(-a1*t) in the function formula
                  P2= unname(coef2), p2= -unname(exp2), 
                  P3= unname(coef3), p3= -unname(exp3),
                  P= 0
    )
  } else {
    pars3[[i]]= c(P1= unname(coef1), p1= -unname(exp1),  # -a1 because we have A1exp(-a1*t) in the function formula
                  P2= unname(coef2), p2= -unname(exp2), 
                  P3= unname(coef3), p3= -unname(exp3)
    )
  }

  } else if (l%/%3 == 2) {
      
      if (type == "metabolite") {P1= -0.1} else {P1=0.1}
      
    if(constant) {
      pars3[[i]] = c(P1= P1 ,p1= 0.1,
                     P2= 0.01,p2= 0.01,
                     P3= 0.001,p3= 0.001
                     # P= 0.0001
      )
      
      initial_est_3exp [[i]]= (pars3[[i]]["P1"]* exp(-pars3[[i]]["p1"]*t)+
                                 pars3[[i]]["P2"]*exp(-pars3[[i]]["p2"]*t)+
                                 pars3[[i]]["P3"]*exp(-pars3[[i]]["p3"]*t)
                               # pars3[[i]]["P"]
      )
    } else {
      pars3[[i]] = c(P1= P1 ,p1= 0.1,
                     P2= 0.01,p2= 0.01,
                     P3= 0.001 #,p3= 0.001
                     )
      
      initial_est_3exp [[i]]= (pars3[[i]]["P1"]* exp(-pars3[[i]]["p1"]*t)+
                                 pars3[[i]]["P2"]*exp(-pars3[[i]]["p2"]*t)+
                                 pars3[[i]]["P3"]*exp(-pars3[[i]]["p3"]*t)#+
                               #pars3[[i]]["P"]
      )
    } 

      #}
      
    } else {
      pars3[[i]]= NULL
      initial_est_3exp[[i]]= NULL
    }
      

    ## for 2 exponentilas: ====================================================
    # one exponential fitted on first 3 points,
    # one exponential fitted on the rest (4th to last points)
    if (l >= 5) {
      
      if(l >= 9) {
        if(type == "metabolite") {
          ymax_indx= which(d == max(d))
          indx_1= 1:(min(ymax_indx, 6))
          indx_2= setdiff(1:l, indx_1)
          rm(ymax_indx)
        } else if (type == "pulse") {
          indx_1= 1:3
          indx_2= 4:l
        }
      } else if (l < 9) {
        if(type == "metabolite") {
          ymax_indx= which(d == max(d))
          indx_1= 1:(min(ymax_indx, 3))
          indx_2= (min(ymax_indx, 3)): l
          rm(ymax_indx)
        } else if (type == "pulse") {
          indx_1= 1:4
          indx_2= 4:l
        }
      }
      
      ### 2nd exp term ######
      t_2= t[indx_2]
      d_2= d[indx_2]
      if (any(d_2 <= 0)) {
          d_2= d_2+ (-min(d_2)) + 1e-6
      }
      
      d_2l= log(d_2)
      
      d_2l_ave= mean(d_2l)
      t_2_ave= mean(t_2)
      num= d_2l_ave - d_2l
      den= t_2 - t_2_ave
      den[den==0]= 1
      k_rng= range(num/den)
      lnA_rng= range(d_2l_ave - k_rng * t_2_ave)
      
      x= t(rbind(rep(1, length(t_2)), t_2))
      Rinv= solve(chol(t(x) %*% x))     # inverse of uppertriangular factor of the choleski decomposition
      Amat= t(cbind(c(1, -1, 0, 0),
                    c(0, 0, 1, -1)))
     
      bvec= c(min(lnA_rng),
              -max(lnA_rng),
              -max(k_rng), 
              min(k_rng))
      
      # bvec= c(min(lnA_rng),
      #         -min(0,max(lnA_rng)),
      #         min(1e-6, -max(k_rng)), 
      #         max(1e-5, min(k_rng)))

      dvec= t(d_2l) %*% x
      fit2_2exp= solve.QP(Dmat= Rinv, factorized = T, dvec, Amat, bvec)
      coef2_2exp= exp(fit2_2exp$solution[1])
      exp2_2exp= fit2_2exp$solution[2]
      
      d_fitted2= coef2_2exp* exp(exp2_2exp * t)
      d_updated2= d - d_fitted2
      
      rm(x, Rinv, Amat, bvec, dvec)
      
      ### 1st exp term ######
      d_1= d_updated2[c(indx_1)]
      t_1= t[c(indx_1)]
      #if (any(d_1 == 0)) {d_1[d_1 == 0]= 1e-5}
      if(any(d_1 <= 0)) {
        d_1= d_1 + (abs(min(d_1))) + 1e-6
      } 
      if (type == "metabolite") {
        d_1_neg= -d_1
        rm(d_1)
        d_1= d_1_neg + abs(min(d_1_neg)) + 1e-5
      }
      d_1l= log(d_1)
       #else stop("inconsistent with type metabolite, first point is not smaller than other points")
      
      if (length(d_1l) == 2) {
        fit1_2exp= list(solution= NA)
        fit1_2exp$solution[2]= diff(d_1l)/diff(t_1) 
        fit1_2exp$solution[1]= d_1l[1]-(fit1_2exp$solution[2]*t_1[1])
      } else if((length(d_1l) > 2)) {
        d_1l_ave= mean(d_1l)
        t_1_ave= mean(t_1)
        num= d_1l_ave - d_1l
        den= t_1 - t_1_ave
        den[den==0]= 1
        k_rng= range(num/den)
        lnA_rng= range(d_1l_ave - k_rng * t_1_ave)
        
        x= t(rbind(rep(1, length(t_1)), t_1))
        Rinv= solve(chol(t(x) %*% x))     # inverse of uppertriangular factor of the choleski decomposition
        Amat= t(cbind(c(1, -1, 0, 0),
                      c(0, 0, 1, -1)))
        # Amat= -diag(ncol(x))
        # Amat[1]=0
        # bvec= c(-max(lnA_rng),
        #         min(k_rng))
        bvec= c(min(lnA_rng),
                -max(lnA_rng),
                -max(k_rng), 
                min(k_rng))
        dvec= t(d_1l) %*% x
        fit1_2exp= solve.QP(Dmat= Rinv, factorized = T, dvec, Amat, bvec)
      } #else stop("inconsistent with type metabolite, first point is not smaller than other points")
      
    
    if (type == "metabolite") {
      exp1_2exp= (fit1_2exp$solution[2])
      if(exp1_2exp == 0) {exp1_2exp= -1e-7}
      tmp0_2exp= abs(min(d_1_neg))/exp(exp1_2exp*t_1)
      tmp1_2exp= median(tmp0_2exp)
      coef1_2exp= -exp(fit1_2exp$solution[1]) - tmp1_2exp  # - abs(min(d_1_neg))/2)
      
      rm(d_1_neg)
      # rm(k_rng, lnA_rng, x, Rinv, Amat, bvec, dvec)
      
      d_fitted1= coef1*exp(exp1 * t)
    } else if(type == "pulse") {
      coef1_2exp= exp(fit1_2exp$solution[1])
      exp1_2exp= (fit1_2exp$solution[2])
    }
    
      # rm(k_rng, lnA_rng, x, Rinv, Amat, bvec, dvec)
      
      d_fitted1= coef1_2exp*exp(exp1_2exp * t)

      
      #### estimates #####
      initial_est_2exp[[i]]= d_fitted1 + d_fitted2
      
      if(constant) {
        pars2[[i]]= c(P1= unname(coef1_2exp), p1= -unname(exp1_2exp), 
                      P2= unname(coef2_2exp), p2= -unname(exp2_2exp)
                      # P= 1e-6
        )
      } else {
        pars2[[i]]= c(P1= unname(coef1_2exp), p1= -unname(exp1_2exp), 
                      P2= unname(coef2_2exp), p2= -unname(exp2_2exp)#, 
                      #P= 1e-6
        )
      }

    } else if(l == 4) {
      if (constant) {
        pars2 [[i]] = c(P1= 0.05,
                        p1= 0.05,
                        P2= 0.005,
                        p2= 0.005
                        # P= 0.0001
        )
      } else {
        pars2 [[i]] = c(P1= 0.05,
                        p1= 0.05,
                        P2= 0.005,
                        p2= 0.005#,
                        #P= 0.0001
        ) 
      }

      
      initial_est_2exp [[i]]= (pars2[[i]]["P1"]* exp(-pars2[[i]]["p1"]*t)+
                                 pars2[[i]]["P2"]*exp(-pars2[[i]]["p2"]*t)#+
                                 #pars2[[i]]["P"]
                               )
      #}
    } else {
      pars2 [[i]]= NULL
      initial_est_2exp[[i]]= NULL
    }
    
    # plot(t,D, main= names(type))
    # points(t, initial_est_2exp[[i]], col="blue")
    # points(t, initial_est_3exp[[i]], col="red")
    
  }
  names(pars2)= names(typ)
  names(pars3)= names(typ)
  names(initial_est_2exp)= names(typ)
  names(initial_est_3exp)= names(typ)
  
  initial= list (pars2= pars2, 
                 pars3= pars3,
                 initial_estimate_2exp= initial_est_2exp,
                 initial_estimate_3exp= initial_est_3exp)
  
  initial
  
}
