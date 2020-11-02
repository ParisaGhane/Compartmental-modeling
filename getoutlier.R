
#Parisa Ghane (pghane@tamu.edu)

###################################################################################################


getoutlier= function (model) 
{
  #tr_val= true_value
  
  indx= vector(mode = "list", length = length(model))
  
  for (i in 2:length(model))
  {
    res= resid(model[[i]])
    std= sqrt(var(res))
    m= mean(res)
    thr1= m - 2* std
    thr2= m + 2*std
      
#     rng= max(res)-min(res)
#     thr= rng*0.8 + min(res) 
    
    indx[[i]]= which(res<thr1 | res > thr2)
    
  }
  indx
#   f= function(pars, t){
#     expr = expression(pars$A* exp(-pars$a*t)+
#                         pars$B*exp(-pars$b*t)+pars$C)
#     eval(expr)
#   }
#   
#   fval= f(as.list(coef(model.2exp)))
#   dist(fval, tr_val)
          
}

rm_outlier= function(data_df, out_indx) 
  {
  data_l= as.list(data_df)
  
  for (i in 2:length(out_indx)) 
  {
    ind= out_indx[[i]]
    if (length(ind) == 0) 
    {
      data_l[[i]]=data_l[[i]]
      }else 
      {
        data_l[[i]][c(ind)]= NA
      }
  }
  data_out= as.data.frame(data_l)
  data_out
}


  