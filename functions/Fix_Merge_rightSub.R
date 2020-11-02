# Here are 3 functions: mydataFix, MyExcelMerge, and rightSub

        # MydataFix Finds and Removes
                    # COLUMNAMES==0 or NA
                    # DUPLICATED COLUMNS 
                    # ROWS WITH ALL NA VALUES

        # MyExcelMerge
                    # takes the name for excel file and spreadsheets 
                    # Returns all spreadsheets merged

        # rightSub
                  # takes all data
                  # returns the data for subjects starting with the character specified in starting_char



#Parisa Ghane (pghane@tamu.edu)

###mydataFix#######################################################################################
mydataFix =  function(data){
  
  dup_index= which(duplicated(colnames(data)))
  if (length(dup_index)!=0){
    new_data_1= data [, -dup_index]
  }else {
    new_data_1= data
  }
  zero_vals_indx= apply(new_data_1, 2, function(x) {   #remove columns with all values == 0
    all(as.numeric(x)==0)== FALSE
  } )
  zero_vals_indx[is.na(zero_vals_indx)]= TRUE
  new_data= new_data_1[, zero_vals_indx]
  
  tmpcoln= substr(colnames(new_data), 
                  start = 1, stop = 1
  )
  zero_col= which(tmpcoln== "0") #columns with name = "0"
  if (length(zero_col)!=0){
    data_nonzero= new_data[,-zero_col]
  }else {
    data_nonzero=new_data
  }
  
  na_col= which(colnames(data_nonzero)== "NA")
  if(length(na_col)!=0) {
    mydata_fixed= data_nonzero[, -na_col]
  }else {
    mydata_fixed= data_nonzero
  }
  
  ind= apply(mydata_fixed, 1, function(x) all(is.na(x)))
  nonNA_data_fixed= mydata_fixed[!ind,]
  rm(ind)
  
    return(nonNA_data_fixed)
}

###myExcelMerge####################################################################################

myExcelMerge= function(file_name, sheet_name){
  library(readxl)
  #source("T:/Pghane/Categorical_Analysis/functions/mydataFix.R")
  
  m= length(sheet_name)
  data_list= vector(mode = "list", length = m)
  
  for (i in 1:length(sheet_name)){
    raw_data= read_excel(file_name,
                         sheet = sheet_name[i],
                         col_names=TRUE
    )
    data_fixed= mydataFix(data = raw_data)
    data_list[[i]]= data_fixed 
    rm(raw_data, data_fixed)
  }
  
  
  all_data= Reduce( function (x,y)
    merge(x,y, 
          by= c("subjectid", "ctrid", "group"), 
          all= TRUE),
    data_list
  )
  return(all_data)
}


###rightSubj#######################################################################################
rightSubj= function(data, num_char, starting_char){
  data= as.data.frame(data)
  tempid= substr(data[ ,"subjectid"], 
                 start = 1, stop = num_char)
  data[tempid== starting_char, ]
}