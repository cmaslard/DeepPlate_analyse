# R Script - multiple mean comparison mor than two mean

# Author: Maslard Corentin
# Email Address: corentin.maslard@gmail.com
# Created Date: July 28, 2023
# Last Updated: July 28, 2023

# Description: two mean comparison with all verification for T.test
# Web source see ...

# Pkg 
library(rstatix)
library(ggpubr) # for the ggqqplot

# ok if length(column_value)==1&&length(category_variables)==1

# Test for equality of variances for consisdiated explanatory variables

## Fisher test for each variable considered for the two factor categorical variables 

verification_of_t_test_assumptions<-function(data_i,column_value,category_variables,variable_transformation_requirements=F,transformation="",force_student=F){
  #data_i$group=paste(sep="_",data_i$genotype, data_i$water_condition, data_i$heat_condition)
  
    list_cat_var=category_variables
  
  if(transformation=="log"){
    data_i[,column_value] = log(data_i[,column_value])
  }
  
  if(transformation=="sqrt"){
    data_i[,column_value] = sqrt(data_i[,column_value])
  }
  
  #Identify outliers by group: ----
  out_group=data_i %>%
    dplyr::group_by_at(list_cat_var) %>%
    rstatix::identify_outliers(column_value)
  n_out_group=nrow(out_group)
  
  if(n_out_group>0){
    cat_col(paste0("- Warning, there is (are): ",n_out_group," outlier (s) that can lead to false results in student test and in tests of equality of variances. (in the returne call 'outliers') \n"),"red")
    #print(out_group)
  }else{
    cat_col("- No outlier found \n","green")
  }
  
       # Creation of the new column
  data_i$group <- "" 
  for (i in list_cat_var) {
    data_i$group <- paste(data_i$group, data_i[, i], sep = "")
  }
  
  #collect the two condition to see normality
  vect_var=levels(as.factor(data_i[,category_variables]))
  ## Test the normality hypothesis by group ----
  
  normality_by_group=data_i %>%
    dplyr::group_by_at(list_cat_var) %>%
    shapiro_test(!!as.symbol(column_value))
  
  n_F_normality_group= normality_by_group %>% filter(p < 0.05) %>% nrow()
  
  if(n_F_normality_group>0){
    cat_col(paste0("- Normality on the two factor.  Warning, there is (are): ",disp_round(n_F_normality_group)," group (s) that do not respect the normality hypothesis (in the returne call 'normality_by_group') \n"),"red")
    
  }else{
    cat_col("- Normality on the model residuals. Shapiro Wilk validate for each groupe \n","green")
  }
  
  pval_fisher= as.numeric(var.test(data_i[,column_value]~data_i[,category_variables])["p.value"])
  if(pval_fisher>0.05){
    cat_col(paste0("- Homogeneity of variances, Fisher's test is not significant (p > 0.05):",disp_round(pval_fisher), ". Consequently, we can assume homogeneity of variances in the two different factor.", "\n"),color = "green")
  }else{
    cat_col(paste0("- Homogeneity of variances, Fisher's test is significant (p < 0.05):",disp_round(pval_fisher), ". Consequently, we can't assume homogeneity of variances in the two different factor. We therefore normally use a non-parametric test or Welch teste", "\n"),color = "red")
  }
  
result_t_test<- t.test(data_i [data_i[,category_variables]==vect_var[1],column_value], data_i [data_i[,category_variables]==vect_var[2],column_value],var.equal = TRUE) #assume equality of variance
result_welch_test<- t.test(data_i[data_i[,category_variables]==vect_var[1],column_value], data_i [data_i[,category_variables]==vect_var[2],column_value])
result_u_test <- wilcox.test(data_i [data_i[,category_variables]==vect_var[1],column_value], data_i [data_i[,category_variables]==vect_var[2],column_value])

  
  if(force_student==T){
    cat_col("Attention you have forced the test to make the student. If you're a biologist, this may be understandable.But be careful with interpretations \n",color = "yellow")
    cat_col("Here are the Student results: \n",color = "red")
    print(result_t_test)
    test_output="t.test" #think to add method.args = list(var.equal = TRUE)
    
  }else{
    if (n_F_normality_group>0) {
      cat_col("You need to use a non-parametric test or transform variables !!! \n",color = "white")
      variable_transformation_requirements=T
      
      # Effectuer le test de Mann-Whitney
      cat_col("Here are the non parametric test of Mann-Whitney or the Wilcoxon rank-sum  results: \n",color = "yellow")
      print(result_u_test)
      test_output="wilcox.test"
      
      }else if (pval_fisher>0.05& n_F_normality_group==0){
        cat_col("Here are the Student results: \n",color = "green")
        print(result_t_test)
        test_output="t.test" #think to add method.args = list(var.equal = TRUE)

      }else{
        cat_col("Here are the Welch results: \n",color = "red")
        print(result_welch_test)
        test_output="welch" #think to add method.args = list(var.equal = FALSE)
    }
  }

  return(list(outliers=out_group,
              normality_by_group=normality_by_group,
              result_test=list(result_t_test=result_t_test,result_u_test=result_u_test,result_welch_test=result_welch_test),
              variable_transformation_requirements=variable_transformation_requirements,
              data_used=data_i,
              test_output=test_output, #all test performe att this end (parametric or not)
              vect_var=vect_var
              ))
}

