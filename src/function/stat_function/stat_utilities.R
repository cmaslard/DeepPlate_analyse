# R Script - stats short tools

# Author: Maslard Corentin
# Email Address: corentin.maslard@gmail.com
# Created Date: July 12, 2023
# Last Updated: July 12, 2023

# Description: Very useful little functions for my different tests
cat_col <- function(text, color) {
  colors <- c("red", "green", "yellow", "blue", "magenta", "cyan", "white")
  color_code <- switch(tolower(color),
                       "red"     = "31",
                       "green"   = "32",
                       "yellow"  = "33",
                       "blue"    = "34",
                       "magenta" = "35",
                       "cyan"    = "36",
                       "white"   = "37")
  
  # Vérifier si la couleur est valide
  if (color_code == "") {
    cat("Couleur non valide.")
    return(invisible())
  }
  
  # Afficher le texte avec la couleur spécifiée
  cat(paste0("\033[", color_code, "m", text, "\033[0m"))
}

# Display dot after comma  ----
disp_round<-function(pval,rnd=0.001){
  if(pval<rnd){
    pval=paste0("<",rnd)
  }else{
    pval=paste0("= ",round(pval,3))
  }
    return(pval)
}

# compile_table function ----
#The compile_table function takes as input a list (my_list) and the name of a table (table_name). Its main goal is to compile the specified table from each element in the list into a single consolidated table.
#To achieve this, the function uses lapply to iterate over each element in my_list and extract the specified table using the provided table_name. The individual tables are then combined using do.call(rbind, tables) to create a single, merged table.
#The function returns the compiled table as the output.


compile_table <- function(my_list, table_name) {
  # Compiling the tables
  tables <- lapply(my_list, function(x) x[[table_name]])
  compiled_table <- do.call(rbind, tables)
  
  # Converting to data.frame
  #compiled_table <- as.data.frame(compiled_table)
  
  return(compiled_table)
}


## calculates the mean, standard deviation and number of occurrences for a given value (column_value) and multiple or one cathegories variables. I can also add a grp value which is a grouping. For example, several genotypes. 
stats_summary <- function(data, column_value_1, category_variables,grp_var="") {
  library(dplyr)
  
  # Select specified columns
  if(grp_var==""){
    colonnes <- c(column_value_1, category_variables)
    donnees <- data[colonnes]
    
    if (anyNA(donnees)) {
      stop("The dataset contains missing values (NA). Use drop_na on the numeric variable before using this function.")
    }
    
    grp=c(category_variables)  
  }else{
    colonnes <- c(column_value_1, category_variables,grp_var)
    donnees <- data[colonnes]
    if (anyNA(donnees)) {
      stop("The dataset contains missing values (NA). Use drop_na on the numeric variable before using this function.")
    }
    grp=c(category_variables,grp_var)  
  }
  # Calculate mean and standard deviation by combining categorical variables
  stats <- donnees %>%
    group_by(across(all_of(grp))) %>%
    dplyr::summarise(Mean = mean(!!as.symbol(column_value_1)),
              Sd = sd(!!as.symbol(column_value_1)),
              N=n(),
              .groups = "drop" # to delete the message
              )%>% 
    mutate(column_value = column_value_1) %>%
    relocate(column_value, .before = 1)
  
  # Return results
  return(stats)
}

# calculates the rate of change. say the cathegoriel column and say in that column in order the control conditions (normally only one per column). Then check the table
calculate_variation_rate <- function(data, grp_var="", category_variables, control_conditions) {
  if(grp_var==""){
    vr=data %>%
      dplyr::mutate(type_cond = {
        condition = rep(TRUE, n())
        for (i in seq_along(category_variables)) {
          condition = condition & (get(category_variables[i]) == control_conditions[i])
        }
        ifelse(condition, "controle", "stress")
      })
  }else{
    if (length(category_variables)!=length(control_conditions)) {
      stop("The number of columns with categorical data does not correspond to the number of control conditions.")
    }
    vr=data %>%
      dplyr::group_by(across(all_of(grp_var))) %>%
      dplyr::mutate(type_cond = {
        condition = rep(TRUE, n())
        for (i in seq_along(category_variables)) {
          condition = condition & (get(category_variables[i]) == control_conditions[i])
        }
        ifelse(condition, "controle", "stress")
      })
  }
  vr$VarRate <- with(vr, ifelse(type_cond == "stress", ((Mean - Mean[type_cond == "controle"]) / Mean[type_cond == "controle"]) * 100, NA))
  
  return(vr)
}

#compile the to function before
stats_summary_tv <-function(data,column_value_1,category_variables,grp_var="",control_conditions=""){
  # Make table resume mean sd and n by grouping by grp_var (ex: genotype) and by category_variables
  
  result <- stats_summary(data,column_value_1,category_variables,grp_var) 
  #print(result)
  # Calculate the variation rates for each group and add to previous table
  
  if(control_conditions[1]!=""){
  result <- calculate_variation_rate(result, grp_var, category_variables, control_conditions)
  }else{
   cat_col(paste0("! We can't do the rate of change for '",column_value_1 ,"' because you haven't entered the control parameters condition \n"),"yellow") 
  }
  # Return results
  return(result)
}

# Same as above but just repeat if there are more than one numerical value.
stats_summary_tv_n<-function(data,column_value,category_variables,grp_var="",control_conditions=""){
  result=data.frame()
  n_column_value=length(column_value)
  if(n_column_value==1){
    #cat("! You can also use the stats_summary_tv function to generate statistics for a single column with numerical values.\n") 
    cat_col("! You can also use the stats_summary_tv function to generate statistics for a single column with numerical values.\n", "green")    
    
    result=stats_summary_tv(data,column_value,category_variables,grp_var,control_conditions) 
  }else{
    for (i in column_value){
      result_i=stats_summary_tv(data,column_value_1 = i,category_variables,grp_var,control_conditions)
      result<-rbind(result,result_i)
    }
  }
  return(result)
}

# 
# ################################# test ############################
# source(here::here("xp1_analyse/script/global_analyse/function/data_importation.R"))
#  df_global_archi=data_import_kinetic(df_ra_i = T,df_height_plant_i = T,df_root_tips_nb_10000_i = F,df_root_tips_angle_i = F,df_root_tips_length_i = F,df_root_nb_ramif_i = F,key = "plant_num_shooting_date") %>%
#   mutate(condition=factor(condition,levels=c("Sto_WW_OT","Sto_WS_OT","Sto_WW_HS","Sto_WS_HS","Wen_WW_OT","Wen_WS_OT","Wen_WW_HS","Wen_WS_HS"))) %>%
#   mutate(climat_condition=factor(climat_condition,levels=c("WW_OT","WS_OT","WW_HS","WS_HS"))) %>%
#   filter(days_after_transplantation!=3) %>% #all data compile by taskid, shooting date and plant num
#   filter(days_after_transplantation ==19) %>%
#   dplyr::select(-one_of(c("nod","profondeur","background","days_after_transplantation"))) %>%
#   dplyr::select(-contains("_3")) %>%
#   mutate(plant_num=as.numeric(as.character(plant_num)))
# 
# grp_var <- c("genotype")
# group_columns <- c("")
# category_variables <- c("water_condition","heat_condition")
# category_variables <- c("water_condition","heat_condition","genotype")
# control_conditions <- c("WW","OT")
# control_conditions <- c("WW","OT","Wendy")
# 
# stats_summary(data=df_global_archi,column_value_1="root",category_variables = c("water_condition","heat_condition"),grp_var = "genotype")
# 
# stats_summary_tv(df_global_archi %>% drop_na("diff_correct_plant_height"),column_value_1 = "diff_correct_plant_height",grp_var = "genotype",category_variables = c("water_condition","heat_condition"),control_conditions = c("WW","OT"))
# result_test <- stats_summary_tv(data=df_global_archi,column_value_1="root",category_variables = c("water_condition","heat_condition"),grp_var="genotype",control_conditions = c("WW","OT"))
# result_test <- stats_summary_tv(data=df_global_archi,column_value_1="root",category_variables,grp_var="",control_conditions) ; result_test
# 
#  
# stats_summary_tv_n(df_global_archi %>% drop_na("diff_correct_plant_height"),column_value = c("diff_correct_plant_height","root"),grp_var = "genotype",category_variables = c("heat_condition"),control_conditions = "OT")
