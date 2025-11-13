# R Script - Statistical_Best_Practices

# Author: Maslard Corentin
# Email Address: corentin.maslard@gmail.com
# Created Date: July 12, 2023
# Last Updated: July 17, 2023

# Description: Clean statistics by following literature and certain courses. Certain values then go into other scripts to create graphs, for example.
# Work close to https://biostatgv.sentiweb.fr/?module=tests and and Abdou's stats courses. 

# Required packages ----
library(tidyverse)
library(readxl)
# Sources: ----
# - Source 1: Just a few useful little functions for scripting
source(here::here("src/function/stat_function/stat_utilities.R"))
source(here::here("src/function/stat_function/theme_plot.R"))


#need to add in this function
# - an export for all graph
# - Source 2: Description of source 2

# Quick summary of data ----

# Explain y with x ? (two numeric data)
# cathegorical variable ? 

# yes 
  # fisher test for each variable categorial (condition)
  # shapiro on anconva ressiduals if ns ancova if not two non parametric test by categorial variable (may be not answer to the question)

  # if interaction in X input (categorial variable plus numeric variable) is significant realise two simple regression by categorial variable

# no
# Shapiro pour voir distribution sur chaque variable.
## si les deux on une distribution normal alors correlation de pearson si l'une des deux est pas normale alors Spearman

# Mean comparison test ----

# Number of mean

## Two mean ---- 
  ### Independent data ?
  # yes 
    # Normal distribution of data for all variables ?   
      # !Shapiro test for all variables
      # yes 
        # Egality of variances  ?
        # ! Fisher test
          # yes
            # t-test !!!!       
        # no 
          # welch test or need to transform variable
      # no 
        # Need to transform or test de Man Whitney !!!
  # no 
    # Normal distribution of data by variable ? 
      # yes
        # paired data t-test !!!!
      # no
        # Wilcoxon test !!!!
## More than two means ----
  ### Independent data ?
  # yes
    # Egality of variances and normal distribution of model residuals? 
      # One variables with more than two modalities (if two is comparaison mean for two)
          # !Bartlett test for variance and !Shapiro test on anova residuals
            # yes
              # anova test 
              # post hoc Tukey
      # two or more modality with multiple factor in each variables  ?
        # Fisher test (for equality of variance) on each explanatory variable
        # !Bartlett test (on all combine variables) for variance and !Shapiro test on anova residuals
          # yes
            # anova test 
            # post hoc Tukey  
      # no (one one or tho variable and error to ressiduals)
        # one or multiple kruskal wallis test
          # post hoc, Test all pairs of means by Mannn-Whitney tests (with Bonferroni correction!)
  # no
    # Egality of variances and normal distribution of model residuals? 
      # !Bartlett test for variance and !Shapiro test on anova residuals
        # yes
        # anovaR test    
        # no
        # Friedman test

####################### Begin of this  function ###############
# normaly i will follow the previous condition
stat_analyse <- function(
    data,                   # the dataset to analyse
    column_value = 0,       # a name of column with numeric value, one or multiple 
    category_variables,     # a name or vector of multiple column name with categorical variable inside
    grp_var = "",           # if you want to make stats or variation rate intragroup, for example two genotypes and we compare each condition for genotype
    control_conditions = "",# if you have control condition and you want to see the variation rate
    independent_data = T,
    biologist_stats = F,
    
    # for plot
    show_plot = F,
    hex_pallet = manu_palett2,
    outlier_show=F,
    label_outlier = "plant_num",
    Ylab_i="y",
    strip_normale="", 
    strip = "",
    strip_fill_vec = "white",
    strip_text_vec = "black"
) {
  source(here::here("src/function/stat_function/test_functions/multiple_mean_comparison.R"))
  source(here::here("src/function/stat_function/test_functions/two_mean_comparison.R"))
  source(here::here("src/function/stat_function/plot_functions/plot_compare_two_independent_groups.R"))
  source(here::here("src/function/stat_function/plot_functions/plot_compare_more_than_two_groups.R"))
  source(here::here("src/function/stat_function/plot_functions/plot_compare_more_than_two_groups_two_factor.R"))
  source(here::here("src/function/stat_function/plot_functions/plot_utility.R"))
  
  if (independent_data != T) {
    cat_col("For the moment I haven't put in this function, tests with dependent data \n", color = "red") ; break
    
  }
  if (length(column_value) == 0) {
    cat_col("(These tests are not currently available. (Make an AFC-ACM) \n", color = "red") ; break
    
  }
  # summary of the stats (mean, sd and variation rate if you add the control value)
  if (length(column_value) == 2) { # it's correlation
    source(here::here("src/function/stat_function/sub_main_analysis/stat_analysis_main_corr.R"))
    
  } else if (length(column_value) > 2) {
    cat_col("(These tests are not currently available. (Make an ACP) \n", color = "red") ; break
    
  } else if (length(column_value) == 1) { # for mean comparison
    if (
      #category_variables != "" &&
      length(category_variables) == 1&& length(levels(as.factor(data[,category_variables]))) == 2 && grp_var == ""
    ) {
      cat_col("Comparison of two averages (non-apariate data)  if ok, T.test\n", color = "magenta")
      summary_result = stats_summary_tv(data, column_value, category_variables, grp_var, control_conditions)
      print(head(summary_result, n = 10))
      if (biologist_stats == T) {
        summary_result_stat = verification_of_t_test_assumptions(
          data_i = data,
          column_value,
          category_variables,
          force_student = T
        )
      } else {
        summary_result_stat = verification_of_t_test_assumptions(
          data_i = data,
          column_value,
          category_variables, 
          force_student = F
        )
      }
      # plot student 
      if (show_plot == T) {
        plot = plot_compare_two_independent_groups(
          data = summary_result_stat[["data_used"]],
          result_test = summary_result_stat[["result_test"]],
          x_val = category_variables,
          y_val = column_value,
          method = summary_result_stat[["test_output"]],
          hex_pallet = hex_pallet,
          Ylab_i=Ylab_i
        )
      }
    } else if (
    #category_variables != "" &&
      length(category_variables) == 1 &&
      grp_var != "" &&
      length(grp_var) == 1
    ) {
      cat_col("(ANOVA 2 factor but I prefer if you use cat variable = 2 for stats)\n", color = "magenta")
      anova_factor=2
      cat_col("(Usefull for one factor multiple and one comparison)\n", color = "green")
      summary_result = stats_summary_tv(data, column_value, category_variables, grp_var, control_conditions)
      print(head(summary_result, n = 10))
      if (biologist_stats == T) {
        summary_result_stat = verification_of_anova_multifactor_assumptions(
          data_i = data,
          column_value=column_value,
          category_variables=category_variables,
          grp_var=grp_var,
          anova_factor=anova_factor,
          force_ANOVA = T
        )
      } else {
        cat_col("(Not available for the moment)\n", color = "red")
      }
      
      
      if (show_plot == T) {
        #View(summary_result_stat[["data_used"]])
        plot = plot_one_variable_one_comparison(
          Z = summary_result_stat[["data_used"]],
          X = summary_result_stat[["data_used"]][, category_variables],#%>% mutate(group_cat_var=factor(group_cat_var,levels=c("Sto_WW_OT","Sto_WS_OT","Sto_WW_HS","Sto_WS_HS","Wen_WW_OT","Wen_WS_OT","Wen_WW_HS","Wen_WS_HS"))),
          Y = summary_result_stat[["data_used"]][, column_value],
          compile_cond=summary_result_stat[["data_used"]][,"compile_cond"],
          group_cat_var=summary_result_stat[["data_used"]][,"group_cat_var"],#%>% mutate(group_cat_var=factor(group_cat_var,levels=c("Sto_WW_OT","Sto_WS_OT","Sto_WW_HS","Sto_WS_HS","Wen_WW_OT","Wen_WS_OT","Wen_WW_HS","Wen_WS_HS"))),
          Xlab = category_variables,
          if(Ylab_i=="y"){
            Ylab =  column_value
          }else{
            Ylab=Ylab_i
          },
          Y_bis = column_value,
          outlier_show = outlier_show,
          label_outlier = label_outlier,
          category_variables = category_variables,
          grp_var=grp_var,
          if(summary_result_stat[["test_output"]]=="anova"|biologist_stats==T){
            Tukey = T  
          }else if (summary_result_stat[["test_output"]]=="kruskal"){
            Tukey = F
          }else{
            Tukey=F
          },hex_pallet=hex_pallet,
          strip_normale=strip_normale
          
        )
      }
      data_stats = plot[["Z"]]
      plot = plot[["p"]]
      
      
    } else if (
      length(category_variables) == 1 &&
      length(levels(as.factor(data[,category_variables])))  > 2 &&
      grp_var == "" | 
      
      length(category_variables) > 1 |
      length(levels(as.factor(category_variables))) > 2 
    ) {
      #cat_col("(This is an Anova of X factor)\n", color = "magenta") # !!!!!!!!!!!!!!!!!!!!
      if (
        length(category_variables) == 1 &&
        length(levels(as.factor(data[,category_variables])))  > 2 &&
        grp_var == ""
      ) {
        cat_col("(ANOVA 1 foctor)\n", color = "magenta") # !!!!!!!!!!!!!!!!!!!!
        anova_factor=1
      }else if (   
        length(category_variables) > 1 |
        length(levels(as.factor(category_variables))) > 2 # peut importe si grp_var ou non car je cumule dans la fonction
      ) {
        
        if(grp_var==""){
          list_cat_var=category_variables
        }else{
          list_cat_var=c(grp_var,category_variables)
        }
        anova_factor=length(list_cat_var)
        cat_col(paste0("(ANOVA of ",anova_factor," categorial factor)\n"), color = "magenta")
      }
        summary_result = stats_summary_tv(data, column_value, category_variables, grp_var, control_conditions)
        print(head(summary_result, n = 10))
        
        if (biologist_stats == T) {
          summary_result_stat = verification_of_anova_multifactor_assumptions(
            data_i = data,
            column_value=column_value,
            category_variables=category_variables,
            grp_var=grp_var,
            anova_factor=anova_factor,
            force_ANOVA = T
          )
        } else {
          summary_result_stat = verification_of_anova_multifactor_assumptions(
            data_i = data,
            column_value=column_value,
            category_variables=category_variables,
            grp_var=grp_var,
            anova_factor=anova_factor,
            force_ANOVA = F
          )
          
          if (summary_result_stat[["variable_transformation_requirements"]] == T) {
            p_input <- as.character(readline(prompt = "Perform a transformation of the numerical variable, if so enter the type of transformation (log, sqrt)?: "))
            # Remove quotation marks if present
            p_input <- gsub("\"|\'", "", p_input)
            print(p_input)
            if (p_input == "log") {
              cat_col("(Currently being created ANOVA2bis transformation log)\n", color = "magenta")
              summary_result_stat = verification_of_anova_multifactor_assumptions(
                data_i = data,
                column_value,
                category_variables,
                grp_var,
                anova_factor=anova_factor,
                force_ANOVA = F,
                transformation = "log"
              )
            } else if (p_input == "sqrt") {
              cat_col("(Currently being created ANOVA2bis transformation square root)\n", color = "magenta")
              summary_result_stat = verification_of_anova_multifactor_assumptions(
                data_i = data,
                column_value=column_value,
                category_variables=category_variables,
                grp_var=grp_var,
                anova_factor=anova_factor,
                force_ANOVA = F,
                transformation = "sqrt"
              )
            }
          }
        }
      
      # plot
      if (show_plot == T) {
        #View(summary_result_stat[["data_used"]])
        plot = plot_multi_mean(
          Z = summary_result_stat[["data_used"]],
          X = summary_result_stat[["data_used"]][, category_variables],#%>% mutate(group_cat_var=factor(group_cat_var,levels=c("Sto_WW_OT","Sto_WS_OT","Sto_WW_HS","Sto_WS_HS","Wen_WW_OT","Wen_WS_OT","Wen_WW_HS","Wen_WS_HS"))),
          Y = summary_result_stat[["data_used"]][, column_value],
          compile_cond=summary_result_stat[["data_used"]][,"compile_cond"],
          group_cat_var=summary_result_stat[["data_used"]][,"group_cat_var"],#%>% mutate(group_cat_var=factor(group_cat_var,levels=c("Sto_WW_OT","Sto_WS_OT","Sto_WW_HS","Sto_WS_HS","Wen_WW_OT","Wen_WS_OT","Wen_WW_HS","Wen_WS_HS"))),
          Xlab = category_variables,
          if(Ylab_i=="y"){
            Ylab =  column_value
            }else{
            Ylab=Ylab_i
          },
          Y_bis = column_value,
          outlier_show = outlier_show,
          label_outlier = label_outlier,
          category_variables = category_variables,
          grp_var=grp_var,
          if(summary_result_stat[["test_output"]]=="anova"|biologist_stats==T){
            Tukey = T  
          }else if (summary_result_stat[["test_output"]]=="kruskal"){
            Tukey = F
          }else{
            Tukey=F
          },hex_pallet=hex_pallet,
          strip = strip,
          strip_fill_vec = strip_fill_vec,
          strip_text_vec = strip_text_vec
          
        )
        
        # MAJ 20250312
        summary_result_stat$data_used = plot$data_used
        plot = plot$p
      }
    
    } else {
      cat_col(
        paste0(
          "I haven't set up this possibility yet. I don't know what test to do :( \n",
          "column_value = ", length(column_value), ":\t", column_value, "\n"  ,    
          "category_variables = ", length(category_variables), ":\t", category_variables, "\n",
          "grp_var = ", length(grp_var), ":\t", grp_var, "\n"
        ),
        color = "red"
      )
      break
    }
  }
  
  summary_result_stat = c(summary_result_stat, list(summary_result = summary_result, plot = plot, data_stats = if(exists("data_stats")){data_stats}else{"NA"}))
  return(summary_result_stat)
}

# ######### data test
# result_test=stat_analyse(
#   data=df_global_archi %>% drop_na("diff_correct_plant_height"),
#   #column_value = c("diff_correct_plant_height"),
#   column_value = c("diff_correct_plant_height","root","area","perimeter","largeur","density","length_skull"),
#   category_variables = c("water_condition","heat_condition"),
#   grp_var = "genotype",biologist_stats = T,control_conditions = c("WW","OT")
# )
# 
# result_test2=stat_analyse(
#   data=df_global_archi %>% drop_na("diff_correct_plant_height"),
#   #column_value = c("diff_correct_plant_height"),
#   column_value = c("root"),
#   category_variables = c("water_condition","heat_condition"),
#   grp_var = "genotype",biologist_stats = T,control_conditions = c("WW","OT")
# )
# 
# result_test=stat_analyse(
#   data=df_global_archi %>% drop_na("diff_correct_plant_height") %>% filter(!plant_num %in% c(1005,1010,1059,1072,1063,1118,1069,1128)),
#   #column_value = c("diff_correct_plant_height"),
#   column_value = c("root"),
#   category_variables = c("water_condition","heat_condition"),
#   grp_var = "genotype"
# )
# 
# result_test=stat_analyse(
#   data=df_global_archi %>% drop_na("diff_correct_plant_height") %>% filter(!plant_num %in% c(82,1005,1010,1069,1128,1072,1063,1118,1104,1021)),
#   #column_value = c("diff_correct_plant_height"),
#   column_value = c("root"),
#   category_variables = c("water_condition","heat_condition"),
#   grp_var = "genotype",control_conditions = c("WW","OT")
# )
# 
# result_test[["result_stat_anova"]]
# 
# test=list(root=result_test, result_test2)
# class(result_test)
# 
# test[["root"]]
# 
# # Utilisation de la fonction pour compiler le tableau "outliers"
# tableau_outliers <- compile_table(result_test, "outliers")
# occurrences <- tableau_outliers %>%
#   dplyr::count(plant_num) %>% 
#   arrange(desc(n))
# 
# print(occurrences)
# 
# tableau_stat_anova <- compile_tableau(result_test, "result_stat_anova")
# 
# tableau_stat_anova["root",]
# 
# 
# 
# # verification of the ionomic data ----
# # importation of the data
# 
# df_physio_iono=data_import_plant_num(df_physio_i = T, df_iono_conc_i = T)
# 
# colnames(df_physio_iono)
# vec_column=colnames(df_physio_iono[8:length(df_physio_iono)])
# 
# vec_column=colnames(df_physio_iono[19:length(df_physio_iono)])
# 
# 
# result_test=stat_analyse(
#   data=df_physio_iono %>% drop_na(vec_column),
#   #column_value = c("diff_correct_plant_height"),
#   column_value = vec_column,
#   category_variables = c("water_condition","heat_condition"),
#   grp_var = "genotype",control_conditions = c("WW","OT")
# )
# 
# tableau_outliers <- compile_table(result_test, "outliers")
# occurrences <- tableau_outliers %>%
#   dplyr::count(plant_num) %>% 
#   arrange(desc(n))
# 
# tableau_stat_anova <- compile_table(result_test, "result_stat_anova")
# test=tableau_stat_anova[,"coefficients"]$Be9_stem_concentration
# test[5]
# 
# tableau_stat_anova[["Be9_stem_concentration"]]$`Pr(>F)`[1]
# 
# occurrences
# 
# length(vec_column)
