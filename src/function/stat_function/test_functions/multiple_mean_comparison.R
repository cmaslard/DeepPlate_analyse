# R Script - multiple mean comparison mor than two mean

# Author: Maslard Corentin
# Email Address: corentin.maslard@gmail.com
# Created Date: July 13, 2023
# Last Updated: august 10, 2023

# Description: multiple mean comparison with all verification for ANOVA test
# Web source see https://www.datanovia.com/en/fr/lessons/anova-dans-r/ 

# Pkg 
library(rstatix)
library(ggpubr) # for the ggqqplot

# ok if length(column_value)==1&&length(category_variables)>1

# Test for equality of variances for consisdiated explanatory variables

## Fisher test for each variable considered and Bartlett test for the combination of all categorical variables 
# je ne fais pas de Bartlette mais leven a la place. Car je dois le faire pour chaque groupe sinon. De plus, . Si les données sont normalement distribuées, le test de Bartlett peut être utilisé, mais dans le cas contraire, le test de Levene est souvent préfére
# tres complexe d'analyser une anova a plus de trois facteur.

#Si un facteur en posthoc fais un tuckey et non parametrique shapiro
# si deux facteur ou plus complexe de faire du post hoc tuckey surtout si on a de l'interaction entre les facteurs. mais possible
# en non parametrique on peut faire du LOES si plus d'un facteur. Sinon a chaque fois on fait un choix, et on réduit pour avoir moins de facteur mais avec plus de variable (combine les variables) 
# pour le moment vue que je suis biologiste je vais compiler dans un même facteur pour faire une anova simple (a 8 variable) et faire un tuket. Normalement pas le droit. :(
# 

verification_of_anova_multifactor_assumptions<-function(data_i,column_value,category_variables,grp_var,variable_transformation_requirements=F,transformation="",anova_factor,force_ANOVA=F){
  data_i_init=data_i
  #data_i$group=paste(sep="_",data_i$genotype, data_i$water_condition, data_i$heat_condition)
  if(grp_var==""){
    list_cat_var=category_variables
  }else{
    list_cat_var=c(grp_var,category_variables)
  }
  
  if(transformation=="log"){
    data_i[,column_value] = log(data_i[,column_value])
  }
  
  if(transformation=="sqrt"){
    data_i[,column_value] = sqrt(data_i[,column_value])
  }
  #Identify outliers by group: ----
  out_group=data_i %>%
    dplyr::group_by_at(list_cat_var) %>%
    rstatix::identify_outliers(column_value) %>% 
    ungroup()
  if(nrow(out_group)==0){
    
  } else{
    out_group=data_i %>%
      dplyr::group_by_at(list_cat_var) %>%
      rstatix::identify_outliers(column_value) %>% 
      ungroup() %>% 
    mutate(var_analyse=column_value)
  }
  n_out_group=nrow(out_group)
  
  if(n_out_group>0){
    cat_col(paste0("- Warning, there is (are): ",n_out_group," outlier (s) that can lead to false results in ANOVA and in tests of equality of variances. (in the returne call 'outliers') \n"),"red")
  }else{
    cat_col("- No outlier found \n","green")
  }
  
  # Creation of the new column
  data_i=as.data.frame(data_i) # si un tibble in input create probleme for new column mais supprimer le level des facteurs 

    # data_i$group <- ""
    # for (i in list_cat_var) {
    #   data_i$group <- paste(data_i$group, data_i[, i], sep = "")
    # }
  
  data_i=data_i %>% 
    rowwise() %>% 
    dplyr::mutate(compile_cond = factor(paste(c_across(all_of(list_cat_var)), collapse = ""), levels = unique(paste(c_across(all_of(list_cat_var)), collapse = ""))))
  
    # data_i$group_cat_var <- ""
    # for (i in category_variables) {
    #   data_i$group_cat_var <- paste(data_i$group_cat_var, data_i[, i], sep = "")
    # }
 
  data_i=data_i %>% 
    rowwise() %>% 
    dplyr::mutate(group_cat_var = factor(paste(c_across(all_of(category_variables)), collapse = "."), levels = unique(paste(c_across(all_of(category_variables)), collapse = "."))))

  
  data_i=as.data.frame(data_i)
  data_i <- data_i %>%
    dplyr::mutate(compile_cond = factor(compile_cond, levels = unique(compile_cond))) %>% 
    dplyr::mutate(group_cat_var = factor(group_cat_var, levels = unique(group_cat_var))) 
  
    #creation of the ANOVA
    ## formula
  f <- paste0(column_value, "~", paste(sep="+",list_cat_var)) #befor
  f<-paste(column_value, "~", paste(list_cat_var, collapse = "+"))
  f<-paste(column_value, "~", paste(list_cat_var, collapse = "*")) # interaction and addition !!!
  Mx <- lm(f, data = data_i, na.action = "na.omit") #anova
  result_stat_anova <- summary(Mx)
  
  # The different test 
  pval_leven=data_i %>% 
    levene_test(as.formula(f)) %>% 
    pull(p) #test of the equality of the variance

  bar <- bartlett.test(data=data_i,data_i[,column_value]~compile_cond)$p.value # Bartlett test for combination of all variable for egality of variance 
  s_test_global=shapiro_test(residuals(Mx))$p.value ### Calculate the Shapiro-Wilk normality test on residuals
  normality_by_group=data_i %>% 
    dplyr::group_by_at(list_cat_var) %>% 
    mutate(Mean=mean(!!as.symbol(column_value))) %>% #dplyr::mutate(Mean=mean(column_value)) %>% #filter(Mean>0) %>%   #mutate(Mean=mean(!!as.symbol(column_value))
    shapiro_test(!!as.symbol(column_value)) #for the the normality of each group
  
  n_F_normality_group= normality_by_group %>% filter(p < 0.05) %>% nrow()
  if(anova_factor==1){
  # # Variance assumption ########################################################################################
    if(bar>0.05|pval_leven>0.05){
      cat_col(paste0("- Homogeneity of variances, Bartlett test (p :",disp_round(bar)," or Leven test (p :",disp_round(pval_leven),"is (are) not significant (p > 0.05). Consequently, we can assume homogeneity of variances.","\n"),color = "green")
    }else if (bar<0.05&&pval_leven<0.05){
      cat_col(paste0("- Homogeneity of variances, Bartlett (p :",disp_round(bar)," and Leven test (p :", disp_round(pval_leven), " are significant (p < 0.05). Consequently, we can't assume homogeneity of variances.", "\n"),color = "red")
    }else{
      cat_col(paste0("- !!! Homogeneity of variances error !!!, \n"),color = "red")
    }
    
    #residus du model
    
    
  }else{
    # # Variance assumption ########################################################################################
    #variance not use bartlet see introduction
    if(pval_leven>0.05){
      cat_col(paste0("- Homogeneity of variances, Levene's test is not significant (p > 0.05), is equal to ",disp_round(pval_leven), ". Consequently, we can assume homogeneity of variances in the different groups.", "\n"),color = "green")
    }else if (bar<0.05&&pval_leven<0.05){
      cat_col(paste0("- Homogeneity of variances, Levene's test is significant (p < 0.05):",disp_round(pval_leven), ". Consequently, we can't assume homogeneity of variances in the different groups. We therefore normally use a non-parametric test", "\n"),color = "red")
    }else{
      cat_col(paste0("- !!! Homogeneity of variances error !!!, \n"),color = "red")
    }
    
  }
  
  
  # Normality assumption ########################################################################################
  ## Check the normality hypothesis by analyzing the model's residuals
  ### Create a QQ plot of residuals (one for anova 1 factor and multiple for multiple factor)
  
  # pour une annova a un facteur il faut faire un shapiro sur le global mais aussi sur chaque groupe.  Idem pour les multiples facteur
  #resultats stats prints pour le global
  
    if(s_test_global<0.05){ 
    cat_col(paste0("- Normality on the model residuals, the Shapiro-Wilk test has not been validated. The p-value is significant (p ",disp_round(s_test_global),"), so we can assume that the anova is invalid because there is no normality on the model residuals.\n"),"red")
  }else{
    cat_col(paste0("- In the QQ plot, as all points lie approximately along the reference line, we can assume normality. This conclusion is supported by the Shapiro-Wilk test. The p-value is not significant (p =",disp_round(s_test_global),"), so we can assume normality. \n"),"green")
  }
  
  #resultats stats prints pour chaque groupe
  if(n_F_normality_group>0){
    cat_col(paste0("- Normality on the model residuals. Warning, there is (are): ",disp_round(n_F_normality_group)," group (s) that do not respect the normality hypothesis (in the returne call 'normality_by_group') \n"),"red")
  }else{
    cat_col("- Normality on the model residuals. Shapiro-Wilk validate for each groupe \n","green")
  }  
  
  # visualisation a l'aide des qqplot  
  print(ggqqplot(residuals(Mx))) # sur le global
  
  # ## Test the normality hypothesis by group visualisation ggqqplot----
  #
  
  if(anova_factor==1){
   print(ggqqplot(data_i,column_value,facet.by="compile_cond")) #visualisation si un seul facteur et multiple variable quali
  }else if (anova_factor==2){
    print(ggqqplot(data_i, column_value) +
            facet_grid(reformulate(list_cat_var[1],list_cat_var[2]),labeller = "label_both"))
  }else if(anova_factor==3){
    
    t1 <- c(list_cat_var[1],list_cat_var[2])
    t2 <- list_cat_var[3]
    #need to be modify for all type of value begin #####
    print(ggqqplot(data_i, column_value, ggtheme = theme_bw()) +
            facet_grid(reformulate(t1,t2),labeller = "label_both"))
    #need to be modify for all type of value end #####
  }else{
    cat_col("- Anova more than 3 factor QQplot not make yet \n","red")
    }
    
  ########### partie non parametrique ----
  Kteste <-kruskal.test(data_i[,column_value] ~ compile_cond,data = data_i) # usefull if equality of variance but not a normale distribution, to do for replace anova 1 factor !!!!
  p_value_k = Kteste$p.value
  
  ####### si on force l'anova
  
  if(force_ANOVA==T){
    cat_col("Attention you have forced the test to make the ANOVA. If you're a biologist, this may be understandable.But be careful with interpretations \n",color = "yellow")
    cat_col("Here are the ANOVA results: \n",color = "red")
    print(result_stat_anova)
    cat_col("- Now perform a posthoc test: \n",color = "magenta")
    cat_col(paste0("- If 1 factor we can do Kruskal-Wallis test, normali note if we performe a test for more thant 1 factor. Here equal to ",anova_factor," factor. Result of ",disp_round(p_value_k)," for the Kruskal-Wallis test\n"),"yellow")
    test_output="anova" # si on force a faire de l'anova ou si hypothese validite ok alors anova sinon non parametrique donc kruskal donc mann-whitney ou pairwise.t.test (voir graph compare more than two groups)
    
  }else{
    if (s_test_global<0.05 | n_F_normality_group > 0 | pval_leven<0.05 | p_value_k <0.05 | bar <0.05 ) {
      cat_col("You need to use a non-parametric test or transform variables if you whant performe ANOVA test !!! \n",color = "red")
      test_output="kruskal" # si on force a faire de l'anova ou si hypothese validite ok alors anova sinon non parametrique donc kruskal donc mann-whitney ou pairwise.t.test (voir graph compare more than two groups) 
      variable_transformation_requirements=T
      
      }else{
        cat_col("Here are the ANOVA results: \n",color = "green")
        print(result_stat_anova)
        test_output="anova" # si on force a faire de l'anova ou si hypothese validite ok alors anova sinon non parametrique donc kruskal donc mann-whitney ou pairwise.t.test (voir graph compare more than two groups) 
        cat_col("- Now perform a posthoc-test: \n",color = "magenta") 
        
        if(anova_factor==1){
          cat_col(paste0("- ANOVA 1 factor so we can do Kruskal-Wallis test. The result is",disp_round(p_value_k)," for the Kruskal-Wallis test\n"),"green")
          
        }else{
          cat_col("- ANOVA more than 1 factor so we can do Kruskal-Wallis test. Need to do LOESS test (not making for the moment) \n","yellow")
        }
        
    }
  }
  return(list(outliers=out_group,
              normality_by_group=normality_by_group,
              result_stat_anova=result_stat_anova,
              variable_transformation_requirements=variable_transformation_requirements,
              data_used=data_i,
              test_output=test_output #all test performe att this end (parametric or not)
              ))
  
}

# ############ non parametrique faires loes pour plus d'un facteur
# #
# # Specify the parameters.
# #
# 
# 
# beta <- rbind("Genotype 1" = c(1, 1/4, 1/8, 1),
#               "Genotype 2" = c(1/2, 3/4, 1/2, 7/4),
#               "Genotype 3" = c(3/2, 1/2, 1/4, 1/4),
#               "Genotype 4" = c(2, 5/4, 1/4, 3/4))
# colnames(beta) <- paste("Treatment", 1:4)
# sigma <- 1/4
# tau <- 1/4
# rate <- 2
# #
# # Simulate.
# #
# ff <- function(x) exp(-abs(x))
# # ff <- function(x) (1 + x^2)^(-1)
# rf <- function(x, beta, tau, sigma, rate) {
#   y <- beta * ff(rate * x)
#   rnorm(length(y), y, sqrt(sigma^2 + y * tau^2))
# }
# 
# n.per.group <- 120*2
# set.seed(17)
# x <- seq(-2, 2, length.out=n.per.group)
# X <- expand.grid(Distance = x,
#                  Genotype = rownames(beta),
#                  Treatment = colnames(beta))
# X$Response <- c(sapply(c(beta), function(b) rf(x, b, tau, sigma, rate)))
# #
# # Plot.
# #
# library(ggplot2)
# ggplot(X, aes(abs(Distance), Response, color = Treatment)) + 
#   geom_hline(yintercept=0) + 
#   geom_point(alpha = 1/4, show.legend = FALSE) + 
#   stat_smooth(method = "loess", formula = "y ~ x", size=1.5, 
#               show.legend = FALSE, se=FALSE) + 
#   labs(x = "Distance", y = "Response") + 
#   theme_classic () +
#   facet_grid(Genotype ~ Treatment) + 
#   ggtitle("Data", "(Randomly Generated)")
# 
# ggplot(X, aes(Distance, Response, color = Treatment)) + 
#   geom_hline(yintercept=0) + 
#   geom_point(alpha = 1/4, show.legend = FALSE) + 
#   stat_smooth(method = "loess", formula = "y ~ x", size=1.5, show.legend = FALSE) + 
#   labs(x = "Distance", y = "Response") + 
#   theme_classic () +
#   facet_grid(Genotype ~ Treatment) + 
#   ggtitle("Data", "(Randomly Generated)")
# #
# #
# # Fit.
# #
# f <- function(theta) {
#   beta <- matrix(exp(theta[1:16]), 4, 
#                  dimnames=list(paste("Genotype", 1:4), 
#                                paste("Treatment", 1:4)))
#   tau <- exp(theta[17])
#   sigma <- exp(theta[18])
#   rate <- exp(theta[19])
#   
#   y <- beta[cbind(X$Genotype, X$Treatment)] * ff(X$Distance * rate)
#   -sum(dnorm(X$Response, y, sqrt(sigma^2 + y * tau^2), log = TRUE))
# }
# 
# theta <- rep(0, 19)
# fit <- nlm(f, theta, hessian=TRUE)
# 
# theta.hat <- fit$estimate
# beta.hat <- matrix(exp(theta.hat[1:16]), 4, 
#                    dimnames=list(paste("Genotype", 1:4), 
#                                  paste("Treatment", 1:4)))
# tau.hat <- exp(theta.hat[17])
# sigma.hat <- exp(theta.hat[18])
# rate.hat <- exp(theta.hat[19])
# #
# # Plot.
# #
# X$Prediction <- c(beta.hat[cbind(X$Genotype, X$Treatment)] * ff(X$Distance * rate.hat))
# ggplot(X, aes(Prediction, Response)) +
#   geom_point(alpha=1/2) + 
#   geom_abline(intercept=0, slope=1, size=1.5, color="#d01010") + 
#   ggtitle("Model Fit")
# 
# ggplot(X, aes(abs(Distance), Prediction, color = Treatment))  + 
#   geom_hline(yintercept=0) + 
#   geom_point(aes(y=Response), alpha=1/4, show.legend=FALSE) + 
#   geom_line(size=1.5, show.legend=FALSE) + 
#   labs(x = "Distance", y = "Response") + 
#   theme_classic () +
#   facet_grid(Genotype ~ Treatment) + 
#   ggtitle("Fitted Values", "(With Original Data)")
# #
# # Extract the covariance matrix.
# #
# V <- solve(fit$hessian)
# s <- rbind(Actual = c(beta, tau, sigma, rate),
#            Fit = exp(fit$estimate),
#            SE = sqrt(diag(V)),
#            t = (fit$estimate - log(c(beta, tau, sigma, rate))) / sqrt(diag(V)))
# colnames(s) <- c(c(outer(rownames(beta), colnames(beta), paste)), "tau", "sigma", "rate")
# print(t(s), digits=3)
# 
# colnames(V) <- rownames(V) <- colnames(s)
# #
# # Check the correlation among the estimates.
# #
# Sigma <- t(V/sqrt(diag(V))) / sqrt(diag(V))
# # Shows little correlation among the betas; negative correlation between tau and sigma
# image(seq_along(colnames(Sigma)), seq_along(rownames(Sigma)), Sigma,
#       main="Variance-Covariance Matrix of Estimates")
# 
# Sigma <- pmin(1, pmax(-1, Sigma))
# h <- hclust(as.dist(matrix(acos(Sigma), nrow(V))), method="median")
# plot(h)
# 


# SLOAP pour faire une sorte d'anova non parametrique (mais tres tres complexe)
#Installation des packages si nécessaire
# install.packages("car")
# install.packages("agricolae")



