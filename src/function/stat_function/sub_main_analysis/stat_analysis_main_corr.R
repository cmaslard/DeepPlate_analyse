if (category_variables == "" && grp_var == "") {
  cat_col("(These tests are not currently available.) !!!!Correlation without categorical variable!!!\n", color = "red")
  break
} else if (category_variables != "" && grp_var == "") {
  cat_col("(These tests are not currently available.) Correlation with one variable categorical (ANCOVA) \n", color = "red")
  break
} else {
  cat_col("(error we haven't this possibility, see the code \n", color = "red")
  break
} # end correlation
