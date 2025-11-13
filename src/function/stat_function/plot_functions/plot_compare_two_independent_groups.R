#script for plot two mean independent groups
library(ggpubr)
library(ggplot2)

plot_compare_two_independent_groups<-function(data,result_test=summary_result_stat[["result_test"]],x_val,y_val,method,hex_pallet=c("#EF5B0C","#003865"),Ylab_i=Ylab_i
){
  p<-ggplot(data, aes_string(x = x_val, y = y_val,col=x_val))+
              geom_boxplot(outlier.shape =NA )+
              geom_jitter()+
              scale_color_manual(values=hex_pallet)
               #              color = x_val, palette = hex_pallet,
               #              add = "jitter")+
  # p<-ggboxplot(data, x = x_val, y = y_val,
  #              color = x_val, palette = hex_pallet,
  #              add = "jitter")+
  
#for captation
if(method=="welch"){
  #captation_x=paste0("Semi-parametric test: ",disp_round(result_test$result_t_test$p.value))
  captation_x <- paste0("Semi-parametric test (Welch Test): ",disp_round(result_test$result_welch_test$p.value),"
                        \n (Parametric test (T.test):  ",disp_round(result_test$result_t_test$p.value),")",
                        "\n (Nonparametric test (Wilcox test): ",disp_round(result_test$result_u_test$p.value),")")
            # stat_compare_means(method = "t.test",method.args = list(var.equal = FALSE),label.x = 1.5, label.y = max(data[,y_val]))
}else if (method=="t.test"){
  captation_x <- paste0("\n Parametric test (T.test):  ",disp_round(result_test$result_t_test$p.value),
                        "\n (Nonparametric test (Wilcox test): ",disp_round(result_test$result_u_test$p.value),")")
}else if (method=="wilcox.test"){
  captation_x <- paste0("\n (Parametric test (T.test):  ",disp_round(result_test$result_t_test$p.value),")",
                        "\n Nonparametric test (Wilcox test): ",disp_round(result_test$result_u_test$p.value))
}
  
  
p=p+my_theme+ labs(caption = captation_x)+ylab(Ylab_i)
print(p)
return(p)
}

