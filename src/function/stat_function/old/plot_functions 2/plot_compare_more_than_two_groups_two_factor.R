# plot for more than two groupe (witouth intergroup)
library(ggplot2)
library(multcompView)
library(ggrepel) #pour bien positioner le nom des outliers
library(ggpattern)

plot_one_variable_one_comparison=function(
    Y,
    Z,
    X,
    Xlab="x",
    Ylab="y",
    Tukey=T,
    Y_bis,
    outlier_show=T,
    label_outlier="plant_num",
    category_variables=category_variables,
    grp_var,
    compile_cond,
    group_cat_var,
    hex_pallet,
    strip_normale=F
    )
{
  if(Tukey==T){
    cat_col("For letter we use Tukey method \n","yellow")
  }else{
    cat_col("For letter we use pairwise.t.test, with correction holm \n","yellow")
  }
    #for verif normalité des résidus et homogénéité des variances
  mod2<-aov(Y~compile_cond)
  
  #anova_for_graph
  #calcule lettre significativiter ---- 
  #normalement ne fonctionne que si le teste de kruscal walis est significatif
  Kteste=kruskal.test(Y ~compile_cond, data=Z) # YES !!!!
  # réalisation des comparaisons deux à deux des moyennes par des tests de Wilcoxon (pool.sd=F), avec ajustement des pvalues par la méthode de Holm
  # et récupération de la matrice triangulaire des p-values
  
  if(Tukey==F){
    cat_col("pairwise.t.test with holm correction \n","green")
    pp <- pairwise.t.test(Y, compile_cond,p.adjust.method ="holm",paired = F,pool.sd=F) #voire ici pourquoi methode holm et pourquoi on a besoins d'ajuster les pvalue (pour diminuer les faux positif autrement dit le risque derreur) https://delladata.fr/comparaisons-multiples-et-ajustement-des-pvalues-avec-le-logiciel-r/  #si vraiment je veux avoir les resultats de student je dois a la place de holm ecrire none et utiliser pool.sd (car je ne pas fait de test sur l'egalité de variance). (Simply put, if you are making multiple comparisons, you are increasing the chances of finding spurious results. A correction adjusts for this.)
    # réalisation des comparaisons deux à deux des moyennes par des tests de Wilcoxon (ou man_witney),(pas besoins de verif l’hypothèse d’homoscedasticité ou d’égalité des variances) avec ajustement des pvalues par la méthode de Holm. Foncionne bien sur les petits echantillon. https://jonathanlenoir.files.wordpress.com/2013/12/tests-de-comparaison-de-moyennes-non-param.pdf / https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html
    
    #voir ici pourquoi j'ai choisis pool.sd=F
    #https://stats.stackexchange.com/questions/337749/in-r-how-can-the-p-value-during-pairwise-t-tests-without-adjustment-be-differen
    
    pp
    # transformation de la matrice triangulaire de p-values en une matrice rectangulaire
    mymat <-tri.to.squ(pp$p.value)
    mymat
    
    # création des lettres correspondant à chaque moyenne
    myletters <- multcompLetters(mymat,compare="<=",threshold=0.05,Letters=letters)
    myletters
    # conserver les lettres dans un data.frame
    myletters_df <- data.frame(group=names(myletters$Letters),letter = myletters$Letters)
  }else{
    cat_col("Tukey processing \n","red")
    tuktuk<-TukeyHSD(mod2, ordered=FALSE, conf.level=0.95, na.rm=T)
    multcompLetters4(mod2, tuktuk)
    myletters=multcompLetters4(mod2, TukeyHSD(mod2,ordered=FALSE, conf.level=0.95, na.rm=T))
    myletters_df=as.data.frame(rownames(myletters$compile_cond$LetterMatrix))
    for (i in 1:length(myletters$compile_cond$Letters)){
      #print(myletters$X$Letters[[i]])
      myletters_df$letter[i]=myletters$compile_cond$Letters[[i]] 
    }
    
    colnames(myletters_df)=c("groupe","letter")
    rownames(myletters_df)=myletters_df$groupe
  }
  #Afficher les résultats des différences à l’aide de lettres   
  # plot
  #faire un merge pour que ça fonctionne
  Z$group=NA
  for(i in 1:length(Y)){
    for(j in 1:length(myletters_df$group)){
      if(compile_cond[i]==myletters_df$group[j]){
        Z$group[i]=as.character(myletters_df$letter[j])
      }else{}
    }
  }

    #avt  
  #######################
  # row.names(Z)=as.data.frame(Z[,label_outlier])[,1] #avant label outlier
  # dat <- Z %>% tibble::rownames_to_column(var="outlier") %>% group_by_at(vars({{category_variables}})) %>% mutate(is_outlier=ifelse(is_outlier(get(Y_bis)), Y, as.numeric(NA)))
  # dat$outlier[which(is.na(dat$is_outlier))] <- as.numeric(NA)
  # Z$outlier=dat$outlier
  # ########################
  # 
  
  #######################
  #print(as.data.frame(Z[,label_outlier])[,1]) ###################
  row.names(Z)=as.data.frame(Z[,label_outlier])[,1] #avant label outlier
  dat <- Z %>% tibble::rownames_to_column(var="outlier") %>% group_by(compile_cond) %>% mutate(is_outlier=ifelse(is_outlier(get(Y_bis)), Y, as.numeric(NA)))
  dat$outlier[which(is.na(dat$is_outlier))] <- as.numeric(NA)
  Z$outlier=dat$outlier
  #Z$column_wrap=Z[,column_wrap]
  
  ########################
  #graph
  if(grp_var==""){
    if(Xlab=="condition"){
      Z=Z%>% mutate(group_cat_var=factor(group_cat_var,levels=c("Sto_WW_OT","Sto_WS_OT","Sto_WW_HS","Sto_WS_HS","Wen_WW_OT","Wen_WS_OT","Wen_WW_HS","Wen_WS_HS")))
    }
    if(Xlab=="Compartment"|Xlab=="compartment_l"){
      Z=Z%>% mutate(group_cat_var=factor(group_cat_var,levels=c("Phylloplane","Leaf endosphere" , "Root endosphere" , "Rhizoplane","Rhizosphere" ,"Bulk soil")))
    } 
    if(Xlab=="climat_condition"){
      Xlab="Treatment"
      Z=Z%>% mutate(group_cat_var=factor(group_cat_var,levels=c("WW_OT","WS_OT" , "WW_HS" , "WS_HS")))
      
    }
    
  p=ggplot(Z, aes(x=group_cat_var, y=Y, colour=group_cat_var, fill=group_cat_var))+
    geom_boxplot(outlier.alpha = 0, alpha=0.25)+
    geom_jitter(alpha=0.5,position = position_jitter(seed = 1),shape=16)+ 
    #geom_text(aes(label=outlier),na.rm=TRUE,nudge_y=0.05)+
    #stat_summary(fun.y=mean, colour="black", geom="point", 
    #            shape=18, size=3) +
    my_theme+
    theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+
    geom_text(data = Z, aes(label = group, y = 1.05*max(Y[!is.na(Y)])))+xlab(Xlab)+ylab(Ylab)+
    labs(caption=paste0("kruskal test=",format.pval(Kteste$p.value,digits = 2)))+
    scale_color_manual(values=hex_pallet)+
    scale_fill_manual(values=hex_pallet)
  }else{
    if(grp_var=="genotype"){
      if(Xlab=="condition"){
        Z=Z%>% mutate(group_cat_var=factor(group_cat_var,levels=c("Sto_WW_OT","Sto_WS_OT","Sto_WW_HS","Sto_WS_HS","Wen_WW_OT","Wen_WS_OT","Wen_WW_HS","Wen_WS_HS")))
      }
      if(Xlab=="Compartment"|Xlab=="compartment_l"){
        Z=Z%>% mutate(group_cat_var=factor(group_cat_var,levels=c("Phylloplane","Leaf endosphere" , "Root endosphere" , "Rhizoplane","Rhizosphere" ,"Bulk soil")))
      } 
      if(Xlab=="climat_condition"){
        Xlab="Treatment"
        Z=Z%>% mutate(group_cat_var=factor(group_cat_var,levels=c("WW_OT","WS_OT" , "WW_HS" , "WS_HS")))
        
      }
    }
    if(strip_normale==F){
    p=ggplot(Z, aes(x=group_cat_var, y=Y, colour=group_cat_var, fill=group_cat_var))+#geom_boxplot()+facet_grid(~genotype)
      geom_boxplot(outlier.alpha = 0, alpha=0.25)+
      geom_jitter(alpha=0.5,position = position_jitter(seed = 1),shape=16)+ 
      #geom_text(aes(label=outlier),na.rm=TRUE,nudge_y=0.05)+
      #stat_summary(fun.y=mean, colour="black", geom="point", 
      #            shape=18, size=3) +
      theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+
      geom_text(data = Z, aes(label = group, y = 1.06*max(Y[!is.na(Y)])))+
      xlab(Xlab)+
      ylab(Ylab)+
      labs(caption=paste0("kruskal test=",format.pval(Kteste$p.value,digits = 2)))+
      scale_fill_manual(values=hex_pallet)+
      scale_color_manual(values=hex_pallet)+
      theme(legend.key.size = unit(1.2, 'cm'))+
      scale_pattern_manual(values=c('none', "stripe")) + # invert here to apply strip to the good condition
      scale_pattern_type_manual(values=c(NA, NA))+
      #facet_grid(~genotype)
      facet_grid(~genotype,scales="free",switch = "x")+
      theme(panel.grid = element_blank(),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.border = element_rect(fill = NA, colour = "white"), 
            axis.line = element_line(),
            strip.background = element_blank(),
            panel.spacing = unit(0.6, "lines"), #peut creee des eerreur je ne sais pas pourquoi. Avant 0.7. maintenant 0.6
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())
      # my_theme+
    }else{
      p=ggplot(Z, aes(x=group_cat_var, y=Y, colour=group_cat_var, fill=group_cat_var))+#geom_boxplot()+facet_grid(~genotype)
        geom_text(data = Z, aes(label = group, y = 1.06*max(Y[!is.na(Y)])))+
        geom_boxplot(outlier.alpha = 0, alpha=0.25)+
        geom_jitter(alpha=0.5,position = position_jitter(seed = 1),shape=16)+ 
        #geom_text(aes(label=outlier),na.rm=TRUE,nudge_y=0.05)+
        #stat_summary(fun.y=mean, colour="black", geom="point", 
        #            shape=18, size=3) +
        theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+
        xlab(Xlab)+
        ylab(Ylab)+
        labs(caption=paste0("kruskal test=",format.pval(Kteste$p.value,digits = 2)),
             color=Xlab,fill=Xlab)+
        scale_fill_manual(values=hex_pallet)+
        scale_color_manual(values=hex_pallet)+
        #theme(legend.key.size = unit(1.2, 'cm'))+
        #scale_pattern_manual(values=c('none', "stripe")) + # invert here to apply strip to the good condition
        #scale_pattern_type_manual(values=c(NA, NA))+
        #facet_grid(~genotype)
        facet_grid(~genotype,scales="free")+
        theme(panel.grid = element_blank(),
              panel.background = element_rect(fill = "white", colour = "white"), 
              panel.border = element_rect(fill = NA, colour = "black"), 
              strip.background = element_rect(color = "black")
              #axis.line = element_line()#,
              #strip.background = element_blank(),
              #panel.spacing = unit(0.6, "lines"), #peut creee des eerreur je ne sais pas pourquoi. Avant 0.7. maintenant 0.6
              #axis.title.x=element_blank(),
              #axis.text.x=element_blank(),
              #axis.ticks.x=element_blank()
              )
    }
  }
  if(outlier_show==T){
    p=p+ggrepel::geom_label_repel(aes(label = outlier),fill = 'white',position = position_jitter(seed = 1),na.rm=TRUE,box.padding   = 0.35, point.padding = 0.5,segment.color = 'grey50',max.overlaps = 100)
    print(p)
    return(p)
  }else{
    print(p)
    
    return(p)
  }
}
