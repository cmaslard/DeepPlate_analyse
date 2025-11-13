Contrasts <- function(model, Target, Replicate, Interaction) {
  
  # ## Subdirectory
  #   if (!I(Project_Name %in% dir(Results_Directory))) dir.create(paste0(Results_Directory,"/",Project_Name), showWarnings=FALSE)
  #   subdirectory <- "DiffAnalysis"
  #   path=paste0(Results_Directory,"/",Project_Name,"/",subdirectory)
  #   if (!I(subdirectory %in% dir(paste0(Results_Directory,"/",Project_Name)))) dir.create(path, showWarnings=FALSE)
  
  ## Nb of biological factors
  NbBioFactors <- ncol(Target)-1
  
  ## Contrasts
  FacBio <- 1:(NbBioFactors)
  coeff.name = colnames(model)
  nl = sapply(Target[, FacBio], function(x)
    length(levels(x)))
  nl1 <- levels(Target[, FacBio[1]])
  contrast.names <- as.character("")
  
  ## Number of biological factor == 1
    if (NbBioFactors == 1)
    {
        contrast.factor1 <- matrix(0, ncol = length(coeff.name), nrow = 1)
        for (h in 1:(length(nl1) - 1))
        {
            for (i in (h + 1):length(nl1))
            {
                f2 <- paste(colnames(Target)[FacBio][1], nl1[i], sep = "")
                contrast.definition <- rep(0, length(coeff.name))
                f1 <- paste(colnames(Target)[FacBio][1], nl1[h], sep = "")
                if (h != 1)
                {
                    contrast.definition[which(f1 == coeff.name)] = 1
                }
                contrast.definition[which(f2 == coeff.name)] = (-1)
                contrast.factor1 = rbind(contrast.factor1, contrast.definition)
                contrast.names = c(contrast.names, paste0("[", nl1[h], "-", nl1[i], "]"))
            }
        }
        contrast.matrix <- contrast.factor1
        colnames(contrast.matrix) = coeff.name
        rownames(contrast.matrix) = contrast.names
    }
  
  ## Number of biological factor == 2
  if (NbBioFactors == 2) {
    if (Interaction == FALSE) {
      ## Factor 1
      contrast.factor1 <- matrix(0, ncol = length(coeff.name), nrow = 1)
      for (h in 1:(length(nl1) - 1))
      {
        for (i in (h + 1):length(nl1))
        {
          f2 <- paste(colnames(Target)[FacBio][1], nl1[i], sep = "")
          contrast.definition <- rep(0, length(coeff.name))
          f1 <- paste(colnames(Target)[FacBio][1], nl1[h], sep = "")
          if (h != 1)
          {
            contrast.definition[which(f1 == coeff.name)] = 1
          }
          contrast.definition[which(f2 == coeff.name)] = (-1)
          contrast.factor1 = rbind(contrast.factor1, contrast.definition)
          contrast.names = c(contrast.names, paste0("[", nl1[h], "-", nl1[i], "]"))
        }
      }
      ## Factor 2
      nl2 <- levels(Target[, FacBio[2]])
      contrast.factor2 <- matrix(0, ncol = length(coeff.name), nrow = 1)
      for (j in 1:(length(nl2) - 1))
      {
        for (k in (j + 1):length(nl2))
        {
          f2 <- paste(colnames(Target)[FacBio][2], nl2[k], sep = "")
          contrast.definition <- rep(0, length(coeff.name))
          f1 <- paste(colnames(Target)[FacBio][2], nl2[j], sep = "")
          if (j != 1)
          {
            contrast.definition[which(f1 == coeff.name)] = 1
          }
          contrast.definition[which(f2 == coeff.name)] = (-1)
          contrast.factor2 = rbind(contrast.factor2, contrast.definition)
          contrast.names = c(contrast.names, paste0("[", nl2[j], "-", nl2[k], "]"))
        }
      }
      contrast.matrix <- unique(rbind(contrast.factor1, contrast.factor2))
      colnames(contrast.matrix) = coeff.name
      rownames(contrast.matrix) = contrast.names
    }
    
    if ( Interaction == TRUE){
      nl2 <- levels(Target[, FacBio[2]])
      contrast.Interaction<-matrix(0,ncol=length(coeff.name),nrow=1)
      for (h in 1:(length(nl1)-1))
      {
        f1<-paste(colnames(Target)[FacBio][1],nl1[h],sep="")
        for (i in (h+1):length(nl1))
        {
          f2<-paste(colnames(Target)[FacBio][1],nl1[i],sep="")
          
          for(j in 1:(length(nl2)-1))
          {
            g1<-paste(colnames(Target)[FacBio][2],nl2[j],sep="")
            for (k in (j+1):length(nl2))
            {
              g2<-paste(colnames(Target)[FacBio][2],nl2[k],sep="")
              contrast.definition<-rep(0,length(coeff.name))
              contrast.definition[grep(paste(f1,g1,sep=":"),coeff.name)]=
                ifelse(is.element(paste(f1,g1,sep=":"),coeff.name),1,0)
              contrast.definition[grep(paste(f1,g2,sep=":"),coeff.name)]=
                ifelse(is.element(paste(f1,g2,sep=":"),coeff.name),-1,0)
              contrast.definition[grep(paste(f2,g1,sep=":"),coeff.name)]=
                ifelse(is.element(paste(f2,g1,sep=":"),coeff.name),-1,0)
              contrast.definition[grep(paste(f2,g2,sep=":"),coeff.name)]=
                ifelse(is.element(paste(f2,g2,sep=":"),coeff.name),1,0)
              contrast.Interaction=rbind(contrast.Interaction,contrast.definition)
              contrast.names=c(contrast.names,paste0("[",nl1[h],"_",nl2[j],"-",nl1[h],"_",nl2[k],"]-[",nl1[i],"_",nl2[j],"-",nl1[i],"_",nl2[k],"]"))
            }
          }
        }
      }
      colnames(contrast.Interaction)=coeff.name
      rownames(contrast.Interaction)=contrast.names
      
      ## effect of biological factor 1 averaged on biological factor 2
      contrast.factor1<-matrix(0,ncol=length(coeff.name),nrow=1)
      contrast.names <- as.character("")
      for (h in 1:(length(nl1)-1))
      {
        for (i in (h+1):length(nl1))
        {
          f2<-paste(colnames(Target)[FacBio][1],nl1[i],sep="")
          contrast.definition<-rep(0,length(coeff.name))
          f1<-paste(colnames(Target)[FacBio][1],nl1[h],sep="")
          if(h!=1)
          {
            contrast.definition[which(f1==coeff.name)]=1
            contrast.definition[grep(paste0(f1,":"),coeff.name)]=(1/nl[2])
          }
          contrast.definition[which(f2==coeff.name)]=(-1)
          contrast.definition[grep(paste0(f2,":"),coeff.name)]=(-1/nl[2])
          contrast.factor1=rbind(contrast.factor1,contrast.definition)
          contrast.names=c(contrast.names,paste0("[",nl1[h],"-",nl1[i],"]"))
        }
      }
      
      colnames(contrast.factor1)=coeff.name
      rownames(contrast.factor1)=contrast.names
      
      
      ## effect of biological factor 2 averaged on biological factor 1
      contrast.factor2<-matrix(0,ncol=length(coeff.name),nrow=1)
      contrast.names <- as.character("")
      for (j in 1:(nl[2]-1))
      {
        for (k in (j+1):nl[2])
        {
          contrast.definition<-rep(0,length(coeff.name))
          g1<-paste(colnames(Target)[FacBio][2],nl2[j],sep="")
          if(j!=1)
          {
            contrast.definition[which(g1==coeff.name)]=1
            contrast.definition[grep(paste0(":",g1),coeff.name)]=(1/nl[1])
          }
          g2<-paste(colnames(Target)[FacBio][2],nl2[k],sep="")
          contrast.definition[which(g2==coeff.name)]=(-1)
          contrast.definition[grep(paste0(":",g2),coeff.name)]=(-1/nl[1])
          contrast.factor2=rbind(contrast.factor2,contrast.definition)
          contrast.names=c(contrast.names,paste0("[",nl2[j],"-",nl2[k],"]"))
        }
      }
      
      colnames(contrast.factor2)=coeff.name
      rownames(contrast.factor2)=contrast.names
      
      
      ## effect of biological factor 2 given one level of biological factor 1
      contrast.factor2sachant1<-matrix(0,ncol=length(coeff.name),nrow=1)
      contrast.names <- as.character("")
      for (h in 1:nl[1])
      {
        h1<-paste(colnames(Target)[FacBio][1],nl1[h],sep="")
        for (j in 1:(nl[2]-1))
        {
            for (k in (j+1):nl[2])
          {
              contrast.definition<-rep(0,length(coeff.name))
              g1<-paste(colnames(Target)[FacBio][2],nl2[j],sep="")
              if(j!=1)
              {
                  contrast.definition[which(g1==coeff.name)]=1
                  if(h!=1)
                      contrast.definition[grep(paste0(h1,":",g1),coeff.name)]=1
              }
              g2<-paste(colnames(Target)[FacBio][2],nl2[k],sep="")
              contrast.definition[which(g2==coeff.name)]=-1
              if(h!=1)
                  contrast.definition[grep(paste0(h1,":",g2),coeff.name)]=-1
              contrast.factor2sachant1=rbind(contrast.factor2sachant1,contrast.definition)
              contrast.names=c(contrast.names,paste0("[",nl1[h],"_",nl2[j],"-",nl1[h],"_",nl2[k],"]"))
          }
        }
      }
      colnames(contrast.factor2sachant1)=coeff.name
      rownames(contrast.factor2sachant1)=contrast.names
      
      ## effect of biological factor 1 given one level of biological factor 2
      contrast.factor1sachant2<-matrix(0,ncol=length(coeff.name),nrow=1)
      contrast.names <- as.character("")
      for (h in 1:nl[2])
      {
          h1<-paste(colnames(Target)[FacBio][2],nl2[h],sep="")
          for (j in 1:(nl[1]-1))
          {
              for (k in (j+1):nl[1])
              {
                  contrast.definition<-rep(0,length(coeff.name))
                  g1<-paste(colnames(Target)[FacBio][1],nl1[j],sep="")
                  if(j!=1)
                  {
                      contrast.definition[which(g1==coeff.name)]=1
                      if(h!=1)
                          contrast.definition[grep(paste0(g1,":",h1),coeff.name)]=1
                  }
                  g2<-paste(colnames(Target)[FacBio][1],nl1[k],sep="")
                  contrast.definition[which(g2==coeff.name)]=(-1)
                  if(h!=1)
                      contrast.definition[grep(paste0(g2,":",h1),coeff.name)]=(-1)
                  contrast.factor1sachant2=rbind(contrast.factor1sachant2,contrast.definition)
                  contrast.names=c(contrast.names,paste0("[",nl2[h],"_",nl1[j],"-",nl2[h],"_",nl1[k],"]"))
              }
          }
      }
      colnames(contrast.factor1sachant2)=coeff.name
      rownames(contrast.factor1sachant2)=contrast.names
      
      contrast.matrix <- unique(rbind(contrast.factor1, contrast.factor2, contrast.factor1sachant2, contrast.factor2sachant1, contrast.Interaction ))
    }
  }
    if(nrow(contrast.matrix)==2)
    {
        contrast.matrix <- as.data.frame(contrast.matrix)
        contrast.matrix = contrast.matrix[-1,]
    } else
    {
        contrast.matrix <- contrast.matrix[-1,]
    }
    
    # Output <- file(paste0(path,"/","GLM_Contrasts.txt"), open="wt")
    # sink(Output)
    # sink(Output, type = "message")
    # cat("################################################\n")
    # cat("Names of the contrasts \n")
    # cat("################################################\n\n")
    # Contrast_Names <- data.frame(row.names(contrast.matrix))
    # colnames(Contrast_Names) <- c("Contrast_Names")
    # print(Contrast_Names)
    # sink(type = "message")
    # sink()
    
    ## Write contrast.matrix in table file
    # if(nrow(contrast.matrix)==1)
    # {
    #     contrast.matrix[1,1]=rownames(contrast.matrix)
    #     fileout<- paste0(path,"/","Contrasts_Matrix.txt")
    #     write.table(contrast.matrix,fileout,row.names=FALSE,col.names=TRUE, sep="\t", quote = FALSE)
    # }
    #     else
    # {
    #     contrast.matrix <- as.data.frame(contrast.matrix[,-1])
    #     contrast.matrix <- tibble::rownames_to_column(contrast.matrix, "Contrasts")
    #     fileout<- paste0(path,"/","Contrasts_Matrix.txt")
    #     write.table(contrast.matrix,fileout,row.names=FALSE,col.names=TRUE, sep="\t", quote = FALSE)
    # }   
   
 # Model <- list("GLM_Model" = glm_model, "Contrasts" = contrast.matrix)
  return(contrast.matrix)
}
