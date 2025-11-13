# Corentin Maslard
# Date: 20250220

# The aim of this scipt is to automatically perform GOs on different groups, and then to compare the different groups with each other. 

# pkg 
library(org.Psativum1c.eg.db) # if not working install it
library(clusterProfiler)
library(GO.db)
library(AnnotationDbi)
library(dplyr)
library(parallel)
library(qvalue)
# function
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

run_enrichment <- function(group) {
  library(tidyverse)
  library(clusterProfiler)
  # Extraire la liste des gènes pour le groupe
  gene_list <- group_info %>% 
    filter(group == !!group) %>% 
    pull(ID)
  
  # Analyse d'enrichissement avec clusterProfiler
  enrichment_results <- enricher(gene = gene_list, 
                                 TERM2GENE = all_GO_type, 
                                 pvalueCutoff = pvalueCutoff_i,
                                 pAdjustMethod = pAdjustMethod_i,
                                 minGSSize = minGSSize_i,
                                 maxGSSize = maxGSSize_i,
                                 universe = universe_i)
  
  # Retourner les résultats s'ils ne sont pas vides
  if (!is.null(enrichment_results) && nrow(as.data.frame(enrichment_results)) > 0) {
    enrichment_df <- as.data.frame(enrichment_results)
    enrichment_df$group <- group  # Ajouter le nom du groupe
    return(enrichment_df)
  }
  
  return(NULL)
}

run_enrichment_manual <- function(group) {
  
  library(dplyr)
  library(qvalue)
  # library(GO.db) # si tu veux ensuite ajouter les descriptions

  # ---------------------------------------------------------------------------- #
  # 3) Extraire la liste de gènes (p.ex. DEG) pour ce groupe
  #    en s'assurant de ne prendre que ceux présents dans all_GO_type_filter
  # ---------------------------------------------------------------------------- #
  gene_list <- group_info %>% 
    filter(group == !!group) %>% 
    filter(ID %in% unique(all_GO_type_filter$gene_id)) %>% 
    pull(ID)
  
  # Si aucun gène, on renvoie un data.frame vide
  if (length(gene_list) == 0) {
    return(NULL)
  }
  
  # ---------------------------------------------------------------------------- #
  # 4) L'ensemble (geneset) = (GO,gène) restreint aux gènes de ce groupe
  # ---------------------------------------------------------------------------- #
  geneset <- all_GO_type_filter %>%
    filter(gene_id %in% gene_list)
  
  # Si aucun GO retenu, on renvoie un data.frame vide
  if (nrow(geneset) == 0) {
    return(NULL)
  }
  
  # ---------------------------------------------------------------------------- #
  # 5) Créer un data.frame de base pour héberger nos résultats
  # ---------------------------------------------------------------------------- #
  result_subgroup <- data.frame(
    ID            = unique(geneset$go_id),
    Description   = unique(geneset$go_id),  # On écrira plus tard le vrai nom si dispo
    GeneRatio     = NA_character_,  # On va stocker du texte "x/y"
    BgRatio       = NA_character_,
    RichFactor    = NA_real_,
    FoldEnrichment= NA_real_,
    zScore        = NA_real_,
    pvalue        = NA_real_,
    p.adjust      = NA_real_,
    qvalue        = NA_real_,
    geneID        = NA_character_,  # si tu veux y mettre la liste des gènes
    Count         = NA_real_,
    group         = group,
    stringsAsFactors = FALSE
  )
  
  # ---------------------------------------------------------------------------- #
  # 6) Boucler sur chaque GO du groupe et calculer Count, ratios, p-value...
  # ---------------------------------------------------------------------------- #
  for (GO_i in seq_along(result_subgroup$ID)) {
    
    current_GO <- result_subgroup$ID[GO_i]
    
    # Nombre de gènes (dans ce groupe) associés à ce GO
    count_g <- sum(geneset$go_id == current_GO)
    result_subgroup[GO_i,"Count"] <- count_g
    
    # GeneRatio = "Count / taille de la liste" (sous forme "x/y")
    gene_ratio_d <- length(gene_list)  # taille de la liste d'intérêt
    result_subgroup[GO_i,"GeneRatio"] <- paste0(count_g, "/", gene_ratio_d)
    
    # BgRatio = "nb gènes (univers) associés à ce GO / taille de l'univers"
    bgRatio_n <- sum(all_GO_type_filter$go_id == current_GO)
    bgRatio_d <- length(unique(all_GO_type_filter$gene_id))
    result_subgroup[GO_i,"BgRatio"] <- paste0(bgRatio_n, "/", bgRatio_d)
    
    # RichFactor = Count / setSize (= nb gènes totaux du GO dans l'ANNO globale)
    # ATTENTION : tu utilises "all_GO_type" vs "all_GO_type_filter" ? 
    # Souvent on fait Count / bgRatio_n, ou Count / setSize (selon la définition)
    # On peut par exemple faire:
    #   setSize = nombre de gènes annotés à ce GO DANS l'UNIVERS
    setSize <- bgRatio_n 
    result_subgroup[GO_i,"RichFactor"] <- count_g / setSize
    
    # FoldEnrichment = (Count / length(gene_list)) / (bgRatio_n / bgRatio_d)
    FE <- (count_g / gene_ratio_d) / (bgRatio_n / bgRatio_d)
    result_subgroup[GO_i,"FoldEnrichment"] <- FE
    
    # p-value hypergéométrique
    # k = taille de la liste d'intérêt
    # q = count_g -1
    # m = bgRatio_n
    # n = bgRatio_d - bgRatio_n
    if (type_calcul == "hypergeometric") {
      # Test hypergéométrique
      # phyper(q, m, n, k, lower.tail=FALSE) => P(X >= q)
      pval <- phyper(
        q = count_g - 1,
        m = bgRatio_n,
        n = bgRatio_d - bgRatio_n,
        k = gene_ratio_d,
        lower.tail = FALSE
      )
      
    } else if (type_calcul == "fisher") {
      # Test de Fisher
      # On construit la matrice 2x2 :
      #         Dans GO    Pas dans GO
      # Liste     a            b
      # Bg-List   c            d
      #
      # a = count_g
      # b = gene_ratio_d - count_g
      # c = bgRatio_n - count_g
      # d = (bgRatio_d - gene_ratio_d) - (bgRatio_n - count_g)
      
      a <- count_g
      b <- gene_ratio_d - count_g
      c <- bgRatio_n - count_g
      d <- (bgRatio_d - gene_ratio_d) - c
      mat_2x2 <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
      
      # alternative="greater" pour tester la sur-représentation
      fisher_res <- fisher.test(mat_2x2, alternative = "greater")
      pval <- fisher_res$p.value
      
    } else {
      # Par défaut, si type_calcul n'est ni "hypergeometric" ni "fisher"
      pval <- NA
    }
    
    # Enregistrer la p-value
    result_subgroup[GO_i,"pvalue"] <- pval
    
    # Optionnel : on pourrait stocker la liste des gènes pour ce GO
    # Par exemple :
    these_genes <- geneset %>%
      filter(go_id == current_GO) %>%
      pull(gene_id)
    result_subgroup[GO_i,"geneID"] <- paste(these_genes, collapse = "/")
  }
  
  # ---------------------------------------------------------------------------- #
  # 7) p.adjust + qvalue
  # ---------------------------------------------------------------------------- #
  # Correction multiple BH (FDR)
  result_subgroup$p.adjust <- p.adjust(result_subgroup$pvalue, method = pAdjustMethod_i)
  
  # qvalue (nécessite {qvalue})
  # On vérifie qu'il y a assez de p-values non-nulles
  if (sum(!is.na(result_subgroup$pvalue)) > 1) {
    result_subgroup$qvalue <- tryCatch({
      qvalue(result_subgroup$pvalue)$qvalues
    }, error = function(e) {
      message("Erreur perso lors du calcul des q-values : ", e$message)
      rep(NA, length(result_subgroup$pvalue))
    })
  } else {
    result_subgroup$qvalue <- NA
  }
  # ---------------------------------------------------------------------------- #
  # 8) Option : filtrer par pvalueCutoff
  # ---------------------------------------------------------------------------- #
  result_subgroup <- result_subgroup %>%
    filter(p.adjust <= pvalueCutoff_i)
  
  # Optionnel : si tu veux la description GO
  # library(GO.db)
  # go_terms <- keys(GO.db)
  # go_descriptions <- Term(go_terms)
  # go_table <- data.frame(ID = go_terms, Description_GO = go_descriptions)
  # result_subgroup <- left_join(result_subgroup, go_table, by = "ID")
  
  # ---------------------------------------------------------------------------- #
  # 9) Retour
  # ---------------------------------------------------------------------------- #
  if (nrow(result_subgroup) == 0) {
    return(NULL)  # si aucun terme ne passe le cutoff
  } else {
    return(result_subgroup)
  }
}

source(here::here("src/function/run_parallel.R")) # This function allow to performe multiprocessing on multiple platforme 

GO_on_different_group <- function(
    functional_roles = "BP", #MF or CC
    group_info, # here is a dataframe with in ID (the name of the genes in database) and group (the group (cluster or category)). 
    group, # column of the different group
    ID, # column of the DF with the name of the different genes
    universe_i = NULL, # vector containing all the genes you consider to be the reference universe (DEG gene) , By default, I take all the genes, which is not perfect.  
    minGSSize_i = 10, 
    maxGSSize_i = 500,
    pvalueCutoff_i = 0.05,
    pAdjustMethod_i = "BH",
    top = FALSE,
    all_GO = AnnotationDbi::select(x = org.Psativum1c.eg.db, keys = keys(x=org.Psativum1c.eg.db, keytype = "GID"), columns =c("GO","ONTOLOGY") ,  keytype = "GID"),
    type_calcul = "clusterProfiler" # fisher # hypergeometric
){ #keytypes(org.Psativum1c.eg.db) 
  # verification of all inputs ####
  if( !functional_roles %in% c("BP", "CC", "MF")){
    stop("You need to add BP, CC or MF")
  }
  
  # add column in the data frame corresponding to ID and group
  group_info = group_info %>% mutate(group = !!sym(group), ID = !!sym(ID)) %>% dplyr::select(group, ID)
  
  # Select Biological Processes
  all_GO_type <- all_GO[all_GO$ONTOLOGY == functional_roles,] %>% 
    na.omit() %>% 
    dplyr::select("GO", "GID") %>% 
    dplyr::rename(go_id = GO, 
                  gene_id = GID)
  
  cat_col("- Performing Gene Ontology enrichment analysis for each group \n", "green")
  
  num_groups <- length(unique(group_info$group))  # Total number of iterations
  #pb <- txtProgressBar(min = 0, max = num_groups, style = 3)  # Create progress bar
  
  all_results <- list()
  group_list <- unique(group_info$group)  # Store unique groups
  
  if(type_calcul=="clusterProfiler"){
    
    # Fonction interne pour l'analyse d'enrichissement par groupe 
    assign("group_info", group_info, envir = parent.frame())
    assign("all_GO_type", all_GO_type, envir = parent.frame())
    assign("pvalueCutoff_i", pvalueCutoff_i, envir = parent.frame())
    assign("pAdjustMethod_i", pAdjustMethod_i, envir = parent.frame())
    assign("minGSSize_i", minGSSize_i, envir = parent.frame())
    assign("maxGSSize_i", maxGSSize_i, envir = parent.frame())
    assign("universe_i", universe_i, envir = parent.frame())
    
    all_results <- run_parallel(
      group_list = group_list,
      varlist =  c("group_info", "all_GO_type", "pvalueCutoff_i", "pAdjustMethod_i", "minGSSize_i", "maxGSSize_i", "universe_i"),
      function_to_run = run_enrichment
    )
    
    # Supprimer les résultats NULL
    all_results <- all_results[!sapply(all_results, is.null)]
    
    # Combiner tous les résultats dans un seul data frame
    combined_results <- bind_rows(all_results)
    
    # Récupérer la description des termes GO
    go_terms <- keys(GO.db)
    go_descriptions <- Term(go_terms)
    go_table <- data.frame(GO_Term = go_terms, Description_GO = go_descriptions)
    
    names(go_table)[1] <- "ID"
    
    # Fusionner les résultats avec les descriptions GO
    combined_results <- left_join(combined_results, go_table, by = "ID")
    
  #### Here a new section based on a hypergeometric or fisher test ##########    
    # see  (not add in https://ifb-elixirfr.github.io/EBAII/2021/ebaiin2/RNASeq/EBAIIn2_RNASeq.html#3_Du_g%C3%A8ne_au_gene-set)
  } else if (type_calcul %in% c("fisher", "hypergeometric")){
    all_GO_type_filter <- all_GO_type %>% filter(gene_id %in% unique(universe_i))
    go_gene_counts <- all_GO_type_filter %>% 
      dplyr::group_by(go_id) %>% 
      dplyr::summarise(gene_count = n_distinct(gene_id), .groups = "drop")
    
    go_id_to_keep <- go_gene_counts %>%
      filter(gene_count >= minGSSize_i, gene_count <= maxGSSize_i) %>%
      pull(go_id)
    
    all_GO_type_filter <- all_GO_type_filter %>%
      filter(go_id %in% go_id_to_keep)
    
    # test multiprocessing
    assign("group_info", group_info, envir = parent.frame())
    assign("top", top, envir = parent.frame())
    assign("all_GO_type", all_GO_type, envir = parent.frame())
    assign("pvalueCutoff_i", pvalueCutoff_i, envir = parent.frame())
    assign("pAdjustMethod_i", pAdjustMethod_i, envir = parent.frame())
    assign("minGSSize_i", minGSSize_i, envir = parent.frame())
    assign("maxGSSize_i", maxGSSize_i, envir = parent.frame())
    assign("universe_i", universe_i, envir = parent.frame())
    assign("type_calcul", type_calcul, envir = parent.frame())
    assign("all_GO_type_filter", all_GO_type_filter, envir = parent.frame())
    assign("go_gene_counts", go_gene_counts, envir = parent.frame())
    
    
    all_results <- run_parallel(
      group_list = group_list,
      varlist =  c("group_info", "all_GO_type", "pvalueCutoff_i", "pAdjustMethod_i", "minGSSize_i", "maxGSSize_i", "universe_i","all_GO_type_filter", "type_calcul", "top","go_gene_counts"),
      function_to_run = run_enrichment_manual
    )
    
    combined_results <- bind_rows(all_results)
    
    # Get the GO description
    go_terms <- keys(GO.db)
    go_descriptions <- Term(go_terms)
    go_table <- data.frame(GO_Term = go_terms, Description_GO = go_descriptions)
    
    names(go_table)[1] <- "ID"
    combined_results <- left_join(combined_results, go_table, by = "ID")
  }
  
  if(top == FALSE){
    cat_col("- Give all go terms \n", "green")
    return(combined_results)
  }else if(top >0){
    combined_results.top <- combined_results %>%
      dplyr::group_by(group) %>%
      top_n(-top, wt = p.adjust) %>%
      arrange(group, p.adjust)
    
    # Define the list of selected significant GO terms
    Signif_GO <- unique(combined_results.top$Description)
    
    # Filter results based on this list
    group_selected_GO <- filter(combined_results, ID %in% Signif_GO)
    
    # Ensure GeneRatio and BgRatio are character columns
    group_selected_GO_2 <- group_selected_GO
    #str(group_selected_GO_2)
    group_selected_GO_2$GeneRatio <- as.character(group_selected_GO_2$GeneRatio)
    group_selected_GO_2$BgRatio <- as.character(group_selected_GO_2$BgRatio)
    
    # Calculate the fold enrichment
    group_selected_GO_2$FoldEnrichment <- with(group_selected_GO_2, 
                                                 (sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))) /
                                                   (sapply(strsplit(BgRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))))
    
    # Calculate the -log10(pvalue)
    group_selected_GO_2$'|-log10(Pval)|' <- -(log10(group_selected_GO_2$p.adjust))
    
    # Ensure that `Term` column has factor levels including all possible terms
    group_selected_GO_2$Description_GO <- factor(group_selected_GO_2$Description_GO, levels = unique(group_selected_GO_2$Description_GO))
    
    # Now, create a data frame with all combinations of `Term` and `filename`
    all_combinations <- expand.grid(Description_GO = levels(group_selected_GO_2$Description_GO), group = unique(group_info$group))
    
    # Merge this with your original data to fill missing combinations
    group_selected_GO_filled <- merge(all_combinations, group_selected_GO_2, by = c("Description_GO", "group"), all.x = TRUE) 
    
    cat_col("- End of the analyse \n", "green")
    
    return(group_selected_GO_filled)
  }
}