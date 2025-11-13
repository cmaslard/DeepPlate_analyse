# The create_presence_matrix function generates a presence-absence matrix from a list of omics data, indicating the unique variables found across the datasets and including a specified sign in the output.
create_presence_matrix <- function(omic_data_list, sign) {
  unique_variable <- c(unique(unlist(omic_data_list)),"") # add "" if only on variable
  result_df <- data.frame(Variable = unique_variable)
  result_df <- cbind(result_df, sapply(omic_data_list, function(g) as.integer(unique_variable %in% g)))
  colnames(result_df)[-1] <- names(omic_data_list)
  result_df$Sign <- sign
  result_df <- result_df[1:(nrow(result_df) - 1), ] # Delete the last line of the dataframe 
  return(result_df)
}
library(colorspace)
# This function upsetplot_condition_merge_sign creates an UpSet plot to visualize the intersections of deregulated conditions based on a dataset. It customizes the plot by applying color palettes, setting annotations for "Up" and "Down" regulations, and defining specific intersections to display.
upsetplot_condition_merge_sign<-function(df_deregulate_condition,title_i,color_palette = c('Down' = "#710A22", 'Up' = "#F0567A")){
  v_condition = colnames(df_deregulate_condition)[2:length(df_deregulate_condition)]
  upsetplot_condition <-ComplexUpset::upset(
    df_deregulate_condition, 
    rev(v_condition),
    name='Condition',
    mode = "exclusive_intersection",# 'inclusive_intersection' #'exclusive_union' #'inclusive_union'
    #min_size = 100,
    width_ratio = 0.15,
    height_ratio = 1.5, #0.65
    #max_degree      = 2,
    #n_intersections = 40,
    
    stripes=ComplexUpset::upset_stripes(
      # mapping=aes(color=climat_condition),geom=geom_segment(size=6),
      # colors=c(
      #   'WS_OT'=as.character(climate_pallet[2]),
      #   'WW_SD'=as.character(climate_pallet[3]),
      #   'WS_SD'=as.character(climate_pallet[4])
      geom=geom_segment(size=10),  # make the stripes larger
      #colors=c('grey90', 'white')
      colors = c(
        lighten(as.character(pallet_edaphic_condition[4]), amount = .8),
        lighten(as.character(pallet_edaphic_condition[3]), amount = .8),
        lighten(as.character(pallet_edaphic_condition[2]), amount = .8)
        )
      #),
      #data=condition_metadata
      # geom=geom_segment(size=10)
      #colors=c(as.character(climate_pallet[2:4]))
    ),
    base_annotations = setNames(list(
      list(
        aes = aes(x = intersection, fill = Sign),
        geom = list(
          geom_bar(stat = 'count', position = position_dodge2(width = .9, reverse = TRUE), na.rm = TRUE),
          geom_text(
            vjust = -0.2,
            aes(label = after_stat(paste0(round(..count.., 1)))),
            stat = 'count',
            position = position_dodge2(width = .9, reverse = TRUE),
            na.rm = TRUE
          ),
          scale_fill_manual(values = color_palette),  # Define colors for each 'Sign'
          scale_y_continuous(
            limits = c(0, NA),  # Starts at 0 and adjusts the upper limit automatically
            expand = expansion(mult = c(0, 0.15))  # Adds 10% space above y-axis max.
          ),
          theme_classic(),
          theme(panel.grid = element_blank(),
                axis.title.x = element_blank(),
                axis.text.x = element_blank()
                )
        )
      )
    ), title_i),
    
    set_sizes = (
      ComplexUpset::upset_set_size(
        geom = geom_bar(
          aes(fill = Sign, x = group),
          width = 0.8,
          just = 0.5
        ),
        position = 'right'
      ) + 
        scale_fill_manual(values = color_palette, guide=FALSE)  # Use the same color palette here
    ),
    guides='over',
    querie=list(
      ComplexUpset::upset_query(set='KAY_WS', fill=pallet_edaphic_condition[2],only_components=c('intersections_matrix')),
      ComplexUpset::upset_query(set='KAY_SD', fill=pallet_edaphic_condition[3],only_components=c('intersections_matrix')),
      ComplexUpset::upset_query(set='KAY_WS_SD', fill=pallet_edaphic_condition[4],only_components=c('intersections_matrix')),
      
      ComplexUpset::upset_query(set='WT1_WS', fill=pallet_edaphic_condition[2],only_components=c('intersections_matrix')),
      ComplexUpset::upset_query(set='WT1_SD', fill=pallet_edaphic_condition[3],only_components=c('intersections_matrix')),
      ComplexUpset::upset_query(set='WT1_WS_SD', fill=pallet_edaphic_condition[4],only_components=c('intersections_matrix')),
      ComplexUpset::upset_query(set='W78*_WS', fill=pallet_edaphic_condition[2],only_components=c('intersections_matrix')),
      ComplexUpset::upset_query(set='W78*_SD', fill=pallet_edaphic_condition[3],only_components=c('intersections_matrix')),
      ComplexUpset::upset_query(set='W78*_WS_SD', fill=pallet_edaphic_condition[4],only_components=c('intersections_matrix')),

      ComplexUpset::upset_query(set='WT2_WS', fill=pallet_edaphic_condition[2],only_components=c('intersections_matrix')),
      ComplexUpset::upset_query(set='WT2_SD', fill=pallet_edaphic_condition[3],only_components=c('intersections_matrix')),
      ComplexUpset::upset_query(set='WT2_WS_SD', fill=pallet_edaphic_condition[4],only_components=c('intersections_matrix')),

      ComplexUpset::upset_query(set='E568K_WS', fill=pallet_edaphic_condition[2],only_components=c('intersections_matrix')),
      ComplexUpset::upset_query(set='E568K_SD', fill=pallet_edaphic_condition[3],only_components=c('intersections_matrix')),
      ComplexUpset::upset_query(set='E568K_WS_SD', fill=pallet_edaphic_condition[4],only_components=c('intersections_matrix'))
    ),
    matrix=ComplexUpset::intersection_matrix(
      #mapping=aes(color=climat_condition),
      geom=geom_point(
        shape='circle filled',
        size=3.5,
        stroke=0.45
      )
    ),
    #sort_intersections_by='degree'
    sort_intersections_by=c('degree', 'cardinality'),
    # sort_intersections=FALSE,
    #sort_intersections="descending",
    sort_sets=FALSE,
    # intersections=list(
    #   c('KAY_WS','KAY_SD','KAY_WS_SD','WT1_WS','WT1_SD','WT1_WS_SD'),
    #   c('KAY_WS','KAY_WS_SD','WT1_WS','WT1_WS_SD'),
    #   c('KAY_SD','KAY_WS_SD','WT1_SD','WT1_WS_SD'),
    #   c('KAY_WS','KAY_SD','WT1_WS','WT1_SD'),
    #   c('KAY_WS','WT1_WS'),
    #   c('KAY_SD','WT1_SD'),
    #   c('KAY_WS_SD','WT1_WS_SD'),
    #   'KAY_WS',
    #   'WT1_WS',
    #   'KAY_SD',
    #   'WT1_SD',
    #   "KAY_WS_SD",
    #   "WT1_WS_SD"
     #)
    
    # intersections=list(
    #   c('KAY_WS','KAY_SD','KAY_WS_SD','WT1_WS','WT1_SD','WT1_WS_SD', 'W78*_WS','W78*_SD','W78*_WS_SD'),
    #   c('KAY_WS','KAY_WS_SD','WT1_WS','WT1_WS_SD', 'W78*_WS','W78*_WS_SD'),
    #   c('KAY_WS','KAY_WS_SD','WT1_WS','WT1_WS_SD'),
    #   c('W78*_WS','W78*_WS_SD','WT1_WS','WT1_WS_SD'),
    #   c('KAY_SD','KAY_WS_SD','WT1_SD','WT1_WS_SD'),
    #   c('KAY_WS','KAY_SD','WT1_WS','WT1_SD'),
    #   c('KAY_WS','WT1_WS'),
    #   c('KAY_SD','WT1_SD'),
    #   c('KAY_WS_SD','WT1_WS_SD'),
    #   'KAY_WS',
    #   'WT1_WS',
    #   'W78*_WS',
    #   'KAY_SD',
    #   'WT1_SD',
    #   'W78*_SD',
    #   "KAY_WS_SD",
    #   "WT1_WS_SD", 
    #   "WT1_WS_SD"
    # )
    
  )
  return(upsetplot_condition)
}

