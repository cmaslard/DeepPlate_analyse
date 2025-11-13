# This function saves a given plot (plot_x) as both a PDF and a high-resolution PNG file at specified dimensions.
# fig_export <- function(path, plot_x, height_i, width_i, res_i = 300){
#   #path = "report/multi_omics/plot/LMM/barplot_cluster"
#   # PNG ####
#   png(here::here(paste0(path,".png")), height = height_i*res_i, width = width_i*res_i, res=res_i)
#   print(plot_x)
#   dev.off()
#   
#   # PDF ####
#   pdf(here::here(paste0(path,".pdf")), height = height_i, width = width_i)
#   print(plot_x)
#   dev.off()
#   
#   # SVG ####
#   ggsave(here::here(paste0(path,".svg")),
#          plot =  plot_x,
#          height = height_i,
#          width = width_i)
# }

fig_export <- function(path, plot_x, height_i, width_i, res_i = 300, format = c("png", "pdf", "svg")){
  
  if("png" %in% format){
    png(here::here(paste0(path, ".png")), height = height_i * res_i, width = width_i * res_i, res = res_i)
    print(plot_x)
    dev.off()
  }
  
  if("pdf" %in% format){
    pdf(here::here(paste0(path, ".pdf")), height = height_i, width = width_i)
    print(plot_x)
    dev.off()
  }
  
  if("svg" %in% format){
    ggsave(here::here(paste0(path, ".svg")),
           plot = plot_x,
           height = height_i,
           width = width_i)
  }
}