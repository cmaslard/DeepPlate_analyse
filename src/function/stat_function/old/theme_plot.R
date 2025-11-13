manu_palett=c("#1E8449",  "#5E5A93","#D35400","#E6AB02","#2E86C1","#A7790C","#88a70c","#666666","#51870C")
manu_palett2=c("#1E8449",  "#5E5A93","#D35400","#E6AB02","#2E86C1","#A7790C","#88a70c","#666666","#51870C")

hex_pallet=c("#EF5B0C","#003865")

temperature_color=c("#274690","#950714")

temperature_color=c("#457B9D","#E63946")


my_theme=theme_classic()+theme(
  plot.title = element_text(size = 14, face = "bold"),
  axis.title = element_text(size = 12),
  axis.text = element_text(size = 10),
  #legend.title = element_blank(),
  legend.text = element_text(size = 10),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()#,
  # axis.ticks.length = unit(-0.1, "cm"),  # Ajuste la longueur des tics (vers l'intérieur)
  # axis.ticks.margin = unit(0.1, "cm")    # Ajuste la marge des tics par rapport à l'axe attention a la superposition
)
