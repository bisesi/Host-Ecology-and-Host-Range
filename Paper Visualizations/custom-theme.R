#ATB
#Figure theme
#theme_bisesi()

#custom theme
theme_bisesi <- function(x){
  theme_bw(base_size=12, base_family="serif") %+replace% 
    theme(
      panel.background  = element_blank(),
      plot.background = element_rect(fill="white", colour=NA), 
      legend.background = element_rect(fill="transparent", colour=NA),
      legend.key = element_rect(fill="transparent", colour=NA),
      strip.background = element_rect(color = "white", fill = "white"),
      strip.text = element_text(size = 12),
      strip.placement = "outside",
      strip.switch.pad.grid = unit(0.05, "in"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank()
    )
}