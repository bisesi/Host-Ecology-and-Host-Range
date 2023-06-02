#ATB
#Generate all png outputs for figures
#Figs 2, 3, 4, 5, 6, Supp 1, 2, 3
#Figure 1 generated using biorender

#load library
library("ggtext")

#load figure 2 
rm(list = ls())
ecoli = "#000000"
senterica = "#BB5566"
eh7 = "#DDAA33"
p22vir = "#004488"
source(here::here("figures", "final-figs", "scripts", "figure-2.R"))

png(here::here("figures", "final-figs", "imgs", "figure-2.png"), res = 300, width = 2750, height = 2350)
fig2
dev.off()

#load figure 3
rm(list = ls())
ecoli = "#000000"
senterica = "#BB5566"
eh7 = "#DDAA33"
p22vir = "#004488"
source(here::here("figures", "final-figs", "scripts", "figure-3.R"))

png(here::here("figures", "final-figs", "imgs", "figure-3.png"), res = 300, width = 3000, height = 2600)
fig3
dev.off()

#load figure 4
rm(list = ls())
ecoli = "#000000"
senterica = "#BB5566"
eh7 = "#DDAA33"
p22vir = "#004488"
source(here::here("figures", "final-figs", "scripts", "figure-4.R"))

png(here::here("figures", "final-figs", "imgs", "alt-figure-4.png"), res = 300, width = 3200, height = 2800)
alt_fig4
dev.off()

#load figure 5
rm(list = ls())
ecoli = "#000000"
senterica = "#BB5566"
eh7 = "#DDAA33"
p22vir = "#004488"
source(here::here("figures", "final-figs", "scripts", "figure-5.R"))

png(here::here("figures", "final-figs", "imgs", "figure-5.png"), res = 300, width = 2750, height = 2350)
fig5
dev.off()

#load figure 6
rm(list = ls())
ecoli = "#000000"
senterica = "#BB5566"
eh7 = "#DDAA33"
p22vir = "#004488"
source(here::here("figures", "final-figs", "scripts", "figure-6.R"))

png(here::here("figures", "final-figs", "imgs", "figure-6.png"), res = 300, width = 1600, height = 1000)
fig6
dev.off()

#load supplemental figure 1
rm(list = ls())
ecoli = "#000000"
senterica = "#BB5566"
eh7 = "#DDAA33"
p22vir = "#004488"
source(here::here("figures", "final-figs", "scripts", "supplemental-figure-1.R"))

png(here::here("paper-supplement", "figs", "supplemental-figure-1.png"), res = 300, width = 3000, height = 2600)
supp_fig1
dev.off()

#load supplemental figure 2
rm(list = ls())
ecoli = "#000000"
senterica = "#BB5566"
eh7 = "#DDAA33"
p22vir = "#004488"
source(here::here("figures", "final-figs", "scripts", "supplemental-figure-2.R"))

png(here::here("paper-supplement", "figs", "supplemental-figure-2.png"), res = 200, width = 2200, height = 1800)
supp_fig2
dev.off()
