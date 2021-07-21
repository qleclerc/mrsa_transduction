
library(ggpubr)
library(cowplot)
library(png)
library(ggplot2)

pa = ggplot() + background_image(readPNG(here::here("Figures", "fig3a.png"))) + theme_void()
pb = ggplot() + background_image(readPNG(here::here("Figures", "fig3b.png"))) + theme_void()
pd = ggplot() + background_image(readPNG(here::here("Figures", "fig3c.png"))) + theme_void()

pc = ggplot() +
  geom_line(aes(x = c(1,1e1,5e1,1e2,5e2,1e3,5e3,1e4,5e4,
                      1e5,5e5,1e6,5e6,1e7,5e7,1e8,5e8,1e9),
                y = 1*(1-c(1,1e1,5e1,1e2,5e2,1e3,5e3,1e4,5e4,
                           1e5,5e5,1e6,5e6,1e7,5e7,1e8,5e8,1e9)/1e9)), size = 0.8) +
  scale_x_continuous(trans=log10_trans(),
                     breaks=trans_breaks("log10", function(x) 10^x, n=6),
                     labels=trans_format("log10", math_format(10^.x))) +
  theme_bw() +
  labs(x = "cfu per mL", y = "Phage parameter multiplier") +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12),
        plot.title = element_text(face = "bold"))

plot_grid(plot_grid(pa,NULL,pc, rel_widths = c(1,0.05,0.7),
                    nrow = 1, labels = c("a", "", "d")),
          NULL,
          plot_grid(pb,NULL,pd, rel_heights = c(1,0.05,1),
                    nrow = 3, labels = c("b", "", "c")),
          nrow = 3,
          rel_heights = c(0.3,0.05,1))

plot_grid(plot_grid(pa,NULL,pc, rel_heights = c(0.7,0.05,0.7),
                    ncol = 1, labels = c("a", "", "c")),
          NULL,
          plot_grid(pb,NULL,pd, rel_heights = c(1,0.05,1),
                    ncol = 1, labels = c("b", "", "d")),
          ncol = 3,
          rel_widths = c(0.5,0.05,1))

ggsave("fig3full.png", height = 15, width = 10, dpi = 600)
