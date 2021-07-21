
library(ggplot2)
library(cowplot)

data = read.csv(here::here("Fitting", "Fitted_params", "params_freq_burst.csv"))
data1 = data[1:50000,]
data2 = data[50001:100000,]

#beta l alpha tau

pbeta = ggplot() +
  geom_line(data = data1, aes(seq(1,50000,1), beta)) +
  geom_line(data = data2, aes(seq(1,50000,1), beta), colour = "red") +
  scale_x_continuous(breaks = c(0,25000,50000)) +
  theme_bw() +
  labs(x = "Step", y = "Adsorption rate") +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12),
        plot.title = element_text(face = "bold"))

pbetahist = ggplot() +
  geom_histogram(data = data, aes(beta)) +
  #geom_histogram(aes(rnorm(50000,40,7)), fill = "blue", alpha = 0.3) +
  theme_bw() +
  labs(x = "Adsorption rate", y = "Count") +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12),
        plot.title = element_text(face = "bold"))

pL = ggplot() +
  geom_line(data = data1, aes(seq(1,50000,1), L)) +
  geom_line(data = data2, aes(seq(1,50000,1), L), colour = "red") +
  scale_x_continuous(breaks = c(0,25000,50000)) +
  theme_bw() +
  labs(x = "Step", y = "Burst size") +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12),
        plot.title = element_text(face = "bold"))

pLhist = ggplot() +
  geom_histogram(data = data, aes(L)) +
  geom_histogram(aes(rnorm(50000,40,7)), fill = "blue", alpha = 0.3) +
  theme_bw() +
  labs(x = "Burst size", y = "Count") +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12),
        plot.title = element_text(face = "bold"))


palpha = ggplot() +
  geom_line(data = data1, aes(seq(1,50000,1), alpha)) +
  geom_line(data = data2, aes(seq(1,50000,1), alpha), colour = "red") +
  scale_x_continuous(breaks = c(0,25000,50000)) +
  theme_bw() +
  labs(x = "Step", y = "Prop. transducing phage") +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12),
        plot.title = element_text(face = "bold"))

palphahist = ggplot() +
  geom_histogram(data = data, aes(alpha)) +
  #geom_histogram(aes(rnorm(50000,40,7)), fill = "blue", alpha = 0.3) +
  theme_bw() +
  labs(x = "Prop. transducing phage", y = "Count") +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12),
        plot.title = element_text(face = "bold"))


ptau = ggplot() +
  geom_line(data = data1, aes(seq(1,50000,1), tau)) +
  geom_line(data = data2, aes(seq(1,50000,1), tau), colour = "red") +
  scale_x_continuous(breaks = c(0,25000,50000)) +
  theme_bw() +
  labs(x = "Step", y = "Latent period") +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12),
        plot.title = element_text(face = "bold"))

ptauhist = ggplot() +
  geom_histogram(data = data, aes(tau)) +
  geom_histogram(aes(rnorm(50000,0.67,0.07)), fill = "blue", alpha = 0.3) +
  theme_bw() +
  labs(x = "Latent period", y = "Count") +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        strip.text.x = element_text(size=12),
        plot.title = element_text(face = "bold"))


plot_grid(pbeta, NULL, pbetahist,
          pL, NULL, pLhist,
          palpha, NULL, palphahist,
          ptau, NULL, ptauhist,
          nrow = 4,
          rel_widths = c(0.8,0.05,1),
          labels)

ggsave("supp_fig_6.png", height = 12, width = 10, dpi = 600)
