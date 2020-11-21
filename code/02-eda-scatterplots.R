# Paper Fig 1, scatterplots of C, N, CN ratio to EM
library(tidyverse)
library(cowplot)


min1.dat <- read_rds("models/Soil model data.rds") %>%
  # filter(lyr %in% c('ff_total', 'min_1', 'min_2')) %>%
  filter(lyr == "min_1") %>%
  mutate(cnr = sc / sn) %>%
  filter(cnr > quantile(cnr, .01), cnr < quantile(cnr, .99)) %>%
  select(lyr, em = myc_e_pct, c = sc, n = sn, cnr)

pts.alp <- .5
smooth.method <- "auto"

# # facet panels
# min1.dat %>%
#   gather(c, n, cnr, key = "var", value = "soil") %>%
#   mutate(var = factor(
#     var, c("c", "n", "cnr"),
#     c("C (g cm-3)", "N (g cm-3)", "C:N ratio")
#   )) %>%
#   ggplot(aes(em, soil)) +
#   geom_point(alpha = pts.alp) +
#   geom_smooth() +
#   facet_grid(var ~ ., scales = "free", switch = "y") +
#   scale_y_log10(labels = scales::comma) +
#   annotation_logticks(sides = "l") +
#   labs(y = "", x = "EM dominance") + # expression('Soil content (g'~cm^-3*')')
#   theme_bw()
# ggsave("figures/C, N, CN ratio - EM.pdf", h = 10, w = 10 * .618)

# separate panels
c.gg <- ggplot(min1.dat, aes(em, c)) +
  geom_point(alpha = pts.alp) +
  geom_smooth(method = smooth.method) +
  scale_y_log10(labels = scales::comma) +
  annotation_logticks(sides = "l") +
  labs(y = expression("Soil C (g" ~ cm^-3 * ")"), x = "EM dominance")

n.gg <- ggplot(min1.dat, aes(em, n)) +
  geom_point(alpha = pts.alp) +
  geom_smooth(method = smooth.method) +
  scale_y_log10(labels = scales::comma) +
  annotation_logticks(sides = "l") +
  labs(y = expression("Soil N (g" ~ cm^-3 * ")"), x = "EM dominance")

cnr.gg <- ggplot(min1.dat, aes(em, cnr)) +
  geom_point(alpha = pts.alp) +
  geom_smooth(method = smooth.method) +
  scale_y_log10(labels = scales::comma, breaks = c(10, 20, 50, 100)) +
  annotation_logticks(sides = "l") +
  labs(y = "Soil C:N", x = "EM dominance")

plot_grid(cnr.gg,
  plot_grid(c.gg, n.gg, labels = c("b", "c"), nrow = 2),
  labels = c("a", ""), nrow = 1
)

ggsave("figures/Fig 2 Scatterplots.pdf", h = 6.18, w = 10)
