# plot coefficients from all layers and all models
library(tidyverse)
library(cowplot)


# read values -------------------------------------------------------------

tag.ls <- c("C only", "N only", "CN ratio")
beta.sum <- vector("list", 3)
names(beta.sum) <- tag.ls

for (tag in tag.ls) {
  in.file <- paste0("Models/", tag, ".rds")

  cov.tag <- read_rds(in.file)$tag$cov
  lyr.tag <- read_rds(in.file)$tag$lyr

  beta.sum[[tag]] <- read_rds(in.file)$summary %>%
    filter(stringr::str_detect(term, "^beta")) %>%
    mutate(tag = stringr::str_extract(term, "\\d")) %>%
    full_join(cov.tag %>% # change beta coef name
      rownames_to_column("tag"), ., by = "tag") %>%
    select(coef = name, grp = group, lyr, mean:upr) %>%
    full_join(lyr.tag, ., by = c("abbr" = "lyr")) %>%
    select(-abbr) %>%
    rename(lyr = name) %>%
    mutate(
      lyr = factor(lyr, lyr.tag$name, c("Forest floor", lyr.tag$name[-1])),
      coef = factor(coef, cov.tag$name),
      grp = factor(grp, unique(cov.tag$group))
    ) %>%
    filter(coef != "Intercept") %>%
    filter(!grepl("^Mineral layer total", lyr)) %>%
    droplevels() %>%
    mutate(
      lyr = recode(lyr,
        `Mineral layer 1` = "0 - 10 cm\nmineral layer",
        `Mineral layer 2` = "10 - 20 cm\nmineral layer"
      ),
      # `Mineral layer total` = '0 - 20 cm\nmineral layer'),
      coef = recode(coef,
        `Leaf nitrogen` = "Leaf N"
      )
    )
}


# make figures ------------------------------------------------------------

coef.lim <- c(-2.7, 2.1)
cov.tag$name <- recode(cov.tag$name, `Leaf nitrogen` = "Leaf N")

c.gg <- beta.sum[["C only"]] %>%
  ggplot(aes(
    x = coef,
    y = med, ymin = lwr, ymax = upr,
    col = lyr
  )) +
  geom_hline(yintercept = 0, alpha = .5, linetype = "dashed") +
  geom_pointrange(position = position_dodge(width = -.5)) + # minus to reverse factor levels, will lead to warning
  scale_x_discrete(limits = rev(cov.tag$name[-1])) +
  ylim(coef.lim) +
  labs(x = "", y = "Coefficient value (median and 95% CI)", col = "", title = "C-only model") +
  coord_flip() +
  theme(legend.position = "none")

n.gg <- beta.sum[["N only"]] %>%
  ggplot(aes(
    x = coef,
    y = med, ymin = lwr, ymax = upr,
    col = lyr
  )) +
  geom_hline(yintercept = 0, alpha = .5, linetype = "dashed") +
  geom_pointrange(position = position_dodge(width = -.5)) + # minus to reverse factor levels, will lead to warning
  scale_x_discrete(limits = rev(cov.tag$name[-1])) +
  ylim(coef.lim) +
  labs(x = "", y = "Coefficient value (median and 95% CI)", col = "", title = "N-only model") +
  coord_flip() +
  theme(legend.position = "none")

cnr.gg <- beta.sum[["CN ratio"]] %>%
  ggplot(aes(
    x = coef,
    y = med, ymin = lwr, ymax = upr,
    col = lyr
  )) +
  geom_hline(yintercept = 0, alpha = .5, linetype = "dashed") +
  geom_pointrange(position = position_dodge(width = -.5)) + # minus to reverse factor levels, will lead to warning
  scale_x_discrete(limits = rev(cov.tag$name[-1])) +
  ylim(coef.lim) +
  labs(x = "", y = "Coefficient value (median and 95% CI)", col = "", title = "C:N model") +
  coord_flip() +
  theme(
    legend.position = c(0.2, 0.9),
    legend.key.height = grid::unit(2, "line")
  )

plot_grid(cnr.gg,
  plot_grid(c.gg, n.gg, labels = c("b", "c"), nrow = 2),
  labels = c("a", ""), nrow = 1
)

ggsave("Figures/Fig 3 Coefplots.pdf", h = 6.18, w = 10)


# make table --------------------------------------------------------------

reshape2::melt(beta.sum) %>%
  select(mod = L1, lyr, coef, stat = variable, val = value) %>%
  filter(stat %in% c("med", "lwr", "upr")) %>%
  mutate(mod = factor(mod, tag.ls), val = round(val, 3)) %>%
  spread(stat, val) %>%
  mutate(val = paste0(med, " (", lwr, ", ", upr, ")")) %>%
  select(Model = mod, Layer = lyr, coef, val) %>%
  spread(coef, val) %>%
  write_csv("Figures/Tab S1 Coef table.csv")
