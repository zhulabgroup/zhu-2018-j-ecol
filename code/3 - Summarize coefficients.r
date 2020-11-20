# Summarize regression coefficients (beta's)
library(tidyverse)

in.path <- "Models/N only.rds"
out.path <- "Figures/N only/"
if (!file.exists(out.path)) {
  dir.create(out.path)
}

# Summarize beta's --------------------------------------------------------

cov.tag <- read_rds(in.path)$tag$cov
lyr.tag <- read_rds(in.path)$tag$lyr
n.plt <- read_rds(in.path)$data %>% # number of plots (sample size)
  group_by(lyr) %>%
  count() %>%
  `[[`("n")

beta.sum <- read_rds(in.path)$summary %>%
  filter(stringr::str_detect(term, "^beta")) %>%
  mutate(tag = stringr::str_extract(term, "\\d")) %>%
  full_join(cov.tag %>% # change beta coef name
    rownames_to_column("tag"), ., by = "tag") %>%
  select(coef = name, grp = group, lyr, mean:upr) %>%
  full_join(lyr.tag, ., by = c("abbr" = "lyr")) %>%
  select(-abbr) %>%
  rename(lyr = name) %>%
  mutate(
    lyr = factor(lyr, lyr.tag$name, paste0(lyr.tag$name, ", n = ", n.plt)),
    coef = factor(coef, cov.tag$name),
    grp = factor(grp, unique(cov.tag$group))
  )

# multi-panel = layers
beta.sum %>%
  filter(coef != "Intercept") %>%
  droplevels() %>%
  ggplot(aes(med, coef, col = grp)) +
  geom_vline(xintercept = 0, alpha = .5, linetype = "dashed") +
  geom_point() +
  geom_errorbarh(aes(xmin = lwr, xmax = upr), height = 0, size = .5) +
  facet_wrap(~lyr, scale = "fixed") +
  scale_y_discrete(limits = rev(cov.tag$name[-1])) +
  labs(x = "Coefficient value (median and 95% CI)", y = "", col = "") +
  theme_bw()
ggsave(paste0(out.path, "Coefs multi-panel.pdf"), width = 10, height = 10 * .618)

# one-panel, group = layers
beta.sum %>%
  filter(coef != "Intercept") %>%
  filter(!grepl("^Mineral layer total", lyr)) %>%
  droplevels() %>%
  ggplot(aes(
    x = coef,
    y = med, ymin = lwr, ymax = upr,
    col = lyr
  )) +
  geom_hline(yintercept = 0, alpha = .5, linetype = "dashed") +
  geom_pointrange(position = position_dodge(width = -.5)) + # minus to reverse factor levels, will lead to warning
  scale_x_discrete(limits = rev(cov.tag$name[-1])) +
  labs(x = "", y = "Coefficient value (median and 95% CI)", col = "", title = "CN ratio") +
  coord_flip() +
  theme_bw()
ggsave(paste0(out.path, "Coefs single-panel.pdf"), width = 7, height = 7 * .618)
