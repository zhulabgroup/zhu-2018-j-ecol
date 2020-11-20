# In response to reviews, explore patterns among ecoregions
# DO NOT PUBLISH
library(tidyverse)
library(cowplot)
library(R2jags)

if (!file.exists("Figures/Ecoregion/")) {
  dir.create("Figures/Ecoregion/")
}

# make scatterplots -------------------------------------------------------

min1.dat <- read_rds("Models/Soil model data.rds") %>%
  # filter(lyr %in% c('ff_total', 'min_1', 'min_2')) %>%
  filter(lyr == "min_1") %>%
  mutate(cnr = sc / sn) %>%
  filter(cnr > quantile(cnr, .01), cnr < quantile(cnr, .99)) %>%
  select(lyr, eco, em = myc_e_pct, c = sc, n = sn, cnr)

# choose only ecoregions with 50 or more plots
min1.dat <- count(min1.dat, eco) %>%
  rename(nn = n) %>%
  filter(nn >= 50) %>% # at least 50 plots
  left_join(min1.dat, by = "eco") %>%
  select(-nn)

pts.alp <- .5
smooth.method <- "lm"

c.gg <- ggplot(min1.dat, aes(em, c, group = eco, color = eco)) +
  geom_point(alpha = pts.alp) +
  geom_smooth(method = smooth.method, se = F) +
  scale_y_log10(labels = scales::comma) +
  annotation_logticks(sides = "l") +
  labs(y = expression("Soil C (g" ~ cm^-3 * ")"), x = "EM dominance")

n.gg <- ggplot(min1.dat, aes(em, n, group = eco, color = eco)) +
  geom_point(alpha = pts.alp) +
  geom_smooth(method = smooth.method, se = F) +
  scale_y_log10(labels = scales::comma) +
  annotation_logticks(sides = "l") +
  labs(y = expression("Soil N (g" ~ cm^-3 * ")"), x = "EM dominance")

cnr.gg <- ggplot(min1.dat, aes(em, cnr, group = eco, color = eco)) +
  geom_point(alpha = pts.alp) +
  geom_smooth(method = smooth.method, se = F) +
  scale_y_log10(labels = scales::comma, breaks = c(10, 20, 50, 100)) +
  annotation_logticks(sides = "l") +
  labs(y = "Soil C:N", x = "EM dominance")

plot_grid(cnr.gg,
  plot_grid(c.gg, n.gg, labels = c("b", "c"), nrow = 2),
  labels = c("a", ""), nrow = 1
)

ggsave("Figures/Ecoregion/Scatterplots linreg.pdf", h = 6.18, w = 10)

# summarize slope coefficients (mean and 95% CI) from linear regressions
min1.dat %>%
  group_by(eco) %>%
  do(linreg = lm(log(cnr) ~ em, data = .)) %>%
  mutate(
    mean = coef(linreg)[2],
    lwr = confint(linreg)[2, 1],
    upr = confint(linreg)[2, 2]
  ) %>%
  ggplot(aes(eco, mean, ymin = lwr, ymax = upr)) +
  geom_pointrange() +
  labs(x = "Ecoregion (province)", y = "Slope (mean and 95% CI)")
ggsave("Figures/Ecoregion/Slope coefs.pdf", h = 6.18, w = 10)


# run HB models -----------------------------------------------------------

# setup model
run_jags <- function(in.dat, sel.mod = "cnr", ...) {
  # in.dat:  a list of (C, N, X)
  # sel.mod: selected model, 'cnr' or 'conly'
  # ...: n.chains, n.iter

  # CN ratio model
  cnr.model <- function() {

    # likelihood
    for (i in 1:n.plt) {
      C[i] ~ dnorm(k[i] * N[i], tau)
      log(k[i]) <- X[i, ] %*% beta # X is a matrix with 1st col intercept
    }

    # prior
    for (j in 1:n.x) {
      beta[j] ~ dnorm(0, 1e-3)
    }
    tau ~ dgamma(1e-3, 1e-3)
  }

  # C only model
  conly.model <- function() {

    # likelihood
    for (i in 1:n.plt) {
      C[i] ~ dnorm(k[i], tau) # no N[i] here
      log(k[i]) <- X[i, ] %*% beta # X is a matrix with 1st col intercept
    }

    # prior
    for (j in 1:n.x) {
      beta[j] ~ dnorm(0, 1e-3)
    }
    tau ~ dgamma(1e-3, 1e-3)
  }

  # N only model
  nonly.model <- function() {

    # likelihood
    for (i in 1:n.plt) {
      N[i] ~ dnorm(k[i], tau) # no C[i] here
      log(k[i]) <- X[i, ] %*% beta # X is a matrix with 1st col intercept
    }

    # prior
    for (j in 1:n.x) {
      beta[j] ~ dnorm(0, 1e-3)
    }
    tau ~ dgamma(1e-3, 1e-3)
  }

  pars <- c("beta", "tau", "k")

  in.dat$n.plt <- nrow(in.dat$X)
  in.dat$n.x <- ncol(in.dat$X)

  out.jags <- jags.parallel(
    data = in.dat, inits = NULL,
    parameters.to.save = pars,
    model.file = switch(sel.mod,
      cnr = cnr.model,
      conly = conly.model,
      nonly = nonly.model
    ), ...
  )

  out.jags
}

# prep data
all.dat <- read_rds("Models/Soil model data.rds")

# all covariate names
cov.tag <- tribble( # in order of appearance
  ~abbr, ~name, ~group,
  "int", "Intercept", "Intercept",
  "myc_e_pct", "EM", "Mycorrhiza",
  "phycat0_G_pct", "Gymnosperm", "Phylogeny",
  "leaf_ntg", "Leaf N", "Trait",
  "clay_pct", "Clay", "Soil texture",
  "tmp", "Temperature", "Climate",
  "ppt", "Precipitation", "Climate"
)

# all layer names
lyr.tag <- tribble(
  ~abbr, ~name,
  "ff_total", "Forest floor total",
  "min_1", "Mineral layer 1",
  "min_2", "Mineral layer 2",
  "min_total", "Mineral layer total"
)

# all available ecoregion names (note leading space)
eco.tag <- tribble(
  ~abbr, ~name,
  " 212", "212 Laurentian Mixed Forest Province",
  " 223", "223 Central Interior Broadleaf Forest Province",
  " 232", "232 Outer Coastal Plain Mixed Forest Province"
  # ' 221', '221 Eastern Broadleaf Forest Province'
)

# covariates to normalize x_nor = (x - min)/(max - min)
nor.tag <- c("leaf_ntg", "tmp", "ppt")

# run models
mcmc.mat <- NULL
in.dat <- NULL
for (lyr.sel in lyr.tag$abbr) {
  for (eco.sel in eco.tag$abbr) {

    # lyr.sel <- 'min_total'
    cnr.dat <- all.dat %>%
      filter(lyr == lyr.sel, eco == eco.sel) %>%
      arrange(plt)

    # 0-1 transform all covariates
    x.mat <- as.matrix(cnr.dat[, cov.tag$abbr[-1]])
    for (cov.nm in nor.tag) {
      x.min <- min(x.mat[, cov.nm])
      x.max <- max(x.mat[, cov.nm])
      x.mat[, cov.nm] <- (x.mat[, cov.nm] - x.min) / (x.max - x.min)
    }
    x.mat <- cbind(int = 1, x.mat)

    # run jags
    cnr.jags <- run_jags(
      in.dat = list(
        C = cnr.dat$sc,
        N = cnr.dat$sn,
        X = x.mat
      ),
      sel.mod = "cnr",
      n.chains = 3, n.iter = 2e3
    )

    # summarize output
    mcmc.mat <- cnr.jags %>%
      coda::as.mcmc() %>%
      ggmcmc::ggs() %>%
      rename(iter = Iteration, chain = Chain, term = Parameter) %>%
      cbind(lyr = lyr.sel, eco = eco.sel, .) %>%
      rbind(mcmc.mat, .)

    in.dat <- cnr.dat %>%
      rownames_to_column(var = "tag") %>%
      bind_rows(in.dat, .)
  }
}

mcmc.dat <- as_tibble(mcmc.mat)

# summarize parameters
mcmc.sum <- mcmc.dat %>%
  group_by(lyr, eco, term) %>%
  summarize(
    mean = mean(value), sd = sd(value), se = plotrix::std.error(value),
    lwr = quantile(value, .025), med = median(value), upr = quantile(value, .975)
  ) %>%
  ungroup() %>%
  mutate(lyr = as.character(lyr), term = as.character(term))

list(data = in.dat, mcmc = mcmc.dat, summary = mcmc.sum, tag = list(cov = cov.tag, lyr = lyr.tag)) %>%
  write_rds("Models/Ecoregion CN ratio.rds")

# make coefplots
beta.sum <- mcmc.sum %>%
  filter(stringr::str_detect(term, "^beta")) %>%
  mutate(tag = stringr::str_extract(term, "\\d")) %>%
  full_join(cov.tag %>% # change beta coef name
    rownames_to_column("tag"), ., by = "tag") %>%
  select(coef = name, grp = group, lyr, eco, mean:upr) %>%
  mutate(eco = as.character(eco)) %>%
  full_join(lyr.tag, ., by = c("abbr" = "lyr")) %>%
  rename(lyr = name) %>%
  full_join(eco.tag, ., by = c("abbr" = "eco")) %>%
  rename(eco = name) %>%
  select(coef, lyr, eco, mean:upr) %>%
  mutate(
    lyr = factor(lyr, lyr.tag$name, c("Forest floor", lyr.tag$name[-1])),
    coef = factor(coef, cov.tag$name)
  ) %>%
  filter(coef != "Intercept") %>%
  filter(!grepl("^Mineral layer total", lyr)) %>%
  droplevels() %>%
  mutate(lyr = recode(lyr,
    `Mineral layer 1` = "0 - 10 cm\nmineral layer",
    `Mineral layer 2` = "10 - 20 cm\nmineral layer"
  ))

ggplot(beta.sum, aes(
  x = coef,
  y = med, ymin = lwr, ymax = upr,
  col = lyr, shape = eco
)) +
  geom_hline(yintercept = 0, alpha = .5, linetype = "dashed") +
  geom_pointrange(position = position_dodge(width = -.5)) + # minus to reverse factor levels, will lead to warning
  scale_x_discrete(limits = rev(cov.tag$name[-1])) +
  labs(x = "", y = "Coefficient value (median and 95% CI)", col = "", title = "C:N model") +
  coord_flip()
ggsave("Figures/Ecoregion/CN ratio coefs.pdf", w = 10, h = 6.18)
