# Run full HB model, can use different responses and covariates
library(tidyverse)
library(R2jags)

in.path <- "Models/Soil model data.rds"
out.path <- "Models/CN ratio.rds"


# Setup model -------------------------------------------------------------

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

  out.jags <- jags(
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


# Prepare data ------------------------------------------------------------

all.dat <- read_rds(in.path)

# all covariate names
cov.tag <- tribble( # in order of appearance
  ~abbr, ~name, ~group,
  "int", "Intercept", "Intercept",
  "myc_e_pct", "EM", "Mycorrhiza",
  "phycat0_G_pct", "Gymnosperm", "Phylogeny",
  "leaf_ntg", "Leaf nitrogen", "Trait",
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

# covariates to normalize x_nor = (x - min)/(max - min)
nor.tag <- c("leaf_ntg", "tmp", "ppt")

# Run model ---------------------------------------------------------------

mcmc.mat <- NULL
in.dat <- NULL
for (lyr.sel in lyr.tag$abbr) {

  # lyr.sel <- 'min_total'
  cnr.dat <- all.dat %>%
    filter(lyr == lyr.sel) %>%
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
    cbind(lyr = lyr.sel, .) %>%
    rbind(mcmc.mat, .)

  in.dat <- cnr.dat %>%
    rownames_to_column(var = "tag") %>%
    bind_rows(in.dat, .)
}

mcmc.dat <- as_tibble(mcmc.mat)

# Summarize parameters ----------------------------------------------------

mcmc.sum <- mcmc.dat %>%
  group_by(lyr, term) %>%
  summarize(
    mean = mean(value), sd = sd(value), se = plotrix::std.error(value),
    lwr = quantile(value, .025), med = median(value), upr = quantile(value, .975)
  ) %>%
  ungroup() %>%
  mutate(lyr = as.character(lyr), term = as.character(term))

list(data = in.dat, mcmc = mcmc.dat, summary = mcmc.sum, tag = list(cov = cov.tag, lyr = lyr.tag)) %>%
  write_rds(out.path)
