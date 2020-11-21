# install packages
list.of.packages <- c(
  "tidyverse", "maptools", "cowplot",
  "R2jags" # also install jags http://mcmc-jags.sourceforge.net
)
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)

# set paths
if (!file.exists("figures/")) {
  dir.create("figures/")
}

if (!file.exists("models/")) {
  dir.create("models")
  # file.symlink("/data/ZHULAB/FIA/Mycorrhiza/Models/", "models")
}
