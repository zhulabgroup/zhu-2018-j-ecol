# Make maps for CN ratio and other variables
library(tidyverse)
library(maptools)

in.path <- "Models/CN ratio.rds"
out.path <- "Figures/CN ratio/"
if (!file.exists(out.path)) {
  dir.create(out.path)
}

# Summarize k's -----------------------------------------------------------

k.tag <- read_rds(in.path)$data %>%
  select(lyr, tag, plt)

k.sum <- read_rds(in.path)$summary %>%
  filter(stringr::str_detect(term, "^k")) %>%
  mutate(
    tag = stringr::str_extract(term, "\\(?[0-9,.]+\\)?"),
    lyr = as.character(lyr)
  ) %>%
  inner_join(k.tag, ., by = c("tag", "lyr")) %>%
  select(lyr, plt, mean:upr)

k.dat <- read_rds(in.path)$data %>%
  select(lyr, plt, sc, sn) %>%
  inner_join(k.sum, by = c("lyr", "plt")) %>%
  full_join(read_rds(in.path)$tag$lyr, ., by = c("abbr" = "lyr")) %>%
  mutate(obs = sc / sn) %>%
  select(lyr = name, plt, sc, sn, obs, mean:upr)

# xy plot
k.dat %>%
  ggplot(aes(obs, mean)) +
  geom_pointrange(aes(ymin = lwr, ymax = upr), alpha = .2) +
  geom_abline(intercept = 0, slope = 1, col = "red") +
  facet_wrap(~lyr) +
  scale_x_log10() +
  scale_y_log10() +
  coord_fixed(ratio = 1)
ggsave(paste0(out.path, "States obs and mod k.pdf"), w = 10, h = 10 * .618)


# Plot maps ---------------------------------------------------------------

# prepare base map
usfs.prj <- "+proj=aea +lat_1=38 +lat_2=42 +lat_0=40 +lon_0=-82 +x_0=0 +y_0=0 +ellps=clrk66 +datum=NAD83 +units=m"
latlon.prj <- "+proj=longlat +datum=NAD83"

us_st100.dat <- readShapePoly("Data/GIS/us_st100/us_st100.shp",
  proj4string = CRS(usfs.prj)
) %>%
  spTransform(CRS(latlon.prj)) %>%
  fortify()

k.map <- k.dat %>%
  inner_join(read_rds("Data/Soil/Plot table.rds") %>%
    select(plt, lon, lat), ., by = "plt")

k.map %>%
  select(lon, lat, lyr, obs, mod = med) %>%
  gather("stat", "val", obs:mod) %>%
  mutate(
    log.val = log(val),
    stat = factor(stat, c("obs", "mod"), c("Observed", "Modeled"))
  ) %>%
  group_by(lon, lat, lyr, stat) %>%
  summarize(log.val = mean(log.val)) %>%
  ggplot() +
  geom_point(aes(lon, lat, col = log.val), alpha = .5) +
  scale_color_gradientn(colors = rainbow(10)) + # fields::tim.colors()) +
  facet_grid(lyr ~ stat) +
  geom_path(data = us_st100.dat, aes(long, lat, group = group)) +
  coord_map()

# bad maps -- maybe this isn't a predictive model? or CN ratio can't be captured?
