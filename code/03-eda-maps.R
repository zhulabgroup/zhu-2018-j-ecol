# make maps to show CN ratio and EM%
library(tidyverse)
library(maptools)
library(cowplot)

usfs.prj <- "+proj=aea +lat_1=38 +lat_2=42 +lat_0=40 +lon_0=-82 +x_0=0 +y_0=0 +ellps=clrk66 +datum=NAD83 +units=m"
latlon.prj <- "+proj=longlat +datum=NAD83"

us_st100.dat <- readShapePoly("data-raw/GIS/us_st100.shp",
  proj4string = CRS(usfs.prj)
) %>%
  spTransform(CRS(latlon.prj)) %>%
  fortify()

min1.dat <- read_rds("models/Soil model data.rds") %>%
  filter(lyr == "min_1") %>%
  mutate(cnr = sc / sn) %>%
  filter(cnr > quantile(cnr, .01), cnr < quantile(cnr, .99))


# make C:N ratio plots ----------------------------------------------------

cnr.brk <- c(
  1e-4, 2e-4, 5e-4,
  1e-3, 2e-3, 5e-3,
  1e-2, 2e-2, 5e-2,
  1e-1, 2e-1, 5e-1,
  1, 2, 5,
  1e1, 2e1, 5e1,
  1e2, 2e2, 5e2
)
cnr.lab <- "Soil C:N"
c.lab <- expression("Soil C (g" ~ cm^-3 * ")")
n.lab <- expression("Soil N (g" ~ cm^-3 * ")")
myc.lab <- "EM dominance"
lon.lab <- "Longitude (deg)"
lat.lab <- "Latitude (deg)"
pts.alp <- .5

c.gg <- min1.dat %>%
  ggplot() +
  geom_point(aes(lon, lat, col = sc), alpha = pts.alp) +
  geom_path(data = us_st100.dat, aes(long, lat, group = group)) +
  coord_map() +
  scale_color_gradientn(
    trans = "log",
    breaks = cnr.brk, labels = cnr.brk,
    colors = topo.colors(6)
  ) +
  labs(x = lon.lab, y = lat.lab, color = c.lab)

n.gg <- min1.dat %>%
  ggplot() +
  geom_point(aes(lon, lat, col = sn), alpha = pts.alp) +
  geom_path(data = us_st100.dat, aes(long, lat, group = group)) +
  coord_map() +
  scale_color_gradientn(
    trans = "log",
    breaks = cnr.brk, labels = cnr.brk,
    colors = topo.colors(6)
  ) +
  labs(x = lon.lab, y = lat.lab, color = n.lab)

cnr.gg <- min1.dat %>%
  ggplot() +
  geom_point(aes(lon, lat, col = cnr), alpha = pts.alp) +
  geom_path(data = us_st100.dat, aes(long, lat, group = group)) +
  coord_map() +
  scale_color_gradientn(
    trans = "log",
    breaks = cnr.brk, labels = cnr.brk,
    colors = topo.colors(6)
  ) +
  labs(x = lon.lab, y = lat.lab, color = cnr.lab)

myc.gg <- min1.dat %>%
  ggplot() +
  geom_point(aes(lon, lat, col = myc_e_pct), alpha = pts.alp) +
  geom_path(data = us_st100.dat, aes(long, lat, group = group)) +
  coord_map() +
  scale_color_gradientn(colors = topo.colors(6)) +
  labs(x = lon.lab, y = lat.lab, col = myc.lab)

plot_grid(plot_grid(myc.gg, cnr.gg, labels = c("a", "b"), ncol = 1),
  plot_grid(c.gg, n.gg, labels = c("c", "d"), ncol = 1),
  labels = c("", "")
)
ggsave("figures/Fig 1 Maps.pdf", h = 6.18, w = 10)
