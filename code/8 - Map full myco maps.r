# map a full EM% map showing all plots, regardless of soil data, for paper Fig S1
library(tidyverse)
library(maptools)


usfs.prj <- "+proj=aea +lat_1=38 +lat_2=42 +lat_0=40 +lon_0=-82 +x_0=0 +y_0=0 +ellps=clrk66 +datum=NAD83 +units=m"
latlon.prj <- "+proj=longlat +datum=NAD83"

us_st100.dat <- readShapePoly("Data/GIS/us_st100/us_st100.shp",
  proj4string = CRS(usfs.prj)
) %>%
  spTransform(CRS(latlon.prj)) %>%
  fortify()


# get data ----------------------------------------------------------------

plot.dat <- read_rds("Data/Soil/Plot full table.rds") %>%
  filter(lon >= -100) # eastern US

spec.dat <- read_csv("Data/Species Traits/Species traits for Kai w LeafN final.csv") %>%
  select(spcd, myc, clade = Clade, phycat0 = phy, phycat1 = Phycat1, phycat2 = Phycat2, leaf_ntg = `Leaf N_mg_g`) %>%
  left_join(read_rds("Data/Soil/Species table.rds"), by = "spcd") %>%
  mutate(sp = tolower(species_symbol)) %>%
  select(spcd, sp, genus, species, common_name, myc:leaf_ntg) %>%
  filter(myc != "er") %>% # get rid of 'er' species (only 1)
  mutate(myc = ifelse(myc == "b", "e", myc)) %>% # change 'both' to 'em'
  filter(!is.na(sp))

tree.dat <- read_rds("Data/Soil/Tree full table.rds") %>%
  filter(statuscd == 1) %>% # only for live trees
  mutate(
    plt.area = 4 * ifelse(dia < 5, 13.5e-4, 168e-4), # plot area for 4 units
    ba = 0.25 * pi * ((dia * 0.0254)^2) / plt.area
  ) %>% # unit = m2/ha
  inner_join(plot.dat %>%
    select(plt, lon, lat), ., by = "plt")

myco.dat <- tree.dat %>%
  left_join(spec.dat, by = "spcd") %>%
  group_by(plt, myc) %>%
  summarize(ba.sum = sum(ba, na.rm = T)) %>%
  spread(myc, ba.sum, fill = 0, sep = "_") %>%
  mutate(
    myc_tot = myc_a + myc_e,
    myc_a_pct = myc_a / myc_tot,
    myc_e_pct = myc_e / myc_tot
  ) %>%
  select(plt, myc_e_pct) %>%
  inner_join(plot.dat %>%
    select(plt, lon, lat), ., by = "plt") %>%
  drop_na()


# plot map ----------------------------------------------------------------

myco.dat %>%
  ggplot() +
  geom_point(aes(lon, lat, col = myc_e_pct), shape = ".") + # make points very small to avoid overplotting
  geom_path(data = us_st100.dat, aes(long, lat, group = group)) +
  coord_map() +
  scale_color_gradientn(colors = topo.colors(6)) +
  labs(x = "Longitude (deg)", y = "Latitude (deg)", col = "EM dominance")
ggsave("Figures/Fig S1 EM map.pdf", h = 6.18 * .75, w = 10 * .75)
