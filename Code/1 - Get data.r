# Organize data to be read by CNR model
# Inherit code from HB model
library(tidyverse)
library(maptools)


# Get plot data -----------------------------------------------------------

plot.dat <- read_rds("Data/Soil/Plot table.rds") %>%
  mutate(eco = stringr::str_sub(ecosubcd, 1, 4)) %>% # ecoregion: province level (highest)
  filter(lon >= -100) # eastern US


# Get climate data --------------------------------------------------------

coords.mat <- as.matrix(plot.dat[, c("lon", "lat")])
worldclim.ras <- raster::getData("worldclim", var = "bio", res = 10, download = T, path = tempdir())
worldclim.mat <- raster::extract(worldclim.ras, coords.mat)
worldclim.dat <- tibble(
  tmp = worldclim.mat[, 1] / 10, # bio1 is mean annual temp (needs to divided by 10, in deg c)
  ppt = worldclim.mat[, 12]
) # bio12 is total annaul precip (in mm)
clim.dat <- plot.dat %>%
  select(plt, lon, lat, eco) %>%
  bind_cols(worldclim.dat)


# Get soil data -----------------------------------------------------------

soil.dat <- read_rds("Data/Soil/Soil table.rds") %>%
  mutate(
    c.org = bulk_density * c_org_pct / 100,
    c.iog = bulk_density * c_inorg_pct / 100,
    c.tot = bulk_density * c_total_pct / 100, # unit = g/(cm3)
    n.tot = bulk_density * n_total_pct / 100,
    lyr = tolower(layer_type)
  ) %>%
  select(plt, lyr, sc = c.tot, sn = n.tot, txtrlyr1, txtrlyr2) %>%
  gather(txtrlyr, txtrcode, txtrlyr1:txtrlyr2) %>% # average clay from 0-4 and 4-8 inch
  mutate(
    clay = recode(txtrcode,
      `0` = 0, `1` = .45, `2` = .6, `3` = .1, `4` = .1, `9` = 999
    ), # Nina's email
    clay = ifelse(clay == 999, NA, clay)
  ) %>%
  group_by(plt, lyr) %>%
  summarize(
    sc = mean(sc, na.rm = T),
    sn = mean(sn, na.rm = T),
    clay = mean(clay, na.rm = T)
  ) %>%
  filter(lyr %in% c("ff_total", "min_1", "min_2")) # only relevant layers; other layers too much missing data

# calculate min_total data and combine with other layers
soil.dat <- soil.dat %>%
  filter(lyr %in% c("min_1", "min_2")) %>%
  complete(plt = unique(plt), lyr = c("min_1", "min_2")) %>%
  group_by(plt) %>%
  summarize(sc = mean(sc, na.rm = F), sn = mean(sn, na.rm = F), clay = mean(clay, na.rm = F)) %>% # average min1 and min2
  mutate(lyr = "min_total") %>%
  select(plt, lyr, sc, sn, clay) %>%
  bind_rows(soil.dat) %>%
  arrange(plt, lyr) %>%
  rename(clay_pct = clay) %>%
  inner_join(plot.dat %>%
    select(plt, lon, lat, eco), ., by = "plt")


# Get species trait data --------------------------------------------------

spec.dat <- read_csv("Data/Species Traits/Species traits for Kai w LeafN final.csv") %>%
  select(spcd, myc, clade = Clade, phycat0 = phy, phycat1 = Phycat1, phycat2 = Phycat2, leaf_ntg = `Leaf N_mg_g`) %>%
  left_join(read_rds("Data/Soil/Species table.rds"), by = "spcd") %>%
  mutate(sp = tolower(species_symbol)) %>%
  select(spcd, sp, genus, species, common_name, myc:leaf_ntg) %>%
  filter(myc != "er") %>% # get rid of 'er' species (only 1)
  mutate(myc = ifelse(myc == "b", "e", myc)) %>% # change 'both' to 'em'
  filter(!is.na(sp))

# # set "Salix" to AM species
# spec2.dat <- spec.dat %>%
#   mutate(myc = ifelse(genus == 'Salix', 'a', myc))

# Get tree basal area (BA) data -------------------------------------------

tree.dat <- read_rds("Data/Soil/Tree table.rds") %>%
  filter(statuscd == 1) %>% # only for live trees
  mutate(
    plt.area = 4 * ifelse(dia < 5, 13.5e-4, 168e-4), # plot area for 4 units
    ba = 0.25 * pi * ((dia * 0.0254)^2) / plt.area
  ) %>% # unit = m2/ha
  inner_join(plot.dat %>%
    select(plt, lon, lat), ., by = "plt")

ba.myc <- tree.dat %>%
  left_join(spec.dat, by = "spcd") %>%
  group_by(plt, myc) %>%
  summarize(ba.sum = sum(ba, na.rm = T)) %>%
  spread(myc, ba.sum, fill = 0, sep = "_") %>%
  mutate(
    myc_tot = myc_a + myc_e,
    myc_a_pct = myc_a / myc_tot,
    myc_e_pct = myc_e / myc_tot
  ) %>%
  select(plt, myc_e_pct) # set myc_a as baseline

ba.clade <- tree.dat %>%
  left_join(spec.dat, by = "spcd") %>%
  group_by(plt, clade) %>%
  summarize(ba.sum = sum(ba, na.rm = T)) %>%
  spread(clade, ba.sum, fill = 0, sep = "_") %>%
  mutate(
    clade_tot = clade_Asterids + clade_Fabids + clade_Magnolids + clade_Malvids + clade_Monocot + clade_Pinales,
    clade_Asterids_pct = clade_Asterids / clade_tot,
    clade_Fabids_pct = clade_Fabids / clade_tot,
    clade_Magnolids_pct = clade_Magnolids / clade_tot,
    clade_Malvids_pct = clade_Malvids / clade_tot,
    clade_Monocot_pct = clade_Monocot / clade_tot,
    clade_Pinales_pct = clade_Pinales / clade_tot
  ) %>%
  select(plt, clade_Asterids_pct, clade_Fabids_pct, clade_Malvids_pct, clade_Monocot_pct, clade_Pinales_pct) # set Magnolids as baseline

ba.phycat0 <- tree.dat %>%
  left_join(spec.dat, by = "spcd") %>%
  group_by(plt, phycat0) %>%
  summarize(ba.sum = sum(ba, na.rm = T)) %>%
  spread(phycat0, ba.sum, fill = 0, sep = "_") %>%
  mutate(
    phycat0_tot = phycat0_A + phycat0_G,
    phycat0_A_pct = phycat0_A / phycat0_tot,
    phycat0_G_pct = phycat0_G / phycat0_tot
  ) %>%
  select(plt, phycat0_G_pct) # set A as baseline

ba.phycat1 <- tree.dat %>%
  left_join(spec.dat, by = "spcd") %>%
  group_by(plt, phycat1) %>%
  summarize(ba.sum = sum(ba)) %>%
  spread(phycat1, ba.sum, fill = 0, sep = "_") %>%
  mutate(
    phycat1_tot = phycat1_Fagaceae + phycat1_Other + phycat1_Pinaceae + phycat1_Sapindaceae,
    phycat1_Fagaceae_pct = phycat1_Fagaceae / phycat1_tot,
    phycat1_Other_pct = phycat1_Other / phycat1_tot,
    phycat1_Pinaceae_pct = phycat1_Pinaceae / phycat1_tot,
    phycat1_Sapindaceae_pct = phycat1_Sapindaceae / phycat1_tot
  ) %>%
  select(plt, phycat1_Fagaceae_pct, phycat1_Pinaceae_pct, phycat1_Sapindaceae_pct) # set Other as basline

ba.phycat2 <- tree.dat %>%
  left_join(spec.dat, by = "spcd") %>%
  group_by(plt, phycat2) %>%
  summarize(ba.sum = sum(ba)) %>%
  spread(phycat2, ba.sum, fill = 0, sep = "_") %>%
  mutate(
    phycat2_tot = phycat2_Betulaceae + phycat2_Cupressaceae + phycat2_Fagaceae + phycat2_Juglandaceae + phycat2_Oleaceae +
      phycat2_Other + phycat2_Pinaceae + phycat2_Salicaceae + phycat2_Sapindaceae,
    phycat2_Betulaceae_pct = phycat2_Betulaceae / phycat2_tot,
    phycat2_Cupressaceae_pct = phycat2_Cupressaceae / phycat2_tot,
    phycat2_Fagaceae_pct = phycat2_Fagaceae / phycat2_tot,
    phycat2_Juglandaceae_pct = phycat2_Juglandaceae / phycat2_tot,
    phycat2_Oleaceae_pct = phycat2_Oleaceae / phycat2_tot,
    phycat2_Other_pct = phycat2_Other / phycat2_tot,
    phycat2_Pinaceae_pct = phycat2_Pinaceae / phycat2_tot,
    phycat2_Salicaceae_pct = phycat2_Salicaceae / phycat2_tot,
    phycat2_Sapindaceae_pct = phycat2_Sapindaceae / phycat2_tot
  ) %>%
  select(plt, phycat2_Betulaceae_pct:phycat2_Sapindaceae_pct, -phycat2_Other_pct) # set Other as baseline

ba.trait <- tree.dat %>%
  left_join(spec.dat, by = "spcd") %>%
  group_by(plt) %>%
  summarize(leaf_ntg = sum(ba * leaf_ntg, na.rm = T) / sum(ba, na.rm = T))

ba.dat <- ba.myc %>%
  # full_join(ba.clade, by = 'plt') %>%
  full_join(ba.phycat0, by = "plt") %>%
  # full_join(ba.phycat1, by = 'plt') %>%
  # full_join(ba.phycat2, by = 'plt') %>%
  full_join(ba.trait, by = "plt") %>%
  inner_join(plot.dat %>%
    select(plt, lon, lat, eco), ., by = "plt")


# Combine all data --------------------------------------------------------

all.dat <- soil.dat %>%
  filter(sc != 0, sn != 0) %>% # can't infer C:N ratio if either C or N is 0
  inner_join(ba.dat, by = c("plt", "lon", "lat", "eco")) %>%
  inner_join(clim.dat, by = c("plt", "lon", "lat", "eco")) %>%
  drop_na() %>%
  arrange(plt)

write_rds(all.dat, "Models/Soil model data.rds")


# EDA plots ---------------------------------------------------------------

# correlation plots
all.dat %>%
  # filter(lyr == 'ff_total') %>%
  select(-plt, -lon, -lat, -lyr, -sc, -sn) %>%
  GGally::ggcorr(label = T)
ggsave("Figures/Correlation matrix.pdf", w = 7, h = 7)

# maps of available plots
# prepare base map
usfs.prj <- "+proj=aea +lat_1=38 +lat_2=42 +lat_0=40 +lon_0=-82 +x_0=0 +y_0=0 +ellps=clrk66 +datum=NAD83 +units=m"
latlon.prj <- "+proj=longlat +datum=NAD83"

us_st100.dat <- readShapePoly("Data/GIS/us_st100/us_st100.shp",
  proj4string = CRS(usfs.prj)
) %>%
  spTransform(CRS(latlon.prj)) %>%
  fortify()

all.dat %>%
  select(lon, lat) %>%
  unique() %>%
  ggplot() +
  geom_point(aes(lon, lat), alpha = .5) +
  geom_path(data = us_st100.dat, aes(long, lat, group = group)) +
  coord_map()
ggsave("Figures/Plot map.pdf", w = 10, h = 10)

# # marginal plot of C ~ N | myc or others
# var.ls <- all.dat %>%
#   select(-plt, -lon, -lat, -lyr, -sc, -sn) %>%
#   colnames()
#
# pdf("Figures/Marginal plots.pdf", w = 10, h = 10 * .618)
# for (grp.var in var.ls) {
#   sub.dat <- all.dat[, c("lyr", "sc", "sn", grp.var)] %>%
#     `colnames<-`(c("lyr", "sc", "sn", "grp.val"))
#
#   sub.gg <- sub.dat %>%
#     mutate(grp.cat = ifelse(grp.val > mean(grp.val),
#       paste0("High (", mean(grp.val) %>% round(2), ", ", max(grp.val) %>% round(2), "]"),
#       paste0("Low [", min(grp.val) %>% round(2), ", ", mean(grp.val) %>% round(2), "]")
#     )) %>%
#     ggplot(aes(sn, sc, col = grp.cat, group = grp.cat)) +
#     geom_point(alpha = .5) +
#     geom_smooth(method = "lm", formula = y ~ x - 1, se = F, fullrange = T) +
#     facet_wrap(~lyr) +
#     scale_x_sqrt() +
#     scale_y_sqrt() +
#     labs(x = "Soil nitrogen (g/cm3)", y = "Soil carbon (g/cm3)", col = grp.var)
#
#   print(sub.gg)
# }
# dev.off()
#
#
# # Summarize for Rich Phillips, 12/5/2017 ----------------------------------
#
# phillips.dat <- all.dat %>%
#   select(em_pct = myc_e_pct, sc, sn) %>%
#   mutate(em_cut = cut(em_pct, breaks = seq(0, 1, by = .1), include.lowest = T)) %>%
#   group_by(em_cut) %>%
#   summarize(
#     c_mean = mean(sc), c_se = plotrix::std.error(sc),
#     n_mean = mean(sn), n_se = plotrix::std.error(sn)
#   )
#
# ggplot(phillips.dat, aes(em_cut)) +
#   geom_point(aes(y = c_mean)) +
#   geom_point(aes(y = n_mean)) +
#   scale_y_log10()
#
# write_csv(phillips.dat, "Models/CN summary.csv")
#
# unique(all.dat$plt) %>% length()
