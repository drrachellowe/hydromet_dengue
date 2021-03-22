# Combined effects of hydrometeorological hazards and urbanisation on dengue risk in Brazil: a spatiotemporal modelling study

# Rachel Lowe (2021)
# https://github.com/drrachellowe/hydromet_dengue

# R script to visualise and explore dengue, climate and socio-economic datasets

# load packages and data
source("00_load_packages_data.R")

# extract geographical and socio-economic variables from first time slice

df <- dplyr::select(data[data$time == 1,], 
             micro_code, region_name, biome_name, main_climate, 
             urban, water_network, water_shortage)

# obtain correlation between percentage of residents living in urban areas and with access to the piped water network
cor.test(df$urban, df$water_network)

# obtain correlation between percentage of residents living in urban areas and frequency of water supply shortages
cor.test(df$urban, df$water_shortage)

# rename code column in the shape object and attach spatial data
map_df <- full_join(map, df, by = c("code" = "micro_code"))

# Manuscript Fig 1
# plot dengue incidence rate (per 100,000 inhabitants) heat maps (month and year) per state 
dir_facet <- 
  data %>% 
  group_by(year, month, state_code) %>%
  # calculate state level incidence rate
  summarise(cases = sum(dengue_cases),
            pop = sum(population)) %>% 
  mutate(var = cases / pop * 10^5) %>% 
  # add the predefined state grid by state code
  left_join(grid, by = c("state_code" = "code_num")) %>% 
  ggplot(aes(x = month, y = year, fill = var)) + 
  geom_raster() +
  ylab("Year") + 
  xlab("Month") +
  scale_fill_gradientn(name = "DIR", colours=brewer.pal(9, "PuRd"), trans = "log1p", breaks = c(0, 10, 100, 300, 1000), labels = c(0, 10, 100, 300, 1000) ) + 
  scale_y_continuous() +
  scale_x_continuous(breaks = c(1,4,7,10), labels = c("Jan", "Apr", "Jul", "Oct")) +
  theme_bw() + 
  # organise by state name in grid file
  facet_geo( ~name, grid = grid)

ggsave("figs/fig_01_dir_facet.eps", height = 30, width = 25, units = "cm")

# Manuscript Fig 2a
# plot map of % residents living in urban area 
urban_map <- ggplot(map_df) + 
  geom_sf(aes(fill = urban), lwd = 0) +
  scale_fill_gradient_tableau("Blue-Teal") +
  labs(fill = "% residents living \n in urban areas") + 
  theme_void()

# ggsave(urban_map, filename = "figs/fig_02a_urban.eps", height = 4, width = 5)

# Manuscript Fig 2b

# water access against level of urbanisation
water_scatter <- ggplot(df) +
  geom_point(aes(x = urban, y = water_network, colour = region_name, 
                 shape = region_name), size = 3) +
  xlab("% residents living in urban areas") +
  ylab("% residents with access to water network") +
  labs(colour = "Region") +
  scale_color_tableau("Nuriel Stone", direction = 1, name = "Region") +
  scale_shape_manual(values = c(15, 17, 18, 19, 20), name = "Region") 

# ggsave(water_scatter, filename = "figs/fig_02b_water_scatter.eps", height = 5, width = 7)

# water supply shortage frequency against level of access to the water network
shortage_scatter <- ggplot(df) +
  geom_point(aes(x = urban, y = water_shortage, colour = region_name, 
                 shape = region_name), size = 3) +
  xlab("% residents living in urban areas") +
  ylab("Frequency of water shortages") +
  ylim(0,1) +
  labs(colour = "Region") +
  scale_color_tableau("Nuriel Stone", direction = 1, name = "Region") +
  scale_shape_manual(values = c(15, 17, 18, 19, 20), name = "Region") 

# ggsave(shortage_scatter, filename = "figs/fig_02c_shortage_scatter.eps", height = 5, width = 7)

# Manuscript Fig 2
# make composite plot of urban map and urban-water scatter
urban_water <- ggarrange(urban_map, water_scatter, shortage_scatter, 
                         ncol = 3, labels = c("a", "b", "c"), 
                         hjust = -1,
                         vjust = 1,
                         font.label = list(size = 14, face = "plain"))

ggsave(urban_water, filename = "figs/fig_02_urban_water.eps", height = 5, width = 20)

# Appendix Fig A1a
# plot maps of 6 biomes
biome_map <-ggplot(map_df) +
  geom_sf(aes(fill = biome_name), lwd = 0) +
  scale_fill_tableau("Summer", direction = -1) +
  labs(fill = "Biome") +
  theme_void()

# ggsave(biome_map, filename = "figs/fig_S01a_biomes.eps", height = 4, width = 5.5)

# Appendix Fig A1b
# plot map of Koppen climate zones
climate_map <-ggplot(map_df) +
  geom_sf(aes(fill = main_climate), lwd = 0) +
  scale_fill_tableau("Green-Orange-Teal", direction = 1) +
  labs(fill = "Climate zone") +
  theme_void()

# ggsave(climate_map, filename = "figs/fig_S01b_climate.eps", height = 4, width = 5.5)

# Appendix Fig A1c
# plot map of 5 geo-political regions 
region_map <- ggplot(map_df) + 
  geom_sf(aes(fill = region_name), lwd = 0) +
  scale_fill_tableau("Nuriel Stone", direction = 1) +
  labs(fill = "Region") + 
  theme_void()

# ggsave(region_map, filename = "figs/fig_S01c_regions.eps", height = 4, width = 5.5)

# Appendix Fig A1
biome_climate_region <- ggarrange(biome_map, climate_map, region_map, 
                          ncol = 3, labels = c("a", "b", "c"), 
                          hjust = -1,
                          vjust = 3,
                          font.label = list(size = 14, face = "plain"))

ggsave(biome_climate_region, filename = "figs/fig_S01_maps.eps", height = 4, width = 16)

# Appendix Fig A2
# plot maps of dengue incidence rate (DIR) per 100,000 population per year 
dir_year <- 
  data %>% 
  group_by(year, micro_code) %>%
  # calculate annual incidence rate
  summarise(cases = sum(dengue_cases),
            pop = sum(population)) %>% 
  mutate(var = cases / pop * 10^5)  %>%
  # add the map
  left_join(map, ., by = c("code" = "micro_code")) %>% 
  ggplot() + 
  geom_sf(aes(fill = var), lwd = 0, color = NA) +
  scale_fill_gradientn(name = "DIR", colours = brewer.pal(9, "PuRd"), 
                       trans = "log1p", breaks = c(0, 10, 100, 300, 1000), 
                       labels = c(0, 10, 100, 300, 1000) ) + 
  theme_void() +
  facet_wrap(~year, ncol = 5)

ggsave("figs/fig_S02_observed_DIR.eps", height = 20, width = 30, units = "cm")

# Appendix Fig A3
# plot Palmer drought severity index heat maps (month and year) per state
pdsi_facet <- 
  data %>% 
  group_by(year, month, state_code) %>%
  # calculate mean PDSI
  summarise( var = mean(pdsi, na.rm = T)) %>% 
  # add the predefined state grid by state code
  left_join(grid, by = c("state_code" = "code_num")) %>% 
  ggplot(aes(x = month, y = year, fill = var)) + 
  geom_raster() +
  ylab("Year") + 
  xlab("Month") +
  scale_fill_gradientn(name = "PDSI", colours = brewer.pal(11, "BrBG")) + 
  scale_y_continuous() +
  scale_x_continuous(breaks = c(1,4,7,10), labels = c("Jan", "Apr", "Jul", "Oct")) +
  theme_bw() + 
  # organise by state name in grid file
  facet_geo( ~name, grid = grid)

ggsave("figs/fig_S03_pdsi_facet.eps", height = 30, width = 25, units = "cm")

# Appendix Fig A4
# plot minimum temperature heat maps (month and year) per state 
tmin_facet <- 
  data %>% 
  group_by(year, month, state_code) %>%
  # calculare mean tmin
  summarise(tmin = mean(tmin)) %>% 
  # add the predefined state grid by state code
  left_join(grid, by = c("state_code" = "code_num")) %>% 
  ggplot(aes(x = month, y = year, fill = tmin)) + 
  geom_raster() +
  ylab("Year") + 
  xlab("Month") +
  scale_fill_gradientn(name = "Tmin", colours = rev(brewer.pal(11, "RdBu"))) + 
  scale_y_continuous() +
  scale_x_continuous(breaks = c(1,4,7,10), labels = c("Jan", "Apr", "Jul", "Oct")) +
  theme_bw() + 
  # organise by state name in grid file
  facet_geo( ~name, grid = grid)

ggsave("figs/fig_S04_tmin_facet.eps", height = 30, width = 25, units = "cm")

# Appendix Fig A5a

# plot % residents with access to the piped water network
access_map <- ggplot(map_df) + 
  geom_sf(aes(fill = water_network), lwd = 0) +
  scale_fill_gradient_tableau("Blue-Teal") +
  labs(fill = "% of residents with \n access to water network") + 
  theme_void()

# ggsave(access_map, filename = "figs/fig_S05a_water_access.eps", height = 4, width = 5)

# Appendix Fig A5b

# plot % residents with access to the piped water network
shortage_map <- ggplot(map_df) + 
  geom_sf(aes(fill = water_shortage), lwd = 0) +
  scale_fill_gradient_tableau("Blue-Teal") +
  labs(fill = "Water shortage frequency") + 
  theme_void()

# ggsave(shortage_map, filename = "figs/fig_S05b_water_shortage.eps", height = 4, width = 5)

# Appendix Fig A5c

# water supply shortage frequency against level of access to the water network
access_shortage_scatter <- ggplot(df) +
  geom_point(aes(x = water_network, y = water_shortage, colour = region_name, 
                 shape = region_name), size = 3) +
  xlab("% residents with access to water network") +
  ylab("Frequency of water shortages") +
  ylim(0,1) +
  labs(colour = "Region") +
  scale_color_tableau("Nuriel Stone", direction = 1, name = "Region") +
  scale_shape_manual(values = c(15, 17, 18, 19, 20), name = "Region") 

# ggsave(access_shortage_scatter, filename = "figs/fig_S05c_access_shortage_scatter.eps", height = 5, width = 7)

# Appendix Fig A5
# make composite plot of water access and water shortage
water <- ggarrange(access_map, shortage_map, access_shortage_scatter, 
                         ncol = 3, labels = c("a", "b", "c"), 
                         hjust = -1,
                         vjust = 1,
                         font.label = list(size = 14, face = "plain"))

ggsave(water, filename = "figs/fig_S05_water.eps", height = 5, width = 20)


