

# Setup -------------------------------------------------------------------

library(geodata)
library(sf)

## Load county data. -----
# https://www.whitenosesyndrome.org/where-is-wns
countyStatus <- read.csv(file.path(wd$data, "wns_county_status", "wns_county_status.csv"))
NoAm0_county <- gadm(country = c("USA", "MEX", "CAN"), level = 2, path = locationGADMdata) %>%
  st_as_sf()

wns_counties <- left_join(NoAm0_county, countyStatus, by = c( "GID_0" = "country",  "NAME_1" = "stateprov", "NAME_2" = "county")) %>%
  dplyr::filter(!is.na(determination), determination  %in% c("WNS Positive", "Pd Positive")) %>%
  st_transform(proj_eqd) %>%
  dplyr::mutate(
    ID = row_number(),
    year = as.numeric(stringr::word(status_season,1,1, sep = "-"))
    )

# Specify point of "jump"
myPt <- wns_counties %>%
  dplyr::filter(status_season == "2015-16", NAME_1 == "Washington") %>%
  st_union() %>%
  st_centroid()

# Load North American map.
NoAm0 <- gadm(country = c("USA", "MEX", "CAN"), level = 1, path = locationGADMdata)
NoAm <- NoAm0 %>%
  st_as_sf() %>%
  st_transform(crs = proj_eqd) %>%
  st_simplify(dTolerance = 1e3) # Simplify to 1km resolution.


# Plot map ----------------------------------------------------------------
lowcol <- "yellow"
midcol <- "red"
highcol <- "darkorchid4"
midpoint <- 2012

myMap <- ggplot() +
  geom_sf(wns_counties, mapping = aes(fill = year)) +
  geom_sf(NoAm, mapping = aes(), fill = NA, color = "grey40", size = 0.5) +
  scale_fill_gradient2(
    low=lowcol, mid = midcol, high = highcol,
    midpoint= midpoint,
    name = "Winter of\nfirst detection",
    breaks = seq(2000,2023, 1),
    labels = paste0(seq(2000,2023, 1), "-", as.numeric(substr(seq(2000,2023, 1), 3,4)) +1)
  ) +
  geom_sf(myPt, mapping = aes(), shape = 1, size = 10 ) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    legend.key.height = unit(1, "in")
  )  +
  coord_sf(
    xlim = c(ext(wns_counties)[1] - 50e3, ext(wns_counties)[2]+ 50e3),
    ylim = c(ext(wns_counties)[3] - 50e3, ext(wns_counties)[4] + 50e3)
  )

ggsave(myMap, filename = file.path(wd$figs, "spreadMap.png"), width = 10, height = 6)



# Plot distances between counties. ----------------------------------------

# Find distances between each county.
distMatrix <- wns_counties %>%
  {st_distance(.,.)} %>%
  data.frame()
colnames(distMatrix) <- wns_counties$ID
distMatrix$to <- wns_counties$ID

# Find minimum distances between countries of different years.
df <- distMatrix %>%
  pivot_longer(cols = -c(to), names_to = "from",values_to = "dist") %>%
  dplyr::mutate(from = as.numeric(from)) %>%
  left_join(., dplyr::select(st_drop_geometry(wns_counties), ID, year), by = c("to"= "ID")) %>%
  rename(year_to = year) %>%
  left_join(., dplyr::select(st_drop_geometry(wns_counties), ID, year), by = c("from"= "ID")) %>%
  rename(year_from = year) %>%
  dplyr::filter(to != from) %>%
  dplyr::select(from, year_from, to, year_to, dist)

minDists <- df %>%
  dplyr::filter(
    year_from < year_to,
    to != from
    ) %>%
  group_by(to) %>%
  arrange(dist) %>%
  slice(1) %>%
  dplyr::mutate(minDist = as.numeric(dist))

myPt2 <-minDists %>%
  ungroup %>%
  arrange(desc(minDist)) %>%
  slice(1)

gainPlot <- ggplot() +
  aes(x = year_to, y = minDist/1e3) +
  geom_point(data = minDists,aes( color = year_to), position = position_jitter(height = 0, width = 0.1)) +
  theme_bw() +
  scale_color_gradient2(
    low=lowcol, mid = midcol, high = highcol,
    midpoint= midpoint,
  ) +
  scale_y_continuous(
    name = "Smallest geographic gain\nfrom prior year (km)",
    limits = c(0,2600),
    breaks = seq(0,3000,500),
    expand = c(0,0)
  ) +
  scale_x_continuous(
    name = "Winter of first detection",
    breaks = seq(2000,2023, 1),
    labels = paste0(seq(2000,2023, 1), "-", as.numeric(substr(seq(2000,2023, 1), 3,4)) +1)
  ) +
  geom_point(
    data = myPt2,
    shape = 1, size = 10
  ) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position="none",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
ggsave(gainPlot, filename = file.path(wd$figs, "gainPlot.png"), width = 10, height = 6)


# Cumulative species confirmed --------------------------------------------
# See https://www.whitenosesyndrome.org/static-page/bats-affected-by-wns for most of the data sourced here

library(tidyverse)
library(viridisLite)

firstDetections <- read.csv("data/species_firstDetections.csv") %>%
  mutate(WNS.Susceptible = factor(WNS.Susceptible)) %>%
  group_by(WNS.Susceptible, Year.Pd.confirmed) %>%
  summarise(numSpp = n()) %>%
  ungroup%>%
  tidyr::complete(Year.Pd.confirmed = 2006:2023, WNS.Susceptible = c("Yes", "No"), fill = list(numSpp = 0)) %>%
  group_by(WNS.Susceptible) %>%
  arrange(Year.Pd.confirmed) %>%
  mutate(cumSpp = cumsum(numSpp))

col_yes <- "#9a031e"
col_no <-  "#fb8b24"

speciesCount <- ggplot(firstDetections) +
  aes(x=Year.Pd.confirmed, y = cumSpp, color = WNS.Susceptible, fill = WNS.Susceptible) +
  geom_col() +
  scale_color_manual(name = "Cumulative number of species confirmed Pd+", breaks  = c("Yes", "No"), values = c(col_yes, col_no), labels = c("Species with\nWNS diagnostic signs", "Species with\nno diagnostic signs")) +
  scale_fill_manual( name = "Cumulative number of species confirmed Pd+", breaks  = c("Yes", "No"), values = c(col_yes, col_no), labels = c("Species with\nWNS diagnostic signs", "Species with\nno diagnostic signs")) +
  scale_y_continuous(name = "Cumulative Species", expand = c(0,0), breaks = seq(0,20, by = 5), limits = c(0,20)) +
  scale_x_continuous(name = "Year",               expand = c(0,0), breaks = seq(2006,2022, by = 1), limits = c(2006, 2023)) +
 # ggtitle("Cumulative number of species confirmed Pd+") +
  theme(legend.position = "bottom") +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(0.37, 0.8),
    legend.background = element_rect(fill = "transparent", color = "transparent")
  )
ggsave(speciesCount, filename = file.path(wd$figs, "speciesCount.png"), width = 10, height = 6)



# Combine -----------------------------------------------------------------

bigPlot <- ggpubr::ggarrange(
  myMap + theme(legend.key.height = unit(1, "cm")),
  gainPlot,
  speciesCount,
  ncol = 1, labels = LETTERS, vjust = c(1.5,0.5, 1)
  )

ggsave(bigPlot, filename = file.path(wd$figs, "bigPlot.png"), width = 5, height = 8)

