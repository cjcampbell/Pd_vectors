
# Setup -------------------------------------------------------------------
source("~/Pd_vectors/R/00_Setup.R")
library(gridExtra)
library(ggpubr)
library(ggtext)
library(ggnewscale)
library(data.table)

myCRS <- proj_eqd

# Load records ------------------------------------------------------------

myRecords <- list.files(recursive = T, pattern = "SI1 - Records of tree bats in subterranean roosts.csv", full.names = T) %>%
  read.csv() %>%
  # Some manual tidying
  dplyr::mutate(
    address = case_when(
      locality == "30.79, -86.22" ~ "30.79, -86.22",
      Location == "3 miles SW Elco, Alexander County, Illinois" ~ "Alexander County, Illinois, US",
      county != "" ~
        paste(sep = ", ",
              # If I have county info, stop at that level.
              municipality, paste0(county, " County"),  stateProvince, countryCode),
      municipality != "" ~
        paste(sep = ", ", municipality,stateProvince, countryCode),
      TRUE ~
        paste(sep = ", ",
              locality, municipality,stateProvince, countryCode)
    ),
    address = gsub(", , ", ", ", address),
    address = gsub(", , ", ", ", address),
    address = gsub(", , ", ", ", address),
    address = gsub("^, ", "", address),

    living_record = case_when(
      living == 1 ~ "Y",
      recently_living == 1 ~ "Y",
      TRUE ~ "N"
    ),

    inside = case_when(
      Outside.entrance == 1 ~ "N",
      TRUE ~ "Y"
    ),
    YearObserved = as.numeric(as.character(YearObserved)),
    countryCode = case_when(stateProvince == "California" & county == "Inyo" ~ "US", TRUE ~ countryCode),
    stateProvince = case_when(Location == "Carlsbad caverns" ~ "New Mexico", TRUE ~ stateProvince)
  )

## Geocoding of records ---------------------------------------------------------------

rerun_geocoding <- FALSE
if(rerun_geocoding != FALSE){

  if(file.exists( file.path(wd$bin, "coded_addresses.Rdata") ) ) {
    load(file.path(wd$bin, "coded_addresses.Rdata"))
  }
  locations_key <- myRecords %>%
    dplyr::select(address) %>% distinct

  to_look_up <- locations_key %>%
    anti_join(., coded_addresses, by = c("address" = "input_address"))

  # Prep for georeferencing
  require(ggmap)
  if(has_google_key() == FALSE){
    myAPIkey <- read.table("~/googleAPIkey.txt") %>% unlist %>% paste
    register_google(key = myAPIkey)
    if( has_google_key() ) message("Google API key registered.")
  }


  coded <- lapply(to_look_up$address, function(i){
    x <- ggmap::geocode(
      i, output = "latlona", messaging = F,
      override_limit = T, force = T, source = "google"
    )
    while(is.na(x$lon)){
      Sys.sleep(10)
      x <- geocode(
        i, output = "latlona", messaging = F,
        override_limit = T, force = T, source = "google"
      )
    }
    Sys.sleep(1)
    tb <- cbind(ADDRESS = i, x)
    return(tb)
  } )

  coded_addresses_new <- coded %>%
    plyr::ldply() %>%
    dplyr::rename(input_address = ADDRESS, address_geocoded = address)

  coded_addresses <- rbind(coded_addresses_new, coded_addresses)

  fwrite(coded_addresses, file = file.path(wd$bin, "coded_addresses.csv"), row.names = F)

} else {
  coded_addresses <- fread(file.path(wd$bin, "coded_addresses.csv"))
}

myRecords <- left_join(myRecords, coded_addresses, by = c("address" = "input_address")) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(crs = myCRS)

# Make as data frame (compatible with point jitter.)
myRecords_df <- myRecords %>%
  st_coordinates() %>%
  bind_cols(myRecords) %>%
  dplyr::mutate(subGenusOrder = case_when(subGenus == "Lasionycteris" ~ 1, subGenus == "Lasiurus" ~ 2, subGenus == "Aorestes" ~ 3, TRUE ~ 4)) %>%
  arrange(subGenusOrder)

## Exploration of records ------

# Approx year observed
myRecords %>%
  ggplot() +
  aes(YearObserved, fill = subGenus) +
  geom_bar()

# Approx numberRecords
myRecords %>%
  ggplot() +
  aes(subGenus, fill = subGenus) +
  geom_bar()

# Other count info
dplyr::filter(myRecords_df, living_record == "Y") %>%
  count(subGenus, inside)

# Load map data -----------------------------------------------------------

# Set analysis boundaries using the records data
myMapExt <- st_bbox(myRecords) %>%
  st_as_sfc() %>%
  st_buffer(dist = 500e3)

if(!file.exists(file.path(wd$bin, "americasMaps.RData"))) {
    # Get country list.
    americas <- rnaturalearth::countries110 %>%
      st_as_sf %>%
      dplyr::filter(region_un == "Americas")

    # Use level 1 when available, if not use level 0.
    # First, download both.
    for(i in  americas$adm0_a3) {
      #print(i)
      try({ gadm(country = i, level = 0, path = locationGADMdata) })
      try({ gadm(country = i, level = 1, path = locationGADMdata) })
    }

    myCountries_1 <- lapply(americas$adm0_a3, function(code) {
      whichFiles <- list.files(locationGADMdata, pattern = code, recursive = T, full.names = T)
      myRDS <- grep(whichFiles, pattern = "_1_pk.rds", value = T)
      if(length(myRDS) == 0) {
        myRDS <- grep(whichFiles, pattern = "_0_pk.rds", value = T)
      }
      out <- myRDS %>%
        terra::vect() %>%
        sf::st_as_sf()
      if(!is.null(out)) return(out)
    })

    americas_sf_1 <- bind_rows(myCountries_1) %>%
      st_transform(myCRS) %>%
      st_crop(myMapExt) %>%
      st_simplify(dTolerance = 5e3)

    ## Do by level 0 only for country borders
    myCountries_0 <- lapply(americas$adm0_a3, function(code) {
      whichFiles <- list.files(locationGADMdata, pattern = code, recursive = T, full.names = T)
      myRDS <- grep(whichFiles, pattern = "_0_pk.rds", value = T)
      out <- myRDS %>%
        terra::vect() %>%
        sf::st_as_sf()
      if(!is.null(out)) return(out)
    })

    americas_sf_0 <- bind_rows(myCountries_0) %>%
      st_transform(myCRS) %>%
      st_crop(myMapExt) %>%
      st_simplify(dTolerance = 5e3)

    save(americas_sf_1, americas_sf_0, file = file.path(wd$bin, "americasMaps.RData"))

  } else {
    load(file.path(wd$bin, "americasMaps.RData"))
}


# Make large compound plot ------------------------------------------------

ggItems <- list(
  scale_color_manual(values = c("#70AE6E", "#E09F3E", "#A43828", "#335C67"), breaks = c("Dasypterus", "Aeorestes", "Lasiurus", "Lasionycteris")) ,
  scale_fill_manual(values = c("#70AE6E", "#E09F3E", "#A43828", "#335C67"), breaks = c("Dasypterus", "Aeorestes", "Lasiurus", "Lasionycteris")) ,
  theme_void() ,
  theme(
    plot.margin = margin(0,0,0,0),
    axis.title = element_blank(),
    plot.background = element_rect(color = "grey50", fill="white", linewidth=0.25)
  )
)

hawaiiBox <- americas_sf_1 %>%
  dplyr::filter(NAME_1 == "Hawaii") %>%
  st_bbox() %>%
  st_as_sfc() %>%
  st_buffer(dist = 100e3) %>%
  st_bbox()

insetHawaii <- function(p, xmin = -3.0e6, xmax = 3.5e6, ymin = -3.5e6, ymax=NA) {
  suppressMessages({
    p +
      ggItems +
      coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
      annotation_custom(
        ggplotGrob({
          p +
            ggItems +
            coord_sf(ylim = hawaiiBox[c(2,4)], xlim = hawaiiBox[c(1,3)]) +
            theme(
              legend.position = "none",
              plot.background = element_rect(color = "black", fill="white", size=0.25)
            )
        }),
        ymin = -3e6, ymax = -2e6, xmin = -3E6, xmax = -1e6
      )
  })
}


## "Alive, inside hibernation site"
set.seed(420); p1 <- ggplot() +
  geom_sf(americas_sf_1, mapping = aes(), size = 0.1, colour = "grey50", fill = "grey98") +
  geom_sf(americas_sf_0, mapping = aes(), size = 0.25, colour = "grey10", fill = NA) +
  geom_jitter(dplyr::filter(myRecords_df, inside == "Y", living_record == "Y"),
              mapping=aes(x=X,y=Y, fill = subGenus, color = subGenus),
              shape = 21, size = 2,
              alpha = 0.8, width = 10E4, height = 10E4)  # 10km jitter
# "Dead, inside hibernation site"
set.seed(420); p2 <- ggplot() +
  geom_sf(americas_sf_1, mapping = aes(), size = 0.1, colour = "grey50", fill = "grey98") +
  geom_sf(americas_sf_0, mapping = aes(), size = 0.25, colour = "grey10", fill = NA) +
  geom_jitter(dplyr::filter(myRecords_df, inside == "N", living_record == "Y"),
              mapping=aes(x=X,y=Y, fill = subGenus, color = subGenus),
              shape = 21, size = 2,
              alpha = 0.8, width = 10E4, height = 10E4)  # 10km jitter
# "Alive, outside hibernation site"
set.seed(420); p3 <- ggplot() +
  geom_sf(americas_sf_1, mapping = aes(), size = 0.1, colour = "grey50", fill = "grey98") +
  geom_sf(americas_sf_0, mapping = aes(), size = 0.25, colour = "grey10", fill = NA) +
  geom_jitter(dplyr::filter(myRecords_df, living_record == "N"),
              mapping=aes(x=X,y=Y, fill = subGenus, color = subGenus),
              shape = 21, size = 2,
              alpha = 0.8, width = 10E4, height = 10E4)  # 10km jitter

# g_legend <- function(a.gplot){
#   tmp <- ggplot_gtable(ggplot_build(a.gplot))
#   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
#   legend <- tmp$grobs[[leg]]
#   return(legend)}
#
# myLegend <- g_legend(p1)

# legendDeets <- list(
#   theme(
#     legend.position = c(0.85,0.85),
#     legend.title = element_blank(),
#     legend.background = element_rect(fill = "white", color = "grey50")
#     )
# )

p1b <- {insetHawaii(p1) + ggtitle("Alive, inside cave or mine")    + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) }
p2b <- {insetHawaii(p2) + ggtitle("Dead, inside cave or mine")     + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) }
p3b <- {insetHawaii(p3) + ggtitle("Alive, outside cave or mine")   + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) }


pp <- ggarrange(
  plotlist = list(
    p1b,
    p2b,
    p3b
  ),
  ncol = 3, common.legend = T
)
# plot(pp)
ggsave(plot = pp, filename = file.path(wd$figs, "bigMap.png"), units = "in", height = 7.5, width = 13)

# Add bar charts---------

xsub1 <- "Silver-haired<br>(<i>Lasionycteris</i>)"
xsub2 <- "Red<br>(<i>Lasiurus</i>)"
xsub3 <- "Hoary<br>(<i>Aeorestes</i>)"
xsub4 <- "Yellow<br>(<i>Dasypterus</i>)"


deets <- list(
  geom_bar() ,
  geom_text(stat='count', aes(label=..count..), vjust=-1) ,
  scale_fill_manual(values = c("#70AE6E", "#E09F3E", "#A43828", "#335C67"), breaks = c("Dasypterus", "Aeorestes", "Lasiurus", "Lasionycteris")) ,
  scale_x_discrete(
    labels = c(xsub1,
               xsub2,
               xsub3,
               xsub4)
  ) ,
  ggpubr::theme_pubr() ,
  scale_y_continuous(name = "Number of records", limits = c(0,NA), expand  = expansion(mult = c(0, .15) ) ) ,
  xlab("Subgenus") ,
  theme(
    legend.position = "none",
    axis.text.x = element_markdown(size = rel(0.9))
  )
)


b1 <- dplyr::filter(myRecords, inside == "Y", living_record == "Y") %>%
  mutate(subGenus = factor(subGenus, levels = c("Lasionycteris", "Lasiurus", "Aeorestes", "Dasypterus"))) %>%
  ggplot() +
  aes(subGenus, fill = subGenus) +
  deets

b2 <- dplyr::filter(myRecords, living_record == "N") %>%
  mutate(subGenus = factor(subGenus, levels = c("Lasionycteris", "Lasiurus", "Aeorestes", "Dasypterus"))) %>%
  ggplot() +
  aes(subGenus, fill = subGenus) +
  deets

b3 <- dplyr::filter(myRecords, inside == "N") %>%
  mutate(subGenus = factor(subGenus, levels = c("Lasionycteris", "Lasiurus", "Aeorestes", "Dasypterus"))) %>%
  ggplot() +
  aes(subGenus, fill = subGenus) +
  deets

ppp <- ggpubr::ggarrange(
 plotlist = list(
   p1b,
   p2b,
   p3b,
   b1,
   b2,
   b3
 ),
  ncol = 3, nrow = 2, heights = c(1,1), labels = LETTERS
)

# plot(ppp)

ggsave(plot = ppp, filename = file.path(wd$figs, "bigMap_bars.png"), units = "in", height = 8, width = 12.5)


ppp2 <- ggpubr::ggarrange(
  plotlist = list(
    p1b,b1,
    p2b,b2,
    p3b,b3
  ),
  ncol = 2, nrow = 3, widths = c(1,1), labels = LETTERS
)
ggsave(plot = ppp2, filename = file.path(wd$figs, "bigMap_bars2.png"), units = "in", height = 12.5, width = 9)



## Plot individually ----------------------------
insetHawaii(p3, xmax = NA, ymin = NA, ymax = NA) +
  ggtitle("Alive, inside cave or mine") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))


{insetHawaii(p1) + ggtitle("Alive, inside cave or mine")  + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))} %>%
  ggsave(filename = file.path(wd$figs, "map_alive_inside.png"), width = 6, height = 6)
insetHawaii(p2) + ggtitle("Dead, inside cave or mine")  + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
insetHawaii(p3,  xmax = NA, ymin = NA, ymax = NA) + ggtitle("Alive, outside cave or mine") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))




# Timing figs -------------------------------------------------------------

timing_records <- dplyr::filter(myRecords, living_record == "Y") %>%
  dplyr::mutate(tidy_month = case_when(
    month == "Two winters" ~ 1,
    month == "Fall"        ~ 10,
    month == "winter"      ~ 1,
    month == "summer"      ~ 6,
    month == "8-Jul"       ~ 8,  # Should read "8-9",
    month == "8-10"        ~ 9,
    month == "7-9"         ~ 7,
    month == "8-9"         ~ 8,
    month == "1-2"         ~ 1,
    TRUE ~ as.numeric(month)
  )
  ) %>%
  dplyr::filter(stateProvince != "Hawaii")

timingPlotStuff <- list(
  geom_density(
    aes(tidy_month, group = subGenus, fill = subGenus),
    alpha = 0.5, adjust = 2
  ) ,
  scale_x_continuous(
    name = "Month", limits = c(5,16), expand = c(0,0),
    breaks = seq(5,16,1),
    labels = c("May", "June", "July", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr")
  ) ,
  ylab("Density") ,
  ggpubr::theme_pubr() ,
  theme(
    panel.grid = element_blank(),
    legend.position = c(0.75,0.85),
    # text = element_text(size = 15),
    # axis.title.x = element_text(size = 20),
    legend.text.align = 0,
    plot.margin = unit(c(0.1,1,0.1,0.15), "cm"),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )
)

# myRecords %>% filter(stateProvince != "Hawaii", living_record == "Y", inside == "Y") %>% group_by(subGenus, living_record, inside) %>% dplyr::summarise(n=n())

xsub1 <- "Silver-haired (<i>Lasionycteris</i>)<br><i>n</i>=157"
xsub2 <- "Red (<i>Lasiurus</i>)<br><i>n</i>=16"
xsub3 <- "Hoary (<i>Aeorestes</i>)<br><i>n</i>=6"
xsub4 <- "Yellow (<i>Dasypterus</i>)"
o1 <- timing_records %>%
  rbind({mutate(timing_records, tidy_month = tidy_month + 12)}) %>%
  dplyr::filter(inside == "Y") %>%
  ggplot() +
  timingPlotStuff +
  ggtitle("Alive, inside hibernation site") +
  scale_fill_manual(
    values = c("#70AE6E", "#E09F3E", "#A43828", "#335C67"), breaks = c("Dasypterus", "Aeorestes", "Lasiurus", "Lasionycteris"),
    labels = rev(c(xsub1, xsub2, xsub3, xsub4))
  ) +
  theme(
    legend.text = element_markdown()
  )

# myRecords %>% filter(stateProvince != "Hawaii", living_record == "Y", inside == "N") %>% group_by(subGenus, living_record, inside) %>% dplyr::summarise(n=n())

xsub1 <- "Silver-haired (<i>Lasionycteris</i>)<br><i>n</i>=6"
xsub2 <- "Red (<i>Lasiurus</i>)<br><i>n</i>=16"
xsub3 <- "Hoary (<i>Aeorestes</i>)<br><i>n</i>=1"
xsub4 <- "Yellow (<i>Dasypterus</i>)<br><i>n</i>=1"
o2 <- timing_records %>%
  rbind({mutate(timing_records, tidy_month = tidy_month + 12)}) %>%
  dplyr::filter(inside == "N") %>%
  ggplot() +
  timingPlotStuff +
  ggtitle("Alive, outside hibernation site") +
  scale_fill_manual("Subgenus",
    values = c("#70AE6E", "#E09F3E", "#A43828", "#335C67"), breaks = c("Dasypterus", "Aeorestes", "Lasiurus", "Lasionycteris"),
    labels = rev(c(xsub1, xsub2, xsub3, xsub4))
  )+
  theme(
    legend.text = element_markdown()
  )

oo <- arrangeGrob(
  grobs = list(o1, o2),
  ncol = 1, heights = c(1,1)
)

plot(oo)

ggsave(plot = oo, filename = file.path(wd$figs, "timings.png"), units = "in", height = 7, width = (13/2)-0.5)



# alternate timing plots --------------------------------------------------
timing_records %>%
  as.data.frame() %>%
  count(inside, subGenus)

myPlotDeets <- list(
  geom_violin(alpha = 0.2, adjust = 2, draw_quantiles = c(0.05, 0.95)),
  geom_jitter(height = 0.1, width = 0.05 ),
  scale_x_continuous(
    name = "Month", limits = c(1,12), breaks = 1:12#,
    #labels = c("Jan", "Feb", "Mar", "Apr", "May", "June", "July", "Aug", "Sep", "Oct", "Nov", "Dec")
  ) ,
  ggpubr::theme_pubr() ,
  theme(
    legend.text = element_markdown(),
    axis.text.y = element_markdown(),
  #  axis.text.y = element_blank(),
   # axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid = element_blank(),
    panel.grid.major.x = element_line(color = "grey70", linewidth = 0.1),
    legend.position = "none",
    legend.text.align = 0,
    plot.margin = unit(c(0.1,1,0.1,0.15), "cm"),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA)
  )

)

## Inside -----
timing_records %>%
  dplyr::filter(inside == "Y", !is.na(tidy_month)) %>%
  as.data.frame() %>%
  count(inside, subGenus)

xsub1_inside <- "Silver-haired<br>(<i>Lasionycteris</i>)<br><i>n</i>=154"
xsub2_inside <- "Red<br>(<i>Lasiurus</i>)<br><i>n</i>=15"
xsub3_inside <- "Hoary<br>(<i>Aeorestes</i>)<br><i>n</i>=6"
xsub4_inside <- "Yellow<br>(<i>Dasypterus</i>)"

t_inside <- timing_records %>%
  dplyr::mutate(subGenus = factor(subGenus, levels = c("Dasypterus", "Aeorestes", "Lasiurus", "Lasionycteris"))) %>%
  dplyr::filter(inside == "Y", !is.na(tidy_month)) %>%
  ggplot() +
  aes(tidy_month, y = subGenus, color = subGenus, fill = subGenus) +
  scale_color_manual(
    values = c("#70AE6E", "#E09F3E", "#A43828", "#335C67"), breaks = c("Dasypterus", "Aeorestes", "Lasiurus", "Lasionycteris") #,
  # labels = c(xsub1_inside, xsub2_inside, xsub3_inside, xsub4_inside),
  ) +
  scale_fill_manual(
    values = c("#70AE6E", "#E09F3E", "#A43828", "#335C67"), breaks = c("Dasypterus", "Aeorestes", "Lasiurus", "Lasionycteris") #,
    # labels = c(xsub1_inside, xsub2_inside, xsub3_inside, xsub4_inside),
  ) +
  scale_y_discrete(
    drop = TRUE,
    labels = c(xsub1, xsub2, xsub3),
    breaks = c("Lasionycteris", "Lasiurus", "Aeorestes")
  ) +
  myPlotDeets

## Outside ----
timing_records %>%
  dplyr::filter(inside == "N", !is.na(tidy_month)) %>%
  as.data.frame() %>%
  count(inside, subGenus)

xsub1_outside <- "Silver-haired<br>(<i>Lasionycteris</i>)<br><i>n</i>=6"
xsub2_outside <- "Red<br>(<i>Lasiurus</i>)<br><i>n</i>=11"
xsub3_outside <- "Hoary<br>(<i>Aeorestes</i>)<br><i>n</i>=1"
xsub4_outside <- "Yellow<br>(<i>Dasypterus</i>)"

t_outside <- timing_records %>%
  dplyr::mutate(subGenus = factor(subGenus, levels = c("Dasypterus", "Aeorestes", "Lasiurus", "Lasionycteris"))) %>%
  dplyr::filter(inside == "N",  !is.na(tidy_month)) %>%
  ggplot() +
  aes(tidy_month, y = subGenus, color = subGenus, fill = subGenus) +
  scale_color_manual(
    values = c("#70AE6E", "#E09F3E", "#A43828", "#335C67"), breaks = c("Dasypterus", "Aeorestes", "Lasiurus", "Lasionycteris")#,
    # labels = c(xsub1_outside, xsub2_outside, xsub3_outside, xsub4_outside),
  ) +
  scale_fill_manual(
    values = c("#70AE6E", "#E09F3E", "#A43828", "#335C67"), breaks = c("Dasypterus", "Aeorestes", "Lasiurus", "Lasionycteris")#,
    # labels = c(xsub1_outside, xsub2_outside, xsub3_outside, xsub4_outside),
  ) +
  scale_y_discrete(
    drop = T,
    labels = c(xsub1, xsub2, xsub3),
    breaks = c("Lasionycteris", "Lasiurus", "Aeorestes")
  ) +
  myPlotDeets

ggarrange(
  t_inside + ggtitle("Alive, inside cave or mine") ,
  t_outside+ ggtitle("Alive, outside cave or mine") + theme(axis.text.y = element_blank()) ,
  labels = LETTERS) %>%
  ggsave(filename = file.path(wd$figs, "timingPlot.png"), width = 10, height = 4)




### Redo timing plots -----

tp <- timing_records %>%
  dplyr::filter(inside %in% c("Y", "N"),  !is.na(tidy_month)) %>%
  dplyr::mutate(
    subGenus = factor(subGenus, levels = c("Dasypterus", "Aeorestes", "Lasiurus", "Lasionycteris")),
    inside_lab = case_when(inside == "Y" ~ "Alive, inside cave or mine", inside == "N" ~ "Alive, outside cave or mine"),
    inside_lab = factor(inside_lab, levels = c("Alive, inside cave or mine", "Alive, outside cave or mine"))
    ) %>%
  ggplot() +
  aes(tidy_month, y = subGenus, color = subGenus, fill = subGenus) +
  scale_color_manual(
    values = c("#70AE6E", "#E09F3E", "#A43828", "#335C67"), breaks = c("Dasypterus", "Aeorestes", "Lasiurus", "Lasionycteris")#,
    # labels = c(xsub1_outside, xsub2_outside, xsub3_outside, xsub4_outside),
  ) +
  scale_fill_manual(
    values = c("#70AE6E", "#E09F3E", "#A43828", "#335C67"), breaks = c("Dasypterus", "Aeorestes", "Lasiurus", "Lasionycteris")#,
    # labels = c(xsub1_outside, xsub2_outside, xsub3_outside, xsub4_outside),
  ) +
  scale_y_discrete(
    drop = T,
    labels = c(xsub1, xsub2, xsub3),
    breaks = c("Lasionycteris", "Lasiurus", "Aeorestes")
  ) +
  myPlotDeets +
  facet_wrap(~inside_lab) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12)
  )

 ggsave(tp, filename = file.path(wd$figs, "timingPlot2.png"), width = 10, height = 4)




## Quantile values ------

timing_records %>%
  as.data.frame %>%
  dplyr::filter( !is.na(tidy_month)) %>%
  dplyr::mutate(month_offset = case_when(tidy_month < 7 ~ tidy_month + 12, TRUE~ tidy_month)) %>%
  group_by(inside, subGenus) %>%
  dplyr::summarise(
    n = n(),
    q05_noOffset = quantile(tidy_month, 0.05),
    q95_noOffset = quantile(tidy_month, 0.95),
    q05_offset = quantile(month_offset, 0.05),
    q95_offset = quantile(month_offset, 0.95)
    )


# Save all of the above ---------------------------------------------------

myPlots <- list(p1b,
                p2b,
                p3b,
                b1,
                b2,
                b3, t_inside, t_outside)

lapply(1:length(myPlots), function(x) {
              ggsave(plot = myPlots[[x]], filename = file.path(wd$figs, paste0("caveRecordsFig", x, ".png")), width = 5, height = 3)
            })


testPlot <- ggarrange(
  plotlist = list(p1b,b1, t_inside,
                  p2b,b2, ggplot(),
                  p3b,b3, t_outside ),
  ncol = 3, nrow = 3,
  widths = c(1,1.2,1.4))

ggsave(testPlot, filename = file.path(wd$figs, "timing_bigcombo.png"), width = 12, height = 10)


# testPlot2 <- ggarrange(
#   plotlist = list(p1b,
#                   p2b,
#                   p3b,
#                   b1,
#                   b2,
#                   b3, t_inside, ggplot(), t_outside),
#   ncol = 3, nrow = 3,
#   heights = c(1,1,1.2))
#
# ggsave(testPlot2, filename = file.path(wd$figs, "timing_bigcombo2.png"), width = 10, height = 10)

# Timing compass plots ----------------------------------------------------
# Sys.unsetenv("GITHUB_PAT")
# devtools::install_github("clauswilke/relayer")
library(relayer)
pos_Pd <- 4
pos_cav <- 3
pos_Mig <- 2

library(tidyverse)
library(lubridate)
library(relayer)
library(ggnewscale)

timing_records <- timing_records %>%
  dplyr::filter(living_record == "Y") %>%
  dplyr::filter(stateProvince != "Hawaii") %>%
  dplyr::filter(subGenus != "Dasypterus") %>%
  mutate(subGenus = factor(subGenus, levels = c("Lasiurus",
                                                "Aeorestes",
                                                "Lasionycteris")))

Pd_detections <- read_csv(file.path("/Users/cjcampbell/cavesNmines", "data", "Pd_detections.csv") ) %>%
  dplyr::mutate(subGenus = case_when(
    Species == "LABO" ~ "Lasiurus",
    Species == "LACI" ~ "Aeorestes",
    Species == "LANO" ~ "Lasionycteris"
  ))

# Define angles.
myAng <- seq(-20, -340, length.out = 12)

# Plot!
p <- data.frame(
  subGenus = c(rep("Aeorestes", 6), rep("Lasionycteris", 6), rep("Lasiurus", 6)),
  tidy_month = rep(c(3,4,5,8,9,10), 3),
  n = rep("Migration",6*3),
  y = rep(pos_Mig,6*3)
) %>%
  ggplot() +
  # Cave record plotting.
  geom_tile( data = {
    timing_records %>%
      dplyr::filter(inside == "Y") %>%
      group_by(subGenus, tidy_month) %>%
      dplyr::summarise(n = n(), y = pos_cav) %>%
      ungroup
  },
  aes(x = tidy_month-0.5, y = y, fill = n)
  ) +
  scale_color_viridis_c(option = "A") +
  # Migration plotting
  geom_tile(aes(x = tidy_month-0.5, y = y, fill2 = n) ) %>%
  rename_geom_aes(new_aes = c("fill" = "fill2")) +
  # Pd detection plotting
  geom_tile(
    data = {
      Pd_detections %>%
        dplyr::rename(tidy_month = Month) %>%
        group_by(subGenus, tidy_month) %>%
        dplyr::summarise(n = n(), y = pos_Pd) %>%
        ungroup %>% mutate(is_Y = if_else(n > 0, 1, 0))
    },
    aes(x = tidy_month-0.5, y = y, fill3 = is_Y)
  ) %>%
  relayer::rename_geom_aes(new_aes = c("fill" = "fill3")) +
  # Color scaling
  scale_colour_viridis_c(aesthetics = "fill", guide = "colorbar", name = "Cave records", option = "viridis", direction = -1, end = 0.8, begin = 0.2, limits= c(1,65), trans = "log", breaks = c(1,2,5,10,50)) +
  scale_colour_viridis_d(aesthetics = "fill2", guide = "legend", name = NULL, option = "magma", begin = 0.9) +
  scale_colour_viridis_c(aesthetics = "fill3", guide = "legend", name = "Pd detections", option = "magma", direction = -1, breaks = c(1,2), begin = 0.5, end = 0.7) +

  # Plotting details.
  facet_grid(~subGenus, as.table = F) +
  geom_hline(yintercept = c(2:5-0.5), colour = "grey90", size = 0.2) +
  theme_classic()+
  coord_polar() +
  scale_x_continuous(labels =  month.abb[1:12], breaks = 1:12-0.5) +
  scale_y_continuous(limits = c(0,4.5), breaks = 2:5-0.5) +
  theme(
    legend.direction = "horizontal",
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    strip.background = element_blank(),
    # Angle axis text
    axis.text.x = element_text(size = 10, angle = myAng),
    strip.text = element_text(size = 14),
    legend.position = "none"
  )
p

ggsave(p, filename = file.path(wd$figs, "compassPlots.png"), width = 10, height = 7)


# Make key quick ----------------------------------------------------------

# Define angles.
myAng <- seq(-20, -340, length.out = 12)

# Plot!
pkey <- data.frame(
  #subGenus = c(rep("Aeorestes", 6), rep("Lasionycteris", 6), rep("Lasiurus", 6)),
  tidy_month = rep(1,1),
  n = rep(1,1),
  y = rep(pos_Mig,1)
) %>%
  ggplot() +
  # Cave record plotting.
  geom_tile( data = { data.frame(tidy_month = 1, n = 10, y = pos_cav) },
             aes(x = tidy_month-0.5, y = y, fill = n)
  ) +
  scale_color_viridis_c(option = "A") +
  # Migration plotting
  geom_tile(aes(x = tidy_month-0.5, y = y, fill2 = n) ) %>%
  rename_geom_aes(new_aes = c("fill" = "fill2")) +
  # Pd detection plotting
  geom_tile(
    data = { data.frame(tidy_month = 1, is_Y = 1, y = pos_Pd) },
    aes(x = tidy_month-0.5, y = y, fill3 = is_Y)
  ) %>%
  rename_geom_aes(new_aes = c("fill" = "fill3")) +
  # Color scaling
  scale_colour_viridis_c(aesthetics = "fill", guide = "colorbar", name = "Cave records", option = "viridis", direction = -1, end = 0.8, begin = 0.2, limits= c(1,65), trans = "log", breaks = c(1,2,10,50)) +
  scale_colour_viridis_c(aesthetics = "fill2", guide = "legend", name = "Migration", option = "magma", begin = 0.9) +
  scale_colour_viridis_c(aesthetics = "fill3", guide = "legend", name = "Pd detections", option = "magma", direction = -1, breaks = c(1,2), begin = 0.5, end = 0.7) +

  #guides(fill = guide_colorbar(direction = "vertical") )+
  # Plotting details.
  geom_hline(yintercept = c(2:5-0.5), colour = "grey90", size = 0.2) +
  theme_classic()+
  coord_polar() +
  scale_x_continuous(labels =  month.abb[1:12], breaks = 1:12-0.5, limits = c(0,12)) +
  scale_y_continuous(limits = c(0,4.5), breaks = 2:5-0.5) +
  theme(
    legend.direction = "horizontal",
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    strip.background = element_blank(),
    # Angle axis text
    axis.text.x = element_text(size = 10, angle = myAng),
    strip.text = element_text(size = 14),
    legend.position = "right", legend.justification = c(0, 1),
    legend.key.width = unit(1, "cm")
  )
pkey

ggsave(pkey, filename = file.path("figs","compassPlots_key.png"), width = 8, height = 3)
