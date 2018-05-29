library(tidyverse)
library(sf)
dl_eviction_data("CA", "cities", "csv")
dl_eviction_data("CA", "cities", "geojson")
ca_city <- readr::read_csv("data/ca-cities.csv")
ca_city_geo <- sf::read_sf("data/ca-cities.geojson")
ca_city_geo <- ca_city_geo %>% select(GEOID, geometry) %>% 
    rename(city_geoid = GEOID) %>% 
    distinct
ca_city %>% filter(year == 2016) %>% 
    anti_join(ca_city_geo,
              by = "GEOID")


ca_tract <- readr::read_csv("data/ca-tracts.csv")
ca_tract_geo <- sf::read_sf("data/ca-tracts.geojson")
ca_tract_geo <- ca_tract_geo %>% select(GEOID, geometry) %>% 
    rename(tract_geoid = GEOID)
    distinct
ca_city <- fix_csv_names(ca_city)
ca_tract <- fix_csv_names(ca_tract)

ca_geo_x <- ca_tract_geo %>% st_join(ca_city_geo, join = st_within)

tract2016 <- ca_tract %>%
    filter(year == 2016) %>% 
    inner_join(ca_geo_x, by = c("GEOID" = "tract_geoid"))

tract2016 <- tract2016 %>% select(-geometry)
tract2016

library(brms)
theme_set(hrbrthemes::theme_ipsum_rc())
get_prior(eviction_filings | trials(renter_occupied_households) ~ 
              (1 | city_geoid), data = tract2016, family = binomial())

mod1 <- brms::brm(
    eviction_filings | trials(renter_occupied_households) ~ (1 | city_geoid),
    data = tract2016, family = binomial())

mod2 <- brm(eviction_filings | trials(renter_occupied_households) ~ 
                (1 + poverty_rate | city_geoid), 
            data = tract2016, family = binomial(), iter = 4000, 
            chains = 2, cores = 2)

library(tidybayes)
bloop <- tract2016 %>% filter(!is.na(eviction_filings), !is.na(poverty_rate)) %>% 
    add_predicted_samples(model = mod2)

bloop2 <- bloop %>% 
    median_hdi(.prob = c(0.73, 0.89, 0.97))

bloop2 %>% 
    ungroup %>% 
    select(eviction_filings, conf.low, conf.high, .prob) %>% 
    mutate(y = eviction_filings >= conf.low & 
               eviction_filings <= conf.high) %>% 
    group_by(.prob) %>% 
    summarise(mean(y))

bloop2 %>%
    ungroup %>% 
    mutate(prob_fac = factor(.prob),
           prob_fac = forcats::fct_reorder(prob_fac, .prob, .desc = TRUE)) %>%
    ggplot(aes(x = renter_occupied_households, y = pred)) + 
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = prob_fac)) + 
    scale_fill_brewer(palette = "Greens", direction = -1) + 
    geom_point(aes(y = eviction_filings), size = .1, alpha = .3) + 
    scale_x_log10()




options(tigris_use_cache=TRUE)
options(tigris_class="sf")

towns <- tigris::county_subdivisions("06", county = "001")
alameda_roads <- tigris::roads("06", "001")
ca_county <- tigris::counties(state = "06")
alameda <- ca_county %>% filter(COUNTYFP == "001")

# create town labels by finding the centroid of each town
# ggplot's label functions work better with X/Y dataframes rather 
# than sf objects
town_labels <- towns %>% select(NAME) %>%
    mutate(center = st_centroid(geometry)) %>%
    as.tibble %>%
    mutate(center = map(center, ~st_coordinates(.) %>%
                            as_data_frame)) %>%
    select(NAME, center) %>% unnest()


ca_tract %>% 
    filter(year %>% between(2006, 2011),
           parent_location == "Alameda County, California") %>% 
    inner_join(ca_tract_geo, by = c(GEOID = "tract_geoid")) %>% 
    st_as_sf %>% 
    ggplot() +
    geom_sf(data = alameda, size = .1, fill = NA) +
    geom_sf(aes(fill = log(eviction_rate + .01), 
                alpha = renter_occupied_households), 
            size = 0) +
    # scale_alpha_continuous(range = c(1, 1)) +
    viridis::scale_fill_viridis() +    
    geom_sf(data = alameda_roads %>% filter(RTTYP %in% c("I", "S")),
            size = .2, colour = "black") +
    # coord_sf(xlim = c(-122.3, -122.1), ylim = c(37.75, 37.9)) +
    ggrepel::geom_label_repel(
        data = town_labels,
        aes(x = X, y = Y, label = NAME),
        size = 3, family = "Roboto Condensed",
        label.padding = unit(.1, "lines"), alpha = .7) +
    theme_minimal() +
    theme(panel.grid.major = element_line(size = 0),
          plot.background = element_rect(fill = "#fdfdfd",
                                         colour = NA),
          axis.title = element_blank(),
          text = element_text(family = "Roboto Condensed"),
          axis.text = element_blank(),
          legend.position = "bottom") +
    facet_wrap(~year)