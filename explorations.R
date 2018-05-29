library(tidyverse)
library(rethinking)
library(tidybayes)
library(tidybayes.rethinking)
# options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


theme_set(
    theme_minimal() +
        theme(
            panel.grid.major = element_line(size = 0.1),
            panel.grid.minor = element_blank(),
            plot.background = element_rect(fill = "#fdfdfd",
                                           colour = NA),
            text = element_text(family = "Roboto Condensed"),
            axis.title = element_text(hjust = 1),
            plot.title = element_text(hjust = 0),
            plot.subtitle = element_text(hjust = 0)
        )
)

if (!dir.exists("data")) dir.create("data")
if (!file.exists("data/ca-block-groups.geojson"))
    download.file("https://eviction-lab-data-downloads.s3.amazonaws.com/CA/block-groups.geojson",
                  destfile = "data/ca-block-groups.geojson")

if (!file.exists("data/ca-tracts.geojson"))
    download.file("https://eviction-lab-data-downloads.s3.amazonaws.com/CA/tracts.geojson",
                  destfile = "data/ca-tracts.geojson")

if (!file.exists("data/ca-tracts.csv"))
    download.file("https://eviction-lab-data-downloads.s3.amazonaws.com/CA/tracts.csv",
                  destfile = "data/ca-tracts.csv")

if (!file.exists("data/ca-cities.csv"))
    download.file("https://eviction-lab-data-downloads.s3.amazonaws.com/CA/cities.csv",
                  destfile = "data/ca-cities.csv")

ca_tract <- sf::read_sf("data/ca-tracts.geojson")

dict <- "population	poverty-rate	pct-renter-occupied	median-gross-rent	median-household-income	median-property-value	rent-burden	pct-white	pct-af-am	pct-hispanic	pct-am-ind	pct-asian	pct-nh-pi	pct-multiple	pct-other	renter-occupied-households	eviction-filings	evictions	eviction-rate	eviction-filing-rate	imputed	subbed" %>% 
    str_split("\\s") %>% pluck(1)

alameda_tract_evictions <- ca_tract %>% filter(str_detect(pl, "Alameda")) %>% 
    select(contains("e.")) %>% 
    gather(year, evictions, -geometry) %>% 
    mutate(year = str_extract(year, "[0-9]{2}$") %>% 
               paste0("20", .) %>% as.integer)

df_lbls <- ca_tract %>% 
    select(contains(".16")) %>% 
    names %>% 
    str_replace_all("\\.16$", "") %>% 
    "["(1:22)

dict_lbls <- str_replace_all(dict, "-", "_")
dict_lbls <- set_names(dict_lbls, df_lbls)

ct <- ca_tract %>% as.data.frame %>% as.tibble %>% 
    select(GEOID, matches("^[a-z]+\\.[0-9]{2}$")) %>% 
    gather(variable, value, -GEOID) %>% 
    separate(variable, into = c("variable", "year"))

ct <- ct %>% 
    mutate(variable = dict_lbls[variable],
           year = as.integer(paste0("20", year)))

ct <- ct %>% spread(variable, value)

ct <- ca_tract %>% select(GEOID) %>% 
    left_join(ct, by = "GEOID")

ct
## actually:
theme_set(hrbrthemes::theme_ipsum_rc())
ct <- readr::read_csv("data/ca-tracts.csv")
ct <- ct %>% 
    set_names(~str_replace_all(., "-", "_"))

ct <- ct %>% mutate_at(vars(ends_with("rate"), starts_with("pct")), ~./100)

# ca_tract %>% select(GEOID) %>% 
#     inner_join(ct, by = "GEOID")


ct_model_data <- ct %>% 
    filter(!is.na(rent_burden), !is.na(renter_occupied_households),
           !is.na(evictions)) %>% 
    mutate(rent_burden = rent_burden / 100,
           log_renter_occupied_households = log(renter_occupied_households + 1)) %>% 
    mutate(pct_renter = renter_occupied_households / population,
           renter_occupied_households = as.integer(renter_occupied_households))

ct_model_data <- ct_model_data %>% filter(year == 2016)

model1 <- 
    map2stan(
        alist(
            evictions ~ dpois(lambda),
            log(lambda) <- a + (b_rb * rent_burden) + (b_pop * log_renter_occupied_households),
            a ~ dnorm(0, 10),
            b_rb ~ dnorm(0, 10),
            b_pop ~ dnorm(0, 10)
        ),
        data = tidybayes::compose_data(ct_model_data),
        iter = 4000, chains = 4
    )

model1 <- tidybayes::recover_types(model1, ct_model_data)

bloop <- tidy_sim(fit = model1)
bloop <- bloop %>% 
    summarise(eviction_hpdi = list(HPDI(evictions.predicted))) %>% 
    mutate(eviction_low = map_dbl(eviction_hpdi, 1), 
           eviction_hi = map_dbl(eviction_hpdi, 2))
bloop <- bloop %>% select(-eviction_hpdi) %>% 
    ungroup

bloop %>% 
    ungroup %>% 
    sample_n(25, weight = renter_occupied_households) %>% 
    mutate(tempid = seq_len(nrow(.))) %>% 
    ggplot(aes(y = tempid)) + 
    geom_point(aes(x = evictions)) + 
    geom_segment(aes(x = eviction_low, 
                     xend = eviction_hi, yend = tempid), 
                 colour = "grey65")

ct %>% 
    mutate(pcta = pct_af_am > .5) %>% 
    ggplot(aes(x = rent_burden, y = eviction_rate)) + 
    geom_smooth() + 
    scale_x_continuous(limits = c(15, 40)) + 
    facet_wrap(~pcta)

ggplot(bloop, aes(x = rent_burden)) + 
    geom_point(aes(y = evictions), size = .2, position = position_jitter()) +
    geom_ribbon(aes(ymin = eviction_low, ymax = eviction_hi),
                alpha = .4) + 
    scale_x_continuous(limits = c(.15, .4)) +
    scale_y_continuous(limits = c(0, 25))

ct_model_data %>% 
    mutate(afam = cut(pct_af_am, breaks = quantile(pct_af_am), 
                      include.lowest = TRUE),
           rntrs = renter_occupied_households / population,
           rntr_quantile = cut(rntrs, breaks = quantile(rntrs),
                               include.lowest = TRUE)) %>% 
    ggplot(aes(x = eviction_rate + .00001)) + 
    geom_histogram() + 
    scale_x_log10() + 
    facet_wrap(~rntr_quantile)

model2 <- 
    map2stan(
        alist(
            evictions ~  dpois(lambda),
            log(lambda) <- a + (b_rb * rent_burden) + 
                (b_pop * log_renter_occupied_households) +
                (b_rntr * pct_renter) +
                (b_af * pct_af_am),
            a ~ dnorm(0, 10),
            b_rb ~ dnorm(0, 10),
            b_pop ~ dnorm(0, 3),
            b_rntr ~ dnorm(0, 6),
            b_af ~ dnorm(0, 6)
        ),
        data = tidybayes::compose_data(ct_model_data),
        iter = 1000, chains = 2
    )
library(tidybayes); library(tidybayes.rethinking)
simdata2 <- tidy_sim(fit = model2)
simdata2 <- simdata2 %>% 
    summarise(eviction_hpdi = list(HPDI(evictions.predicted))) %>% 
    mutate(eviction_low = map_dbl(eviction_hpdi, 1), 
           eviction_hi = map_dbl(eviction_hpdi, 2))
simdata2 <- simdata2 %>% select(-eviction_hpdi) %>% 
    ungroup

simdata2 %>% 
    ungroup %>% 
    sample_n(25, weight = renter_occupied_households) %>% 
    mutate(tempid = seq_len(nrow(.))) %>% 
    ggplot(aes(y = tempid)) + 
    geom_point(aes(x = evictions)) + 
    geom_segment(aes(x = eviction_low, 
                     xend = eviction_hi, yend = tempid), 
                 colour = "grey65")

model3 <-
    map2stan(
        alist(
            evictions ~  dzipois(p, lambda),
            log(lambda) <- a + (b_rb * rent_burden) +
                (b_pop * log_renter_occupied_households) +
                (b_rntr * pct_renter) +
                (b_af * pct_af_am),
            logit(p) ~ zp,
            zp <- a_zp + (b_zp_af * pct_af_am) + (b_zp_rntr * pct_renter),
            a ~ dnorm(0, 10),
            a_zp ~ dnorm(0, 1),
            b_rb ~ dnorm(0, 10),
            b_pop ~ dnorm(0, 3),
            b_rntr ~ dnorm(0, 6),
            b_af ~ dnorm(0, 6),
            b_zp_af ~ dnorm(0, 1),
            b_zp_rntr ~ dnorm(0, 1)
        ),
        data = ct_model_data_tidy,
        iter = 4000, chains = 2,
        start = list(a = -5, b_rb = 2, b_pop = 1, b_rntr = -2, b_af = 2)
    )

model3 <- recover_types(model3, ct_model_data)
simdata3 <- tidy_sim(fit = model3)

simdata3 %>% 
    ungroup %>% filter(.iteration < 10) %>% 
    mutate(tempid = paste0(.chain, "-", .iteration)) %>% 
    ggplot(aes(y = tempid)) + 
    ggridges::geom_density_ridges(aes(x = evictions.predicted))
    
simdata3 <- simdata3 %>% 
    summarise(eviction_hpdi = list(HPDI(evictions.predicted))) %>% 
    mutate(eviction_low = map_dbl(eviction_hpdi, 1), 
           eviction_hi = map_dbl(eviction_hpdi, 2))

simdata3 <- simdata3 %>% select(-eviction_hpdi) %>% 
    ungroup

simdata3 %>% 
    ungroup %>% 
    sample_n(100, weight = renter_occupied_households) %>% 
    mutate(tempid = seq_len(nrow(.))) %>% 
    ggplot(aes(y = tempid)) + 
    geom_point(aes(x = evictions)) + 
    geom_segment(aes(x = eviction_low, 
                     xend = eviction_hi, yend = tempid), 
                 colour = "grey65")


ct_model_data_tidy <- tidybayes::compose_data(ct_model_data)
ct_model_data_tidy$renter_occupied_households <- ct_model_data$renter_occupied_households

model4 <- map2stan(
    alist(
        evictions ~ dbinom(renter_occupied_households, p),
        logit(p) <- a + (b_af * pct_af_am) + 
            (b_rb * rent_burden) +
            (b_rntr * pct_renter),
        a ~ dnorm(0, 10),
        b_rb ~ dnorm(0, 10),
        b_rntr ~ dnorm(0, 10),
        b_af ~ dnorm(0, 10)
    ), data = ct_model_data_tidy,
    iter = 1000, chains = 2
)

simdata4 <- tidy_sim(fit = model4)

simdata4 %>% 
    ungroup %>% filter(.iteration < 10) %>% 
    mutate(tempid = paste0(.chain, "-", .iteration)) %>% 
    ggplot(aes(y = tempid)) + 
    ggridges::geom_density_ridges(aes(x = evictions.predicted))

simdata4 <- simdata4 %>% 
    summarise(eviction_hpdi = list(HPDI(evictions.predicted))) %>% 
    mutate(eviction_low = map_dbl(eviction_hpdi, 1), 
           eviction_hi = map_dbl(eviction_hpdi, 2)) %>% 
    select(-eviction_hpdi) %>% ungroup

simdata4 %>% 
    sample_n(100, weight = renter_occupied_households) %>% 
    mutate(tempid = seq_len(nrow(.))) %>% 
    ggplot(aes(y = tempid)) + 
    geom_point(aes(x = evictions), size = .3) + 
    geom_segment(aes(x = eviction_low, 
                     xend = eviction_hi, yend = tempid), 
                 colour = "grey65")
##############

model5 <- map2stan(
    alist(
        evictions ~ dbinom(renter_occupied_households, p),
        logit(p) <- a + a_county[parent_location] + (b_af * pct_af_am) + 
            (b_rb * rent_burden) +
            (b_rntr * pct_renter),
        a ~ dnorm(0, 10),
        a_county[parent_location] ~ dnorm(0, sigma_county),
        sigma_county ~ dcauchy(0, 2),
        b_rb ~ dnorm(0, 10),
        b_rntr ~ dnorm(0, 10),
        b_af ~ dnorm(0, 10)
    ), data = ct_model_data_tidy,
    iter = 2000, chains = 4
)
model5 <- recover_types(model5, ct_model_data)
simdata5 <- tidy_sim(fit = model5)

# simdata5 %>% 
#     ungroup %>% filter(.iteration < 10) %>% 
#     mutate(tempid = paste0(.chain, "-", .iteration)) %>% 
#     ggplot(aes(y = tempid)) + 
#     ggridges::geom_density_ridges(aes(x = evictions.predicted))

simdata5 <- simdata5 %>% 
    summarise(eviction_hpdi = list(HPDI(evictions.predicted))) %>% 
    mutate(eviction_low = map_dbl(eviction_hpdi, 1), 
           eviction_hi = map_dbl(eviction_hpdi, 2)) %>% 
    select(-eviction_hpdi) %>% ungroup

simdata5 %>% 
    sample_n(100, weight = renter_occupied_households) %>% 
    mutate(tempid = seq_len(nrow(.))) %>% 
    ggplot(aes(y = tempid)) + 
    geom_segment(aes(x = eviction_low, 
                     xend = eviction_hi, yend = tempid), 
                 colour = "grey65") +
    geom_point(aes(x = evictions), size = .3)

simdata5 %>% 
    ggplot(aes(x = renter_occupied_households)) +
    geom_segment(aes(xend = renter_occupied_households,
                     y = eviction_low, yend = eviction_hi), 
                 colour = "grey65", size = .1, 
                 position = position_jitter(height = 0)) + 
    geom_point(aes(y = evictions), size = .1, position = "jitter") + 
    scale_x_log10() +
    scale_y_log10()

###################
###############

model6 <- map2stan(
    alist(
        evictions ~ dbinom(renter_occupied_households, p),
        logit(p) <- a + 
            a_county[parent_location] + 
            (b_af + b_afc[parent_location]) * pct_af_am + 
            (b_rb + b_rbc[parent_location]) * rent_burden +
            (b_rntr + b_rntrc[parent_location]) * pct_renter,
        a ~ dnorm(0, 10),
        c(a_county, b_afc, b_rbc, b_rntrc)[parent_location] ~ dmvnorm2(c(a_county0, b_afc0, b_rbc0, b_rntrc0), sigma_county, rho),
        a_county0 ~ dnorm(0, 10),
        b_afc0 ~ dnorm(0, 10),
        b_rbc0 ~ dnorm(0, 10),
        b_rntrc0 ~ dnorm(0, 10),
        sigma_county ~ dcauchy(0, 2),
        b_rb ~ dnorm(0, 10),
        b_rntr ~ dnorm(0, 10),
        b_af ~ dnorm(0, 10),
        rho ~ dlkjcorr(2)
    ), data = ct_model_data_tidy,
    iter = 2000, chains = 4
)

model6 <- recover_types(model6, ct_model_data)
simdata6 <- tidy_sim(fit = model6)

simdata6 <- simdata6 %>% 
    summarise(eviction_hpdi = list(HPDI(evictions.predicted))) %>% 
    mutate(eviction_low = map_dbl(eviction_hpdi, 1), 
           eviction_hi = map_dbl(eviction_hpdi, 2)) %>% 
    select(-eviction_hpdi) %>% ungroup

simdata6 %>% 
    sample_n(100, weight = renter_occupied_households) %>% 
    mutate(tempid = seq_len(nrow(.))) %>% 
    ggplot(aes(y = tempid)) + 
    geom_segment(aes(x = eviction_low, 
                     xend = eviction_hi, yend = tempid), 
                 colour = "grey65") +
    geom_point(aes(x = evictions), size = .3)

simdata6 %>% 
    ggplot(aes(x = renter_occupied_households)) +
    geom_segment(aes(xend = renter_occupied_households,
                     y = eviction_low, yend = eviction_hi), 
                 colour = "grey65", size = .1, 
                 position = position_jitter(height = 0)) + 
    geom_point(aes(y = evictions), size = .1, position = "jitter") + 
    scale_x_log10() +
    scale_y_log10()

