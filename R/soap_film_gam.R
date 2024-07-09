## ----setup, include=FALSE--------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
library(mgcv)
library(tidyverse)
library(purrr)
library(data.table)
library(terra)
library(sjPlot)
require(sf)
require(sp)
require(raster)
require(gdistance)
require(rgeos)
library(rmapshaper)
library(gratia)
library(viridis)
library(spbabel)
library(soapcheckr)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
n_fish = fread(here::here("jack-help", "gam", "n_fish.csv"), header = TRUE)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
n_fish <- n_fish %>%
  dplyr::select(new_receiver_id, longitude,
                latitude, study_year, origin,
                year, week, new_season, month, max_fish, number_of_fish)

n_fish$origin <- as.factor(n_fish$origin)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
n_fish_spring <- n_fish %>% filter(new_season == 'spring')


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
#Convert to UTM
n_fish_spring_utm <- n_fish_spring %>%
  st_as_sf(., coords=c("longitude", "latitude"), crs=4326) %>%
  st_transform(32633) %>%
  as(., "Spatial") %>%
  as_tibble() %>%
  dplyr::rename(x=coords.x1, y=coords.x2)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
mod1 <- gam(number_of_fish ~ origin + s(x, y, by = origin, k = 40) +
              s(study_year, k = 3) + origin, offset = log(max_fish),
            data = n_fish_spring_utm,
            family=nb(link = "log"),
            method = "REML")


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
# gam.check(mod1) #Not great, but going with this for now
appraise(mod1)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
summary(mod1)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
gratia::draw(mod1, contour = FALSE, scales = 'fixed', select = c(1,2))


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
shape <- st_read(here::here("jack-help",
                            "gam",
                            "shapefile_new",
                            "siljan_new.shp"))
shap_spat <- shape %>%
  st_transform(32633) %>%
  as(., "Spatial") %>%
  sp::geometry(.)

plot(shape, col = "blue")

ggplot() +
  geom_sf(data = shape)

## --------------------------------------------------------------------------------------------------------------------------------------------------------------
shape_simp <- ms_simplify(shape, keep = 0.01, keep_shapes = FALSE, explode = TRUE)

plot(shape_simp, col = "blue")

ggplot() +
  geom_sf(data = shape_simp)
## --------------------------------------------------------------------------------------------------------------------------------------------------------------
shape.xy.aut <- sptable(shape_simp)

shape.xy.aut <- shape.xy.aut %>%
  dplyr::rename(x = x_, y = y_)

# bens' edits
#
shape_simp_ed <- shape_simp %>%
  dplyr::select(geometry) %>%
  st_cast("MULTIPOINT") %>%
  mutate(
    id = 1:nrow(.)
  )

bnd_pt <- shape_simp_ed %>%
  split(.$id) %>%
  purrr::map(~ st_cast(.x, "POINT") %>%
               mutate(
                 x = st_coordinates(.)[,"X"],
                 y = st_coordinates(.)[,"Y"]
               ) %>%
               st_drop_geometry() %>%
               dplyr::select(-id)
  )


nr <- 1:length(bnd_pt)

shap_bnd_ls <- lapply(nr, function(n) as.list.data.frame(bnd_pt[[n]]))


soap_check(bnd = shap_bnd_ls)




# #Make knots from the geographical extent of the observations
# N <- floor((abs((max(n_fish_spring_utm$x)-min(n_fish_spring_utm$x)))/1000)) #every 1000 m
# gx <- seq(min(n_fish_spring_utm$x), max(n_fish_spring_utm$x), length.out = N)
# gy <- seq(min(n_fish_spring_utm$y), max(n_fish_spring_utm$y), length.out = N)
# gp <- expand.grid(gx, gy)
# names(gp) <- c("x","y")
# plot(gp$x, gp$y)
#
# #The GAM needs the border coordinates as a list of lists,
# #where each list describes one border segment or island:
# coords <- shape.xy.aut %>%
#   dplyr::select(x,y,branch_)
# names(coords) <- c("x", "y", "branch_")
# borderlist <- split(coords, coords$branch_)
# names(borderlist)
#
# border.aut <- lapply(borderlist, `[`, c(1,2))
# nr <- seq(1,3)
# border.aut <- lapply(nr, function(n) as.list.data.frame(border.aut[[n]]))

#Make knots
lake_grid <- st_as_sf(shape_simp) %>%
  st_make_grid(cellsize = 1000, square = TRUE, what = "centers") %>%
  st_as_sf()

st_geometry(lake_grid) <- "geometry"

#remove knots that fall outside the boundry
lake_intesects <- st_intersection(shape_simp, lake_grid)

ggplot() +
  geom_sf(data = shape_simp) +
  geom_sf(data = lake_intesects)

lake_knots <- lake_intesects %>%
  mutate(
    lon = st_coordinates(.)[,"X"],
    lat = st_coordinates(.)[,"Y"]
  ) %>%
  st_drop_geometry() %>%
  as.data.frame() %>%
  dplyr::select(lon, lat)

#And then check that the border and knots are in order with the soap_check function
#devtools::install_github("dill/soapcheckr")
library(soapcheckr)

#vignette("how_to_use_soapcheckr")

lake_knots <- lake_knots %>%
  dplyr::mutate(x = lon,
                y = lat) %>%
  dplyr::select(x,y)

#check knots
soap_check(bnd = shap_bnd_ls, knots = lake_knots)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
crunch_ind <- autocruncher(bnd = shap_bnd_ls, knots = lake_knots)
crunch_ind
# remove knots that are problematic
lake_knots <- lake_knots[-crunch_ind, ]
soap_check(shap_bnd_ls, knots = lake_knots)


glimpse(shap_bnd_ls)



shap_bnd_ls <- lapply(nr,
                             function(n)
                               shap_bnd_ls[[n]] <- c(
                                 shap_bnd_ls[[n]],
                                 list(f = rep(0, length(shap_bnd_ls[[n]]$x))
                                 )
                               )
)

glimpse(shap_bnd_ls)
## --------------------------------------------------------------------------------------------------------------------------------------------------------------
mod2 <- gam(number_of_fish ~ origin +
              s(x, y, by = origin, k = 40, bs = "so",
                xt = list(bnd = shap_bnd_ls, nmax = 1500)) +
              s(study_year, k = 3),
            offset = log(max_fish),
            data = n_fish_spring_utm,
            family = nb(link = 'log'),
            method = "REML",
            knots = lake_knots)


## --------------------------------------------------------------------------------------------------------------------------------------------------------------
# gam.check(mod2) #Not great

appraise(mod2)

## --------------------------------------------------------------------------------------------------------------------------------------------------------------
summary(mod2)





## --------------------------------------------------------------------------------------------------------------------------------------------------------------
gratia::draw(mod2, contour = FALSE, scales = 'fixed', select = c(1,2))

## --------------------------------------------------------------------------------------------------------------------------------------------------------------
plot(mod2)


lake_pred <- shape_simp %>%
  st_make_grid(cellsize = 500, square = TRUE, what = "centers") %>%
  st_as_sf()
st_geometry(lake_pred) <- "geometry"


lake_pred <- st_intersection(lake_pred, shape_simp) %>%
  dplyr::select(geometry)

lake_pred_df <- lake_pred %>%
  mutate(
    x = st_coordinates(.)[,"X"],
    y = st_coordinates(.)[,"Y"],
  ) %>%
  st_drop_geometry()


lake_pred_df <- tidyr::expand_grid(
  lake_pred_df,
  tibble(
    origin = unique(n_fish_spring_utm$origin)
  ),
  tibble(
    study_year = unique(n_fish_spring_utm$study_year)
  )
)
lake_pred_df_sf <- st_as_sf(lake_pred_df, coords = c("x", "y"),
                            crs = 32633)

ggplot() +
  # geom_sf(data = shape_simp) +
  geom_sf(data = lake_pred_df_sf) +
  facet_grid(origin ~ study_year)



pred <- broom.mixed::augment(mod2, newdata = lake_pred_df)
beepr::beep()

pred_1 <- pred %>%
  mutate(

    # lower = exp(1) ^ (fit - 1.96 * se.fit),
    #
    # upper = exp(1) ^ (fit + 1.96 * se.fit),
    # fit = exp(1) ^ fit,
    lower = exp(1) ^ (.fitted - 1.96 * .se.fit),
    higher = exp(1) ^ (.fitted + 1.96 * .se.fit),
    .fitted = exp(1) ^ .fitted
  ) %>%
  filter(.fitted<= 9)

summary(n_fish_spring_utm)


ggplot() +
  geom_raster(data = pred_1, aes(x = x, y = y, fill = .fitted)) +
  geom_sf(data = shape_simp, fill = NA, colour = "black") +
  facet_grid(origin ~ study_year) +
  scale_fill_viridis_c(name = "Number of Observation",
                       trans = "reverse",
                       # breaks = rev(seq(0, 60, 15)
                                    # )
  ) +
  theme_void(
    base_size = 15
  ) +
  theme(
    legend.background = element_blank(),

    # legend.position.inside =  = c(0.98, 0.82),
  ) +
  guides(fill = guide_colourbar(
    frame.linewidth = 0.3,
    ticks.colour = 'black',
    frame.colour = 'black')) +
  labs(x = "Longitude",
       y = "Latitude")
#> Warning: A numeric `legend.position` argument in `theme()` was deprecated in ggplot2
#> 3.5.0.
#> ℹ Please use the `legend.position.inside` argument of `theme()` instead.
