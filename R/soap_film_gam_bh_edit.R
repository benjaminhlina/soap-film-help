# ---- bring in packages ----
{
  library(broom.mixed)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(gratia)
  library(here)
  library(mgcv)
  library(purrr)
  library(rmapshaper)
  library(sf)
  library(soapcheckr)
}
# ---- bring in data ----
n_fish <- fread(here::here("data",
                          "n_fish.csv"),
               header = TRUE)

# ---- bring in shapefile ----

shape <- st_read(here::here("shapefile_new",
                            "siljan_new.shp"))

# ---- select the columns we want -----
n_fish <- n_fish %>%
  dplyr::select(new_receiver_id, longitude,
                latitude, study_year, origin,
                year, week, new_season, month, max_fish, number_of_fish) %>%
  mutate(
    origin = factor(origin)
  )

# ---- filter out just spring ----
n_fish_spring <- n_fish %>%
  filter(new_season == 'spring')

# ---- convert to utms for model ----

n_fish_spring_utm <- n_fish_spring %>%
  st_as_sf(., coords = c("longitude", "latitude"),
           crs = 4326) %>%
  st_transform(crs = 32633) %>%
  as(., "Spatial") %>%
  as_tibble() %>%
  dplyr::rename(x = coords.x1,
                y = coords.x2)

# ---- first model without soap-film ---

# moved where the fixed effect was in the code as I usally do fixed effects
# then smoothers
mod1 <- gam(number_of_fish ~ origin +
              s(x, y, by = origin, k = 40) +
              s(study_year, k = 3),
            offset = log(max_fish),
            data = n_fish_spring_utm,
            family = nb(link = "log"),
            method = "REML"
            )

# I like using gratia to check the model vs. gam.check()
appraise(mod1)
# yep you're right, not great,

# ---- assess main effects ----
anova.gam(mod1)

# ---- draw model ----

draw(mod1, contour = FALSE, scales = 'fixed',
     select = c(1,2)
     )


# ------------------ SOAP-FILM -------------------------

# ---- check our shape ----
ggplot() +
  geom_sf(data = shape)

# shape is too complex need to simplify

shape_simp <- ms_simplify(shape,
                          keep = 0.01,
                          keep_shapes = FALSE,
                          explode = TRUE)

# check simplified shape

ggplot() +
  geom_sf(data = shape_simp)



# ---- we now need to convert to boundary list -----
# this code is from the following post and is in the vigentte for
# soapcheckr
#
# https://blog.benjaminhlina.com/posts/post-with-code/soapcheckr/

# ---- split into multipoin objects ----
shape_simp_ed <- shape_simp %>%
  dplyr::select(geometry) %>%
  st_cast("MULTIPOINT") %>%
  mutate(
    id = 1:nrow(.)
  )

# ---- split into points for each multipoint object ----
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

# create a vector for the length of the list
nr <- 1:length(bnd_pt)

# ---- create our boundary list -----
shap_bnd_ls <- lapply(nr, function(n) as.list.data.frame(bnd_pt[[n]]))

# ---- check with soap checker ----
soap_check(bnd = shap_bnd_ls)

# boundary looks good

# ---- Make knots ----
lake_grid <- shape_simp %>%
  st_make_grid(cellsize = 1000, square = TRUE, what = "centers") %>%
  st_as_sf()

st_geometry(lake_grid) <- "geometry"

# remove knots that fall outside the boundry
lake_intesects <- st_intersection(shape_simp, lake_grid)

# check using ggplot
ggplot() +
  geom_sf(data = shape_simp) +
  geom_sf(data = lake_intesects)

# some knots look close we will use soapcheckr to remove
# but first we need to convert to data.frame
lake_knots <- lake_intesects %>%
  mutate(
    x = st_coordinates(.)[,"X"],
    y = st_coordinates(.)[,"Y"]
  ) %>%
  st_drop_geometry() %>%
  as.data.frame() %>%
  dplyr::select(x, y)

# ---- check knots ----
soap_check(bnd = shap_bnd_ls, knots = lake_knots)

# ---- multiple knots are going to cause issues, id them using autocrunch ----
crunch_ind <- autocruncher(bnd = shap_bnd_ls, knots = lake_knots)
crunch_ind
# remove knots that are problematic
lake_knots <- lake_knots[-crunch_ind, ]

# check knots agian
soap_check(shap_bnd_ls, knots = lake_knots)

# they look great

# ---- check our boundry ----
glimpse(shap_bnd_ls)

# we are missiing one key thing and that is f
# f tells the model what to do when it hits the boundary, we do not
# expect fish to be there so f is 0


# ---- add f ----

shap_bnd_ls <- lapply(nr,
                      function(n)
                        shap_bnd_ls[[n]] <- c(
                          shap_bnd_ls[[n]],
                          list(f = rep(0, length(shap_bnd_ls[[n]]$x))
                          )
                        )
)

glimpse(shap_bnd_ls)
# f is now part of the boundary


# ---- we can now create a soap-film model ----
mod2 <- gam(number_of_fish ~ origin +
              s(x, y, by = origin, k = 40, bs = "so",
                xt = list(bnd = shap_bnd_ls, nmax = 1500)) +
              s(study_year, k = 3),
            offset = log(max_fish),
            data = n_fish_spring_utm,
            family = nb(link = 'log'),
            method = "REML",
            knots = lake_knots)


# ---- lets check the model fit ----
appraise(mod2)

# again model fit isn't great

# ---- check main effects ----
anova.gam(mod2)
# origin as fixed effect with soap-film is not significant
# where with a thin-plate smoother it is
anova.gam(mod1)

# ---- draw using gratia ----
draw(mod2, contour = FALSE, scales = 'fixed',
     select = c(1,2))

# this creats an error that gavin has commented on and will update

# ---- we can still preidct

# ---- create grid to predict -----
lake_pred <- shape_simp %>%
  st_make_grid(cellsize = 500, square = TRUE, what = "centers") %>%
  st_as_sf()

# name geometry
st_geometry(lake_pred) <- "geometry"

# trim predicted points to be inside lake
lake_pred <- st_intersection(lake_pred, shape_simp) %>%
  dplyr::select(geometry)

# convert to tibble
lake_pred_df <- lake_pred %>%
  mutate(
    x = st_coordinates(.)[,"X"],
    y = st_coordinates(.)[,"Y"],
  ) %>%
  st_drop_geometry()

# ---- exapand points that we are going to predict to include origin & yr ----
lake_pred_df <- tidyr::expand_grid(
  lake_pred_df,
  tibble(
    origin = unique(n_fish_spring_utm$origin)
  ),
  tibble(
    study_year = unique(n_fish_spring_utm$study_year)
  )
)

# ---- convert back to sf to check if this makes sense ----
lake_pred_df_sf <- st_as_sf(lake_pred_df, coords = c("x", "y"),
                            crs = 32633)
# ---- plot using ggplot ----
ggplot() +
  geom_sf(data = lake_pred_df_sf) +
  facet_grid(origin ~ study_year)

# points make sense

# ---- predict nubmer of observations using augment ----
pred <- broom.mixed::augment(mod2, newdata = lake_pred_df)
beepr::beep()

# ---- check the number of observations ----
summary(n_fish_spring_utm)
# ---- need to back transform because we used a log scale ----
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
  filter(.fitted <= 9) # added in filter to remove predicted values > 9 since
# max nubmerof obs is 9.


# ---- plot predicted results from the model ----

ggplot() +
  geom_tile(data = pred_1, aes(x = x, y = y, fill = .fitted)) +
  geom_sf(data = shape_simp, fill = NA, colour = "black") +
  facet_grid(origin ~ study_year) +
  scale_fill_viridis_c(name = "Number of\nObservation",
                       trans = "reverse",
  ) +
  theme_void(
    base_size = 15
  ) +
  theme(
    legend.background = element_blank(),
  ) +
  guides(fill = guide_colourbar(
    frame.linewidth = 0.3,
    ticks.colour = 'black',
    frame.colour = 'black')) +
  labs(x = "Longitude",
       y = "Latitude")

