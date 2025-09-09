## ------------------------------------------------------------------------------------------------------------------------------------
#| label: load_pkg
#| include: false
#source("~/.Rprofile")
require(NinaR)
require(Norimon) #devtools::install_github("NINAnor/Norimon)
require(DBI)
require(RPostgres)
require(ggplot2)
require(xtable)
require(sf)
require(ggforce)
# require(bookdown)
require(knitr)
require(cowplot)
require(grid)
require(gridExtra)
require(landscapemetrics)#
require(fasterize)
require(sf)
require(tmap)
require(betapart)
require(vegan)
require(PerformanceAnalytics)
require(lme4) 
#install.packages(reformulas)
# require(lmerTest)
require(MuMIn)
require(iNEXT)
require(latex2exp)
require(ggvegan)
require(bookdown)
require(Polychrome)
require(ggthemes)
require(sjPlot)
require(knitr)
require(datawizard)
require(iNEXT)

require(terra)
require(rnaturalearth)
require(mgcv)
require(tidyterra)
require("RColorBrewer")
require(ggeffects)
require(gam.hp)
require(effects)
require(gllvm)
# remotes::install_github("davidsjoberg/ggsankey") 
library(ggsankey)
library(ggpattern)

require(rnaturalearthhires)
suppressPackageStartupMessages(require(tidyverse))


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: setup
#| echo: false
#| message: false
#| warning: false

opts_chunk$set(
  tidy = T,
  tidy.opts = list(width.cutoff = 50),
  echo = F,
  message = F,
  warning = F,
  # dpi = 600,
  # dev = c("pdf", "png", "svg"),
  # fig.path = "figure/",
  cache.extra = 12 # rename this to clear cache
)
options(
  width = 60,
  "xtable.comment" = F,
  "xtable.include.rownames" = F
)

theme_set(theme_few())

# rmarkdown::find_pandoc(version = "2.5")
# regex to change ref format from html to latex (\\@ref)(\()(.*?)(\)) --to--- \\ref{\3}
#


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: cache_setting
#| comment: set up rerun cache for faster renders
cache_var <- TRUE
recalc <- FALSE


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: db_connect
#| comment: requires .pgpass file
connect_to_insect_db(timezone = "Europe/Oslo",
                     bigint = "character")



## ------------------------------------------------------------------------------------------------------------------------------------
#| label: load_localities
#| cache: !expr cache_var
norge <- get_map()

localities <- get_localities(dataset = "NorIns") %>%
  mutate(habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type)) |> 
  filter(year < 2025)

n_localities <- localities %>%
  st_drop_geometry() %>%
  summarize(n_distinct(locality)) %>%
  pull()

localities_point <- localities %>%
  mutate(geom_point = st_centroid(geom)) %>%
  st_set_geometry("geom_point")


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: load_sampling_info
#| cache: !expr cache_var

year_locality <- tbl(
  con,
  Id(
    schema = "events",
    table = "year_locality"
  )
)

locality_sampling <- tbl(
  con,
  Id(
    schema = "events",
    table = "locality_sampling"
  )
)




## ------------------------------------------------------------------------------------------------------------------------------------
#| label: load_cache_biomass

load("out/biomass.RData")




## ------------------------------------------------------------------------------------------------------------------------------------
#| label: load_coarse_richness_year_data
load("out/richness.RData")


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-lokaler20202024
#| fig-cap: Overvåkingslokaliteter i perioden 2020-2024. Hvert år overvåkes 10 lokaliteter per habitatstype og region i et forskjøvet design, for totalt 250 semi-naturlige lokaliteter spredt over hele landet, og 50 lokaliteter i skog på Østlandet.

ggplot(norge) +
  geom_sf() +
  geom_sf(
    data = localities_point,
    aes(color = habitat_type),
    lwd = 2
  ) +
  scale_color_nina(
    name = "Habitat",
    palette = "orange-green"
  )


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: !expr cache_var
nas_ins_lok_1000m_agg_q <- "
SELECT ar5.*, ar5_kode.egenskapsverdi as arealtype
FROM backgrounds.locality_1000m_buffer_arealtype_agg ar5,
lookup.ar5_arealtype ar5_kode,
events.year_locality yl
WHERE ar5.locality_id = yl.locality_id
AND yl.project_short_name = 'NorIns'
AND yl.year < 2025
AND ar5.\"arealressursArealtype\" = ar5_kode.kodeverdi

"
nas_ins_lok_1000m_agg <- read_sf(con,
  query = nas_ins_lok_1000m_agg_q
)




## ------------------------------------------------------------------------------------------------------------------------------------
#| label:  calc_ar5_shannon_index,
#| cache: !expr cache_var
load(file = "out/nas_ins_lok_1000m_agg_brick.Rdata")
nas_ins_lok_1000m_arealtype_shannon_index <- lsm_l_shdi(nas_ins_lok_1000m_agg_brick)
nas_ins_lok_1000m_arealtype_shannon_index$layer <- names(nas_ins_lok_1000m_agg_brick)

locality_ar5_shdi <- nas_ins_lok_1000m_arealtype_shannon_index %>%
  dplyr::select(
    locality = layer,
    lsm_l_shdi = value
  )


## ------------------------------------------------------------------------------------------------------------------------------------
#| label:  calc_ar5_perimeter_length
#| cache: !expr cache_var
locality_ar5_tot_perimeter <- nas_ins_lok_1000m_agg %>%
  mutate(perimeter = lwgeom::st_perimeter_lwgeom(geom)) %>%
  st_drop_geometry() %>%
  group_by(locality) %>%
  summarize(tot_perimeter_length = sum(perimeter) - units::as_units(2 * pi * 1000, "meter")) %>%
  select(
    locality,
    tot_perimeter_length
  ) %>%
  mutate(tot_perimeter_length = units::drop_units(tot_perimeter_length))

# locality_ar5_tot_perimeter %>%
#  print(n = Inf)


## ------------------------------------------------------------------------------------------------------------------------------------
#| label:  calc_ar5_area_sum
#| cache: !expr cache_var

locality_ar5_area <- nas_ins_lok_1000m_agg %>%
  mutate(area = st_area(geom)) %>%
  st_drop_geometry() %>%
  group_by(
    locality,
    arealtype
  ) %>%
  summarize(tot_area = sum(area)) %>%
  pivot_wider(
    id_cols = locality,
    names_from = arealtype,
    values_from = tot_area
  ) %>%
  mutate_if(is.numeric, units::drop_units) %>%
  mutate_if(is.numeric, replace_na, replace = 0)


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: !expr cache_var
loc_types <- localities %>%
  select(
    locality,
    habitat_type,
    year
  )

locality_ar5_area_data_2024 <- locality_ar5_area %>%
  left_join(loc_types,
    by = c("locality" = "locality"),
    copy = T
  ) %>%
  filter(year == 2024) %>% # to look at only 2023, comment out
  ungroup() %>%
  select(-c(
    "locality",
    "habitat_type",
    "Ikke kartlagt",
    "Hav",
    "geom",
    "year"
  ))

locality_ar5_area_group_2024 <- locality_ar5_area %>%
  left_join(loc_types,
    by = c("locality" = "locality"),
    copy = T
  ) %>%
  filter(year == 2024) %>% # to look at only 2023, comment out
  ungroup() %>%
  select(c("locality", "habitat_type")) %>%
  mutate(habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type))


## ------------------------------------------------------------------------------------------------------------------------------------
#| message: false
#| cache: !expr cache_var

locality_ar5_area_data <- locality_ar5_area %>%
  left_join(loc_types,
    by = c("locality" = "locality"),
    copy = T
  ) %>%
  ungroup() %>%
  select(-c(
    "locality",
    "habitat_type",
    "Ikke kartlagt",
    "Hav",
    "geom",
    "year"
  ))

locality_ar5_area_group <- locality_ar5_area %>%
  left_join(loc_types,
    by = c("locality" = "locality"),
    copy = T
  ) %>%
  ungroup() %>%
  mutate(habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type)) |> 
  mutate(habitat_type_year = paste0(habitat_type, "_", year)) |> 
  select(c("locality", "habitat_type_year"))
 


## ------------------------------------------------------------------------------------------------------------------------------------
#| message: false
#| results: hide
example_NMDS <- metaMDS(locality_ar5_area_data_2024, # Our community-by-species matrix
  k = 2
)


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-ar5nmds
#| fig-height: 3.5
#| fig-cap: NMDS over variasjonen i landskapskomposisjon (blå tekst i figuren) i en sirkel med 1km radius rundt lokalitetene besøkt i 2024. "Samferdsel" og "Fulldyrka jord" ligger for nære "Ferskvann" for å skrives ut og og vises som en sirkelsymbol. Det samme gjelder "Åpen fastmark" som ligger nære "Overflatedyrka jord". Nummer i grønt refererer til lokalitsnavnene i semi-naturlig mark og nummer i lilla til lokalitetsnavnene i skog.

vegan::ordiplot(example_NMDS, type = "n")


fixed_loc_names <- locality_ar5_area_group_2024$locality %>%
  gsub("(.*_)(.)", "\\2", .)

orditorp(example_NMDS,
  display = "sites",
  labels = fixed_loc_names,
  col = as.character(factor(locality_ar5_area_group_2024$habitat_type,
    labels = ninaColors()[5:4]
  )),
  air = 0.05,
  cex = 0.7
)

ordihull(example_NMDS,
  groups = locality_ar5_area_group_2024$habitat_type,
  draw = "polygon",
  col = "grey90",
  label = T,
  cex = 1
)

orditorp(example_NMDS,
  display = "species",
  col = ninaColors()[2],
  air = 1.1,
  cex = 0.8
)


## ------------------------------------------------------------------------------------------------------------------------------------
#| message: false
#| results: hide
NMDS_all_years <- metaMDS(locality_ar5_area_data, # Our community-by-species matrix
  k = 2
)


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-ar5nmds_all_years
#| fig-height: 3.5
#| fig-cap: Tilsvarende som @fig-ar5nmds, men for alle lokaliteter i 2020-2024. Skogslokaliteter er tegnet i grønt.

vegan::ordiplot(NMDS_all_years, type = "n")

fixed_loc_names <- locality_ar5_area_group$locality %>%
  gsub("(.*_)(.)", "\\2", .)

ordihull(NMDS_all_years,
  groups = locality_ar5_area_group$habitat_type_year,
  draw = "polygon",
  col = c(rep("grey90", 5), rep("green", 5)),
  label = F,
  cex = 1
)

orditorp(NMDS_all_years,
  display = "species",
  col = ninaColors()[2],
  air = 1.1,
  cex = 0.8
)


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-semi-nat-ar5-plot
#| fig-pos: '!hb'
#| fig-cap: Eksempel på variasjon i landskapskomposisjon (arealtype i AR5) i semi-naturlig mark innen en 1 km buffer. SSB 500x500m rutene er markert i blått og felleplassering i rødt. Semi-nat_143 har lite fulldyrka jord, og mye skog og myr. Semi-nat_126 har mye fulldyrka jord. Semi-nat_112 har mye veier og bebyggelse og er svært heterogent, mens Semi-nat_107 har mye innmarksbeite og overflatedyrka jord.
#| fig-height: 8
#| fig-width: 5
p1 <- locality_ar5_plot("Semi-nat_143",
  legend = F,
  strip = T
)
p2 <- locality_ar5_plot("Semi-nat_126",
  legend = F,
  strip = T
)

p3 <- locality_ar5_plot("Semi-nat_112",
  legend = F,
  strip = T
)

p4 <- locality_ar5_plot("Semi-nat_107",
  legend = F,
  strip = T
)


leg <- locality_ar5_plot("Semi-nat_107",
  legend = T,
  strip = T
) +
  theme(
    legend.position = "bottom",
    # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 7),
    legend.key.width = unit(0.5, "cm")
  ) +
  guides(fill = guide_legend(
    nrow = 3,
    byrow = TRUE
  ))

legend <- as_grob(ggpubr::get_legend(leg))

blank <- ggplot() +
  theme_void()

# grid.arrange(grobs = c(glist_1), ncol = 3, as.table = T)

lm <- matrix(c(1:4, 5, 5),
  ncol = 2,
  byrow = T
)

marrangeGrob(
  list(
    p1,
    p2,
    p3,
    p4,
    legend
  ),
  top = "",
  layout_matrix = lm
)


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-skog-ar5-plot
#| fig-cap: Eksempel på variasjon i landskapskomposisjon (arealtype i AR5) i skogsmark innen en 1 km buffer. SSB 500x500m ruter vises i blått, felleplassering i rødt og de kartlagte ANO-sirklene i lilla (små sirkler). Skog_41 og Skog_44 er relativt homogene og inneholder mest skog, Skog_49 er heterogen med mye dyrket mark, mens Skog_50 består mest av skog og lite bevokst fjell, her klassifisert som 'Åpen fastmark'.
#| fig-height: 8
#| fig-width: 5
p1 <- locality_ar5_plot("Skog_44",
  legend = F,
  strip = T
)

p2 <- locality_ar5_plot("Skog_49",
  legend = F,
  strip = T
)

p3 <- locality_ar5_plot("Skog_41",
  legend = F,
  strip = T
)

p4 <- locality_ar5_plot("Skog_50",
  legend = F,
  strip = T
)

leg <- locality_ar5_plot("Skog_49",
  legend = T,
  strip = T
) +
  theme(
    legend.position = "bottom",
    # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 7),
    legend.key.width = unit(0.5, "cm")
  ) +
  guides(fill = guide_legend(
    nrow = 3,
    byrow = TRUE
  ))

legend <- as_grob(ggpubr::get_legend(leg))

blank <- ggplot() +
  theme_void()

# grid.arrange(grobs = c(glist_1), ncol = 3, as.table = T)

lm <- matrix(c(1:4, 5, 5),
  ncol = 2,
  byrow = T
)

marrangeGrob(
  list(
    p1,
    p2,
    p3,
    p4,
    legend
  ),
  top = "",
  layout_matrix = lm
)


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: get_tot_biomass
tot_biomass <- dbGetQuery(
  con,
  "
                          SELECT sum(st.wet_weight) as sum_vekt,
                          avg(st.wet_weight) as avg_vekt
                          FROM events.sampling_trap st,
                          events.locality_sampling ls,
                          events.year_locality yl
                          WHERE st.locality_sampling_id = ls.id
                          AND ls.year_locality_id = yl.id
                          AND yl.project_short_name = 'NorIns'

                          "
)

biomass_per_loc_year_raw <- get_biomass(
  trap_type = "All",
  agg_level = "year_locality",
  dataset = "NorIns",
  subset_region = NULL
)

# tot_biomass <- get_biomass(trap_type = "All",
#                                          agg_level = "total",
#                                          dataset = "NorIns"
#                                           #  ,subset_region = "Østlandet"
#                                         ) %>%
#   select(sum_wet_weight)

biomass_per_loc_year_mean <- biomass_per_loc_year_raw %>%
  summarise(
    mean_sum_wet_weight = mean(sum_wet_weight),
    sd_sum_wet_weight = sd(sum_wet_weight)
  )

biomass_per_year_q <- "

  SELECT year,  count(distinct st.id)::integer as no_samples,
  count(distinct yl.locality_id) as no_lok,
  min(ls.start_date)::date, max(ls.end_date)::date, round((sum(st.wet_weight) / 1000)::numeric, 2) as kg
                          FROM events.sampling_trap st,
                          events.locality_sampling ls,
                          events.year_locality yl
                          WHERE st.locality_sampling_id = ls.id
                          AND st.received_at_lab IS TRUE
                          AND ls.year_locality_id = yl.id
                          AND yl.project_short_name = 'NorIns'
                          GROUP BY year
                          ORDER BY year


"

biomass_per_year <- dbGetQuery(
  con,
  biomass_per_year_q
) %>%
  select(
    `År` = year,
    `Start` = min,
    `Slutt` = max,
    `Antall lok.` = no_lok,
    `Antall prøver` = no_samples,
    `Våtvekt (kg)` = kg
  ) %>%
  janitor::adorn_totals("row",
    name = "Sum"
  )


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: tbl-biomass-year
#| results: asis
#| tbl-cap: "Total biomasse insekter innsamlet i prosjektet så langt."


if (names(rmarkdown::metadata$format) == "ninareport-pdf") {
  print(xtable(biomass_per_year))
} else {
  kable(biomass_per_year)
}


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: !expr cache_var

no_samplings_per_lok_q <- "
SELECT l.region_name, l.locality,
l.habitat_type,
count(distinct ls.sampling_name)::integer,
min(ls.start_date)::date as start_date,
max(ls.end_date)::date as end_date
FROM events.locality_sampling ls,
events.year_locality yl,
locations.localities l,
events.sampling_trap st
WHERE ls.year_locality_id = yl.id
AND yl.locality_id = l.id
AND yl.year = 2024
AND yl.project_short_name = 'NorIns'
AND st.locality_sampling_id = ls.id
AND st.received_at_lab IS TRUE
GROUP BY region_name, habitat_type,  locality
ORDER BY region_name, locality
"
no_samplings_per_lok <- dbGetQuery(
  con,
  no_samplings_per_lok_q
) %>%
  as_tibble() %>%
  mutate(
    start_date = as.character(start_date),
    end_date = as.character(end_date),
    lok_within_hab = as.integer(gsub("(.*_)(.)", "\\2", locality))
  ) %>%
  arrange(
    habitat_type,
    lok_within_hab
  ) %>%
  select(
    Region = region_name,
    Lokalitet = locality,
    `Antall felletømminger` = count,
    `Start prøvetaking` = start_date,
    `Slutt prøvetaking` = end_date
  )


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: no_samples
#| cache: !expr cache_var
no_samples_q <- "

SELECT trap_model as \"Felletype\",
count(distinct st.sample_name)::integer as \"Antall prosesserte felleprøver\"
FROM events.locality_sampling ls,
events.year_locality yl,
locations.traps,
events.sampling_trap st
WHERE ls.year_locality_id = yl.id
AND st.trap_id = traps.id
AND yl.year = 2024
AND yl.project_short_name = 'NorIns'
AND st.locality_sampling_id = ls.id
AND st.received_at_lab IS TRUE
GROUP BY trap_model
"

no_samples <- dbGetQuery(
  con,
  no_samples_q
) %>%
  mutate(Felletype = ifelse(Felletype == "Malaise_Watkins_Black", "Malaisefelle", "Vindusfelle"))


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: tbl-antall-prover
#| results: asis
#| tbl-cap: Antall prosesserte prøver i 2024. De fire vindusfelleprøvene i hver tømmerunde på hver plass ble slått sammmen og metastrekkodet sammen.
if (names(rmarkdown::metadata$format) == "ninareport-pdf") {
  xtable(no_samples)
} else {
  kable(no_samples)
}

## ------------------------------------------------------------------------------------------------------------------------------------
#| label: tbl-nr-samplings-per-lok
#| results: asis
#| tbl-cap: "Antall felletømminger per lokalitet i 2024."

if (names(rmarkdown::metadata$format) == "ninareport-pdf") {
  print(xtable(no_samplings_per_lok), size = "\\footnotesize")
} else if (names(rmarkdown::metadata$format) == "ninareport-html") {
   kable(no_samplings_per_lok)
  # DT::datatable(no_samplings_per_lok,
  #   style = "bootstrap",
  #   options = list(
  #     pageLength = 10,
  #     autoWidth = TRUE,
  #     columnDefs = list(list(
  #       width = "100px",
  #       targets = 1,
  #       caption = "Antall felletømminger per lokalitet i 2023.",
  #       elementID = "tbl-nr-samplings-per-lok"
  #     ))
  #   )
  # )
} else {
  kable(no_samplings_per_lok)
}


## ------------------------------------------------------------------------------------------------------------------------------------
ano_herb_agg <- tbl(
  con,
  Id(
    schema = "views",
    table = "ano_herb_agg"
  )
) %>%
  filter(project_short_name == "NorIns")


## ------------------------------------------------------------------------------------------------------------------------------------
loc_reg <- localities %>%
  st_drop_geometry() %>%
  select(
    locality,
    region_name
  ) %>%
  mutate(region_name = factor(region_name, levels = c("Østlandet", "Sørlandet", "Trøndelag", "Nord-Norge", "Vestlandet")))


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-ano-herb-spec-nor
#| fig-pos: '!b'
#| out-width: '80%'
#| out-height: '30%'
#| fig-cap: "Dekningsgrad (y-akse) og artsantall av karplanter (størrelse) fra ANO-kartlegging ved overvåkingslokalitetene. For skogshabitat med full ANO-kartlegging viser vi medianverdier fra de 18 rutene som ble målt. Dekningsgrader på over 100 \\% er mulig hvor forskjellige arter vokser over og under hverandre med overlapp, da man registrerer dekningsgrad hvor hver enkelt art. Punktene er noe spredt ut for å kunne skille mellom region (x-akse) og overlappende verdier (y-akse)."

plot_ano_herb_sum(2020:2024)


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: !expr cache_var
year_locality_biodiv_ano <- year_locality_biodiv %>%
  left_join(ano_herb_agg,
    by = c("locality" = "locality"),
    copy = T
  )


## ------------------------------------------------------------------------------------------------------------------------------------
tree_agg <- tbl(
  con,
  Id(
    schema = "views",
    table = "landsskog_agg"
  )
) |> 
filter(year < 2025) %>%
  collect() %>%
  mutate(dom_tree = ifelse(dom_tree == "Dunbjerk", "Dunbjørk", dom_tree)) %>%
  mutate(dom_tree = ifelse(dom_tree == "Hengebjerk", "Hengebjørk", dom_tree)) %>%
  mutate(dom_tree = factor(dom_tree, levels = c(
    "Dunbjørk",
    "Hengebjørk",
    "Bøk",
    "Furu",
    "Gran",
    "Selje"
  ))) %>%
  left_join(loc_reg,
    by = c("locality" = "locality"),
    copy = TRUE
  ) %>%
  select(-c(region_name, project_short_name))


## ------------------------------------------------------------------------------------------------------------------------------------
tree_agg_to_plot <- tree_agg %>%
  group_by(
    year,
    dom_tree
  ) %>%
  summarise(
    avg_avg_age = mean(avg_age),
    sd_avg_age = sd(avg_age, na.rm = TRUE)
  )


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-tree-age-type
#| fig-pos: '!hb'
#| fig-cap: "Dominerende treslag og skogsalder (eldste trær) i den forenklete landsskogstakseringeren gjennom årene 2020-2024. Punktene er noe spredt ut for å kunne skille mellom treslag (x-akse) og overlappende verdier (y-akse)."

ggplot(tree_agg) +
  geom_point(
    aes(
      x = year,
      y = avg_age,
      col = dom_tree
    ),
    alpha = 0.7,
    size = 3,
    position = position_jitterdodge(
      dodge.width = 0.4,
      jitter.height = 2,
      jitter.width = 0
    )
  ) +
  # geom_line(aes(x = year,
  #                y = avg_avg_age,
  #                col = dom_tree),
  #           data = tree_agg_to_plot,
  #           linewidth = 2) +
  scale_x_continuous(breaks = seq(2020, 2024)) +
  scale_color_nina(name = "Dominerende\ntreslag") +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  xlab("År") +
  ylab("Skogsalder")


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-clim-trend-sor
#| fig-cap: "Klimatrend juni-juli på Sørlandet. Øverste figur viser summen av daglig gjennomsnittstemperatur, mens nedre figur viser nedbørsum. Feillinjer viser +- 1 standardavvik basert på alle lokaliteter i regionen."
#| fig-pos: "!h"

clim_trend_plot(
  region = "Sørlandet",
  from_month = 6,
  to_month = 8,
  from_year = 2013,
  to_year = 2024
)


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-clim-trend-vest
#| fig-cap: "Klimatrend juni-juli på Vestlandet. Øverste figur viser summen av daglig gjennomsnittstemperatur, mens nedre figur viser summen av daglig nedbørsum. Feillinjer viser +- 1 standardavvik basert på alle lokaliteter i regionen."

clim_trend_plot(
  region = "Vestlandet",
  from_month = 6,
  to_month = 8,
  from_year = 2013,
  to_year = 2024
)


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-clim-trend-ost
#| fig-cap: "Klimatrend juni-juli på Østlandet. Øverste figur viser summen av daglig gjennomsnittstemperatur, mens nedre figur viser summen av daglig nedbørsum. Feillinjer viser +- 1 standardavvik basert på alle lokaliteter i regionen."

clim_trend_plot(
  region = "Østlandet",
  from_month = 6,
  to_month = 8,
  from_year = 2013,
  to_year = 2024
)


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-clim-trend-trond
#| fig-pos: '!hb'
#| fig-cap: "Klimatrend juni-juli i Trøndelag. Øverste figur viser summen av daglig gjennomsnittstemperatur, mens nedre figur viser summen av daglig nedbørsum. Feillinjer viser +- 1 standardavvik basert på alle lokaliteter i regionen."

clim_trend_plot(
  region = "Trøndelag",
  from_month = 6,
  to_month = 8,
  from_year = 2013,
  to_year = 2024
)


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-clim-trend-nord
#| fig-cap: "Klimatrend juni-juli i Nord-Norge. Øverste figur viser summen av daglig gjennomsnittstemperatur, mens nedre figur viser summen av daglig nedbørsum. Feillinjer viser +- 1 standardavvik basert på alle lokaliteter i regionen."

clim_trend_plot(
  region = "Nord-Norge",
  from_month = 6,
  to_month = 8,
  from_year = 2013,
  to_year = 2024
)


## ------------------------------------------------------------------------------------------------------------------------------------
climdata_semi_nat_105 <- get_climate_data(locality = "Semi-nat_105")


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-sor-temp-2024
#| fig-cap: "Temperatur på Sørlandet 2024, eksemplifisert av lokalitet Semi-nat_105, nær Sandnes."

plot_climate_comparison(climdata_semi_nat_105,
  main_title = FALSE,
  focal_year = 2024,
  rolling_mean = TRUE,
  y_high_limit = 70
)


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-sor-precip-2024
#| fig-cap: "Nedbør på Sørlandet 2024, eksemplifisert av lokalitet Semi-nat_105, nær Sandnes."
plot_climate_comparison(climdata_semi_nat_105,
  main_title = FALSE,
  variable = "precipitation",
  rolling_mean = TRUE,
  focal_year = 2024,
  y_low_limit = -10,
  y_high_limit = 80
)


## ------------------------------------------------------------------------------------------------------------------------------------
climdata_semi_nat_125 <- get_climate_data(locality = "Semi-nat_125")


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-ost-temp-2024
#| fig-pos: '!hb'
#| fig-cap: "Temperatur på Østlandet 2024, eksemplifisert av lokalitet Semi-nat_125, nær Drammen."

plot_climate_comparison(climdata_semi_nat_125,
  main_title = FALSE,
  focal_year = 2024,
  rolling_mean = TRUE
)


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-ost-precip-2024
#| fig-pos: '!hb'
#| fig-cap: "Nedbør på Østlandet 2024, eksemplifisert av lokalitet Semi-nat_125, nær Drammen."
#|
plot_climate_comparison(climdata_semi_nat_125,
  main_title = FALSE,
  variable = "precipitation",
  rolling_mean = TRUE,
  focal_year = 2024,
  y_low_limit = -10
)


## ------------------------------------------------------------------------------------------------------------------------------------
climdata_semi_nat_117 <- get_climate_data(locality = "Semi-nat_117")


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-vest-temp-2024
#| fig-cap: "Temperatur på Vestlandet 2024, eksemplifisert av lokalitet Semi-nat_117, nær Stryn."

plot_climate_comparison(climdata_semi_nat_117,
  main_title = FALSE,
  focal_year = 2024,
  rolling_mean = TRUE
)


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-vest-precip-2024
#| fig-cap: "Nedbør på Vestlandet 2024, eksemplifisert av lokalitet Semi-nat_117, nær Stryn."
plot_climate_comparison(climdata_semi_nat_117,
  main_title = FALSE,
  variable = "precipitation",
  rolling_mean = TRUE,
  focal_year = 2024,
  y_low_limit = -10,
  y_high_limit = 80
)


## ------------------------------------------------------------------------------------------------------------------------------------
climdata_semi_nat_135 <- get_climate_data(locality = "Semi-nat_135")


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-trond-temp-2024
#| fig-pos: '!hb'
#| fig-cap: "Temperatur i Trøndelag 2024, eksemplifisert av lokalitet Semi-nat_135, nær Stjørdal."

plot_climate_comparison(climdata_semi_nat_135,
  main_title = FALSE,
  focal_year = 2024,
  rolling_mean = TRUE
)


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-trond-precip-2024
#| fig-pos: '!hb'
#| fig-cap: "Nedbør i Trøndelag 2024, eksemplifisert av lokalitet Semi-nat_135, nær Stjørdal."

plot_climate_comparison(climdata_semi_nat_135,
  main_title = FALSE,
  variable = "precipitation",
  rolling_mean = TRUE,
  focal_year = 2024,
  y_low_limit = -10
)


## ------------------------------------------------------------------------------------------------------------------------------------
climdata_semi_nat_147 <- get_climate_data(locality = "Semi-nat_147")


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-nord-temp-2024
#| fig-cap: "Temperatur i Nord-Norge 2024, eksemplifisert av lokalitet Semi-nat_147, nær Lakselv."

plot_climate_comparison(climdata_semi_nat_147,
  main_title = FALSE,
  focal_year = 2024,
  y_high_limit = 70,
  rolling_mean = TRUE
)


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-nord-precip-2024
#| fig-cap: "Nedbør i Nord-Norge 2024, eksemplifisert av lokalitet Semi-nat_147, nær Lakselv."

plot_climate_comparison(climdata_semi_nat_147,
  main_title = FALSE,
  variable = "precipitation",
  rolling_mean = TRUE,
  focus_year = 2024,
  y_low_limit = -10
)


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-templogger
#| fig-pos: '!b'
#| out-height: 3in
#| out-width: 2in
#| fig-show: hold
#| fig-cap: "Temperatur- og luftfuktighetslogger av typen Hobo MX2301A. Loggeren er plassert i skygge og data kan hentes ut med blåtann og deretter til en skyløsning."

include_graphics("figures/templogger.jpg")


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-lyslogger
#| fig-pos: '!b'
#| out-height: 3in
#| out-width: 2in
#| fig-show: hold
#| fig-cap: "Temperatur- og lysmåler av typen Hobo MX2202 (grønn knapp). Loggeren er plassert i sol og data kan hentes ut med blåtann og deretter til en skyløsning."
include_graphics("figures/lyslogger_small.png")




## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: !expr cache_var
load(file = "out/locality_sampling_loggerdata.Rdata")


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-logger-corr
#| eval: false
#| fig-cap: "Korrelasjonsplot over loggerdata. Relativ luftfuktighet (%) måles av loggermodell MX2301A og lys (lux) av loggermodell MX2202. MX2301A er montert under en lysskjerm og reagerer ikke like hyppig på sollys som en uskjermete loggeren MX2202 som også måler lys."

# locality_sampling_loggerdata %>%
#   select(
#     Temp_MX2301A = avg_values_MX2301A_temperature,
#     Temp_MX2202 = avg_values_MX2202_temperature,
#     RL = avg_values_MX2301A_rh,
#     Lys = avg_values_MX2202_light
#   ) %>%
#   filter(Lys >= 100) %>%
#   Norimon::chart_correlation(.,
#     histogram = TRUE,
#     method = "pearson",
#     color = ninaColors()[1]
#   )


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: !expr cache_var
locality_sampling_weight_loggerdata <- biomass_mf_locality_sampling_time %>%
  left_join(locality_sampling_loggerdata,
    by = c(
      "sampling_name" = "sampling_name",
      "year" = "year",
      "locality" = "locality",
      "habitat_type" = "habitat_type"
    )
  ) 


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-locality-sampling-biomass-loggerdata
#| fig-height: 3
#| fig-width: 5
#| fig-cap: "Sammenheng mellom temperatur og innsamlet biomasse per enkelt fangstperiode i malaisefeller (2 uker), her vist ved en GAM-modell. Biomasse er avrendt våtvekt delt på antall felledager (gram insekter per dag). Dataene antyder en logaritmisk sammenheng, særlig i semi-naturlig mark."
#|
locality_sampling_weight_loggerdata %>%
  group_by(habitat_type) %>%
  ggplot(aes( # y = log_bm_per_day,
    y = avg_wet_weight / no_trap_days,
    x = avg_values_MX2301A_temperature
  )) +
  geom_point(aes(color = habitat_type),
             alpha = 0.5) +
  scale_color_nina(
    name = "Habitatstype",
    palette = "darkblue-orange"
  ) +
  geom_smooth(aes(color = habitat_type),
    method = "gam",
    se = T,
    linewidth = 2
  ) +
  ylab("Biomasse per felledag") +
  xlab("Middeltemperatur ved felle i fangstperioden")


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: !expr cache_var
locality_sampling_biodiv_loggerdata <- locality_sampling_biodiv %>%
  left_join(locality_sampling_loggerdata,
    by = c(
      "sampling_name" = "sampling_name",
      "year" = "year",
      "locality" = "locality",
      "habitat_type" = "habitat_type"
    )
  ) 


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-locality-sampling-biodiv-loggerdata
#| fig-height: 3
#| fig-width: 5
#| fig-cap: "Sammenheng mellom temperatur og antall arter innsamlet per fangstperiode i malaisefeller, her vist som en lineær funksjon."

locality_sampling_biodiv_loggerdata %>%
  group_by(habitat_type) %>%
  ggplot(aes(
    y = no_species,
    x = avg_values_MX2301A_temperature
  )) +
  geom_point(aes(color = habitat_type),
             alpha = 0.5) +
  scale_color_nina(
    name = "Habitatstype",
    palette = "darkblue-orange"
  ) +
  geom_smooth(aes(color = habitat_type),
    method = "glm",
    se = T,
    linewidth = 2
  ) +
  ylab("Antall arter") +
  xlab("Middeltemperatur ved felle i fangstperioden")


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: !expr cache_var
nas_ins_lok_1000m_skog_q <- "
SELECT skog.*
FROM backgrounds.locality_1000m_buffer_skogtype_agg skog,
events.year_locality yl
WHERE skog.locality_id = yl.locality_id
AND yl.project_short_name = 'NorIns'

"

nas_ins_lok_1000m_skog <- read_sf(con,
  query = nas_ins_lok_1000m_skog_q
)


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: !expr cache_var
locality_skog_area <- nas_ins_lok_1000m_skog %>%
  mutate(area = st_area(geom)) %>%
  st_drop_geometry() %>%
  group_by(
    locality,
    treslag,
    bonitet
  ) %>%
  summarize(tot_area = sum(area)) %>%
  pivot_wider(
    id_cols = locality,
    names_from = c(treslag, bonitet),
    values_from = tot_area
  ) %>%
  mutate_if(is.numeric, units::drop_units) %>%
  mutate_if(is.numeric, replace_na, replace = 0)


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: !expr cache_var
locality_skog_area_data <- locality_skog_area %>%
  ungroup() %>%
  select(-"locality")

skog_rda <- decorana(locality_skog_area_data)

# summary(skog_rda) #hm, too large axis length?


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: !expr cache_var
skog_PCA <- rda(locality_skog_area_data, scale = F)

# summary(skog_PCA)


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: !expr cache_var
skog_scores <- scores(skog_PCA)$sites %>%
  as_tibble() %>%
  cbind(locality_skog_area["locality"]) %>%
  select(locality,
    skog_PC1 = PC1,
    skog_PC2 = PC2
  )

locality_skog_area <- locality_skog_area %>%
  left_join(skog_scores,
    by = c("locality" = "locality")
  )


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: !expr cache_var
locality_ar5_area_data <- locality_ar5_area %>%
  left_join(loc_types,
    by = c("locality" = "locality"),
    copy = T
  ) %>%
  ungroup() %>%
  select(-c(
    "locality",
    "habitat_type",
    "Ikke kartlagt",
    "Hav",
    "geom"
  ))

locality_ar5_area_group <- locality_ar5_area %>%
  left_join(loc_types,
    by = c("locality" = "locality"),
    copy = T
  ) %>%
  ungroup() %>%
  select(c("locality", "habitat_type")) %>%
  mutate(habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type))


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: !expr cache_var
#| label: bind_ar5_scores
# decorana(locality_ar5_area_data) #check axis length

# ar5.var <- apply (locality_ar5_area_data, 2, var)
# ar5.var <- ar5.var/sum(ar5.var)

# sort(ar5.var, decreasing = T)

ar5_PCA_norm <- rda(locality_ar5_area_data, scale = T)
# summary(ar5_PCA_norm)

ar5_PCA <- rda(locality_ar5_area_data, scale = F)
# summary(ar5_PCA)

ar5_scores <- scores(ar5_PCA)$sites %>%
  as_tibble() %>%
  cbind(locality_ar5_area["locality"]) %>%
  select(locality,
    ar5_PC1 = PC1,
    ar5_PC2 = PC2
  )


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: make_ar5_pca_plot
ar5_PCA_plot <- autoplot(ar5_PCA, legend.position = "none") +
  scale_color_nina(palette = "orange-green") +
  xlim(c(-2000, 3000))


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: make_skog_pca_plot
skog_PCA_plot <- autoplot(skog_PCA,
  legend.position = "none"
) +
  scale_color_nina(palette = "orange-green") +
  xlim(c(-2500, 3000)) +
  ylim(c(-1500, 1500))


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-ar5-PCA
#| fig-pos: '!bh'
#| fig-cap: "PCA-analyse over landskapskomposisjon i arealtype (AR5). Verdiene er ikke normalisert, da alle er fra samme datakilde og har samme skala (km2). 'Bebygd', 'Overflatedyrka jord', 'Samferdsel', 'Innmarksbeite', 'Myr' og 'Ferskvann' overlapper i figuren fordi de bidrar med relativt like areal i de undersøkte lokalitetene."
ar5_PCA_plot


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-skog-PCA
#| fig-pos: '!hb'
#| fig-cap: "PCA-analyse over landskapskomposisjon i skogtype og bonitet ('skogtype'_'bonitet') etter AR5 (basert på kartdata fra SR16). Verdiene er ikke normaliserte, da alle er fra samme datakilde og har samme skala (km2). Første aksen gjenspeiler først og fremst variasjon i treslag (fra barskog til ikke tresatt), mens den andre aksen gjenspeiler bonitet (fra impediment til høy). De fleste kombinasjoner har lik mengde og overlapper derfor i figuren. Treslag og bonitet klassifiseres kun på skog, myr og tresatt innmarksbeite, på øvrige arealtyper brukes verdiet 'Ikke_relevant'."
skog_PCA_plot


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: !expr cache_var
loc_reg_year <- localities %>%
  st_drop_geometry() %>%
  select(year, locality, region_name, habitat_type) %>%
  mutate(region_name = factor(region_name, levels = c("Østlandet", "Sørlandet", "Trøndelag", "Nord-Norge", "Vestlandet")))


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-ar5-shannon-plot
#| fig-pos: '!hb'
#| fig-cap: "Fordeling av Shannon indeks av arealtype fra kartlaget AR5 (omtalt som AR5-arealtype_Shannon_indeks i de statistiske modellene nedenfor), innen en buffer med radius 1000m rundt lokalitetene."

locality_ar5_shdi %>%
  left_join(loc_reg_year,
    by = c("locality" = "locality")
  ) %>%
  mutate(habitat_type = ifelse(grepl("Semi-nat", locality), "Semi-nat", "Skog")) %>%
  mutate(habitat_type = factor(habitat_type, levels = c("Semi-nat", "Skog"))) %>%
  group_by(region_name) %>%
  ggplot(.) +
  geom_point(
    aes(
      x = year,
      y = lsm_l_shdi,
      color = region_name,
      pch = habitat_type
    ),
    size = 3,
    alpha = 0.7,
    position = position_jitterdodge(
      jitter.height = 0.1,
      jitter.width = 0,
      dodge.width = 0.3
    )
  ) +
  scale_color_nina(name = "Region") +
  scale_x_continuous(breaks = seq(2020, 2024)) +
  guides(
    color = guide_legend(override.aes = list(size = 3)),
    pch = guide_legend(
      title = "Habitattype",
      override.aes = list(size = 3)
    )
  ) +
  xlab("År") +
  ylab("Shannon index")


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-ar5-perimeter-plot
#| fig-pos: '!h'
#| fig-cap: "Fordeling av total lengde av grenser mellom arealtyper, innen et buffetområde på 1 km rundt lokalitetene. Lanskapet med desidert mest liniære element etter denne definisjon er 'Semi-nat_112', vist i @fig-semi-nat-ar5-plot."
locality_ar5_tot_perimeter %>%
  left_join(loc_reg_year,
    by = c("locality" = "locality")
  ) %>%
  mutate(habitat_type = ifelse(grepl("Semi-nat", locality), "Semi-nat", "Skog")) %>%
  mutate(habitat_type = factor(habitat_type, levels = c("Semi-nat", "Skog"))) %>%
  group_by(region_name) %>%
  ggplot(.) +
  geom_point(
    aes(
      x = year,
      y = tot_perimeter_length,
      color = region_name,
      pch = habitat_type
    ),
    size = 3,
    alpha = 0.7,
    position = position_jitterdodge(
      jitter.height = 0.1,
      jitter.width = 0,
      dodge.width = 0.3
    )
  ) +
  scale_color_nina(name = "Region") +
  scale_x_continuous(breaks = seq(2020, 2024)) +
  guides(
    color = guide_legend(override.aes = list(size = 3)),
    pch = guide_legend(
      title = "Habitattype",
      override.aes = list(size = 3)
    )
  ) +
  xlab("År") +
  ylab("Tot. lengde lineære element")




## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: !expr cache_var
sampling_week_q <- "
SELECT sampling_name, date_part('week', end_date) as ukenummer
FROM events.locality_sampling ls
"

sampling_week <- dbGetQuery(
  con,
  sampling_week_q
) %>%
  select(
    Ukenummer = ukenummer,
    everything()
  )






## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: !expr cache_var

locality_sampling_biomass_temp_compl_light_week <- read_csv(file = "out/locality_sampling_biomass_temp_compl_light_week.csv") %>%
  mutate(
    year_num = Year - 2020,
    region_name = factor(region_name, levels = c("Østlandet", "Sørlandet", "Trøndelag", "Nord-Norge", "Vestlandet")),
    habitat_type = factor(habitat_type, levels = c("Skog", "Semi-nat"))
  )



## ------------------------------------------------------------------------------------------------------------------------------------
locality_ar5_area <- locality_ar5_area %>%
  left_join(ar5_scores,
    by = c("locality" = "locality")
  )




## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: !expr cache_var
locality_sampling_biomass_big_join <- read_csv(file = "out/locality_sampling_biomass_big_join.csv")

locality_sampling_biomass_big_join <- locality_sampling_biomass_big_join %>%
  mutate(region_name = factor(region_name, levels = c("Østlandet", "Sørlandet", "Trøndelag", "Nord-Norge", "Vestlandet")))

locality_sampling_biomass_big_join_forest <- locality_sampling_biomass_big_join %>%
  filter(habitat_type == "Skog") 

locality_sampling_biomass_big_join_semi_nat <- locality_sampling_biomass_big_join %>%
  filter(habitat_type == "Semi-nat")


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: tbl-refbase
#| results: asis
#| tbl-cap: "Antall arter i referansebasene som er blitt brukt i prosjektet så langt. I 2020 tok vi utgangspunkt i en base fra Porter, men i 2021 og 2022 har vi utvidet denne med en betydlig mengde norske arter. Analyser i denne rapport bruker samme versjon av referansedatabasen som i 2022."
#|

refbase_tab <- tibble(
  Klasse = c(
    "Arachnida (edderkoppdyr)",
    "Chilopoda (skolependere)",
    "Collembola (spretthaler)",
    "Diplopoda (tusenbein)",
    "Insecta (insekter)",
    "Ikke målgruppe"
  ),
  `COIClassifier v4` = c(
    4778,
    135,
    268,
    246,
    65660,
    39693
  ),
  `NorInvert v1.1` = c(
    4952,
    147,
    377,
    246,
    69426,
    39695
  ),
  `NorInvert v1.2` = c(
    5124,
    147,
    388,
    246,
    73211,
    39698
  ),
  `NorInvert v1.3` = c(
    5124,
    147,
    388,
    246,
    72151,
    39698
  ),
    `NorInvert v2.0` = c(
    10682,
    251,
    570,
    565,
    110729,
    62583
  )
)

if (names(rmarkdown::metadata$format) == "ninareport-pdf") {
  print(xtable(refbase_tab, digits = 0), size = "footnotesize")
} else {
  kable(refbase_tab, digits = 0)
}


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: make_coverage_table
#| cache: !expr cache_var
#| comment: db_coverage

nt.inv = read.table("../../rawData/norins_metabarcoding/Nortaxa_13022025_inverts.csv", 
                  header = T,
                  sep = ";",
                  comment.char = "")


norwegian_spec <- nt.inv |> 
  filter(ExistsInNorway) |> 
  group_by(Class, Order) |> 
  summarize(`Arter påvist i Norge` = n()) 
  
in_db <- nt.inv |> 
  filter(ExistsInNorway & NorInv2) |> 
  group_by(Class, Order) |> 
  summarize(`av hvilke i database` = n()) 

sumtab_2 <- norwegian_spec |> 
  left_join(in_db,
             by = c("Class" = "Class",
                    "Order" = "Order")) |> 
  mutate(across(where(is.integer), ~ replace_na(.x, 0))) |> 
  mutate(`% i db` = round(`av hvilke i database` / `Arter påvist i Norge` * 100, 2))



## ------------------------------------------------------------------------------------------------------------------------------------
#| label: tbl-coverage_table
#| results: asis
#| tbl-cap: "Antall kjente arter fra Norge (fra offentlige kilder), samt antall arter med strekkoder i referansebasen NorInvert v2.0."
#| tbl-pos: 'ht'
#| layout-ncol: 1

if (names(rmarkdown::metadata$format) == "ninareport-pdf") {
  print(xtable(sumtab_2, digits = 1), size = "small")
} else {
  kable(sumtab_2, digits = 1)
}


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: load_sankey_data
#| cache: !expr cache_var

track21data <- read.table("../../../../store_filer_for_rapport/sannkey_data_tracking2021seqs_13052025.txt",
  header = T,
  sep = "\t",
  fileEncoding = "latin1"
)

figdat.taxonomic_level <- track21data |>
  ggsankey::make_long(classification_level21,
    classification_level22,
    classification_level23,
    classification_level24,
    value = abundance
  )

figdat.taxonomic_level$node <- factor(figdat.taxonomic_level$node,
  levels = c(
    "non-target",
    "class",
    "order",
    "family",
    "genus",
    "species"
  )
)


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: mk_sankey_fig
#| cache: !expr cache_var

taxonomy.seq_level <- ggplot(figdat.taxonomic_level, aes(
  x = x,
  next_x = next_x,
  node = node,
  next_node = next_node,
  fill = factor(node), value = value
)) +
  geom_alluvial(
    flow.alpha = .6,
    node.color = "grey30",
    width = 0.2
  ) +
  scale_fill_manual(
    drop = FALSE,
    values = c(
      "#A2AAAD",
      "#004F71",
      "#008C95",
      "#E57200",
      "#93328E",
      "#7A9A01",
      "#2DCCD3",
      "#FFB25B"
    ),
    name = "Taksonomisk\nrank",
    labels = c(
      "IiFG",
      "klasse",
      "orden",
      "familie",
      "slekt",
      "art"
    )
  ) +
  scale_x_discrete(labels = c(
    "2021\nNorInvert v1.0",
    "2022\nNorInvert v1.2",
    "2023\nNorInvert v1.3",
    "2024\nNorInvert v2.0"
  )) +
  theme(
    axis.title.x = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "grey"),
    axis.text.x = element_text(size = 7)
  ) +
  ylab("Antall sekvenser")



## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-seq_level
#| cache: !expr cache_var
#| fig-cap: "Sankey diagram som viser hvordan oppgradering av referansedatabasen har påvirket taksonomisk klassifisiering av dataen fra 2021.  Kolonnene viser andelen av sekvensene som ble klassifisert til forskjellige taksonomisk nivåer og fargete strømninger viser endringer år over år. IiFG = ikke i fokus gruppen."
#cache force
taxonomy.seq_level



## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-dna-flow-1
#| fig-cap: "Flytskjema for identifisering av målgruppe sekvensene og kvalitetssikring av vår referansedatabase for kontaminasjon og referansesekvensmangel."
#| out-width: '100%'
#| out-height: '90%'

if (names(rmarkdown::metadata$format) == "ninareport-pdf") {
  include_graphics("figures/flyt_2024_1.pdf")
} else {
  include_graphics("figures/flyt_2024_1.png")
}


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-dna-flow-2
#| fig-cap: "Flytskjema for klassifisering og kvalitetssikring av arter basert på DNA-metastrekkoding. Kun arter med høy konfidens er tatt med i de analyser som avhenger av artsidentitet i denne rapporten."
#| out-width: '100%'
#| out-height: '90%'

if (names(rmarkdown::metadata$format) == "ninareport-pdf") {
  include_graphics("figures/flyt_2024_2.pdf")
} else {
  include_graphics("figures/flyt_2024_2.png")
}


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-dna-flow-3
#| fig-cap: "Flytskjema for horisontskanning for fremmede arter og nye arter for Norge. *I tillegg til å kryssjekke arter mot den norske fremmedartslista har vi også kryssjekket artsfunn i dette prosjektet mot databasene Daisie, GRIIS, Saul Pathways, Roy og List of Union Concern. Arter som er listet i en av disse databasene, er antatt å være fremmede arter i et eller flere land i Europa."
#| out-width: '90%'
#| out-height: '90%'
#| fig-pos: '!h'

if (names(rmarkdown::metadata$format) == "ninareport-pdf") {
  include_graphics("figures/flyt_3.pdf")
} else {
  include_graphics("figures/flyt_3.png")
}


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: true
rank_abund_raw_q <- "
SELECT (row_number() over(order by abund DESC))::integer as row_number,
(rank() over(order by abund DESC))::integer as rank,
taxonomic_level,
id_order as order,
species_latin,
abund
FROM

(SELECT
taxonomic_level, id_order, species_latin, sum(no_reads)::integer as abund
FROM occurrences.observations o,
events.identifications i,
events.sampling_trap st,
events.locality_sampling ls,
events.year_locality yl,
locations.traps
WHERE o.identification_id = i.id
AND i.sampling_trap_id = st.id
AND st.locality_sampling_id = ls.id
AND ls.year_locality_id = yl.id
AND st.trap_id = traps.id
AND yl.project_short_name = 'NorIns'
AND trap_short_name != 'MF2'
GROUP BY taxonomic_level, id_order, species_latin
) foo
ORDER BY row_number, rank

"

rank_abund_raw <- dbGetQuery(
  con,
  rank_abund_raw_q
) %>%
  as_tibble()


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: true
rank_abund <- rank_abund_raw %>%
  arrange(row_number, rank) %>%
  select(
    row_number,
    abund,
    taxonomic_level,
    species_latin
  ) %>%
  pivot_wider(
    names_from = taxonomic_level,
    values_from = abund,
    values_fill = 0
  ) %>%
  mutate(
    count_species = cumsum(species > 0),
    count_genus = cumsum(genus > 0),
    count_family = cumsum(family > 0),
    count_order = cumsum(order > 0),
    count_class = cumsum(class > 0),
    all_abund = species + family + order + genus + class
  ) %>%
  mutate(cumsum_abund = cumsum(as.numeric(all_abund))) %>%
  mutate(pct_all_sequences = cumsum_abund / sum(all_abund))

rank_abund_long <- rank_abund %>%
  select(
    rank_abundance = row_number,
    species_latin,
    Art = count_species,
    Slekt = count_genus,
    Familje = count_family,
    Orden = count_order,
    Klasse = count_class
  ) %>%
  pivot_longer(
    cols = c(
      "Art",
      "Slekt",
      "Familje",
      "Orden",
      "Klasse"
    ),
    names_to = "taxonomic_level"
  ) %>%
  mutate(taxonomic_level = factor(taxonomic_level, levels = c("Klasse", "Orden", "Familje", "Slekt", "Art"))) %>%
  arrange(
    taxonomic_level,
    value
  )


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: true
rank_abund_order <- rank_abund_raw %>%
  arrange(row_number, rank) %>%
  group_by(order) %>%
  mutate(order_rank = rank(rank, ties.method = "first")) %>%
  mutate(to_species = as.integer(taxonomic_level == "species")) %>%
  filter(order %in% c("Diptera", "Hymenoptera", "Coleoptera", "Lepidoptera", "Hemiptera")) %>%
  select(
    rank,
    order,
    species_latin,
    order_rank,
    to_species
  )


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: true
order_spec_mod <- glm(to_species ~ order_rank * order, family = binomial(), data = rank_abund_order)
# summary(order_spec_mod)


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: true
new_data_limits <- rank_abund_order %>%
  group_by(order) %>%
  summarise(no_otu = max(order_rank))


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: true

# Ugly, could make an apply version

new_data <- tibble(
  order = "Diptera",
  order_rank = seq(1, new_data_limits$no_otu[new_data_limits$order == "Diptera"], by = 50)
) %>%
  rbind(tibble(
    order = "Hymenoptera",
    order_rank = seq(1, new_data_limits$no_otu[new_data_limits$order == "Hymenoptera"], by = 50)
  )) %>%
  rbind(tibble(
    order = "Coleoptera",
    order_rank = seq(1, new_data_limits$no_otu[new_data_limits$order == "Coleoptera"], by = 50)
  )) %>%
  rbind(tibble(
    order = "Lepidoptera",
    order_rank = seq(1, new_data_limits$no_otu[new_data_limits$order == "Lepidoptera"], by = 50)
  )) %>%
  rbind(tibble(
    order = "Hemiptera",
    order_rank = seq(1, new_data_limits$no_otu[new_data_limits$order == "Hemiptera"], by = 50)
  ))


pred_order_id <- predict(order_spec_mod, newdata = new_data, type = "response")
new_data$pred <- pred_order_id

new_data <- new_data %>%
  mutate(order = factor(order, levels = c("Lepidoptera", "Hymenoptera", "Diptera", "Coleoptera", "Hemiptera")))

new_data <- new_data %>%
  group_by(order) %>%
  mutate(rel_order_rank = rank(rev(order_rank)) / max(rank(rev(order_rank))))


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-artsniva
#| fig-pos: '!ht'
#| fig-cap: "**A)** Fordeling over hvor nøye vi kan sortere de operasjonelle 'artene' (OTUs) til ulike taksonomiske nivåer, fra de som avgitt mye DNA til venstre (vanlige arter) til de som avgitt lite DNA til høyre (sjeldne arter). Sett over hele datasettet kan vi sette et artsnavn på litt unner 50 % av artene, mens andelen er over 75 % for de arter som tilsammen bidratt med over 90 % av den totale mengden DNA (striplete linien). **B)** sannsynligheten for å kunne koble et artsnavn til en operasjonell 'art' for fem utvalgte ordener. For biller (Coleoptera) kan vi sette et artsnavn på over 75 % av alle operasjonelle 'arter', mens andelen synker bratt for tovinger og veps (Diptera, Hymenoptera) jo mindre vanlig 'arten' er."
#| #| fig-format: png

Afig <- rank_abund_long %>%
  # filter(rank_abundance <= 13200) %>%
  ggplot(., aes(x = rank_abundance, y = value, fill = taxonomic_level)) +
  geom_bar(stat = "identity", position = "fill", width = 100) +
  ylab("Andel arter identifisert") +
  xlab("Arter i rekkefølge av mengde DNA") +
  scale_fill_nina() +
    scale_x_continuous(
    breaks = c(1, 5000, 10000, 15000, 20000, max(rank_abund_long$rank_abundance)),
    labels = c("Mest\nDNA", 5000, 10000, 15000, 20000, "Minst\nDNA")
    ) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = rel(0.6)),
    axis.text.x = element_text(size = rel(0.8)),
    axis.title.x = element_text(size = rel(0.8))
  ) +
  guides(fill = guide_legend(override.aes = list(size = 5), nrow = 2, byrow = T)) +
  geom_vline(xintercept = min(rank_abund$row_number[rank_abund$pct_all_sequences > 0.899999999]), linetype = "dashed") +
  ggtitle("A")

# logit_perc <- scales::trans_new("logit perc",
# transform = function(x)qlogis(x/100),
# inverse = function(x)100*plogis(x)
# )

Bfig <- ggplot(new_data) +
  geom_line(aes(x = rel_order_rank, y = pred, color = order),
    linewidth = 2.5
  ) +
  theme(
    # axis.text.x = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.background = element_blank(),
    # axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = rel(0.4)), #[Sjekk at all tekst vises]
    axis.title.y = element_text(size = rel(0.7)),
    axis.text.x = element_text(size = rel(0.8)),
    axis.title.x = element_text(size = rel(0.8))
  ) +
  guides(color = guide_legend(override.aes = list(size = 5), nrow = 2, byrow = T)) +
  ylab("Sannsynlighet for identifikasjon til artsnivå") +
  scale_color_nina() +
  xlab("Arter i rekkefølge av mengde DNA") +
  scale_x_reverse(
    breaks = c(1, 0.75, 0.5, 0.25, 0),
    labels = c("Mest\nDNA", "3/4", "1/2", "1/4", "Minst\nDNA")
  ) +
  ggtitle("B")

#+
# xlim(c(0, 6000))
# +  scale_y_continuous(trans = "logit")

grid.arrange(Afig, Bfig, ncol = 2)


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-barcode-map
#| fig-cap: Lokaliteter for insektinnsamling for strekkoding til referansebibliotek.
#| fig-pos: '!h'


barcodes_raw <- readxl::read_excel("../../rawData/strekkoding_lok/Sampling points NorIns reference library 2024.xlsx",
                                   skip = 1)

barcodes <- barcodes_raw |> 
  mutate(Lat = ifelse(!is.na(`Exact Site`) & (`Exact Site` == 'Seiloy' | `Exact Site` == 'SeilOy' ), 59.12496, Lat)) |> 
  mutate(Lat = ifelse(!is.na(`Exact Site`) & (`Exact Site` == 'Rosnesbukta'), 59.34515, Lat)) |> 
  mutate(Lon = ifelse(!is.na(Sector) & (Sector == 'Lovoya'), 10.44156, Lon)) |> 
  mutate(Lat = ifelse(!is.na(`Exact Site`) & (`Exact Site` == 'Borselvnes' & Elev == '10'), 70.3163, Lat))
  
barcodes <- barcodes |> 
  filter(`Country/Ocean` != 'Sweden')

barcodes <- sf::st_as_sf(barcodes,
                     coords = c("Lon", "Lat"),
                     crs = 4326
                     ) |> 
  st_transform(crs = 25833)

ggplot(norge) +
  geom_sf() +
  geom_sf(aes(),
          color = nina_colors["blue"],
          data = barcodes) +
  scale_color_nina(name = "Fylke") +
  theme(legend.position = "")



## ------------------------------------------------------------------------------------------------------------------------------------
#| label: get_barcode_data
#| cache: !expr cache_var
sk_figdat = read.table("../../../../store_filer_for_rapport/strekkode_db_coverage_data.txt", header = T, sep = "\t")


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: mk_added_barcode_fig
#| cache: !expr cache_var
strekkod_fig <- ggplot(sk_figdat,
                       aes(x = year, fill = in_NI2.0)) + 
  geom_bar(stat = "count", position = "stack", colour = "black") + 
  scale_fill_manual(values = c("white", "#004F71"), 
                    name = "Finnes i\nNorInvert 2.0", 
                    labels = c("Nei","Ja")) + 
  theme_classic() + 
  ylab("Antall arter") + 
  theme(axis.title.x = element_blank())


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-added-barcodes
#| fig-cap: Antall taksa strekkodet gjennom prosjektet mellom 2019-2023 og andelen av disse som er inkludert i sekvensdatabasen (NorInvert 2.0) og dermed brukes til artsbestemmelse. Arter strekkodet i 2024 vil bli inkludert i neste runde.
#| cache: !expr cache_var
strekkod_fig


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-physiphora_alceae
#| fig-pos: '!h'
#| out-height: 5in
#| out-width: 4in
#| fig-show: hold
#| fig-cap: "Fluen Physiphora alceae ble i 2024 verifisert av Jørn Gustad, som også har funnet arten uavhengig. Arten er ikke tidligere kjent fra Norge og er også første representant for sin slekt i Norge. Foto: Arnstein Staverløkk."

include_graphics("figures/Physiphora_alceae_small.jpg")








## ------------------------------------------------------------------------------------------------------------------------------------
#| label: load_hym_family_comp
#| cache: !expr cache_var
load(file = "out/hym_family_comp.Rdata")

hym_family_comp <- hym_family_comp |> 
 mutate(type = factor(type,
                      levels = c("Metastrekkoding",
                                 "Manuell")
                      )
        )


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-ver_hym_family
#| fig-cap: Sammenligning av antall observerte arter veps (Hymenoptera) for ulike identifiseringsmetoder, fordelt på familie. Noter at flere arter er manuelt bestemt bare til morfoart, og identifisert via strekkoding som unike arter selv om de ikke er blitt navngitt. 
#| fig-pos: "!h"
#| fig-height: 7

ggplot(hym_family_comp) +
  geom_bar(aes(x = richness, 
               y = id_family,
               fill = type),
           position = "dodge",
           stat = "identity") +
  facet_wrap(facets = "sample_name_short",
             ncol = 1) +
  coord_flip() +
  scale_fill_nina(name = "Identifiseringstype",
                  palette = "contrast") +
  xlab("Antall arter/morfoarter") +
  ylab("Familie") +
  theme(axis.text.x = element_text(size = 5, angle = 90, hjust = 0.5))




## ------------------------------------------------------------------------------------------------------------------------------------
avg_beetles_met <- beetles_methods_richn |> 
  group_by(identification_name) |> 
  summarise(tot_richn = mean(no_species),
            sd_tot_rich = sd(no_species))

mult <- qnorm(1 - (1 - 0.95) / 2)

avg_beetles_met <- avg_beetles_met |> 
    mutate(
    ci_lower = tot_richn - mult * sd_tot_rich,
    ci_upper = tot_richn + mult * sd_tot_rich
  )


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-ident_coleopt_eff_plot
#| fig-pos: "!h"
#| fig-cap: Estimert antall arter biller (Coleoptera) per prøve fordelt på identifiseringsmetodikk. Feilmarginene er basert på modellert effekt av identifiseringsmetodikkene, en glm med felleprøve og radnummer som tilfeldig effekt, og 'Poisson' som fordeling på modellfeilene.

beetle_eff_plot +
  theme_few()
   


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: get_hym_crushed
#| cache: !expr cache_var

hym_metabar_methods_richn_q <- "
SELECT st.sample_name,  i.identification_name, count(distinct species_latin_fixed) as no_species
FROM occurrences.observations o,
events.identifications i,
events.sampling_trap st,
(SELECT distinct i.sampling_trap_id
FROM occurrences.observations o,
events.identifications i 
WHERE o.identification_id = i.id
AND i.identification_name = 'crushed_novaseq') foo
WHERE i.sampling_trap_id = foo.sampling_trap_id
AND o.identification_id = i.id
AND st.id = foo.sampling_trap_id
--AND i.identification_name = 'metabarcoding_novaseq'
AND o.id_order = 'Hymenoptera'
GROUP BY st.sample_name, i.sampling_trap_id, i.identification_name
--AND i.sampling_trap_id = '1945f610-44fc-48ab-8b52-56b7fcc7c7ad'
ORDER BY st.sample_name, i.identification_name
"
hym_metabar_crushed_richn <- dbGetQuery(con,
                                        hym_metabar_methods_richn_q)  |> 
  as_tibble()



## ------------------------------------------------------------------------------------------------------------------------------------
#| label: arrange_hym_crushed
#| cache: !expr cache_var

hym_methods_richn <- hym_metabar_crushed_richn |> 
  mutate(no_species = as.integer(no_species)) |> 
  mutate(sample_name_short = gsub("(Ost_2024_)(.*)(_VF4)$", "\\2", sample_name)) |> 
  mutate(identification_name = case_when(
    identification_name == 'metabarcoding_novaseq' ~ 'Lysert metastrekkodet',
    identification_name == 'crushed_novaseq' ~ 'Knust metastrekkodet',
    .default = identification_name)
    ) |>
  mutate(identification_name = factor(identification_name, 
                                      levels = c("Knust metastrekkodet",
                                                 "Manuell",
                                                 "Lysert metastrekkodet")
                                      )
         ) |> 
  arrange(sample_name, 
          identification_name)

## ------------------------------------------------------------------------------------------------------------------------------------
#| label: mk_hym_crushed_fig
#| cache: !expr cache_var

hym_crushed_fig <- ggplot(hym_methods_richn) +
  geom_bar(aes(y = no_species, 
               x = sample_name_short,
               fill = identification_name),
           position = "dodge",
           stat = "identity") +
  #coord_flip() +
  scale_fill_nina(name = "Identifiseringstype",
                  palette = "contrast") +
  ylab("Vindusfelle") +
  xlab("Antall arter") +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5))



## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-hym_crushed_fig
#| fig-cap: Sammenligning mellom antall arter veps (Hymenoptera) i lyserte respektive knuste vindusfelleprøver (etter lysering). Helt blå stabler betyr at ikke noen art ble funnet etter kun lysering.

hym_crushed_fig 


## ------------------------------------------------------------------------------------------------------------------------------------
overdisp_fun <- function(model) {
  ## number of variance parameters in 
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}


## ------------------------------------------------------------------------------------------------------------------------------------
hym_methods_richn <- hym_methods_richn |> 
  as_tibble() |> 
  mutate(rownumber = as.character(row_number()),
         identification_name = factor(identification_name),
         no_species = as.integer(no_species))

hym_methods_richn_mod <- glmer(no_species ~ identification_name + (1 | sample_name_short), data = hym_methods_richn, family = "poisson")
#summary(hym_methods_richn_mod)
#overdisp_fun(hym_methods_richn_mod)

hym_pred <- ggeffects::predict_response(hym_methods_richn_mod, 
                                        terms = "identification_name",
                                        bias_correction = TRUE)

hym_eff_plot <- plot(hym_pred,
                     show_title = FALSE,
                     colors = "#004F71") +
  ylab("Antall arter") +
  xlab("Identifiseringstype")




## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-ident_hym_eff_plot
#| fig-pos: "!h"
#| fig-cap: Estimert antall arter vepser (Hymenoptera) per prøve fordelt på identifiseringsmetodikk. Feilmarginene er basert på modellert effekt av identifiseringsmetodikkene.

hym_eff_plot +
  theme_few()


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-skog-ost-biomass
#| fig-cap: "Variasjon i biomasse over tid for skogslokaliteter på Østlandet."
#| fig-height: 4
#| fig-pos: '!h'

xlims <- biomass_mf_locality_sampling_time %>%
  filter(
    habitat_type == "Skog",
    region_name == "Østlandet"
  ) %>%
  select(julian_day) %>%
  pull()

xlims <- c(
  min(xlims - 1, na.rm = T),
  max(xlims, na.rm = T)
)

ylims <- biomass_mf_locality_sampling_time %>%
  filter(
    habitat_type == "Skog",
    region_name == "Østlandet",
    year == 2024
  ) %>%
  mutate(amount = sum_wet_weight / no_trap_days) %>%
  pull()

ylims <- c(
  min(ylims, na.rm = T) * 1.1,
  max(ylims, na.rm = T) * 1.1
)

p4 <- biomass_mf_locality_sampling_time %>%
  filter(
    habitat_type == "Skog",
    region_name == "Østlandet",
    year == 2024
  ) %>%
  ggplot(.) +
  geom_line(
    aes(
      y = sum_wet_weight / no_trap_days,
      x = as.Date(end_date),
      color = locality
    ),
    linewidth = 0.75
  ) +
  geom_point(
    aes(
      y = sum_wet_weight / no_trap_days,
      x = as.Date(end_date),
      color = locality
    ),
    size = 2.5
  ) +
  scale_x_date(
    date_labels = "%d.%m.%y",
    limits = as.Date(xlims, origin = as.Date("2024-01-01"))
  ) +
  scale_y_continuous(limits = ylims) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ylab("Biomasse per\nfelledøgn (g/døgn)") +
  xlab("Tømmedato") +
  scale_color_nina(name = "Lokalitet") +
  theme(
    legend.position = "bottom",
    # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 5),
    legend.key.width = unit(0.5, "cm")
  )


p4


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-semi-nat-sor-biomass
#| fig-pos: '!b'
#| fig-cap: "Variasjon i biomasse over tid for semi-naturlig mark i Sørlandet."
#| fig-height: 4

xlims <- biomass_mf_locality_sampling_time %>%
  filter(
    habitat_type == "Semi-nat",
    region_name == "Sørlandet"
  ) %>%
  select(julian_day) %>%
  pull()

xlims <- c(
  min(xlims - 1, na.rm = T), # -1 because julian days seem to start at 0 in R...
  max(xlims, na.rm = T)
)

ylims <- biomass_mf_locality_sampling_time %>%
  filter(
    habitat_type == "Skog",
    region_name == "Østlandet",
    year == 2024
  ) %>%
  mutate(amount = sum_wet_weight / no_trap_days) %>%
  pull()

ylims <- c(
  min(ylims, na.rm = T) * 1.1,
  max(ylims, na.rm = T) * 1.1
)

p2 <- biomass_mf_locality_sampling_time %>%
  filter(
    habitat_type == "Semi-nat",
    region_name == "Sørlandet",
    year == 2024
  ) %>%
  ggplot(.) +
  geom_line(
    aes(
      y = sum_wet_weight / no_trap_days,
      x = as.Date(end_date),
      color = locality
    ),
    linewidth = 0.75
  ) +
  geom_point(
    aes(
      y = sum_wet_weight / no_trap_days,
      x = as.Date(end_date),
      color = locality
    ),
    size = 2.5
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ylab("Biomasse per felledøgn (g/døgn)") +
  scale_x_date(
    date_labels = "%d.%m.%y",
    limits = as.Date(xlims, origin = as.Date("2024-01-01"))
  ) +
  # scale_y_continuous(limits = ylims) +
  xlab("Tømmedato") +
  scale_color_nina(name = "Lokalitet") +
  theme(
    legend.position = "bottom",
    # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 5),
    legend.key.width = unit(0.5, "cm")
  )


p2


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-semi-nat-vest-biomass
#| fig-pos: '!b'
#| fig-cap: "Variasjon i biomasse over tid for semi-naturlig mark på Vestlandet"
#| fig-height: 4

xlims <- biomass_mf_locality_sampling_time %>%
  filter(
    habitat_type == "Semi-nat",
    region_name == "Vestlandet"
  ) %>%
  select(julian_day) %>%
  pull()

xlims <- c(
  min(xlims - 1, na.rm = T),
  max(xlims, na.rm = T)
)

ylims <- biomass_mf_locality_sampling_time %>%
  filter(
    habitat_type == "Semi-nat",
    region_name == "Vestlandet",
    year == 2024
  ) %>%
  mutate(amount = sum_wet_weight / no_trap_days) %>%
  pull()

ylims <- c(
  min(ylims, na.rm = T) * 1.1,
  max(ylims, na.rm = T) * 1.1
)

p4 <- biomass_mf_locality_sampling_time %>%
  filter(
    habitat_type == "Semi-nat",
    region_name == "Vestlandet",
    year == 2024
  ) %>%
  ggplot(.) +
  geom_line(
    aes(
      y = sum_wet_weight / no_trap_days,
      x = as.Date(end_date),
      color = locality
    ),
    linewidth = 0.75
  ) +
  geom_point(
    aes(
      y = sum_wet_weight / no_trap_days,
      x = as.Date(end_date),
      color = locality
    ),
    size = 2.5
  ) +
  scale_x_date(
    date_labels = "%d.%m.%y",
    limits = as.Date(xlims, origin = as.Date("2024-01-01"))
  ) +
  scale_y_continuous(limits = ylims) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ylab("Biomasse per\nfelledøgn (g/døgn)") +
  xlab("Tømmedato") +
  scale_color_nina(name = "Lokalitet") +
  theme(
    legend.position = "bottom",
    # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 5),
    legend.key.width = unit(0.5, "cm")
  )


p4


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-semi-nat-ost-biomass
#| fig-pos: '!b'
#| fig-cap: "Variasjon i biomasse over tid for semi-naturlig mark på Østlandet."
#| fig-height: 4

xlims <- biomass_mf_locality_sampling_time %>%
  filter(
    habitat_type == "Semi-nat",
    region_name == "Østlandet"
  ) %>%
  select(julian_day) %>%
  pull()

xlims <- c(
  min(xlims - 1, na.rm = T),
  max(xlims, na.rm = T)
)

ylims <- biomass_mf_locality_sampling_time %>%
  filter(
    habitat_type == "Semi-nat",
    region_name == "Østlandet",
    year == 2024
  ) %>%
  mutate(amount = sum_wet_weight / no_trap_days) %>%
  pull()

ylims <- c(
  min(ylims, na.rm = T) * 1.1,
  max(ylims, na.rm = T) * 1.1
)

p4 <- biomass_mf_locality_sampling_time %>%
  filter(
    habitat_type == "Semi-nat",
    region_name == "Østlandet",
    year == 2024
  ) %>%
  ggplot(.) +
  geom_line(
    aes(
      y = sum_wet_weight / no_trap_days,
      x = as.Date(end_date),
      color = locality
    ),
    linewidth = 0.75
  ) +
  geom_point(
    aes(
      y = sum_wet_weight / no_trap_days,
      x = as.Date(end_date),
      color = locality
    ),
    size = 2.5
  ) +
  scale_x_date(
    date_labels = "%d.%m.%y",
    limits = as.Date(xlims, origin = as.Date("2024-01-01"))
  ) +
  scale_y_continuous(limits = ylims) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ylab("Biomasse per\nfelledøgn (g/døgn)") +
  xlab("Tømmedato") +
  scale_color_nina(name = "Lokalitet") +
  theme(
    legend.position = "bottom",
    # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 5),
    legend.key.width = unit(0.5, "cm")
  )


p4


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-semi-nat-trond-biomass
#| fig-cap: "Variasjon i biomasse over tid for semi-naturlig mark i Trøndelag."
#| fig-height: 4
xlims <- biomass_mf_locality_sampling_time %>%
  filter(
    habitat_type == "Semi-nat",
    region_name == "Trøndelag",
    year == 2024
  ) %>%
  select(julian_day) %>%
  pull()

xlims <- c(
  min(xlims - 1, na.rm = T),
  max(xlims, na.rm = T)
)

ylims <- biomass_mf_locality_sampling_time %>%
  filter(
    habitat_type == "Semi-nat",
    region_name == "Trøndelag",
    year == 2024
  ) %>%
  mutate(amount = sum_wet_weight / no_trap_days) %>%
  pull()

ylims <- c(
  min(ylims, na.rm = T) * 1.1,
  max(ylims, na.rm = T) * 1.1
)

p3 <- biomass_mf_locality_sampling_time %>%
  filter(
    habitat_type == "Semi-nat",
    region_name == "Trøndelag",
    year == 2024
  ) %>%
  ggplot(.) +
  geom_line(
    aes(
      y = sum_wet_weight / no_trap_days,
      x = as.Date(end_date),
      color = locality
    ),
    linewidth = 0.75
  ) +
  geom_point(
    aes(
      y = sum_wet_weight / no_trap_days,
      x = as.Date(end_date),
      color = locality
    ),
    size = 2.5
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ylab("Biomasse per felledøgn (g/døgn)") +
  scale_x_date(
    date_labels = "%d.%m.%y",
    limits = as.Date(xlims, origin = as.Date("2024-01-01"))
  ) +
  scale_y_continuous(limits = ylims) +
  xlab("Tømmedato") +
  scale_color_nina(name = "Lokalitet") +
  theme(
    legend.position = "bottom",
    # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 5),
    legend.key.width = unit(0.5, "cm")
  )


p3


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-semi-nat-nord-biomass
#| fig-cap: "Variasjon i biomasse over tid for semi-naturlig mark i Nord-Norge."
#| fig-height: 4
xlims <- biomass_mf_locality_sampling_time %>%
  filter(
    habitat_type == "Semi-nat",
    region_name == "Nord-Norge"
  ) %>%
  select(julian_day) %>%
  pull()

xlims <- c(
  min(xlims - 1, na.rm = T),
  max(xlims, na.rm = T)
)

#Hard set the mi
if(xlims[1] > 151){xlims[1] <- 151}

ylims <- biomass_mf_locality_sampling_time %>%
  filter(
    habitat_type == "Semi-nat",
    region_name == "Nord-Norge",
    year == 2024
  ) %>%
  mutate(amount = sum_wet_weight / no_trap_days) %>%
  pull()

ylims <- c(
  min(ylims, na.rm = T) * 1.1,
  max(ylims, na.rm = T) * 1.1
)

p1 <- biomass_mf_locality_sampling_time %>%
  filter(
    habitat_type == "Semi-nat",
    region_name == "Nord-Norge",
    year == 2024
  ) %>%
  ggplot(.) +
  geom_line(
    aes(
      y = sum_wet_weight / no_trap_days,
      x = as.Date(end_date),
      color = locality
    ),
    linewidth = 0.75
  ) +
  geom_point(
    aes(
      y = sum_wet_weight / no_trap_days,
      x = as.Date(end_date),
      color = locality
    ),
    size = 2.5
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ylab("Biomasse per felledøgn (g/døgn)") +
  scale_x_date(
    date_labels = "%d.%m.%y",
    limits = as.Date(xlims, origin = as.Date("2024-01-01"))
  ) +
  scale_y_continuous(limits = ylims) +
  xlab("Tømmedato") +
  scale_color_nina(name = "Lokalitet") +
  theme(
    legend.position = "bottom",
    # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 5),
    legend.key.width = unit(0.5, "cm")
  )


p1


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: get_order_phen
#| cache: !expr cache_var 
øst_order_phen <- get_phenology(taxonomic_level = "Order",
                                subset_region = "Østlandet")


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: get_phen_temps
#| cache: !expr cache_var

locs <- øst_order_phen |> 
  select(locality) |>
  unique()

clim_tbl <- dplyr::tbl(con, DBI::Id(schema = "climate_data",
                                    table = "se_norge"))

#Create cumulative sum, not done
loc_clim <- clim_tbl |> 
  filter(date < '2025-01-01',
         date > '2019-12-31',
         locality %in% locs$locality) |> 
  mutate(locality_tg = ifelse(locality_tg >=0, locality_tg, 0)) |>  #set negative temps to 0
  arrange(locality, date) |> 
  collect() |> 
  mutate(year = lubridate::year(date),
         jday = lubridate::yday(date)) |> 
  group_by(locality, year) |> 
  mutate(csum = cumsum(locality_tg))



## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: !expr cache_var
loc_clim_strip <- loc_clim |> 
  ungroup() |> 
  select(locality,
         date,
         csum)

øst_order_phen_temp <- øst_order_phen |> 
  left_join(loc_clim_strip,
            by = c("end_date_obs" = "date",
                   "locality" = "locality")
  )



## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: !expr cache_var
#| label: fig-ost_order_phen_dato_tempsum
#| fig-pos: '!h'
#| fig-height: 8
#| fig-width: 7
#| fig-cap: Fenologi for de største insektordener i forhold til innsamlingsdato (juliansk dato) og temperatursum for årene 2000-20024. Y-aksen viser logaritmert antall arter per flasketømming gjennom sesongen

pp1 <- plot(øst_order_phen_temp,
            type = "no_species",
            x_axis_type = "date",
            scale_to_max = FALSE) +
    theme(legend.position = "bottom",
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 7),
    legend.key.width = unit(0.5, "cm")
    ) +
  guides(fill = guide_legend(
    nrow = 2,
    byrow = TRUE
  ))

pp2 <- plot(øst_order_phen_temp,
            type = "no_species",
            x_axis_type = "temperature_sum",
            scale_to_max = FALSE)

leg <- cowplot::as_grob(ggpubr::get_legend(pp1))

pp1 <- pp1 + theme_bw() + 
  theme(legend.position = "none",
        plot.margin = unit(c(2, 5.5, 5.5, 0.1), "pt")
        ) + 
  ylab("Antall arter (logaritmert)")

pp2 <- pp2 + 
  theme_bw() + 
  theme(legend.position = "none",
        plot.margin = unit(c(2, 5.5, 5.5, 0.1), "pt")
        ) + 
  ylab("")

# grid.arrange(pp1,pp2,leg,
#              layout_matrix = rbind(c(1,1,1, 2,2,2, 3) )  )

lm <- matrix(c(1, 1, 1, 3, 2, 2, 2, 3),
  nrow = 4,
  byrow = F
)

gridExtra::marrangeGrob(
  list(
    pp1,
    pp2,
    leg
  ),
  top = "",
  layout_matrix = lm,
  padding = unit(0, "line")
)





## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: !expr cache_var
load(file = "out/gllvm_div_biom.Rdata")






## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: !expr cache_var
load(file = "out/fit_multi_div.Rdata")




## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: !expr cache_var
load(file = "out/fit_multi_div_csum.Rdata")



## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: !expr cache_var
#| label: fig-phenology_div
#| fig-cap: Modellert aktivitetstopp i diversitet for de 5 største insektordnene i overvåkingsprogrammet. Y-aksen viser  logaritmert artsantall for hver orden per prøvetilfelle.

phen_peak_plot(fit_multi_div) +
  scale_color_nina(name = "Orden") +
  ylab("Log. artsantall") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-biomass-year-hab
#| fig-pos: '!h'
#| fig-cap: "Gjennomsnitt og 95\\% konfidensintervall for insektbiomasse  i malaisefeller per dag og lokalitet, fra juli-august (+- k * se der k er 97.5 persentilen i en t-studentfordeling). Merk at populasjonstopper kan skille seg noe i tid mellom år, og at disse måneder ikke representerer hele sesongen."

to_plot <- biomass_mf_locality_sampling_time %>%
  # filter(habitat_type == 'Semi-nat') %>%
  mutate(habitat_no = ifelse(habitat_type == "Semi-nat", "Gressmark", "Skog")) %>%
  mutate(habitat_no = factor(habitat_no, levels = c("Gressmark", "Skog"))) %>%
  mutate(year = as.factor(year)) %>%
  filter(
    start_month >= 7,
    start_month <= 8
  ) %>%
  group_by(
    year,
    region_name,
    habitat_type
  ) %>%
  summarise(
    mean_biomass = mean(avg_wet_weight / no_trap_days, na.rm = TRUE),
    sd_biomass = sd(avg_wet_weight / no_trap_days, na.rm = TRUE),
    mean_log_bm_per_day = mean(exp(log_bm_per_day), na.rm = TRUE),
    sd_log_bm_per_day = sd(exp(log_bm_per_day), na.rm = TRUE),
    n = n(),
    se_biomass = sd_biomass / sqrt(n),
    se_log_bm_per_day = sd_log_bm_per_day / sqrt(n)
  ) %>%
  mutate(
    year = as.factor(year),
    region_name = as.factor(region_name),
    habitat_type = as.factor(habitat_type)
  ) %>%
  ungroup()

region_cols <- c("Nord-Norge" =  "#004F71",
                 "Trøndelag" =  "#93328E",
                 "Sørlandet" =  "#E57200",
                 "Østlandet" =   "#008C95",
                 "Vestlandet" =  "#7A9A01") 

ciMult <- qt(.975, to_plot$n - 1)
dodge <- 0.4

ggplot(
  to_plot,
  aes(
    x = year,
    y = mean_biomass,
    group = region_name:habitat_type
  )
) +
  geom_point(
    aes(
      x = year,
      y = mean_biomass,
      color = region_name,
      pch = habitat_type
    ),
    cex = 4,
    position = position_dodge(width = dodge)
  ) +
  geom_errorbar(
    aes(
      x = year,
      ymin = mean_biomass - ciMult * se_biomass,
      ymax = mean_biomass + ciMult * se_biomass,
      color = region_name,
      # group = habitat_type
    ),
    width = .2,
    linewidth = 1,
    position = position_dodge(width = dodge)
  ) +
  geom_line(
    aes(
      x = year,
      y = mean_biomass,
      # group = region_name:habitat_type
    ),
    position = position_dodge(width = dodge)
  ) +
  scale_color_manual(
    name = "Region",
    values = region_cols
  ) +
  scale_shape_discrete(name = "Habitattype") +
  ylab("Middelverdi biomasse per dag og lokalitet") +
  xlab("År") +
  theme(
    legend.position = "right",
    axis.title.y = element_text(size = 10)
  )
# guides(color=guide_legend(nrow=2 , byrow=TRUE))
# ylim(c(0, max(ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range))) +
# ggtitle("B")



## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: !expr cache_var
pollinator_table <- tibble(
  "Norsk navn" = c(
    "Gravebier",
    "Langtungebier",
    "Korttungebier",
    "Markbier",
    "Buksamlerbier",
    "Blomsterbier",
    "Blomfluer",
    "Smygere",
    "Glansvinger",
    "Nymfevinger",
    "Svalestjerter",
    "Hvitvinger",
    "Metallmerker (uoffisiell)"
  ),
  Familie = c(
    "Andrenidae",
    "Apidae",
    "Colletidae",
    "Halictidae",
    "Megachilidae",
    "Melittidae",
    "Syrphidae",
    "Hesperiidae",
    "Lycaenidae",
    "Nymphalidae",
    "Papilionidae",
    "Pieridae",
    "Riodinidae"
  )
)

pollinators <- pollinator_table %>%
  select(Familie) %>%
  pull()


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: !expr cache_var
year_locality_pollinators <- get_observations(
  subset_families = pollinators,
  subset_region = "Østlandet",
  agg_level = "year_locality",
  trap_type = "MF"
) %>%
  mutate(habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type))


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-div-year-hab
#| fig-pos: '!h'
#| fig-cap: "Gjennomsnitt av antall arter fanget per lokalitet og år i skog og semi-naturlig mark. Merk at i 2021 ble novaseq-maskinen tatt i bruk for metastrekkoding, som påviser flere arter grunnet større sekvenseringsdybde. Figuren viser data kun fra malaisefeller. Usikkerhetsstolper viser 95 % konfidensintervall. Kun data fra Østlandet er tatt med for sammenlignbare data."

year_mean_div <- year_locality_mf_ost_div %>%
  group_by(
    year,
    habitat_type
  ) %>%
  summarise(
    mean_richness = mean(no_species),
    sd_richness = sd(no_species),
    n = n(),
    se_richness = sd_richness / sqrt(n)
  ) %>%
  mutate(
    year = as.factor(year),
    habitat_type = as.factor(habitat_type)
  ) %>%
  ungroup()

dodge <- 0.2

ciMult <- qt(.975, year_mean_div$n - 1)

p1 <- ggplot(
  year_mean_div,
  aes(
    x = year,
    y = mean_richness,
    group = habitat_type
  )
) +
  geom_point(
    aes(
      x = year,
      y = mean_richness,
      color = habitat_type
    ),
    cex = 3,
    position = position_dodge(dodge)
  ) +
  geom_errorbar(
    aes(
      x = year,
      ymin = mean_richness - ciMult * se_richness,
      ymax = mean_richness + ciMult * se_richness,
      color = habitat_type
    ),
    width = .2,
    linewidth = 1,
    position = position_dodge(dodge)
  ) +
  geom_line(
    aes(
      x = year,
      y = mean_richness,
      group = habitat_type
    ),
    position = position_dodge(dodge)
  ) +
  scale_color_nina(
    name = "Habitatstype",
    palette = "darkblue-orange"
  ) +
  ylab("Antall arter per lokalitet") +
  xlab("År") +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) 

p1


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-div-year-hab-poll
#| fig-pos: '!h'
#| fig-cap: "Gjennomsnitt av antall pollinatorer fanget per lokalitet og år i skog og semi-naturlig mark (alle arter i 13 insektfamilier kjente som pollinatorer). Merk at i 2021 ble novaseq-maskinen tatt i bruk for metastrekkoding, som påviser flere arter grunnet større sekvenseringsdybde. Figuren viser data kun fra malaisefeller. Usikkerhetsstolper viser 95 % konfidensintervall. Kun data fra Østlandet er tatt med for sammenlignbare data."
year_mean_pollinators <- year_locality_pollinators %>%
  group_by(
    year,
    habitat_type
  ) %>%
  summarise(
    mean_richness = mean(no_species),
    sd_richness = sd(no_species),
    n = n(),
    se_richness = sd_richness / sqrt(n)
  ) %>%
  mutate(
    year = as.factor(year),
    habitat_type = as.factor(habitat_type)
  ) %>%
  ungroup()

p2 <- ggplot(
  year_mean_pollinators,
  aes(
    x = year,
    y = mean_richness,
    group = habitat_type
  )
) +
  geom_point(
    aes(
      x = year,
      y = mean_richness,
      color = habitat_type
    ),
    cex = 3,
    position = position_dodge(dodge)
  ) +
  geom_errorbar(
    aes(
      x = year,
      ymin = mean_richness - ciMult * se_richness,
      ymax = mean_richness + ciMult * se_richness,
      color = habitat_type
    ),
    width = .2,
    linewidth = 1,
    position = position_dodge(dodge)
  ) +
  geom_line(
    aes(
      x = year,
      y = mean_richness,
      group = habitat_type
    ),
    position = position_dodge(dodge)
  ) +
  scale_color_nina(
    name = "Habitatstype",
    palette = "darkblue-orange"
  ) +
  ylab("Antall arter pollinatorer per lokalitet") +
  xlab("År") +
  theme(
    legend.position = "bottom",
    axis.title.y = element_text(size = 8)
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) 

p2


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fit_year_model

year_model <- lmer(
  log_bm_per_day ~
    region_name +
    habitat_type +
    year_num +
    # (1|year_num) +
    (1 | Ukenummer) +
    (1 | locality),
  data = subset(
    locality_sampling_biomass_temp_compl_light_week,
    # habitat_type == "Semi-nat" &
    # region_name == "Østlandet"
  ),
  na.action = "na.fail",
  REML = F
)


## ------------------------------------------------------------------------------------------------------------------------------------
#| eval: false
# summary(year_model)


## ------------------------------------------------------------------------------------------------------------------------------------
regional_year_model <- lmer(
  log_bm_per_day ~
    habitat_type +
    region_name*year_num +
    (1 | Ukenummer) +
    (1 | locality),
  data = subset(
    locality_sampling_biomass_temp_compl_light_week,
    # habitat_type == "Semi-nat" &
    # region_name == "Østlandet"
  ),
  na.action = "na.fail",
  REML = F
)


## ------------------------------------------------------------------------------------------------------------------------------------
#| eval: false
# summary(regional_year_model)


## ------------------------------------------------------------------------------------------------------------------------------------
year_model_rand_year <- lmer(
  log_bm_per_day ~
    region_name +
    habitat_type +
    year_num +
    (1 | year_num) +
    (1 | Ukenummer) +
    (1 | locality),
  data = subset(
    locality_sampling_biomass_temp_compl_light_week,
    # habitat_type == "Semi-nat" &
    # region_name == "Østlandet"
  ),
  na.action = "na.fail",
  REML = F
)


## ------------------------------------------------------------------------------------------------------------------------------------
tt <- summary(year_model)
year_to_tab <- tt$coefficients %>%
  as_tibble(rownames = "Parameter") %>%
  mutate(Estimate.norm = exp(Estimate))

year_to_tab <- year_to_tab %>%
  mutate(Parameter = case_match(
    Parameter,
    "(Intercept)" ~ "Intercept",
    "year_num" ~ "Årstrend",
    "region_nameNord-Norge" ~ "Region: Nord-Norge",
    "region_nameSørlandet" ~ "Region: Sørlandet",
    "region_nameTrøndelag" ~ "Region: Trøndelag",
    "region_nameVestlandet" ~ "Region: Vestlandet",
    "habitat_typeSemi-nat" ~ "Habitat: Semi-nat"
    
  )) %>%
  select(
    Parameter,
    Estimate.norm,
    Estimate,
    `Std. Error`,
    # df,
    `t value` # ,
    # `Pr(>|t|)`
  )


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: tbl-year-mod
#| results: asis
#| tbl-cap: "Parameterestimat fra modell av årstrend i biomasse. Parameterestimatene vises på tilbaketransformert (Estimate.norm) og lineær skala."

if (names(rmarkdown::metadata$format) == "ninareport-pdf") {
  xtable(year_to_tab, digits = 3)
} else {
  kable(year_to_tab, digits = 3)
}


## ------------------------------------------------------------------------------------------------------------------------------------
tt <- summary(regional_year_model)
regional_year_to_tab <- tt$coefficients %>%
  as_tibble(rownames = "Parameter") %>%
  mutate(Estimate.norm = exp(Estimate))

regional_year_to_tab <- regional_year_to_tab %>%
  mutate(Parameter = case_match(
    Parameter,
    "(Intercept)" ~ "Intercept",
    "year_num" ~ "Årstrend Østlandet",
    "region_nameNord-Norge" ~ "Region: Nord-Norge",
    "region_nameSørlandet" ~ "Region: Sørlandet",
    "region_nameTrøndelag" ~ "Region: Trøndelag",
    "region_nameVestlandet" ~ "Region: Vestlandet",
    "habitat_typeSemi-nat" ~ "Habitat: Semi-nat",
    "region_nameSørlandet:year_num" ~ "Forskjell årstrend Sørlandet",
    "region_nameTrøndelag:year_num" ~ "Forskjell årstrend Trøndelag",
    "region_nameNord-Norge:year_num" ~ "Forskjell årstrend Nord-Norge",
    
    
  )) %>%
  select(
    Parameter,
    Estimate.norm,
    Estimate,
    `Std. Error`,
    # df,
    `t value` # ,
    # `Pr(>|t|)`
  )


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: tbl-regional-year-mod
#| results: asis
#| tbl-cap: "Parameterestimat fra modell av regionale biomasse-årstrender. Parameterestimatene vises på tilbaketransformert (Estimate.norm) og lineær skala. Referansenivået er Østlandet. Derfor vises årstrender for Sørlandet, Trøndelag og Nord-Norge som hvor mye disse trender skiller seg fra trenden på Østlandet. Vestlandet har ikke noen estimert trend da vi her har kun ett år med data."

if (names(rmarkdown::metadata$format) == "ninareport-pdf") {
  xtable(regional_year_to_tab, digits = 3)
} else {
  kable(regional_year_to_tab, digits = 3)
}


## ------------------------------------------------------------------------------------------------------------------------------------
r2_verdi_biomasse_year_model <- r.squaredGLMM(year_model)
r2_verdi_biomasse_year_model_rand_year <- r.squaredGLMM(year_model_rand_year)


## ------------------------------------------------------------------------------------------------------------------------------------
fe <- fixef(year_model)

to_mult <- locality_sampling_biomass_temp_compl_light_week %>%
  select(
    sampling_name,
    region_name,
    year_num,
    habitat_type
  ) %>%
  mutate("TRUE" = TRUE) %>%
  pivot_wider(
    names_from = region_name,
    values_from = "TRUE",
    values_fill = FALSE
  ) %>%
  mutate(
    `(Intercept)` = 1,
    `semi-nat` = ifelse(habitat_type == "Semi-nat", 1, 0)
  ) %>%
  select(
    `(Intercept)`,
    Sørlandet,
    Trøndelag,
    `Nord-Norge`,
    Vestlandet,
    `semi-nat`,
    year_num
  )


prediction <- as.matrix(to_mult) %*% fe

to_cor_test <- locality_sampling_biomass_temp_compl_light_week %>%
  select(
    Year,
    log_bm_per_day
  ) %>%
  cbind(prediction)


## ------------------------------------------------------------------------------------------------------------------------------------
#| eval: true
rho <- cor.test(~ prediction + log_bm_per_day, data = to_cor_test)$estimate



## ------------------------------------------------------------------------------------------------------------------------------------
#| label: plot_model1
#| cache: !expr cache_var

p1 <- plot_model(year_model,
  sort.est = TRUE,
  axis.labels = c(
    "year_num" = "År",
    "habitat_typeSemi-nat" = "Habitat: Semi-nat",
    "region_nameNord-Norge" = "Region: Nord-Norge",
    "region_nameTrøndelag" = "Region-Trøndelag",
    "region_nameSørlandet" = "Region-Sørlandet",
    "region_nameVestlandet" = "Region-Vestlandet"
  ),
  title = ""
) +
  xlab("") +
  ggtitle("A") +
  theme(axis.title.y = element_text(size = 9))
# theme_sjplot()


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: plot_model2
#| cache: !expr cache_var

p2 <- plot_model(regional_year_model,
  sort.est = TRUE,
  axis.labels = c(
    "year_num" = "År",
    "habitat_typeSemi-nat" = "Habitat: Semi-nat",
    "region_nameNord-Norge" = "Region: Nord-Norge",
    "region_nameTrøndelag" = "Region-Trøndelag",
    "region_nameSørlandet" = "Region-Sørlandet",
    "region_nameVestlandet" = "Region-Vestlandet",
    "region_nameSørlandet:year_num" = "Forskjell årstrend Sørlandet",
    "region_nameTrøndelag:year_num" = "Forskjell årstrend Trøndelag",
    "region_nameNord-Norge:year_num" = "Forskjell årstrend Nord-Norge"
  ),
  title = ""
) +
  xlab("") +
  ggtitle("B") +
  theme(axis.text.y = element_text(size = 7)) 



## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-year-mod
#| fig-cap: "Estimert tidstrend og regionale forskjeller i insektbiomasse. (A) viser den estimerte felles tidstrenden for alle regioner (Årstrend), og (B) viser den estimerte årstrenden for Østlandet (Årstrend) og forskjellen i trend for øvrige regioner. Et negativt estimat for år betyr at insektbiomassen minker over tid. Estimatene for habitat og region er i forhold til referansenivåene skog og Østlandet."
#| fig-pos: '!hb'
#| fig-height: 6

# theme_sjplot()

marrangeGrob(
  list(
    p1,
    p2
  ),
  top = "",
  ncol = 1,
  nrow = 2
)


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: vær-mod
global_model <- lmer(
  log_bm_per_day ~
    scale_avg_values_MX2301A_temperature +
    scale_avg_values_MX2301A_rh +
    scale_avg_values_MX2202_light +
    scale_stddev_values_MX2301A_temperature +
    scale_stddev_values_MX2301A_rh +
    scale_stddev_values_MX2202_light +
    scale_avg_values_MX2301A_temperature:scale_avg_values_MX2301A_rh +
    scale_avg_values_MX2301A_temperature:scale_avg_values_MX2202_light +
    scale_avg_values_MX2301A_rh:scale_avg_values_MX2202_light +
    scale_stddev_values_MX2301A_temperature:scale_stddev_values_MX2301A_rh +
    scale_stddev_values_MX2301A_temperature:scale_stddev_values_MX2202_light +
    scale_stddev_values_MX2301A_rh:scale_stddev_values_MX2202_light +
    region_name +
    habitat_type +
    year_num +
    # (1|year_num) +
    (1 | Ukenummer) +
    (1 | locality),
  data = locality_sampling_biomass_temp_compl_light_week,
  na.action = "na.fail",
  REML = F
)

## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-var-mod
#| fig-pos: '!hb'
#| fig-cap: "Parameterestimat fra modellen over klimaloggerdata. Avhengig variabel er logaritmen av biomasse insekter per dag (avrendt våtvekt (g/dag)). Forklaringsvariabler er standardiserte (omgjort til Z-verdier). Feillinjer viser 95 % konfidensintervall, og variabler med estimater som ikke overlapper med null effekt (0) er i fet stil (eller som har en interaksjon som ikke overlapper med 0). Reginale estimat er i forhold til referanse-regionen Østlandet. Kontinuerlige variabler er estimert på tvers av region og habitat."

p <- plot_model(global_model,
  sort.est = TRUE,
  # rm.terms = c("habitat_typeSemi-nat"),
  axis.labels = c(
    "scale_avg_values_MX2301A_temperature" = "<b style='color:black'>Temperatur</b>",
    "scale_avg_values_MX2301A_rh" = "<b style='color:black'>Luftfuktighet</b>",
    "scale_avg_values_MX2202_light" = "<b style='color:black'>Lys</b>",
    "scale_stddev_values_MX2301A_temperature" = "<b style='color:black'>Variasjon_temp.</b>",
    "scale_stddev_values_MX2301A_rh" = "<b style='color:black'>Variasjon_fukt.</b>",
    "scale_stddev_values_MX2202_light" = "<b style='color:black'>Variasjon_lys</b>",
    "year_num" = "<b style='color:black'>År</b>",
    "habitat_typeSemi-nat" = "<b style='color:black'>Habitat_Semi-nat</b>",
    "region_nameNord-Norge" = "Region_Nord-Norge",
    "region_nameTrøndelag" = "<b style='color:black'>Region_Trøndelag</b>",
    "region_nameSørlandet" = "<b style='color:black'>Region_Sørlandet</b>",
    "region_nameVestlandet" = "Region_Vestlandet",
    "scale_avg_values_MX2301A_temperature:scale_avg_values_MX2301A_rh" = "Temp:Fukt.",
    "scale_avg_values_MX2301A_temperature:scale_avg_values_MX2202_light" = "<b style='color:black'>Temp:Lys</b>",
    "scale_avg_values_MX2301A_rh:scale_avg_values_MX2202_light" = "Fukt:Lys",
    "scale_stddev_values_MX2301A_temperature:scale_stddev_values_MX2301A_rh" = "<b style='color:black'>Variasjon_temp:Variasjon_fukt</b>",
    "scale_stddev_values_MX2301A_temperature:scale_stddev_values_MX2202_light" = "<b style='color:black'>Variasjon_temp:Variasjon_lys</b>",
    "scale_stddev_values_MX2301A_rh:scale_stddev_values_MX2202_light" = "Variasjon_fukt:Variasjon_lys"
  ),
  title = ""
) +
  xlab("Forklaringsvariabel") +
  # theme_sjplot() +
  ggplot2::theme(axis.text.y = ggtext::element_markdown())

p


## ------------------------------------------------------------------------------------------------------------------------------------
tt <- summary(global_model)
to_tab <- tt$coefficients %>%
  as_tibble(rownames = "Parameter") %>%
  mutate(Estimate.norm = exp(Estimate))

to_tab <- to_tab %>%
  mutate(Parameter = case_match(
    Parameter,
    "(Intercept)" ~ "Intercept",
    "scale_avg_values_MX2301A_temperature" ~ "Temperatur",
    "scale_avg_values_MX2301A_rh" ~ "Luftfuktighet",
    "scale_avg_values_MX2202_light" ~ "Lys",
    "scale_stddev_values_MX2301A_temperature" ~ "Variasjon_temp.",
    "scale_stddev_values_MX2301A_rh" ~ "Variasjon_fukt",
    "scale_stddev_values_MX2202_light" ~ "Variasjon_lys",
    "habitat_typeSemi-nat" ~ "Habitat_Semi-nat",
    "year_num" ~ "År",
    "region_nameNord-Norge" ~ "Region_Nord-Norge",
    "region_nameSørlandet" ~ "Region_Sørlandet",
    "region_nameTrøndelag" ~ "Region_Trøndelag",
    "region_nameVestlandet" ~ "Region_Vestlandet",
    "habitat_typeSemi-nat" ~ "Habitat_Semi-nat",
    "scale_avg_values_MX2301A_temperature:scale_avg_values_MX2301A_rh" ~ "Temperatur:Fukt",
    "scale_avg_values_MX2301A_temperature:scale_avg_values_MX2202_light" ~ "Temperatur:Lys",
    "scale_avg_values_MX2301A_rh:scale_avg_values_MX2202_light" ~ "Fukt:Lys",
    "scale_stddev_values_MX2301A_temperature:scale_stddev_values_MX2301A_rh" ~ "Variasjon_temp.:Variasjon_fukt",
    "scale_stddev_values_MX2301A_temperature:scale_stddev_values_MX2202_light" ~ "Variasjon_temp.:Variasjon_lys",
    "scale_stddev_values_MX2301A_rh:scale_stddev_values_MX2202_light" ~ "Variasjon_fukt:Variasjon_lys"
  )) %>%
  select(
    Parameter,
    Estimate.norm,
    Estimate,
    `Std. Error`,
    # df,
    `t value` # ,
    # `Pr(>|t|)`
  )


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: tbl-var-mod
#| results: asis
#| tbl-cap: "Parameterestimat fra klimaloggermodell. Parameterestimatene vises på tilbaketransformert (Estimate.norm) og linejær skala."

if (names(rmarkdown::metadata$format) == "ninareport-pdf") {
  xtable(to_tab, digits = 3)
} else {
  kable(to_tab, digits = 3)
}


## ------------------------------------------------------------------------------------------------------------------------------------
# sw(model_averages)
r2_verdi_biomasse_temp_glob <- r.squaredGLMM(global_model)


## ----sub_model-----------------------------------------------------------------------------------------------------------------------
# obs jeg tar bort region også for enklere beregning av prediksjon til figur

sub_model_form <- update.formula(global_model, . ~ . - scale_avg_values_MX2301A_temperature:scale_avg_values_MX2301A_rh
  - scale_avg_values_MX2301A_temperature:scale_avg_values_MX2202_light
  - scale_avg_values_MX2301A_rh:scale_avg_values_MX2202_light
  - scale_stddev_values_MX2301A_temperature:scale_stddev_values_MX2301A_rh
  - scale_stddev_values_MX2301A_rh:scale_stddev_values_MX2202_light
  - region_name)

sub_model <- lmer(sub_model_form,
  data = locality_sampling_biomass_temp_compl_light_week,
  na.action = "na.fail",
  REML = F
)


## ----eval = F------------------------------------------------------------------------------------------------------------------------
# summary(sub_model)


## ------------------------------------------------------------------------------------------------------------------------------------
r2_verdi_biomasse_sub <- r.squaredGLMM(sub_model)



## ------------------------------------------------------------------------------------------------------------------------------------
#| eval: true
#| #| eval: false
rho <- cor.test(~ prediction + log_bm_per_day, data = to_cor_test)$estimate
# cor.test(~ prediction + log_bm_per_day, data = to_cor_test)
# 44% vs 73% R2 vs Spearman-Rho



## ------------------------------------------------------------------------------------------------------------------------------------
weather_summary <- tbl(
  con,
  Id(
    schema = "climate_data",
    table = "weather_summary_raw"
  )
) %>%
  collect() %>%
  mutate(
    year_prev = as.integer(year_prev),
    month_prev = as.integer(month_prev)
  )

weather_summary_scaled <- weather_summary %>%
  mutate(across(where(is.double), scale)) %>%
  mutate(across(where(is.double), as.vector))

# names(weather_summary_scaled) <- names(weather_summary)


## ------------------------------------------------------------------------------------------------------------------------------------
locality_sampling_biomass_temp_compl_light_week_weather_anom <- locality_sampling_biomass_temp_compl_light_week %>%
  left_join(weather_summary_scaled,
    by = c(
      "sampling_name" = "sampling_name",
      "locality" = "locality"
    )
  )


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: sub-model-anom

anom_model_form <- update.formula(
  sub_model,
  . ~ . -
    region_name +
    temp_ano_winter +
    prec_ano_winter +
    temp_ano_winter:prec_ano_winter +
    temp_may_cur +
    prec_may_cur +
    temp_may_cur:prec_may_cur +
    temp_may_prev +
    prec_may_prev +
    temp_may_prev:prec_may_prev +
    temp_ano_month_prev +
    prec_ano_month_prev +
    temp_ano_month_prev:prec_ano_month_prev
)

anom_model <- lmer(anom_model_form,
  data = locality_sampling_biomass_temp_compl_light_week_weather_anom,
  #na.action = "na.fail",
  na.action = "na.omit",
  REML = F
)

anom_model_form_april <- update.formula(
  sub_model,
  . ~ . -
    region_name +
    temp_ano_winter +
    prec_ano_winter +
    temp_ano_winter:prec_ano_winter +
    temp_april_cur +
    prec_april_cur +
    temp_april_cur:prec_april_cur +
    temp_april_prev +
    prec_april_prev +
    temp_april_prev:prec_april_prev +
    temp_ano_month_prev +
    prec_ano_month_prev +
    temp_ano_month_prev:prec_ano_month_prev
)

anom_model_april <- lmer(anom_model_form_april,
  data = locality_sampling_biomass_temp_compl_light_week_weather_anom,
  #na.action = "na.fail",
  na.action = "na.omit",
  REML = F
)


## ------------------------------------------------------------------------------------------------------------------------------------
r2_verdi_biomasse_anom <- r.squaredGLMM(anom_model)


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-anom-mod
#| fig-cap: "Parameterestimat fra modell over værforhold utenfor innsamlingsperioden. Avhengig variabel er biomasse insekter per dag (avrendt våtvekt (g/dag)). Feillinjer viser 95 % konfidensintervall, og parametere i fet stil angir de som ikke overlapper med null effekt (0). Estimatene for variablene for værforhold under innsamlingsperioden var lite påvirket av de nye variablene og vises ikke her, for å forenkle figuren. Estimatet for årstrenden skal forstås som tidstrenden, når effektene av alle andre parametre har blitt tatt høyde for."

plot_model(anom_model,
  sort.est = TRUE,
  rm.terms = c(
    "scale_avg_values_MX2301A_temperature",
    "scale_avg_values_MX2301A_rh",
    "scale_avg_values_MX2202_light",
    "scale_stddev_values_MX2301A_temperature",
    "scale_stddev_values_MX2202_light",
    "scale_stddev_values_MX2301A_rh",
    "scale_avg_values_MX2301A_temperature:scale_avg_values_MX2301A_rh",
    "scale_avg_values_MX2301A_temperature:scale_avg_values_MX2202_light",
    "scale_avg_values_MX2301A_rh:scale_avg_values_MX2202_light",
    "scale_stddev_values_MX2301A_temperature:scale_stddev_values_MX2202_light"
  ),
  axis.labels = c(
    "scale_avg_values_MX2301A_temperature" = "<b style='color:black'>Temperatur</b>",
    "scale_avg_values_MX2301A_rh" = "<b style='color:black'>Luftfuktighet</b>",
    "scale_avg_values_MX2202_light" = "<b style='color:black'>Lys</b>",
    "scale_stddev_values_MX2301A_temperature" = "<b style='color:black'>Variasjon_temp.</b>",
    "scale_stddev_values_MX2202_light" = "<b style='color:black'>Variasjon_lys</b>",
    "scale_avg_values_MX2301A_temperature:scale_avg_values_MX2301A_rh" = "Temp:Fukt",
    "scale_avg_values_MX2301A_temperature:scale_avg_values_MX2202_light" = "Temp:Lys",
    "scale_avg_values_MX2301A_rh:scale_avg_values_MX2202_light" = "Fukt:Lys",
    "scale_stddev_values_MX2301A_temperature:scale_stddev_values_MX2202_light" = "<b style='color:black'>Variasjon_temp:Variasjon_lys</b>",
    "scale_stddev_values_MX2301A_rh" = "<b style='color:black'>Variasjon_fukt</b>",
    "year_num" = "<b style='color:black'>År</b>",
    "habitat_typeSemi-nat" = "<b style='color:black'>Habitat_Semi-nat</b>",
    "prec_ano_winter" = "<b style='color:black'>Nedbør_vinter</b>",
    "temp_ano_month_prev" = "Temperatur_-1år",
    "temp_may_prev" = "<b style='color:black'>Temperatur_mai_-1år</b>",
    "temp_may_cur:prec_may_cur" = "Temperatur_mai:Nedbør_mai",
    "temp_ano_month_prev:prec_ano_month_prev" = "Temperatur_-1år:Nedbør_-1år",
    "prec_ano_month_prev" = "Nedbør_-1år",
    "temp_may_prev:prec_may_prev" = "<b style='color:black'>Temperatur_mai_-1år:Nedbør_mai_-1år</b>",
    "temp_ano_winter:prec_ano_winter" = "Temperatur_vinter:Nedbør_vinter",
    "prec_may_cur" = "<b style='color:black'>Nedbør_mai</b>",
    "prec_may_prev" = "Nedbør_mai_-1år",
    "temp_ano_winter" = "Temperatur_vinter",
    "temp_may_cur" = "<b style='color:black'>Temperatur_mai</b>"
  ),
  title = ""
) +
  ggplot2::theme(axis.text.y = ggtext::element_markdown())


## ------------------------------------------------------------------------------------------------------------------------------------
##No tree age measurements was taken in Skog_49 (Sweco) 

global_model_big_forest <- lmer(
  log_bm_per_day ~
    scale_avg_values_MX2301A_temperature +
    scale_avg_values_MX2301A_rh +
    scale_avg_values_MX2202_light +
    scale_stddev_values_MX2301A_temperature +
    scale_stddev_values_MX2301A_rh +
    scale_stddev_values_MX2202_light +
    # habitat_type +
    scale_stddev_values_MX2301A_temperature:scale_stddev_values_MX2202_light +
    scale(ar5_PC1) +
    scale(ar5_PC2) +
    scale(lsm_l_shdi) +
    scale(tot_perimeter_length) +
    scale(skog_PC1) +
    scale(skog_PC2) +
    scale(ano_mean_no_spec) +
    scale(ano_mean_cover) +
    scale(avg_age) +
    dom_tree +
    year_num +
    (1 | locality) +
    (1 | Ukenummer),
  data = locality_sampling_biomass_big_join_forest,
  #na.action = "na.fail",
  na.action = "na.omit",
  REML = F,
  control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-4))
)


## ------------------------------------------------------------------------------------------------------------------------------------
tt <- summary(global_model_big_forest)
to_tab <- tt$coefficients %>%
  as_tibble(rownames = "Parameter") %>%
  mutate(`Estimate.norm` = exp(Estimate))

to_tab <- to_tab %>%
  mutate(Parameter = case_match(
    Parameter,
    "(Intercept)" ~ "Intercept",
    "scale_avg_values_MX2301A_temperature" ~ "Temperatur",
    "scale_avg_values_MX2301A_rh" ~ "Luftfuktighet",
    "scale_avg_values_MX2202_light" ~ "Lys",
    "scale_stddev_values_MX2301A_temperature" ~ "Variasjon_temp.",
    "scale_stddev_values_MX2301A_rh" ~ "Variasjon_fukt",
    "scale_stddev_values_MX2202_light" ~ "Variasjon_lys",
    "habitat_typeSemi-nat" ~ "Habitat_semi-nat",
    "year_num" ~ "År",
    "scale-avg_values_MX2301A_temperature:scale_avg_values_MX2301A_rh" ~ "Temperatur:Fukt",
    "scale_avg_values_MX2301A_temperature:scale_avg_values_MX2202_light" ~ "Temperatur:Lys",
    "scale_avg_values_MX2301A_rh:scale_avg_values_MX2202_light" ~ "Fukt:Lys",
    "scale_stddev_values_MX2301A_temperature:scale_stddev_values_MX2202_light" ~ "Variasjon_temp:Variasjon_lys",
    "dom_treeDunbjørk" ~ "Dom._treslag_Dunbjørk ",
    "dom_treeFuru" ~ "Dom._treslag_Furu",
    "dom_treeGran" ~ "Dom._treslag_Gran",
    "dom_treeHengebjørk" ~ "Dom._treslag_Hengebjørk",
    "dom_treeSelje" ~ "Dom._treslag_Selje",
    "scale(avg_age)" ~ "Landsskog_skogalder",
    "scale(ano_mean_cover)" ~ "ANO_dekningsgrad_karpl.",
    "scale(ano_mean_no_spec)" ~ "ANO_artsantall_karpl.",
    "scale(ar5_PC1)" ~ "AR5-arealtype_PCA-akse_1",
    "scale(ar5_PC2)" ~ "AR5-arealtype_PCA-akse_2",
    "scale(skog_PC1)" ~ "AR5-skog_PCA-akse_1",
    "scale(skog_PC2)" ~ "AR5-skog_PCA-akse_2",
    "scale(lsm_l_shdi)" ~ "AR5-arealtype_Shannon_indeks",
    "scale(tot_perimeter_length)" ~ "AR5-arealtype_lineære_elementer"
  )) %>%
  filter(!(Parameter %in% c(
    "Intercept",
    "Luftfuktighet",
    "Variasjon_fukt",
    "Temperatur",
    "Lys",
    "Variasjon_temp.",
    "Temperatur:Lys",
    "Fuktighet:Lys",
    "Variasjon_lys",
    "Variasjon_temp:Variasjon_lys"
  ))) %>%
  select(
    Parameter,
    Estimate.norm,
    Estimate,
    `Std. Error`,
    # df,
    `t value` # ,
    # `Pr(>|t|)`
  )


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: tbl-forest-mod
#| results: asis
#| tbl-cap: "Parameterestimat fra modell på insektbiomasse i skogslokaliteter. Værvariabler er ikke vist."
#| eval: true

if (names(rmarkdown::metadata$format) == "ninareport-pdf") {
  xtable(to_tab, digits = 3)
} else {
  kable(to_tab, digits = 3)
}


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-forest-mod
#| fig-cap: "Parameterestimat fra modell på insektbiomasse i skogslokaliteter. Værvariabler er ikke vist."

plot_model(global_model_big_forest,
  sort.est = TRUE,
  rm.terms = c(
    "scale_avg_values_MX2301A_temperature",
    "scale_avg_values_MX2301A_rh",
    "scale_avg_values_MX2202_light",
    "scale_stddev_values_MX2301A_temperature",
    "scale_avg_values_MX2301A_temperature:scale_avg_values_MX2202_light",
    "scale_stddev_values_MX2301A_temperature:scale_stddev_values_MX2202_light",
    "scale_avg_values_MX2301A_rh:scale_avg_values_MX2202_light",
    "scale_stddev_values_MX2202_light",
    "scale_stddev_values_MX2301A_rh"
  ),
  axis.labels = c(
    "year_num" = "År",
    "dom_treeHengebjørk" = "Dom._treslag_Hengebjørk",
    "dom_treeGran" = "Dom._treslag_Gran",
    "dom_treeFuru" = "Dom._treslag_Furu",
    "dom_treeSelje" = "Dom._treslag_Selje",
    "dom_treeDunbjørk" = "Dom._treslag_Dunbjørk",
    "scale(avg_age)" = "Landsskog_skogalder",
    "scale(ano_mean_cover)" = "ANO_dekningsgrad_karpl.",
    "scale(ano_mean_no_spec)" = "ANO_artsantall_karpl.",
    "scale(ar5_PC1)" = "AR5-arealtype_PCA-akse_1",
    "scale(ar5_PC2)" = "AR5-arealtype_PCA-akse_2",
    "scale(skog_PC1)" = "AR5-skog_PCA-akse_1",
    "scale(skog_PC2)" = "AR5-skog_PCA-akse_2",
    "scale(lsm_l_shdi)" = "AR5-arealtype_Shannon indeks",
    "scale(tot_perimeter_length)" = "AR5-arealtype_lineære elementer"
  ),
  title = ""
)


## ------------------------------------------------------------------------------------------------------------------------------------
r2_verdi_biomasse_big_forest <- r.squaredGLMM(global_model_big_forest)


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: !expr cache_var

global_model_big_semi_nat <- lmer(
  log_bm_per_day ~
    scale_avg_values_MX2301A_temperature +
    scale_avg_values_MX2301A_rh +
    scale_avg_values_MX2202_light +
    scale_stddev_values_MX2301A_temperature +
    scale_stddev_values_MX2202_light +
    scale_avg_values_MX2301A_temperature:scale_avg_values_MX2202_light +
    scale_avg_values_MX2301A_rh:scale_avg_values_MX2202_light +
    scale_stddev_values_MX2301A_temperature:scale_stddev_values_MX2202_light +
    scale(ar5_PC1) +
    scale(ar5_PC2) +
    scale(lsm_l_shdi) +
    scale(tot_perimeter_length) +
    scale(skog_PC1) +
    scale(skog_PC2) +
    scale(ano_mean_no_spec) +
    scale(ano_mean_cover) +
    region_name +
    year_num +
    (1 | locality) +
    (1 | Ukenummer),
  data = locality_sampling_biomass_big_join_semi_nat,
  na.action = "na.fail",
  REML = F,
  control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-4))
) # Suppress warnings due to near zero estimate for locality random effect


## ------------------------------------------------------------------------------------------------------------------------------------
#| eval: false
# summary(global_model_big_semi_nat)


## ------------------------------------------------------------------------------------------------------------------------------------
tt <- summary(global_model_big_semi_nat)
to_tab <- tt$coefficients %>%
  as_tibble(rownames = "Parameter") %>%
  mutate(`Estimate.norm` = exp(Estimate))

to_tab <- to_tab %>%
  mutate(Parameter = case_match(
    Parameter,
    "(Intercept)" ~ "Intercept",
    "scale_avg_values_MX2301A_temperature" ~ "Temperatur",
    "scale_avg_values_MX2301A_rh" ~ "Luftfuktighet",
    "scale_avg_values_MX2202_light" ~ "Lys",
    "scale_stddev_values_MX2301A_temperature" ~ "Variasjon_temp.",
    "scale_stddev_values_MX2301A_rh" ~ "Variasjon_fukt",
    "scale_stddev_values_MX2202_light" ~ "Variasjon_lys",
    "habitat_typeSemi-nat" ~ "Habitat_Semi-nat",
    "year_num" ~ "År",
    "scale_avg_values_MX2301A_temperature:scale_avg_values_MX2301A_rh" ~ "Temperatur:Fukt",
    "scale_avg_values_MX2301A_temperature:scale_avg_values_MX2202_light" ~ "Temperatur:Lys",
    "scale_avg_values_MX2301A_rh:scale_avg_values_MX2202_light" ~ "Fukt:Lys",
    "scale_stddev_values_MX2301A_temperature:scale_stddev_values_MX2202_light" ~ "Variasjon_temp:Variasjon_lys",
    "dom_treeDunbjørk" ~ "Dom._treslag_Dunbjørk ",
    "dom_treeFuru" ~ "Dom._treslag_Furu",
    "dom_treeGran" ~ "Dom._treslag_Gran",
    "dom_treeHengebjørk" ~ "Dom._treslag_Hengebjørk",
    "dom_treeSelje" ~ "Dom._treslag_Selje",
    "scale(avg_age)" ~ "Landsskog_skogalder",
    "scale(ano_mean_cover)" ~ "ANO_dekningsgrad_karpl.",
    "scale(ano_mean_no_spec)" ~ "ANO_artsantall_karpl.",
    "scale(ar5_PC1)" ~ "AR5-arealtype_PCA-akse_1",
    "scale(ar5_PC2)" ~ "AR5-arealtype_PCA-akse_2",
    "scale(skog_PC1)" ~ "AR5-skog_PCA-akse_1",
    "scale(skog_PC2)" ~ "AR5-skog_PCA-akse_2",
    "scale(lsm_l_shdi)" ~ "AR5-arealtype_Shannon_indeks",
    "scale(tot_perimeter_length)" ~ "AR5-arealtype_lineære_elementer",
    "region_nameSørlandet" ~ "Region_Sørlandet",
    "region_nameTrøndelag" ~ "Region_Trøndelag",
    "region_nameNord-Norge" ~ "Region_Nord-Norge",
    "region_nameVestlandet" ~ "Region_Vestlandet"
  )) %>%
  filter(!(Parameter %in% c(
    "Intercept",
    "Temperatur",
    "Lys",
    "Variasjon_temp.",
    "Tempertur:Lys",
    "Variasjon_temp:Variasjon_lys",
    "Fukt:Lys",
    "Luftfuktighet",
    "Variasjon_lys"
  ))) %>%
  select(
    Parameter,
    Estimate.norm,
    Estimate,
    `Std. Error`,
    # df,
    `t value` # ,
    # `Pr(>|t|)`
  )


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: tbl-seminat-mod
#| results: asis
#| tbl-cap: "Parameterestimat fra modell over insektbiomasse i semi-naturlig mark. Værvariabler er ikke vist."

if (names(rmarkdown::metadata$format) == "ninareport-pdf") {
  xtable(to_tab, digits = 3)
} else {
  kable(to_tab, digits = 3)
}


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-seminat-mod
#| fig-cap: "Parameterestimat fra modell over insektbiomasse i semi-naturlig mark. Værvariabler er ikke vist."
#|
plot_model(global_model_big_semi_nat,
  sort.est = TRUE,
  rm.terms = c(
    "scale_avg_values_MX2301A_temperature",
    "scale_avg_values_MX2301A_rh",
    "scale_avg_values_MX2202_light",
    "scale_stddev_values_MX2301A_temperature",
    "scale_avg_values_MX2301A_temperature:scale_avg_values_MX2202_light",
    "scale_stddev_values_MX2202_light",
    "scale_stddev_values_MX2301A_temperature:scale_stddev_values_MX2202_light",
    "scale_avg_values_MX2301A_rh:scale_avg_values_MX2202_light",
    "scale_stddev_values_MX2301A_rh"
  ),
  axis.labels = c(
    "year_num" = "År",
    "dom_treeHengebjørk" = "Dom._treslag_Hengebjørk",
    "dom_treeGran" = "Dom._treslag_Gran",
    "dom_treeFuru" = "Dom._treslagLFuru",
    "dom_treeSelje" = "Dom._treslag_Selje",
    "dom_treeDunbjørk" = "Dom._treslag_Dunbjørk",
    "scale(avg_age)" = "Landsskog_skogalder",
    "scale(ano_mean_cover)" = "ANO_dekningsgrad_karpl.",
    "scale(ano_mean_no_spec)" = "ANO_artsantall_karpl.",
    "scale(ar5_PC1)" = "AR5-arealtype_PCA-akse_1",
    "scale(ar5_PC2)" = "AR5-arealtype_PCA-akse_2",
    "scale(skog_PC1)" = "AR5-skog_PCA-akse_1",
    "scale(skog_PC2)" = "AR5-skog_PCA-akse_2",
    "scale(lsm_l_shdi)" = "AR5-arealtype_Shannon indeks",
    "scale(tot_perimeter_length)" = "AR5-arealtype_lineære_elementer",
    "region_nameTrøndelag" = "Region_Trøndelag",
    "region_nameNord-Norge" = "Region_Nord-Norge",
    "region_nameSørlandet" = "Region_Sørlandet",
    "region_nameVestlandet" = "Region_Vestlandet"
  ),
  title = ""
)


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: !expr cache_var
r2_verdi_biomasse_big_semi_nat <- r.squaredGLMM(global_model_big_semi_nat)


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: true
# Here summing all the no_reads per species_latin_fixed and sample

otu_loc_q <- "
SELECT l.locality, l.region_name, l.habitat_type, ls.sampling_name, ls.sampling_number, date_part('week', ls.end_date) as end_week, st.sample_name, o.species_latin_fixed as otu, sum(o.no_reads)::integer no_reads
FROM occurrences.observations o,
events.identifications i,
events.sampling_trap st,
events.locality_sampling ls,
events.year_locality yl,
locations.traps,
locations.localities l
WHERE o.identification_id = i.id
AND i.sampling_trap_id = st.id
AND st.locality_sampling_id = ls.id
AND ls.year_locality_id = yl.id
AND yl.locality_id = l.id
AND st.trap_id = traps.id
AND yl.year = '2024'
AND yl.project_short_name = 'NorIns'
AND traps.trap_short_name = 'MF1'
AND (o.taxonomic_level = 'species'
OR o.taxonomic_level = 'genus'
OR o.taxonomic_level = 'family'
OR o.taxonomic_level = 'order'
)
GROUP BY l.locality, l.region_name, l.habitat_type, ls.end_date, ls.sampling_name, ls.sampling_number, st.sample_name, o.species_latin_fixed

"

otu_loc_raw <- dbGetQuery(
  con,
  otu_loc_q
) %>%
  as_tibble() %>%
  mutate(habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type))


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: separate_ost_comm
otu_loc_ost <- otu_loc_raw %>%
  filter(region_name == "Østlandet") %>%
  select(-"sampling_number")

otu_loc_ost_wider <- otu_loc_ost %>%
  pivot_wider(
    names_from = otu,
    values_from = no_reads,
    values_fill = 0,
  ) %>%
  as.data.frame()

rownames(otu_loc_ost_wider) <- otu_loc_ost_wider$sampling_name
otu_loc_ost <- otu_loc_ost_wider[-c(1:6)]
# dim(otu_loc)


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: separate_seminat_comm

otu_loc_semi <- otu_loc_raw %>%
  filter(habitat_type == "Semi-nat") %>%
  filter(region_name %in% c("Trøndelag", "Nord-Norge", "Østlandet")) %>%
  select(-"sampling_number")

otu_loc_semi_wider <- otu_loc_semi %>%
  pivot_wider(
    names_from = otu,
    values_from = no_reads,
    values_fill = 0,
  ) %>%
  as.data.frame()

rownames(otu_loc_semi_wider) <- otu_loc_semi_wider$sampling_name
otu_loc_semi <- otu_loc_semi_wider[-c(1:6)]


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: mds_ost
#| message: false
#| results: hide
#| cache: true
ord_malaise_ost <- vegan::metaMDS(otu_loc_ost,
  distance = "bray",
  # trymax = 1000
  trymax = 200
)

ord_malaise_ost_aug <- cbind(otu_loc_ost_wider[c(1:6)], vegan::scores(ord_malaise_ost)[1])
# %>%
# mutate(sampling_number = factor(sampling_number),
#       habitat_type = factor(habitat_type))


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: mds_semi
#| message: false
#| results: hide
#| cache: true
ord_malaise_semi <- vegan::metaMDS(otu_loc_semi,
  distance = "bray",
  # trymax = 1000
  trymax = 400
)

ord_malaise_semi_aug <- cbind(otu_loc_semi_wider[c(1:6)], vegan::scores(ord_malaise_semi)[1])
# %>%
# mutate(sampling_number = factor(sampling_number),
#       habitat_type = factor(habitat_type))


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: true
fig_habitats_ost <-
  ggplot(
    ord_malaise_ost_aug,
    aes(
      x = sites.NMDS1,
      y = sites.NMDS2,
      shape = habitat_type,
      color = end_week
    )
  ) +
  geom_point() +
  stat_ellipse(
    aes(
      x = sites.NMDS1,
      y = sites.NMDS2,
      group = habitat_type
    ),
    show.legend = F,
    # data = ord_malaise_ost_aug,
    inherit.aes = F
  ) +
  scale_color_nina(
    name = "Tømmeuke",
    discrete = F, palette = "blue-orange"
  ) +
  scale_shape_discrete(name = "Habitattype") +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  )


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: true
fig_seminat_mds <-
  ggplot(
    ord_malaise_semi_aug,
    aes(
      x = sites.NMDS1,
      y = sites.NMDS2,
      shape = region_name,
      color = end_week
    )
  ) +
  geom_point() +
  stat_ellipse(
    aes(
      x = sites.NMDS1,
      y = sites.NMDS2,
      group = region_name
    ),
    show.legend = F,
    data = ord_malaise_semi_aug,
    inherit.aes = F
  ) +
  scale_color_nina(
    name = "Tømmeuke",
    discrete = F, palette = "blue-orange"
  ) +
  scale_shape_discrete(name = "Region") +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  )


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: true
#| label: get_otu_year_loc
# Here summing all the no_reads per species_latin_fixed and locality

otu_year_loc_q <- "
SELECT l.locality, l.region_name, l.habitat_type, o.species_latin_fixed as otu, sum(o.no_reads)::integer no_reads
FROM occurrences.observations o,
events.identifications i,
events.sampling_trap st,
events.locality_sampling ls,
events.year_locality yl,
locations.traps,
locations.localities l
WHERE o.identification_id = i.id
AND i.sampling_trap_id = st.id
AND st.locality_sampling_id = ls.id
AND ls.year_locality_id = yl.id
AND yl.locality_id = l.id
AND st.trap_id = traps.id
AND yl.year = '2024'
AND yl.project_short_name = 'NorIns'
AND traps.trap_short_name = 'MF1'
AND (o.taxonomic_level = 'species'
OR o.taxonomic_level = 'genus'
OR o.taxonomic_level = 'family'
OR o.taxonomic_level = 'order'
)
GROUP BY l.locality, l.region_name, l.habitat_type, o.species_latin_fixed

"

otu_year_loc_raw <- dbGetQuery(
  con,
  otu_year_loc_q
) %>%
  as_tibble() %>%
  mutate(habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type))


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: separate_ost_comm_year_loc
otu_year_loc_ost <- otu_year_loc_raw %>%
  filter(region_name == "Østlandet")

otu_year_loc_ost_wider <- otu_year_loc_ost %>%
  pivot_wider(
    names_from = otu,
    values_from = no_reads,
    values_fill = 0,
  ) %>%
  as.data.frame()

rownames(otu_year_loc_ost_wider) <- otu_year_loc_ost_wider$locality
otu_year_loc_ost <- otu_year_loc_ost_wider[-c(1:4)]
# dim(otu_loc)


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: separate_seminat_comm_year_loc

otu_year_loc_semi <- otu_year_loc_raw %>%
  filter(habitat_type == "Semi-nat") %>%
  filter(region_name %in% c("Trøndelag", "Nord-Norge", "Østlandet", "Sørlandet", "Vestlandet")) 

otu_year_loc_semi_wider <- otu_year_loc_semi %>%
  pivot_wider(
    names_from = otu,
    values_from = no_reads,
    values_fill = 0,
  ) %>%
  as.data.frame()

rownames(otu_year_loc_semi_wider) <- otu_year_loc_semi_wider$locality
otu_year_loc_semi <- otu_year_loc_semi_wider[-c(1:4)]


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: mds_ost_year_loc
#| message: false
#| results: hide
#| cache: true
ord_malaise_ost_year_loc <- vegan::metaMDS(otu_year_loc_ost,
  distance = "bray",
  # trymax = 1000
  trymax = 200
)

ord_malaise_ost_year_loc_aug <- cbind(otu_year_loc_ost_wider[c(1:3)], vegan::scores(ord_malaise_ost_year_loc)[1])
# %>%
# mutate(sampling_number = factor(sampling_number),
#       habitat_type = factor(habitat_type))


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: mds_semi_year_loc
#| message: false
#| results: hide
#| cache: true
ord_malaise_semi_year_loc <- vegan::metaMDS(otu_year_loc_semi,
  distance = "bray",
  # trymax = 1000
  trymax = 800
)

ord_malaise_semi_year_loc_aug <- cbind(otu_year_loc_semi_wider[c(1:3)], vegan::scores(ord_malaise_semi_year_loc)[1])
# %>%
# mutate(sampling_number = factor(sampling_number),
#       habitat_type = factor(habitat_type))


## ------------------------------------------------------------------------------------------------------------------------------------
#| message: false
#| results: hide
#| cache: true

loc_expl_var <- locality_sampling_biomass_big_join_semi_nat |> 
  group_by(locality) |> 
  summarize(mean_temp = mean(scale_avg_values_MX2301A_temperature, na.rm = T),
            mean_light = mean(scale_avg_values_MX2202_light, na.rm = T) ,
            mean_rh = mean(scale_avg_values_MX2301A_rh, na.rm = T),
            lsm_l_shdi = mean(lsm_l_shdi),
            tot_perimeter_length = mean(tot_perimeter_length),
            skog_PC1 = mean(skog_PC1),
            skog_PC2 = mean(skog_PC1),
            ar5_PC1 = mean(ar5_PC1),
            ar5_PC2 = mean(ar5_PC2),
            ano_mean_no_spec = mean(ano_mean_no_spec),
            ano_mean_cover = mean(ano_mean_cover))



points <- ord_malaise_semi_year_loc_aug %>%
  left_join(loc_expl_var,
    by = c("locality" = "locality")
  ) 

# fit <- vegan::envfit(ord_malaise_semi ~ avg_values_MX2301A_temperature + avg_values_MX2202_light + ar5_PC1 + ar5_PC2,
#                     data = points,
#                     na.rm = TRUE)

fit <- vegan::envfit(ord_malaise_semi_year_loc ~ mean_temp + mean_light + mean_rh + lsm_l_shdi + tot_perimeter_length + skog_PC1 + skog_PC2 + ar5_PC1 + ar5_PC2 + ano_mean_no_spec + ano_mean_cover,
  data = points,
  na.rm = TRUE
)

fit_scores <- data.frame(vegan::scores(fit, display = "vectors"))
fit_scores <- cbind(fit_scores, vars = rownames(fit_scores))


## ------------------------------------------------------------------------------------------------------------------------------------
#| results: hide
#| warning: false
#| message: false
#| cache: !expr cache_var
fig_seminat <- ggplot(points) +
  geom_point(
    data = points,
    aes(
      x = sites.NMDS1, y = sites.NMDS2,
      color = region_name
    )
  ) +
  # scale_color_manual(name = "Lokalitet",
  #                    values=paste0(Polychrome::palette36.colors())) +
  scale_color_nina(name = "Region") +
  coord_fixed() +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 8),#[sjekk størrelse]
    legend.key.width = unit(0.5, "cm")
      ) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  guides(color = guide_legend(override.aes = list(size = 5), nrow = 2, byrow = T)) +
  coord_flip()

fig_seminat <- fig_seminat +
  geom_segment(
    data = fit_scores[c(1:3, 5, 7, 9, 10, 11), ],
    aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
    show.legend = F,
    arrow = arrow(length = unit(0.25, "cm"))
  ) +
  geom_segment(
    data = fit_scores[c(4, 6, 8), ],
    aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
    linetype = 2,
    show.legend = F,
    arrow = arrow(length = unit(0.25, "cm"))
  ) +
  geom_text(
    data = fit_scores,
    # aes(x = NMDS1-0.05, y = NMDS2-0.05, label = c("temp","lys","ar5_PC1","ar5_PC2")),
    aes(x = NMDS1 - 0.05, y = NMDS2 - 0.05, label = c("temp", "lys", "rh", "shannon", "perim", "skog_PC1", "skog_PC2", "ar5_PC1", "ar5_PC2", "ano_spec", "ano_cov")),
    size = 3
  ) +
  ggtitle("Semi-Nat")


## ------------------------------------------------------------------------------------------------------------------------------------
#| message: false
#| results: hide
#| cache: true


points <- ord_malaise_ost_year_loc_aug %>%
  left_join(loc_expl_var,
    by = c("locality" = "locality")
  ) 

# fit <- vegan::envfit(ord_malaise_ost ~ avg_values_MX2301A_temperature + avg_values_MX2202_light + ar5_PC1 + ar5_PC2,
#                      data = points,
#                      na.rm = TRUE)

fit <- vegan::envfit(ord_malaise_ost_year_loc ~ mean_temp + mean_light + mean_rh + lsm_l_shdi + tot_perimeter_length + skog_PC1 + skog_PC2 + ar5_PC1 + ar5_PC2 + ano_mean_no_spec + ano_mean_cover,
  data = points,
  na.rm = TRUE
)



fit_scores <- data.frame(vegan::scores(fit, display = "vectors"))
fit_scores <- cbind(fit_scores, vars = rownames(fit_scores))


## ------------------------------------------------------------------------------------------------------------------------------------
#| results: hide
#| warning: false
#| message: false
fig_ost <- ggplot(points) +
  geom_point(
    data = points,
    aes(
      x = sites.NMDS1, y = sites.NMDS2,
      color = habitat_type
    )
  ) +
  # scale_color_manual(name = "Lokalitet",
  #                    values=paste0(Polychrome::palette36.colors())) +
  scale_color_nina(name = "Lokalitet") +
  coord_fixed() +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 8), #[sjekk størrelse]
    legend.key.width = unit(0.5, "cm")
  ) +
  xlab("NMDS1") +
  ylab("NMDS2")

fig_ost <- fig_ost +
  geom_segment(
    data = fit_scores[c(1:4, 6:11), ],
    aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
    show.legend = F,
    arrow = arrow(length = unit(0.25, "cm"))
  ) +
  geom_segment(
    data = fit_scores[5, ],
    aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
    linetype = 2,
    show.legend = F,
    arrow = arrow(length = unit(0.25, "cm"))
  ) +
  geom_text(
    data = fit_scores,
    aes(x = NMDS1 - 0.05, y = NMDS2 - 0.05, label = c("temp", "lys", "rh", "shannon", "perim", "skog_PC1", "skog_PC2", "ar5_PC1", "ar5_PC2", "ano_spec", "ano_cov")),
    size = 3
  ) +
  ggtitle("Østlandet")


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-beta-div-patterns-ost
#| fig-pos: '!ht'
#| fig-cap: "Ordinasjon som viser hovedmønstrene i artssammensetning i insektprøvene innsamlet på Østlandet i 2024. Punktene er farget etter innsamlingsdato og ellipsene viser et 95% konfidens-intervall for hvert habitat."
fig_habitats_ost


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-beta-div-patterns-semi
#| fig-pos: '!ht'
#| fig-cap: "Ordinasjon som viser hovedmønstrene i artssammensetning i insektprøvene innsamlet på i semi-naturlig mark i 2024. Punktene er farget etter innsamlingsdato og ellipsene viser et 95% konfidens-intervall for hver region. Data fra Sørlandet og Vestlandet er ikke tatt med da rutinen ikke ville finne en løsning med disse inkludert."
fig_seminat_mds


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-betadiv-pattern-semi-nat
#| fig-pos: '!hb'
#| fig-height: 7
#| fig-width: 4
#| fig-cap: "Ordinasjon som viser sammenheng mellom forklaringsvariabler og mønster i samfunnskomposisjon i semi-naturlig mark. Aksene har bytt posisjon og y aksen retning for å etterligne den geografiske plasseringen av regionene. Solide linjer viser signifikante sammenhenger med forklaringsvariabler. 'rh', 'lys', og 'temp' står henholdsvis for middelverdiene av luftfuktighet, lysintensitet og temperatur over hele fellesesongen. 'ar5_PC1' står for PCA-akse 1 i arealtyper der positive verdier er angir med større mengde skog. 'ar5_PC2' er på samme måte angir større mengde åpen fastmark og mindre mengde fulldyrka jord. 'Skog_PC2' står for PCA-akse 2 for skogstype der positive verdier angir høyere bonitet. 'ano_cov' angir høyere dekningsgrad av karplanter og 'ano_spec' et høyere artsantall karplanter fra ANO-kartleggingen. 'perim' står for mer lineære element i landskapet og 'shannon' angir landskapsheterogenitet (shannon indeks for arealtyper)."
fig_seminat


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-betadiv-pattern-ost
#| fig-pos: '!hb'
#| fig-height: 5
#| fig-cap: "Ordinasjon som viser sammenheng mellom forklaringsvariabler og mønster i samfunnskomposisjon på Østlandet. Solide linjer viser signifikante sammenhenger med forklaringsvariabler. 'rh', 'lys', og 'temp' står henholdsvis for luftfuktighet, lysintensitet og temperatur. 'ar5_PC1' (mengde skog), 'ar5_PC2' (mengde fulldyrka jord vs fastmark), 'ano_spec' (artsantall karplanter), 'ano_cov' (dekningsgrad karplanter) ligger tett på hverandre oppe til høyre. 'skog_PC1' (mengde barskog), 'skog_PC2' (bonitet), 'shannon' (shannon indeks for arealtyper) og 'perim' (lengde lineære element) ligger tett på hverandre nede til venstre."
fig_ost


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: get_biodiv_data_to_gam
#| cache: !expr cache_var
# Get biodiv info on locality & year level
year_loc_biodiv <- get_observations(dataset = "NorIns",
                                    agg_level = "year_locality",
                                    trap_type = "MF",
                                    subset_year = 2020:2024)
  

#year_loc_biodiv


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: get_biomass_data_to_gam
#| cache: !expr cache_var

year_loc_biomass <- get_biomass(trap_type = "MF",
                                dataset = "NorIns",
                                agg_level = "year_locality",
                                subset_year = 2020:2024) 
#year_loc_biomass


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: get_loc_data_to_gam
#| 
# Get locale info
#Here we get the lon, lat for the centroid of the locality squares. 
#loc_info <- Norimon::get_localities(dataset = "NorIns")

loc_info <- localities |> 
  mutate(centroid_lon = sf::st_coordinates(st_transform(st_centroid(geom), 4326))[, 1],
         centroid_lat = sf::st_coordinates(st_transform(st_centroid(geom), 4326))[, 2])
#loc_info


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: get_trap_data_to_gam
# Get more detailed trap position info. 
traps <- read_sf(con,
             Id(schema = "locations",
                table = "traps")) |> 
  filter(trap_short_name == "MF1") |> #only malaisetraps
  mutate(trap_lon = sf::st_coordinates(st_transform(geom, 4326))[, 1],
         trap_lat = sf::st_coordinates(st_transform(geom, 4326))[, 2]) |> 
  select(locality,
         coordinate_precision_m,
         elev_m,
         trap_lon,
         trap_lat) |> 
  st_drop_geometry() |> 
  collect()


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: join_for_gam
# Join
loc_stripped <- loc_info |> 
  st_drop_geometry() |> 
  select(locality,
         centroid_lon,
         centroid_lat)

loc_sum_stats <- year_loc_biodiv |> 
  inner_join(year_loc_biomass) |> 
  left_join(loc_stripped) |> 
  left_join(traps)
  
loc_sum_stats_ost <- loc_sum_stats |> 
  filter(region_name == "Østlandet")
#loc_sum_stats_ost


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: load_pal_for_gam

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: get_hoyde_rast
### get altitude data ###
rastlist <- list.files(path = "../../../../store_filer_for_rapport/hoyde_ostland", pattern='*.tif',
                       all.files=TRUE, full.names=T)
allrasters <- sprc(lapply(rastlist, rast))
hoyde <- mosaic(allrasters)
#plot(hoyde)


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: mk_elev_rast
#| cache: true
### make rasters with x and y coordinates ###
coords <- c(init(hoyde, 'x'), init(hoyde, 'y'))
names(coords) <- c('x', 'y') # I will give the names only to keep it organized
coords[[1]] -> coords_x 
coords[[2]] -> coords_y 

### get map of ostlandet ####
map <- ne_countries(country = "norway", type = "countries", scale="large")
map2 <- ne_states(country = "norway")
map2[which(map2$gns_name %in% c("Hedmark", "Akershus", "Ostfold", "Oslo", "Buskerud", "Vestfold", "Telemark", "Oppland")),] -> map3

# make habitatmap - only seminat, for predicting
hoyde -> habitat
as.factor(habitat) -> habitat
values(habitat) <- as.factor("Semi-nat")
#ggplot() +  geom_spatraster(data=habitat)

## prep map for predictions ##
single_raster <- c(hoyde, coords_x, coords_y, habitat)
names(single_raster) <- c("h", "x", "y", "habitat_type")
# crop on ostlandet
map_ <- st_transform(map, crs = crs(single_raster, proj=T))
map2_ <- st_transform(map2, crs = crs(single_raster, proj=T))
map3_ <- st_transform(map3, crs = crs(single_raster, proj=T))
st_union(map3_) -> map33_

#### data massage ####
dat <- as.data.frame(loc_sum_stats_ost)
for (i in 5:length(dat)) as.numeric(dat[,i]) -> dat[,i]
projcrs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
dat <- st_as_sf(x = dat, coords = c("trap_lon", "trap_lat"), crs = projcrs)
dat <- st_transform(dat, crs = crs(hoyde, proj=T))
sf::st_coordinates(dat)[,1] -> dat$x
sf::st_coordinates(dat)[,2] -> dat$y
dat$elev_m -> dat$h
as.factor(dat$habitat_type) -> dat$habitat_type
dat$habitat_type <- factor(dat$habitat_type, levels = c("Semi-nat","Forest"))


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: mk_gam_plots
#| cache: false
### plot data ###


p1 <- ggplot() + geom_sf(data=map2, fill = "#FAFAFA") +
  geom_sf(data=map3, fill = "grey", alpha=.4) + 
  xlim(4.5,13) + ylim(58,63.4)  +
  geom_sf(data=dat, aes(color = no_species)) +
  scale_color_gradientn(colours = myPalette(100), name = "", limits=c(min(dat$no_species)-10, max(dat$no_species)+10)) +
  geom_sf(data=dat, colour="black", pch=21) +
  theme_bw() + ggtitle("Antall arter")
p2 <- ggplot() + geom_sf(data=map2, fill = "#FAFAFA") +
  geom_sf(data=map3, fill = "grey", alpha=.4) + 
  xlim(4.5,13) + ylim(58,63.4)  +
  geom_sf(data=dat, aes(color = shannon_div)) +
  scale_color_gradientn(colours = myPalette(100), name = "", limits=c(min(dat$shannon_div)-10, max(dat$shannon_div)+10)) +
  geom_sf(data=dat, colour="black", pch=21) +
  theme_bw() + ggtitle("Shannon diversitet")
p3 <- ggplot() + geom_sf(data=map2, fill = "#FAFAFA") +
  geom_sf(data=map3, fill = "grey", alpha=.4) + 
  xlim(4.5,13) + ylim(58,63.4)  +
  geom_sf(data=dat, aes(color = avg_wet_weight)) +
  scale_color_gradientn(colours = myPalette(100), name = "", limits=c(min(dat$avg_wet_weight)-3, max(dat$avg_wet_weight)+3)) +
  geom_sf(data=dat, colour="black", pch=21) +
  theme_bw() + ggtitle("Biomasse (g)")
p4 <- ggplot() + geom_sf(data=map2, fill = "#FAFAFA") +
  geom_sf(data=map3, fill = "grey", alpha=.4) + 
  xlim(4.5,13) + ylim(58,63.4)  +
  geom_sf(data=dat, aes(color = GDE_by_asv)) +
  scale_color_gradientn(colours = myPalette(100), name = "", limits=c(min(dat$GDE_by_asv)-.3, max(dat$GDE_by_asv)+.3)) +
  geom_sf(data=dat, colour="black", pch=21) +
  theme_bw() + ggtitle("GDE asv")


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: no_spec_gam
#| cache: false
# model antall arter #
#gg1 <- gam(no_species~s(x)+s(y)+s(h)+s(x,y)+s(x,h)+s(y,h)+habitat_type, data=dat)
gg1 <- gam(no_species~s(x)+s(y)+s(h)+habitat_type, data=dat)
#summary(gg1)
#print(paste("Rsqr:", round(summary(gg1)$r.sq, 2)))# Rsqr: 0.42
# residual-plot:
dat$resid_quant <- rep("[10. - 90.] persentil",length(dat$year))
dat$resid_quant[which(resid(gg1) < as.numeric(quantile(resid(gg1), 0.1)))] <- "10. persentil"
dat$resid_quant[which(resid(gg1) > as.numeric(quantile(resid(gg1), 0.9)))] <- "90. persentil"
as.factor(dat$resid_quant) -> dat$resid_quant
dat$resid_quant <- factor(dat$resid_quant, levels = c("90. persentil","[10. - 90.] persentil","10. persentil"))
p111 <- ggplot() + geom_sf(data=map2, fill = "#FAFAFA") +
  geom_sf(data=map3, fill = "lightgrey", alpha=.4) + 
  xlim(4.5,13) + ylim(58,63.4)  +
  geom_sf(data=dat, aes(size = no_species,color=resid_quant), alpha=0.7) +
  scale_colour_manual(name="Residualer", values = c("darkred", "darkgrey", "darkgreen")) +
  geom_sf(data=dat, aes(size = no_species), colour="black", pch=21) +
  theme_bw() + ggtitle("Antall arter") +
  labs(size="Antall arter") +
  scale_size_continuous(range = c(0.6, 5))
#
ef <- ggemmeans(gg1, term="y")
p1y <- plot(ef, show_data=T, dot_size=1) + theme_bw() + ggtitle("") + xlab("Breddegrad") + ylab("Antall arter") +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(axis.text.y=element_text(size=5)) +
  theme(axis.text.x=element_text(size=5))
ef <- ggemmeans(gg1, term="x")
p1x <- plot(ef, show_data=T, dot_size=1) + theme_bw() + ggtitle("") + xlab("Lengdegrad") + ylab("") +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(axis.text.y=element_text(size=5)) +
  theme(axis.text.x=element_text(size=5))
ef <- ggemmeans(gg1, term="h")
p1h <- plot(ef, show_data=T, dot_size=1) + theme_bw() + ggtitle("") + xlab("moh") + ylab("") +
  theme(axis.text.y=element_text(size=5)) +
  theme(axis.text.x=element_text(size=5))
ef <- ggemmeans(gg1, term="habitat_type")
p1ht <- plot(ef, show_data=T, dot_size=1) + theme_bw()


## ------------------------------------------------------------------------------------------------------------------------------------
### read raster ###
pred_r1 <- rast('../../../../store_filer_for_rapport/pred_r1.tif')


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: plot_no_spec_gam
#| cache: false
### plot antall arter - main figure ###
p11 <- ggplot() + geom_sf(data=map2_, fill = "#f5f8f5") +
  geom_spatraster(data=pred_r1) +   xlim(56875,387075) + ylim(6529675,6961825)  +
  scale_fill_gradientn(colours = myPalette(100), limits=c(minmax(pred_r1)[1,], minmax(pred_r1)[2,]),
                       name = "", na.value = "transparent") + 
  geom_sf(data=map33_, fill = NA, color="#03841c", linewidth=0.9) + 
  theme_bw() + ggtitle("Antall arter - romlig prediksjon")
leg1 <- ggpubr::get_legend(p11)
p1_ <- p1 + theme(legend.position = "none")
p11_ <- p11 + theme(legend.position = "none")


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: mk_shann_gam
#| cache: false
### model shannon ###
#gg2 <- gam(shannon_div~s(x)+s(y)+s(h)+s(x,y)+s(x,h)+s(y,h)+habitat_type, data=dat)
gg2 <- gam(shannon_div~s(x)+s(y)+s(h)+habitat_type, data=dat)
#summary(gg2)
#print(paste("Rsqr:", round(summary(gg2)$r.sq, 2))) #Rsqr: 0.47
gg2p <- gam.hp::gam.hp(gg2, type="adjR2")$hierarchical.partitioning
#
dat$resid_quant2 <- rep("[10. - 90.] persentil",length(dat$year))
dat$resid_quant2[which(resid(gg2) < as.numeric(quantile(resid(gg2), 0.1)))] <- "10. persentil"
dat$resid_quant2[which(resid(gg2) > as.numeric(quantile(resid(gg2), 0.9)))] <- "90. persentil"
as.factor(dat$resid_quant2) -> dat$resid_quant2
dat$resid_quant2 <- factor(dat$resid_quant2, levels = c("90. persentil","[10. - 90.] persentil","10. persentil"))
p222 <- ggplot() + geom_sf(data=map2, fill = "#FAFAFA") +
  geom_sf(data=map3, fill = "grey", alpha=.4) + 
  xlim(4.5,13) + ylim(58,63.4)  +
  geom_sf(data=dat, aes(size = shannon_div,color=resid_quant2), alpha=0.7) +
  scale_colour_manual(name="Residualer", values = c("darkred", "darkgrey", "darkgreen")) +
  geom_sf(data=dat, aes(size = shannon_div), colour="black", pch=21) +
  theme_bw() + ggtitle("Shannon diversitet") +
  labs(size="Shannon") +
  scale_size_continuous(range = c(0.6, 5))
#
ef <- ggemmeans(gg2, term="y")
p2y <- plot(ef, show_data=T, dot_size=1) + theme_bw() + ggtitle("") + xlab("Breddegrad") + ylab("Shannon diversitet") +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(axis.text.y=element_text(size=5)) +
  theme(axis.text.x=element_text(size=5))
ef <- ggemmeans(gg2, term="x")
p2x <- plot(ef, show_data=T, dot_size=1) + theme_bw() + ggtitle("") + xlab("Lengdegrad") + ylab("") +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(axis.text.y=element_text(size=5)) +
  theme(axis.text.x=element_text(size=5))
ef <- ggemmeans(gg2, term="h")
p2h <- plot(ef, show_data=T, dot_size=1) + theme_bw() + ggtitle("") + xlab("moh") + ylab("") +
  theme(axis.text.y=element_text(size=5)) +
  theme(axis.text.x=element_text(size=5))
ef <- ggemmeans(gg2, term="habitat_type")
p2ht <- plot(ef, show_data=T, dot_size=1) + theme_bw()


## ------------------------------------------------------------------------------------------------------------------------------------
### read raster ###
pred_r2 <- rast('../../../../store_filer_for_rapport/pred_r2.tif')


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: plot_shann_gam
#| cache: false
# plot
p22 <- ggplot() + geom_sf(data=map2_, fill = "#f5f8f5") +
  geom_spatraster(data=pred_r2) +   xlim(56875,387075) + ylim(6529675,6961825)  +
  scale_fill_gradientn(colours = myPalette(100), limits=c(minmax(pred_r2)[1,], minmax(pred_r2)[2,]),
                       name = "", na.value = "transparent") + 
  geom_sf(data=map33_, fill = NA, color="#03841c", linewidth=0.9) + 
  theme_bw() + ggtitle("Shannon diversitet - romlig prediksjon")
leg2 <- ggpubr::get_legend(p22)
p2_ <- p2 + theme(legend.position = "none")
p22_ <- p22 + theme(legend.position = "none")


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: mk_resid_gam
#| cache: false
### model biomass ###
#gg3 <- gam(avg_wet_weight~s(x)+s(y)+s(h)+s(x,y)+s(x,h)+s(y,h)+habitat_type, data=dat)
gg3 <- gam(avg_wet_weight~s(x)+s(y)+s(h)+habitat_type, data=dat)
#summary(gg3)
#print(paste("Rsqr:", round(summary(gg3)$r.sq, 2))) # Rsqr: 0.06
gg3p <- gam.hp::gam.hp(gg3, type="adjR2")$hierarchical.partitioning
#gg3p[which(rownames(gg3p) == "habitat_type"), which(colnames(gg3p) == "I.perc(%)")]
#
dat$resid_quant3 <- rep("[10. - 90.] persentil",length(dat$year))
dat$resid_quant3[which(resid(gg3) < as.numeric(quantile(resid(gg3), 0.1)))] <- "10. persentil"
dat$resid_quant3[which(resid(gg3) > as.numeric(quantile(resid(gg3), 0.9)))] <- "90. persentil"
as.factor(dat$resid_quant3) -> dat$resid_quant3
dat$resid_quant3 <- factor(dat$resid_quant3, levels = c("90. persentil","[10. - 90.] persentil","10. persentil"))
p333 <- ggplot() + geom_sf(data=map2, fill = "#FAFAFA") +
  geom_sf(data=map3, fill = "grey", alpha=.4) + 
  xlim(4.5,13) + ylim(58,63.4)  +
  geom_sf(data=dat, aes(size = avg_wet_weight,color=resid_quant3), alpha=0.7) +
  scale_colour_manual(name="Residualer", values = c("darkred", "darkgrey", "darkgreen")) +
  geom_sf(data=dat, aes(size = avg_wet_weight), colour="black", pch=21) +
  theme_bw() + ggtitle("Biomasse") +
  labs(size="Biomasse (g)") +
  scale_size_continuous(range = c(0.6, 5))
#
ef <- ggemmeans(gg3, term="y")
p3y <- plot(ef, show_data=T, dot_size=1) + theme_bw() + ggtitle("") + xlab("Breddegrad") + ylab("Biomasse (g)") +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(axis.text.y=element_text(size=5)) +
  theme(axis.text.x=element_text(size=5))
ef <- ggemmeans(gg3, term="x")
p3x <- plot(ef, show_data=T, dot_size=1) + theme_bw() + ggtitle("") + xlab("Lengdegrad") + ylab("") +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(axis.text.y=element_text(size=5)) +
  theme(axis.text.x=element_text(size=5))
ef <- ggemmeans(gg3, term="h")
p3h <- plot(ef, show_data=T, dot_size=1) + theme_bw() + ggtitle("") + xlab("moh") + ylab("") +
  theme(axis.text.y=element_text(size=5)) +
  theme(axis.text.x=element_text(size=5))
ef <- ggemmeans(gg3, term="habitat_type")
p3ht <- plot(ef, show_data=T, dot_size=1) + theme_bw()


## ------------------------------------------------------------------------------------------------------------------------------------
### read raster ###
pred_r3 <- rast('../../../../store_filer_for_rapport/pred_r3.tif')


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: plot_biomass_gam
### plot ###
p33 <- ggplot() + geom_sf(data=map2_, fill = "#f5f8f5") +
  geom_spatraster(data=pred_r3) +   xlim(56875,387075) + ylim(6529675,6961825)  +
  scale_fill_gradientn(colours = myPalette(100), limits=c(minmax(pred_r3)[1,], minmax(pred_r3)[2,]),
                       name = "", na.value = "transparent") + 
  geom_sf(data=map33_, fill = NA, color="#03841c", linewidth=0.9) + 
  theme_bw() + ggtitle("Biomasse - romlig prediksjon")
leg3 <- ggpubr::get_legend(p33)
p3_ <- p3 + theme(legend.position = "none")
p33_ <- p33 + theme(legend.position = "none")


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: plot_resid_gam
### model GDE ###
#gg4 <- gam(GDE_by_asv~s(x)+s(y)+s(h)+s(x,y)+s(x,h)+s(y,h)+habitat_type, data=dat)
gg4 <- gam(GDE_by_asv~s(x)+s(y)+s(h)+habitat_type, data=dat)
#summary(gg4)
#print(paste("Rsqr:", round(summary(gg4)$r.sq, 2))) # Rsqr: 0.06
#
dat$resid_quant4 <- rep("[10. - 90.] persentil",length(dat$year))
dat$resid_quant4[which(resid(gg4) < as.numeric(quantile(resid(gg4), 0.1)))] <- "10. persentil"
dat$resid_quant4[which(resid(gg4) > as.numeric(quantile(resid(gg4), 0.9)))] <- "90. persentil"
as.factor(dat$resid_quant4) -> dat$resid_quant4
dat$resid_quant4 <- factor(dat$resid_quant4, levels = c("90. persentil","[10. - 90.] persentil","10. persentil"))
p444 <- ggplot() + geom_sf(data=map2, fill = "#FAFAFA") +
  geom_sf(data=map3, fill = "grey", alpha=.4) + 
  xlim(4.5,13) + ylim(58,63.4)  +
  geom_sf(data=dat, aes(size = GDE_by_asv,color=resid_quant4), alpha=0.7) +
  scale_colour_manual(name="Residualer", values = c("darkred", "darkgrey", "darkgreen")) +
  geom_sf(data=dat, aes(size = GDE_by_asv), colour="black", pch=21) +
  theme_bw() + ggtitle("GDE asv") +
  labs(size="GDE asv") +
  scale_size_continuous(range = c(0.6, 5))
#
ef <- ggemmeans(gg4, term="y")
p4y <- plot(ef, show_data=T, dot_size=1) + theme_bw() + ggtitle("") + xlab("Breddegrad") + ylab("Biomasse (g)") +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(axis.text.y=element_text(size=5)) +
  theme(axis.text.x=element_text(size=5))
ef <- ggemmeans(gg4, term="x")
p4x <- plot(ef, show_data=T, dot_size=1) + theme_bw() + ggtitle("") + xlab("Lengdegrad") + ylab("") +
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
  theme(axis.text.y=element_text(size=5)) +
  theme(axis.text.x=element_text(size=5))
ef <- ggemmeans(gg4, term="h")
p4h <- plot(ef, show_data=T, dot_size=1) + theme_bw() + ggtitle("") + xlab("moh") + ylab("") +
  theme(axis.text.y=element_text(size=5)) +
  theme(axis.text.x=element_text(size=5))
ef <- ggemmeans(gg4, term="habitat_type")
p4ht <- plot(ef, show_data=T, dot_size=1) + theme_bw()


## ------------------------------------------------------------------------------------------------------------------------------------
### read raster ###
pred_r4 <- rast('../../../../store_filer_for_rapport/pred_r4.tif')


## ------------------------------------------------------------------------------------------------------------------------------------
### plot ###
p44 <- ggplot() + geom_sf(data=map2_, fill = "#f5f8f5") +
  geom_spatraster(data=pred_r4) +   xlim(56875,387075) + ylim(6529675,6961825)  +
  scale_fill_gradientn(colours = myPalette(100), limits=c(minmax(pred_r4)[1,], minmax(pred_r4)[2,]),
                       name = "", na.value = "transparent") + 
  geom_sf(data=map33_, fill = NA, color="#03841c", linewidth=0.9) + 
  theme_bw() + ggtitle("GDE asv - romlig prediksjon")

leg4 <- ggpubr::get_legend(p44)
p4_ <- p4 + theme(legend.position = "none")
p44_ <- p44 + theme(legend.position = "none")

## ------------------------------------------------------------------------------------------------------------------------------------
#| label: plot_weight_shann_eff
lmm1 <- lm(avg_wet_weight~shannon_div+habitat_type, data=dat)
#summary(lmm1)
eff <- ggemmeans(lmm1, term="shannon_div")
pyy <- plot(eff, show_residuals=T, dot_size=1) + theme_bw() + ggtitle("") + xlab("Shannon index") + ylab("Biomasse (g)") +
  theme(axis.text.y=element_text(size=8)) +
  theme(axis.text.x=element_text(size=8))


## ----eval=F--------------------------------------------------------------------------------------------------------------------------
#| label: fig-ostland-antallarter
#| fig-pos: '!h'
#| fig-height: 7
#| fig-width: 8
#| fig-cap: Figuren viser antall arter ved østlandslokalitetene (venstre), og en romlig prediksjon av antall arter over hele østlandsområdet (høyre; under maksimal moh verdi i datasettet). Antall arter er modellert på variablene breddegrad, lengdegrad og moh. Effekten av disse på antall arter vises i de tre panelene nederst i figuren.  
# grid.arrange(p1_,p11_,leg1,
#              p1y,p1x,p1h,
#              layout_matrix = rbind(c(1,1,1, 2,2,2, 3),
#                                    c(1,1,1, 2,2,2, 3),
#                                    c(4,4, 5,5, 6,6 ,NA))  )


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-ostland-shannon
#| fig-pos: '!h'
#| fig-height: 7
#| fig-width: 7
#| fig-cap: Figuren viser Shannon indeks ved østlandslokalitetene (venstre), og en romlig prediksjon av denne indeksen over hele østlandsområdet (høyre; under maksimal moh verdi i datasettet). Shannon indeks er modellert på variablene breddegrad, lengdegrad og moh. Effekten av disse på indeksen vises i de tre panelene nederst i figuren. 
grid.arrange(p2_,p22_,leg2,
             p2y,p2x,p2h,
             layout_matrix = rbind(c(1,1,1, 2,2,2, 3),
                                   c(1,1,1, 2,2,2, 3),
                                   c(4,4, 5,5, 6,6 ,NA))  )


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-ostland-biomasse
#| fig-pos: '!h'
#| fig-height: 7
#| fig-width: 7
#| fig-cap: Figuren viser insektbiomasse ved østlandslokalitetene (venstre), og en romlig prediksjon av biomasse over hele østlandsområdet (høyre; under maksimal moh verdi i datasettet). Insektbiomasse er modellert på variablene breddegrad, lengdegrad og moh. Effekten av disse på biomassen vises i de tre panelene nederst i figuren. 
grid.arrange(p3_,p33_,leg3,
             p3y,p3x,p3h,
             layout_matrix = rbind(c(1,1,1, 2,2,2, 3),
                                   c(1,1,1, 2,2,2, 3),
                                   c(4,4, 5,5, 6,6 ,NA))  )


## ----eval=F--------------------------------------------------------------------------------------------------------------------------
#| label: fig-ostland-GDE
#| fig-pos: '!h'
#| fig-height: 7
#| fig-width: 7
#| fig-cap: Figuren viser GDE ved østlandslokalitetene (venstre), og en romlig prediksjon av GDE over hele østlandsområdet (høyre; under maksimal moh verdi i datasettet). GDE er modellert på variablene breddegrad, lengdegrad og moh. Effekten av disse på GDE vises i de tre panelene nederst i figuren. 
# grid.arrange(p4_,p44_,leg4,
#              p4y,p4x,p4h,
#              layout_matrix = rbind(c(1,1,1, 2,2,2, 3),
#                                    c(1,1,1, 2,2,2, 3),
#                                    c(4,4, 5,5, 6,6 ,NA))  )


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-ostland-shannon-v-biomasse
#| fig-pos: '!h'
#| fig-height: 5
#| fig-width: 4
#| fig-cap: Korrelasjon mellom Shannon diversitet og insektbiomasse. Punktene i figuren er residualer.  
########## residual plots ############
pyy


## ------------------------------------------------------------------------------------------------------------------------------------
leggy <- ggpubr::get_legend(p111 + guides(size = "none"))
p222 <- p222 + theme(legend.title=element_blank())
p333 <- p333 + theme(legend.title=element_blank())


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-ostland-residualer
#| fig-pos: '!h'
#| fig-height: 5
#| fig-width: 7
#| fig-cap: Residualverdier fra de romlige modellene som beskriver Shannon diversitet og insektbiomasse. De lokaliteten som har størst positivt eller negativt avvik fra modellen (dvs. høyere eller lavere observerte verdier enn forventet av modellen), vises i henholdsvis rødt og grønt. De grå verdiene samsvarer til sammenligningen med modellens prediksjoner.  
########## residual plots ############
grid.arrange(p222 + guides(color = "none") + theme(legend.position = "bottom"),
             p333 + guides(color = "none") + theme(legend.position = "bottom"),
             leggy,
             layout_matrix = rbind(c(1,1,1,1, 2,2,2,2,  3,3)))


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: true
redlisted_obs_1 <- "
SELECT
vitenskapelig_navn_id,
\"populærnavn\",
kategori_2021,
\"påvirkningsfaktorer\" as expl_fact,
obs.id_order,
obs.species_latin_fixed,
yl.year,
ls.sampling_name,
l.locality,
ST_Centroid(l.geom) as geom
FROM occurrences.observations obs LEFT JOIN
lookup.redlist_2021 redlist
ON (obs.id_genus = redlist.genus
AND obs.id_species = redlist.species),
events.identifications,
events.year_locality yl,
events.locality_sampling ls,
locations.localities l,
events.sampling_trap st
WHERE redlist.kategori_2021 IN ('DD', 'NT', 'VU', 'EN', 'CR', 'RE')
AND obs.identification_id = identifications.id
AND identifications.sampling_trap_id = st.id
AND st.locality_sampling_id = ls.id
AND ls.year_locality_id = yl.id
AND yl.locality_id = l.id
AND yl.project_short_name = 'NorIns'
AND obs.identification_confidence = 'HIGH'
"

redlisted_obs_2021 <- read_sf(con,
  query = redlisted_obs_1
)

redlisted_obs_2021 <- redlisted_obs_2021 %>%
  # mutate(expl_fact = gsub("([a-zA-Z]*)( >)(.*)", "\\1", expl_fact)) %>%
  separate(expl_fact,
    into = c(
      "expl_1",
      "expl_2",
      "expl_3",
      "expl_4",
      "expl_5"
    ),
    sep = ">"
  ) %>%
  separate(expl_3,
    into = "expl_3_main",
    sep = "_"
  ) %>%
  mutate(expl_3_main = ifelse(is.na(expl_3_main), expl_1, expl_3_main)) %>%
  mutate(expl_3_main = stringr::str_trim(expl_3_main)) %>%
  mutate(
    expl_3_main = ifelse(expl_3_main == " ", NA, expl_3_main),
    expl_3_main = ifelse(expl_3_main == "", NA, expl_3_main),
    expl_3_main = ifelse(expl_3_main == "Ukjent ", NA, expl_3_main),
    expl_3_main = ifelse(expl_3_main == "Ukjent", NA, expl_3_main)
  )

redlisted_obs_2021 <- redlisted_obs_2021 %>%
  group_by(kategori_2021) %>%
  mutate(no_spec_per_kat = n_distinct(species_latin_fixed)) %>%
  ungroup() %>%
  filter(kategori_2021 != "DD")


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: true
redlisted_obs_1_mf <- "
SELECT
vitenskapelig_navn_id,
\"populærnavn\",
kategori_2021,
\"påvirkningsfaktorer\" as expl_fact,
obs.id_order,
obs.species_latin_fixed,
yl.year,
ls.sampling_name,
l.locality,
ST_Centroid(l.geom) as geom
FROM occurrences.observations obs LEFT JOIN
lookup.redlist_2021 redlist
ON (obs.id_genus = redlist.genus
AND obs.id_species = redlist.species),
events.identifications,
events.year_locality yl,
events.locality_sampling ls,
locations.localities l,
events.sampling_trap st,
locations.traps,
lookup.trap_types
WHERE redlist.kategori_2021 IN ('DD', 'NT', 'VU', 'EN', 'CR', 'RE')
AND obs.identification_id = identifications.id
AND identifications.sampling_trap_id = st.id
AND st.locality_sampling_id = ls.id
AND ls.year_locality_id = yl.id
AND yl.locality_id = l.id
AND yl.project_short_name = 'NorIns'
AND st.trap_id = traps.id
AND trap_types.trap_model = traps.trap_model
AND trap_types.trap_type = 'Malaise'
AND obs.identification_confidence = 'HIGH'
"

redlisted_obs_2021_mf <- read_sf(con,
  query = redlisted_obs_1_mf
)

redlisted_obs_2021_mf <- redlisted_obs_2021_mf %>%
  # mutate(expl_fact = gsub("([a-zA-Z]*)( >)(.*)", "\\1", expl_fact)) %>%
  separate(expl_fact,
    into = c(
      "expl_1",
      "expl_2",
      "expl_3",
      "expl_4",
      "expl_5"
    ),
    sep = ">"
  ) %>%
  separate(expl_3,
    into = "expl_3_main",
    sep = "_"
  ) %>%
  mutate(expl_3_main = ifelse(is.na(expl_3_main), expl_1, expl_3_main)) %>%
  mutate(expl_3_main = stringr::str_trim(expl_3_main)) %>%
  mutate(
    expl_3_main = ifelse(expl_3_main == " ", NA, expl_3_main),
    expl_3_main = ifelse(expl_3_main == "", NA, expl_3_main),
    expl_3_main = ifelse(expl_3_main == "Ukjent ", NA, expl_3_main),
    expl_3_main = ifelse(expl_3_main == "Ukjent", NA, expl_3_main)
  ) %>%
  filter(kategori_2021 != "DD")


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: true

no_detected_redlist_species_per_hab_type <- redlisted_obs_2021 %>%
  st_drop_geometry() %>%
  left_join(localities,
            by = c("locality" = "locality")) |> 
  filter(kategori_2021 %in% c(
    "NT",
    "VU",
    "EN",
    "CR",
    "RE"
  )) %>%
  group_by(habitat_type) |> 
  summarise(no_spec = n_distinct(vitenskapelig_navn_id)) 

no_detected_redlist_species <- no_detected_redlist_species_per_hab_type |> 
  summarise(sum(no_spec)) |> 
  pull()

no_detected_redlist_species_skog <- no_detected_redlist_species_per_hab_type |> 
  filter(habitat_type == "Skog") |> 
  select(no_spec) |> 
  pull()

no_detected_redlist_species_semi_nat <- no_detected_redlist_species_per_hab_type |> 
  filter(habitat_type == "Semi-nat") |> 
  select(no_spec) |> 
  pull()



## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: true

redlisted_obs_2021_agg <- redlisted_obs_2021 %>%
  group_by(
    locality,
    kategori_2021
  ) %>%
  summarise(no_spec = n_distinct(species_latin_fixed))

kat_order <- tibble(
  kategori_2021 = c(
    "NT",
    "VU",
    "EN",
    "CR",
    "RE"
  ),
  kat_order = 1:5
)

no_spec_per_kat_df <- redlisted_obs_2021 %>%
  st_drop_geometry() %>%
  select(
    kategori_2021,
    no_spec_per_kat
  ) %>%
  distinct()

redlisted_obs_2021_agg <- redlisted_obs_2021_agg %>%
  left_join(kat_order,
    by = c("kategori_2021" = "kategori_2021")
  ) %>%
  ungroup() %>%
  left_join(no_spec_per_kat_df,
    by = c("kategori_2021" = "kategori_2021")
  ) %>%
  mutate(kategori_2021 = factor(paste0(kategori_2021, " (", no_spec_per_kat, " stk.)"))) %>%
  mutate(kategori_2021 = reorder(kategori_2021, kat_order))

redlisted_obs_2021_label <- redlisted_obs_2021_agg %>%
  select(locality) %>%
  distinct()

redlisted_obs_2021_agg <- redlisted_obs_2021_agg %>%
  st_jitter(redlisted_obs_2021_agg, amount = 7000)


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: true
no_locales_with_redlisted_species <- redlisted_obs_2021_agg %>%
  st_drop_geometry() %>%
  ungroup() %>%
  summarise(n_distinct(locality)) %>%
  pull()


## ------------------------------------------------------------------------------------------------------------------------------------
# redlist_colors <- c("#5a5b5d","#990013", "#d80024" , "#ff4c29", "#f36c21", "#f9a61a")
# names(redlist_colors) <- c("RE", "CR", "EN", "VU", "NT", "DD")

# redlist_colors <- c("#f9a61a", "#f36c21", "#ff4c29", "#d80024") #2022 colors
redlist_colors <- c("#ed6c26", "#e94f33", "#d80f27", "#99191a", "#5a5b5d")

# redlist_colors <- c("#f7a520", "#ed6c26", "#d80f27", "#5a5b5d") #2023 powertoys color picker tool

names(redlist_colors) <- levels(redlisted_obs_2021_agg$kategori_2021)


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-redlist-map
#| fig-cap: "Fordeling av funn av truede arter (VU, EN, CR, RE), navngitt med 'høy' sikkerhet. Sirkelstørrelsen viser antall arter funnet ved hver lokalitet, og fargen viser rødlistekategori. Posisjonene er noe justert på kartet for at de ikke skal skjule hverandre. Prikkene med 'RE' består av funn av 'Elachista quadripunctella', som er bekreftet morfologisk av Kai Berggren. Ingen CR arter ble identifisert med 'høy' sikkerhet."
#| fig-height: 8
set.seed(12345)
norge %>%
  # filter(fylke == "Innlandet" |
  #        fylke =="Vestfold og Telemark" |
  #        fylke == "Oslo" |
  #        fylke == "Viken" |
  #        fylke == "Trøndelag" |
  #        fylke == "Rogaland" |
  #        fylke == "Agder") %>%
  ggplot(.) +
  geom_sf() +
  geom_sf(
    aes(
      size = no_spec,
      color = kategori_2021
    ),
    alpha = 0.8,
    data = redlisted_obs_2021_agg %>% arrange(kategori_2021)
  ) +
  scale_color_manual(
    name = "Rødlistekategori",
    values = redlist_colors
  ) +
  scale_size(
    breaks = c(1, seq(
      from = 0,
      to = max(redlisted_obs_2021_agg$no_spec),
      by = 5
    )),
    name = "Antall arter per\nkategori og lokalitet"
  ) +
  xlab("") +
  ylab("") +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 7), #[sjekk size]
    legend.key.width = unit(0.5, "cm")
  )



## ------------------------------------------------------------------------------------------------------------------------------------

redlisted_obs_2021_to_fig <- redlisted_obs_2021 %>%
  st_drop_geometry() |> 
  left_join(localities,
            by = c("locality" = "locality")) |> 
  mutate(kategori_2021 = factor(kategori_2021, 
                                levels = c("DD",
                                           "NT",
                                           "VU",
                                           "EN",
                                           "CR",
                                           "RE")),
         region_name = factor(region_name,
                              levels = c("Sørlandet",
                                         "Østlandet",
                                         "Vestlandet",
                                         "Trøndelag",
                                         "Nord-Norge")
                              )) |>
  group_by(region_name, kategori_2021) |> 
  summarise(no_species = n_distinct(species_latin_fixed),
            no_lok = n_distinct(locality)) |> 
  arrange(region_name,
          kategori_2021) 



## ------------------------------------------------------------------------------------------------------------------------------------
region_cols <- c("Nord-Norge" =  "#004F71",
                 "Trøndelag" =  "#93328E",
                 "Sørlandet" =  "#E57200",
                 "Østlandet" =   "#008C95",
                 "Vestlandet" =  "#7A9A01") 

 
 
fig_redlist_detection <- ggplot(redlisted_obs_2021_to_fig) +
  geom_bar(aes(y = region_name,
               x = no_species,
               fill = region_name),
           stat = "identity") +
 facet_grid(cols = vars(kategori_2021),
            switch = "x") +
  coord_flip() +
  xlab("Antall arter") +
  ylab("Rødlistekategori") +
  scale_fill_manual(name = "Region",
                     values = region_cols) +
  theme(axis.text.x = element_blank())

  


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-redlist-detection
#| fig-cap: "Antall truede arter funnet i hver region, fordelt på rødlistekategori. Kun arter med 'høy' sikkerhet er tatt med." 

fig_redlist_detection


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: redlisted-obs-agg
#| cache: true

kat_order <- tibble(
  kategori_2021 = c(
    "DD",
    "NT",
    "VU",
    "EN",
    "CR",
    "RE"
  ),
  kat_order = 1:6
)

redlisted_obs_2021_less_agg <- redlisted_obs_2021 %>%
  st_drop_geometry() %>%
  dplyr:::select(
    species_latin = species_latin_fixed,
    id_order,
    popular_name = populærnavn,
    everything()
  ) %>%
  group_by(
    locality,
    kategori_2021,
    id_order,
    sampling_name,
    species_latin,
    popular_name,
    no_spec_per_kat
  ) %>%
  distinct() %>%
  group_by(
    locality,
    kategori_2021,
    id_order,
    species_latin,
    popular_name,
    no_spec_per_kat
  ) %>%
  summarise(times_detected = n()) %>%
  ungroup() %>%
  left_join(kat_order,
    by = c("kategori_2021" = "kategori_2021")
  ) %>%
  mutate(kategori_2021_append = factor(paste0(kategori_2021, " (", no_spec_per_kat, " stk.)"))) %>%
  mutate(kategori_2021_append = reorder(kategori_2021_append, kat_order)) %>%
  arrange(
    locality,
    kategori_2021_append,
    species_latin
  )

#Add regions
loc_reg <- localities %>%
  st_drop_geometry() %>%
  mutate(region_name = factor(region_name, levels = c("Nord-Norge", "Trøndelag", "Vestlandet", "Østlandet", "Sørlandet"))) %>%
  select(locality,
         year,
         region_name,
         habitat_type)

redlisted_obs_2021_less_agg <- redlisted_obs_2021_less_agg %>%
  left_join(loc_reg,
            by = c("locality" = "locality")) %>%
  arrange(region_name,
          habitat_type,
          year) %>%
  mutate(y_order = row_number(),
         y_colors = recode(region_name,
                           'Nord-Norge' = "#004F71",
                           'Trøndelag' = "#93328E",
                           'Sørlandet' = "#E57200",
                           'Østlandet' = "#008C95",
                           'Vestlandet' = "#7A9A01")) %>%
  mutate(y_label = paste("<span style = 'color: ",
                          y_colors,
                         ";'>",
                         locality,
                         "</span>", sep = ""))



## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: true

spec_no_loc <- redlisted_obs_2021_less_agg %>%
  group_by(
    kategori_2021,
    id_order,
    species_latin,
    popular_name
  ) %>%
  summarise(
    no_lok = n_distinct(locality),
    no_detected = sum(times_detected)
  ) %>%
  arrange(
    #kategori_2021,
    id_order,
    species_latin
    # desc(no_lok),
    # desc(no_detected)
  ) %>%
  dplyr::select(
    `Rødlistekat.` = kategori_2021,
    `Orden` = id_order,
    `Latinsk navn` = species_latin,
    `Norsk navn` = popular_name,
    `Antall lokaliteter` = no_lok,
    `Antall tilfeller` = no_detected
  )

spec_no_loc <- spec_no_loc %>%
  mutate(`Norsk navn` = str_to_title(`Norsk navn`))


spec_no_loc1 <- spec_no_loc %>%
  filter(`Rødlistekat.` == "NT")

spec_no_loc2 <- spec_no_loc %>%
  filter(`Rødlistekat.` == "VU" |
    `Rødlistekat.` == "EN" |
    `Rødlistekat.` == "RE" |
    `Rødlistekat.` == "CR")

## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: !expr cache_var
no_alien_q <- "SELECT obs.id_phylum, obs.id_class, obs.id_order, count(distinct obs.species_latin_fixed)
FROM occurrences.observations obs,
lookup.fremmedartslista_2023_artsdatabanken alien,
events.identifications,
events.year_locality yl,
events.locality_sampling ls,
locations.localities l,
events.sampling_trap st
WHERE obs.identification_id = identifications.id
--AND obs.species_latin_fixed = alien.vitenskapelig_navn
AND obs.latinsknavnid = alien.vitenskapelig_navn_id
AND identifications.sampling_trap_id = st.id
AND st.locality_sampling_id = ls.id
AND ls.year_locality_id = yl.id
AND yl.locality_id = l.id
AND yl.project_short_name = 'NorIns'
AND alien.risikokategori_2023 IN ('SE', 'HI', 'PH', 'LO')
AND obs.identification_confidence = 'HIGH'
GROUP BY id_phylum, id_class, id_order
ORDER BY id_phylum, id_class, id_order

"

no_alien <- dbGetQuery(
  con,
  no_alien_q
) %>%
  # arrange(id_class, id_order, species_latin_fixed) %>%
  mutate(count = as.integer(count)) %>%
  select(everything(),
    no_alien = count
  )


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: !expr cache_var
no_fennoscand_q <- "
SELECT obs.id_phylum, obs.id_class, obs.id_order, count(distinct obs.species_latin_fixed)
FROM occurrences.observations obs,
lookup.horizon_scan_2025 h,
events.identifications,
events.year_locality yl,
events.locality_sampling ls,
locations.localities l,
events.sampling_trap st
WHERE obs.identification_id = identifications.id
AND obs.species_latin_fixed = h.species_latin_fixed
AND identifications.sampling_trap_id = st.id
AND st.locality_sampling_id = ls.id
AND ls.year_locality_id = yl.id
AND yl.locality_id = l.id
AND yl.project_short_name = 'NorIns'
AND obs.identification_confidence = 'HIGH'
AND h.fennoscandinavian IS TRUE
GROUP BY id_phylum, id_class, id_order
ORDER BY id_phylum, id_class, id_order
"

no_fennoscand <- dbGetQuery(
  con,
  no_fennoscand_q
) %>%
  mutate(count = as.integer(count)) %>%
  select(everything(),
    no_fennoscand_alien = count
  )


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: !expr cache_var
no_pot_alien_q <- "
SELECT obs.id_phylum, obs.id_class, obs.id_order, count(distinct obs.species_latin_fixed)
FROM occurrences.observations obs,
lookup.horizon_scan_2025 alien,
events.identifications,
events.year_locality yl,
events.locality_sampling ls,
locations.localities l,
events.sampling_trap st
WHERE obs.identification_id = identifications.id
AND obs.species_latin_fixed = alien.species_latin_fixed
AND identifications.sampling_trap_id = st.id
AND st.locality_sampling_id = ls.id
AND ls.year_locality_id = yl.id
AND yl.locality_id = l.id
AND yl.project_short_name = 'NorIns'
AND obs.identification_confidence = 'HIGH'
AND alien.potential_alien_species IS TRUE
GROUP BY id_phylum, id_class, id_order
ORDER BY id_phylum, id_class, id_order
"

no_pot_alien <- dbGetQuery(
  con,
  no_pot_alien_q
) %>%
  mutate(count = as.integer(count)) %>%
  select(everything(),
    no_pot_alien = count
  ) |> 
  as_tibble()


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: !expr cache_var
no_all_species_q <- "
SELECT obs.id_phylum, obs.id_class, obs.id_order, count(distinct obs.species_latin_fixed)
FROM occurrences.observations obs,
events.identifications,
events.year_locality yl,
events.locality_sampling ls,
locations.localities l,
events.sampling_trap st,
locations.traps
WHERE obs.identification_id = identifications.id
AND identifications.sampling_trap_id = st.id
AND st.locality_sampling_id = ls.id
AND ls.year_locality_id = yl.id
AND yl.locality_id = l.id
AND st.trap_id = traps.id
AND yl.project_short_name = 'NorIns'
AND ((yl.year = 2020 AND traps.trap_short_name IN ('MF1', 'MF3', 'VF1',  'VF3'))
OR yl.year != 2020)
--AND obs.identification_confidence = 'HIGH'
GROUP BY id_phylum, id_class, id_order
ORDER BY id_phylum, id_class, id_order
"

no_all_species <- dbGetQuery(
  con,
  no_all_species_q
)

# This is a hacky 'fix' to get the total species numbers to align. Should be worked out properly!
no_all_species <- no_all_species %>%
  mutate(count = ifelse(id_class == "Insecta" & is.na(id_order), as.integer(count) - 2, as.integer(count))) %>%
  mutate(count = as.integer(count)) %>%
  select(everything(),
    no_all_species = count
  )


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: !expr cache_var
no_all_species_to_species_q <- "
SELECT obs.id_phylum, obs.id_class, obs.id_order, count(distinct obs.species_latin_fixed)
FROM occurrences.observations obs,
events.identifications,
events.year_locality yl,
events.locality_sampling ls,
locations.localities l,
events.sampling_trap st,
locations.traps
WHERE obs.identification_id = identifications.id
AND identifications.sampling_trap_id = st.id
AND st.locality_sampling_id = ls.id
AND ls.year_locality_id = yl.id
AND yl.locality_id = l.id
AND st.trap_id = traps.id
AND yl.project_short_name = 'NorIns'
AND ((yl.year = 2020 AND traps.trap_short_name IN ('MF1', 'MF3', 'VF1',  'VF3'))
OR yl.year != 2020)
AND obs.taxonomic_level = 'species'
AND obs.identification_confidence = 'HIGH'
GROUP BY id_phylum, id_class, id_order
ORDER BY id_phylum, id_class, id_order
"

no_all_species_to_species <- dbGetQuery(
  con,
  no_all_species_to_species_q
) %>%
  mutate(count = as.integer(count)) %>%
  select(everything(),
    no_all_species_to_species = count
  )


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: !expr cache_var
no_all_native_species_to_species_q <- "
SELECT obs.id_phylum, obs.id_class, obs.id_order, count(distinct obs.species_latin_fixed)
FROM occurrences.observations obs,
events.identifications,
events.year_locality yl,
events.locality_sampling ls,
locations.localities l,
events.sampling_trap st,
locations.traps
WHERE obs.identification_id = identifications.id
AND identifications.sampling_trap_id = st.id
AND st.locality_sampling_id = ls.id
AND ls.year_locality_id = yl.id
AND yl.locality_id = l.id
AND yl.project_short_name = 'NorIns'
AND st.trap_id = traps.id
AND ((yl.year = 2020 AND traps.trap_short_name IN ('MF1', 'MF3', 'VF1',  'VF3'))
OR yl.year != 2020)
AND obs.taxonomic_level = 'species'
AND obs.native IS TRUE
AND obs.identification_confidence = 'HIGH'
GROUP BY id_phylum, id_class, id_order
ORDER BY id_phylum, id_class, id_order
"


no_all_native_species_to_species <- dbGetQuery(
  con,
  no_all_native_species_to_species_q
) %>%
  mutate(count = as.integer(count)) %>%
  select(everything(),
    no_all_native_species_to_species = count
  )


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: !expr cache_var

no_all_aliens <- no_fennoscand %>%
  full_join(no_pot_alien,
    by = c(
      "id_phylum" = "id_phylum",
      "id_class" = "id_class",
      "id_order" = "id_order"
    )
  ) %>%
  full_join(no_alien,
    by = c(
      "id_phylum" = "id_phylum",
      "id_class" = "id_class",
      "id_order" = "id_order"
    )
  ) %>%
  full_join(no_all_species,
    by = c(
      "id_phylum" = "id_phylum",
      "id_class" = "id_class",
      "id_order" = "id_order"
    )
  ) %>%
  full_join(no_all_species_to_species,
    by = c(
      "id_phylum" = "id_phylum",
      "id_class" = "id_class",
      "id_order" = "id_order"
    )
  ) %>%
  full_join(no_all_native_species_to_species,
    by = c(
      "id_phylum" = "id_phylum",
      "id_class" = "id_class",
      "id_order" = "id_order"
    )
  ) %>%
  arrange(
    id_phylum,
    id_class,
    id_order
  ) %>%
  replace(is.na(.), 0) %>%
  mutate(id_order = ifelse(id_order == 0, "Ukjent", id_order),
         id_phylum = ifelse(id_phylum == 0, "Ukjent", id_phylum),
         id_class = ifelse(id_class == 0, "Ukjent", id_class))


no_all_aliens <- no_all_aliens %>%
  as_tibble() %>%
  mutate(
    no_all_species = as.integer(no_all_species),
    no_all_species_to_species = as.integer(no_all_species_to_species),
    no_all_native_species_to_species = as.integer(no_all_native_species_to_species),
    no_alien = as.integer(no_alien),
    no_fennoscand_alien = as.integer(no_fennoscand_alien),
    no_pot_alien = as.integer(no_pot_alien)
  ) %>%
  mutate(id_order = ifelse(id_order == "Plecoptera_Insecta", "Plecoptera", id_order)) %>%
  select(
    Fylum = "id_phylum",
    Klasse = "id_class",
    Orden = "id_order",
    "Antall takson" = no_all_species,
    "Takson bestemt til art" = no_all_species_to_species,
    "Finnes i Norge" = no_all_native_species_to_species,
    "Fremmede arter" = no_alien,
    "Fennoskandiske arter" = no_fennoscand_alien,
    "Potensielt fremmede arter" = no_pot_alien
  )

no_all_aliens <- no_all_aliens %>%
  add_row(
    Fylum = "",
    Klasse = "",
    Orden = "Sum",
    `Antall takson` = sum(.$`Antall takson`),
    `Takson bestemt til art` = sum(.$`Takson bestemt til art`),
    `Finnes i Norge` = sum(.$`Finnes i Norge`),
    `Fremmede arter` = sum(.$`Fremmede arter`),
    `Fennoskandiske arter` = sum(.$`Fennoskandiske arter`),
    `Potensielt fremmede arter` = sum(.$`Potensielt fremmede arter`)
  )


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: tbl-no-species-order
#| results: asis
#| tbl-cap: "Antall arter og høyere taksa påvist med DNA-metastrekkoding i alle felleprøver (2020 - 2024) sortert på orden. Antallet av de observerte arter som finnes i Norge (registrert i artsnavnebasen), som er kjente fremmede arter for Norge (registrert i fremmedartslista), som har tilstedeværelse i Fennoskandia (i GBIF), og som er potensielt fremmede arter (se figur 8) gis i separate kolonner. Kun artsidentifiseringer med 'høy' konfidens er tatt med (unntatt 'Antall takson')."

xtab <- xtable(no_all_aliens)
align(xtab) <- "p{1.5cm}p{1.2cm}p{1.2cm}p{2.2cm}p{1cm}p{1cm}p{1cm}p{1cm}p{1.7cm}p{1cm}"

if (names(rmarkdown::metadata$format) == "ninareport-pdf") {
  print(xtab, sanitize.text.function = function(x) {
    x
  }, size = "scriptsize")
} else {
  kable(xtab)
}


## ------------------------------------------------------------------------------------------------------------------------------------
pot_alien_spec <- tbl(
  con,
  Id(
    schema = "lookup",
    table = "horizon_scan_2025"
  )
) %>%
  filter(potential_alien_species) %>%
  select(species_latin_fixed) %>%
  pull()

fennoscand_spec <- tbl(
  con,
  Id(
    schema = "lookup",
    table = "horizon_scan_2025"
  )
) %>%
  filter(fennoscandinavian) %>%
  select(species_latin_fixed) %>%
  pull()

alien_spec <- tbl(
  con,
  Id(
    schema = "lookup",
    table = "fremmedartslista_2023_artsdatabanken"
  )
) %>%
  filter(ekspertkomite == "Terrestriske invertebrater") %>%
  select(vitenskapelig_navn) %>%
  pull()


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: true
fennoskand_obs_q <- "
SELECT DISTINCT ON(species_latin_fixed, year, sampling_name, locality, geom)
obs.species_latin_fixed,
yl.year,
ls.sampling_name,
l.locality,
ST_Centroid(l.geom) as geom
FROM occurrences.observations obs,
lookup.horizon_scan_2025 fennoscand,
events.identifications,
events.year_locality yl,
events.locality_sampling ls,
locations.localities l,
events.sampling_trap st
WHERE obs.identification_id = identifications.id
AND obs.species_latin_fixed = fennoscand.species_latin_fixed
AND identifications.sampling_trap_id = st.id
AND st.locality_sampling_id = ls.id
AND ls.year_locality_id = yl.id
AND yl.locality_id = l.id
AND yl.project_short_name = 'NorIns'
AND obs.identification_confidence = 'HIGH'
AND fennoscand.fennoscandinavian IS TRUE
"

fennoskand_obs <- read_sf(con,
  query = fennoskand_obs_q
) %>%
  mutate(alien_kat = "Fennoskandisk forek.")


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: true
pot_alien_obs_q <- "
SELECT DISTINCT ON (obs.species_latin_fixed, year, sampling_name, locality, geom)
obs.species_latin_fixed,
yl.year,
ls.sampling_name,
l.locality,
ST_Centroid(l.geom) as geom
FROM occurrences.observations obs,
lookup.horizon_scan_2025 pot_alien,
events.identifications,
events.year_locality yl,
events.locality_sampling ls,
locations.localities l,
events.sampling_trap st
WHERE obs.identification_id = identifications.id
AND obs.species_latin_fixed = pot_alien.species_latin_fixed
AND identifications.sampling_trap_id = st.id
AND st.locality_sampling_id = ls.id
AND ls.year_locality_id = yl.id
AND yl.locality_id = l.id
AND yl.project_short_name = 'NorIns'
AND obs.identification_confidence = 'HIGH'
AND pot_alien.potential_alien_species IS TRUE
"

pot_alien_obs <- read_sf(con,
  query = pot_alien_obs_q
) %>%
  mutate(alien_kat = "Potensielt fremmede arter")


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: true
alien_obs_q <- "
SELECT
obs.species_latin_fixed,
yl.year,
ls.sampling_name,
l.locality,
ST_Centroid(l.geom) as geom
FROM occurrences.observations obs,
lookup.horizon_scan_2025 alien,
events.identifications,
events.year_locality yl,
events.locality_sampling ls,
locations.localities l,
events.sampling_trap st
WHERE obs.identification_id = identifications.id
--AND obs.species_latin_fixed = alien.vitenskapelig_navn
AND obs.species_latin_fixed = alien.species_latin_fixed
AND identifications.sampling_trap_id = st.id
AND st.locality_sampling_id = ls.id
AND ls.year_locality_id = yl.id
AND yl.locality_id = l.id
AND yl.project_short_name = 'NorIns'
AND alien.risikokategori_2023 IN ('SE', 'HI', 'PH', 'LO')
AND obs.identification_confidence = 'HIGH'
"

alien_obs <- read_sf(con,
  query = alien_obs_q
) %>%
  mutate(alien_kat = "På fremmedartslista")


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: true

all_alien_obs <- fennoskand_obs %>%
  rbind(pot_alien_obs) %>%
  rbind(alien_obs)


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: true
kat_order <- tibble(
  alien_kat = c(
    "Fennoskandisk forek.",
    "Potensielt fremmede arter",
    "På fremmedartslista"
  ),
  kat_order = 1:3
)

all_alien_obs_agg <- all_alien_obs %>%
  group_by(alien_kat) %>%
  mutate(no_spec_per_kat = n_distinct(species_latin_fixed)) %>%
  group_by(
    locality,
    alien_kat,
    no_spec_per_kat
  ) %>%
  summarise(no_spec = n_distinct(species_latin_fixed)) %>%
  ungroup() %>%
  mutate(alien_kat_append = factor(paste0(alien_kat, " (", no_spec_per_kat, " stk.)"))) %>%
  left_join(kat_order,
    by = c("alien_kat" = "alien_kat")
  ) %>%
  mutate(alien_kat_append = reorder(alien_kat_append, kat_order))

all_alien_obs_agg <- all_alien_obs_agg %>%
  st_jitter(all_alien_obs_agg, amount = 7000)


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-non-native-species-kat
#| fig-cap: "Insektarter uten tidligere kjent forekomst i Norge, basert på offentlige digitale kilder, registrert i overvåkingen i 2020-2024. Fennoskandisk forekomst baseres på registreringer i GBIF. Potensielt fremmede arter er arter som ikke er registrert observert i Norge eller Fennoskandia. Observasjonene er basert på DNA-metastrekkoding og er ikke morfologisk verifisert. Arter vurdert til NR (Ikke risikovurdert) på fremmedartslista er ikke tatt med. Posisjonene er spredt ut noe på kartet for at punktene ikke skal overlappe."
#| fig-height: 8
set.seed(12345)
norge %>%
  ggplot(.) +
  geom_sf() +
  geom_sf(
    aes(
      size = no_spec,
      color = alien_kat_append
    ),
    alpha = 0.8,
    data = all_alien_obs_agg
  ) +
  scale_color_nina(
    name = "Kategori",
    palette = "contrast"
  ) +
  scale_size(
    breaks = c(
      1, seq(
        from = 10,
        to = max(all_alien_obs_agg$no_spec),
        by = 20
      ),
      max(all_alien_obs_agg$no_spec)
    ),
    name = "Antall arter per\nkategori og lokalitet"
  ) +
  xlab("") +
  ylab("") +
  guides(color = guide_legend(override.aes = list(size = 5), nrow = 2, byrow = T)) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 7),
    legend.key.width = unit(0.5, "cm")
  )


## ------------------------------------------------------------------------------------------------------------------------------------

all_alien_obs_to_fig <- all_alien_obs %>%
  st_drop_geometry() |> 
  left_join(localities,
            by = c("locality" = "locality")) |> 
  mutate(alien_kat = factor(alien_kat, 
                                levels = c("Fennoskandisk forek.",
                                           "Potensielt fremmede arter",
                                           "På fremmedartslista")
                            ),
         region_name = factor(region_name,
                              levels = c("Sørlandet",
                                         "Østlandet",
                                         "Vestlandet",
                                         "Trøndelag",
                                         "Nord-Norge")
                              )) |>
  group_by(region_name, alien_kat) |> 
  summarise(no_species = n_distinct(species_latin_fixed),
            no_lok = n_distinct(locality)) |> 
  arrange(region_name,
          alien_kat) 



## ------------------------------------------------------------------------------------------------------------------------------------
region_cols <- c("Nord-Norge" =  "#004F71",
                 "Trøndelag" =  "#93328E",
                 "Sørlandet" =  "#E57200",
                 "Østlandet" =   "#008C95",
                 "Vestlandet" =  "#7A9A01") 

fig_alien_detection <- ggplot(all_alien_obs_to_fig) +
  geom_bar(aes(y = region_name,
               x = no_species,
               fill = region_name),
           stat = "identity") +
 facet_grid(cols = vars(alien_kat),
            switch = "x",
            ) +
  coord_flip() +
  xlab("Antall arter") +
  ylab("Kategori") +
  scale_fill_manual(name = "Region",
                     values = region_cols) +
  theme(axis.text.x = element_blank(),
        strip.text.x = element_text(size = 7),
        axis.ticks.x = element_blank())

  


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-alien-detection
#| fig-cap: Antall potensielt fremmede arter funnet i hver region, fordelt på kategori.

fig_alien_detection



## ------------------------------------------------------------------------------------------------------------------------------------
#| label: alien-obs-agg
#| cache: true

kat_order <- tibble(
  alien_kat = c(
    "Fennoskandisk forek.",
    "Ikke fennoskandisk forek.",
    "På fremmedartslista"
  ),
  kat_order = 1:3
)

all_alien_obs_less_agg <- all_alien_obs %>%
  st_drop_geometry() %>%
  group_by(alien_kat) %>%
  mutate(no_spec_per_kat = n_distinct(species_latin_fixed)) %>%
  select(
    species_latin = species_latin_fixed,
    everything()
  ) %>%
  group_by(
    locality,
    alien_kat,
    sampling_name,
    species_latin
  ) %>%
  distinct() %>%
  group_by(
    locality,
    alien_kat,
    no_spec_per_kat,
    species_latin
  ) %>%
  summarise(times_detected = n()) %>%
  mutate(alien_kat_append = factor(paste0(alien_kat, " (", no_spec_per_kat, " stk.)"))) %>%
  left_join(kat_order,
    by = c("alien_kat" = "alien_kat")
  ) %>%
  mutate(alien_kat_append = reorder(alien_kat_append, kat_order)) %>%
  arrange(
    locality,
    alien_kat,
    species_latin
  ) %>%
  left_join(loc_reg,
            by = c("locality" = "locality")) %>%
  ungroup() %>%
  arrange(region_name,
          habitat_type,
          year) %>%
  mutate(y_order = row_number(),
         y_colors = recode(region_name,
                           'Nord-Norge' = "#004F71",
                           'Trøndelag' = "#93328E",
                           'Sørlandet' = "#E57200",
                           'Østlandet' = "#008C95",
                           'Vestlandet' = "#7A9A01")) %>%
  mutate(y_label = paste("<span style = 'color: ",
                          y_colors,
                         ";'>",
                         locality,
                         "</span>", sep = ""))



## ------------------------------------------------------------------------------------------------------------------------------------
#| label: load_unkn_biodiv
#| cache: !expr cache_var

figdat_raw <- read.table("../../../../store_filer_for_rapport/unknown_biodiversity_figure_data.txt",
  header = T,
  sep = "\t",
  fileEncoding = "latin1"
)


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: mk_unkn_biodiv_fig
#| cache: !expr cache_var

figdat <- figdat_raw |> 
  mutate(taxa = ifelse(taxa == "andre Insecta", "Øvrige", taxa))

#collate all known norwegian species

figdat <- figdat |> 
  mutate(variable = ifelse(xcat == "Norge" &
                variable == "antall_norskarter_NasIns", "tabort_antall_norskarter_NasIns", variable))


norwegian_NorIns <- figdat |>
  filter(xcat == 'Norge',
         variable == 'tabort_antall_norskarter_NasIns') |> 
  select(taxa,
         value,
         xcat)

norwegian_not_NorIns <- figdat |>
  filter(xcat == 'Norge',
         variable == 'antall_norsk_arter_ikke_i_NasIns') |> 
  select(taxa,
         value,
         xcat)
  
norske_arter <- norwegian_NorIns |> 
  inner_join(norwegian_not_NorIns,
             by = c("taxa" = "taxa",
                    "xcat" = "xcat")) |> 
  mutate(value = value.x + value.y,
         variable = 'all_known_norske_arter') |>
  select(taxa,
         variable,
         value,
         xcat)
  
figdat <- figdat |> 
  rbind(norske_arter) |> 
  filter(variable != "tabort_antall_norskarter_NasIns")



figdat_to_plot <- figdat[!figdat$variable %in% c(
    "antall_norsk_arter_ikke_i_NasIns",
    "antall_fremmede_NasIns",
    "Antatt.antall.arter.i.Norge",
    "antall_norsk_arter",
    "antall_arter_NasIns",
    "antall_norsk_taxa_i_NorInvert2.0"
  ) &
    figdat$taxa %in% c(
      "Øvrige",
      "Coleoptera",
      "Lepidoptera",
      "Hymenoptera",
      "Hemiptera",
      "Diptera",
      "Arachnida"
    ), ]

figdat_to_plot$variable <- factor(figdat_to_plot$variable,
  levels = rev(
    c(
      "antall_norskarter_NasIns",
      "all_known_norske_arter",
      "estimated_overlooked",
      "antall_fennoskandiske_NasIns",
      "darkmatter"
    )
  )
)

figdat$taxa <- factor(figdat$taxa,
  levels = c(
    "Arachnida",
    "Chilopoda",
    "Diplopoda",
    "Entognatha",
    "Øvrige",
    "Coleoptera",
    "Diptera",
    "Hemiptera",
    "Hymenoptera",
    "Lepidoptera"
  )
)

darkmatter <- ggplot(figdat_to_plot,
  aes(
    x = xcat,
    y = value,
    fill = variable
  ),
  pattern = "stripe"
) +
  geom_col_pattern(
    aes(
      pattern = variable,
      fill = variable
    ),
    pattern_fill = "white",
    pattern_color = "white",
    pattern_spacing = .2,
    pattern_key_scale_factor = 0.25,
    color = "black"
  ) +
  facet_grid(~taxa,
    switch = "x"
  ) +
  scale_pattern_manual(
    values = c(
      "darkmatter" = "none",
      "antall_fennoskandiske_NasIns" = "none",
      "estimated_overlooked" = "stripe",
      "all_known_norske_arter" = "none",
      "antall_norskarter_NasIns" = "none"
    ),
    labels = c(
      "Observert\nukjent biomangfold",
      "Observert\nfennoskandiske",
      "Estimert ukjent\nbiomangfold",
      "Kjente norske arter",
      "Observerte\nnorske arter"
    ),
    name = "Kategori"
  ) +
 # theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
    strip.text = element_text(face = "bold", size = 5),
    legend.text = element_text(size = 6)
    ) +
  scale_fill_manual(
    values = c(
      "grey10",
      "#E57200",
      "grey50",
      "#7FC5CA",
      "#008C95"
    ),
    labels = c(
      "Observert\nukjent biomangfold",
      "Observert\nfennoskandiske",
      "Estimert ukjent\nbiomangfold",
      "Kjente norske arter",
      "Observerte\nnorske arter"
    ),
    name = "Kategori"
  ) +
  ylab("Antall arter")


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-darkmatter
#| cache: !expr cache_var
#| fig-cap: "Dekning i NorIns av Norges biomangfold av insekter og edderkoppdyr. For hver taksonomisk gruppe vises det totale estimerte biomangfoldet for Norge i første kolonne. Denne kolonnen er delt opp i de artene som er blitt oppdaget i NorIns (mørkeblått), de som ikke er blitt observert i NorIns (lysblått), og de som aldri er blitt observert, men antas å finnes i Norge (i grått). Dette estimerte ukjente biomangfoldet for Norge er tatt fra @elven_kunnskapsstatus_2021. Andre kolonne for hver artsgruppe viser det som er blitt observert i NorIns. Denne viser kjente arter fra tidligere (mørkeblått, samme som i venstre kolonne), arter som er kjent fra naboland (oransje), og observert ukjent biomangfold (svart). Observert ukjent biomangfold består av arter som enten er nye for Norge eller for vitenskapen, samt automatisert grupperte arter som ikke har kunnet kobles til et artsnavn med identifiseringsmetodikken brukt i insektovervåkingsprogrammet."

darkmatter



## ------------------------------------------------------------------------------------------------------------------------------------
#| label: tbl-div-tab
#| results: asis
#| tbl-cap: "Antall observerte arter leddyr fra overvåkingen i 2020-2024. Merk at artsantallet delvis er basert på en automatisk oppdeling av genetiske strenger, og skal ses som et estimat. Test-feller med dobbel fangstperiode fra 2020 er ikke tatt med, derfor skiller seg tallet noe fra @tbl-no-species-order."

to_tab <- region_habitat_div %>%
  arrange(
    År,
    Region,
    `Habitat-type`
  )

if (names(rmarkdown::metadata$format) == "ninareport-pdf") {
  print(xtable(to_tab))
} else {
  kable(to_tab)
}


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-richness-locality
#| fig-cap: "Totalt observert artsantall per lokalitet. Merk at DNA-metastrekkodingsrutinen ble oppdatert i 2021 og at dette kan forklare det lavere artsantallet i 2020."

year_locality_biodiv %>%
  ggplot(
    .,
    aes(
      x = year,
      y = no_species
    )
  ) +
  geom_boxplot(
    aes(fill = region_name)
    #        stat = "identity"
  ) +
  scale_fill_nina(
    name = "År",
    palette = "main"
  ) +
  ylab("Totalt antall observerte arter per lokalitet") +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1,
    size = 8
  )) +
  facet_wrap("habitat_type",
    nrow = 2
  )




## ------------------------------------------------------------------------------------------------------------------------------------
load(file = "out/ost_iNEXT.rda")


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-div-est
#| fig-height: 7
#| fig-cap: "Observert og estimert totalt artsantall (Chao) i malaisefeller i de undersøkte regionene. Merk at ekstrapoleringen forutsetter en perfekt deteksjonsevne. Usikkerhetene rundt estimatene er så små at de dekkes av linjene."

lines <- c(
  "Rarefaction" = "solid",
  "Extrapolation" = "dotted"
)



p <- ggiNEXT(ost_iNEXT) 

p + 
  ylab("Artsantall") +
  xlab("Antall lokaliteter") +
  theme(
    legend.text = element_text(size = 8),
    legend.key.width = unit(1, "cm"),
    axis.title = element_text(size = 14)
  ) +
  scale_color_manual(values = c("Semi-nat_Nord"= "#004F71", 
                                "Semi-nat_Øst" = "#008C95", 
                                "Semi-nat_Sør" = "#E57200",
                                "Semi-nat_Trond" = "#93328E", 
                                "Skog_Øst" = "#008C95", 
                                "Semi-nat_Vest" = "#7A9A01")) +
  guides(
    color = guide_legend("",
                         nrow = 2, byrow = TRUE),
    shape = guide_legend("", nrow = 2, byrow = TRUE),
    fill = guide_legend("", nrow = 2, byrow = TRUE)
  ) +
  scale_linetype_manual(
    values = lines,
    name = "Method",
    labels = c("Interpolert", "Ekstrapolert")
  )


## ------------------------------------------------------------------------------------------------------------------------------------
#| cache: true
spec_occ <- tbl(
  con,
  Id(
    schema = "views",
    table = "species_occurrences"
  )
) %>%
  filter(project_short_name == "NorIns")

spec_no_occ <- spec_occ %>%
  filter(!grepl("MF2", sample_name)) %>%
  group_by(species_latin_fixed) %>%
  summarise(no_occ = as.integer(n())) %>%
  arrange(desc(no_occ)) %>%
  collect() %>%
  mutate(art_no = row_number()) %>%
  group_by(species_latin_fixed)

no_samples <- spec_occ %>%
  filter(!grepl("MF2", sample_name)) %>%
  summarise(no_samples = n_distinct(sample_name)) %>%
  pull()

no_singletons <- spec_no_occ |> 
  filter(no_occ < 2) |> 
  nrow() 

no_doubletons <- spec_no_occ |> 
  filter(no_occ == 2) |> 
  nrow() 

no_less_than_10 <- spec_no_occ |> 
  filter(no_occ < 10) |> 
  nrow() 

more_than_100 <- spec_no_occ |> 
   filter(no_occ > 100) |> 
   nrow() 


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-spec-occ
#| fig-cap: "Rangering av i hvor mange felleprøver de ulike insektartene har blitt funnet."
#| fig-width: 9
spec_no_occ %>%
  select(no_occ, art_no) %>%
  ggplot(
    .,
    aes(x = art_no, y = no_occ)
  ) +
  geom_bar(
    stat = "identity",
    aes(fill = no_occ),
    width = 100
  ) +
  ylab(paste0("Observert antall ganger av ", no_samples, " mulige")) +
  xlab("Art") +
  scale_fill_nina(discrete = FALSE) +
  scale_x_continuous(
    position = "bottom",
    breaks = c(1, 5000, 10000, 15000, 20000)
  ) +
  theme(legend.position = "")
#+
# coord_flip() +
# scale_y_reverse() +
# theme(axis.title.x = element_text(angle = 180),
#      axis.title.y = element_text(angle = 180))
#+
# geom_point(aes(x = 1, y = max(no_occ)), col = ninaColors()[3])





## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-haplotype-kart
#| fig-cap: "Frekvenser av haplotyper (varianter av ASV/DNA-kopier) ved ulike lokaliteter som viser genetisk variasjon innen hver enkeltart. Som et oversiktlig eksempel vises de tre artene som er fanget i flest tilfeller innenfor hver av ordenene sommerfugler (Lepidoptera), biller (Coleoptera), tovinger (Diptera) og vepser (Hymenoptera). Hver farge representerer en unik genetisk variant. Størrelsen på sirklene er logaritmisk skalert etter mengden DNA av arten i hver lokalitet. Besøk programmets nettside for å se valgfri art i større kart."
#| out-width: '100%'
#| out-height: '90%'

include_graphics("figures/asv_composite.png")



## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-poll-haplotype-kart
#| fig-cap: "Frekvenser av haplotyper (varianter av ASV/DNA-kopier) for et utvalg av pollinatorer-familier, langtungebier (Apidae), nymfevinger (Nymphalidae), glansvinger (Lycaenidae), og blomsterfluer (Syrphidae). Hver farge representerer en unik genetisk variant. Størrelsen på sirklene er logaritmisk skalert etter mengden DNA av arten ved hver lokalitet. Besøk programmets nettside for å se valgfri art i større kart."
#| out-width: '100%'
#| out-height: '90%'

include_graphics("figures/asv_poll_composite.png")


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: fig-frysekontainer
#| fig-cap: "Frysekontainer i innflyttingsfase. Kontaineren beregnes har kapasitet for ca. 10 års drift av programmet."
#| out-height: 5in
#| out-width: 4in

include_graphics("figures/Frysekontainer_2_rot_small.jpg")


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: join_for_app1_tab
yl_catch_tab <- year_locality_biodiv %>%
  #mutate(year = as.integer(as.character(year))) %>%
  left_join(year_locality_biomass) %>%
  arrange(year, region_name, habitat_type, locality) %>%
  select(
    År = year,
    Lokalitet = locality,
    Økosystemtype = habitat_type,
    Region = region_name,
    `Antall arter` = no_species,
    `Shannon div.` = shannon_div,
    `Våtvekt (middelv.)` = avg_wet_weight
  )


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: ylcatchtab
#| results: asis

if (names(rmarkdown::metadata$format) == "ninareport-pdf") {
  print(xtable(yl_catch_tab),
    size = "footnotesize",
    floating = F,
    floating.environment = "table*",
    tabular.environment = "longtable"
  )
} else if (names(rmarkdown::metadata$format) == "ninareport-html") {
  DT::datatable(yl_catch_tab,
    style = "bootstrap",
    options = list(
      pageLength = 10,
      autoWidth = TRUE,
      columnDefs = list(list(
        width = "100px",
        targets = 1
      ))
    )
  )
} else {
  kable(yl_catch_tab)
}


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: table_redlist_spec
#| results: asis
#| eval: true


if (names(rmarkdown::metadata$format) == "ninareport-pdf") {
  print(
    xtable(spec_no_loc1,
      align = rep("l", 7)
    ),
    size = "footnotesize",
    sanitize.colnames.function = function(x) {
      x
    },
    # table.placement = "p",
    floating.environment = "table*",
    tabular.environment = "longtable",
    floating = F
  )
} else if (names(rmarkdown::metadata$format) == "ninareport-html") {
  DT::datatable(spec_no_loc1,
    style = "bootstrap",
    options = list(
      pageLength = 10,
      autoWidth = TRUE,
      columnDefs = list(list(
        width = "100px",
        targets = 1
      ))
    )
  )
} else {
  kable(spec_no_loc1)
}


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: table_end_spec
#| results: asis
#| eval: true


if (names(rmarkdown::metadata$format) == "ninareport-pdf") {
  print(
    xtable(spec_no_loc2,
      align = rep("l", 7)
    ),
    size = "footnotesize",
    sanitize.colnames.function = function(x) {
      x
    },
    # table.placement = "p",
    floating.environment = "table*",
    tabular.environment = "longtable",
    floating = F
  )
} else if (names(rmarkdown::metadata$format) == "ninareport-html") {
  DT::datatable(spec_no_loc2,
    style = "bootstrap",
    options = list(
      pageLength = 10,
      autoWidth = TRUE,
      columnDefs = list(list(
        width = "100px",
        targets = 1
      ))
    )
  )
} else {
  kable(spec_no_loc2)
}


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: get_obs_for_appendix
#| cache: !expr cache_var
obs <- tbl(
  con,
  Id(
    schema = "occurrences",
    table = "observations"
  )
) %>%
  filter(identification_confidence == "HIGH") |> 
  mutate(identification_id = as.character(identification_id))

project_ids_q <- "
SELECT project_short_name, i.id::text as identification_id
FROM events.identifications i,
events.sampling_trap st,
events.locality_sampling ls,
events.year_locality yl
WHERE i.sampling_trap_id = st.id
AND st.locality_sampling_id = ls.id
AND ls.year_locality_id = yl.id
"

project_ids <- dbGetQuery(con,
                         project_ids_q) |> 
  as_tibble()

obs <- obs |> 
  inner_join(project_ids,
             by = c("identification_id" = "identification_id"),
             copy = TRUE) |> 
  filter(project_short_name == 'NorIns')


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: get_pot_alien_new
#| cache: !expr cache_var

horizon_scan_2025 <- tbl(
  con,
  Id(
    schema = "lookup",
    table = "horizon_scan_2025"
  )
)

pot_alien <- horizon_scan_2025 |> 
  filter(potential_alien_species) |> 
  mutate(
    saul_boolean = ifelse(is.na(saul_boolean), FALSE, saul_boolean),
    daisie_boolean = ifelse(is.na(daisie_boolean), FALSE, daisie_boolean)
  ) |> 
  inner_join(obs,
             by = c("species_latin_fixed" = "species_latin_fixed")) |> 
  select(
    Order = "id_order",
    Artsnavn = species_latin_fixed,
    "Listet i Saul" = saul_boolean,
    "Listet i Daisie" = daisie_boolean
  ) |> 
  arrange(Artsnavn) |> 
  distinct() |> 
  collect()


fennoscand_spec <- horizon_scan_2025 |> 
  filter(fennoscandinavian) |> 
  mutate(
    saul_boolean = ifelse(is.na(saul_boolean), FALSE, saul_boolean),
    daisie_boolean = ifelse(is.na(daisie_boolean), FALSE, daisie_boolean)
  ) |> 
  inner_join(obs,
             by = c("species_latin_fixed" = "species_latin_fixed")) |> 
  select(
    Order = "id_order",
    Artsnavn = species_latin_fixed,
    "Listet i Saul" = saul_boolean,
    "Listet i Daisie" = daisie_boolean
  ) |> 
  arrange(Artsnavn) |> 
  distinct() |> 
  collect()




## ------------------------------------------------------------------------------------------------------------------------------------
#| label: get_alien_obs_new
#| cache: !expr cache_var

alien_obs_q <- "
SELECT
foo.id_order,
foo.species_latin_fixed,
foo.\"risikokategori_2023\",
foo.habitat_type,
count(distinct foo.sampling_name)::integer as antall_tilfeller,
count(distinct foo.locality)::integer as antall_lokaliteter
FROM
(SELECT
obs.id_order,
obs.species_latin_fixed,
alien.\"risikokategori_2023\",
l.habitat_type,
ls.sampling_name,
l.locality
FROM occurrences.observations obs,
lookup.horizon_scan_2025 alien,
events.identifications,
events.year_locality yl,
events.locality_sampling ls,
locations.localities l,
events.sampling_trap st
WHERE obs.identification_id = identifications.id
AND obs.species_latin_fixed = alien.species_latin_fixed
AND identifications.sampling_trap_id = st.id
AND st.locality_sampling_id = ls.id
AND ls.year_locality_id = yl.id
AND yl.locality_id = l.id
AND yl.project_short_name = 'NorIns'
AND alien.\"risikokategori_2023\" IN ('SE', 'HI', 'PH', 'LO')
AND obs.identification_confidence = 'HIGH') foo
GROUP BY id_order, species_latin_fixed, \"risikokategori_2023\", habitat_type
ORDER BY  \"risikokategori_2023\", species_latin_fixed, habitat_type
"

alien_obs <- dbGetQuery(
  con,
  alien_obs_q
) %>%
  mutate(habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type)) %>%
  select(
    Orden = id_order,
    Artsnavn = species_latin_fixed,
    Risikovurdering = risikokategori_2023,
    Habitatstype = habitat_type,
    `Antall lokaliteter` = antall_lokaliteter,
    `Antall tilfeller` = antall_tilfeller
  ) %>%
  arrange(Orden,
          Artsnavn)


## ------------------------------------------------------------------------------------------------------------------------------------
#| label: table_alien_obs
#| results: asis

if (names(rmarkdown::metadata$format) == "ninareport-pdf") {
  print(xtable(alien_obs),
        size = "footnotesize")
} else if (names(rmarkdown::metadata$format) == "ninareport-html") {
  DT::datatable(alien_obs,
    style = "bootstrap",
    options = list(
      pageLength = 10,
      autoWidth = TRUE,
      columnDefs = list(list(
        width = "100px",
        targets = 1
      ))
    )
  )
} else {
  kable(alien_obs)
}

