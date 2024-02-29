## ----load_packages, include = F-------------------------------------------------------------------------------------------------------
require(NinaR)
require(Norimon)
require(DBI)
require(RPostgres)
require(ggplot2)
require(xtable)
require(sf)
require(ggplot2)
require(ggforce)
# require(bookdown)
require(knitr)
require(cowplot)
require(grid)
require(gridExtra)
require(landscapemetrics)
require(fasterize)
require(sf)
require(tmap)
require(betapart)
require(vegan)
library(PerformanceAnalytics)
require(lme4)
# require(lmerTest)
require(MuMIn)
require(iNEXT)
require(latex2exp)
require(ggvegan)
require(bookdown)
require(phyloseq)
require(Polychrome)
require(ggthemes)
require(sjPlot)
suppressPackageStartupMessages(require(tidyverse))


## -------------------------------------------------------------------------------------------------------------------------------------
opts_chunk$set(
  tidy = T,
  tidy.opts = list(width.cutoff = 50),
  echo = F,
  message = F,
  warning = F,
  # dpi = 600,
  # dev = c("pdf", "png", "svg"),
  # fig.path = "figure/",
  cache.extra = 4 # rename this to clear cache
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


## -------------------------------------------------------------------------------------------------------------------------------------
cache_var <- TRUE
recalc <- FALSE


## -------------------------------------------------------------------------------------------------------------------------------------
connect_to_insect_db()




## -------------------------------------------------------------------------------------------------------------------------------------
load("out/richness.RData")


## -------------------------------------------------------------------------------------------------------------------------------------
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


## -------------------------------------------------------------------------------------------------------------------------------------
localities <- get_localities(dataset = "NasIns") %>%
  mutate(habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type))

n_localities <- localities %>%
  st_drop_geometry() %>%
  summarize(n_distinct(locality)) %>%
  pull()

year_locality_biodiv <- get_observations(
  agg_level = "year_locality",
  dataset = "NasIns",
  as_tibble = T
) %>%
  mutate(
    habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type),
    year = as.factor(year),
    region_name = factor(region_name, levels = c(
      "Østlandet",
      "Sørlandet",
      "Trøndelag",
      "Nord-Norge"
    ))
  )


locality_sampling_biodiv <- get_observations(
  agg_level = "locality_sampling",
  dataset = "NasIns",
  as_tibble = T
)
# locality_sampling_biodiv

biomass_mf_locality_sampling <- get_biomass(
  agg_level = "locality_sampling",
  dataset = "NasIns",
  trap_type = "MF",
  subset_region = NULL
) %>%
  mutate(log_bm_per_day = log(sum_wet_weight / no_trap_days))

# biomass_mf_locality_sampling

year_locality_biomass <- get_biomass(
  dataset = "NasIns",
  trap_type = "MF",
  agg_level = "year_locality",
  subset_region = NULL
)


## -------------------------------------------------------------------------------------------------------------------------------------
localities_point <- localities %>%
  mutate(geom_point = st_centroid(geom)) %>%
  st_set_geometry("geom_point")


## -------------------------------------------------------------------------------------------------------------------------------------
norge <- get_map()


## -------------------------------------------------------------------------------------------------------------------------------------
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


## -------------------------------------------------------------------------------------------------------------------------------------
nas_ins_lok_1000m_agg_q <- "
SELECT ar5.*, ar5_kode.egenskapsverdi as arealtype
FROM backgrounds.locality_1000m_buffer_arealtype_agg ar5,
lookup.ar5_arealtype ar5_kode,
events.year_locality yl
WHERE ar5.locality_id = yl.locality_id
AND yl.project_short_name = 'NasIns'
AND ar5.\"arealressursArealtype\" = ar5_kode.kodeverdi

"
nas_ins_lok_1000m_agg <- read_sf(con,
  query = nas_ins_lok_1000m_agg_q
)




## -------------------------------------------------------------------------------------------------------------------------------------
load(file = "out/nas_ins_lok_1000m_agg_brick.Rdata")
nas_ins_lok_1000m_arealtype_shannon_index <- lsm_l_shdi(nas_ins_lok_1000m_agg_brick)
nas_ins_lok_1000m_arealtype_shannon_index$layer <- names(nas_ins_lok_1000m_agg_brick)

locality_ar5_shdi <- nas_ins_lok_1000m_arealtype_shannon_index %>%
  dplyr::select(
    locality = layer,
    lsm_l_shdi = value
  )


## -------------------------------------------------------------------------------------------------------------------------------------
locality_ar5_tot_perimeter <- nas_ins_lok_1000m_agg %>%
  mutate(perimeter = lwgeom::st_perimeter(geom)) %>%
  st_drop_geometry() %>%
  group_by(locality) %>%
  summarize(tot_perimeter_length = sum(perimeter)) %>%
  select(
    locality,
    tot_perimeter_length
  ) %>%
  mutate(tot_perimeter_length = units::drop_units(tot_perimeter_length))

# locality_ar5_tot_perimeter %>%
#  print(n = Inf)


## -------------------------------------------------------------------------------------------------------------------------------------
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


## -------------------------------------------------------------------------------------------------------------------------------------
loc_types <- localities %>%
  select(
    locality,
    habitat_type,
    year
  )

locality_ar5_area_data_2023 <- locality_ar5_area %>%
  left_join(loc_types,
    by = c("locality" = "locality"),
    copy = T
  ) %>%
  filter(year == 2023) %>% # to look at only 2023, comment out
  ungroup() %>%
  select(-c(
    "locality",
    "habitat_type",
    "Ikke kartlagt",
    "Hav",
    "geom"
  ))

locality_ar5_area_group_2023 <- locality_ar5_area %>%
  left_join(loc_types,
    by = c("locality" = "locality"),
    copy = T
  ) %>%
  filter(year == 2023) %>% # to look at only 2023, comment out
  ungroup() %>%
  select(c("locality", "habitat_type")) %>%
  mutate(habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type))


## -------------------------------------------------------------------------------------------------------------------------------------
example_NMDS <- metaMDS(locality_ar5_area_data_2023, # Our community-by-species matrix
  k = 2
)


## -------------------------------------------------------------------------------------------------------------------------------------
ordiplot(example_NMDS, type = "n")

ordihull(example_NMDS,
  groups = locality_ar5_area_group_2023$habitat_type,
  draw = "polygon",
  col = "grey90",
  label = T,
  cex = 1
)

fixed_loc_names <- locality_ar5_area_group_2023$locality %>%
  gsub("(.*_)(.)", "\\2", .)

orditorp(example_NMDS,
  display = "sites",
  labels = fixed_loc_names,
  col = as.character(factor(locality_ar5_area_group_2023$habitat_type,
    labels = ninaColors()[5:4]
  )),
  air = 0.05,
  cex = 0.7
)

orditorp(example_NMDS,
  display = "species",
  col = ninaColors()[2],
  air = 0.01,
  cex = 0.8
)


## -------------------------------------------------------------------------------------------------------------------------------------
p1 <- locality_ar5_plot("Semi-nat_87",
  legend = F,
  strip = T
)
p2 <- locality_ar5_plot("Semi-nat_71",
  legend = F,
  strip = T
)

p3 <- locality_ar5_plot("Semi-nat_81",
  legend = F,
  strip = T
)

p4 <- locality_ar5_plot("Semi-nat_100",
  legend = F,
  strip = T
)


leg <- locality_ar5_plot("Semi-nat_36",
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

legend <- as_grob(get_legend(leg))

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


## -------------------------------------------------------------------------------------------------------------------------------------
p1 <- locality_ar5_plot("Skog_34",
  legend = F,
  strip = T
)

p2 <- locality_ar5_plot("Skog_31",
  legend = F,
  strip = T
)

p3 <- locality_ar5_plot("Skog_38",
  legend = F,
  strip = T
)

p4 <- locality_ar5_plot("Skog_33",
  legend = F,
  strip = T
)

leg <- locality_ar5_plot("Skog_12",
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

legend <- as_grob(get_legend(leg))

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


## -------------------------------------------------------------------------------------------------------------------------------------
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
                          AND yl.project_short_name = 'NasIns'

                          "
)

biomass_per_loc_year_raw <- get_biomass(
  trap_type = "All",
  agg_level = "year_locality",
  dataset = "NasIns",
  subset_region = NULL
)

# tot_biomass <- get_biomass(trap_type = "All",
#                                          agg_level = "total",
#                                          dataset = "NasIns"
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
  min(ls.start_date), max(ls.end_date), round((sum(st.wet_weight) / 1000)::numeric, 2) as kg
                          FROM events.sampling_trap st,
                          events.locality_sampling ls,
                          events.year_locality yl
                          WHERE st.locality_sampling_id = ls.id
                          AND st.received_at_lab IS TRUE
                          AND ls.year_locality_id = yl.id
                          AND yl.project_short_name = 'NasIns'
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
    `Antall prøver` = no_samples,
    `Våtvekt (kg)` = kg
  ) %>%
  janitor::adorn_totals("row",
    name = "Sum"
  )


## -------------------------------------------------------------------------------------------------------------------------------------
if (names(rmarkdown::metadata$format) == "ninareport-pdf") {
  xtable(biomass_per_year)
} else {
  kable(biomass_per_year)
}


## -------------------------------------------------------------------------------------------------------------------------------------
no_samplings_per_lok_q <- "
SELECT l.region_name, l.locality,
l.habitat_type,
count(distinct ls.sampling_name)::integer,
min(ls.start_date) as  start_date,
max(ls.end_date) as end_date
FROM events.locality_sampling ls,
events.year_locality yl,
locations.localities l,
events.sampling_trap st
WHERE ls.year_locality_id = yl.id
AND yl.locality_id = l.id
AND yl.year = 2023
AND yl.project_short_name = 'NasIns'
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


## -------------------------------------------------------------------------------------------------------------------------------------
if (names(rmarkdown::metadata$format) == "ninareport-pdf") {
  print(xtable(no_samplings_per_lok), size = "\\footnotesize")
} else if (names(rmarkdown::metadata$format) == "ninareport-html") {
  # kable(no_samplings_per_lok)
  DT::datatable(no_samplings_per_lok,
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
  kable(no_samplings_per_lok)
}


## -------------------------------------------------------------------------------------------------------------------------------------
no_samples_q <- "

SELECT trap_model as \"Felletype\",
count(distinct st.sample_name)::integer as \"Antall prosesserte felleprøver\"
FROM events.locality_sampling ls,
events.year_locality yl,
locations.traps,
events.sampling_trap st
WHERE ls.year_locality_id = yl.id
AND st.trap_id = traps.id
AND yl.year = 2023
AND yl.project_short_name = 'NasIns'
AND st.locality_sampling_id = ls.id
AND st.received_at_lab IS TRUE
GROUP BY trap_model
"

no_samples <- dbGetQuery(
  con,
  no_samples_q
) %>%
  mutate(Felletype = ifelse(Felletype == "Malaise_Watkins_Black", "Malaisefelle", "Vindusfelle"))


## -------------------------------------------------------------------------------------------------------------------------------------
if (names(rmarkdown::metadata$format) == "ninareport-pdf") {
  xtable(no_samples)
} else {
  kable(no_samples)
}


## -------------------------------------------------------------------------------------------------------------------------------------
include_graphics("figures/hans_merged.jpg")


## -------------------------------------------------------------------------------------------------------------------------------------
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
  Porter = c(
    4778,
    135,
    268,
    246,
    65660,
    39693
  ),
  `Modifisert 2021` = c(
    4952,
    147,
    377,
    246,
    69426,
    39695
  ),
  `Modifisert 2022` = c(
    5124,
    147,
    388,
    246,
    72151,
    39698
  )
)

if (names(rmarkdown::metadata$format) == "ninareport-pdf") {
  xtable(refbase_tab, digits = 0)
} else {
  kable(refbase_tab, digits = 0)
}


## -------------------------------------------------------------------------------------------------------------------------------------
if (names(rmarkdown::metadata$format) == "ninareport-pdf") {
  include_graphics("figures/flyt_2023_1.pdf")
} else {
  include_graphics("figures/flyt_2023_1.png")
}


## -------------------------------------------------------------------------------------------------------------------------------------
if (names(rmarkdown::metadata$format) == "ninareport-pdf") {
  include_graphics("figures/flyt_2023_2.pdf")
} else {
  include_graphics("figures/flyt_2023_2.png")
}


## -------------------------------------------------------------------------------------------------------------------------------------
if (names(rmarkdown::metadata$format) == "ninareport-pdf") {
  include_graphics("figures/flyt_3.pdf")
} else {
  include_graphics("figures/flyt_3.png")
}


## -------------------------------------------------------------------------------------------------------------------------------------
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
AND yl.project_short_name = 'NasIns'
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


## -------------------------------------------------------------------------------------------------------------------------------------
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


## -------------------------------------------------------------------------------------------------------------------------------------
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


## -------------------------------------------------------------------------------------------------------------------------------------
order_spec_mod <- glm(to_species ~ order_rank * order, family = binomial(), data = rank_abund_order)
# summary(order_spec_mod)


## -------------------------------------------------------------------------------------------------------------------------------------
new_data_limits <- rank_abund_order %>%
  group_by(order) %>%
  summarise(no_otu = max(order_rank))


## -------------------------------------------------------------------------------------------------------------------------------------
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


## -------------------------------------------------------------------------------------------------------------------------------------
Afig <- rank_abund_long %>%
  # filter(rank_abundance <= 13200) %>%
  ggplot(., aes(x = rank_abundance, y = value, fill = taxonomic_level)) +
  geom_bar(stat = "identity", position = "fill", width = 100) +
  ylab("Andel arter identifisert") +
  xlab("Arter i rekkefølge av mengde DNA") +
  scale_fill_nina() +
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
  geom_vline(xintercept = min(rank_abund$row_number[rank_abund$pct_all_sequences > 0.899999999]), linetype = "dashed")

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
    legend.text = element_text(size = rel(0.45)),
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
  )

#+
# xlim(c(0, 6000))
# +  scale_y_continuous(trans = "logit")

grid.arrange(Afig, Bfig, ncol = 2)


## -------------------------------------------------------------------------------------------------------------------------------------
spec_seq_depth_q <- "
SELECT yl.year,
l.habitat_type, st.sample_name, count(distinct o.species_latin)::numeric as species_richn, i.read_abundance::numeric
FROM occurrences.observations o,
events.identifications i,
events.sampling_trap st,
events.locality_sampling ls,
events.year_locality yl,
locations.localities l,
locations.traps,
lookup.trap_types tt
WHERE o.identification_id = i.id
AND i.sampling_trap_id = st.id
AND st.locality_sampling_id = ls.id
AND ls.year_locality_id = yl.id
AND yl.locality_id = l.id
AND st.trap_id = traps.id
AND traps.trap_model = tt.trap_model
AND tt.trap_type = 'Malaise'
--AND yl.year = 2023
--OR yl.year = 2021) --could include 2021 as well
GROUP BY st.sample_name, yl.year, i.read_abundance, l.habitat_type
"

spec_seq_depth <- dbGetQuery(
  con,
  spec_seq_depth_q
) %>%
  as_tibble() %>%
  mutate(
    habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type),
    year = as.factor(year),
    scaled_species_richn = sjmisc::center(species_richn),
    scaled_read_abundance = sjmisc::center(read_abundance)
  )

# Trond_2022_Semi-nat_53_sampling-01_MF1 looks to be an outlier. Typing error?

# spec_seq_depth <- spec_seq_depth %>%
#  filter(sample_name != "Trond_2022_Semi-nat_53_sampling-01_MF1")


## -------------------------------------------------------------------------------------------------------------------------------------
## seq_depth_mod <- glm(scaled_species_richn ~ scaled_read_abundance +
##   habitat_type +
##   # read_abundance : habitat_type +
##   year +
##   scaled_read_abundance:year, data = spec_seq_depth)
## 
## summary(seq_depth_mod)
## 
## mydf <- ggeffects::ggpredict(
##   terms = c("scaled_read_abundance", "year"),
##   model = seq_depth_mod
## )
## ggplot(
##   mydf,
##   aes(
##     x = x,
##     y = predicted,
##     colour = group
##   )
## ) +
##   geom_line() +
##   geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.1, )


## -------------------------------------------------------------------------------------------------------------------------------------
spec_seq_depth %>%
  filter(
    read_abundance < 1e+06,
    year == 2023
  ) %>%
  ggplot(., aes(x = read_abundance, y = species_richn, col = habitat_type)) +
  geom_point() +
  theme(panel.background = element_blank(), axis.line = element_line(color = "black")) +
  xlab("Sekvenseringsdybde") +
  ylab("Antall arter") +
  scale_color_manual(values = c("#004F71", "#008C95"), name = "Habitat")


## -------------------------------------------------------------------------------------------------------------------------------------
physeq.species <- readRDS("figures/ord_physeq.rds")


## -------------------------------------------------------------------------------------------------------------------------------------
# col_arne <- openxlsx::read.xlsx("../../../../../2022/analysis/insectsurvey_2022/rawData/validering/Coleoptera_man_ident_2022_Arne.xlsx") %>%
#   as_tibble()

col_arne <- openxlsx::read.xlsx("../../rawData/validering/Coleoptera_morf_ident_2023_Arne.xlsx") %>%
  as_tibble()

col_man <- col_arne %>%
  mutate(
    species_latin = merged.name,
    sample_name = gsub("(.*)(-4$)", "\\1", sample_name)
  ) %>%
  select(
    sample_name,
    species_latin,
    amount,
    id_certainty,
    id_comment
  ) %>%
  filter(!is.na(sample_name)) %>%
  distinct()

# col_man
col_val_samples <- col_man %>%
  select(sample_name) %>%
  distinct() %>%
  pull()

col_val_samples <- paste0(
  c(
    "('",
    paste0(col_val_samples,
      collapse = "','"
    ),
    "')"
  ),
  collapse = ""
)


## -------------------------------------------------------------------------------------------------------------------------------------
col_dna_q <- paste0(
  "
SELECT distinct on(st.sample_name, o.species_latin_fixed)
st.sample_name,
o.species_latin_fixed as species_latin,
taxonomic_level
FROM occurrences.observations o,
events.identifications i,
events.sampling_trap st
WHERE o.identification_id = i.id
AND i.sampling_trap_id = st.id
AND id_order = 'Coleoptera'
AND st.sample_name in ",
  col_val_samples
)

col_dna <- dbGetQuery(
  con,
  col_dna_q
)

# col_dna %>% filter(is.na(taxonomic_level))


## -------------------------------------------------------------------------------------------------------------------------------------
col_merged <- col_man %>%
  full_join(col_dna,
    by = c(
      "sample_name" = "sample_name",
      "species_latin" = "species_latin"
    ),
    keep = FALSE,
    suffix = c("_man", "_dna")
  )

# col_man %>% filter(is.na(id_certainty))

col_merged <- col_merged %>%
  mutate(
    man_to_species = !grepl("(.p\\.)", species_latin) & !is.na(id_certainty),
    man_to_genus = grepl("(.p\\.)", species_latin) & !is.na(id_certainty),
    dna_to_species = ifelse(taxonomic_level == "species", TRUE, FALSE),
    dna_to_genus = ifelse(taxonomic_level == "genus", TRUE, FALSE)
  )

openxlsx::write.xlsx(col_merged,
  file = "../../rawData/validering/col_merged.xlsx"
)


## -------------------------------------------------------------------------------------------------------------------------------------
# Got col_merged back after manual revision

# col_merged <- openxlsx::read.xlsx("../../rawData/validering/col_merged_back.xlsx",
#                                    startRow = 2) %>%
#   as_tibble()
# #col_merged

# using this as placeholder until I get it back
# tax_unc should have values both for manual and metabarcoding.
# like tax_unc_man, tax_unc_dna, talk to Arne of which id_certainty to set at uncertain.

col_merged <- openxlsx::read.xlsx("../../rawData/validering/col_merged_jens_mod.xlsx") %>%
  as_tibble() %>%
  mutate(tax_unc = ifelse(id_certainty %in% c("3", "2", "1"), TRUE, FALSE))


## -------------------------------------------------------------------------------------------------------------------------------------
col_agg <- col_merged %>%
  group_by(
    sample_name,
    species_latin
  ) %>%
  mutate(n_dup = n()) %>%
  filter(n_dup < 2) %>%
  group_by(sample_name) %>%
  summarise(
    tot_spec = n(),
    no_man_to_species = sum(man_to_species, na.rm = TRUE),
    no_man_to_genus = sum(man_to_genus, na.rm = TRUE),
    no_dna_to_species = sum(dna_to_species, na.rm = TRUE),
    no_dna_to_genus = sum(dna_to_genus, na.rm = TRUE)
    # no_tax_uncertain = sum(tax_unc, na.rm = TRUE)
  ) %>%
  mutate(
    no_dna_missed = tot_spec - (no_dna_to_species + no_dna_to_genus),
    no_man_missed = tot_spec - (no_man_to_species + no_man_to_genus)
  )

col_longer <- col_agg %>%
  select(-tot_spec) %>%
  pivot_longer(cols = !"sample_name") %>%
  mutate(
    method = ifelse(grepl("man", name), "Morfologi", "Metastrekkoding"),
    kategori = ifelse(grepl("to_species", name), "Bestemt til art", NA)
  ) %>%
  mutate(kategori = ifelse(grepl("to_genus", name), "Bestemt til slekt", kategori)) %>%
  mutate(kategori = ifelse(grepl("missed", name), "Ikke oppdaget", kategori)) %>%
  # mutate(kategori = ifelse(grepl("uncertain", name), "Usikker manuell best.", kategori)) %>%
  mutate(kategori = factor(kategori, levels = c("Ikke oppdaget", "Bestemt til slekt", "Bestemt til art"))) %>%
  mutate(locality = gsub("(Ost_2023_)(.*)(_sampling.*)", "\\2", sample_name))


## -------------------------------------------------------------------------------------------------------------------------------------
col_val_plot <- ggplot(
  col_longer,
  aes(
    x = method,
    y = value,
    fill = kategori
  )
) +
  geom_bar(stat = "identity") +
  facet_wrap(~locality,
    ncol = 5
  ) +
  # scale_fill_nina() +
  scale_fill_manual(values = c("grey85", "#008C95", "#004F71")) +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 90),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title.x = element_blank(),
    strip.text = element_text(size = 5.5)
  ) +
  scale_x_discrete(labels = c("DNA", "Morf")) +
  ylab("Antall arter") +
  ggtitle("Coleoptera")


## -------------------------------------------------------------------------------------------------------------------------------------
col_man_matrix <- col_merged %>%
  filter(!is.na(species_latin)) %>%
  mutate(
    locality = gsub("(Ost_2023_)(.*)(_sampling.*)", "\\2", sample_name),
    present = 1
  ) %>%
  select(
    locality,
    species_latin,
    present
  ) %>%
  distinct() %>%
  pivot_wider(
    names_from = species_latin,
    values_from = present,
    values_fill = 0
  )

col_dna_matrix <- col_merged %>%
  filter(!is.na(taxonomic_level)) %>%
  mutate(
    locality = gsub("(Ost_2023_)(.*)(_sampling.*)", "\\2", sample_name),
    present = 1
  ) %>%
  select(
    locality,
    species_latin,
    present
  ) %>%
  distinct() %>%
  pivot_wider(
    names_from = species_latin,
    values_from = present,
    values_fill = 0
  )

col_to_iNEXT <- list(
  t(col_man_matrix[-1]),
  t(col_dna_matrix[-1])
)

names(col_to_iNEXT) <- c("Morf", "DNA")


## -------------------------------------------------------------------------------------------------------------------------------------
col_inext <- iNEXT(col_to_iNEXT,
  datatype = "incidence_raw",
  endpoint = 50
)

col_fort <- fortify(col_inext)


## -------------------------------------------------------------------------------------------------------------------------------------
syrph_sondre <- openxlsx::read.xlsx("../../rawData/validering/Syrphidae_morf_ident_2023.xlsx") %>%
  as_tibble()

syrph_man <- syrph_sondre %>%
  mutate(
    species_latin = paste0(str_trim(genus), " ", str_trim(species)),
    sample_name = gsub("(.*)(_MF1$)", "\\1", sample_name)
  ) %>%
  mutate(sample_name = paste0(sample_name, "_MF1")) %>%
  select(
    sample_name,
    genus,
    species,
    species_latin
  ) %>%
  filter(!is.na(genus)) %>%
  # mutate(species_latin = paste0(genus, " ", species)) %>%
  distinct()

# syrph_man
syrph_val_samples <- syrph_man %>%
  select(sample_name) %>%
  distinct() %>%
  pull()

syrph_val_samples <- paste0(
  c(
    "('",
    paste0(syrph_val_samples,
      collapse = "','"
    ),
    "')"
  ),
  collapse = ""
)


## -------------------------------------------------------------------------------------------------------------------------------------
# The Sør samples had single digits as sample events, but the rest had 2 digit values. (I updated the trigger after the Sør locations was initiated). Doing this manually.
syrph_dna_q <- paste0(
  "
SELECT distinct on(st.sample_name, o.species_latin_fixed)
st.sample_name,
o.species_latin_fixed as species_latin,
taxonomic_level
FROM occurrences.observations o,
events.identifications i,
events.sampling_trap st
WHERE o.identification_id = i.id
AND i.sampling_trap_id = st.id
AND id_family = 'Syrphidae'
AND st.sample_name in ",
  syrph_val_samples
)

syrph_dna <- dbGetQuery(
  con,
  syrph_dna_q
)


## -------------------------------------------------------------------------------------------------------------------------------------
syrph_merged <- syrph_man %>%
  full_join(syrph_dna,
    by = c(
      "sample_name" = "sample_name",
      "species_latin" = "species_latin"
    ),
    keep = FALSE,
    suffix = c("_man", "_dna")
  )

syrph_merged <- syrph_merged %>%
  mutate(
    man_to_family = grepl("Ukjent", genus),
    man_to_species = !is.na(species),
    man_to_genus = !grepl("Ukjent", genus) & grepl("NA", species_latin),
    dna_to_species = taxonomic_level == "species",
    dna_to_genus = taxonomic_level == "genus"
  )

openxlsx::write.xlsx(syrph_merged,
  file = "../../rawData/validering/Syrphidae_merged_2023.xlsx"
)


## -------------------------------------------------------------------------------------------------------------------------------------
# Got syrph_merged back after manual check
syrph_merged <- openxlsx::read.xlsx("../../rawData/validering/Syrphidae_merged_2023_JÅ_mod.xlsx") %>%
  as_tibble()

# syrph_merged


## -------------------------------------------------------------------------------------------------------------------------------------
syrph_agg <- syrph_merged %>%
  select(
    sample_name,
    species_latin,
    man_to_species,
    man_to_genus,
    man_likely_error,
    dna_to_species,
    dna_to_genus,
    dna_likely_error,
    dna_ambig
  ) %>%
  distinct() %>%
  group_by(sample_name) %>%
  summarise(
    tot_spec = n(),
    # no_man_to_family = sum(man_to_family, na.rm = TRUE),
    no_man_to_species = sum(man_to_species, na.rm = TRUE),
    no_man_to_genus = sum(man_to_genus, na.rm = TRUE),
    no_man_error = sum(man_likely_error, na.rm = TRUE),
    no_dna_to_species = sum(dna_to_species, na.rm = TRUE),
    no_dna_to_genus = sum(dna_to_genus, na.rm = TRUE),
    no_dna_ambig = sum(dna_ambig, na.rm = TRUE),
    no_dna_error = sum(dna_likely_error, na.rm = TRUE)
  ) %>%
  mutate(
    no_dna_missed = tot_spec - (no_dna_to_species + no_dna_to_genus + no_dna_error + no_dna_ambig),
    no_man_missed = tot_spec - (no_man_to_species + no_man_to_genus + no_man_error)
  )

syrph_longer <- syrph_agg %>%
  select(-tot_spec) %>%
  pivot_longer(cols = !"sample_name") %>%
  mutate(
    method = ifelse(grepl("man", name), "Morfologi", "Metastrekkoding"),
    kategori = ifelse(grepl("to_species", name), "Bestemt til art", NA)
  ) %>%
  mutate(
    kategori = ifelse(grepl("to_genus", name), "Bestemt til slekt", kategori),
    kategori = ifelse(grepl("to_family", name), "Bestemt til familie", kategori),
    kategori = ifelse(grepl("no_dna_ambig", name), "DNA kan ikke skille\nmellom arter", kategori),
    kategori = ifelse(grepl("error", name), "Feilbestemt", kategori)
  ) %>%
  mutate(kategori = ifelse(grepl("missed", name), "Ikke oppdaget", kategori)) %>%
  mutate(kategori = factor(kategori, levels = c("Ikke oppdaget", "DNA kan ikke skille\nmellom arter", "Bestemt til familie", "Bestemt til slekt", "Bestemt til art", "Feilbestemt"))) %>%
  mutate(locality = gsub("(.*)(_2023_)(.*)(_sampling.*)", "\\3", sample_name))


## -------------------------------------------------------------------------------------------------------------------------------------
syrph_val_plot <- ggplot(
  syrph_longer,
  aes(
    x = method,
    y = value,
    fill = kategori
  )
) +
  geom_bar(stat = "identity") +
  facet_wrap(~locality,
    ncol = 10
  ) +
  # scale_fill_nina() +
  scale_fill_manual(values = c("grey85", "#E57200", "#008C95", "#004F71", "#93328E")) +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 90),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title.x = element_blank(),
    strip.text = element_text(
      size = 5.5,
      angle = 90
    )
  ) +
  scale_x_discrete(labels = c("DNA", "Morf")) +
  ylab("Antall arter") +
  ggtitle("Syrphidae")


## -------------------------------------------------------------------------------------------------------------------------------------
syrph_man_matrix <- syrph_merged %>%
  filter(!is.na(species)) %>%
  mutate(
    locality = gsub("(.*)(_2023_)(.*)(_sampling.*)", "\\3", sample_name),
    present = as.integer(1)
  ) %>%
  select(
    locality,
    species_latin,
    present
  ) %>%
  distinct() %>%
  pivot_wider(
    names_from = species_latin,
    values_from = present,
    values_fill = 0
  )

syrph_dna_matrix <- syrph_merged %>%
  filter(!is.na(taxonomic_level)) %>%
  mutate(
    locality = gsub("(.*)(_2023_)(.*)(_sampling.*)", "\\3", sample_name),
    present = as.integer(1)
  ) %>%
  select(
    locality,
    species_latin,
    present
  ) %>%
  distinct() %>%
  pivot_wider(
    names_from = species_latin,
    values_from = present,
    values_fill = 0
  )

syrph_to_iNEXT <- list(
  t(syrph_man_matrix[-1]),
  t(syrph_dna_matrix[-1])
)

names(syrph_to_iNEXT) <- c("Morf", "DNA")


## -------------------------------------------------------------------------------------------------------------------------------------
syrph_inext <- iNEXT(syrph_to_iNEXT,
  datatype = "incidence_raw",
  endpoint = 50
)

syrph_fort <- fortify(syrph_inext)


## -------------------------------------------------------------------------------------------------------------------------------------
col_fort <- col_fort %>%
  mutate(color = "Col_C")

syrph_fort <- syrph_fort %>%
  mutate(color = "Col_S")

df <- col_fort %>%
  rbind(syrph_fort)

df.point <- df[which(df$Method == "Observed"), ]
df.point$shape[grepl("Morf", df.point$Assemblage)] <- "Morfologi"
df.point$shape[grepl("DNA", df.point$Assemblage)] <- "DNA"
# df.point$color=gsub("Morf|DNA","",df.point$site)
# df.point$color=gsub("DNA","",df.point$color)

df.line <- df[which(df$Method != "Observed"), ]
df.line$method <- factor(
  df.line$Method,
  c("Rarefaction", "Extrapolation"),
  c("Rarefaction", "Extrapolation")
)
# df.line$color=gsub("Morph|DNA","",df.line$site)
df.line$shape[grepl("Morf", df.line$Assemblage)] <- "Morfologi"
df.line$shape[grepl("DNA", df.line$Assemblage)] <- "DNA"

df$shape[grepl("Morph", df$Assemblage)] <- "Morfologi"
df$shape[grepl("DNA", df$Assemblage)] <- "DNA"
# df$color=gsub("Morph|DNA","",df$site)


## -------------------------------------------------------------------------------------------------------------------------------------
lines <- c(
  "Rarefaction" = "solid",
  "Extrapolation" = "dotted"
)

ggplot(df, aes(x = x, y = y, colour = color, shape = shape, fill = color)) +
  geom_point(aes(shape = shape, colour = color), size = 5, data = df.point) +
  geom_line(aes(linetype = method), linewidth = 1.5, data = df.line) +
  geom_ribbon(aes(
    ymin = y.lwr, ymax = y.upr,
    fill = color, colour = NULL
  ), alpha = 0.2, data = df) +
  labs(x = "Antall feller", y = "Antall arter") +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    text = element_text(size = 10),
    legend.box = "vertical"
  ) +
  scale_fill_manual(
    values = c("#004F71", "#E57200"),
    labels = c("Coleoptera vindusf.", "Syrphidae malaisef.")
  ) +
  scale_color_manual(
    values = c("#004F71", "#E57200"),
    labels = c("Coleoptera vindusf.", "Syrphidae malaisef.")
  ) +
  xlim(c(0, 20)) +
  scale_linetype_manual(
    values = lines,
    name = "Method",
    labels = c("Interpolert", "Ekstrapolert")
  )


## -------------------------------------------------------------------------------------------------------------------------------------
grid.arrange(col_val_plot,
  syrph_val_plot,
  ncol = 1
)


## -------------------------------------------------------------------------------------------------------------------------------------
include_graphics("figures/syrphidae_sortering.jpg")


## -------------------------------------------------------------------------------------------------------------------------------------
include_graphics("figures/frysereol1.jpg")


## -------------------------------------------------------------------------------------------------------------------------------------
include_graphics("figures/container_merged.jpg")


## -------------------------------------------------------------------------------------------------------------------------------------
boxes_per_shelf <- 6 * 3 * 2
shelves_per_section <- 4
bottles_per_box <- 24

bottles_per_seminat_year <- 10 * 14
bottles_per_forest_year <- 10 * 2 * 14

no_seminat_per_year <- 5
no_forest_per_year <- 1
no_forest_per_year_full <- 5

no_bottles_per_year <- bottles_per_seminat_year * no_seminat_per_year + bottles_per_forest_year * no_forest_per_year
no_bottles_per_year_full <- bottles_per_seminat_year * no_seminat_per_year + bottles_per_forest_year * no_forest_per_year_full

bottles_per_section <- boxes_per_shelf * bottles_per_box * shelves_per_section

year_per_section <- bottles_per_section / no_bottles_per_year
year_per_section_full <- bottles_per_section / no_bottles_per_year_full


## -------------------------------------------------------------------------------------------------------------------------------------
# AJ-kasse er 60 * 40 * 23 cm
# 40 ft fryskontainere har innvendig mål:
# 1160 * 229 * 261 cm

# 40ft
# 11 på høyden
# 2 på bredden + 1 meter gang
# 28 på lengden
# 616 bokser totalt

# 20ft frysekontainer har innvendig mål:
# 546 * 229 * 230

# 20ft
# 9 på høyden
# 2 på bredden + 1 meter gang
# 14 på lengden
# 252 bokser totalt


no_bokses_height_40ft <- 11
no_boxes_width_40ft <- 2 # incl 1 meter walkway
no_boxes_length_40ft <- 28

no_bokses_height_20ft <- 9
no_boxes_width_20ft <- 2 # incl 1 meter walkway
no_boxes_length_20ft <- 14


no_boxes_in_cont_40ft <- no_bokses_height_40ft * no_boxes_width_40ft * no_boxes_length_40ft
no_bottles_in_cont_40ft <- no_boxes_in_cont_40ft * bottles_per_box
year_per_cont_40ft_full <- no_bottles_in_cont_40ft / no_bottles_per_year_full # 7.5
year_per_cont_40ft_as_2024 <- no_bottles_in_cont_40ft / no_bottles_per_year # 15


no_boxes_in_cont_20ft <- no_bokses_height_20ft * no_boxes_width_20ft * no_boxes_length_20ft
no_bottles_in_cont_20ft <- no_boxes_in_cont_20ft * bottles_per_box
year_per_cont_20ft_full <- no_bottles_in_cont_20ft / no_bottles_per_year_full # 2.88
year_per_cont_20ft_as_2024 <- no_bottles_in_cont_20ft / no_bottles_per_year # 6


## -------------------------------------------------------------------------------------------------------------------------------------
ano_herb_agg <- tbl(
  con,
  Id(
    schema = "views",
    table = "ano_herb_agg"
  )
) %>%
  filter(project_short_name == "NasIns")


## -------------------------------------------------------------------------------------------------------------------------------------
loc_reg <- localities %>%
  st_drop_geometry() %>%
  select(
    locality,
    region_name
  ) %>%
  mutate(region_name = factor(region_name, levels = c("Østlandet", "Sørlandet", "Trøndelag", "Nord-Norge")))


## -------------------------------------------------------------------------------------------------------------------------------------
plot_ano_herb_sum()


## ---- cache = cache_var---------------------------------------------------------------------------------------------------------------
year_locality_biodiv_ano <- year_locality_biodiv %>%
  left_join(ano_herb_agg,
    by = c("locality" = "locality"),
    copy = T
  )


## -------------------------------------------------------------------------------------------------------------------------------------
tree_agg <- tbl(
  con,
  Id(
    schema = "views",
    table = "landsskog_agg"
  )
) %>%
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


## -------------------------------------------------------------------------------------------------------------------------------------
tree_agg_to_plot <- tree_agg %>%
  group_by(
    year,
    dom_tree
  ) %>%
  summarise(
    avg_avg_age = mean(avg_age),
    sd_avg_age = sd(avg_age, na.rm = TRUE)
  )


## -------------------------------------------------------------------------------------------------------------------------------------
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
  scale_x_continuous(breaks = seq(2020, 2023)) +
  scale_color_nina(name = "Dominerende\ntreslag") +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  xlab("År") +
  ylab("Skogsalder")


## -------------------------------------------------------------------------------------------------------------------------------------
include_graphics("figures/skog_38_00023.jpg")


## -------------------------------------------------------------------------------------------------------------------------------------
clim_trend_plot(
  region = "Trøndelag",
  from_month = 6,
  to_month = 8,
  from_year = 2013,
  to_year = 2023
)


## -------------------------------------------------------------------------------------------------------------------------------------
clim_trend_plot(
  region = "Østlandet",
  from_month = 6,
  to_month = 8,
  from_year = 2013,
  to_year = 2023
)


## -------------------------------------------------------------------------------------------------------------------------------------
clim_trend_plot(
  region = "Sørlandet",
  from_month = 6,
  to_month = 8,
  from_year = 2013,
  to_year = 2023
)


## -------------------------------------------------------------------------------------------------------------------------------------
clim_trend_plot(
  region = "Nord-Norge",
  from_month = 6,
  to_month = 8,
  from_year = 2013,
  to_year = 2023
)


## -------------------------------------------------------------------------------------------------------------------------------------
climdata_semi_nat_86 <- get_climate_data(locality = "Semi-nat_86")


## -------------------------------------------------------------------------------------------------------------------------------------
plot_climate_comparison(climdata_semi_nat_86,
  main_title = FALSE,
  focus_year = 2023,
  rolling_mean = TRUE
)


## -------------------------------------------------------------------------------------------------------------------------------------
plot_climate_comparison(climdata_semi_nat_86,
  main_title = FALSE,
  variable = "precipitation",
  rolling_mean = TRUE,
  focus_year = 2023,
  y_low_limit = -20
)


## -------------------------------------------------------------------------------------------------------------------------------------
climdata_semi_nat_73 <- get_climate_data(locality = "Semi-nat_73")


## -------------------------------------------------------------------------------------------------------------------------------------
plot_climate_comparison(climdata_semi_nat_73,
  main_title = FALSE,
  focus_year = 2023,
  rolling_mean = TRUE
)


## -------------------------------------------------------------------------------------------------------------------------------------
#|
plot_climate_comparison(climdata_semi_nat_73,
  main_title = FALSE,
  variable = "precipitation",
  rolling_mean = TRUE,
  focus_year = 2023,
  y_low_limit = -10
)


## -------------------------------------------------------------------------------------------------------------------------------------
climdata_semi_nat_61 <- get_climate_data(locality = "Semi-nat_61")


## -------------------------------------------------------------------------------------------------------------------------------------
plot_climate_comparison(climdata_semi_nat_61,
  main_title = FALSE,
  focus_year = 2023,
  rolling_mean = TRUE
)


## -------------------------------------------------------------------------------------------------------------------------------------
plot_climate_comparison(climdata_semi_nat_61,
  main_title = FALSE,
  variable = "precipitation",
  rolling_mean = TRUE,
  focus_year = 2023,
  y_low_limit = -10,
  y_high_limit = 80
)


## -------------------------------------------------------------------------------------------------------------------------------------
climdata_semi_nat_96 <- get_climate_data(locality = "Semi-nat_96")


## -------------------------------------------------------------------------------------------------------------------------------------
plot_climate_comparison(climdata_semi_nat_96,
  main_title = FALSE,
  focus_year = 2023,
  y_high_limit = 70,
  rolling_mean = TRUE
)


## -------------------------------------------------------------------------------------------------------------------------------------
plot_climate_comparison(climdata_semi_nat_96,
  main_title = FALSE,
  variable = "precipitation",
  rolling_mean = TRUE,
  focus_year = 2023,
  y_low_limit = -10
)


## -------------------------------------------------------------------------------------------------------------------------------------
include_graphics("figures/templogger.jpg")


## -------------------------------------------------------------------------------------------------------------------------------------
include_graphics("figures/lyslogger_small.png")




## -------------------------------------------------------------------------------------------------------------------------------------
load(file = "out/locality_sampling_loggerdata.Rdata")


## -------------------------------------------------------------------------------------------------------------------------------------
## locality_sampling_loggerdata %>%
##   select(
##     Temp_MX2301A = avg_values_MX2301A_temperature,
##     Temp_MX2202 = avg_values_MX2202_temperature,
##     RL = avg_values_MX2301A_rh,
##     Lys = avg_values_MX2202_light
##   ) %>%
##   filter(Lys >= 100) %>%
##   Norimon::chart_correlation(.,
##     histogram = TRUE,
##     method = "pearson",
##     color = ninaColors()[1]
##   )


## -------------------------------------------------------------------------------------------------------------------------------------
locality_sampling_weight_loggerdata <- biomass_mf_locality_sampling %>%
  left_join(locality_sampling_loggerdata,
    by = c(
      "sampling_name" = "sampling_name",
      "year" = "year",
      "locality" = "locality",
      "habitat_type" = "habitat_type"
    )
  ) %>%
  mutate(habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type))


## -------------------------------------------------------------------------------------------------------------------------------------
#|
locality_sampling_weight_loggerdata %>%
  group_by(habitat_type) %>%
  ggplot(aes( # y = log_bm_per_day,
    y = avg_wet_weight / no_trap_days,
    x = avg_values_MX2301A_temperature
  )) +
  geom_point(aes(color = habitat_type)) +
  scale_color_nina(
    name = "Habitatstype",
    palette = "darkblue-orange"
  ) +
  geom_smooth(aes(color = habitat_type),
    method = "gam",
    # fullrange = T,
    se = T,
    linewidth = 2
  ) +
  ylab("Avrunnen våtvekt (g.) per felledag") +
  xlab("Middeltemperatur ved felle i fangstperioden")


## -------------------------------------------------------------------------------------------------------------------------------------
locality_sampling_biodiv_loggerdata <- locality_sampling_biodiv %>%
  left_join(locality_sampling_loggerdata,
    by = c(
      "sampling_name" = "sampling_name",
      "year" = "year",
      "locality" = "locality",
      "habitat_type" = "habitat_type"
    )
  ) %>%
  mutate(habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type))


## -------------------------------------------------------------------------------------------------------------------------------------
locality_sampling_biodiv_loggerdata %>%
  group_by(habitat_type) %>%
  ggplot(aes(
    y = no_species,
    x = avg_values_MX2301A_temperature
  )) +
  geom_point(aes(color = habitat_type)) +
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


## -------------------------------------------------------------------------------------------------------------------------------------
nas_ins_lok_1000m_skog_q <- "
SELECT skog.*
FROM backgrounds.locality_1000m_buffer_skogtype_agg skog,
events.year_locality yl
WHERE skog.locality_id = yl.locality_id
AND yl.project_short_name = 'NasIns'

"

nas_ins_lok_1000m_skog <- read_sf(con,
  query = nas_ins_lok_1000m_skog_q
)


## -------------------------------------------------------------------------------------------------------------------------------------
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

## -------------------------------------------------------------------------------------------------------------------------------------
locality_skog_area_data <- locality_skog_area %>%
  ungroup() %>%
  select(-"locality")

skog_rda <- decorana(locality_skog_area_data)

# summary(skog_rda) #hm, too large axis length?


## -------------------------------------------------------------------------------------------------------------------------------------
skog_PCA <- rda(locality_skog_area_data, scale = F)

# summary(skog_PCA)


## -------------------------------------------------------------------------------------------------------------------------------------
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


## -------------------------------------------------------------------------------------------------------------------------------------
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


## -------------------------------------------------------------------------------------------------------------------------------------
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


## -------------------------------------------------------------------------------------------------------------------------------------
hm <- autoplot(ar5_PCA, legend.position = "none") +
  scale_color_nina(palette = "orange-green")


## -------------------------------------------------------------------------------------------------------------------------------------
hm


## -------------------------------------------------------------------------------------------------------------------------------------
tt <- autoplot(skog_PCA,
  legend.position = "none"
) +
  scale_color_nina(palette = "orange-green")


## -------------------------------------------------------------------------------------------------------------------------------------
tt


## -------------------------------------------------------------------------------------------------------------------------------------
loc_reg_year <- localities %>%
  st_drop_geometry() %>%
  select(year, locality, region_name, habitat_type) %>%
  mutate(region_name = factor(region_name, levels = c("Østlandet", "Sørlandet", "Trøndelag", "Nord-Norge")))


## -------------------------------------------------------------------------------------------------------------------------------------
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
  scale_x_continuous(breaks = seq(2020, 2023)) +
  guides(
    color = guide_legend(override.aes = list(size = 3)),
    pch = guide_legend(
      title = "Habitattype",
      override.aes = list(size = 3)
    )
  ) +
  xlab("År") +
  ylab("Shannon index")


## -------------------------------------------------------------------------------------------------------------------------------------
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
  scale_x_continuous(breaks = seq(2020, 2023)) +
  guides(
    color = guide_legend(override.aes = list(size = 3)),
    pch = guide_legend(
      title = "Habitattype",
      override.aes = list(size = 3)
    )
  ) +
  xlab("År") +
  ylab("Tot. lengde lineære element")




## -------------------------------------------------------------------------------------------------------------------------------------
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





## -------------------------------------------------------------------------------------------------------------------------------------
locality_sampling_biomass_temp_compl_light_week <- read_csv(file = "out/locality_sampling_biomass_temp_compl_light_week.csv") %>%
  mutate(
    year_num = Year - 2020,
    region_name = factor(region_name, levels = c("Østlandet", "Sørlandet", "Trøndelag", "Nord-Norge"))
  )


## -------------------------------------------------------------------------------------------------------------------------------------
## #|
## 
## # sanity check of data (is there something wrong with the logger data in 2023?)
## locality_sampling_biomass_temp_compl_light_week %>%
##   mutate(year = as_factor(Year)) %>%
##   ggplot(.) +
##   geom_boxplot(aes(
##     y = no_trap_days,
##     x = year
##   ))
## 
## locality_sampling_biomass_temp_compl_light_week %>%
##   mutate(year = as_factor(Year)) %>%
##   ggplot(.) +
##   geom_boxplot(aes(
##     y = avg_wet_weight,
##     x = year
##   ))
## 
## locality_sampling_biomass_temp_compl_light_week %>%
##   mutate(year = as_factor(Year)) %>%
##   ggplot(.) +
##   geom_boxplot(aes(
##     y = avg_values_MX2301A_temperature,
##     x = year
##   ))
## 
## locality_sampling_biomass_temp_compl_light_week %>%
##   mutate(year = as_factor(Year)) %>%
##   ggplot(.) +
##   geom_boxplot(aes(
##     y = avg_values_MX2202_light,
##     x = year,
##     fill = region_name
##   ))
## 
## locality_sampling_biomass_temp_compl_light_week %>%
##   mutate(year = as_factor(Year)) %>%
##   ggplot(.) +
##   geom_boxplot(aes(
##     y = avg_values_MX2301A_rh,
##     x = year,
##     fill = region_name
##   ))
## 
## locality_sampling_biomass_temp_compl_light_week %>%
##   mutate(year = as_factor(Year)) %>%
##   ggplot(.) +
##   geom_boxplot(aes(
##     y = stddev_values_MX2301A_temperature,
##     x = year,
##     fill = region_name
##   ))
## 
## locality_sampling_biomass_temp_compl_light_week %>%
##   mutate(year = as_factor(Year)) %>%
##   ggplot(.) +
##   geom_boxplot(aes(
##     y = stddev_values_MX2301A_rh,
##     x = year,
##     fill = region_name
##   ))
## 
## locality_sampling_biomass_temp_compl_light_week %>%
##   mutate(year = as_factor(Year)) %>%
##   ggplot(.) +
##   geom_boxplot(aes(
##     y = stddev_values_MX2202_light,
##     x = year,
##     fill = region_name
##   ))


## -------------------------------------------------------------------------------------------------------------------------------------
locality_ar5_area <- locality_ar5_area %>%
  left_join(ar5_scores,
    by = c("locality" = "locality")
  )




## -------------------------------------------------------------------------------------------------------------------------------------
locality_sampling_biomass_big_join <- read_csv(file = "out/locality_sampling_biomass_big_join.csv")

locality_sampling_biomass_big_join <- locality_sampling_biomass_big_join %>%
  mutate(region_name = factor(region_name, levels = c("Østlandet", "Sørlandet", "Trøndelag", "Nord-Norge")))

locality_sampling_biomass_big_join_forest <- locality_sampling_biomass_big_join %>%
  filter(habitat_type == "Forest") %>%
  mutate(habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type))

locality_sampling_biomass_big_join_semi_nat <- locality_sampling_biomass_big_join %>%
  filter(habitat_type == "Semi-nat")


## -------------------------------------------------------------------------------------------------------------------------------------
locality_sampling_biomass_big_join_forest %>%
  ggplot(.) +
  geom_point(
    aes(
      x = lsm_l_shdi,
      y = tot_perimeter_length
    ),
    color = ninaColors()[1],
    size = 3
  ) +
  ylab("Total lengde lineære element") +
  xlab("Artype Shannon indeks")


## -------------------------------------------------------------------------------------------------------------------------------------
locality_sampling <- tbl(
  con,
  Id(
    schema = "events",
    table = "locality_sampling"
  )
) %>%
  select(
    sampling_name,
    start_date,
    end_date
  ) %>%
  collect()


## -------------------------------------------------------------------------------------------------------------------------------------
biomass_raw <- get_biomass(
  agg_level = "locality_sampling",
  trap_type = "MF",
  subset_region = NULL
)

biomass_mf_locality_sampling_time <- biomass_raw %>%
  left_join(locality_sampling,
    by = c("sampling_name" = "sampling_name"),
    copy = T
  ) %>%
  mutate(
    julian_day = lubridate::yday(end_date),
    habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type),
    start_month = lubridate::month(start_date)
  ) %>%
  filter(sum_wet_weight > 0)


## -------------------------------------------------------------------------------------------------------------------------------------
## biomass_mf_locality_sampling_time %>%
##   mutate(series = factor(
##     interaction(habitat_type,
##       region_name,
##       sep = "_"
##     ),
##     levels = c(
##       "Semi-nat_Østlandet",
##       "Skog_Østlandet",
##       "Semi-nat_Trøndelag",
##       "Semi-nat_Sørlandet",
##       "Semi-nat_Nord-Norge"
##     )
##   )) %>%
##   ggplot(
##     .,
##     aes(
##       y = sum_wet_weight / no_trap_days,
##       x = as.Date(julian_day, origin = as.Date("2010-01-01")),
##       group = series,
##       color = series
##     )
##   ) +
##   geom_line(
##     stat = "smooth",
##     formula = y ~ x,
##     cex = 2,
##     alpha = 0.8
##   ) +
##   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
##   scale_color_nina(name = "") +
##   facet_wrap(~year,
##     ncol = 1
##   ) +
##   ylab("Biomasse per felledøgn (gram)") +
##   xlab("Tømmedato") +
##   scale_x_date(date_labels = "%d.%b") +
##   theme(
##     legend.position = "bottom",
##     legend.text = element_text(size = 7)
##   ) +
##   guides(color = guide_legend(nrow = 2, byrow = TRUE))


## -------------------------------------------------------------------------------------------------------------------------------------
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
    year == 2023
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
    year == 2023
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
    limits = as.Date(xlims, origin = as.Date("2023-01-01"))
  ) +
  scale_y_continuous(limits = ylims) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ylab("Biomasse per\nfelledøgn (gram)") +
  xlab("Tømmedato") +
  scale_color_nina(name = "Lokalitet") +
  theme(
    legend.position = "bottom",
    # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 5),
    legend.key.width = unit(0.5, "cm")
  )
#+
# annotate("text", x = as.Date("2023-10-1"), y = 4, label = "bold('2023')", parse = TRUE)

# marrangeGrob(list(p1,
#              p2,
#              p3,
#              p4),
#              top = "",
#              ncol = 1,
#              nrow = 4)

p4


## -------------------------------------------------------------------------------------------------------------------------------------
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
    year == 2023
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
    year == 2023
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
  ylab("Biomasse per felledøgn (gram)") +
  scale_x_date(
    date_labels = "%d.%m.%y",
    limits = as.Date(xlims, origin = as.Date("2023-01-01"))
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
#+
# annotate("text", x = as.Date("2023-10-1"), y = 6, label = "bold('2023')", parse = TRUE)

# marrangeGrob(list(p1,
#                   p2 ),
#              top = "",
#              ncol = 1,
#              nrow = 2)

p2


## -------------------------------------------------------------------------------------------------------------------------------------
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
    year == 2023
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
    year == 2023
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
    limits = as.Date(xlims, origin = as.Date("2023-01-01"))
  ) +
  scale_y_continuous(limits = ylims) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ylab("Biomasse per\nfelledøgn (gram)") +
  xlab("Tømmedato") +
  scale_color_nina(name = "Lokalitet") +
  theme(
    legend.position = "bottom",
    # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 5),
    legend.key.width = unit(0.5, "cm")
  )
#+
# annotate("text", x = as.Date("2023-10-1"), y = 8, label = "bold('2023')", parse = TRUE)

# marrangeGrob(list(p1,
#              p2,
#              p3,
#              p4),
#              top = "",
#              ncol = 1,
#              nrow = 4)

p4


## -------------------------------------------------------------------------------------------------------------------------------------
xlims <- biomass_mf_locality_sampling_time %>%
  filter(
    habitat_type == "Semi-nat",
    region_name == "Trøndelag"
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
    year == 2023
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
    year == 2023
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
  ylab("Biomasse per felledøgn (gram)") +
  scale_x_date(
    date_labels = "%d.%m.%y",
    limits = as.Date(xlims, origin = as.Date("2023-01-01"))
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
#+
# annotate("text", x = as.Date("2023-10-1"), y = 12, label = "bold('2023')", parse = TRUE)

# marrangeGrob(list(p1,
#                   p2,
#                   p3),
#              top = "",
#              ncol = 1,
#              nrow = 3)

p3


## -------------------------------------------------------------------------------------------------------------------------------------
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

ylims <- biomass_mf_locality_sampling_time %>%
  filter(
    habitat_type == "Semi-nat",
    region_name == "Nord-Norge",
    year == 2023
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
    year == 2023
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
  ylab("Biomasse per felledøgn (gram)") +
  scale_x_date(
    date_labels = "%d.%m.%y",
    limits = as.Date(xlims, origin = as.Date("2023-01-01"))
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
#+
# annotate("text", x = as.Date("2023-10-1"), y = 12, label = "bold('2023')", parse = TRUE)

# marrangeGrob(list(p1,
#                   p2,
#                   p3),
#              top = "",
#              ncol = 1,
#              nrow = 3)

p1


## -------------------------------------------------------------------------------------------------------------------------------------
locality_sampling_biodiv_sampling_time <- locality_sampling_biodiv %>%
  left_join(locality_sampling,
    by = c("sampling_name" = "sampling_name"),
    copy = T
  ) %>%
  mutate(
    julian_day = lubridate::yday(end_date),
    habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type)
  )


## -------------------------------------------------------------------------------------------------------------------------------------
## locality_sampling_biodiv_sampling_time %>%
##   mutate(series = factor(
##     interaction(habitat_type,
##       region_name,
##       sep = "_"
##     ),
##     levels = c(
##       "Semi-nat_Østlandet",
##       "Skog_Østlandet",
##       "Semi-nat_Trøndelag",
##       "Semi-nat_Sørlandet"
##     )
##   )) %>%
##   ggplot(
##     .,
##     aes(
##       y = no_species / no_trap_days,
##       x = as.Date(julian_day, origin = as.Date("2010-01-01")),
##       group = series,
##       color = series
##     )
##   ) +
##   geom_smooth(
##     formula = y ~ x,
##     cex = 2
##   ) +
##   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
##   scale_color_nina(name = "") +
##   facet_wrap(~year,
##     ncol = 1
##   ) +
##   ylab("Artsantall per felledøgn") +
##   xlab("Tømmedato") +
##   scale_x_date(date_labels = "%d.%b") +
##   theme(
##     legend.position = "bottom",
##     legend.text = element_text(size = 7)
##   ) +
##   guides(color = guide_legend(
##     nrow = 2,
##     byrow = TRUE
##   ))


## -------------------------------------------------------------------------------------------------------------------------------------
year_locality_div <- get_observations(
  subset_region = "Østlandet",
  agg_level = "year_locality",
  trap_type = "MF"
) %>%
  mutate(habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type))


## -------------------------------------------------------------------------------------------------------------------------------------
xlims <- locality_sampling_biodiv_sampling_time %>%
  filter(
    habitat_type == "Skog",
    region_name == "Østlandet"
  ) %>%
  select(julian_day) %>%
  pull()

xlims <- c(
  min(xlims - 1, na.rm = T),
  max(xlims + 1, na.rm = T)
)

ylims <- locality_sampling_biodiv_sampling_time %>%
  filter(
    habitat_type == "Skog",
    region_name == "Østlandet"
  ) %>%
  mutate(amount = no_species / no_trap_days) %>%
  pull()

ylims <- c(
  min(ylims, na.rm = T) * 1.1,
  max(ylims, na.rm = T) * 1.1
)

p2 <- locality_sampling_biodiv_sampling_time %>%
  filter(
    habitat_type == "Semi-nat",
    region_name == "Sørlandet",
    year == 2023
  ) %>%
  ggplot(.) +
  geom_line(
    aes(
      y = no_species / no_trap_days,
      x = as.Date(end_date),
      color = locality
    ),
    linewidth = 0.75
  ) +
  geom_point(
    aes(
      y = no_species / no_trap_days,
      x = as.Date(end_date),
      color = locality
    ),
    size = 2.5
  ) +
  scale_x_date(
    date_labels = "%d.%m.%y",
    limits = as.Date(xlims, origin = as.Date("2023-01-01"))
  ) +
  scale_y_continuous(limits = ylims) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ylab("Artsantall per\nfelledøgn") +
  xlab("Tømmedato") +
  scale_color_nina(name = "Lokalitet") +
  theme(
    legend.position = "bottom",
    # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 5),
    legend.key.width = unit(0.5, "cm")
  )
#+
# annotate("text", x = as.Date("2023-10-1"), y = 100, label = "bold('2023')", parse = TRUE)

# marrangeGrob(list(p1,
#                   p2),
#              top = "",
#              ncol = 1,
#              nrow = 2)

p2 # replace with p2


## -------------------------------------------------------------------------------------------------------------------------------------
p4 <- locality_sampling_biodiv_sampling_time %>%
  filter(
    habitat_type == "Skog",
    region_name == "Østlandet",
    year == 2023
  ) %>%
  ggplot(.) +
  geom_line(
    aes(
      y = no_species / no_trap_days,
      x = as.Date(end_date),
      color = locality
    ),
    linewidth = 0.75
  ) +
  geom_point(
    aes(
      y = no_species / no_trap_days,
      x = as.Date(end_date),
      color = locality
    ),
    size = 2.5
  ) +
  scale_x_date(
    date_labels = "%d.%m.%y",
    limits = as.Date(xlims, origin = as.Date("2023-01-01"))
  ) +
  scale_y_continuous(limits = ylims) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ylab("Artsantall per felledøgn") +
  xlab("Tømmedato") +
  scale_color_nina(name = "Lokalitet") +
  theme(
    legend.position = "bottom",
    # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 5),
    legend.key.width = unit(0.5, "cm")
  )
#+
# annotate("text", x = as.Date("2023-10-1"), y = 100, label = "bold('2023')", parse = TRUE)

# marrangeGrob(list(p1,
#                   p2,
#                   p3),
#              top = "",
#              ncol = 1,
#              nrow = 3)

p4 # change to p4


## -------------------------------------------------------------------------------------------------------------------------------------
xlims <- locality_sampling_biodiv_sampling_time %>%
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

ylims <- locality_sampling_biodiv_sampling_time %>%
  filter(
    habitat_type == "Semi-nat",
    region_name == "Østlandet"
  ) %>%
  mutate(amount = no_species / no_trap_days) %>%
  pull()

ylims <- c(
  min(ylims, na.rm = T) * 1.1,
  max(ylims, na.rm = T) * 1.1
)

p4 <- locality_sampling_biodiv_sampling_time %>%
  filter(
    habitat_type == "Semi-nat",
    region_name == "Østlandet",
    year == 2023
  ) %>%
  ggplot(.) +
  geom_line(
    aes(
      y = no_species / no_trap_days,
      x = as.Date(end_date),
      color = locality
    ),
    linewidth = 0.75
  ) +
  geom_point(
    aes(
      y = no_species / no_trap_days,
      x = as.Date(end_date),
      color = locality
    ),
    size = 2.5
  ) +
  scale_x_date(
    date_labels = "%d.%m.%y",
    limits = as.Date(xlims, origin = as.Date("2023-01-01"))
  ) +
  scale_y_continuous(limits = ylims) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ylab("Artsantall per felledøgn") +
  xlab("Tømmedato") +
  scale_color_nina(name = "Lokalitet") +
  theme(
    legend.position = "bottom",
    # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 5),
    legend.key.width = unit(0.5, "cm")
  )
#+
# annotate("text", x = as.Date("2023-10-1"), y = 100, label = "bold('2023')", parse = TRUE)

# marrangeGrob(list(p1,
#                   p2,
#                   p3,
#                   p4),
#              top = "",
#              ncol = 1,
#              nrow = 4)
p4 # change to p4


## -------------------------------------------------------------------------------------------------------------------------------------
p3 <- locality_sampling_biodiv_sampling_time %>%
  filter(
    habitat_type == "Semi-nat",
    region_name == "Trøndelag",
    year == 2023
  ) %>%
  ggplot(.) +
  geom_line(
    aes(
      y = no_species / no_trap_days,
      x = as.Date(end_date),
      color = locality
    ),
    linewidth = 0.75
  ) +
  geom_point(
    aes(
      y = no_species / no_trap_days,
      x = as.Date(end_date),
      color = locality
    ),
    size = 2.5
  ) +
  scale_x_date(
    date_labels = "%d.%m.%y",
    limits = as.Date(xlims, origin = as.Date("2023-01-01"))
  ) +
  scale_y_continuous(limits = ylims) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ylab("Artsantall per\nfelledøgn") +
  xlab("Tømmedato") +
  scale_color_nina(name = "Lokalitet") +
  theme(
    legend.position = "bottom",
    # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 5),
    legend.key.width = unit(0.5, "cm")
  )
#+
# annotate("text", x = as.Date("2023-10-1"), y = 100, label = "bold('2023')", parse = TRUE)

# marrangeGrob(list(p1,
#                   p2,
#                   p3),
#              top = "",
#              ncol = 1,
#              nrow = 3)

p3 # change to p3


## -------------------------------------------------------------------------------------------------------------------------------------
p1 <- locality_sampling_biodiv_sampling_time %>%
  filter(
    habitat_type == "Semi-nat",
    region_name == "Nord-Norge",
    year == 2023
  ) %>%
  ggplot(.) +
  geom_line(
    aes(
      y = no_species / no_trap_days,
      x = as.Date(end_date),
      color = locality
    ),
    linewidth = 0.75
  ) +
  geom_point(
    aes(
      y = no_species / no_trap_days,
      x = as.Date(end_date),
      color = locality
    ),
    size = 2.5
  ) +
  scale_x_date(
    date_labels = "%d.%m.%y",
    limits = as.Date(xlims, origin = as.Date("2023-01-01"))
  ) +
  scale_y_continuous(limits = ylims) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ylab("Artsantall per\nfelledøgn") +
  xlab("Tømmedato") +
  scale_color_nina(name = "Lokalitet") +
  theme(
    legend.position = "bottom",
    # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 5),
    legend.key.width = unit(0.5, "cm")
  )
#+
# annotate("text", x = as.Date("2023-10-1"), y = 100, label = "bold('2023')", parse = TRUE)

# marrangeGrob(list(p1,
#                   p2,
#                   p3),
#              top = "",
#              ncol = 1,
#              nrow = 3)

p1


## -------------------------------------------------------------------------------------------------------------------------------------
biomass_mf_year_locality <- get_biomass(
  trap_type = "MF",
  agg_level = "year_locality",
  subset_region = NULL
)


## -------------------------------------------------------------------------------------------------------------------------------------
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
    mean_biomass = mean(avg_wet_weight / no_trap_days),
    sd_biomass = sd(avg_wet_weight / no_trap_days),
    n = n(),
    se_biomass = sd_biomass / sqrt(n)
  ) %>%
  mutate(
    year = as.factor(year),
    region_name = as.factor(region_name),
    habitat_type = as.factor(habitat_type)
  ) %>%
  ungroup()

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
  scale_color_nina(
    name = "Region",
    palette = "contrast"
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


## -------------------------------------------------------------------------------------------------------------------------------------
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


## -------------------------------------------------------------------------------------------------------------------------------------
year_locality_pollinators <- get_observations(
  subset_families = pollinators,
  subset_region = "Østlandet",
  agg_level = "year_locality",
  trap_type = "MF"
) %>%
  mutate(habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type))


## -------------------------------------------------------------------------------------------------------------------------------------
year_mean_div <- year_locality_div %>%
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

marrangeGrob(
  list(
    p1,
    p2
  ),
  ncol = 2,
  nrow = 1,
  top = ""
)


## -------------------------------------------------------------------------------------------------------------------------------------
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


## -------------------------------------------------------------------------------------------------------------------------------------
## summary(year_model)


## -------------------------------------------------------------------------------------------------------------------------------------
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


## -------------------------------------------------------------------------------------------------------------------------------------
tt <- summary(year_model)
year_to_tab <- tt$coefficients %>%
  as_tibble(rownames = "Parameter") %>%
  mutate(Estimate.norm = exp(Estimate))

year_to_tab <- year_to_tab %>%
  mutate(Parameter = case_match(
    Parameter,
    "(Intercept)" ~ "Intercept",
    "year_num" ~ "År",
    "region_nameNord-Norge" ~ "Region: Nord-Norge",
    "region_nameSørlandet" ~ "Region: Sørlandet",
    "region_nameTrøndelag" ~ "Region: Trøndelag",
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


## -------------------------------------------------------------------------------------------------------------------------------------
if (names(rmarkdown::metadata$format) == "ninareport-pdf") {
  xtable(year_to_tab, digits = 3)
} else {
  kable(year_to_tab, digits = 3)
}


## -------------------------------------------------------------------------------------------------------------------------------------
r2_verdi_biomasse_year_model <- r.squaredGLMM(year_model)
r2_verdi_biomasse_year_model_rand_year <- r.squaredGLMM(year_model_rand_year)


## -------------------------------------------------------------------------------------------------------------------------------------
fe <- fixef(year_model)
# > paste0(names(fe), collapse = "', '")
# "(Intercept)', 'region_nameSørlandet', 'region_nameTrøndelag', 'region_nameNord-Norge', 'habitat_typeSemi-nat', 'year_num"

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


## -------------------------------------------------------------------------------------------------------------------------------------
rho <- cor.test(~ prediction + log_bm_per_day, data = to_cor_test)$estimate
# cor.test(~ prediction + log_bm_per_day, data = to_cor_test)
# 28% vs 11% Rho


## -------------------------------------------------------------------------------------------------------------------------------------
p1 <- plot_model(year_model,
  sort.est = TRUE,
  axis.labels = c(
    "year_num" = "År",
    "habitat_typeSemi-nat" = "Habitat: Semi-nat",
    "region_nameNord-Norge" = "Region: Nord-Norge",
    "region_nameTrøndelag" = "Region-Trøndelag",
    "region_nameSørlandet" = "Region-Sørlandet"
  ),
  title = ""
) +
  xlab("År ikke som tilfeldig effekt") +
  ggtitle("A") +
  theme(axis.title.y = element_text(size = 9))
# theme_sjplot()

p2 <- plot_model(year_model_rand_year,
  sort.est = TRUE,
  axis.labels = c(
    "year_num" = "År",
    "habitat_typeSemi-nat" = "Habitat: Semi-nat",
    "region_nameNord-Norge" = "Region: Nord-Norge",
    "region_nameTrøndelag" = "Region-Trøndelag",
    "region_nameSørlandet" = "Region-Sørlandet"
  ),
  title = ""
) +
  xlab("År som tilfeldig effekt") +
  ggtitle("B") +
  theme(axis.title.y = element_text(size = 9))

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


## ----vær-mod--------------------------------------------------------------------------------------------------------------------------
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


## -------------------------------------------------------------------------------------------------------------------------------------
## plot(global_model)
## 
## locality_sampling_biomass_temp_compl_light_week %>%
##   cbind(resid(global_model)) %>%
##   ggplot(.) +
##   geom_boxplot(aes(
##     x = `resid(global_model)`,
##     fill = locality,
##     group = locality
##   )) +
##   theme(legend.position = "none") +
##   coord_flip()


## -------------------------------------------------------------------------------------------------------------------------------------
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
    "habitat_typeSemi-nat" = "<b style='color:black'>Habitat: Semi-nat</b>",
    "region_nameNord-Norge" = "Region: Nord-Norge",
    "region_nameTrøndelag" = "<b style='color:black'>Region: Trøndelag</b>",
    "region_nameSørlandet" = "Region: Sørlandet",
    "scale_avg_values_MX2301A_temperature:scale_avg_values_MX2301A_rh" = "Temp:Fukt",
    "scale_avg_values_MX2301A_temperature:scale_avg_values_MX2202_light" = "Temp:Lys",
    "scale_avg_values_MX2301A_rh:scale_avg_values_MX2202_light" = "Fukt:Lys",
    "scale_stddev_values_MX2301A_temperature:scale_stddev_values_MX2301A_rh" = "Variasjon_temp:Variasjon_ fukt",
    "scale_stddev_values_MX2301A_temperature:scale_stddev_values_MX2202_light" = "<b style='color:black'>Variasjon_temp:Variasjon_lys</b>",
    "scale_stddev_values_MX2301A_rh:scale_stddev_values_MX2202_light" = "Variasjon_fukt:Variasjon_lys"
  ),
  title = ""
) +
  xlab("Forklaringsvariabel") +
  # theme_sjplot() +
  ggplot2::theme(axis.text.y = ggtext::element_markdown())

p


## -------------------------------------------------------------------------------------------------------------------------------------
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
    "habitat_typeSemi-nat" ~ "Habitat_semi-nat",
    "year_num" ~ "År",
    "region_nameNord-Norge" ~ "Region: Nord-Norge",
    "region_nameSørlandet" ~ "Region: Sørlandet",
    "region_nameTrøndelag" ~ "Region: Trøndelag",
    "habitat_typeSemi-nat" ~ "Habitat: Semi-nat",
    "scale_avg_values_MX2301A_temperature:scale_avg_values_MX2301A_rh" ~ "Temp.:fukt",
    "scale_avg_values_MX2301A_temperature:scale_avg_values_MX2202_light" ~ "Temp.:lys",
    "scale_avg_values_MX2301A_rh:scale_avg_values_MX2202_light" ~ "Fukt:lys",
    "scale_stddev_values_MX2301A_temperature:scale_stddev_values_MX2301A_rh" ~ "Var.temp.:Var.fukt",
    "scale_stddev_values_MX2301A_temperature:scale_stddev_values_MX2202_light" ~ "Var.temp.:Var.lys",
    "scale_stddev_values_MX2301A_rh:scale_stddev_values_MX2202_light" ~ "Var.fukt:Var.lys"
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


## -------------------------------------------------------------------------------------------------------------------------------------
if (names(rmarkdown::metadata$format) == "ninareport-pdf") {
  xtable(to_tab, digits = 3)
} else {
  kable(to_tab, digits = 3)
}


## -------------------------------------------------------------------------------------------------------------------------------------
# sw(model_averages)
r2_verdi_biomasse_temp_glob <- r.squaredGLMM(global_model)


## ----sub_model------------------------------------------------------------------------------------------------------------------------
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


## ---- eval = F------------------------------------------------------------------------------------------------------------------------
## summary(sub_model)


## -------------------------------------------------------------------------------------------------------------------------------------
r2_verdi_biomasse_sub <- r.squaredGLMM(sub_model)


## -------------------------------------------------------------------------------------------------------------------------------------
fe <- fixef(sub_model)
# > paste0(names(fe), collapse = "', '")
# "(Intercept)', 'scale_avg_values_MX2301A_temperature', 'scale_avg_values_MX2301A_rh', 'scale_avg_values_MX2202_light', 'scale_stddev_values_MX2301A_temperature', 'scale_stddev_values_MX2301A_rh', 'scale_stddev_values_MX2202_light', 'region_nameSørlandet', 'region_nameTrøndelag', 'region_nameNord-Norge', 'habitat_typeSemi-nat', 'year_num', 'scale_avg_values_MX2301A_temperature:scale_avg_values_MX2202_light', 'scale_stddev_values_MX2301A_temperature:scale_stddev_values_MX2202_light"

to_mult <- locality_sampling_biomass_temp_compl_light_week %>%
  mutate(
    `(Intercept)` = 1,
    `scale_avg_values_MX2301A_temperature:scale_avg_values_MX2202_light` = scale_avg_values_MX2301A_temperature * scale_avg_values_MX2202_light,
    `scale_avg_values_MX2202_light:scale_avg_values_MX2301A_rh` = scale_avg_values_MX2202_light * scale_avg_values_MX2301A_rh,
    `scale_stddev_values_MX2301A_temperature:scale_stddev_values_MX2202_light` = scale_stddev_values_MX2301A_temperature * scale_stddev_values_MX2202_light,
    # `temp_ano_winter:prec_ano_winter` = temp_ano_winter * prec_ano_winter,
    # `temp_april_cur:prec_april_cur` = temp_april_cur * prec_april_cur,
    # `temp_april_prev:prec_april_prev` = temp_april_prev * prec_april_prev,
    # `temp_ano_month_prev:prec_ano_month_prev` = temp_ano_month_prev * prec_ano_month_prev,
    `habitat_typeSemi-nat` = as.numeric(habitat_type == "Semi-nat")
  ) %>%
  select(
    `(Intercept)`,
    scale_avg_values_MX2301A_temperature,
    scale_avg_values_MX2301A_rh,
    scale_avg_values_MX2202_light,
    scale_stddev_values_MX2301A_temperature,
    scale_stddev_values_MX2301A_rh,
    scale_stddev_values_MX2202_light,
    `habitat_typeSemi-nat`,
    year_num,
    `scale_stddev_values_MX2301A_temperature:scale_stddev_values_MX2202_light`
  )


prediction <- as.matrix(to_mult) %*% fe

to_cor_test <- locality_sampling_biomass_temp_compl_light_week %>%
  select(
    Year,
    log_bm_per_day
  ) %>%
  cbind(prediction)


## -------------------------------------------------------------------------------------------------------------------------------------
rho <- cor.test(~ prediction + log_bm_per_day, data = to_cor_test)$estimate
# cor.test(~ prediction + log_bm_per_day, data = to_cor_test)
# 44% vs 73% R2 vs Spearman-Rho


## -------------------------------------------------------------------------------------------------------------------------------------
ggplot(to_cor_test) +
  geom_point(aes(
    x = prediction,
    y = log_bm_per_day,
    col = as.factor(Year)
  )) +
  # scale_y_log10() +
  scale_color_nina(name = "År") +
  annotate("text",
    x = -.8,
    y = 2,
    label = as.expression(bquote(R2 == .(round(r2_verdi_biomasse_sub[1], 2)))),
    hjust = 1
  ) +
  annotate("text",
    x = -0.8,
    y = 1.5,
    label = as.expression(bquote(rho == .(round(rho, 2)))),
    hjust = 1
  ) +
  xlab("Modellprediksjon") +
  ylab("Log biomasse per dag")


## -------------------------------------------------------------------------------------------------------------------------------------
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


## -------------------------------------------------------------------------------------------------------------------------------------
## biomass_mf_locality_sampling_weather_sum <- biomass_mf_locality_sampling %>%
##   left_join(weather_summary_scaled,
##     by = c("sampling_name" = "sampling_name")
##   )
## 
## ggplot(biomass_mf_locality_sampling_weather_sum) +
##   geom_point(aes(
##     x = scale(temp_ano_winter),
##     y = log_bm_per_day,
##     colour = as.factor(year.x)
##   )) +
##   scale_color_nina()


## -------------------------------------------------------------------------------------------------------------------------------------
locality_sampling_biomass_temp_compl_light_week_weather_anom <- locality_sampling_biomass_temp_compl_light_week %>%
  left_join(weather_summary_scaled,
    by = c(
      "sampling_name" = "sampling_name",
      "locality" = "locality"
    )
  )


## ----sub-model-anom-------------------------------------------------------------------------------------------------------------------
# anom_model <- lmer(log_bm_per_day ~
#                        scale_avg_values_MX2301A_temperature +
#                        scale_avg_values_MX2301A_rh +
#                        scale_avg_values_MX2202_light +
#                        scale_stddev_values_MX2301A_temperature +
#                        scale_stddev_values_MX2301A_rh +
#                        scale_stddev_values_MX2202_light +
#                        #scale_avg_values_MX2301A_temperature:scale_avg_values_MX2301A_rh +
#                        #scale_avg_values_MX2301A_temperature:scale_avg_values_MX2202_light +
#                        #scale_avg_values_MX2301A_rh:scale_avg_values_MX2202_light +
#                        #scale_stddev_values_MX2301A_temperature:scale_stddev_values_MX2301A_rh +
#                        scale_stddev_values_MX2301A_temperature:scale_stddev_values_MX2202_light +
#                        #scale_stddev_values_MX2301A_rh:scale_stddev_values_MX2202_light +
#                     #mean_temp +
#                     #sum_prec +
#                     temp_ano_winter +
#                     prec_ano_winter +
#                     temp_ano_winter:prec_ano_winter +
#                     #temp_april_cur +
#                     #prec_april_cur +
#                     temp_may_cur +
#                     prec_may_cur +
#                     #temp_april_cur:prec_april_cur +
#                     temp_may_cur:prec_may_cur +
#                     #temp_april_prev +
#                     #prec_april_prev +
#                     temp_may_prev +
#                     prec_may_prev +
#                     #temp_april_prev:prec_april_prev +
#                     temp_may_prev:prec_may_prev +
#                     temp_ano_month_prev +
#                     prec_ano_month_prev +
#                     temp_ano_month_prev:prec_ano_month_prev +
#                     #region_name +
#                     habitat_type +
#                     year_num +
#                     (1|Ukenummer) +
#                     (1|locality),
#                     data = locality_sampling_biomass_temp_compl_light_week_weather_anom,
#                     na.action = "na.fail",
#                     REML = F)

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
  na.action = "na.fail",
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
  na.action = "na.fail",
  REML = F
)


## -------------------------------------------------------------------------------------------------------------------------------------
## summary(anom_model)
## summary(anom_model_april)
## 
## AIC(
##   sub_model,
##   anom_model,
##   anom_model_april
## )


## -------------------------------------------------------------------------------------------------------------------------------------
r2_verdi_biomasse_anom <- r.squaredGLMM(anom_model)


## -------------------------------------------------------------------------------------------------------------------------------------
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
    "scale_stddev_values_MX2301A_temperature" = "<b style='color:black'>Variasjon temp.</b>",
    "scale_stddev_values_MX2202_light" = "<b style='color:black'>Variasjon lys</b>",
    "scale_avg_values_MX2301A_temperature:scale_avg_values_MX2301A_rh" = "Temp:Fukt",
    "scale_avg_values_MX2301A_temperature:scale_avg_values_MX2202_light" = "Temp:Lys",
    "scale_avg_values_MX2301A_rh:scale_avg_values_MX2202_light" = "Fukt:Lys",
    "scale_stddev_values_MX2301A_temperature:scale_stddev_values_MX2202_light" = "<b style='color:black'>Variasjon temp:Variasjon lys</b>",
    "scale_stddev_values_MX2301A_rh" = "<b style='color:black'>Variasjon fukt</b>",
    "year_num" = "År",
    "habitat_typeSemi-nat" = "<b style='color:black'>Habitat: Semi-nat</b>",
    "prec_ano_winter" = "<b style='color:black'>Nedbør vinter</b>",
    "temp_ano_month_prev" = "<b style='color:black'>Temp.-1år</b>",
    "temp_may_prev" = "Temp. mai-1år",
    "temp_may_cur:prec_may_cur" = "Temp.mai:Nedb.mai",
    "temp_ano_month_prev:prec_ano_month_prev" = "<b style='color:black'>Temp.-1år:Nedb.-1år</b>",
    "prec_ano_month_prev" = "<b style='color:black'>Nedb.-1år</b>",
    "temp_may_prev:prec_may_prev" = "Temp.mai-1år:Nedb.mai-1år",
    "temp_ano_winter:prec_ano_winter" = "<b style='color:black'>Temp.vinter:Nedb.vinter</b>",
    "prec_may_cur" = "<b style='color:black'>Nedbør mai</b>",
    "prec_may_prev" = "Nedbør mai-1år",
    "temp_ano_winter" = "Temp.vinter",
    "temp_may_cur" = "Temp.mai"
  ),
  title = ""
) +
  ggplot2::theme(axis.text.y = ggtext::element_markdown())


## -------------------------------------------------------------------------------------------------------------------------------------
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
  na.action = "na.fail",
  REML = F,
  control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-4))
)


## -------------------------------------------------------------------------------------------------------------------------------------
## locality_sampling_biomass_big_join_forest_2022 <- locality_sampling_biomass_big_join_forest %>%
##   filter(Year < 2023)
## 
## global_model_big_forest_2022 <- lmer(
##   log_bm_per_day ~
##     scale_avg_values_MX2301A_temperature +
##     scale_avg_values_MX2301A_rh +
##     scale_avg_values_MX2202_light +
##     scale_stddev_values_MX2301A_temperature +
##     scale_stddev_values_MX2301A_rh +
##     scale_stddev_values_MX2202_light +
##     # habitat_type +
##     scale_stddev_values_MX2301A_temperature:scale_stddev_values_MX2202_light +
##     scale(ar5_PC1) +
##     scale(ar5_PC2) +
##     scale(lsm_l_shdi) +
##     scale(tot_perimeter_length) +
##     scale(skog_PC1) +
##     scale(skog_PC2) +
##     scale(ano_mean_no_spec) +
##     scale(ano_mean_cover) +
##     scale(avg_age) +
##     dom_tree +
##     year_num +
##     (1 | locality) +
##     (1 | Ukenummer),
##   data = locality_sampling_biomass_big_join_forest,
##   na.action = "na.fail",
##   REML = F,
##   control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = 1e-4))
## )
## 
## cand_models_big_forest <- dredge(global_model_big_forest_2022,
##   fixed = c(
##     "scale(avg_values_MX2301A_temperature)",
##     # "scale(avg_values_MX2301A_rh)",
##     "scale(avg_values_MX2202_light)",
##     "scale(stddev_values_MX2301A_temperature)",
##     "scale(stddev_values_MX2301A_rh)"
##     # ,"scale(stddev_values_MX2202_light)"
##   )
## )
## 
## model_averages_big_forest <- model.avg(cand_models_big_forest, subset = delta <= 2)
## 
## par_weights_forest <- MuMIn::sw(get.models(cand_models_big_forest, subset = delta <= 4))


## -------------------------------------------------------------------------------------------------------------------------------------
## summary(global_model_big_forest)


## -------------------------------------------------------------------------------------------------------------------------------------
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
    "scale-avg_values_MX2301A_temperature:scale_avg_values_MX2301A_rh" ~ "Temp.:fukt",
    "scale_avg_values_MX2301A_temperature:scale_avg_values_MX2202_light" ~ "Temp.:lys",
    "scale_avg_values_MX2301A_rh:scale_avg_values_MX2202_light" ~ "Fukt:lys",
    "scale_stddev_values_MX2301A_temperature:scale_stddev_values_MX2202_light" ~ "Var.temp:Var.lys",
    "dom_treeDunbjørk" ~ "Dom. treslag: Dunbjørk ",
    "dom_treeFuru" ~ "Dom. treslag: Furu",
    "dom_treeGran" ~ "Dom. treslag: Gran",
    "dom_treeHengebjørk" ~ "Dom. treslag: Hengebjørk",
    "dom_treeSelje" ~ "Dom. treslag: Selje",
    "scale(avg_age)" ~ "Landsskog: skogalder",
    "scale(ano_mean_cover)" ~ "ANO: dekningsgrad karpl.",
    "scale(ano_mean_no_spec)" ~ "ANO: artsantall karpl.",
    "scale(ar5_PC1)" ~ "AR5-arealtype: PCA-akse 1",
    "scale(ar5_PC2)" ~ "AR5-arealtype: PCA-akse 2",
    "scale(skog_PC1)" ~ "AR5-skog: PCA-akse 1",
    "scale(skog_PC2)" ~ "AR5-skog: PCA-akse 2",
    "scale(lsm_l_shdi)" ~ "AR5-arealtype: Shannon indeks",
    "scale(tot_perimeter_length)" ~ "AR5-arealtype: lineære elementer"
  )) %>%
  filter(!(Parameter %in% c(
    "Intercept",
    "Luftfuktighet",
    "Variasjon_fukt",
    "Temperatur",
    "Lys",
    "Variasjon_temp.",
    "Temp.:lys",
    "Fukt:lys",
    "Variasjon_lys",
    "Var.temp:Var.lys"
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


## -------------------------------------------------------------------------------------------------------------------------------------
if (names(rmarkdown::metadata$format) == "ninareport-pdf") {
  xtable(to_tab, digits = 3)
} else {
  kable(to_tab, digits = 3)
}


## -------------------------------------------------------------------------------------------------------------------------------------
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
    "dom_treeHengebjørk" = "Dom. treslag: Hengebjørk",
    "dom_treeGran" = "Dom. treslag: Gran",
    "dom_treeFuru" = "Dom. treslag: Furu",
    "dom_treeSelje" = "Dom. treslag: Selje",
    "dom_treeDunbjørk" = "Dom. treslag: Dunbjørk",
    "scale(avg_age)" = "Landsskog: skogalder",
    "scale(ano_mean_cover)" = "ANO: dekningsgrad karpl.",
    "scale(ano_mean_no_spec)" = "ANO: artsantall karpl.",
    "scale(ar5_PC1)" = "AR5-arealtype: PCA-akse 1",
    "scale(ar5_PC2)" = "AR5-arealtype: PCA-akse 2",
    "scale(skog_PC1)" = "AR5-skog: PCA-akse 1",
    "scale(skog_PC2)" = "AR5-skog: PCA-akse 2",
    "scale(lsm_l_shdi)" = "AR5-arealtype: Shannon indeks",
    "scale(tot_perimeter_length)" = "AR5-arealtype: lineære elementer"
  ),
  title = ""
)


## -------------------------------------------------------------------------------------------------------------------------------------
r2_verdi_biomasse_big_forest <- r.squaredGLMM(global_model_big_forest)


## -------------------------------------------------------------------------------------------------------------------------------------
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


## -------------------------------------------------------------------------------------------------------------------------------------
## summary(global_model_big_semi_nat)


## -------------------------------------------------------------------------------------------------------------------------------------
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
    "habitat_typeSemi-nat" ~ "Habitat_semi-nat",
    "year_num" ~ "År",
    "scale_avg_values_MX2301A_temperature:scale_avg_values_MX2301A_rh" ~ "Temp.:fukt",
    "scale_avg_values_MX2301A_temperature:scale_avg_values_MX2202_light" ~ "Temp.:lys",
    "scale_avg_values_MX2301A_rh:scale_avg_values_MX2202_light" ~ "Fukt:lys",
    "scale_stddev_values_MX2301A_temperature:scale_stddev_values_MX2202_light" ~ "Var.temp:Var.lys",
    "dom_treeDunbjørk" ~ "Dom. treslag: Dunbjørk ",
    "dom_treeFuru" ~ "Dom. treslag: Furu",
    "dom_treeGran" ~ "Dom. treslag: Gran",
    "dom_treeHengebjørk" ~ "Dom. treslag: Hengebjørk",
    "dom_treeSelje" ~ "Dom. treslag: Selje",
    "scale(avg_age)" ~ "Landsskog: skogalder",
    "scale(ano_mean_cover)" ~ "ANO: dekningsgrad karpl.",
    "scale(ano_mean_no_spec)" ~ "ANO: artsantall karpl.",
    "scale(ar5_PC1)" ~ "AR5-arealtype: PCA-akse 1",
    "scale(ar5_PC2)" ~ "AR5-arealtype: PCA-akse 2",
    "scale(skog_PC1)" ~ "AR5-skog: PCA-akse 1",
    "scale(skog_PC2)" ~ "AR5-skog: PCA-akse 2",
    "scale(lsm_l_shdi)" ~ "AR5-arealtype: Shannon indeks",
    "scale(tot_perimeter_length)" ~ "AR5-arealtype: lineære elementer",
    "region_nameSørlandet" ~ "Sørlandet",
    "region_nameTrøndelag" ~ "Trøndelag",
    "region_nameNord-Norge" ~ "Nord-Norge"
  )) %>%
  filter(!(Parameter %in% c(
    "Intercept",
    "Temperatur",
    "Lys",
    "Variasjon_temp.",
    "Temp.:lys",
    "Var.temp:Var.lys",
    "Fukt:lys",
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


## ----tab:seminat-mod, results = "asis"------------------------------------------------------------------------------------------------
if (names(rmarkdown::metadata$format) == "ninareport-pdf") {
  xtable(to_tab, digits = 3)
} else {
  kable(to_tab, digits = 3)
}


## -------------------------------------------------------------------------------------------------------------------------------------
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
    "dom_treeHengebjørk" = "Dom. treslag: Hengebjørk",
    "dom_treeGran" = "Dom. treslag: Gran",
    "dom_treeFuru" = "Dom. treslag: Furu",
    "dom_treeSelje" = "Dom. treslag: Selje",
    "dom_treeDunbjørk" = "Dom. treslag: Dunbjørk",
    "scale(avg_age)" = "Landsskog: skogalder",
    "scale(ano_mean_cover)" = "ANO: dekningsgrad karpl.",
    "scale(ano_mean_no_spec)" = "ANO: artsantall karpl.",
    "scale(ar5_PC1)" = "AR5-arealtype: PCA-akse 1",
    "scale(ar5_PC2)" = "AR5-arealtype: PCA-akse 2",
    "scale(skog_PC1)" = "AR5-skog: PCA-akse 1",
    "scale(skog_PC2)" = "AR5-skog: PCA-akse 2",
    "scale(lsm_l_shdi)" = "AR5-arealtype: Shannon indeks",
    "scale(tot_perimeter_length)" = "AR5-arealtype: lineære elementer",
    "region_nameTrøndelag" = "Trøndelag",
    "region_nameNord-Norge" = "Nord-Norge",
    "region_nameSørlandet" = "Sørlandet"
  ),
  title = ""
)


## -------------------------------------------------------------------------------------------------------------------------------------
r2_verdi_biomasse_big_semi_nat <- r.squaredGLMM(global_model_big_semi_nat)


## -------------------------------------------------------------------------------------------------------------------------------------
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
AND yl.year = '2023'
AND yl.project_short_name = 'NasIns'
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


## -------------------------------------------------------------------------------------------------------------------------------------
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


## -------------------------------------------------------------------------------------------------------------------------------------
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


## -------------------------------------------------------------------------------------------------------------------------------------
ord_malaise_ost <- vegan::metaMDS(otu_loc_ost,
  distance = "bray",
  # trymax = 1000
  trymax = 200
)

ord_malaise_ost_aug <- cbind(otu_loc_ost_wider[c(1:6)], vegan::scores(ord_malaise_ost)[1])
# %>%
# mutate(sampling_number = factor(sampling_number),
#       habitat_type = factor(habitat_type))


## -------------------------------------------------------------------------------------------------------------------------------------
ord_malaise_semi <- vegan::metaMDS(otu_loc_semi,
  distance = "bray",
  # trymax = 1000
  trymax = 800
)

ord_malaise_semi_aug <- cbind(otu_loc_semi_wider[c(1:6)], vegan::scores(ord_malaise_semi)[1])
# %>%
# mutate(sampling_number = factor(sampling_number),
#       habitat_type = factor(habitat_type))


## -------------------------------------------------------------------------------------------------------------------------------------
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


## -------------------------------------------------------------------------------------------------------------------------------------
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


## -------------------------------------------------------------------------------------------------------------------------------------
# otu_loc_wider_semi_nat <- otu_loc_raw %>%
#   filter(habitat_type == "Semi-nat") %>%
#   select(-"sampling_number") %>%
#   pivot_wider(names_from = otu,
#               values_from = no_reads,
#               values_fill = 0,) %>%
#   as.data.frame()
#
# rownames(otu_loc_wider_semi_nat) <- otu_loc_wider_semi_nat$sampling_name
# otu_loc_semi_nat <- otu_loc_wider_semi_nat[-c(1:6)]
#
# ord_semi_nat <- vegan::metaMDS(otu_loc_semi_nat,
#                               distance = "bray",
#                               trymax = 100)
#
# ord_semi_nat_aug <- cbind(otu_loc_wider_semi_nat[c(1:6)], vegan::scores(ord_semi_nat)[1])
#
# points <- data.frame(ord_semi_nat$points)

points <- ord_malaise_semi_aug %>%
  left_join(locality_sampling_biomass_big_join_semi_nat,
    by = c(
      "sampling_name" = "sampling_name",
      "locality" = "locality",
      "region_name" = "region_name",
      "habitat_type" = "habitat_type"
    )
  ) %>%
  mutate(region_habitat = paste0(region_name, " - ", habitat_type))

# fit <- vegan::envfit(ord_malaise_semi ~ avg_values_MX2301A_temperature + avg_values_MX2202_light + ar5_PC1 + ar5_PC2,
#                     data = points,
#                     na.rm = TRUE)

fit <- vegan::envfit(ord_malaise_semi ~ scale_avg_values_MX2301A_temperature + scale_avg_values_MX2202_light + scale_avg_values_MX2301A_rh + lsm_l_shdi + tot_perimeter_length + skog_PC1 + skog_PC2 + ar5_PC1 + ar5_PC2 + ano_mean_no_spec + ano_mean_cover,
  data = points,
  na.rm = TRUE
)

fit_scores <- data.frame(vegan::scores(fit, display = "vectors"))
fit_scores <- cbind(fit_scores, vars = rownames(fit_scores))


## -------------------------------------------------------------------------------------------------------------------------------------
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
    legend.text = element_text(size = 7),
    legend.key.width = unit(0.5, "cm")
  ) +
  xlab("NMDS1") +
  ylab("NMDS2")

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


## -------------------------------------------------------------------------------------------------------------------------------------
# otu_loc_wider_skog <- otu_loc_raw %>%
#   filter(habitat_type == "Skog") %>%
#   select(-"sampling_number") %>%
#   pivot_wider(names_from = otu,
#               values_from = no_reads,
#               values_fill = 0,) %>%
#   as.data.frame()
#
# rownames(otu_loc_wider_skog) <- otu_loc_wider_skog$sampling_name
# otu_loc_skog <- otu_loc_wider_skog[-c(1:6)]
#
# ord_skog <- vegan::metaMDS(otu_loc_skog,
#                               distance = "bray",
#                               trymax = 100)
#
# ord_skog_aug <- cbind(otu_loc_wider_skog[c(1:6)], vegan::scores(ord_skog)[1])
#
# points <- data.frame(ord_skog$points)



points <- ord_malaise_ost_aug %>%
  left_join(locality_sampling_biomass_big_join_forest,
    by = c(
      "sampling_name" = "sampling_name",
      "locality" = "locality",
      "region_name" = "region_name",
      "habitat_type" = "habitat_type"
    )
  ) %>%
  mutate(region_habitat = paste0(region_name, " - ", habitat_type))

# fit <- vegan::envfit(ord_malaise_ost ~ avg_values_MX2301A_temperature + avg_values_MX2202_light + ar5_PC1 + ar5_PC2,
#                      data = points,
#                      na.rm = TRUE)

fit <- vegan::envfit(ord_malaise_ost ~ scale_avg_values_MX2301A_temperature + scale_avg_values_MX2202_light + scale_avg_values_MX2301A_rh + lsm_l_shdi + tot_perimeter_length + skog_PC1 + skog_PC2 + ar5_PC1 + ar5_PC2 + ano_mean_no_spec + ano_mean_cover,
  data = points,
  na.rm = TRUE
)



fit_scores <- data.frame(vegan::scores(fit, display = "vectors"))
fit_scores <- cbind(fit_scores, vars = rownames(fit_scores))


## -------------------------------------------------------------------------------------------------------------------------------------
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
    legend.text = element_text(size = 7),
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


## -------------------------------------------------------------------------------------------------------------------------------------
fig_habitats_ost


## -------------------------------------------------------------------------------------------------------------------------------------
fig_seminat_mds


## -------------------------------------------------------------------------------------------------------------------------------------
fig_seminat


## -------------------------------------------------------------------------------------------------------------------------------------
fig_ost


## -------------------------------------------------------------------------------------------------------------------------------------
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
AND yl.project_short_name = 'NasIns'
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
  mutate(expl_3_main = str_trim(expl_3_main)) %>%
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


## -------------------------------------------------------------------------------------------------------------------------------------
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
AND yl.project_short_name = 'NasIns'
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
  mutate(expl_3_main = str_trim(expl_3_main)) %>%
  mutate(
    expl_3_main = ifelse(expl_3_main == " ", NA, expl_3_main),
    expl_3_main = ifelse(expl_3_main == "", NA, expl_3_main),
    expl_3_main = ifelse(expl_3_main == "Ukjent ", NA, expl_3_main),
    expl_3_main = ifelse(expl_3_main == "Ukjent", NA, expl_3_main)
  ) %>%
  filter(kategori_2021 != "DD")


## -------------------------------------------------------------------------------------------------------------------------------------
no_detected_redlist_species <- redlisted_obs_2021 %>%
  st_drop_geometry() %>%
  filter(kategori_2021 %in% c(
    "NT",
    "VU",
    "EN",
    "CR",
    "RE"
  )) %>%
  summarise(no_spec = n_distinct(vitenskapelig_navn_id)) %>%
  pull()


## -------------------------------------------------------------------------------------------------------------------------------------
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


## -------------------------------------------------------------------------------------------------------------------------------------
no_locales_with_redlisted_species <- redlisted_obs_2021_agg %>%
  st_drop_geometry() %>%
  ungroup() %>%
  summarise(n_distinct(locality)) %>%
  pull()


## -------------------------------------------------------------------------------------------------------------------------------------
# redlist_colors <- c("#5a5b5d","#990013", "#d80024" , "#ff4c29", "#f36c21", "#f9a61a")
# names(redlist_colors) <- c("RE", "CR", "EN", "VU", "NT", "DD")

# redlist_colors <- c("#f9a61a", "#f36c21", "#ff4c29", "#d80024") #2022 colors
redlist_colors <- c("#ed6c26", "#e94f33", "#d80f27", "#5a5b5d")

# redlist_colors <- c("#f7a520", "#ed6c26", "#d80f27", "#5a5b5d") #2023 powertoys color picker tool

names(redlist_colors) <- levels(redlisted_obs_2021_agg$kategori_2021)


## -------------------------------------------------------------------------------------------------------------------------------------
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
    legend.text = element_text(size = 6),
    legend.key.width = unit(0.5, "cm")
  )


## -------------------------------------------------------------------------------------------------------------------------------------
redlisted_obs_2021_agg_loc_expl <- redlisted_obs_2021 %>%
  group_by(
    locality,
    expl_3_main
  ) %>%
  summarise(no_spec_loc_expl = n_distinct(species_latin_fixed)) %>%
  ungroup() %>%
  mutate(expl_3_main = ifelse(is.na(expl_3_main), "Ukjent/Ikke angitt", expl_3_main))

redlisted_obs_2021_agg_expl <- redlisted_obs_2021 %>%
  st_drop_geometry() %>%
  group_by(expl_3_main) %>%
  summarise(no_spec_expl = n_distinct(species_latin_fixed)) %>%
  ungroup() %>%
  mutate(expl_3_main = ifelse(is.na(expl_3_main), "Ukjent/Ikke angitt", expl_3_main)) %>%
  arrange(desc(no_spec_expl))

kat_order <- tibble(
  expl_3_main = c(
    "Skogbruk (kommersielt)",
    "Jordbruk",
    "Ukjent/Ikke angitt",
    "Utbygging/utvinning",
    "Opphørt/redusert drift",
    "Buskap/dyrehold",
    "Fremmede arter",
    "Opphørt drift",
    "Påvirker habitatet (beite tråkk mm.)"
  ),
  kat_order = 1:9
)

redlisted_obs_2021_agg_expl_big <- redlisted_obs_2021_agg_loc_expl %>%
  left_join(redlisted_obs_2021_agg_expl,
    by = c("expl_3_main" = "expl_3_main")
  ) %>%
  left_join(kat_order,
    by = c("expl_3_main" = "expl_3_main")
  ) %>%
  ungroup() %>%
  mutate(expl_3_main = factor(paste0(expl_3_main, " (", no_spec_expl, " stk.)"))) %>%
  mutate(expl_3_main = fct_reorder(expl_3_main, kat_order, .na_rm = TRUE))

# levels(redlisted_obs_2021_agg_expl$expl_3_main)

redlisted_obs_2021_agg_expl_big <- redlisted_obs_2021_agg_expl_big %>%
  st_jitter(redlisted_obs_2021_agg_expl_big, amount = 7000)


## -------------------------------------------------------------------------------------------------------------------------------------
redlist_expl_colors <- c(ninaColors()[c(1, 6, 2, 3, 5, 4)], "yellow", "red", "green")
names(redlist_expl_colors) <- levels(redlisted_obs_2021_agg_expl_big$expl_3_main)

set.seed(12345)
norge %>%
  ggplot(.) +
  geom_sf() +
  geom_sf(
    aes(
      size = no_spec_loc_expl,
      color = expl_3_main
    ),
    alpha = 0.7,
    data = redlisted_obs_2021_agg_expl_big %>% arrange(expl_3_main)
  ) +
  # scale_color_nina(name = "Påvirkning") +
  scale_color_manual(
    name = "Påvirkning",
    values = redlist_expl_colors
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
  guides(color = guide_legend(override.aes = list(size = 7), nrow = 3, byrow = T)) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.margin = unit(x = c(0, 40, 0, 40), units = "mm"),
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 5),
    legend.key.width = unit(0.2, "cm")
  )


## -------------------------------------------------------------------------------------------------------------------------------------
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


## -------------------------------------------------------------------------------------------------------------------------------------
ggplot(
  redlisted_obs_2021_less_agg,
  aes(
    y = locality,
    x = kategori_2021_append
  )
) +
  geom_point(
    aes(
      size = times_detected,
      color = kategori_2021_append
    ),
    alpha = 0.7,
    position = position_jitter(height = 0)
  ) +
  scale_y_discrete(limits = rev) +
  # guides(color = "none") +
  scale_color_manual(
    name = "Rødlistekategori",
    values = redlist_colors
  ) +
  scale_size(
    name = "Antall ganger\ndetektert",
    breaks = c(1, 2, 4, 6, 8, 10)
  ) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 6)
  ) +
  xlab("Rødlistekategori 2021") +
  ylab("Lokalitet")


## -------------------------------------------------------------------------------------------------------------------------------------

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
    kategori_2021,
    desc(no_lok),
    desc(no_detected)
  ) %>%
  dplyr::select(
    `Rødlistekat.` = kategori_2021,
    `Orden` = id_order,
    `Latinsk navn` = species_latin,
    `Norsk navn` = popular_name,
    `Antall lokaler` = no_lok,
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

# names(spec_no_loc1)=c("\\multicolumn{1}{p{0.4in}}{\\raggedright Kategori}",
#  "\\multicolumn{1}{p{0.75in}}{\\raggedright Orden}",
# "\\multicolumn{1}{p{0.75in}}{\\raggedright Latinsk navn}",
# "\\multicolumn{1}{p{0.75in}}{\\raggedright Norsk navn}",
# "\\multicolumn{1}{p{0.75in}}{\\raggedright Funnet i antall lokaliteter}",
# "\\multicolumn{1}{p{0.75in}}{\\raggedright Antall ganger observert}")

# names(spec_no_loc2)=c("\\multicolumn{1}{p{0.4in}}{\\raggedright Kategori}",
#  "\\multicolumn{1}{p{0.75in}}{\\raggedright Orden}",
# "\\multicolumn{1}{p{0.75in}}{\\raggedright Latinsk navn}",
# "\\multicolumn{1}{p{0.75in}}{\\raggedright Norsk navn}",
# "\\multicolumn{1}{p{0.75in}}{\\raggedright Funnet i antall lokaliteter}",
# "\\multicolumn{1}{p{0.75in}}{\\raggedright Antall ganger observert}")


## -------------------------------------------------------------------------------------------------------------------------------------
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
AND yl.project_short_name = 'NasIns'
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


## -------------------------------------------------------------------------------------------------------------------------------------
no_fennoscand_q <- "
SELECT obs.id_phylum, obs.id_class, obs.id_order, count(distinct obs.species_latin_fixed)
FROM occurrences.observations obs,
lookup.horizon_scan_2023 h,
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
AND yl.project_short_name = 'NasIns'
AND obs.identification_confidence = 'HIGH'
AND h.fennoscandian IS TRUE
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


## -------------------------------------------------------------------------------------------------------------------------------------
no_pot_alien_q <- "
SELECT obs.id_phylum, obs.id_class, obs.id_order, count(distinct obs.species_latin_fixed)
FROM occurrences.observations obs,
lookup.horizon_scan_2023 alien,
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
AND yl.project_short_name = 'NasIns'
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
  )


## -------------------------------------------------------------------------------------------------------------------------------------
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
AND yl.project_short_name = 'NasIns'
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


## -------------------------------------------------------------------------------------------------------------------------------------
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
AND yl.project_short_name = 'NasIns'
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


## -------------------------------------------------------------------------------------------------------------------------------------
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
AND yl.project_short_name = 'NasIns'
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


## -------------------------------------------------------------------------------------------------------------------------------------
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
  mutate(id_order = ifelse(id_order == 0, "Ukjent", id_order))


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


## -------------------------------------------------------------------------------------------------------------------------------------
xtab <- xtable(no_all_aliens)
align(xtab) <- "p{1.5cm}p{1.2cm}p{1.2cm}p{2.2cm}p{1cm}p{1cm}p{1cm}p{1cm}p{1.7cm}p{1cm}"

if (names(rmarkdown::metadata$format) == "ninareport-pdf") {
  print(xtab, sanitize.text.function = function(x) {
    x
  }, size = "scriptsize")
} else {
  kable(xtab)
}


## -------------------------------------------------------------------------------------------------------------------------------------
pot_alien_spec <- tbl(
  con,
  Id(
    schema = "lookup",
    table = "horizon_scan_2023"
  )
) %>%
  filter(potential_alien_species) %>%
  select(species_latin_fixed) %>%
  pull()

fennoscand_spec <- tbl(
  con,
  Id(
    schema = "lookup",
    table = "horizon_scan_2023"
  )
) %>%
  filter(fennoscandian) %>%
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


## -------------------------------------------------------------------------------------------------------------------------------------
fennoskand_obs_q <- "
SELECT DISTINCT ON(species_latin_fixed, year, sampling_name, locality, geom)
obs.species_latin_fixed,
yl.year,
ls.sampling_name,
l.locality,
ST_Centroid(l.geom) as geom
FROM occurrences.observations obs,
lookup.horizon_scan_2023 fennoscand,
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
AND yl.project_short_name = 'NasIns'
AND obs.identification_confidence = 'HIGH'
AND fennoscand.fennoscandian IS TRUE
"

fennoskand_obs <- read_sf(con,
  query = fennoskand_obs_q
) %>%
  mutate(alien_kat = "Fennoskandisk forek.")

## -------------------------------------------------------------------------------------------------------------------------------------
pot_alien_obs_q <- "
SELECT DISTINCT ON (obs.species_latin_fixed, year, sampling_name, locality, geom)
obs.species_latin_fixed,
yl.year,
ls.sampling_name,
l.locality,
ST_Centroid(l.geom) as geom
FROM occurrences.observations obs,
lookup.horizon_scan_2023 pot_alien,
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
AND yl.project_short_name = 'NasIns'
AND obs.identification_confidence = 'HIGH'
AND pot_alien.potential_alien_species IS TRUE
"

pot_alien_obs <- read_sf(con,
  query = pot_alien_obs_q
) %>%
  mutate(alien_kat = "Potensielt fremmede arter")


## -------------------------------------------------------------------------------------------------------------------------------------
alien_obs_q <- "
SELECT
obs.species_latin_fixed,
yl.year,
ls.sampling_name,
l.locality,
ST_Centroid(l.geom) as geom
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
AND yl.project_short_name = 'NasIns'
AND alien.risikokategori_2023 IN ('SE', 'HI', 'PH', 'LO')
AND obs.identification_confidence = 'HIGH'
"

alien_obs <- read_sf(con,
  query = alien_obs_q
) %>%
  mutate(alien_kat = "På fremmedartslista")


## -------------------------------------------------------------------------------------------------------------------------------------
all_alien_obs <- fennoskand_obs %>%
  rbind(pot_alien_obs) %>%
  rbind(alien_obs)


## -------------------------------------------------------------------------------------------------------------------------------------
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


## -------------------------------------------------------------------------------------------------------------------------------------
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


## -------------------------------------------------------------------------------------------------------------------------------------
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
  )

# %>%
# mutate(kategori_2021 = factor(kategori_2021, levels = c("DD", "NT", "VU", "EN")))

# redlisted_obs_2021 %>%
#   filter(locality == "Skog_10",
#          species_latin_fixed == "Strophosoma fulvicorne")
#

# redlisted_obs_2021_less_agg %>%
# print(n = Inf)


## -------------------------------------------------------------------------------------------------------------------------------------
ggplot(
  all_alien_obs_less_agg,
  aes(
    y = locality,
    x = alien_kat
  )
) +
  geom_point(
    aes(
      size = times_detected,
      color = alien_kat_append
    ),
    alpha = 0.5,
    position = position_jitter(height = 0)
  ) +
  scale_y_discrete(limits = rev) +
  # guides(color = "none") +
  scale_color_nina(
    name = "Kategori",
    palette = "contrast"
  ) +
  scale_size(
    name = "Antall ganger\ndetektert",
    breaks = c(1, 2, 4, 6, 8, 10)
  ) +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 6),
    axis.text.x = element_text(
      angle = 20,
      vjust = 0.6
    )
  ) +
  xlab("") +
  ylab("Lokalitet")




## -------------------------------------------------------------------------------------------------------------------------------------
load(file = "out/region_habitat_div.rda")


## -------------------------------------------------------------------------------------------------------------------------------------
to_tab <- region_habitat_div %>%
  arrange(
    År,
    Region,
    `Habitat-type`
  )

if (names(rmarkdown::metadata$format) == "ninareport-pdf") {
  xtable(to_tab)
} else {
  kable(to_tab)
}


## -------------------------------------------------------------------------------------------------------------------------------------
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




## -------------------------------------------------------------------------------------------------------------------------------------
load(file = "out/ost_iNEXT.rda")


## -------------------------------------------------------------------------------------------------------------------------------------
lines <- c(
  "Rarefaction" = "solid",
  "Extrapolation" = "dotted"
)

ggiNEXT(ost_iNEXT) +
  ylab("Artsantall") +
  xlab("Antall lokaliteter") +
  theme(
    legend.text = element_text(size = 8),
    legend.key.width = unit(1, "cm"),
    axis.title = element_text(size = 14)
  ) +
  scale_color_manual(values = c("#004F71", "#93328E", "#7A9A01", "#E57200", "#008C95")) +
  guides(
    color = guide_legend("", nrow = 2, byrow = TRUE),
    shape = guide_legend("", nrow = 2, byrow = TRUE),
    fill = guide_legend("", nrow = 2, byrow = TRUE)
  ) +
  scale_linetype_manual(
    values = lines,
    name = "Method",
    labels = c("Interpolert", "Ekstrapolert")
  )


## ---- cache = cache_var---------------------------------------------------------------------------------------------------------------
spec_occ <- tbl(
  con,
  Id(
    schema = "views",
    table = "species_occurrences"
  )
) %>%
  filter(project_short_name == "NasIns")

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


## -------------------------------------------------------------------------------------------------------------------------------------
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


## -------------------------------------------------------------------------------------------------------------------------------------
all_comm_mat <- get_community_matrix(
  trap_type = "All",
  dataset = "NasIns",
  as_tibble = T,
  subset_region = NULL
)


## -------------------------------------------------------------------------------------------------------------------------------------
all_comm_mat_mat <- all_comm_mat %>%
  select(-c(year, locality))

levels <- tibble(
  l1 = seq(1:nrow(all_comm_mat_mat)),
  l2 = rep(1, nrow(all_comm_mat_mat))
)




## -------------------------------------------------------------------------------------------------------------------------------------
load(file = "out/div_part.rds")


## -------------------------------------------------------------------------------------------------------------------------------------
div_part_tab <- tibble(
  Nivå = c(
    "$\\alpha$",
    "$\\beta$",
    "$\\gamma$"
  ),
  "Målt verdi" = div_part$oecosimu$statistic[c(1, 3, 2)],
  "Stand. Effekt" = div_part$oecosimu$z[c(1, 3, 2)],
  "Forventet verdi" = div_part$oecosimu$means[c(1, 3, 2)],
  "Prob (sim)" = div_part$oecosimu$pval[c(1, 3, 2)]
)

if (names(rmarkdown::metadata$format) == "ninareport-pdf") {
  print(xtable(div_part_tab, digits = 3),
    sanitize.text.function = function(x) as.character(x),
    sanitise.rownames.function = function(x) as.character(x),
    include.rownames = F
  )
} else {
  kable(div_part_tab,
    digits = 3
  )
}


## -------------------------------------------------------------------------------------------------------------------------------------
div_betacore <- betapart.core(all_comm_mat_mat)


## -------------------------------------------------------------------------------------------------------------------------------------
all_beta_part_freq <- beta.multi(div_betacore) %>%
  unlist()

all_beta_part_freq <- as.data.frame(t(round(all_beta_part_freq / all_beta_part_freq["beta.SOR"] * 100, 2)))
names(all_beta_part_freq) <- c("$\\beta_{SIM}\\%$", "$\\beta_{NES}\\%$", "$\\beta_{SOR}\\%$")


## -------------------------------------------------------------------------------------------------------------------------------------
if (names(rmarkdown::metadata$format) == "ninareport-pdf") {
  print(xtable(all_beta_part_freq), sanitize.colnames.function = function(x) as.character(x))
} else {
  kable(all_beta_part_freq)
}


## -------------------------------------------------------------------------------------------------------------------------------------
skog_ost_comm_mat <- get_community_matrix( # trap_type = "MF",
  dataset = "NasIns",
  subset_habitat = "Forest",
  subset_region = "Østlandet",
  as_tibble = T,
  transposed_matrix = F
)
# %>%
#  select(-c(year, locality))

semi_nat_comm_mat <- get_community_matrix(
  trap_type = "MF",
  dataset = "NasIns",
  subset_habitat = "Semi-nat",
  subset_region = NULL,
  as_tibble = T,
  transposed_matrix = F
)
# %>%
#  select(-c(year, locality))


## -------------------------------------------------------------------------------------------------------------------------------------
skog_ost_dist_beta <- combine_dist_to_comm_mat(
  comm_mat = skog_ost_comm_mat,
  region_name = "('Østlandet')",
  habitat_type = "Forest"
)


## -------------------------------------------------------------------------------------------------------------------------------------
## beta_sne_lm <- glm(beta_sne ~ distance, data = skog_ost_dist_beta)
## summary(beta_sne_lm)
## beta_sim_lm <- glm(beta_sim ~ distance, data = skog_ost_dist_beta)
## summary(beta_sim_lm)
## beta_sor_lm <- glm(beta_sor ~ distance, data = skog_ost_dist_beta)
## summary(beta_sor_lm)
## 
## with(summary(beta_sor_lm), 1 - deviance / null.deviance)


## -------------------------------------------------------------------------------------------------------------------------------------
plot_beta_part(skog_ost_dist_beta)


## ---- cache = cache_var---------------------------------------------------------------------------------------------------------------
semi_nat_dist_beta <- combine_dist_to_comm_mat(
  comm_mat = semi_nat_comm_mat,
  region_name = "('Østlandet', 'Trøndelag', 'Sørlandet', 'Nord-Norge')",
  habitat_type = "Semi-nat"
)


## -------------------------------------------------------------------------------------------------------------------------------------
plot_beta_part(semi_nat_dist_beta)


## -------------------------------------------------------------------------------------------------------------------------------------
## skog_ost_comm_mat_lepi <- get_community_matrix(
##   trap_type = "MF",
##   dataset = "NasIns",
##   subset_habitat = "Forest",
##   subset_region = "Østlandet",
##   subset_orders = "Lepidoptera",
##   as_tibble = T,
##   transposed_matrix = F,
##   obs_at_least_x_times_somewhere = 2
## )
## # %>%
## #  select(-c(year, locality))
## 
## semi_nat_ost_comm_mat_lepi <- get_community_matrix(
##   trap_type = "MF",
##   dataset = "NasIns",
##   subset_habitat = "Semi-nat",
##   subset_region = "Østlandet",
##   subset_orders = "Lepidoptera",
##   as_tibble = T,
##   transposed_matrix = F,
##   obs_at_least_x_times_somewhere = 2
## )
## # %>%
## #  select(-c(year, locality))


## -------------------------------------------------------------------------------------------------------------------------------------
## skog_ost_dist_beta_lepi <- combine_dist_to_comm_mat(
##   comm_mat = skog_ost_comm_mat_lepi,
##   region_name = "('Østlandet')",
##   habitat_type = "Forest"
## )
## 
## semi_nat_ost_dist_beta_lepi <- combine_dist_to_comm_mat(
##   comm_mat = semi_nat_ost_comm_mat_lepi,
##   region_name = "('Østlandet')",
##   habitat_type = "Semi-nat"
## )


## -------------------------------------------------------------------------------------------------------------------------------------
## skog_ost_dist_beta_lepi_mat <- skog_ost_comm_mat_lepi %>%
##   select(-c(year, locality))
## 
## levels <- tibble(
##   l1 = seq(1:nrow(skog_ost_dist_beta_lepi_mat)),
##   l2 = rep(1, nrow(skog_ost_dist_beta_lepi_mat))
## )


## -------------------------------------------------------------------------------------------------------------------------------------
## div_part_skog_lepi <- multipart(
##   y = skog_ost_dist_beta_lepi_mat,
##   x = levels,
##   global = T,
##   nsimul = 99,
##   scales = 0
## )
## 
## save(div_part_skog_lepi,
##   file = "out/div_part_skog_lepi.rds"
## )


## -------------------------------------------------------------------------------------------------------------------------------------
## div_betacore_skog_lepi <- betapart.core(skog_ost_dist_beta_lepi_mat)


## -------------------------------------------------------------------------------------------------------------------------------------
## all_beta_part_freq_skog_lepi <- beta.multi(div_betacore_skog_lepi) %>%
##   unlist()
## 
## all_beta_part_freq_skog_lepi <- as.data.frame(t(round(all_beta_part_freq_skog_lepi / all_beta_part_freq_skog_lepi["beta.SOR"] * 100, 2)))
## names(all_beta_part_freq) <- c("$\\beta_{SIM}\\%$", "$\\beta_{NES}\\%$", "$\\beta_{SOR}\\%$")


## -------------------------------------------------------------------------------------------------------------------------------------
## plot_beta_part(skog_ost_dist_beta_lepi)


## -------------------------------------------------------------------------------------------------------------------------------------
## plot_beta_part(semi_nat_ost_dist_beta_lepi)


## -------------------------------------------------------------------------------------------------------------------------------------
## most_abund_dna_q <- "
## 
## SELECT *
## FROM
## (
## SELECT id_order, sum, species_latin_fixed, rank()
## OVER (PARTITION BY id_order ORDER BY sum DESC) as pos
## FROM
## 
## (
## SELECT id_order, o.species_latin_fixed, sum(o.no_reads)
## 
## FROM occurrences.observations o,
## events.identifications i,
## events.sampling_trap st,
## events.locality_sampling ls,
## events.year_locality yl
## WHERE o.identification_id = i.id
## AND i.sampling_trap_id = st.id
## AND st.locality_sampling_id = ls.id
## AND ls.year_locality_id = yl.id
## AND project_short_name = 'NasIns'
## AND id_order IN ('Lepidoptera','Coleoptera', 'Diptera', 'Hymenoptera')
## GROUP BY id_order, species_latin_fixed
## ) foo
## )bar
## 
## WHERE pos <4
## 
## "
## 
## most_abund_dna <- dbGetQuery(
##   con,
##   most_abund_dna_q
## )


## -------------------------------------------------------------------------------------------------------------------------------------
## most_common_q <- "
## SELECT * FROM
## (SELECT id_order, species_latin_fixed, n_occ, rank() OVER (PARTITION BY id_order ORDER BY n_occ DESC) as pos
## FROM
## (
## SELECT id_order, species_latin_fixed, count(*) as n_occ
## FROM
## (SELECT distinct on(id_order, species_latin_fixed, st.id) id_order, st.sample_name, o.species_latin_fixed
## FROM occurrences.observations o,
## events.identifications i,
## events.sampling_trap st,
## events.locality_sampling ls,
## events.year_locality yl
## WHERE o.identification_id = i.id
## AND i.sampling_trap_id = st.id
## AND st.locality_sampling_id = ls.id
## AND ls.year_locality_id = yl.id
## AND project_short_name = 'NasIns'
## AND o.identification_confidence = 'HIGH'
## AND id_order IN ('Lepidoptera','Coleoptera', 'Diptera', 'Hymenoptera')
## ORDER BY id_order, species_latin_fixed, st.id) foo
## GROUP BY id_order, species_latin_fixed) bar
## ) last
## WHERE pos < 6
## 
## "
## 
## most_common <- dbGetQuery(
##   con,
##   most_common_q
## )






## -------------------------------------------------------------------------------------------------------------------------------------
include_graphics("figures/asv_composite.png")


## -------------------------------------------------------------------------------------------------------------------------------------
## most_common_poll_q <- "
## 
## SELECT * FROM
## (SELECT id_family, species_latin_fixed, n_occ, rank() OVER (PARTITION BY id_family ORDER BY n_occ DESC) as pos
## FROM
## (
## SELECT id_family, species_latin_fixed, count(*) as n_occ
## FROM
## (SELECT distinct on(id_family, species_latin_fixed, st.id) id_family, st.sample_name, o.species_latin_fixed
## FROM occurrences.observations o,
## events.identifications i,
## events.sampling_trap st,
## events.locality_sampling ls,
## events.year_locality yl
## WHERE o.identification_id = i.id
## AND i.sampling_trap_id = st.id
## AND st.locality_sampling_id = ls.id
## AND ls.year_locality_id = yl.id
## AND project_short_name = 'NasIns'
## AND o.identification_confidence = 'HIGH'
## AND id_family IN ('Apidae','Nymphalidae','Lycaenidae', 'Syrphidae')
## ORDER BY id_family, species_latin_fixed, st.id) foo
## GROUP BY id_family, species_latin_fixed) bar
## ) last
## WHERE pos < 6
## 
## "
## 
## most_common_poll <- dbGetQuery(
##   con,
##   most_common_poll_q
## )






## -------------------------------------------------------------------------------------------------------------------------------------
include_graphics("figures/asv_poll_composite.png")


## -------------------------------------------------------------------------------------------------------------------------------------
yl_catch_tab <- year_locality_biodiv %>%
  mutate(year = as.integer(year)) %>%
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


## -------------------------------------------------------------------------------------------------------------------------------------
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


## -------------------------------------------------------------------------------------------------------------------------------------
obs <- tbl(
  con,
  Id(
    schema = "occurrences",
    table = "observations"
  )
) %>%
  filter(identification_confidence == "HIGH")


## -------------------------------------------------------------------------------------------------------------------------------------
pot_alien_spec <- tbl(
  con,
  Id(
    schema = "lookup",
    table = "pot_alien_species2"
  )
) %>%
  mutate(
    saul = as.character(saul),
    daisie = as.character(daisie)
  ) %>%
  mutate(
    saul = ifelse(saul == "false", "", "X"),
    daisie = ifelse(daisie == "false", "", "X")
  ) %>%
  select(
    Artsnavn = species_latin_fixed,
    "Listet i Saul" = saul,
    "Listet i Daisie" = daisie
  ) %>%
  arrange(Artsnavn) %>%
  collect()

pot_alien_spec <- pot_alien_spec %>%
  inner_join(obs,
    by = c("Artsnavn" = "species_latin_fixed"),
    copy = T
  ) %>%
  select(
    Order = "id_order",
    Artsnavn,
    "Listet i Saul",
    "Listet i Daisie"
  ) %>%
  distinct()

fennoscand_spec <- tbl(
  con,
  Id(
    schema = "lookup",
    table = "fennoscand_species2"
  )
) %>%
  mutate(
    saul = as.character(saul),
    daisie = as.character(daisie)
  ) %>%
  mutate(
    saul = ifelse(saul == "false", "", "X"),
    daisie = ifelse(daisie == "false", "", "X")
  ) %>%
  select(
    Artsnavn = species_latin_fixed,
    "Listet i Saul" = saul,
    "Listet i Daisie" = daisie
  ) %>%
  arrange(Artsnavn) %>%
  collect()

fennoscand_spec <- fennoscand_spec %>%
  inner_join(obs,
    by = c("Artsnavn" = "species_latin_fixed"),
    copy = T
  ) %>%
  select(
    Order = "id_order",
    Artsnavn,
    "Listet i Saul",
    "Listet i Daisie"
  ) %>%
  distinct()


## -------------------------------------------------------------------------------------------------------------------------------------
alien_obs_q <- "
SELECT
obs.id_order,
obs.species_latin_fixed,
alien.\"riskCategory\",
l.habitat_type
FROM occurrences.observations obs,
lookup.fremmedartslista_2018_artsdatabanken alien,
events.identifications,
events.year_locality yl,
events.locality_sampling ls,
locations.localities l,
events.sampling_trap st
WHERE obs.identification_id = identifications.id
AND obs.species_latin_fixed = alien.\"scientificName\"
AND identifications.sampling_trap_id = st.id
AND st.locality_sampling_id = ls.id
AND ls.year_locality_id = yl.id
AND yl.locality_id = l.id
AND yl.project_short_name = 'NasIns'
AND alien.\"riskCategory\" IN ('SE', 'HI', 'PH', 'LO')
AND obs.identification_confidence = 'HIGH'
GROUP BY obs.id_order, species_latin_fixed, \"riskCategory\", habitat_type
ORDER BY  \"riskCategory\", species_latin_fixed, habitat_type
"

alien_obs <- dbGetQuery(
  con,
  alien_obs_q
) %>%
  mutate(habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type)) %>%
  select(
    Orden = id_order,
    Artsnavn = species_latin_fixed,
    Risikovurdering = riskCategory,
    Habitatstype = habitat_type
  ) %>%
  arrange(Artsnavn)


## -------------------------------------------------------------------------------------------------------------------------------------
# names(spec_no_loc1) <- c("Kategori",
#                          "Orden",
#                          "Latinsk navn",
#                          "Norsk navn",
#                          "Funnet i antall lokaliteter",
#                          "Antall ganger observert")
# spec_no_loc1 <- head(spec_no_loc1)
# xtable(spec_no_loc1)

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


## -------------------------------------------------------------------------------------------------------------------------------------
# names(spec_no_loc2) <- c("Kategori",
#                          "Orden",
#                          "Latinsk navn",
#                          "Norsk navn",
#                          "Funnet i antall lokaliteter",
#                          "Antall ganger observert")
# spec_no_loc2 <- head(spec_no_loc2)
#
# xtable(spec_no_loc2)

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


## -------------------------------------------------------------------------------------------------------------------------------------
if (names(rmarkdown::metadata$format) == "ninareport-pdf") {
  xtable(alien_obs)
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

