## ----load_packages, include = F-----------------------------------------------------------------------------------------------------------------
require(NinaR)
require(Norimon)
require(DBI)
require(RPostgres)
require(ggplot2)
require(xtable)
require(sf)
require(ggplot2)
require(ggforce)
#require(bookdown)
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
require(MuMIn)
require(iNEXT)
require(latex2exp)
require(ggvegan)
require(bookdown)
require(phyloseq)
require(Polychrome)
require(ggthemes)
suppressPackageStartupMessages(require(tidyverse))



## ----setup, echo=F, message=F, warning=F--------------------------------------------------------------------------------------------------------

opts_chunk$set(tidy=T, 
               tidy.opts=list(width.cutoff=50),
               echo = F,
               message = F,
               warning = F,
               dpi = 600,
               dev = c("pdf", "png", "svg"),
               fig.path = "figure/",
               cache.extra = 3 #rename this to clear cache
               )
options(width = 60,
        "xtable.comment" = F,
        "xtable.include.rownames" = F)


theme_set(theme_few())



## -----------------------------------------------------------------------------------------------------------------------------------------------
#Set up rerun cache for faster renders
cache_var = TRUE



## -----------------------------------------------------------------------------------------------------------------------------------------------
connect_to_insect_db()


## -----------------------------------------------------------------------------------------------------------------------------------------------
#typo, update Norimon
#c("identification_name", "identification_name") ---> c("identification_name" = "identification_name")
get_observations <- function(id_type = c("metabarcoding"),
         subset_orders = NULL,
         subset_families = NULL,
         subset_genus = NULL,
         subset_species = NULL,
         subset_year = NULL,
         subset_region = NULL,
         subset_habitat = NULL,
         trap_type = "All",
         limit = NULL,
         dataset = "NasIns",
         agg_level = "year_locality",
         digits = 2,
         as_tibble = F){
  
  #Bind these variables to stop R CMD check complaints
  if(!exists("con")) {con <- NULL}
  
  Norimon:::checkCon()
  
  if(!is.null(subset_region)){
    subset_region <- match.arg(subset_region, choices = c("Østlandet", "Trøndelag"))
  }
  
  if(!is.null(subset_habitat)){
    subset_habitat <- match.arg(subset_habitat, choices = c("Forest", "Semi-nat"))
  }
  
  id_type <- match.arg(id_type, choices = c("metabarcoding"))
  dataset <- match.arg(dataset, choices = c("NasIns",
                                            "OkoTrond",
                                            "TidVar",
                                            "Nerlandsøya"))
  
  agg_level <- match.arg(agg_level, choices = c("year_locality",
                                                "locality_sampling",
                                                "region_habitat",
                                                "region_habitat_year",
                                                "total",
                                                "none"))
  
  trap_type <- match.arg(trap_type,
                         choices = c("All", "MF", "VF", NULL))
  
  
  ##Set up table sources
  ##Probably needs updating after new batch of data. Also need to test filtering of different identification types
  observations <- dplyr::tbl(con, dbplyr::in_schema("occurrences", "observations"))
  identifications <- dplyr::tbl(con, dbplyr::in_schema("events", "identifications"))
  sampling_trap <- dplyr::tbl(con, dbplyr::in_schema("events", "sampling_trap"))
  locality_sampling <- dplyr::tbl(con, dbplyr::in_schema("events", "locality_sampling"))
  year_locality <- dplyr::tbl(con, dbplyr::in_schema("events", "year_locality"))
  localities <- dplyr::tbl(con, dbplyr::in_schema("locations", "localities"))
  identification_techniques <- dplyr::tbl(con, dbplyr::in_schema("lookup", "identification_techniques"))
  traps <- dplyr::tbl(con, dbplyr::in_schema("locations", "traps"))
  
  ##Join the tables
  
  joined <- observations %>%
    left_join(identifications,
              by = c("identification_id" = "id"),
              suffix = c("_obs", "_ids")) %>%
    left_join(identification_techniques,
              by = c("identification_name" = "identification_name"),
              suffix = c("_obs", "_idtechn")) %>%
    left_join(sampling_trap,
              by = c("sampling_trap_id" = "id"),
              suffix = c("_obs", "_st")) %>%
    left_join(locality_sampling,
              by = c("locality_sampling_id" = "id"),
              suffix = c("_obs", "_ls")) %>%
    left_join(year_locality,
              by = c("year_locality_id" = "id"),
              suffix = c("_obs", "_yl")) %>%
    left_join(localities,
              by = c("locality_id" = "id"),
              suffix = c("_obs", "_loc"))  %>%
    left_join(traps,
              by = c("trap_id" = "id",
                     "year" = "year",
                     "locality" = "locality")
    ) %>%
    mutate(year = as.character(year))
  
  
  
  ##Exclude 2020 4 week samplings
  
  joined <-  joined %>%
    mutate(weeks_sampled = ifelse(grepl("2020", year) & (grepl("1", .data$trap_short_name) | grepl("3", trap_short_name)), 2, 4)) %>%
    mutate(weeks_sampled = ifelse(grepl("2020", year), .data$weeks_sampled, 2))
  
  joined <- joined %>%
    filter(.data$weeks_sampled == 2)
  
  
  if(id_type == "metabarcoding"){
    joined <- joined %>%
      filter(.data$identification_type == "metabarcoding")
  }
  
  #Filter on region name
  if(!is.null(subset_region)){
    subset_region <- c("", subset_region)
    joined <- joined %>%
      filter(.data$region_name %in% subset_region)
  }
  
  #Filter on region name
  if(!is.null(subset_habitat)){
    subset_habitat <- c("", subset_habitat)
    joined <- joined %>%
      filter(.data$habitat_type %in% subset_habitat)
  }
  
  if(!is.null(subset_orders)){
    subset_orders <- c("", subset_orders) #To allow one-length subsets
    joined <- joined %>%
      filter(.data$id_order %IN% subset_orders)
  }
  
  if(!is.null(subset_families)){
    subset_families <- c("", subset_families)
    joined <- joined %>%
      filter(.data$id_family %in% subset_families)
  }
  
  if(!is.null(subset_species)){
    subset_species <- c("", subset_species)
    joined <- joined %>%
      filter(.data$species_latin_fixed %in% subset_species)
  }
  
  if(!is.null(subset_year)){
    subset_year <- c("", subset_year)
    joined <- joined %>%
      filter(.data$year %in% subset_year)
  }
  
  if(!is.null(subset_genus)){
    subset_genus <- c("", subset_genus)
    joined <- joined %>%
      filter(.data$id_genus %in% subset_genus)
  }
  
  #filter on dataset
  
  if(!is.null(dataset)){
    joined <- joined %>%
      filter(.data$project_short_name == dataset)
  }
  
  #filter on trap type (recommended to only take MF)
  if(!is.null(trap_type) & trap_type != "All"){
    joined <- joined %>%
      filter(grepl((trap_type), .data$sample_name))
  }
  
  
  ##Aggregate data to choosen level
  ##Add more choices?
  
  res <- joined
  
  
  ##This is slow because we have to collect the data before we calculate Shannon index.
  ##Best would be to do the Shannon calc on the database side. Seems harder than I first thought.
  if(agg_level == "year_locality"){
    
    res <- res %>%
      collect() %>%
      group_by(year_locality_id, locality_id, species_latin_fixed) %>%
      summarise(no_asv_per_species = n_distinct(.data$sequence_id)) %>%
      group_by(year_locality_id, locality_id) %>%
      summarise(no_species = n_distinct(.data$species_latin_fixed),
                shannon_div = round(calc_shannon(.data$species_latin_fixed), digits),
                mean_asv_per_species = round(mean(.data$no_asv_per_species), digits)) %>%
      left_join(localities,
                by = c("locality_id" = "id"),
                copy = T) %>%
      left_join(year_locality,
                by = c("year_locality_id" = "id",
                       "locality_id" = "locality_id",
                       "ano_flate_id" = "ano_flate_id",
                       "ssbid" = "ssbid"),
                copy = T) %>%
      ungroup() %>%
      select(year,
             locality,
             habitat_type,
             region_name,
             no_species,
             shannon_div,
             mean_asv_per_species) %>%
      arrange(year,
              region_name,
              habitat_type,
              locality)
    
  }
  
  
  if(agg_level == "locality_sampling"){
    
    res <- res %>%
      collect() %>%
      group_by(start_date_obs, end_date_obs, sampling_name, year_locality_id, locality_id, species_latin_fixed) %>%
      summarise(no_asv_per_species = n_distinct(.data$sequence_id)) %>%
      group_by(sampling_name, year_locality_id, locality_id) %>%
      summarise(no_trap_days = mean(as.numeric(.data$end_date_obs - .data$start_date_obs)), ##to get the mean trap days from all traps within the sampling event (should be the same for all traps)
                no_species = n_distinct(.data$species_latin_fixed),
                shannon_div = round(calc_shannon(.data$species_latin_fixed), digits),
                mean_asv_per_species = round(mean(.data$no_asv_per_species), digits)) %>%
      left_join(localities,
                by = c("locality_id" = "id"),
                copy = T) %>%
      left_join(year_locality,
                by = c("year_locality_id" = "id"),
                copy = T) %>%
      ungroup() %>%
      select(year,
             locality,
             sampling_name,
             habitat_type,
             region_name,
             no_trap_days,
             no_species,
             shannon_div,
             mean_asv_per_species) %>%
      arrange(year,
              region_name,
              habitat_type,
              locality,
              sampling_name)
    
  }
  
  if(agg_level == "region_habitat"){
    
    res <- res %>%
      collect() %>%
      group_by(region_name,
               habitat_type,
               species_latin_fixed) %>%
      summarise(no_asv_per_species = n_distinct(.data$sequence_id)) %>%
      group_by(region_name,
               habitat_type) %>%
      summarise(no_species = n_distinct(.data$species_latin_fixed),
                shannon_div = round(calc_shannon(.data$species_latin_fixed), digits),
                mean_asv_per_species = round(mean(.data$no_asv_per_species), digits)) %>%
      ungroup() %>%
      select(habitat_type,
             region_name,
             no_species,
             shannon_div,
             mean_asv_per_species) %>%
      arrange(habitat_type,
              region_name)
    
  }
  
  
  if(agg_level == "region_habitat_year"){
    
    res <- res %>%
      collect() %>%
      group_by(region_name,
               habitat_type,
               year,
               species_latin_fixed) %>%
      summarise(no_asv_per_species = n_distinct(.data$sequence_id)) %>%
      group_by(region_name,
               habitat_type,
               year) %>%
      summarise(no_species = n_distinct(.data$species_latin_fixed),
                shannon_div = round(calc_shannon(.data$species_latin_fixed), digits),
                mean_asv_per_species = round(mean(.data$no_asv_per_species), digits)) %>%
      ungroup() %>%
      select(year,
             habitat_type,
             region_name,
             no_species,
             shannon_div,
             mean_asv_per_species) %>%
      arrange(year,
              habitat_type,
              region_name
      )
    
  }
  
  if(agg_level == "total"){
    
    res <- res %>%
      collect() %>%
      group_by(species_latin_fixed) %>%
      summarise(no_asv_per_species = n_distinct(.data$sequence_id)) %>%
      summarise(no_species = n_distinct(.data$species_latin_fixed),
                shannon_div = round(calc_shannon(.data$species_latin_fixed), digits),
                mean_asv_per_species = round(mean(.data$no_asv_per_species), digits)) %>%
      ungroup() %>%
      select(no_species,
             shannon_div,
             mean_asv_per_species)
    
    
  }
  
  
  
  if(!is.null(limit)){
    res <- joined %>%
      head(limit)
  }
  
  if(as_tibble){
    res <- res %>%
      as_tibble()
  }
  
  
  
  return(res)
  
}



## -----------------------------------------------------------------------------------------------------------------------------------------------
 get_community_matrix <- function(limit = NULL,
                                     id_type = c("metabarcoding"),
                                     trap_type = "MF",
                                     dataset = "NasIns",
                                     subset_years = NULL,
                                     subset_orders = NULL,
                                     subset_families = NULL,
                                     subset_species = NULL,
                                     subset_habitat = NULL,
                                     subset_region = c(NULL, "Østlandet", "Trøndelag"),
                                     exclude_singletons = F,
                                     transposed_matrix = F,
                                     as_tibble = F){

  Norimon:::checkCon()


  dataset <- match.arg(dataset,
                       choices = c("NasIns",
                                   "OkoTrond",
                                   "TidVar",
                                   "Nerlandsøya"))

  trap_type <- match.arg(trap_type,
                         choices = c("MF", "VF", "All", NULL))

  if(!is.null(subset_years)){
  subset_years <- as.numeric(subset_years)
  }


  ##Set up table sources
  ##Probably needs updating after new batch of data. Also need to test filtering of different identification types

  observations <- dplyr::tbl(con, dbplyr::in_schema("occurrences", "observations"))
  identifications <- dplyr::tbl(con, dbplyr::in_schema("events", "identifications"))
  sampling_trap <- dplyr::tbl(con, dbplyr::in_schema("events", "sampling_trap"))
  locality_sampling <- dplyr::tbl(con, dbplyr::in_schema("events", "locality_sampling"))
  year_locality <- dplyr::tbl(con, dbplyr::in_schema("events", "year_locality"))
  localities <- dplyr::tbl(con, dbplyr::in_schema("locations", "localities"))
  identification_techniques <- dplyr::tbl(con, dbplyr::in_schema("lookup", "identification_techniques"))
  traps <- dplyr::tbl(con, dbplyr::in_schema("locations", "traps"))

  ##Join the tables
  joined <- observations %>%
    left_join(identifications,
              by = c("identification_id" = "id"),
              suffix = c("_obs", "_ids")) %>%
    left_join(identification_techniques,
              by = c("identification_name" = "identification_name"),
              suffix = c("_obs", "_idtechn")) %>%
    left_join(sampling_trap,
              by = c("sampling_trap_id" = "id"),
              suffix = c("_obs", "_st")) %>%
    left_join(locality_sampling,
              by = c("locality_sampling_id" = "id"),
              suffix = c("_obs", "_ls")) %>%
    left_join(year_locality,
              by = c("year_locality_id" = "id"),
              suffix = c("_obs", "_yl")) %>%
    left_join(localities,
              by = c("locality_id" = "id"),
              suffix = c("_obs", "_loc")) %>%
    left_join(traps,
              by = c("trap_id" = "id",
                     "year" = "year",
                     "locality" = "locality"))



  if(id_type == "metabarcoding"){
    joined <- joined %>%
      filter(identification_type == "metabarcoding")
  }


  if(!is.null(subset_years)){
    subset_years <- c(NA, subset_years) #To allow one-length subsets
    joined <- joined %>%
      filter(year %IN% subset_years)
  }

  if(!is.null(subset_orders)){
    subset_orders <- c("", subset_orders) #To allow one-length subsets
    joined <- joined %>%
      filter(id_order %IN% subset_orders)
  }

  if(!is.null(subset_families)){
    subset_families <- c("", subset_families)
    joined <- joined %>%
      filter(id_family %IN% subset_families)
  }

  if(!is.null(subset_species)){
    subset_species <- c("", subset_species)
    joined <- joined %>%
      filter(species_latin_fixed %IN% subset_species)
  }

  #Filter on region name
  if(!is.null(subset_region)){
    subset_region <- c("", subset_region)
    joined <- joined %>%
      filter(region_name %IN% subset_region)
  }

  #Filter on habitat type

  if(!is.null(subset_habitat)){
    subset_habitat <- c("", subset_habitat)
    joined <- joined %>%
      filter(habitat_type %IN% subset_habitat)
  }



  #filter on dataset

  if(!is.null(dataset)){
    joined <- joined %>%
      filter(project_short_name == dataset)
  }



  ##Aggregate data to choosen level


  ##Exclude 2020 4 week samplings

  joined <-  joined %>%
    mutate(year = as.character(year)) %>%
    mutate(weeks_sampled = ifelse(grepl("2020", year) & (grepl("1", trap_short_name) | grepl("3", trap_short_name)), 2, 4)) %>%
    mutate(weeks_sampled = ifelse(grepl("2020", year), weeks_sampled, 2))

  joined <- joined %>%
    filter(weeks_sampled == 2)

  #filter on dataset

  #filter on trap type (recommended to only take MF)
  if(!is.null(trap_type) & trap_type != "All"){
    joined <- joined %>%
      filter(grepl((trap_type), sample_name))
  }

  ##Aggregate data to choosen level
  ##Add more choices!

  res <- joined

  res <- res %>%
    collect() %>%
    select(year,
           locality,
           species_latin_fixed) %>%
    group_by(year,
             locality,
             species_latin_fixed) %>%
    summarise(count = n()) %>%
    mutate(present = 1)

  if(exclude_singletons){ ##exlude only species that when observed, was observed more than once each time
    to_exclude <- res %>%
      filter(count == 1) %>%
      select(species_latin_fixed) %>%
      pull()

    res <- res %>%
      filter(!(species_latin_fixed %in% to_exclude))
  }



  res <- res %>%
    select(-count) %>%
    pivot_wider(names_from = species_latin_fixed,
                values_from = present,
                values_fill = 0) %>%
    arrange(year,
            locality)




  if(!is.null(limit)){
    res <- joined %>%
      head(limit)
  }

  if(as_tibble){
    res <- res %>%
      as_tibble()
  }


  if(transposed_matrix){
    res <- res %>%
      select(-c(1:2)) %>%
      t()

  }



  return(res)

}


## -----------------------------------------------------------------------------------------------------------------------------------------------
richn_2020 <- get_observations(subset_year = 2020,
                          agg_level = "total") %>% 
  select(no_species) %>% 
  pull()

richn_2021 <- get_observations(subset_year = 2021,
                          agg_level = "total") %>% 
  select(no_species) %>% 
  pull()

richn_2022 <- get_observations(subset_year = 2022,
                          agg_level = "total") %>% 
  select(no_species) %>% 
  pull()

richn_tot <- get_observations(agg_level = "total") %>% 
  select(no_species) %>% 
  pull()



## ----forskj, fig.cap = "Skisse over en forskjøvet design der hver lokalitet gjenbesøkes etter 5 år.", out_width = "0.7\\textwidth"--------------
dat <- tibble(expand.grid(lok = 1:10, "År" = 2020:2030)) %>% 
  mutate(Lokalitet = factor(paste0("Lokalitet ", lok), levels = rev(c("Lokalitet 1",
                                                                 "Lokalitet 2",
                                                                 "Lokalitet 3",
                                                                 "Lokalitet 4",
                                                                 "Lokalitet 5",
                                                                 "Lokalitet 6",
                                                                 "Lokalitet 7",
                                                                 "Lokalitet 8",
                                                                 "Lokalitet 9",
                                                                 "Lokalitet 10"))),
         År = factor(År),
         Visit = as.integer(0)) %>% 
  arrange(lok, År)

dat$Visit[dat$lok == 1] <- rep(c(1, 0, 0, 0, 0), 2)
dat$Visit[dat$lok == 2] <- rep(c(0, 1, 0, 0, 0), 2)
dat$Visit[dat$lok == 3] <- rep(c(0, 0, 1, 0, 0), 2)
dat$Visit[dat$lok == 4] <- rep(c(0, 0, 0, 1, 0), 2)
dat$Visit[dat$lok == 5] <- rep(c(0, 0, 0, 0, 1), 2)
dat$Visit[dat$lok == 6] <- rep(c(1, 0, 0, 0, 0), 2)
dat$Visit[dat$lok == 7] <- rep(c(0, 1, 0, 0, 0), 2)
dat$Visit[dat$lok == 8] <- rep(c(0, 0, 1, 0, 0), 2)
dat$Visit[dat$lok == 9] <- rep(c(0, 0, 0, 1, 0), 2)
dat$Visit[dat$lok == 10] <- rep(c(0, 0, 0, 0, 1), 2)

dat <- dat %>%
  mutate(Visit = factor(Visit))

ggplot(dat, aes(y = Lokalitet, x = År, fill = Visit)) + 
  geom_tile(colour="white", 
            width=.9, height=.9) + 
  scale_fill_manual(values = c("white", "#008C95" )) +
  theme(legend.position = "none") +
  ylab("")



## -----------------------------------------------------------------------------------------------------------------------------------------------
year_locality <- tbl(con,
                     Id(schema = "events",
                        table = "year_locality")) 

locality_sampling <- tbl(con,
                         Id(schema = "events",
                            table = "locality_sampling"))



## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
localities <- get_localities(dataset = "NasIns") %>% 
  mutate(habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type))




year_locality_biodiv <- get_observations(agg_level = "year_locality",
                                 dataset = "NasIns",
                                 as_tibble = T)

locality_sampling_biodiv <- get_observations(agg_level = "locality_sampling",
                                 dataset = "NasIns",
                                 as_tibble = T)
#locality_sampling_biodiv


biomass_mf_locality_sampling <- get_biomass(agg_level = "locality_sampling",
                                                dataset = "NasIns",
                                                trap_type = "MF")

#biomass_mf_locality_sampling

year_locality_biomass <- get_biomass(dataset = "NasIns",
                                             trap_type = "MF",
                                             agg_level = "year_locality")


## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------

norge <- get_map()



## ----lokaler20202022, fig.cap = "Overvåkingslokaliteter i perioden 2020-2022.", out.width = "\\textwidth"---------------------------------------
ggplot(norge) +
  geom_sf() +
  geom_sf(data = localities,
          aes(fill = habitat_type,
              color = habitat_type),
          lwd = 2) +
  scale_color_nina(name = "Habitat",
                   palette = "orange-green") +
  scale_fill_nina(name = "Habitat",
                   palette = "orange-green")


## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
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
                                query = nas_ins_lok_1000m_agg_q)


## ----calc_ar5_shannon_index_save, eval = FALSE--------------------------------------------------------------------------------------------------
## nas_ins_lok_1000m_agg_brick <- Norimon:::rast_list(nas_ins_lok_1000m_agg)
## save(nas_ins_lok_1000m_agg_brick,
##      file = "out/nas_ins_lok_1000m_agg_brick.Rdata")


## ----calc_ar5_shannon_index, cache = cache_var--------------------------------------------------------------------------------------------------
load(file = "out/nas_ins_lok_1000m_agg_brick.Rdata")
nas_ins_lok_1000m_arealtype_shannon_index <-lsm_l_shdi(nas_ins_lok_1000m_agg_brick)
nas_ins_lok_1000m_arealtype_shannon_index$layer <- names(nas_ins_lok_1000m_agg_brick)

locality_ar5_shdi <- nas_ins_lok_1000m_arealtype_shannon_index %>% 
  select(locality = layer,
         lsm_l_shdi = value) 



## ----calc_ar5_perimeter_length, cache = cache_var-----------------------------------------------------------------------------------------------
locality_ar5_tot_perimeter <- nas_ins_lok_1000m_agg %>% 
  mutate(perimeter = lwgeom::st_perimeter(geom)) %>% 
  st_drop_geometry() %>% 
  group_by(locality) %>% 
  summarize(tot_perimeter_length = sum(perimeter)) %>% 
  select(locality,
         tot_perimeter_length) %>% 
  mutate(tot_perimeter_length = units::drop_units(tot_perimeter_length))

#locality_ar5_tot_perimeter %>% 
#  print(n = Inf)


## ----calc_ar5_area_sum, cache = cache_var-------------------------------------------------------------------------------------------------------
locality_ar5_area <- nas_ins_lok_1000m_agg %>% 
  mutate(area = st_area(geom)) %>% 
  st_drop_geometry() %>% 
  group_by(locality,
           arealtype) %>% 
  summarize(tot_area = sum(area)) %>% 
  pivot_wider(id_cols = locality,
              names_from = arealtype, 
              values_from = tot_area
              ) %>% 
    mutate_if(is.numeric , units::drop_units) %>% 
  mutate_if(is.numeric , replace_na, replace = 0)




## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
loc_types <- localities %>% 
  select(locality,
         habitat_type) 


locality_ar5_area_data <- locality_ar5_area %>% 
  left_join(loc_types,
            by = c("locality" = "locality"),
            copy = T) %>% 
  ungroup() %>% 
  select(-c("locality", 
            "habitat_type",
            "Ikke kartlagt",
            "Hav",
            "geom"))

locality_ar5_area_group <- locality_ar5_area %>% 
  left_join(loc_types,
            by = c("locality" = "locality"),
            copy = T) %>% 
  ungroup() %>% 
  select(c("locality", "habitat_type")) %>% 
  mutate(habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type))




## ---- message= F, results = "hide"--------------------------------------------------------------------------------------------------------------
example_NMDS=metaMDS(locality_ar5_area_data, # Our community-by-species matrix
                     k=2)


## ----ar5nmds, fig.pos = "!h", fig.cap = "NMDS over variasjonen i landskapskomposisjon (arealtyper i AR5, blå tekst i figuren) i en sirkel med 1km radius rundt lokalitetene. Nummer i grønt viser semi-naturlige lokaliteter og nummer i lilla viser lokaliteter i skog.", out.width = "\\textwidth"----
ordiplot(example_NMDS, type="n")

ordihull(example_NMDS, 
         groups = locality_ar5_area_group$habitat_type, 
         draw = "polygon",
         col = "grey90",
         label = T,
         cex = 1)

fixed_loc_names <- locality_ar5_area_group$locality %>% 
  gsub("(.*)([0-9]{2})", "\\2", .)

orditorp(example_NMDS, 
         display = "sites",
         labels = fixed_loc_names,
         col = as.character(factor(locality_ar5_area_group$habitat_type, 
                      labels =  ninaColors()[5:4])),
         air = 0.05,
         cex = 0.7)

orditorp(example_NMDS, 
         display = "species", 
         col = ninaColors()[2], 
         air = 0.01,
         cex = 0.8)


## ----semi-nat-ar5-plot, fig.pos = "!hb", fig.cap = "Eksempel på variasjonen i landskapskomposisjon (arealtype i AR5) i semi-naturlig mark innen en 1 km buffer, SSB 500x500m rute i blått, felleplassering i rødt og de kartlagte ANO-sirklene i lilla (én liten sirkel). Semi-nat_07 har lite fulldyrka jord og mye skog og myr,  Semi-nat_36 har stor andel fulldyrka jord, Semi-nat_20 har mye veier og bebyggelse, mens Semi-nat_36 har mye overflatedyrka mark. \\label{fig:semi-nat-ar5-plot}", fig.height = 8, fig.width = 5----
 
p1 <- locality_ar5_plot("Semi-nat_07",
                        legend = F,
                        strip = T)
p2 <- locality_ar5_plot("Semi-nat_37",
                        legend = F,
                        strip = T) 

p3 <- locality_ar5_plot("Semi-nat_20",
                        legend = F,
                        strip = T) 

p4 <- locality_ar5_plot("Semi-nat_36",
                        legend = F,
                        strip = T) 
 

leg <- locality_ar5_plot("Semi-nat_36",
                         legend = T,
                    strip = T) +
  theme(legend.position = "bottom",
       # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
       legend.key.size = unit(0.5, 'cm'),
       legend.text = element_text(size=7),
       legend.key.width = unit(0.5, 'cm')) +
  guides(fill = guide_legend(nrow = 3, 
                             byrow = TRUE))

legend <- as_grob(get_legend(leg))

blank <- ggplot() + theme_void()


#grid.arrange(grobs = c(glist_1), ncol = 3, as.table = T)

lm <-  matrix(c(1:4, 5, 5),
              ncol = 2,
              byrow = T)


marrangeGrob(list(p1,
             p2,
             p3,
             p4,
             legend), 
             top = "",
             layout_matrix = lm)




## ----skog-ar5-plot, fig.cap = "Eksempel på variasjonen i landskapskomposisjon (arealtype i AR5) i skogsmark innen en 1 km buffer, SSB 500x500m rute i blått, felleplassering i rødt og de kartlagte ANO-sirklene i lilla (små sirkler). Skog_20 inneholder mye myr, Skog_12 er heterogent med  mye dyrket mark av flere typer, Skog_14 er mer homogent med mye fulldyrka mark, og Skog_06 har mye åpen fastmark\\label{fig:skog-ar5-plot}", fig.height = 8, fig.width = 5----

p1 <- locality_ar5_plot("Skog_20",
                        legend = F,
                        strip = T) 

p2 <- locality_ar5_plot("Skog_12",
                        legend = F,
                        strip = T) 
  

p3 <- locality_ar5_plot("Skog_14",
                        legend = F,
                        strip = T)


p4 <- locality_ar5_plot("Skog_06",
                        legend = F,
                        strip = T) 

leg <- locality_ar5_plot("Skog_12",
                         legend = T,
                    strip = T) +
  theme(legend.position = "bottom",
       # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
       legend.key.size = unit(0.5, 'cm'),
       legend.text = element_text(size=7),
       legend.key.width = unit(0.5, 'cm')) +
  guides(fill = guide_legend(nrow = 3, 
                             byrow = TRUE))

legend <- as_grob(get_legend(leg))

blank <- ggplot() + theme_void()


#grid.arrange(grobs = c(glist_1), ncol = 3, as.table = T)

lm <-  matrix(c(1:4, 5, 5),
              ncol = 2,
              byrow = T)


marrangeGrob(list(p1,
             p2,
             p3,
             p4,
             legend), 
             top = "",
             layout_matrix = lm)




## -----------------------------------------------------------------------------------------------------------------------------------------------
get_biomass <- function(limit = NULL,
                            trap_type = "MF",
                            subset_region = c(NULL, "Østlandet", "Trøndelag", "Sørlandet"),
                            dataset = "NasIns",
                            agg_level = "year_locality",
                            as_tibble = F){

  Norimon:::checkCon()



  dataset <- match.arg(dataset,
                       choices = c("NasIns",
                                   "OkoTrond",
                                   "TidVar",
                                   "Nerlandsøya"))

  agg_level <- match.arg(agg_level,
                         choices = c("year_locality",
                                     "locality_sampling",
                                     "total",
                                     "none"))

  trap_type <- match.arg(trap_type,
                         choices = c("All", "MF", "VF", NULL))


  ##Set up table sources
  ##Probably needs updating after new batch of data. Also need to test filtering of different identification types
  sampling_trap <- dplyr::tbl(con, dbplyr::in_schema("events", "sampling_trap"))
  locality_sampling <- dplyr::tbl(con, dbplyr::in_schema("events", "locality_sampling"))
  year_locality <- dplyr::tbl(con, dbplyr::in_schema("events", "year_locality"))
  localities <- dplyr::tbl(con, dbplyr::in_schema("locations", "localities"))
  traps <- dplyr::tbl(con, dbplyr::in_schema("locations", "traps"))


  ##Join the tables

  joined <- sampling_trap %>%
    left_join(locality_sampling,
              by = c("locality_sampling_id" = "id"),
              suffix = c("_obs", "_ls")) %>%
    left_join(year_locality,
              by = c("year_locality_id" = "id"),
              suffix = c("_obs", "_yl")) %>%
    left_join(localities,
              by = c("locality_id" = "id"),
              suffix = c("_obs", "_loc")) %>%
    left_join(traps,
              by = c("trap_id" = "id",
                     "year" = "year")
    ) %>%
    mutate(year = as.character(year))


  ##Exclude 2020 4 week samplings

  joined <-  joined %>%
    mutate(weeks_sampled = ifelse(grepl("2020", year) & (grepl("1", trap_short_name) | grepl("3", trap_short_name)), 2, 4)) %>%
    mutate(weeks_sampled = ifelse(grepl("2020", year), weeks_sampled, 2))

  joined <- joined %>%
    filter(weeks_sampled == 2)

  #filter on dataset

  if(!is.null(dataset)){
    joined <- joined %>%
      filter(project_short_name == dataset)
  }

  #filter on trap type (recommended to only take MF)
  if(!is.null(trap_type) & trap_type != "All"){
    joined <- joined %>%
      filter(grepl((trap_type), sample_name))
  }

  #Filter on region name
  if(!is.null(subset_region)){
    subset_region <- c("", subset_region)
    joined <- joined %>%
      filter(region_name %IN% subset_region)
  }

  ##Aggregate data to choosen level
  ##Add more choices?

  res <- joined



  ##This is slow because we have to collect the data before we calculate Shannon index.
  ##Best would be to do the Shannon calc on the database side. Seems harder than I first thought.
  if(agg_level == "year_locality"){

    res <- res %>%
      collect() %>%
      group_by(year_locality_id, locality_id) %>%
      summarise(sum_wet_weight = sum(wet_weight_bottle - weight_empty_bottle, na.rm = T),
                avg_wet_weight = mean(wet_weight_bottle - weight_empty_bottle, na.rm = T)) %>%
      left_join(localities,
                by = c("locality_id" = "id"),
                copy = T) %>%
      left_join(year_locality,
                by = c("year_locality_id" = "id",
                       "locality_id" = "locality_id",
                       "ano_flate_id" = "ano_flate_id",
                       "ssbid" = "ssbid"),
                copy = T) %>%
      ungroup() %>%
      select(year,
             locality,
             habitat_type,
             region_name,
             sum_wet_weight,
             avg_wet_weight) %>%
      arrange(year,
              region_name,
              habitat_type,
              locality)

  }


  if(agg_level == "locality_sampling"){

    res <- res %>%
      collect() %>%
      group_by(sampling_name, year_locality_id, locality_id) %>%
      summarise(no_trap_days = mean(as.numeric(end_date_obs - start_date_obs)),
                sum_wet_weight = sum(wet_weight_bottle - weight_empty_bottle, na.rm = T),
                avg_wet_weight = mean(wet_weight_bottle - weight_empty_bottle, na.rm = T)) %>%
      left_join(localities,
                by = c("locality_id" = "id"),
                copy = T) %>%
      left_join(year_locality,
                by = c("year_locality_id" = "id"),
                copy = T) %>%
      ungroup() %>%
      select(year,
             locality,
             sampling_name,
             habitat_type,
             region_name,
             no_trap_days,
             sum_wet_weight,
             avg_wet_weight) %>%
      arrange(year,
              region_name,
              habitat_type,
              locality,
              sampling_name)

  }



  if(agg_level == "total"){

    res <- res %>%
      collect() %>%
      summarise(sum_wet_weight = sum(wet_weight_bottle - weight_empty_bottle, na.rm = T),
                avg_wet_weight = mean(wet_weight_bottle - weight_empty_bottle, na.rm = T)) %>%
      select(sum_wet_weight,
             avg_wet_weight)

  }

  if(!is.null(limit)){
    res <- joined %>%
      head(limit)
  }


  res <- res %>%
    mutate(sum_wet_weight = round(sum_wet_weight, 2),
           avg_wet_weight = round(avg_wet_weight, 2))

  if(as_tibble){
    res <- res %>%
      as_tibble()
  }



  return(res)

}


## -----------------------------------------------------------------------------------------------------------------------------------------------
tot_biomass <- dbGetQuery(con,
                          "
                          SELECT sum(st.wet_weight)
                          FROM events.sampling_trap st,
                          events.locality_sampling ls,
                          events.year_locality yl
                          WHERE st.locality_sampling_id = ls.id
                          AND ls.year_locality_id = yl.id
                          AND yl.project_short_name = 'NasIns'
                          
                          ")


biomass_per_loc_year_raw <- get_biomass(trap_type = "All",
                                         agg_level = "year_locality",
                                         dataset = "NasIns"
                                          #  ,subset_region = "Østlandet"
                                        )

# tot_biomass <- get_biomass(trap_type = "All",
#                                          agg_level = "total",
#                                          dataset = "NasIns"
#                                           #  ,subset_region = "Østlandet"
#                                         ) %>% 
#   select(sum_wet_weight)

biomass_per_loc_year_mean <- biomass_per_loc_year_raw %>% 
  summarise(mean_sum_wet_weight = mean(sum_wet_weight),
            sd_sum_wet_weight = sd(sum_wet_weight))


biomass_per_year_q <-"

  SELECT year, min(ls.start_date), max(ls.end_date), round((sum(st.wet_weight) / 1000)::numeric, 2) as kg
                          FROM events.sampling_trap st,
                          events.locality_sampling ls,
                          events.year_locality yl
                          WHERE st.locality_sampling_id = ls.id
                          AND ls.year_locality_id = yl.id
                          AND yl.project_short_name = 'NasIns'
                          GROUP BY year
                          

" 

biomass_per_year <- dbGetQuery(con,
                               biomass_per_year_q) %>% 
  select(`År` = year,
        `Start` = min,
        `Slutt` = max,
        `Våtvekt (kg)` = kg) %>% 
  janitor::adorn_totals("row",
                        name = "Sum")



## ----biomass_year, results = "asis"-------------------------------------------------------------------------------------------------------------
print(xtable(biomass_per_year, caption = "Total biomasse insekter innsamlet i prosjektet så langt. \\label{tab:biomass_year}"))



## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------

no_samplings_per_lok_q <- "
SELECT l.locality ,  
count(distinct ls.sampling_name)::integer,
min(ls.start_date) as  start_date,
max(ls.end_date) as end_date 
FROM events.locality_sampling ls,
events.year_locality yl,
locations.localities l,
events.sampling_trap st
WHERE ls.year_locality_id = yl.id
AND yl.locality_id = l.id
AND yl.year = 2022
AND yl.project_short_name = 'NasIns'
AND st.locality_sampling_id = ls.id
AND st.received_at_lab IS TRUE
GROUP BY locality
ORDER BY locality
" 
no_samplings_per_lok <- dbGetQuery(con,
                                   no_samplings_per_lok_q) %>% 
  as_tibble %>% 
  mutate(start_date  = as.character(start_date),
         end_date  = as.character(end_date)) %>% 
  select(Lokalitet = locality,
         `Antall felletømminger` = count,
         `Start prøvetaking` = start_date,
         `Slutt prøvetaking` = end_date)



## ----no-samplings-per-lok, results = "asis"-----------------------------------------------------------------------------------------------------
print(xtable(no_samplings_per_lok, caption = "Antall felletømminger per lokalitet i 2022. \\label{tab:no-samplings-per-lok}"))


## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
no_samples_q <-  "

SELECT trap_model as \"Felletype\", 
count(distinct st.sample_name)::integer as \"Antall prosesserte felleprøver\"
FROM events.locality_sampling ls,
events.year_locality yl,
locations.traps,
events.sampling_trap st
WHERE ls.year_locality_id = yl.id
AND st.trap_id = traps.id
AND yl.year = 2022
AND yl.project_short_name = 'NasIns'
AND st.locality_sampling_id = ls.id
AND st.received_at_lab IS TRUE
GROUP BY trap_model
"

no_samples <- dbGetQuery(con,
                           no_samples_q) %>% 
  mutate(Felletype = ifelse(Felletype == "Malaise_Watkins_Black", "Malaisefelle", "Vindusfelle"))


## ---- antall-prover, results = "asis"-----------------------------------------------------------------------------------------------------------
print(xtable(no_samples, caption = "Antall prosesserte prøver i 2022. De 4 vindusfelleprøvene i hver tømmerunde på en plass ble slått sammmen og metastrekkodet sammen. \\label{tab:antall-prover}"))


## ----refbase, results = "asis"------------------------------------------------------------------------------------------------------------------

refbase_tab <- tibble(Klasse = c("Arachnida (edderkoppdyr)",
                                "Chilopoda (skolependere)",
                                "Collembola (spretthaler)",
                                "Diplopoda (tusenbein)",
                                "Insecta (insekter)",
                                "Ikke målgruppe"),
                      Porter = c(4778,
                                 135,
                                 268,
                                 246,
                                 65660,
                                 39693),
                      `Modifisert 2021` = c(4952,
                                         147,
                                         377,
                                         246,
                                         69426,
                                         39695),
                      `Modifisert 2022` = c(5124,
                                         147,
                                         388,
                                         246,
                                         72151,
                                         39698))


xtable(refbase_tab, digits = 0,
       caption = "Antall arter i referansebasene som er blitt brukt i prosjektet så langt. I 2020 tok vi utgangspunkt i en base fra Porter, men i 2021 og 2022 har vi utvided denne med en betydende mengde norske arter.\\label{tab:refbase}")



## ---- fig.cap = "Flytskjema for identifisering av målgruppe sekvensene og kvalitetssikring av vår referansedatabase for kontaminasjon og referansesekvensmangel. \\label{fig:dna-flow-1}", out.width = "\\textwidth"----
include_graphics("figures/flyt_1.pdf")

## ---- fig.cap = "Flytskjema for klassifisering og kvalitetssikring av arter basert på DNA-metastrekkoding. Kun arter med høy arts-konfidens er tatt med i den videre arbeidsflyten som avhenger av artsidentitet. \\label{fig:dna-flow-2}", out.width = "\\textwidth"----
include_graphics("figures/flyt_2.pdf")


## ---- fig.cap = "Flytskjema for horizon scanning for fremmede arter og nye arter for Norge. *I tillegg til å kryssjekke arter mot den norske fremmedartslista har vi også kryssjekket artsfunn i dette prosjektet mot databasene Daisie, Saul Pathways, Roy og List of Union Concern. Arter som er listet i en av disse databasene, er antatt å være fremmede arter i et eller flere land i Europa. Vi fant ingen arter lister i Roy et al. eller List of Union Concern i dette prosjektet. \\label{fig:dna-flow-2}", out.width = "0.9\\textwidth", fig.pos = "H"----
include_graphics("figures/flyt_3.pdf")


## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
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


rank_abund_raw <- dbGetQuery(con,
                         rank_abund_raw_q) %>% 
  as_tibble()


## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
rank_abund <- rank_abund_raw %>% 
  arrange(row_number, rank) %>% 
  select(row_number,
         abund,
         taxonomic_level,
         species_latin) %>% 
  pivot_wider(names_from = taxonomic_level,
              values_from = abund,
              values_fill = 0) %>% 
  mutate(count_species = cumsum(species > 0),
         count_genus = cumsum(genus > 0),
         count_family = cumsum(family > 0),
         count_order = cumsum(order > 0),
         count_class = cumsum(class > 0),
         all_abund = species + family + order + genus + class) %>% 
  mutate(cumsum_abund = cumsum(as.numeric(all_abund)))%>% 
  mutate(pct_all_sequences = cumsum_abund / sum(all_abund)) 


rank_abund_long <- rank_abund %>% 
  select(rank_abundance = row_number,
         species_latin,
         Art = count_species,
         Slekt = count_genus,
         Familje = count_family,
         Orden = count_order,
         Klasse = count_class) %>% 
  pivot_longer(cols = c("Art", 
                        "Slekt",
                        "Familje",
                        "Orden",
                        "Klasse"),
               names_to = "taxonomic_level") %>% 
  mutate(taxonomic_level = factor(taxonomic_level, levels = c("Klasse", "Orden", "Familje", "Slekt", "Art"))) %>% 
  arrange(taxonomic_level,
          value)




## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------

rank_abund_order <- rank_abund_raw %>% 
  arrange(row_number, rank) %>%
  group_by(order) %>% 
  mutate(order_rank = rank(rank, ties.method = "first")) %>% 
  mutate(to_species = as.integer(taxonomic_level == "species")) %>%
   filter(order %in% c("Diptera", "Hymenoptera", "Coleoptera", "Lepidoptera", "Hemiptera")) %>% 
  select(rank,
         order,
         species_latin,
         order_rank,
         to_species)




## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
order_spec_mod <- glm(to_species ~ order_rank * order, family = binomial() , data = rank_abund_order)
#summary(order_spec_mod)


## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
new_data_limits <- rank_abund_order %>% 
  group_by(order) %>% 
  summarise(no_otu = max(order_rank))


## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------

#Ugly, could make an apply version

new_data <- tibble(order = "Diptera",
                   order_rank = seq(1, new_data_limits$no_otu[new_data_limits$order == "Diptera"], by = 50)
) %>% 
  rbind(tibble(order = "Hymenoptera",
                   order_rank = seq(1, new_data_limits$no_otu[new_data_limits$order == "Hymenoptera"], by = 50)
)) %>% 
  rbind(tibble(order = "Coleoptera",
                   order_rank = seq(1, new_data_limits$no_otu[new_data_limits$order == "Coleoptera"], by = 50)
)) %>% 
  rbind(tibble(order = "Lepidoptera",
                   order_rank = seq(1, new_data_limits$no_otu[new_data_limits$order == "Lepidoptera"], by = 50)
)) %>% 
  rbind(tibble(order = "Hemiptera",
                   order_rank = seq(1, new_data_limits$no_otu[new_data_limits$order == "Hemiptera"], by = 50)
))
                  


pred_order_id <- predict(order_spec_mod, newdata = new_data, type = "response")
new_data$pred <- pred_order_id

new_data <- new_data %>% 
   mutate(order = factor(order, levels = c("Lepidoptera", "Hymenoptera", "Diptera", "Coleoptera", "Hemiptera")))



new_data <- new_data %>% 
  group_by(order) %>% 
  mutate(rel_order_rank =  rank(rev(order_rank)) / max(rank(rev(order_rank)))) 
  


## ----artsniva, fig.pos = '!b', fig.cap = "Modellering av sannsynligheten for klassifisering til artsnivå (y-aksen) i forhold til hvor mange DNA-kopier vi finner fra arten (x-aksen). Til venstre på figurene er arter med veldig mange DNA-kopier (vanlige arter), mens til høyre er arter med få DNA-kopier (sjeldne arter). A) Totalt antall arter påvist i 2022, der x-aksen viser artsantall og B) fem utvalgte ordener som viser forskjellen mellom ulike grupper av insekter. \\label{fig:artsniva}"----

Afig <- rank_abund_long %>% 
  #filter(rank_abundance <= 13200) %>% 
ggplot(., aes(x = rank_abundance, y = value, fill = taxonomic_level)) + 
  geom_bar(stat = "identity", position = "fill", width = 100) + 
  ylab("Andel arter identifisert") + 
  xlab("Arter i rekkefølge av mengde DNA") + 
  scale_fill_nina() +
   theme(panel.background = element_blank(),
        axis.line=element_line(colour="black"),
        legend.title = element_blank(),
        legend.position="bottom") +
  guides(fill = guide_legend(override.aes = list(size=5), nrow = 2, byrow = T)) +
  geom_vline(xintercept = min(rank_abund$row_number[rank_abund$pct_all_sequences > 0.899999999]), linetype = "dashed")



# logit_perc <- scales::trans_new("logit perc",
# transform = function(x)qlogis(x/100),
# inverse = function(x)100*plogis(x)
# )

Bfig <- ggplot(new_data) +
  geom_line(aes(x = rel_order_rank, y = pred, color = order),
            size = 2.5) +
  theme(
    #axis.text.x = element_blank(),
        axis.line=element_line(colour="black"),
        panel.background=element_blank(),
       # axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.position="bottom",
       legend.text=element_text(size = 6)) + 
   guides(color = guide_legend(override.aes = list(size =5), nrow = 2, byrow = T)) + 
  ylab("Sannsynlighet for identifikasjon til artsnivå") +
  scale_color_nina() +
  xlab("Arter i rekkefølge av mengde DNA") +
  scale_x_reverse(breaks = c(1, 0.75, 0.5, 0.25, 0),
                  labels = c("Mest\nDNA", "3/4", "1/2", "1/4", "Minst\nDNA"))


#+
 # xlim(c(0, 6000))
# +  scale_y_continuous(trans = "logit")


grid.arrange(Afig, Bfig, ncol = 2)


## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
spec_seq_depth_q <- "
SELECT l.habitat_type, st.sample_name, count(distinct o.species_latin)::numeric as species_richn, i.read_abundance::numeric
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
AND yl.year = 2022
--OR yl.year = 2021) --could include 2021 as well
GROUP BY st.sample_name, i.read_abundance, l.habitat_type
"

spec_seq_depth <- dbGetQuery(con,
                             spec_seq_depth_q) %>% 
  as_tibble()

#Trond_2022_Semi-nat_53_sampling-01_MF1 looks to be an outlier. Typing error?

spec_seq_depth <- spec_seq_depth %>% 
  filter(sample_name != "Trond_2022_Semi-nat_53_sampling-01_MF1")



## ---- eval = F----------------------------------------------------------------------------------------------------------------------------------
## seq_depth_mod <- glm(species_richn ~ read_abundance * habitat_type, data = spec_seq_depth)
## 
## summary(seq_depth_mod)
## 
## plot(allEffects(seq_depth_mod)) #Tja,


## ----sekvenseringsdybde, fig.pos = '!h', fig.cap = "Sammenhengen mellom antall DNA-kopier (sekvenseringsdybde) per prøve på x-aksen og antall arter påvist i den samme prøven på y-aksen.\\label{fig:sekveringsdybde}", out.width = "0.6\\textwidth"----
ggplot(spec_seq_depth, aes(x = read_abundance, y = species_richn, col = habitat_type)) + 
  geom_point() + 
  theme(panel.background = element_blank(), axis.line = element_line(color="black")) + 
  xlab("Sekvenseringsdybde") + 
  ylab("Antall arter") + 
  scale_color_manual(values = c("#004F71", "#008C95"), name = "Habitat")


## -----------------------------------------------------------------------------------------------------------------------------------------------
physeq.species=readRDS("figures/ord_physeq.rds")



## -----------------------------------------------------------------------------------------------------------------------------------------------
col_arne <- openxlsx::read.xlsx("../../rawData/validering/Coleoptera_man_ident_2022_Arne.xlsx") %>% 
  as_tibble

col_man <- col_arne %>% 
  mutate(species_latin = paste0(str_trim(genus), " ", str_trim(species)),
         sample_name = gsub("(.*)(-4$)", "\\1", sample_name)) %>% 
  select(sample_name,
         species_latin,
         species) %>% 
  distinct()

#col_man
col_val_samples <- col_man %>% 
  select(sample_name) %>% 
  distinct() %>% 
  pull()
  
col_val_samples <- paste0(c("('", 
                         paste0(col_val_samples,
                                collapse = "','"),
                         "')"),
                         collapse = "")


## -----------------------------------------------------------------------------------------------------------------------------------------------
col_dna_q <- paste0("
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


col_dna <- dbGetQuery(con,
                      col_dna_q)


## -----------------------------------------------------------------------------------------------------------------------------------------------
col_merged <- col_man %>% 
  full_join(col_dna,
            by = c("sample_name" = "sample_name",
                   "species_latin" = "species_latin"),
            keep = FALSE,
            suffix = c("_man", "_dna"))

col_merged <- col_merged %>% 
  mutate(man_to_species = !grepl("(.p\\.)", species_latin) & !is.na(species),
         man_to_genus = grepl("(.p\\.)", species_latin)  & !is.na(species),
         dna_to_species = taxonomic_level == "species",
         dna_to_genus = taxonomic_level == "genus")


## -----------------------------------------------------------------------------------------------------------------------------------------------
#Got col_merged back after manual revision

col_merged <- openxlsx::read.xlsx("../../rawData/validering/col_merged.xlsx",
                                   startRow = 2) %>% 
  as_tibble()
#col_merged



## -----------------------------------------------------------------------------------------------------------------------------------------------
col_agg <- col_merged %>% 
  group_by(sample_name) %>% 
  summarise(tot_spec = n(),
            no_man_to_species = sum(man_to_species, na.rm = TRUE),
            no_man_to_genus = sum(man_to_genus, na.rm = TRUE),
            no_dna_to_species = sum(dna_to_species, na.rm = TRUE),
            no_dna_to_genus = sum(dna_to_genus, na.rm = TRUE),
            no_tax_uncertain = sum(tax_unc, na.rm = TRUE)) %>% 
  mutate(no_dna_missed = tot_spec - (no_dna_to_species + no_dna_to_genus + no_tax_uncertain),
         no_man_missed = tot_spec - (no_man_to_species + no_man_to_genus)
  )
  
col_longer <- col_agg %>% 
  select(-tot_spec) %>% 
  pivot_longer(cols = !"sample_name") %>% 
  mutate(method = ifelse(grepl("man", name), "Morfologi", "Metastrekkoding"),
         kategori = ifelse(grepl("to_species", name), "Bestemt til art", NA)) %>% 
  mutate(kategori = ifelse(grepl("to_genus", name), "Bestemt til slekt", kategori)) %>% 
  mutate(kategori = ifelse(grepl("missed", name), "Ikke oppdaget", kategori)) %>% 
  mutate(kategori = ifelse(grepl("uncertain", name), "Usikker taksonomi", kategori)) %>% 
  mutate(kategori = factor(kategori, levels = c("Ikke oppdaget", "Usikker taksonomi", "Bestemt til slekt", "Bestemt til art"))) %>% 
  mutate(locality = gsub("(Ost_2022_)(.*)(_sampling.*)", "\\2", sample_name))



## -----------------------------------------------------------------------------------------------------------------------------------------------
col_val_plot <- ggplot(col_longer,
               aes(x=method,
                   y=value,
                   fill=kategori)) +
    geom_bar(stat="identity") +
    facet_wrap(~ locality,
               ncol = 5)+
    #scale_fill_nina() +
    scale_fill_manual(values=c("grey85","#E57200", "#008C95", "#004F71")) +
    theme(legend.title = element_blank(),
          axis.text.x = element_text(angle=90),
          panel.background = element_blank(),
          axis.line = element_line(colour="black"),
          axis.title.x=element_blank(),
          strip.text=element_text(size=5.5))+
    scale_x_discrete(labels=c("DNA","Morf"))+
    ylab("Antall arter")+
    ggtitle("Coleoptera")




## -----------------------------------------------------------------------------------------------------------------------------------------------


col_man_matrix <- col_merged %>% 
  filter(!is.na(species)) %>% 
  mutate(locality = gsub("(Ost_2022_)(.*)(_sampling.*)", "\\2", sample_name),
         present = 1) %>% 
  select(locality,
         species_latin,
         present) %>% 
   pivot_wider(names_from = species_latin, 
              values_from = present,
              values_fill = 0
              )
  
col_dna_matrix <- col_merged %>% 
  filter(!is.na(taxonomic_level)) %>% 
  mutate(locality = gsub("(Ost_2022_)(.*)(_sampling.*)", "\\2", sample_name),
         present = 1) %>% 
  select(locality,
         species_latin,
         present) %>% 
   pivot_wider(names_from = species_latin, 
              values_from = present,
              values_fill = 0
              ) 
  
col_to_iNEXT <- list(t(col_man_matrix[-1]),
                     t(col_dna_matrix[-1]))
  
names(col_to_iNEXT) <- c("Morf", "DNA")



## -----------------------------------------------------------------------------------------------------------------------------------------------
col_inext <- iNEXT(col_to_iNEXT,
                   datatype = "incidence_raw",
                   endpoint = 50)

col_fort <- fortify(col_inext)
 


## -----------------------------------------------------------------------------------------------------------------------------------------------
syrph_sondre <- openxlsx::read.xlsx("../../rawData/validering/Syrphidae_man_ident_2022.xlsx",
                                  sheet = 2) %>% 
  as_tibble

syrph_man <- syrph_sondre %>% 
  mutate(species_latin = paste0(str_trim(genus), " ", str_trim(species)),
         sample_name = gsub("(.*)(-4$)", "\\1", sample_name)) %>% 
  select(sample_name,
         genus,
         species) %>% 
  filter(!is.na(genus)) %>% 
  mutate(species_latin = paste0(genus, " ", species)) %>% 
  distinct()

#syrph_man
syrph_val_samples <- syrph_man %>% 
  select(sample_name) %>% 
  distinct() %>% 
  pull()
  
syrph_val_samples <- paste0(c("('", 
                         paste0(syrph_val_samples,
                                collapse = "','"),
                         "')"),
                         collapse = "")


## -----------------------------------------------------------------------------------------------------------------------------------------------
#The Sør samples had single digits as sample events, but the rest had 2 digit values. (I updated the trigger after the Sør locations was initiated). Doing this manually.
syrph_dna_q <- paste0("
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
AND st.sample_name in ('Sor_2022_Semi-nat_33_sampling-5_MF1',
					   'Sor_2022_Semi-nat_36_sampling-5_MF1'
					   'Sor_2022_Semi-nat_39_sampling-5_MF1',	
					   'Ost_2022_Semi-nat_41_sampling-03_MF1',	
					   'Ost_2022_Semi-nat_43_sampling-03_MF1',
					   'Ost_2022_Semi-nat_44_sampling-03_MF1',
					   'Ost_2022_Semi-nat_46_sampling-03_MF1',
					   'Trond_2022_Semi-nat_51_sampling-04_MF1',
					   'Trond_2022_Semi-nat_55_sampling-04_MF1',
					   'Trond_2022_Semi-nat_56_sampling-04_MF1'
)

")


syrph_dna <- dbGetQuery(con,
                      syrph_dna_q)


## -----------------------------------------------------------------------------------------------------------------------------------------------
syrph_merged <- syrph_man %>% 
  full_join(syrph_dna,
            by = c("sample_name" = "sample_name",
                   "species_latin" = "species_latin"),
            keep = FALSE,
            suffix = c("_man", "_dna"))

syrph_merged <- syrph_merged %>% 
  mutate(man_to_family = grepl("Ukjent", genus) ,
         man_to_species = !is.na(species),
         man_to_genus = !grepl("Ukjent", genus) & grepl("NA", species_latin),
         dna_to_species = taxonomic_level == "species",
         dna_to_genus = taxonomic_level == "genus")


## -----------------------------------------------------------------------------------------------------------------------------------------------
#Got syrph_merged back after manual check
syrph_merged <- openxlsx::read.xlsx("../../rawData/validering/syrph_merged.xlsx",
                                    startRow = 2) %>% 
  as_tibble()

#syrph_merged



## -----------------------------------------------------------------------------------------------------------------------------------------------
syrph_agg <- syrph_merged %>% 
  group_by(sample_name) %>% 
  summarise(tot_spec = n(),
            #no_man_to_family = sum(man_to_family, na.rm = TRUE),
            no_man_to_species = sum(man_to_species, na.rm = TRUE),
            no_man_to_genus = sum(man_to_genus, na.rm = TRUE),
            no_dna_to_species = sum(dna_to_species, na.rm = TRUE),
            no_dna_to_genus = sum(dna_to_genus, na.rm = TRUE),
            no_dna_ambig = sum(dna_ambig, na.rm = TRUE)) %>% 
  mutate(no_dna_missed = tot_spec - (no_dna_to_species + no_dna_to_genus + no_dna_ambig),
         no_man_missed = tot_spec - (no_man_to_species + no_man_to_genus))
  
syrph_longer <- syrph_agg %>% 
  select(-tot_spec) %>% 
  pivot_longer(cols = !"sample_name") %>% 
  mutate(method = ifelse(grepl("man", name), "Morfologi", "Metastrekkoding"),
         kategori = ifelse(grepl("to_species", name), "Bestemt til art", NA)) %>% 
  mutate(kategori = ifelse(grepl("to_genus", name), "Bestemt til slekt", kategori),
         kategori = ifelse(grepl("to_family", name), "Bestemt til familie", kategori),
         kategori = ifelse(grepl("no_dna_ambig", name), "DNA kan ikke skille\nmellom arter", kategori)) %>%
  mutate(kategori = ifelse(grepl("missed", name), "Ikke oppdaget", kategori)) %>% 
  mutate(kategori = factor(kategori, levels = c("Ikke oppdaget","DNA kan ikke skille\nmellom arter",  "Bestemt til familie", "Bestemt til slekt", "Bestemt til art"))) %>% 
  mutate(locality = gsub("(.*)(_2022_)(.*)(_sampling.*)", "\\3", sample_name))



## -----------------------------------------------------------------------------------------------------------------------------------------------
syrph_val_plot <- ggplot(syrph_longer,
               aes(x = method,
                   y = value,
                   fill = kategori)) +
    geom_bar(stat="identity") +
    facet_wrap(~ locality,
               ncol = 10)+
    #scale_fill_nina() +
    scale_fill_manual(values=c("grey85","#E57200", "#008C95", "#004F71")) +
    theme(legend.title = element_blank(),
          axis.text.x = element_text(angle = 90),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.title.x=element_blank(),
          strip.text = element_text(size = 5.5,
                                    angle = 90)) +
    scale_x_discrete(labels=c("DNA","Morf")) +
    ylab("Antall arter") +
    ggtitle("Syrphidae")




## -----------------------------------------------------------------------------------------------------------------------------------------------

syrph_man_matrix <- syrph_merged %>% 
  filter(!is.na(species)) %>% 
  mutate(locality =  gsub("(.*)(_2022_)(.*)(_sampling.*)", "\\3", sample_name),
         present = as.integer(1)) %>% 
  select(locality,
         species_latin,
         present) %>% 
  #distinct() %>% 
   pivot_wider(names_from = species_latin, 
              values_from = present,
              values_fill = 0
              )
  
syrph_dna_matrix <- syrph_merged %>% 
  filter(!is.na(taxonomic_level)) %>% 
  mutate(locality =  gsub("(.*)(_2022_)(.*)(_sampling.*)", "\\3", sample_name),
         present = 1) %>% 
  select(locality,
         species_latin,
         present) %>% 
   pivot_wider(names_from = species_latin, 
              values_from = present,
              values_fill = 0
              ) 
  
syrph_to_iNEXT <- list(t(syrph_man_matrix[-1]),
                     t(syrph_dna_matrix[-1]))
  
names(syrph_to_iNEXT) <- c("Morf", "DNA")



## -----------------------------------------------------------------------------------------------------------------------------------------------
syrph_inext <- iNEXT(syrph_to_iNEXT,
                   datatype = "incidence_raw",
                   endpoint = 50)

syrph_fort <- fortify(syrph_inext)
 


## -----------------------------------------------------------------------------------------------------------------------------------------------
col_fort <- col_fort %>% 
  mutate(color = 'Col_C')

syrph_fort <- syrph_fort %>% 
  mutate(color = 'Col_S')

df <- col_fort %>% 
  rbind(syrph_fort)


df.point <- df[which(df$Method == "Observed"),]
df.point$shape[grepl("Morf",df.point$Assemblage)] = "Morfologi"
df.point$shape[grepl("DNA",df.point$Assemblage)]="DNA"
#df.point$color=gsub("Morf|DNA","",df.point$site)
#df.point$color=gsub("DNA","",df.point$color)

df.line <- df[which(df$Method!="Observed"),]
df.line$method <- factor(df.line$Method, 
                         c("Rarefaction", "Extrapolation"),
                         c("Rarefaction", "Extrapolation"))
#df.line$color=gsub("Morph|DNA","",df.line$site)
df.line$shape[grepl("Morf",df.line$Assemblage)] = "Morfologi"
df.line$shape[grepl("DNA",df.line$Assemblage)] = "DNA"

df$shape[grepl("Morph",df$Assemblage)] = "Morfologi"
df$shape[grepl("DNA",df$Assemblage)] = "DNA"
#df$color=gsub("Morph|DNA","",df$site)



## ----col-lepi-curves, fig.pos = "!hb", fig.cap = "Artsakkumulasjonskurver for biller og sommerfugler fra malaise og vindusfeller ved bruk av morfologiske og DNA bestemmelsemetoder. Noter at hver av vindusfellene egentlig er 4 feller slått sammen, slik at de tilsammen er 20 stk. \\label{fig:col-lepi-curves}", fig.height = 6, out.width = "0.6\\textwidth"----
lines <- c("Rarefaction" = "solid", 
           "Extrapolation" = "dotted")
ggplot(df, aes(x = x, y = y, colour = color, shape = shape, fill = color)) + 
  geom_point(aes(shape = shape, colour = color), size = 5, data = df.point) +
  geom_line(aes(linetype = method), lwd = 1.5, data = df.line) +
  geom_ribbon(aes(ymin = y.lwr, ymax = y.upr,
                  fill = color, colour = NULL), alpha = 0.2, data = df) +
  labs(x = "Antall feller", y = "Antall arter") +
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        text = element_text(size = 10)
        ,
       legend.box = "vertical"
       ) +
  scale_fill_manual(values=c("#004F71","#E57200"),
                    labels=c("Coleoptera vindusf.","Syrphidae malaisef.")) +
  scale_color_manual(values=c("#004F71","#E57200"),
                     labels=c("Coleoptera vindusf.","Syrphidae malaisef.")) +
  xlim(c(0, 20)) +
     scale_linetype_manual(values = lines,
                       name = "Method", 
                       labels = c("Interpolert", "Ekstrapolert"))



## ----dna-comp,fig.pos = "!hb", fig.cap = "Sammenligning av antall arter påvist for DNA-basert og morfologisk artsbestemmelse for biller (Coleoptera) fanget med 4 vindusfeller og blomsterfluer (Syrphidae) fanget med én malaisefelle. De grå søylene representerer det totale artsantallet som begge metodikker tilsammen oppdaget i hver prøve, og de fargete søylene representerer de arter som respektive metodikk oppdaget. For billene er det én art som har usikker taksonomi. Denne er markert kun i DNA-søylen. For blomsterfluer kan DNA-markøren ikke brukes for å skille mellom flere arter. Disse artsbestemminger kan iblant revideres basert på kjente forekomster. Merk at de manuelle artsbestemmelsene heller ikke er 100 % sikre, og at særlig noen blomfluger har vært vanskelige å skille mellom. \\label{fig:dna-comp}", fig.height = 6----
grid.arrange(col_val_plot,
             syrph_val_plot,
             ncol = 1)


## ----syrph-sort, fig.pos = "!b", fig.show="hold", fig.cap = "Blomsterfluer sortert ut for validering etter lysering, fra tredje samplingen i lokalitet Semi\\textemdash nat 43.\\label{fig:syrph-sort}"----
include_graphics("figures/syrphidae_sortering.jpg")


## -----------------------------------------------------------------------------------------------------------------------------------------------
ano_herb_agg <- tbl(con,
                    Id(schema = "views", 
                       table = "ano_herb_agg")) %>% 
  filter(project_short_name == "NasIns")



## -----------------------------------------------------------------------------------------------------------------------------------------------
loc_reg <- localities %>% 
  st_drop_geometry() %>% 
  select(locality,
         region_name)




## ----ano-herb-spec-nor, eval = T, echo = F, fig.pos = "!b", fig.cap = "Dekningsgrad (y-akse) og artsantall av karplanter (størrelse) fra ANO-kartlegging i overvåkingslokalitetene. For skogshabitat med full ANO-kartlegging viser vi medianverdier fra de 18 rutene som ble målt. Dekningsgrader på over 100 \\% er mulig når forskjellige arter vokser over og under hverandre med overlapp, da man registrerer dekningsgrad hvor hver enkelt art. Punktene er noe spredt ut for å kunne skille mellom region (x-akse) og overlappende verdier (y-akse).", out.width = "0.9\\textwidth"----
#Legg til, oppdatere Norimon


ano_herb_agg %>% 
  left_join(loc_reg,
            by = c("locality" = "locality"),
            copy = TRUE) %>% 
  group_by(region_name) %>% 
  ggplot(.) +
  geom_point(aes(x = year,
                 y = ano_median_cover,
                 size = ano_median_no_spec,
                 col = region_name),
             alpha = 0.7,
             position = position_jitterdodge(dodge.width = 0.4, 
                                             jitter.height = 2,
                                             jitter.width = 0)) +
  scale_x_continuous(breaks= seq(2020, 2022)) +
  scale_color_nina(name = "Region") +
  scale_size(name = "Artantall\nkarplanter") +
   guides(color = guide_legend(override.aes = list(size=5))) +
  xlab("År") +
  ylab("Dekningsgrad %")





## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
year_locality_biodiv_ano <- year_locality_biodiv %>%
  left_join(ano_herb_agg,
            by = c("locality" = "locality"),
            copy = T)
  


## -----------------------------------------------------------------------------------------------------------------------------------------------
tree_agg <- tbl(con,
                Id(schema = "views",
                   table = "landsskog_agg"))


## ----tree-age-type, eval = T, echo = F, fig.pos = "!b", fig.cap = "Dominerende treslag og skogsalder (eldste trer) i den forenklete landsskogstakseringeren gjennom årene 2020-2022. Punktene er noe spredt ut for å kunne skille mellom region (x-akse) og overlappende verdier (y-akse).", out.width = "0.9\\textwidth"----


tree_agg %>% 
  collect() %>% 
  mutate(dom_tree = ifelse(dom_tree == "Dunbjerk", "Dunbjørk", dom_tree)) %>% 
  mutate(dom_tree = ifelse(dom_tree == "Hengebjerk", "Hengebjørk", dom_tree)) %>% 
  mutate(dom_tree = factor(dom_tree, levels = c("Dunbjørk",
                                                "Hengebjørk",
                                                "Bøk",
                                                "Furu",
                                                "Gran",
                                                "Selje"))) %>% 
  left_join(loc_reg,
            by = c("locality" = "locality"),
            copy = TRUE) %>% 
  group_by(region_name) %>% 
  ggplot(.) +
  geom_point(aes(x = year,
                 y = avg_age,
                 col = dom_tree),
             alpha = 0.7,
             size = 3,
             position = position_jitterdodge(dodge.width = 0.4, 
                                             jitter.height = 2,
                                             jitter.width = 0)) +
  scale_x_continuous(breaks= seq(2020, 2022)) +
  scale_color_nina(name = "Dominerende\ntreslag") +
   guides(color = guide_legend(override.aes = list(size=5))) +
  xlab("År") +
  ylab("Skogsalder")




## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
locality_sampling <- tbl(con,
          Id(schema = "events",
             table = "locality_sampling")) %>% 
  select(sampling_name,
         start_date,
         end_date) %>% 
  collect()



## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
biomass_raw <- get_biomass(agg_level = "locality_sampling",
                           trap_type = "MF"
                           ) 


biomass_mf_locality_sampling_time <- biomass_raw %>% 
  left_join(locality_sampling,
            by = c("sampling_name" = "sampling_name"),
            copy = T) %>% 
  mutate(julian_day = lubridate::yday(end_date),
         habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type),
         start_month = lubridate::month(start_date))  %>% 
  filter(sum_wet_weight > 0)




## ----biomass-mf-julian-day-nor, fig.pos = "!h" , fig.cap = "Tidstrend over innsamlet biomasse (i malaisefeller). \\label{fig:biomass-mf-julian-day-nor}", out.width = "0.7\\textwidth", fig.height = 6, fig.width = 4----
biomass_mf_locality_sampling_time %>% 
  mutate(series = factor(interaction(habitat_type, 
                               region_name,
                               sep = "_"),
                         levels = c("Semi-nat_Østlandet",
                                    "Skog_Østlandet",
                                    "Semi-nat_Trøndelag",
                                    "Semi-nat_Sørlandet"))) %>% 
ggplot(.,
       aes(y = sum_wet_weight/no_trap_days,
           x = as.Date(julian_day, origin = as.Date("2010-01-01")),
           group = series,
           color = series)) +
  geom_line(stat = "smooth",
            formula = y ~ x,
            cex = 2,
            alpha = 0.8) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_nina(name = "Tidsserie") +
  facet_wrap(~year,
             ncol = 1) +
  ylab("Biomasse per felledøgn (gram)") + 
  xlab("Tømmedato") +
  scale_x_date(date_labels = "%d.%b") +
  theme(legend.position="bottom",
        legend.text = element_text(size = 7)) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) 
  


## ----semi-nat-sor-biomass, fig.pos = "!b", fig.cap = "Variasjon i biomasse over tid for semi-naturlig mark i Sørlandet \\label{fig:semi-nat-sor-biomass}", fig.height = 4----

xlims <- biomass_mf_locality_sampling_time %>% 
  filter(habitat_type == "Skog",
         region_name == "Østlandet") %>% 
  select(julian_day)  %>% 
  pull()

xlims <- c(min(xlims, na.rm = T), 
           max(xlims, na.rm = T))
  
ylims <- biomass_mf_locality_sampling_time %>% 
    filter(habitat_type == "Skog",
         region_name == "Østlandet") %>% 
  mutate(amount = sum_wet_weight/no_trap_days) %>% 
  pull()

ylims <- c(min(ylims, na.rm = T)*1.1, 
           max(ylims, na.rm = T)*1.1)  



p1 <- biomass_mf_locality_sampling_time %>% 
  filter(habitat_type == "Semi-nat",
         region_name == "Sørlandet",
         year == 2022) %>% 
  ggplot(.)+
  geom_line(aes(y = sum_wet_weight/no_trap_days,
           x = as.Date(end_date),
           color = locality
           ),
           lwd = 0.75) +
  geom_point(aes(y = sum_wet_weight/no_trap_days,
           x = as.Date(end_date),
           color = locality
           ),
           size = 2.5) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Biomasse per felledøgn (gram))") + 
  scale_x_date(date_labels = "%d.%m.%y",
               limits = as.Date(xlims, origin = as.Date("2022-01-01"))) +
  #scale_y_continuous(limits = ylims) +
  xlab("Tømmedato") +
  scale_color_nina(name = "Lokalitet") +
theme(legend.position = "bottom",
       # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
       legend.key.size = unit(0.5, 'cm'),
       legend.text = element_text(size=7),
       legend.key.width = unit(0.5, 'cm'))  +
  annotate("text", x = as.Date("2022-10-1"), y = 6, label = "bold('2022')", parse = TRUE) 

p1



## ----skog-ost-biomass, fig.cap = "Variasjon i biomasse over tid for skogslokaliteter på Østlandet. \\label{fig:skog-ost-biomass}", fig.height = 8----

p1 <- biomass_mf_locality_sampling_time %>% 
  filter(habitat_type == "Skog",
         region_name == "Østlandet",
         year == 2020) %>% 
  ggplot(.)+
  geom_line(aes(y = sum_wet_weight/no_trap_days,
           x = as.Date(end_date),
           color = locality
           ),
           lwd = 0.75) +
  geom_point(aes(y = sum_wet_weight/no_trap_days,
           x = as.Date(end_date),
           color = locality
           ),
           size = 2.5) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Biomasse per\nfelledøgn (gram)") + 
  scale_x_date(date_labels = "%d.%m.%y",
               limits = as.Date(xlims, origin = as.Date("2020-01-01"))) +
  scale_y_continuous(limits = ylims) +
  xlab("Tømmedato") +
  scale_color_nina(name = "Lokalitet") +
theme(legend.position = "bottom",
       # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
       legend.key.size = unit(0.5, 'cm'),
       legend.text = element_text(size=7),
       legend.key.width = unit(0.5, 'cm'))  +
  annotate("text", x = as.Date("2020-10-1"), y = 4, label = "bold('2020')", parse = TRUE)



p2 <- biomass_mf_locality_sampling_time %>% 
  filter(habitat_type == "Skog",
         region_name == "Østlandet",
         year == 2021) %>% 
  ggplot(.) +
  geom_line(aes(y = sum_wet_weight/no_trap_days,
           x = as.Date(end_date),
           color = locality),
           lwd = 0.75) +
  geom_point(aes(y = sum_wet_weight/no_trap_days,
           x = as.Date(end_date),
           color = locality),
           size = 2.5) +
    scale_x_date(date_labels = "%d.%m.%y",
               limits = as.Date(xlims, origin = as.Date("2021-01-01"))
               ) + 
  scale_y_continuous(limits = ylims) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Biomasse per\nfelledøgn (gram)") + 
  xlab("Tømmedato") +
  scale_color_nina(name = "Lokalitet") +
theme(legend.position = "bottom",
       # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
       legend.key.size = unit(0.5, 'cm'),
       legend.text = element_text(size=7),
       legend.key.width = unit(0.5, 'cm')) +
  annotate("text", x = as.Date("2021-10-1"), y = 4, label = "bold('2021')", parse = TRUE)

p3 <- biomass_mf_locality_sampling_time %>% 
  filter(habitat_type == "Skog",
         region_name == "Østlandet",
         year == 2022) %>% 
  ggplot(.) +
  geom_line(aes(y = sum_wet_weight/no_trap_days,
           x = as.Date(end_date),
           color = locality),
           lwd = 0.75) +
  geom_point(aes(y = sum_wet_weight/no_trap_days,
           x = as.Date(end_date),
           color = locality),
           size = 2.5) +
    scale_x_date(date_labels = "%d.%m.%y",
               limits = as.Date(xlims, origin = as.Date("2022-01-01"))
               ) + 
  scale_y_continuous(limits = ylims) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Biomasse per\nfelledøgn (gram)") + 
  xlab("Tømmedato") +
  scale_color_nina(name = "Lokalitet") +
theme(legend.position = "bottom",
       # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
       legend.key.size = unit(0.5, 'cm'),
       legend.text = element_text(size=7),
       legend.key.width = unit(0.5, 'cm')) +
  annotate("text", x = as.Date("2022-10-1"), y = 4, label = "bold('2022')", parse = TRUE)


marrangeGrob(list(p1,
             p2,
             p3),
             top = "",
             ncol = 1,
             nrow = 3)



## ----semi-nat-ost-biomass,fig.pos = "!b",  fig.cap = "Variasjon i biomasse over tid for semi-naturlig mark på Østlandet. \\label{fig:semi-nat-ost-biomass}", fig.height = 8----

xlims <- biomass_mf_locality_sampling_time %>% 
  filter(habitat_type == "Semi-nat",
         region_name == "Østlandet") %>% 
  select(julian_day)  %>% 
  pull()

xlims <- c(min(xlims, na.rm = T), 
           max(xlims, na.rm = T))
  
ylims <- biomass_mf_locality_sampling_time %>% 
    filter(habitat_type == "Semi-nat",
         region_name == "Østlandet") %>% 
  mutate(amount = sum_wet_weight/no_trap_days) %>% 
  pull()

ylims <- c(min(ylims, na.rm = T)*1.1, 
           max(ylims, na.rm = T)*1.1)  

p1 <- biomass_mf_locality_sampling_time %>% 
  filter(habitat_type == "Semi-nat",
         region_name == "Østlandet",
         year == 2020) %>% 
  ggplot(.)+
  geom_line(aes(y = sum_wet_weight/no_trap_days,
           x = as.Date(end_date),
           color = locality
           ),
           lwd = 0.75) +
  geom_point(aes(y = sum_wet_weight/no_trap_days,
           x = as.Date(end_date),
           color = locality
           ),
           size = 2.5) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Biomasse per\nfelledøgn (gram)") + 
  scale_x_date(date_labels = "%d.%m.%y",
               limits = as.Date(xlims, origin = as.Date("2020-01-01"))
               ) + 
  scale_y_continuous(limits = ylims) +
  xlab("Tømmedato") +
  scale_color_nina(name = "Lokalitet") +
theme(legend.position = "bottom",
       # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
       legend.key.size = unit(0.5, 'cm'),
       legend.text = element_text(size=7),
       legend.key.width = unit(0.5, 'cm')) +
  annotate("text", x = as.Date("2020-10-1"), y = 8, label = "bold('2020')", parse = TRUE)



p2 <- biomass_mf_locality_sampling_time %>% 
  filter(habitat_type == "Semi-nat",
         region_name == "Østlandet",
         year == 2021) %>% 
  ggplot(.) +
  geom_line(aes(y = sum_wet_weight/no_trap_days,
           x = as.Date(end_date),
           color = locality),
           lwd = 0.75) +
  geom_point(aes(y = sum_wet_weight/no_trap_days,
           x = as.Date(end_date),
           color = locality),
           size = 2.5) +
    scale_x_date(date_labels = "%d.%m.%y",
               limits = as.Date(xlims, origin = as.Date("2021-01-01"))
               ) + 
  scale_y_continuous(limits = ylims) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Biomasse per\nfelledøgn (gram)") + 

  xlab("Tømmedato") +
  scale_color_nina(name = "Lokalitet") +
theme(legend.position = "bottom",
       # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
       legend.key.size = unit(0.5, 'cm'),
       legend.text = element_text(size=7),
       legend.key.width = unit(0.5, 'cm'))  +
  annotate("text", x = as.Date("2021-10-1"), y = 8, label = "bold('2021')", parse = TRUE)

p3 <- biomass_mf_locality_sampling_time %>% 
  filter(habitat_type == "Semi-nat",
         region_name == "Østlandet",
         year == 2022) %>% 
  ggplot(.) +
  geom_line(aes(y = sum_wet_weight/no_trap_days,
           x = as.Date(end_date),
           color = locality),
           lwd = 0.75) +
  geom_point(aes(y = sum_wet_weight/no_trap_days,
           x = as.Date(end_date),
           color = locality),
           size = 2.5) +
    scale_x_date(date_labels = "%d.%m.%y",
               limits = as.Date(xlims, origin = as.Date("2022-01-01"))
               ) + 
  scale_y_continuous(limits = ylims) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Biomasse per\nfelledøgn (gram)") + 

  xlab("Tømmedato") +
  scale_color_nina(name = "Lokalitet") +
theme(legend.position = "bottom",
       # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
       legend.key.size = unit(0.5, 'cm'),
       legend.text = element_text(size=7),
       legend.key.width = unit(0.5, 'cm'))  +
  annotate("text", x = as.Date("2022-10-1"), y = 8, label = "bold('2022')", parse = TRUE)



marrangeGrob(list(p1,
             p2,
             p3),
             top = "",
             ncol = 1,
             nrow = 3)



## ----semi-nat-trond-biomass, fig.cap = "Variasjon i biomasse over tid for semi-naturlig mark i Trøndelag. \\label{fig:semi-nat-trond-biomass}", fig.height = 8----
xlims <- biomass_mf_locality_sampling_time %>% 
  filter(habitat_type == "Semi-nat",
         region_name == "Trøndelag") %>% 
  select(julian_day)  %>% 
  pull()

xlims <- c(min(xlims, na.rm = T), 
           max(xlims, na.rm = T))
  
ylims <- biomass_mf_locality_sampling_time %>% 
    filter(habitat_type == "Semi-nat",
         region_name == "Trøndelag") %>% 
  mutate(amount = sum_wet_weight/no_trap_days) %>% 
  pull()

ylims <- c(min(ylims, na.rm = T)*1.1, 
           max(ylims, na.rm = T)*1.1)  

p1 <- biomass_mf_locality_sampling_time %>% 
  filter(habitat_type == "Semi-nat",
         region_name == "Trøndelag",
         year == 2021) %>% 
  ggplot(.)+
  geom_line(aes(y = sum_wet_weight/no_trap_days,
           x = as.Date(end_date),
           color = locality
           ),
           lwd = 0.75) +
  geom_point(aes(y = sum_wet_weight/no_trap_days,
           x = as.Date(end_date),
           color = locality
           ),
           size = 2.5) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Biomasse per felledøgn (gram))") + 
  scale_x_date(date_labels = "%d.%m.%y",
               limits = as.Date(xlims, origin = as.Date("2021-01-01"))) +
  scale_y_continuous(limits = ylims) +
  xlab("Tømmedato") +
  scale_color_nina(name = "Lokalitet") +
theme(legend.position = "bottom",
       # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
       legend.key.size = unit(0.5, 'cm'),
       legend.text = element_text(size=7),
       legend.key.width = unit(0.5, 'cm'))  +
  annotate("text", x = as.Date("2021-10-1"), y = 12, label = "bold('2021')", parse = TRUE) 

p2 <- biomass_mf_locality_sampling_time %>% 
  filter(habitat_type == "Semi-nat",
         region_name == "Trøndelag",
         year == 2022) %>% 
  ggplot(.)+
  geom_line(aes(y = sum_wet_weight/no_trap_days,
           x = as.Date(end_date),
           color = locality
           ),
           lwd = 0.75) +
  geom_point(aes(y = sum_wet_weight/no_trap_days,
           x = as.Date(end_date),
           color = locality
           ),
           size = 2.5) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Biomasse per felledøgn (gram))") + 
  scale_x_date(date_labels = "%d.%m.%y",
               limits = as.Date(xlims, origin = as.Date("2022-01-01"))) +
  scale_y_continuous(limits = ylims) +
  xlab("Tømmedato") +
  scale_color_nina(name = "Lokalitet") +
theme(legend.position = "bottom",
       # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
       legend.key.size = unit(0.5, 'cm'),
       legend.text = element_text(size=7),
       legend.key.width = unit(0.5, 'cm'))  +
  annotate("text", x = as.Date("2022-10-1"), y = 12, label = "bold('2022')", parse = TRUE) 

marrangeGrob(list(p1,
             p2),
             top = "",
             ncol = 1,
             nrow = 2)



## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
locality_sampling_biodiv_sampling_time <- locality_sampling_biodiv %>% 
    left_join(locality_sampling,
            by = c("sampling_name" = "sampling_name"),
            copy = T) %>% 
  mutate(julian_day = lubridate::yday(end_date),
         habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type)) 
  


## ----diversity-mf-julian-day-nor, fig.pos = "!h", fig.cap = "Tidstrend over artsantall (i malaisefeller). \\label{fig:diversity-mf-julian-day-nor}", out.width = "0.7\\textwidth", fig.height = 6, fig.width = 4----
locality_sampling_biodiv_sampling_time %>% 
    mutate(series = factor(interaction(habitat_type, 
                               region_name,
                               sep = "_"),
                         levels = c("Semi-nat_Østlandet",
                                    "Skog_Østlandet",
                                    "Semi-nat_Trøndelag",
                                    "Semi-nat_Sørlandet"))) %>% 
ggplot(.,
       aes(y = no_species/no_trap_days,
           x = as.Date(julian_day, origin = as.Date("2010-01-01")),
           group = series,
           color = series)) +
  geom_smooth(formula = y ~ x,
             cex = 2) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_nina(name = "Tidsserie") +
  facet_wrap(~year,
             ncol = 1) +
  ylab("Artsantall per felledøgn") + 
  xlab("Tømmedato") +
  scale_x_date(date_labels = "%d.%b") +
  theme(legend.position="bottom",
        legend.text = element_text(size = 7)) +
  guides(color = guide_legend(nrow = 2, 
                            byrow = TRUE)
         ) 
  


## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
year_locality_div <- get_observations(subset_region = "Østlandet",
                                 agg_level = "year_locality",
                                 trap_type = "MF") %>% 
  mutate(habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type))



## ----semi_nat_sor_div, fig.pos = "!b", fig.cap = "Variasjon i artsantall over tid for semi-naturlig mark i Sørlandet \\label{fig:semi_nat_sor_div}", fig.height = 4----

xlims <- locality_sampling_biodiv_sampling_time %>% 
  filter(habitat_type == "Skog",
         region_name == "Østlandet") %>% 
  select(julian_day)  %>% 
  pull()

xlims <- c(min(xlims, na.rm = T), 
           max(xlims, na.rm = T))
  
ylims <- locality_sampling_biodiv_sampling_time %>% 
    filter(habitat_type == "Skog",
         region_name == "Østlandet") %>% 
  mutate(amount = no_species/no_trap_days) %>% 
  pull()

ylims <- c(min(ylims, na.rm = T)*1.1, 
           max(ylims, na.rm = T)*1.1)  

 locality_sampling_biodiv_sampling_time %>% 
  filter(habitat_type == "Semi-nat",
         region_name == "Sørlandet",
         year == 2022) %>% 
  ggplot(.) +
  geom_line(aes(y = no_species/no_trap_days,
           x = as.Date(end_date),
           color = locality),
           lwd = 0.75) +
  geom_point(aes(y = no_species/no_trap_days,
           x = as.Date(end_date),
           color = locality),
           size = 2.5) +
    scale_x_date(date_labels = "%d.%m.%y",
               limits = as.Date(xlims, origin = as.Date("2022-01-01"))
               ) + 
  scale_y_continuous(limits = ylims) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Artsantall per\nfelledøgn") + 

  xlab("Tømmedato") +
  scale_color_nina(name = "Lokalitet") +
theme(legend.position = "bottom",
       # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
       legend.key.size = unit(0.5, 'cm'),
       legend.text = element_text(size=7),
       legend.key.width = unit(0.5, 'cm')) +
  annotate("text", x = as.Date("2022-10-1"), y = 100, label = "bold('2022')", parse = TRUE)


## ----skog-ost-div, fig.pos = "!h", fig.cap = "Variasjon i artsantall over tid for skogslokaliteter på Østlandet. \\label{fig:skog-ost-div}", fig.height = 8----

p1 <- locality_sampling_biodiv_sampling_time %>% 
  filter(habitat_type == "Skog",
         region_name == "Østlandet",
         year == 2020) %>% 
  ggplot(.)+
  geom_line(aes(y = no_species/no_trap_days,
           x = as.Date(end_date),
           color = locality
           ),
           lwd = 0.75) +
  geom_point(aes(y = no_species/no_trap_days,
           x = as.Date(end_date),
           color = locality
           ),
           size = 2.5) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Artsantall per felledøgn") + 
  scale_x_date(date_labels = "%d.%m.%y",
               limits = as.Date(xlims, origin = as.Date("2020-01-01"))
               ) + 
  scale_y_continuous(limits = ylims) +
  xlab("Tømmedato") +
  scale_color_nina(name = "Lokalitet") +
theme(legend.position = "bottom",
       # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
       legend.key.size = unit(0.5, 'cm'),
       legend.text = element_text(size=7),
       legend.key.width = unit(0.5, 'cm')) +
  annotate("text", x = as.Date("2020-10-1"), y = 100, label = "bold('2020')", parse = TRUE) 



p2 <- locality_sampling_biodiv_sampling_time %>% 
  filter(habitat_type == "Skog",
         region_name == "Østlandet",
         year == 2021) %>% 
  ggplot(.) +
  geom_line(aes(y = no_species/no_trap_days,
           x = as.Date(end_date),
           color = locality),
           lwd = 0.75) +
  geom_point(aes(y = no_species/no_trap_days,
           x = as.Date(end_date),
           color = locality),
           size = 2.5) +
    scale_x_date(date_labels = "%d.%m.%y",
               limits = as.Date(xlims, origin = as.Date("2021-01-01"))
               ) + 
  scale_y_continuous(limits = ylims) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Artsantall per felledøgn") + 

  xlab("Tømmedato") +
  scale_color_nina(name = "Lokalitet") +
theme(legend.position = "bottom",
       # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
       legend.key.size = unit(0.5, 'cm'),
       legend.text = element_text(size=7),
       legend.key.width = unit(0.5, 'cm')) +
  annotate("text", x = as.Date("2021-10-1"), y = 100, label = "bold('2021')", parse = TRUE) 

p3 <- locality_sampling_biodiv_sampling_time %>% 
  filter(habitat_type == "Skog",
         region_name == "Østlandet",
         year == 2022) %>% 
  ggplot(.) +
  geom_line(aes(y = no_species/no_trap_days,
           x = as.Date(end_date),
           color = locality),
           lwd = 0.75) +
  geom_point(aes(y = no_species/no_trap_days,
           x = as.Date(end_date),
           color = locality),
           size = 2.5) +
    scale_x_date(date_labels = "%d.%m.%y",
               limits = as.Date(xlims, origin = as.Date("2022-01-01"))
               ) + 
  scale_y_continuous(limits = ylims) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Artsantall per felledøgn") + 

  xlab("Tømmedato") +
  scale_color_nina(name = "Lokalitet") +
theme(legend.position = "bottom",
       # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
       legend.key.size = unit(0.5, 'cm'),
       legend.text = element_text(size=7),
       legend.key.width = unit(0.5, 'cm')) +
  annotate("text", x = as.Date("2022-10-1"), y = 100, label = "bold('2022')", parse = TRUE) 


marrangeGrob(list(p1,
                  p2,
                  p3),
             top = "",
             ncol = 1,
             nrow = 3)



## ----semi-nat-div, fig.cap = "Variasjon i artsantall over tid for semi-naturlig mark på Østlandet. \\label{fig:semi-nat-div}", fig.height = 8----
xlims <- locality_sampling_biodiv_sampling_time %>% 
  filter(habitat_type == "Semi-nat",
         region_name == "Østlandet") %>% 
  select(julian_day)  %>% 
  pull()

xlims <- c(min(xlims, na.rm = T), 
           max(xlims, na.rm = T))
  
ylims <- locality_sampling_biodiv_sampling_time %>% 
    filter(habitat_type == "Semi-nat",
         region_name == "Østlandet") %>% 
  mutate(amount = no_species/no_trap_days) %>% 
  pull()

ylims <- c(min(ylims, na.rm = T)*1.1, 
           max(ylims, na.rm = T)*1.1)  

p1 <- locality_sampling_biodiv_sampling_time %>% 
  filter(habitat_type == "Semi-nat",
         region_name == "Østlandet",
         year == 2020) %>% 
  ggplot(.)+
  geom_line(aes(y = no_species/no_trap_days,
           x = as.Date(end_date),
           color = locality
           ),
           lwd = 0.75) +
  geom_point(aes(y = no_species/no_trap_days,
           x = as.Date(end_date),
           color = locality
           ),
           size = 2.5) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Artsantall per felledøgn") + 
  scale_x_date(date_labels = "%d.%m.%y",
               limits = as.Date(xlims, origin = as.Date("2020-01-01"))
               ) + 
  scale_y_continuous(limits = ylims) +
  xlab("Tømmedato") +
  scale_color_nina(name = "Lokalitet") +
theme(legend.position = "bottom",
       # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
       legend.key.size = unit(0.5, 'cm'),
       legend.text = element_text(size=7),
       legend.key.width = unit(0.5, 'cm')) +
  annotate("text", x = as.Date("2020-10-1"), y = 100, label = "bold('2020')", parse = TRUE)  



p2 <- locality_sampling_biodiv_sampling_time %>% 
  filter(habitat_type == "Semi-nat",
         region_name == "Østlandet",
         year == 2021) %>% 
  ggplot(.) +
  geom_line(aes(y = no_species/no_trap_days,
           x = as.Date(end_date),
           color = locality),
           lwd = 0.75) +
  geom_point(aes(y = no_species/no_trap_days,
           x = as.Date(end_date),
           color = locality),
           size = 2.5) +
    scale_x_date(date_labels = "%d.%m.%y",
               limits = as.Date(xlims, origin = as.Date("2021-01-01"))
               ) + 
  scale_y_continuous(limits = ylims) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Artsantall per felledøgn") + 

  xlab("Tømmedato") +
  scale_color_nina(name = "Lokalitet") +
theme(legend.position = "bottom",
       # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
       legend.key.size = unit(0.5, 'cm'),
       legend.text = element_text(size=7),
       legend.key.width = unit(0.5, 'cm')) +
  annotate("text", x = as.Date("2021-10-1"), y = 100, label = "bold('2021')", parse = TRUE) 


p3 <- locality_sampling_biodiv_sampling_time %>% 
  filter(habitat_type == "Semi-nat",
         region_name == "Østlandet",
         year == 2022) %>% 
  ggplot(.) +
  geom_line(aes(y = no_species/no_trap_days,
           x = as.Date(end_date),
           color = locality),
           lwd = 0.75) +
  geom_point(aes(y = no_species/no_trap_days,
           x = as.Date(end_date),
           color = locality),
           size = 2.5) +
    scale_x_date(date_labels = "%d.%m.%y",
               limits = as.Date(xlims, origin = as.Date("2022-01-01"))
               ) + 
  scale_y_continuous(limits = ylims) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Artsantall per felledøgn") + 

  xlab("Tømmedato") +
  scale_color_nina(name = "Lokalitet") +
theme(legend.position = "bottom",
       # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
       legend.key.size = unit(0.5, 'cm'),
       legend.text = element_text(size=7),
       legend.key.width = unit(0.5, 'cm')) +
  annotate("text", x = as.Date("2022-10-1"), y = 100, label = "bold('2022')", parse = TRUE) 

marrangeGrob(list(p1,
                  p2,
                  p3),
             top = "",
             ncol = 1,
             nrow = 3)


## ----semi-nat-trond-div, fig.cap = "Variasjon i artsantall over tid for semi-naturlig mark i Trøndelag. \\label{fig:semi-nat-trond-div}", fig.height = 8----
p1 <- locality_sampling_biodiv_sampling_time %>% 
  filter(habitat_type == "Semi-nat",
         region_name == "Trøndelag",
         year == 2021) %>% 
  ggplot(.) +
  geom_line(aes(y = no_species/no_trap_days,
           x = as.Date(end_date),
           color = locality),
           lwd = 0.75) +
  geom_point(aes(y = no_species/no_trap_days,
           x = as.Date(end_date),
           color = locality),
           size = 2.5) +
    scale_x_date(date_labels = "%d.%m.%y",
               limits = as.Date(xlims, origin = as.Date("2021-01-01"))
               ) + 
  scale_y_continuous(limits = ylims) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Artsantall per\nfelledøgn") + 

  xlab("Tømmedato") +
  scale_color_nina(name = "Lokalitet") +
theme(legend.position = "bottom",
       # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
       legend.key.size = unit(0.5, 'cm'),
       legend.text = element_text(size=7),
       legend.key.width = unit(0.5, 'cm')) +
  annotate("text", x = as.Date("2021-10-1"), y = 100, label = "bold('2021')", parse = TRUE) 


p2 <- locality_sampling_biodiv_sampling_time %>% 
  filter(habitat_type == "Semi-nat",
         region_name == "Trøndelag",
         year == 2022) %>% 
  ggplot(.) +
  geom_line(aes(y = no_species/no_trap_days,
           x = as.Date(end_date),
           color = locality),
           lwd = 0.75) +
  geom_point(aes(y = no_species/no_trap_days,
           x = as.Date(end_date),
           color = locality),
           size = 2.5) +
    scale_x_date(date_labels = "%d.%m.%y",
               limits = as.Date(xlims, origin = as.Date("2022-01-01"))
               ) + 
  scale_y_continuous(limits = ylims) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Artsantall per\nfelledøgn") + 

  xlab("Tømmedato") +
  scale_color_nina(name = "Lokalitet") +
theme(legend.position = "bottom",
       # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
       legend.key.size = unit(0.5, 'cm'),
       legend.text = element_text(size=7),
       legend.key.width = unit(0.5, 'cm')) +
  annotate("text", x = as.Date("2022-10-1"), y = 100, label = "bold('2022')", parse = TRUE)

marrangeGrob(list(p1,
                  p2),
             top = "",
             ncol = 1,
             nrow = 2)




## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
pollinator_table <- tibble("Norsk navn" = c("Gravebier",
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
                                    "Metallmerker (uoffisiell)"),
                           Familie = c("Andrenidae", 
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
                                      "Riodinidae")
                  )


pollinators <- pollinator_table %>% 
  select(Familie) %>% 
  pull()



## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
locality_sampling_pollinators <- get_observations(subset_families = pollinators,
                                             agg_level = "locality_sampling",
                                             trap_type = "MF")


## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
locality_sampling_pollinators_sampling_time <- locality_sampling_pollinators %>% 
    left_join(locality_sampling,
            by = c("sampling_name" = "sampling_name"),
            copy = T) %>% 
  mutate(julian_day = lubridate::yday(end_date),
         habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type)) 
  


## ----pollinators-mf-julian-day-nor, fig.pos = "!h", fig.cap = "Tidstrend over artsantall pollinatorer (i malaisefeller). Pollinatorer regnes her som arter i familiene Gravebier, Langtungebier,  Korttungebier, Markbier, Buksamlerbier, Blomsterbier, Blomfluer, Smygere, Glansvinger, Nymfevinger, Svalestjerter, Hvitvinger, Metallmerker (uoffisiell). \\label{fig:pollinators-mf-julian-day-nor}", out.width = "0.7\\textwidth", fig.height = 6, fig.width = 4----
locality_sampling_pollinators_sampling_time %>% 
      mutate(series = factor(interaction(habitat_type, 
                               region_name,
                               sep = "_"),
                         levels = c("Semi-nat_Østlandet",
                                    "Skog_Østlandet",
                                    "Semi-nat_Trøndelag",
                                    "Semi-nat_Sørlandet"))) %>% 
ggplot(.,
       aes(y = no_species/no_trap_days,
           x = as.Date(julian_day, origin = as.Date("2010-01-01")),
           group = series,
           color = series)) +
  geom_smooth(formula = y ~ x,
             cex = 2) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_nina(name = "Tidsserie") +
  facet_wrap(~year,
             ncol = 1) +
  ylab("Artsantall per felledøgn") + 
  xlab("Tømmedato") +
  scale_x_date(date_labels = "%d.%b") +
  theme(legend.position="bottom",
        legend.text = element_text(size = 7)) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) 
  


## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
year_locality_pollinators <- get_observations(subset_families = pollinators,
                                         subset_region = "Østlandet",
                                         agg_level = "year_locality",
                                         trap_type = "MF") %>% 
  mutate(habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type))



## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
biomass_mf_year_locality <- get_biomass(trap_type = "MF",
                                            agg_level = "year_locality"
                                          #  ,subset_region = "Østlandet"
                                        )


## ----biomass-year-hab, fig.pos = "!h", fig.cap = "A) Middelverdi av biomasse insekter i malaisefeller per lokalitet og år i skog og semi-naturlig mark. Notere at disse data er innsamlet unner delvis ulike tidsperioder og fra ulike regioner og derfor ikke er direkt sammenlignbare. B) Biomasse insekter på Østlandet i juli-august (i malaisefeller). Innsamlingsinnsatsen for denne periode og region er sammenlignbar, men populasjonstopper kan skille seg noe i tid mellom de ulike årene, og disse måneder representerer ikke hele sesongen. \\label{fig:biomass-year-hab}", out.width = "\\textwidth"----
year_mean_biomass <- biomass_mf_year_locality %>% 
    filter(region_name == 'Østlandet') %>% 
  mutate(habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type)) %>% 
  group_by(year,
           habitat_type) %>% 
  summarise(mean_biomass = mean(avg_wet_weight),
            sd_biomass = sd(avg_wet_weight)) %>% 
  mutate(year = as.factor(year),
         habitat_type = as.factor(habitat_type)) %>% 
  ungroup()

dodge = 0.2

p1 <-  ggplot(year_mean_biomass,
         aes(x = year,
             y = mean_biomass,
             group = habitat_type)) +
    geom_point(aes(x = year,
               y = mean_biomass,
               color = habitat_type),
               cex = 3,
               position=position_dodge(dodge)) +
    geom_errorbar(aes(x = year,
                      ymin = mean_biomass - sd_biomass,
                      ymax = mean_biomass + sd_biomass,
                      color = habitat_type),
                   width = .2,
                  lwd = 1,
                  position=position_dodge(dodge)) +
    geom_line(aes(x = year,
                  y = mean_biomass,
                  group = habitat_type),
              position=position_dodge(dodge)) +
    scale_color_nina(name = "Habitat",
                     palette = "darkblue-orange") +
    ylab("Middelverdi biomasse per sampletilfelle og lokalitet") +
    xlab("År") +
  theme(legend.position = "bottom") +
  ggtitle("A")

 
dodge = 0.2

p2 <- biomass_mf_locality_sampling_time %>%   
  filter(region_name == 'Østlandet') %>% 
  mutate(habitat_no = ifelse(habitat_type == "Semi-nat", "Gressmark", "Skog")) %>% 
  mutate(habitat_no = factor(habitat_no, levels = c("Gressmark", "Skog"))) %>% 
  mutate(year = as.factor(year)) %>% 
  filter(start_month >=7,
         start_month <=8) %>% 
  group_by(year, habitat_type) %>% 
  summarise(mean_biomass = mean(avg_wet_weight),
            sd_biomass = sd(avg_wet_weight)) %>% 
  mutate(year = as.factor(year),
         habitat_type = as.factor(habitat_type)) %>% 
  ungroup() %>% 

  ggplot(.,
         aes(x = year,
             y = mean_biomass,
             group = habitat_type)) +
    geom_point(aes(x = year,
               y = mean_biomass,
               color = habitat_type),
               cex = 3,
               position=position_dodge(dodge)) +
    geom_errorbar(aes(x = year,
                      ymin = mean_biomass - sd_biomass,
                      ymax = mean_biomass + sd_biomass,
                      color = habitat_type),
                   width = .2,
                  lwd = 1,
                  position=position_dodge(dodge)) +
    geom_line(aes(x = year,
                  y = mean_biomass,
                  group = habitat_type),
              position=position_dodge(dodge)) +
    scale_color_nina(name = "Habitat",
                     palette = "darkblue-orange") +
    ylab("") +
    xlab("År") +
    theme(legend.position = "bottom") +
    #ylim(c(0, max(ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range))) +
    ggtitle("B")

marrangeGrob(list(p1, 
                  p2),
             ncol = 2, 
             nrow = 1,
             top = "")



## ----div_year_hab, fig.pos = "!h", fig.cap = "A) Middelverdi av antall arter fanget per lokalitet og år i skog og semi-naturlig mark i alle regioner. B) Middelverdi av antall pollinatorer fanget per lokalitet og år i skog og semi-naturlig mark i alle regioner. Merk at i 2021 ble novaseq-maskinen tatt i bruk for metastrekkoding. Begge figurer viser data kun fra malaisefeller. \\label{fig:div_year_hab}", out.width = "\\textwidth"----
year_mean_div <- year_locality_div %>% 
  group_by(year,
           habitat_type) %>% 
  summarise(mean_richness = mean(no_species),
            sd_richness = sd(no_species)) %>% 
  mutate(year = as.factor(year),
         habitat_type = as.factor(habitat_type)) %>% 
  ungroup()

dodge = 0.2

p1 <-  ggplot(year_mean_div,
         aes(x = year,
             y = mean_richness,
             group = habitat_type)) +
    geom_point(aes(x = year,
               y = mean_richness,
               color = habitat_type),
               cex = 3,
                  position=position_dodge(dodge)) +
     geom_errorbar(aes(x = year,
                      ymin = mean_richness - sd_richness,
                      ymax = mean_richness + sd_richness,
                      color = habitat_type),
                   width = .2,
                  lwd = 1,
                  position=position_dodge(dodge)) +
    geom_line(aes(x = year,
                  y = mean_richness,
                  group = habitat_type),
                  position=position_dodge(dodge)) +
    scale_color_nina(name = "Habitatstype",
                     palette = "darkblue-orange") +
    ylab("Antall arter per lokalitet (middelv.)") +
    xlab("År") +
    theme(legend.position = "bottom")



 
year_mean_pollinators <- year_locality_pollinators %>% 
  group_by(year,
           habitat_type) %>% 
  summarise(mean_richness = mean(no_species),
            sd_richness = sd(no_species)) %>% 
  mutate(year = as.factor(year),
         habitat_type = as.factor(habitat_type)) %>% 
  ungroup()

p2 <- ggplot(year_mean_pollinators,
         aes(x = year,
             y = mean_richness,
             group = habitat_type)) +
    geom_point(aes(x = year,
               y = mean_richness,
               color = habitat_type),
               cex = 3,
               position=position_dodge(dodge)) +
    geom_errorbar(aes(x = year,
                      ymin = mean_richness - sd_richness,
                      ymax = mean_richness + sd_richness,
                      color = habitat_type),
                      width = .2,
                      lwd = 1,
                  position=position_dodge(dodge)) +
    geom_line(aes(x = year,
                  y = mean_richness,
                  group = habitat_type),
                  position=position_dodge(dodge)) +
    scale_color_nina(name = "Habitatstype",
                     palette = "darkblue-orange") +
    ylab("Antall arter pollinatorer per lokalitet (middelv.)") +
    xlab("År") +
    theme(legend.position = "bottom")

marrangeGrob(list(p1, 
                  p2),
             ncol = 2, 
             nrow = 1,
             top = "")




## ----clim-trend-trond, fig.pos = "!hb", fig.cap = "Klimatrend juni-juli i Trøndelag. Øverste figur viser summen av daglig gjennomsnittstemperatur, nedre figur viser summen av daglig nedbørsum. Feillinjer viser +- 1 standardavvik basert på alle lokaliteter i regionen. \\label{fig:clim-trend-trond}", out.width = "0.7\\textwidth"----
clim_trend_plot(region = "Trøndelag",
                from_month = 6,
                to_month = 8) 



## ----clim-trend-ost, fig.cap = "Klimatrend juni-juli på Østlandet. Øverste figur viser summen av daglig gjennomsnittstemperatur, nedre figur viser summen av daglig nedbørsum. Feillinjer viser +- 1 standardavvik basert på alle lokaliteter i regionen.\\label{fig:clim-trend-ost}", out.width = "0.7\\textwidth"----
clim_trend_plot(region = "Østlandet",
                from_month = 6,
                to_month = 8) 



## ----clim-trend-sor, fig.cap = "Klimatrend juni-juli på Sørlandet. Øverste figur viser summen av daglig gjennomsnittstemperatur, nedre figur viser summen av daglig nedbørsum. Feillinjer viser +- 1 standardavvik basert på alle lokaliteter i regionen. ", out.width = "0.7\\textwidth"----
clim_trend_plot(region = "Sørlandet",
                from_month = 6,
                to_month = 8) 



## -----------------------------------------------------------------------------------------------------------------------------------------------
climdata_semi_nat_55 <- get_climate_data(locality = "Semi-nat_55")


## ----semi-nat-55-temp, fig.pos = "!hb", fig.cap = "Temperatur i Trondheim 2022, eksemplifisert av lokalitet Semi-nat\\_55. \\label{fig:trond-temp-2022}", out.width = "0.9\\textwidth"----
plot_climate_comparison(climdata_semi_nat_55,
                        main_title = FALSE) 



## ----semi-nat-55-precip,fig.pos = "!hb", fig.cap = "Nedbør i Trondheim 2022, eksemplifisert av lokalitet Semi-nat\\_55. \\label{fig:trond-precip-2022}", out.width = "0.9\\textwidth"----
plot_climate_comparison(climdata_semi_nat_55,
                        main_title = FALSE,
                        variable = "precipitation",
                        rolling_mean = TRUE) 


## -----------------------------------------------------------------------------------------------------------------------------------------------
climdata_semi_nat_42 <- get_climate_data(locality = "Semi-nat_42")


## ----semi-nat-42-temp, fig.pos = "!hb",  fig.cap = "Temperatur i Lillehammer 2022, eksemplifisert av lokalitet Semi-nat\\_42. \\label{fig:ost-temp-2022}", out.width = "0.9\\textwidth"----
plot_climate_comparison(climdata_semi_nat_42,
                        main_title = FALSE) 



## ----semi-nat-42-precip, fig.pos = "!hb",  fig.cap = "Nedbør i Lillehammer 2022, eksemplifisert av lokalitet Semi-nat\\_42. \\label{fig:ost-precip-2022}", out.width = "0.9\\textwidth"----
plot_climate_comparison(climdata_semi_nat_42,
                        main_title = FALSE,
                        variable = "precipitation",
                        rolling_mean = TRUE) 



## -----------------------------------------------------------------------------------------------------------------------------------------------
climdata_semi_nat_32 <- get_climate_data(locality = "Semi-nat_32")


## ----semi-nat-32-temp, fig.cap = "Temperatur i Kristiansand 2022, eksemplifisert av lokalitet Semi-nat\\_32. \\label{fig:sor-temp-2022}", out.width = "0.9\\textwidth"----
plot_climate_comparison(climdata_semi_nat_32,
                        main_title = FALSE) 



## ----semi-nat-32-precip, fig.cap = "Temperatur i Kristiansand 2022, eksemplifisert av lokalitet Semi-nat\\_32. \\label{fig:sor-precip-2022}", out.width = "0.9\\textwidth"----
plot_climate_comparison(climdata_semi_nat_32,
                        main_title = FALSE,
                        variable = "precipitation",
                        rolling_mean = TRUE) 



## ---- figures-side, fig.pos = "!b", fig.show="hold", out.width="40%", fig.cap = "Temperatur og luftfuktighetsmåler av typen Hobo MX2301A (venstre) og temperatur og lysmåler av typen Hobo MX2202 (høyre). MX2301A er plassert i skygge og begge tømmes med blåtann til en skyløsning. \\label{fig:templogger}"----
include_graphics("figures/templogger.jpg")
include_graphics("figures/lyslogger_small.png")


## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
locality_sampling_loggerdata <- get_logger_data(dataset = "NasIns",
                                                   agg_level = "locality_sampling") 



## ----logger-corr, eval = F, fig.cap = "Korrelasjonsplot over loggerdata. Relativ luftfuktighet (%) måles av loggermodell MX2301A og lys (lux) av loggermodell MX2202. MX2301A er montert unner en lysskjerm og reagerer ikke like hyppig på sollys som en uskjermete loggeren MX2202 som også måler lys."----
## 
## locality_sampling_loggerdata %>%
##   select(Temp_MX2301A = avg_values_MX2301A_temperature,
##          Temp_MX2202 = avg_values_MX2202_temperature,
##          RL = avg_values_MX2301A_rh,
##          Lys = avg_values_MX2202_light
##          ) %>%
##   filter(Lys >= 100) %>%
##   Norimon::chart_correlation(.,
##                   histogram = TRUE,
##                   method = "pearson",
##                   color = ninaColors()[1])
## 


## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
locality_sampling_weight_loggerdata <- biomass_mf_locality_sampling %>% 
  left_join(locality_sampling_loggerdata,
            by = c("sampling_name" = "sampling_name",
                   "year" = "year",
                   "locality" = "locality",
                   "habitat_type" = "habitat_type")) %>% 
  mutate(habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type))




## ----locality-sampling-biomass-loggerdata, fig.cap = "Sammenheng mellom temperatur og innsamlet biomasse per fangstperiode i malaisefeller, her vist ved en enkel linjær funksjon, selv om andre funksjoner er tenkbare.\\label{fig:locality-sampling-biomass-loggerdata}", out.width = "0.8\\textwidth"----

locality_sampling_weight_loggerdata %>% 
  group_by(habitat_type) %>% 
  ggplot(aes(y = avg_wet_weight/no_trap_days, 
             x = avg_values_MX2301A_temperature)) +
  geom_point(aes(color = habitat_type)) +
  scale_color_nina(name = "Habitatstype",
                   palette = "darkblue-orange") +
     geom_smooth(aes(color = habitat_type),
              method = "glm",
              #fullrange = T,
              se = T,
              lwd = 2) +
  ylab("Avrunnen våtvekt (g.) per felledag") +
  xlab("Middeltemperatur ved felle unner fangstperioden")


## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
locality_sampling_biodiv_loggerdata <- locality_sampling_biodiv %>% 
  left_join(locality_sampling_loggerdata,
            by = c("sampling_name" = "sampling_name",
                   "year" = "year",
                   "locality" = "locality",
                   "habitat_type" = "habitat_type")) %>% 
  mutate(habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type))



## ----locality-sampling-biodiv-loggerdata, fig.cap = "Sammenheng mellom temperatur og antall arter innsamlet per fangstperiode i malaisefeller.\\label{fig:locality-sampling-biodiv-loggerdata}", out.width = "0.8\\textwidth"----
locality_sampling_biodiv_loggerdata %>% 
  group_by(habitat_type) %>% 
  ggplot(aes(y = no_species, 
             x = avg_values_MX2301A_temperature)) +
  geom_point(aes(color = habitat_type)) +
  scale_color_nina(name = "Habitatstype",
                   palette = "darkblue-orange") +
     geom_smooth(aes(color = habitat_type),
              method = "lm",
              se = T,
              lwd = 2) +
  ylab("Antall arter") +
  xlab("Middeltemperatur ved felle unner fangstperioden")


## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
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
                                  query = redlisted_obs_1)



redlisted_obs_2021 <- redlisted_obs_2021 %>% 
  #mutate(expl_fact = gsub("([a-zA-Z]*)( >)(.*)", "\\1", expl_fact)) %>% 
  separate(expl_fact, 
           into = c("expl_1",
                    "expl_2",
                    "expl_3",
                    "expl_4",
                    "expl_5"),
           sep = ">") %>% 
  separate(expl_3,
           into = "expl_3_main",
           sep = "_") %>%
  mutate(expl_3_main = ifelse(is.na(expl_3_main), expl_1, expl_3_main)) %>% 
  mutate(expl_3_main = str_trim(expl_3_main)) %>% 
  mutate(expl_3_main = ifelse(expl_3_main == " ", NA, expl_3_main),
         expl_3_main = ifelse(expl_3_main == "", NA, expl_3_main),
         expl_3_main = ifelse(expl_3_main == "Ukjent ", NA, expl_3_main),
         expl_3_main = ifelse(expl_3_main == "Ukjent", NA, expl_3_main)) 
 
redlisted_obs_2021 <- redlisted_obs_2021 %>% 
  group_by(kategori_2021) %>% 
  mutate(no_spec_per_kat = n_distinct(species_latin_fixed)) %>% 
  ungroup() %>% 
  filter(kategori_2021 != "DD")



## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
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
                                  query = redlisted_obs_1_mf)



redlisted_obs_2021_mf <- redlisted_obs_2021_mf %>% 
  #mutate(expl_fact = gsub("([a-zA-Z]*)( >)(.*)", "\\1", expl_fact)) %>% 
  separate(expl_fact, 
           into = c("expl_1",
                    "expl_2",
                    "expl_3",
                    "expl_4",
                    "expl_5"),
           sep = ">") %>% 
  separate(expl_3,
           into = "expl_3_main",
           sep = "_") %>%
  mutate(expl_3_main = ifelse(is.na(expl_3_main), expl_1, expl_3_main)) %>% 
  mutate(expl_3_main = str_trim(expl_3_main)) %>% 
  mutate(expl_3_main = ifelse(expl_3_main == " ", NA, expl_3_main),
         expl_3_main = ifelse(expl_3_main == "", NA, expl_3_main),
         expl_3_main = ifelse(expl_3_main == "Ukjent ", NA, expl_3_main),
         expl_3_main = ifelse(expl_3_main == "Ukjent", NA, expl_3_main)) %>% 
  filter(kategori_2021 != "DD")
 



## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
no_detected_redlist_species <- redlisted_obs_2021 %>% 
  st_drop_geometry() %>% 
  filter(kategori_2021 %in% c("NT",
                              "VU",
                              "EN",
                              "CR",
                              "RE")) %>%
  summarise(no_spec = n_distinct(vitenskapelig_navn_id)) %>% 
  pull()


## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
redlisted_obs_2021_agg <- redlisted_obs_2021 %>% 
  group_by(locality,
           kategori_2021) %>% 
  summarise(no_spec = n_distinct(species_latin_fixed)) 

kat_order <- tibble(kategori_2021 = c("NT",
                                      "VU",
                                      "EN",
                                      "CR",
                                      "RE"),
                    kat_order = 1:5)


no_spec_per_kat_df <- redlisted_obs_2021 %>% 
  st_drop_geometry() %>% 
  select(kategori_2021, 
         no_spec_per_kat) %>% 
  distinct()

redlisted_obs_2021_agg <- redlisted_obs_2021_agg %>% 
  left_join(kat_order,
            by = c("kategori_2021" = "kategori_2021")) %>% 
  ungroup() %>% 
  left_join(no_spec_per_kat_df,
            by = c("kategori_2021" = "kategori_2021")) %>% 
  mutate(kategori_2021 = factor(paste0(kategori_2021, " (", no_spec_per_kat, " st.)"))) %>% 
  mutate(kategori_2021 = reorder(kategori_2021, kat_order))


redlisted_obs_2021_label <- redlisted_obs_2021_agg %>% 
  select(locality) %>% 
  distinct()

redlisted_obs_2021_agg <- redlisted_obs_2021_agg %>% 
  st_jitter(redlisted_obs_2021_agg, amount = 7000) 


## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
no_locales_with_redlisted_species <- redlisted_obs_2021_agg %>% 
  st_drop_geometry() %>% 
  ungroup() %>% 
  summarise(n_distinct(locality)) %>% 
  pull()


## -----------------------------------------------------------------------------------------------------------------------------------------------
#redlist_colors <- c("#5a5b5d","#990013", "#d80024" , "#ff4c29", "#f36c21", "#f9a61a")
#names(redlist_colors) <- c("RE", "CR", "EN", "VU", "NT", "DD")

#redlist_colors <- c("#f9a61a", "#f36c21", "#ff4c29", "#d80024") #2022 colors
redlist_colors <- c("#ed6c26", "#e94f33", "#d80f27", "#5a5b5d")

#redlist_colors <- c("#f7a520", "#ed6c26", "#d80f27", "#5a5b5d") #2023 powertoys color picker tool

names(redlist_colors) <- levels(redlisted_obs_2021_agg$kategori_2021)


## ----redlist-map, fig.cap = "Fordeling av funn av rødlistete arter (VU-RE), navngitt med 'høy' sikkerhet. Sirkelstørrelsen viser antall arter funnet på hver lokalitet, og fargen viser rødlistekategori. Posisjonene er noe justert på kartet for at de ikke skal skjule hverandre. Prikkene med 'RE' består av funn av 'Elachista quadripunctella', som vil bli morfologisk sjekket. \\label{fig:redlist-map}", fig.height = 8----
set.seed(12345)
norge %>% 
  filter(fylke == "Innlandet" |
         fylke =="Vestfold og Telemark" |
         fylke == "Oslo" |
         fylke == "Viken" | 
         fylke == "Trøndelag" |
         fylke == "Rogaland" |
         fylke == "Agder") %>% 
ggplot(.) +
  geom_sf() +
  geom_sf(aes(size = no_spec,
              color = kategori_2021),
              alpha = 0.8,
          data = redlisted_obs_2021_agg) +
  scale_color_manual(name = "Rødlistekategori",
                     values = redlist_colors
                     ) +
  scale_size(breaks = c(1, seq(from = 0,
                          to = max(redlisted_obs_2021_agg$no_spec), 
                          by = 5)),
             name = "Antall arter per\nkategori og lokalitet") +
  xlab("") +
  ylab("") +
  guides(color = guide_legend(override.aes = list(size=5))) +
     theme(legend.position = "bottom",
           legend.box = "vertical",
       # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
       legend.key.size = unit(0.5, 'cm'),
       legend.text = element_text(size=7),
       legend.key.width = unit(0.5, 'cm'))
  
  


## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------

redlisted_obs_2021_agg_loc_expl <- redlisted_obs_2021 %>% 
  group_by(locality,
           expl_3_main) %>% 
  summarise(no_spec_loc_expl = n_distinct(species_latin_fixed))  %>% 
  ungroup() %>% 
  mutate(expl_3_main = ifelse(is.na(expl_3_main), "Ukjent/Ikke angitt", expl_3_main))


redlisted_obs_2021_agg_expl <- redlisted_obs_2021 %>% 
  st_drop_geometry() %>% 
  group_by(expl_3_main) %>% 
  summarise(no_spec_expl = n_distinct(species_latin_fixed))  %>% 
  ungroup() %>% 
  mutate(expl_3_main = ifelse(is.na(expl_3_main), "Ukjent/Ikke angitt", expl_3_main))

kat_order <- tibble(expl_3_main = c("Ukjent/Ikke angitt",
                                         "Buskap/dyrehold",
                                         "Forurensing",
                                         "Fremmede arter",
                                         "Jordbruk",
                                         "Opphørt/redusert drift",
                                         "Påvirker habitatet (beite tråkk mm.)",
                                         "Skogbruk (kommersielt)",
                                         "Utbygging/utvinning"),
                    kat_order = 1:9)

redlisted_obs_2021_agg_expl_big <- redlisted_obs_2021_agg_loc_expl %>% 
  left_join(redlisted_obs_2021_agg_expl,
             by = c("expl_3_main" = "expl_3_main")) %>% 
  left_join(kat_order,
            by = c("expl_3_main" = "expl_3_main")) %>% 
  ungroup() %>% 
  mutate(expl_3_main = factor(paste0(expl_3_main, " (", no_spec_expl, " st.)"))) %>% 
  mutate(expl_3_main = fct_reorder(expl_3_main, kat_order))

#levels(redlisted_obs_2021_agg_expl$expl_3_main)

redlisted_obs_2021_agg_expl_big <- redlisted_obs_2021_agg_expl_big %>% 
  st_jitter(redlisted_obs_2021_agg_expl_big, amount = 7000) 


## ----redlist-map-expl, fig.cap = "Fordeling av funn av rødlistete arter, navngitt med 'høy' sikkerhet, etter påvirkningsfaktor. Sirkelstørrelsen viser antall arter funnet på hver lokalitet, og fargen viser påvirkningsfaktor (flere nivåer) angitt i rødlisten for 2021. Posisjonene er noe justert på kartet for at de ikke skal skjule hverandre.\\label{fig:redlist-map-expl}", fig.height = 8----

redlist_expl_colors <- c(ninaColors()[c(1, 6, 2, 3, 5, 4)], "yellow")
names(redlist_expl_colors) <- levels(redlisted_obs_2021_agg_expl_big$expl_3_main) 

set.seed(12345)
norge %>% 
  filter(fylke == "Innlandet" |
         fylke =="Vestfold og Telemark" |
         fylke == "Oslo" |
         fylke == "Viken" | 
        fylke == "Trøndelag" |
          fylke == "Rogaland" |
          fylke == "Agder" ) %>% 
ggplot(.) +
  geom_sf() +
  geom_sf(aes(size = no_spec_loc_expl,
              color = expl_3_main),
              alpha = 0.7,
          data = redlisted_obs_2021_agg_expl_big) +
  #scale_color_nina(name = "Påvirkning") +
  scale_color_manual(name = "Påvirkning",
                     values = redlist_expl_colors) +
  scale_size(breaks = c(1, seq(from = 0,
                          to = max(redlisted_obs_2021_agg$no_spec), 
                          by = 5)),
             name = "Antall arter per\nkategori og lokalitet") +

  xlab("") +
  ylab("") +
  guides(color = guide_legend(override.aes = list(size=7), nrow = 3, byrow = T)) +
     theme(legend.position = "bottom",
           legend.box = "vertical",
        legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
       legend.key.size = unit(0.4, 'cm'),
       legend.text = element_text(size = 6),
       legend.key.width = unit(0.2, 'cm'))
  
  


## ----redlisted_obs_agg, cache = cache_var-------------------------------------------------------------------------------------------------------
kat_order <- tibble(kategori_2021 = c("DD",
                                      "NT",
                                      "VU",
                                      "EN",
                                      "CR"),
                    kat_order = 1:5)

redlisted_obs_2021_less_agg <-  redlisted_obs_2021 %>% 
  st_drop_geometry() %>% 
  dplyr:::select(species_latin = species_latin_fixed,
         id_order,
         popular_name = populærnavn,
         everything()) %>% 
  group_by(locality,
           kategori_2021,
           id_order,
           sampling_name,
           species_latin,
           popular_name,
           no_spec_per_kat) %>% 
  distinct() %>% 
  group_by(locality,
           kategori_2021,
           id_order,
           species_latin,
           popular_name,
           no_spec_per_kat) %>% 
  summarise(times_detected = n()) %>% 
   ungroup() %>% 
  left_join(kat_order,
            by = c("kategori_2021" = "kategori_2021")) %>% 
  mutate(kategori_2021_append = factor(paste0(kategori_2021, " (", no_spec_per_kat, " st.)"))) %>% 
  mutate(kategori_2021_append = reorder(kategori_2021_append, kat_order)) %>% 
  arrange(locality,
          kategori_2021_append,
          species_latin)




## ----redlist-detection, fig.cap = "Deteksjon av rødlistede arter i hver lokalitet. Hver sirkel viser funn av en enkelt art i en lokalitet. Størrelsen på sirklene viser antall ganger arten ble detektert på lokaliteten (enkelte felleobservasjoner). Prikkene med 'RE' består av funn av 'Elachista quadripunctella', som vil bli morfologisk sjekket. \\label{fig:redlist-detection}", fig.height = 7----
ggplot(redlisted_obs_2021_less_agg,
       aes(y = locality,
           x = kategori_2021_append)) +
  geom_point(aes(size = times_detected,
                 color = kategori_2021_append),
             alpha = 0.7,
             position = position_jitter(height = 0)) +
  scale_y_discrete(limits = rev) +
  #guides(color = "none") +
    scale_color_manual(name = "Rødlistekategori",
                     values = redlist_colors) +
  scale_size(name = "Antall ganger\ndetektert",
             breaks = c(1, 2, 4, 6, 8, 10)) +
  theme(legend.position="right",
        axis.text.y = element_text(size = 6)) +
  xlab("Rødlistekategori 2021") +
  ylab("Lokalitet")



## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------

spec_no_loc <- redlisted_obs_2021_less_agg %>% 
  group_by(kategori_2021,
           id_order,
           species_latin,
           popular_name) %>% 
  summarise(no_lok = n_distinct(locality),
            no_detected = sum(times_detected)) %>% 
  arrange(kategori_2021,
          desc(no_lok),
          desc(no_detected)) %>% 
  dplyr::select(`Rødlistekategori 2021` = kategori_2021,
         `Orden` = id_order,
         `Latinsk navn` = species_latin,
         `Norsk navn` = popular_name,
         `Funnet i antall lokaler` = no_lok,
         `Funnet antall ganger` = no_detected)

spec_no_loc <- spec_no_loc %>% 
  mutate(`Norsk navn` = str_to_title(`Norsk navn`))
  

spec_no_loc1 <- spec_no_loc %>% 
  filter(`Rødlistekategori 2021` == "NT")

spec_no_loc2 <- spec_no_loc %>% 
  filter(`Rødlistekategori 2021` == "VU" |
         `Rødlistekategori 2021` == "EN" | 
         `Rødlistekategori 2021` == "RE" |
         `Rødlistekategori 2021` == "CR" )

names(spec_no_loc1)=c("\\multicolumn{1}{p{0.4in}}{\\raggedright Kategori}",
                      "\\multicolumn{1}{p{0.75in}}{\\raggedright Orden}",
                     "\\multicolumn{1}{p{0.75in}}{\\raggedright Latinsk navn}",
                     "\\multicolumn{1}{p{0.75in}}{\\raggedright Norsk navn}",
                     "\\multicolumn{1}{p{0.75in}}{\\raggedright Funnet i antall lokaliteter}",
                     "\\multicolumn{1}{p{0.75in}}{\\raggedright Antall ganger observert}")

names(spec_no_loc2)=c("\\multicolumn{1}{p{0.4in}}{\\raggedright Kategori}",
                      "\\multicolumn{1}{p{0.75in}}{\\raggedright Orden}",
                     "\\multicolumn{1}{p{0.75in}}{\\raggedright Latinsk navn}",
                     "\\multicolumn{1}{p{0.75in}}{\\raggedright Norsk navn}",
                     "\\multicolumn{1}{p{0.75in}}{\\raggedright Funnet i antall lokaliteter}",
                     "\\multicolumn{1}{p{0.75in}}{\\raggedright Antall ganger observert}")




## ----redlist-table1, results = "asis"-----------------------------------------------------------------------------------------------------------
print(xtable(spec_no_loc1, 
             caption = "Arter i rødlistekategori nær truet (NT) i prosjektet så langt. Funnene er basert på DNA-metastrekkoding og er ikke morfologisk verifiserte. \\label{tab:redlist-table1}",
             align = rep("l", 7)), 
      size = "footnotesize",
      sanitize.colnames.function=function(x){x},
     # table.placement = "p",
      floating.environment = "table*",
      tabular.environment = "longtable",
      floating = F)



## ----redlist-table2, results = "asis"-----------------------------------------------------------------------------------------------------------
print(xtable(spec_no_loc2, 
             caption = "Rødlistede arter funne i prosjektet så langt, fordelt på kategoriene sterkt truet (EN) sårbar (VU), og regionalt utdødd (RE). Funnene er basert på DNA-metastrekkoding og er ikke morfologisk verifiserte. NB at man ikke kan skille Ampedus sanguinolentus fra Ampedus pomonae med de primers vi bruker. Elachista quadripunctella er overraskende og vil bli morfologisk sjekket. \\label{tab:redlist-table2}",
             align = rep("l", 7)), 
      size = "footnotesize",
      sanitize.colnames.function=function(x){x},
     # table.placement = "p",
      floating.environment = "table*",
      tabular.environment = "longtable",
      floating = F)



## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------

no_alien_q <- "
SELECT obs.id_phylum, obs.id_class, obs.id_order, count(distinct obs.species_latin_fixed)
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
GROUP BY id_phylum, id_class, id_order
ORDER BY id_phylum, id_class, id_order

"  

no_alien <- dbGetQuery(con,
                       no_alien_q) %>% 
  mutate(count = as.integer(count)) %>%
  select(everything(),
         no_alien = count)



## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------

no_fennoscand_q <- "
SELECT obs.id_phylum, obs.id_class, obs.id_order, count(distinct obs.species_latin_fixed)
FROM occurrences.observations obs,
lookup.fennoscand_species2 alien,
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
GROUP BY id_phylum, id_class, id_order
ORDER BY id_phylum, id_class, id_order
"  

no_fennoscand <- dbGetQuery(con,
                       no_fennoscand_q) %>%
  mutate(count = as.integer(count)) %>%
  select(everything(),
         no_fennoscand_alien = count)



## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------

no_pot_alien_q <- "
SELECT obs.id_phylum, obs.id_class, obs.id_order, count(distinct obs.species_latin_fixed)
FROM occurrences.observations obs,
lookup.pot_alien_species2 alien,
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
GROUP BY id_phylum, id_class, id_order
ORDER BY id_phylum, id_class, id_order
"  

no_pot_alien <- dbGetQuery(con,
                       no_pot_alien_q) %>%
  mutate(count = as.integer(count)) %>%
  select(everything(),
         no_pot_alien = count)



## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------

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

no_all_species <- dbGetQuery(con,
                       no_all_species_q) %>% 
  mutate(count = as.integer(count)) %>% 
  select(everything(),
         no_all_species = count)



## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------

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

no_all_species_to_species <- dbGetQuery(con,
                       no_all_species_to_species_q) %>% 
  mutate(count = as.integer(count)) %>% 
  select(everything(),
         no_all_species_to_species = count)



## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------

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

no_all_native_species_to_species <- dbGetQuery(con,
                       no_all_native_species_to_species_q) %>% 
  mutate(count = as.integer(count)) %>% 
  select(everything(),
         no_all_native_species_to_species = count)



## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
no_all_aliens <- no_fennoscand %>% 
  full_join(no_pot_alien,
            by = c("id_phylum" = "id_phylum",
                   "id_class" = "id_class",
                   "id_order" = "id_order")) %>% 
   full_join(no_alien,
            by = c("id_phylum" = "id_phylum",
                   "id_class" = "id_class",
                   "id_order" = "id_order")) %>% 
  full_join(no_all_species,
            by = c("id_phylum" = "id_phylum",
                   "id_class" = "id_class",
                   "id_order" = "id_order")) %>% 
  full_join(no_all_species_to_species,
            by = c("id_phylum" = "id_phylum",
                   "id_class" = "id_class",
                   "id_order" = "id_order")) %>% 
  full_join(no_all_native_species_to_species,
            by = c("id_phylum" = "id_phylum",
                   "id_class" = "id_class",
                   "id_order" = "id_order")) %>% 
  arrange(id_phylum,
          id_class,
          id_order) %>% 
  replace(is.na(.), 0) %>% 
  mutate(id_order = ifelse(id_order == 0, "Ukjent", id_order))
  
  
no_all_aliens <- no_all_aliens %>% 
  as_tibble() %>% 
  mutate(no_all_species = as.integer(no_all_species),
         no_all_species_to_species = as.integer(no_all_species_to_species),
         no_all_native_species_to_species = as.integer(no_all_native_species_to_species),
         no_alien = as.integer(no_alien),
         no_fennoscand_alien = as.integer(no_fennoscand_alien),
         no_pot_alien = as.integer(no_pot_alien)) %>% 
  mutate(id_order = ifelse(id_order == "Plecoptera_Insecta", "Plecoptera", id_order)) %>% 
  select(Fylum = "id_phylum",
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
  add_row(Fylum = "", 
          Klasse = "", 
          Orden = "Sum", 
          `Antall takson` = sum(.$`Antall takson`) ,
          `Takson bestemt til art`= sum(.$`Takson bestemt til art`),
          `Finnes i Norge` = sum(.$`Finnes i Norge`),
          `Fremmede arter` = sum(.$`Fremmede arter`),
          `Fennoskandiske arter` = sum(.$`Fennoskandiske arter`),
          `Potensielt fremmede arter` = sum(.$`Potensielt fremmede arter`) 
           )




## ---- results = "asis"--------------------------------------------------------------------------------------------------------------------------
xtab <- xtable(no_all_aliens, caption = "Antall taksa og arter påvist med DNA-metastrekkoding av alle felleprøver i 2020 - 2022 sortert på orden, finnes i Norge (registrert i artsnavnebasen), kjente fremmede arter for Norge (registrert i fremmedartslista), tilstedeværelse i Fennoskandia (i GBIF), og potentielt fremmede arter (se figur 8). Kun artsidentifiseringer med 'høy' konfidanse er tatt med (unntatt 'Antall takson').\\label{tab:no-species-order}")

align(xtab) <- "p{1.5cm}p{1.2cm}p{1.2cm}p{2.2cm}p{1cm}p{1cm}p{1cm}p{1cm}p{1.7cm}p{1cm}"

print.xtable(xtab, sanitize.text.function=function(x){x}, size = "scriptsize")



## -----------------------------------------------------------------------------------------------------------------------------------------------
pot_alien_spec <- tbl(con,
                  Id(schema = "lookup",
                     table = "pot_alien_species2")) %>% 
  select(species_latin_fixed) %>% 
  pull()

fennoscand_spec <- tbl(con,
                  Id(schema = "lookup",
                     table = "fennoscand_species2")) %>% 
  select(species_latin_fixed) %>% 
  pull()

alien_spec <- tbl(con,
                  Id(schema = "lookup",
                     table = "fremmedartslista_2018_artsdatabanken")) %>% 
  filter(expertGroup == "Ikke-marine invertebrater") %>% 
  select(scientificName) %>%
  pull()


## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------

fennoskand_obs_q <- "
SELECT 
obs.species_latin_fixed,
yl.year,
ls.sampling_name,
l.locality,
ST_Centroid(l.geom) as geom
FROM occurrences.observations obs,
lookup.fennoscand_species2 fennoscand,
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
"

fennoskand_obs <- read_sf(con,
                          query = fennoskand_obs_q) %>% 
  mutate(alien_kat = "Fennoskandisk forek.")


## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------

pot_alien_obs_q <- "
SELECT 
obs.species_latin_fixed,
yl.year,
ls.sampling_name,
l.locality,
ST_Centroid(l.geom) as geom
FROM occurrences.observations obs,
lookup.pot_alien_species2 pot_alien,
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
"

pot_alien_obs <- read_sf(con,
                          query = pot_alien_obs_q) %>% 
  mutate(alien_kat = "Potensielt fremmede arter")


## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------

alien_obs_q <- "
SELECT 
obs.species_latin_fixed,
yl.year,
ls.sampling_name,
l.locality,
ST_Centroid(l.geom) as geom
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
AND alien.\"riskCategory\" IN ('NK', 'SE', 'HI', 'PH', 'LO')
AND obs.identification_confidence = 'HIGH'
"  

alien_obs <- read_sf(con,
                     query = alien_obs_q) %>% 
  mutate(alien_kat = "På fremmedartslista")


## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
all_alien_obs <- fennoskand_obs %>% 
  rbind(pot_alien_obs)  %>% 
  rbind(alien_obs)



## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
kat_order <- tibble(alien_kat = c("Fennoskandisk forek.",
                                  "Potensielt fremmede arter",
                                  "På fremmedartslista"),
                    kat_order = 1:3)

all_alien_obs_agg <- all_alien_obs %>% 
  group_by(alien_kat) %>% 
  mutate(no_spec_per_kat = n_distinct(species_latin_fixed)) %>% 
  group_by(locality,
           alien_kat,
           no_spec_per_kat) %>% 
  summarise(no_spec = n_distinct(species_latin_fixed)) %>% 
  ungroup() %>% 
  mutate(alien_kat_append = factor(paste0(alien_kat, " (", no_spec_per_kat, " st.)"))) %>% 
  left_join(kat_order,
            by = c("alien_kat" = "alien_kat")) %>% 
  mutate(alien_kat_append = reorder(alien_kat_append, kat_order))


all_alien_obs_agg <- all_alien_obs_agg %>% 
  st_jitter(all_alien_obs_agg, amount = 7000) 


## ----non-native-species-kat, fig.cap = "Observasjoner av insektarter uten tidligere kjent forekomst i Norge i offentlige digitale kilder. Fennoskandisk forekomst baseres på registreringer i GBIF. Potensielt fremmede arter er basert på en sammenveing av flere datakilder, se figur 8. Observasjonene er basert på DNA-metastrekkoding og er ikke morfologisk verifisert. Arter vurdert til NR (Ikke risikovurdert) på fremmedartslista er ikke tatt med. Posisjonene er spredt ut noe på kartet for at de ikke skal skjule hverandre.  \\label{fig:non-native-species-kat}", fig.height = 8----
set.seed(12345)
norge %>% 
  filter(fylke == "Innlandet" |
         fylke =="Vestfold og Telemark" |
         fylke == "Oslo" |
         fylke == "Viken" | 
        fylke == "Trøndelag" |
          fylke == "Rogaland" |
          fylke == "Agder") %>% 
ggplot(.) +
  geom_sf() +
  geom_sf(aes(size = no_spec,
              color = alien_kat_append),
              alpha = 0.7,
          data = all_alien_obs_agg) +
  scale_color_nina(name = "Kategori") +
  scale_size(breaks = c(1, seq(from = 10,
                          to = max(all_alien_obs_agg$no_spec), 
                          by = 10),
                        max(all_alien_obs_agg$no_spec)),
             name = "Antall arter per\nkategori og lokalitet") +

  xlab("") +
  ylab("") +
  guides(color = guide_legend(override.aes = list(size=5), nrow = 2, byrow = T)) +
     theme(legend.position = "bottom",
           legend.box = "vertical",
       # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
       legend.key.size = unit(0.5, 'cm'),
       legend.text = element_text(size=7),
       legend.key.width = unit(0.5, 'cm')
     )



## ----alien_obs_agg, cache = cache_var-----------------------------------------------------------------------------------------------------------
kat_order <- tibble(alien_kat = c("Fennoskandisk forek.",
                                  "Ikke fennoskandisk forek.",
                                  "På fremmedartslista"),
                    kat_order = 1:3)

all_alien_obs_less_agg <- all_alien_obs %>% 
  st_drop_geometry() %>% 
  group_by(alien_kat) %>% 
  mutate(no_spec_per_kat = n_distinct(species_latin_fixed)) %>% 
  select(species_latin = species_latin_fixed,
         everything()) %>% 
  group_by(locality,
           alien_kat,
           sampling_name,
           species_latin) %>% 
  distinct() %>% 
  group_by(locality,
           alien_kat,
           no_spec_per_kat,
           species_latin) %>% 
  summarise(times_detected = n()) %>% 
  mutate(alien_kat_append = factor(paste0(alien_kat, " (", no_spec_per_kat, " st.)"))) %>% 
  left_join(kat_order,
            by = c("alien_kat" = "alien_kat")) %>% 
  mutate(alien_kat_append = reorder(alien_kat_append, kat_order)) %>% 
  arrange(locality,
          alien_kat,
          species_latin)

#%>% 
 # mutate(kategori_2021 = factor(kategori_2021, levels = c("DD", "NT", "VU", "EN")))

# redlisted_obs_2021 %>%
#   filter(locality == "Skog_10",
#          species_latin_fixed == "Strophosoma fulvicorne")
# 

#redlisted_obs_2021_less_agg %>% 
 # print(n = Inf)



## ----alien-detection, fig.cap = "Deteksjon av arter uten tidligere har kjent forekomst i Norge. Hver sirkel viser funn av en enkelt art i en lokalitet. Størrelsen på sirklene viser antall ganger arten ble detektert på lokalet (enkelte felleobservasjoner). Fennoskandisk forekomst baseres på registreringer i GBIF. Potensielt fremmede arter er basert på en sammenveing av flere datakilder, se figur 8. Observasjonene er basert på DNA-metastrekkoding og er ikke morfologisk verifisert. Arter vurdert til NR (Ikke risikovurdert) på fremmedartslista er ikke tatt med. \\label{fig:alien-detection}", fig.height = 7----
ggplot(all_alien_obs_less_agg,
       aes(y = locality,
           x = alien_kat)) +
  geom_point(aes(size = times_detected,
                 color = alien_kat_append),
             alpha = 0.5,
             position = position_jitter(height = 0)) +
  scale_y_discrete(limits = rev) +
  #guides(color = "none") +
  scale_color_nina(name = "Kategori") +
  scale_size(name = "Antall ganger\ndetektert",
             breaks = c(1, 2, 4, 6, 8, 10)) +
  theme(legend.position="right",
        axis.text.y = element_text(size = 6),
       axis.text.x = element_text(angle = 20,
                                  vjust = 0.6)) +
  xlab("") +
  ylab("Lokalitet")





## -----------------------------------------------------------------------------------------------------------------------------------------------
load(file = "out/region_habitat_div.rda")


## ----div-tab, results = "asis"------------------------------------------------------------------------------------------------------------------
xtable(region_habitat_div, caption = "Totalt antall identifiserte taksa fra overvåkingen i 2020-2022. \\label{tab:div-tab}")


## ----richness-locality, fig.cap = "Totalt observert artsantall per lokalitet. Merk at DNA-metastrekkodingsrutinen ble oppdatert i 2021 og at dette kan forklare der lavere artsantallet i 2020.", out.width = "\\textwidth"----
year_locality_biodiv %>% 
  mutate(habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type),
         year = as.factor(year)) %>%
  ggplot(., 
         aes(x = year,
             y = no_species)
         ) +
  geom_boxplot(aes(fill = region_name)
   #        stat = "identity"
           ) +
  scale_fill_nina(name = "År",
                  palette = "darkblue-orange") +
  ylab("Totalt antall observerte arter per lokalitet") +
    theme(axis.text.x = element_text(angle = 90, 
                                     vjust = 0.5, 
                                     hjust = 1,
                                     size = 8)) +
  facet_wrap("habitat_type",
             nrow = 2)


## ---- eval = F----------------------------------------------------------------------------------------------------------------------------------
## semi_nat_ost <- get_community_matrix(trap_type = "MF",
##                                          subset_habitat = "Semi-nat",
##                                          subset_region = "Østlandet",
##                                          transposed_matrix = T)
## 
## skog_ost <- get_community_matrix(trap_type = "MF",
##                                          subset_habitat = "Forest",
##                                          subset_region = "Østlandet",
##                                          transposed_matrix = T)
## 
## semi_nat_trond <- get_community_matrix(trap_type = "MF",
##                                          subset_habitat = "Semi-nat",
##                                          subset_region = "Trøndelag",
##                                          transposed_matrix = T)
## 
## semi_nat_sor <- get_community_matrix(trap_type = "MF",
##                                          subset_habitat = "Semi-nat",
##                                          subset_region = "Sørlandet",
##                                          transposed_matrix = T)
## 
## 
## 
## 
## ost_iNEXT_data <- Norimon:::toiNEXT(input = list(semi_nat_ost,
##                                                  skog_ost,
##                                                  semi_nat_trond,
##                                                  semi_nat_sor))
## 
## names(ost_iNEXT_data) <- c("Semi-nat_Øst",
##                            "Skog_Øst",
##                            "Semi-nat_Trond",
##                            "Semi-nat_Sør")
## 
## 
## ost_iNEXT <- iNEXT(ost_iNEXT_data,
##                    datatype = "incidence_raw",
##                    endpoint = 50)
## 
## save(ost_iNEXT,
##     file = "out/ost_iNEXT.rda")
## 


## -----------------------------------------------------------------------------------------------------------------------------------------------
load(file = "out/ost_iNEXT.rda")



## ----div-est, fig.cap = "Oppmålt og estimert total artsantall (Chao) i malaisefeller i de undersøkte regionene. Merk at ekstrapoleringen antar en perfekt deteksjonsevne.\\label{fig:div-est}", out.width = "0.9\\textwidth"----
lines <- c("Rarefaction" = "solid", 
           "Extrapolation" = "dotted")

ggiNEXT(ost_iNEXT) +
  ylab("Artsantall") +
  xlab("Antall lokaliteter") +
  theme(legend.text = element_text(size = 10),
        legend.key.width = unit(1, "cm")) +
  scale_color_manual(values = c("#004F71","#93328E", "#7A9A01", "#E57200")) +
  guides(color = guide_legend("", nrow = 2,byrow = TRUE),
         shape = guide_legend("", nrow = 2,byrow = TRUE),
         fill = guide_legend("", nrow = 2,byrow = TRUE)) +
  scale_linetype_manual(values = lines,
                      name = "Method", 
                      labels = c("Interpolert", "Ekstrapolert"))
  



## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------

spec_occ <- tbl(con,
                Id(schema = "views",
                   table = "species_occurrences")) %>% 
  filter(project_short_name == "NasIns")

spec_no_occ <- spec_occ %>% 
  filter(!grepl("MF2", sample_name)) %>% 
  group_by(species_latin_fixed) %>% 
  summarise(no_occ = as.integer(n())) %>% 
  arrange(desc(no_occ)) %>% 
  collect() %>% 
  mutate(art_no = row_number())

no_samples <- spec_occ %>%   
  filter(!grepl("MF2", sample_name)) %>% 
  summarise(no_samples = n_distinct(sample_name))



## ----spec-occ, fig.cap = "Rangering av i hvor mange felleprøver de ulike insektartene er blitt funnet. \\label{fig:spec-occ}", fig.width =  9----
spec_no_occ %>% 
ggplot(.,
       aes(x = art_no, y = no_occ)) +
  geom_bar(stat = "identity", 
           aes(fill = no_occ),
           width = 100) +
  ylab("Observert antall ganger av 1114 mulige") +
  xlab("Art") +
  scale_fill_nina(discrete = FALSE) +
  scale_x_continuous(position = "bottom",
                     breaks = c(1, 5000, 10000, 15000, 20000)) +
  theme(legend.position = "") 
#+
  #coord_flip() +
  #scale_y_reverse() +
  #theme(axis.title.x = element_text(angle = 180),
  #      axis.title.y = element_text(angle = 180))
#+
 # geom_point(aes(x = 1, y = max(no_occ)), col = ninaColors()[3])
  



## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
all_comm_mat <- get_community_matrix(trap_type = "All",
                                         dataset = "NasIns",
                                         as_tibble = T)




## -----------------------------------------------------------------------------------------------------------------------------------------------
all_comm_mat_mat <- all_comm_mat %>% 
  select(-c(year, locality))

levels <- tibble(l1 = seq(1:nrow(all_comm_mat_mat)),
                 l2 = rep(1, nrow(all_comm_mat_mat))
)


## ---- eval = F----------------------------------------------------------------------------------------------------------------------------------
## div_part <- multipart(y = all_comm_mat_mat,
##           x = levels,
##           global = T,
##           nsimul = 999,
##           scales = 0)
## 
## save(div_part,
##      file = "out/div_part.rds")


## -----------------------------------------------------------------------------------------------------------------------------------------------
load(file = "out/div_part.rds")


## ----div-part-tab, results = 'asis'-------------------------------------------------------------------------------------------------------------
div_part_tab <- tibble(Nivå = c("$\\alpha$", 
                                "$\\beta$",
                                "$\\gamma$"),
                       "Målt verdi" = div_part$oecosimu$statistic[c(1, 3, 2)],
                       "Stand. Effekt" = div_part$oecosimu$z[c(1, 3, 2)],
                       "Forventet verdi" = div_part$oecosimu$means[c(1, 3, 2)],
                       "Prob (sim)" = div_part$oecosimu$pval[c(1, 3, 2)])

print(xtable(div_part_tab, caption = "Randomiseringstest av alpha, beta og gamma-komponenter av diversitet (N = 999).  \\label{tab:div-part-tab}",
             digits = 3),
       sanitize.text.function = function(x) as.character(x),
       sanitise.rownames.function = function(x) as.character(x),
       include.rownames = F)


## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
div_betacore <- betapart.core(all_comm_mat_mat)


## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
all_beta_part_freq <- beta.multi(div_betacore) %>% 
  unlist()

all_beta_part_freq <- as.data.frame(t(round(all_beta_part_freq / all_beta_part_freq["beta.SOR"] * 100, 2)))
names(all_beta_part_freq) <- c("$\\beta_{SIM}\\%$", "$\\beta_{NES}\\%$", "$\\beta_{SOR}\\%$")

  

## ----beta-tab, results = "asis"-----------------------------------------------------------------------------------------------------------------
print(xtable(all_beta_part_freq, caption = "Komponenter i oppmålt beta-diversitet mellom lokalitetetene i overvåkingen. $\\beta_{SIM}\\%$ står for 'species turnover', $\\beta_{NES}\\%$ for 'nestedness' og $\\beta_{SOR}\\%$ for den totale $\\beta-$diversiteten.\\label{tab:beta-tab}"), sanitize.colnames.function = function(x) as.character(x))


## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
skog_ost_comm_mat <- get_community_matrix(#trap_type = "MF",
                                              dataset = "NasIns",
                                              subset_habitat = "Forest",
                                              subset_region = "Østlandet",
                                              as_tibble = T,
                                              transposed_matrix = F
                                              ) 
#%>% 
#  select(-c(year, locality))

semi_nat_comm_mat <- get_community_matrix(trap_type = "MF",
                                              dataset = "NasIns",
                                              subset_habitat = "Semi-nat",
                                              as_tibble = T,
                                              transposed_matrix = F
                                              ) 
#%>% 
#  select(-c(year, locality))





## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
skog_ost_dist_beta <- combine_dist_to_comm_mat(comm_mat = skog_ost_comm_mat,
                                               region_name = "('Østlandet')",
                                               habitat_type = "Forest")


## ---- eval = F----------------------------------------------------------------------------------------------------------------------------------
## beta_sne_lm <- glm(beta_sne~distance, data = skog_ost_dist_beta)
## summary(beta_sne_lm)
## beta_sim_lm <- glm(beta_sim~distance, data = skog_ost_dist_beta)
## summary(beta_sim_lm)
## beta_sor_lm <- glm(beta_sor~distance, data = skog_ost_dist_beta)
## summary(beta_sor_lm)
## 
## 
## with(summary(beta_sor_lm), 1 - deviance/null.deviance)


## ----skog-ost-dist-beta, fig.cap = "Oppdeling av $\\beta-$diversitet i turnover og nestedness, og dess sammenheng med avstand mellom lokaliteter. Skogslokaliteter på Østlandet, alle felletyper. \\label{fig:skog_ost_dist_beta}", fig.height = 8----
plot_beta_part(skog_ost_dist_beta)


## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
semi_nat_dist_beta <- combine_dist_to_comm_mat(comm_mat = semi_nat_comm_mat,
                                               region_name = "('Østlandet', 'Trøndelag', 'Sørlandet')",
                                               habitat_type = "Semi-nat")


## ----semi-nat-ost-dist-beta, fig.cap = " Oppdeling av $\\beta-$diversitet i turnover og nestedness, og dess sammenheng med avstand mellom lokaliteter. Semi-naturlige lokaliteter på Østlandet, i Trøndelag og på Sørlandet, kun malaisefeller. \\label{fig:semi-nat-ost-dist-beta}", fig.height = 9----
plot_beta_part(semi_nat_dist_beta)


## ---- eval = F----------------------------------------------------------------------------------------------------------------------------------
## 
## skog_ost_comm_mat_lepi <- get_community_matrix(trap_type = "MF",
##                                               dataset = "NasIns",
##                                               subset_habitat = "Forest",
##                                               subset_region = "Østlandet",
##                                               subset_orders = "Lepidoptera",
##                                               as_tibble = T,
##                                               transposed_matrix = F,
##                                               obs_at_least_x_times_somewhere = 2
##                                               )
## #%>%
## #  select(-c(year, locality))
## 
## semi_nat_ost_comm_mat_lepi <- get_community_matrix(trap_type = "MF",
##                                               dataset = "NasIns",
##                                               subset_habitat = "Semi-nat",
##                                               subset_region = "Østlandet",
##                                               subset_orders = "Lepidoptera",
##                                               as_tibble = T,
##                                               transposed_matrix = F,
##                                               obs_at_least_x_times_somewhere = 2
##                                               )
## #%>%
## #  select(-c(year, locality))
## 
## 


## ---- eval = F----------------------------------------------------------------------------------------------------------------------------------
## skog_ost_dist_beta_lepi <- combine_dist_to_comm_mat(comm_mat = skog_ost_comm_mat_lepi,
##                                                region_name = "('Østlandet')",
##                                                habitat_type = "Forest")
## 
## semi_nat_ost_dist_beta_lepi <- combine_dist_to_comm_mat(comm_mat = semi_nat_ost_comm_mat_lepi,
##                                                region_name = "('Østlandet')",
##                                                habitat_type = "Semi-nat")


## ---- eval = F----------------------------------------------------------------------------------------------------------------------------------
## skog_ost_dist_beta_lepi_mat <- skog_ost_comm_mat_lepi %>%
##   select(-c(year, locality))
## 
## levels <- tibble(l1 = seq(1:nrow(skog_ost_dist_beta_lepi_mat)),
##                  l2 = rep(1, nrow(skog_ost_dist_beta_lepi_mat))
## )


## ---- eval = F----------------------------------------------------------------------------------------------------------------------------------
## div_part_skog_lepi <- multipart(y = skog_ost_dist_beta_lepi_mat,
##           x = levels,
##           global = T,
##           nsimul = 99,
##           scales = 0)
## 
## save(div_part_skog_lepi,
##      file = "out/div_part_skog_lepi.rds")


## ---- eval = F----------------------------------------------------------------------------------------------------------------------------------
## div_betacore_skog_lepi <- betapart.core(skog_ost_dist_beta_lepi_mat)


## ---- eval = F----------------------------------------------------------------------------------------------------------------------------------
## all_beta_part_freq_skog_lepi <- beta.multi(div_betacore_skog_lepi) %>%
##   unlist()
## 
## all_beta_part_freq_skog_lepi <- as.data.frame(t(round(all_beta_part_freq_skog_lepi / all_beta_part_freq_skog_lepi["beta.SOR"] * 100, 2)))
## names(all_beta_part_freq) <- c("$\\beta_{SIM}\\%$", "$\\beta_{NES}\\%$", "$\\beta_{SOR}\\%$")
## 
## 


## ----skog-ost-dist-beta-lepi, fig.cap = "Skogslokaliteter på Østlandet. Oppdeling av $\\beta-$diversitet i turnover og nestedness i Lepidoptera, og dess sammenheng med avstand mellom lokaliteter. Alle arter som er blitt funnet mindre en 5 ganger i alle lokaliteter exkludert. \\label{fig:skog-ost-dist-beta-lepi}", fig.height = 9, eval = F----
## plot_beta_part(skog_ost_dist_beta_lepi)


## ----semi-nat-ost-dist-beta-lepi, fig.cap = "Semi-naturlige lokaliteter på Østlandet. Oppdeling av $\\beta-$diversitet i turnover og nestedness hos Lepidoptera, og dess sammenheng med avstand mellom lokaliteter. Alle arter som er blitt funnet mindre en 5 ganger i alle lokaliteter exkludert. \\label{fig:semi_nat_ost_dist_beta_lepi}", fig.height = 9, eval = F----
## plot_beta_part(semi_nat_ost_dist_beta_lepi)


## -----------------------------------------------------------------------------------------------------------------------------------------------

get_asv_loc <- function(species = NULL,
                        dataset = c("NasIns"),
                        subset_year = NULL
                        ){
  
  asv_perc_reads <- tbl(con,
                        Id(schema = "views",
                           table = "asv_perc_reads"))
  
  asv_wider <- asv_perc_reads %>% 
    filter(species_latin_fixed %in% species) %>% 
   # mutate(value = round(perc_reads * 100, 4)) %>% 
    collect() %>% 
    as.data.frame()
  #%>% 
    # pivot_wider(names_from = "sequence_id",
    #             values_from = "perc_reads",
    #             values_fill = 0) 
    #   
  return(asv_wider)          

}




## ---- eval = F----------------------------------------------------------------------------------------------------------------------------------
## 
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
## most_abund_dna <- dbGetQuery(con,
##                           most_abund_dna_q)
## 
## 


## ---- eval = F----------------------------------------------------------------------------------------------------------------------------------
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
## most_common <- dbGetQuery(con,
##                           most_common_q)
## 
## 
## 


## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------

plot_asv <- function(species = NULL,
                     background = "ortofoto/norge2.png",
                     pie_scale = 0.5,
                     size = 0.1,
                     caption = TRUE,
                     title = NULL,
                     scale_to_sum_reads = TRUE,
                     ...){

  jon_asv <- get_asv_loc(species = species) %>% 
    mutate(scale_sum_reads = scale(sum_reads, center = min(sum_reads), scale = diff(range(sum_reads)))) 
  # %>% 
    #  filter(locality %in% c("Skog_02", "Semi-nat_11"))

  if(!file.exists(background)) stop("Background image file not found")
  
  tt <- terra::rast(background)
  ext_background <- as.vector(terra::ext(tt))
  img <- png::readPNG(background)
  g <- grid::rasterGrob(img, interpolate=TRUE)


  if(scale_to_sum_reads){
        jon_asv$r <- log(jon_asv$sum_reads) * diff(ext_background[1:2]) / 100 * pie_scale
  } else {
    jon_asv$r <- diff(ext_background[1:2]) / 100 * pie_scale 
  }
  
  jon_asv <- jon_asv %>% 
    mutate(locality = factor(locality, levels = unique(locality[order(desc(r))]))) %>% 
    arrange(desc(r)) 
  
p1 <- ggplot() +
    annotation_custom(g, xmin = ext_background[1], 
                      xmax = ext_background[2], 
                      ymin = ext_background[3],
                      ymax = ext_background[4]) + 
    lapply(split(jon_asv, jon_asv$locality), 
         function(d) {
       geom_arc_bar(aes(x0 = x_25833, 
                         y0 = y_25833,
                         r = r,
                         r0 = 0,
                         amount = perc_reads,
                         fill = sequence_id
                         ), 
                     size = size,
                     data = d, 
                     stat = "pie",
                     inherit.aes = TRUE)
         })  +
  coord_fixed() +
    xlim(ext_background[1:2]) +
    ylim(ext_background[3:4]) +
     theme(legend.position="none",
           panel.border = element_blank(),
           #axis.line=element_line(color = "black"),
           axis.text.x=element_blank(), #remove x axis labels
           axis.ticks.x=element_blank(), #remove x axis ticks
           axis.text.y=element_blank(),  #remove y axis labels
           axis.ticks.y=element_blank(),
           plot.margin = unit(c(0, 0, 0, 0), "pt"),
           plot.caption.position = "panel",
           plot.caption = element_text(vjust = 10,
                                       hjust = 0.5,
                                       size = 14,
                                       margin = margin(0, 0, 0, 0)),
           plot.title = element_text(hjust = 0.5,
                                     vjust = 0)

           ) +
      xlab("") +
      ylab("") +
  scale_fill_nina() 

  if(caption){
    p1 <- p1 +
     labs(caption = species)
  }

  if(!is.null(title)){
    p1 <- p1 +
      labs(title = title)
  }

  p1
  
  }






## ----haplotype-kart, fig.cap = "Frekvenser av haplotyper (varianter av ASV/DNA-kopier) for ulike lokaliteter som viser genetisk variasjon innen hver enkeltart. Her er de tre artene som er fanget i flest tilfeller innenfor hver av ordenene sommerfugler (Lepidoptera), biller (Coleoptera), tovinger (Diptera) og vepser (Hymenoptera) valgt ut som eksempler. Størrelsen på sirklene er logaritmisk skalert etter mengden DNA av arten i hver lokalitet. \\label{fig:haplotype-kart}", eval = T, out.width = "\\textwidth"----
include_graphics("figures/asv_composite.png")


## ---- eval = F----------------------------------------------------------------------------------------------------------------------------------
## most_common_poll_q <- "
## 
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
## most_common_poll <- dbGetQuery(con,
##                                most_common_poll_q)






## ----poll-haplotype-kart, fig.cap = "Frekvenser av haplotyper (varianter av ASV/DNA-kopier) for et utvalg av pollinatorer-familier, langtungebier (Apidae), nymfevinger (Nymphalidae), glansvinger (Lycaenidae), og blomsterfluer (Syrphidae). Størrelsen på sirklene er logaritmisk skalert etter mengden DNA av arten i hver lokalitet. \\label{fig:poll-haplotype-kart}", eval = T, out.height = "\\textwidth"----
include_graphics("figures/asv_poll_composite.png")


## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
locality_sampling_biomass_temp <- biomass_mf_locality_sampling %>% 
  left_join(locality_sampling_loggerdata,
            by = c("sampling_name" = "sampling_name",
                   "year" = "year",
                   "locality" = "locality",
                   "habitat_type" = "habitat_type",
                   "region_name" = "region_name")) %>% 
  select(Year = year,
         everything())

locality_sampling_biomass_temp_compl_light <- locality_sampling_biomass_temp %>% 
  filter(sum_wet_weight > 0,
         !is.na(avg_values_MX2301A_temperature),
         !is.na(avg_values_MX2301A_rh),
         !is.na(avg_values_MX2202_light),
         !is.na(stddev_values_MX2301A_temperature),
         !is.na(stddev_values_MX2301A_rh),
         !is.na(stddev_values_MX2202_light))




## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------

sampling_week_q <- "
SELECT sampling_name, date_part('week', end_date) as ukenummer
FROM events.locality_sampling ls
"

sampling_week <- dbGetQuery(con,
                            sampling_week_q) %>% 
  select(Ukenummer = ukenummer,
         everything())


## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
locality_sampling_biomass_temp_compl_light_week <- locality_sampling_biomass_temp_compl_light %>% 
  left_join(sampling_week,
            by = c("sampling_name" = "sampling_name"))


## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
global_model <- lmer(sum_wet_weight/no_trap_days ~ 
                       scale(avg_values_MX2301A_temperature) +
                       scale(avg_values_MX2301A_rh) + 
                       scale(avg_values_MX2202_light) +
                       scale(stddev_values_MX2301A_temperature) +
                       scale(stddev_values_MX2301A_rh) + 
                       scale(stddev_values_MX2202_light) +
                       habitat_type + 
                       Year +
                       (1|Year:Ukenummer) +
                       (1|locality), 
                     data = locality_sampling_biomass_temp_compl_light_week,
                     na.action = "na.fail",
                     REML = F)

cand_models <- dredge(global_model,
                      fixed = c("habitat_type",
                                "Year"))


## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
model_averages <- model.avg(cand_models, subset = delta <= 4)
#summary(model_averages)


## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
top_models <- get.models(cand_models, subset = delta <= 4)
biomass_temp_AIC <- AICcmodavg::aictab(top_models) %>% 
  as_tibble() %>% 
  # mutate(Modell = c("Full modell", 
  #                   "-avg(RH)")
  #        ) %>% 
  select(Modnames,
         `Antall par.` = K,
         AICc,
         `Delta AICc` = Delta_AICc,
         `AIC-vekt` = AICcWt,
         `Kum.vekt` = Cum.Wt)



## ----temp_AICtab, results = "asis"--------------------------------------------------------------------------------------------------------------
print(xtable(biomass_temp_AIC, caption = "AICc-vekting av de beste modellene av biomasse fra loggerdata. Kun de modeller med max 4 delta-AIC fra den beste modellen er tatt med. \\label{tab:temp_AICtab}" ))



## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
par_weights <- MuMIn::sw(get.models(cand_models, subset = delta <= 4)) 

par_weights <- tibble("Parameter" = attr(par_weights, "names"),
                      "Total vekt" = as.numeric(par_weights),
                      "Med i antall modeller" = round(attr(par_weights, "n.models"))
  )

par_vals_full <- round(model_averages$coefficients[1,], 4) %>% 
  data.frame() 

par_vals_full$Parameter <- rownames(par_vals_full) 
colnames(par_vals_full) <- c("Parameterestimat", "Parameter")



par_weights <- par_weights %>% 
  full_join(par_vals_full,
            by = c("Parameter" = "Parameter")) %>% 
  filter(Parameter != 'habitat_type',
         Parameter != 'habitat_typeSemi-nat',
         Parameter != 'offset(log(no_trap_days))',
         Parameter != '(Intercept)'
         ) %>% 
  mutate(Parameter = ifelse(Parameter == "scale(avg_values_MX2202_light)", "Lys", Parameter),
         Parameter = ifelse(Parameter == "scale(avg_values_MX2301A_temperature)", "Temperatur", Parameter),
         Parameter = ifelse(Parameter == "scale(stddev_values_MX2202_light)", "Variasjon-Lys", Parameter),
         Parameter = ifelse(Parameter == "scale(stddev_values_MX2301A_rh)", "Variasjon-Luftfuktighet", Parameter),
         Parameter = ifelse(Parameter == "scale(avg_values_MX2301A_rh)", "Luftfuktighet", Parameter),
         Parameter = ifelse(Parameter == "scale(stddev_values_MX2301A_temperature)", "Variasjon-Temperatur", Parameter)
         ) %>% 
  arrange(desc(`Total vekt`), Parameter)
  



## ----logger-model-par, results = "asis"---------------------------------------------------------------------------------------------------------
print(xtable(par_weights, digits = 3, caption = "Vekter av parametere i de beste modellene (delta-AICc < 4) av innsamlet biomasse og deres vektede parameterestimat.\\label{tab:logger-model-par}"), size = "footnotesize")


## -----------------------------------------------------------------------------------------------------------------------------------------------
#sw(model_averages)

r2_verdi_biomasse_temp_glob <- r.squaredGLMM(global_model)


## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
nas_ins_lok_1000m_skog_q <- "
SELECT skog.*
FROM backgrounds.locality_1000m_buffer_skogtype_agg skog,
events.year_locality yl
WHERE skog.locality_id = yl.locality_id
AND yl.project_short_name = 'NasIns'

"


nas_ins_lok_1000m_skog <- read_sf(con,
                                query = nas_ins_lok_1000m_skog_q)



## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------

locality_skog_area <- nas_ins_lok_1000m_skog %>% 
  mutate(area = st_area(geom)) %>% 
  st_drop_geometry() %>% 
  group_by(locality,
           treslag,
           bonitet) %>% 
  summarize(tot_area = sum(area)) %>% 
  pivot_wider(id_cols = locality,
              names_from = c(treslag, bonitet), 
              values_from = tot_area
              ) %>% 
    mutate_if(is.numeric , units::drop_units) %>% 
  mutate_if(is.numeric , replace_na, replace = 0)



## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
locality_skog_area_data <- locality_skog_area %>% 
  ungroup() %>% 
  select(-"locality")

skog_rda <- decorana(locality_skog_area_data)

#summary(skog_rda) #hm, too large axis length?


## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
skog_PCA <- rda(locality_skog_area_data, scale = F)

#summary(skog_PCA)


## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------

skog_scores <- scores(skog_PCA)$sites %>% 
  as_tibble() %>% 
  cbind(locality_skog_area["locality"]) %>% 
  select(locality,
         skog_PC1 = PC1,
         skog_PC2 = PC2)


locality_skog_area <- locality_skog_area %>% 
  left_join(skog_scores,
            by = c("locality" = "locality"))


## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
#decorana(locality_ar5_area_data) #check axis length

#ar5.var <- apply (locality_ar5_area_data, 2, var)
#ar5.var <- ar5.var/sum(ar5.var)

#sort(ar5.var, decreasing = T)

ar5_PCA_norm <- rda(locality_ar5_area_data, scale = T)
#summary(ar5_PCA_norm)

ar5_PCA <- rda(locality_ar5_area_data, scale = F)
#summary(ar5_PCA)

ar5_scores <- scores(ar5_PCA)$sites %>% 
  as_tibble() %>% 
  cbind(locality_ar5_area["locality"]) %>% 
  select(locality,
         ar5_PC1 = PC1,
         ar5_PC2 = PC2)


## -----------------------------------------------------------------------------------------------------------------------------------------------
hm <- autoplot(ar5_PCA, legend.position = "none") +
  scale_color_nina(palette = "orange-green") 


## ----ar5-PCA, fig.pos = "!bh",  fig.cap = "PCA-analyse over landskapskomposisjon i arealtype etter AR5. Verdiene er ikke normalisert, da alle er fra samme datakilde og samme skala (km2).\\label{fig:ar5-PCA}", out.width = "0.8\\textwidth"----
hm


## -----------------------------------------------------------------------------------------------------------------------------------------------
tt <- autoplot(skog_PCA,
               legend.position = "none") +
  scale_color_nina(palette = "orange-green") 


## ----skog-PCA, fig.pos = "!hb", fig.cap = "PCA-analyse over landskapskomposisjon i skogtype og bonitet etter AR5 (basert på kartdata fra SR16). Verdiene er ikke normalisert, da alle er fra samme datakilde og samme skala (km2).", fig.dim = c(10, 8), out.width = "0.7\\textwidth"----
tt


## ----ar5-shannon-plot, fig.pos = "!h", fig.cap = "Fordeling av Shannon index av AR5:arealtype innen en buffer med radius 1000m kring lokalitetene. \\label{fig:ar5-shannon-plot}", out.width = "\\textwidth"----
loc_reg_year <- localities %>% 
  st_drop_geometry() %>% 
    select(year, locality, region_name, habitat_type)


locality_ar5_shdi %>% 
  left_join(loc_reg_year,
            by = c("locality" = "locality")
  ) %>% 
  mutate(habitat_type = ifelse(grepl("Semi-nat", locality), "Semi-nat", "Skog")) %>% 
  mutate(habitat_type = factor(habitat_type,levels = c("Semi-nat", "Skog"))) %>% 
    group_by(region_name) %>% 
  ggplot(.) +
  geom_point(aes(x = year,
                 y = lsm_l_shdi,
                 color = region_name,
                 pch = habitat_type),
             size = 3,
             alpha = 0.7,
             position = position_jitterdodge(jitter.height = 0.1,
                                             jitter.width = 0,
                                             dodge.width = 0.3)) +
  scale_color_nina(name = "Region") +
  scale_x_continuous(breaks = seq(2020, 2022)) +
  guides(color = guide_legend(override.aes = list(size = 3)),
         pch = guide_legend(title = "Habitattype",
                            override.aes = list(size = 3))) +
  xlab("År") +
  ylab("Shannon index")




## ----ar5-perimeter-plot, fig.pos = "!h", fig.cap = "Fordeling av total lengde av grenser mellom arealtyper, innen et buffertområde på  1 km kring lokalitetene.", out.width = "\\textwidth"----


locality_ar5_tot_perimeter %>% 
  left_join(loc_reg_year,
            by = c("locality" = "locality")
  ) %>% 
  mutate(habitat_type = ifelse(grepl("Semi-nat", locality), "Semi-nat", "Skog")) %>% 
  mutate(habitat_type = factor(habitat_type,levels = c("Semi-nat", "Skog"))) %>% 
    group_by(region_name) %>% 
  ggplot(.) +
  geom_point(aes(x = year,
                 y = tot_perimeter_length,
                 color = region_name,
                 pch = habitat_type),
             size = 3,
             alpha = 0.7,
             position = position_jitterdodge(jitter.height = 0.1,
                                             jitter.width = 0,
                                             dodge.width = 0.3)) +
  scale_color_nina(name = "Region") +
  scale_x_continuous(breaks = seq(2020, 2022)) +
  guides(color = guide_legend(override.aes = list(size = 3)),
         pch = guide_legend(title = "Habitattype",
                            override.aes = list(size = 3))) +
  xlab("År") +
  ylab("Tot. lengde linjære element")




## -----------------------------------------------------------------------------------------------------------------------------------------------
locality_ar5_area <- locality_ar5_area %>% 
  left_join(ar5_scores,
            by = c("locality" = "locality"))



## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
locality_sampling_biomass_big_join <- locality_sampling_biomass_temp_compl_light_week %>% 
  left_join(locality_ar5_area,
            by = c("locality" = "locality")) %>% 
  left_join(locality_ar5_shdi,
            by = c("locality" = "locality")) %>% 
  left_join(locality_ar5_tot_perimeter,
            by = c("locality" = "locality")) %>% 
  left_join(locality_skog_area,
            by = c("locality" = "locality")) %>% 
  left_join(ano_herb_agg,
            by = c("locality" = "locality",
                   "Year" = "year"),
            copy = T) %>% 
  left_join(tree_agg,
            by = c("locality" = "locality",
                   "Year" = "year"),
            copy = T)

write_csv(locality_sampling_biomass_big_join,
          file = "out/locality_sampling_biomass_big_join.csv")



## -----------------------------------------------------------------------------------------------------------------------------------------------
locality_sampling_biomass_big_join <- read_csv(file = "out/locality_sampling_biomass_big_join.csv")

locality_sampling_biomass_big_join_forest <- locality_sampling_biomass_big_join %>% 
  filter(habitat_type == "Forest")

locality_sampling_biomass_big_join_semi_nat <- locality_sampling_biomass_big_join %>% 
  filter(habitat_type == "Semi-nat")



## ----perim-shann-plot, fig.pos = "!h", fig.cap = "Sammenheng mellom diversiteten av arealtyper (Shannon-indeks) og mengden linjære element i landskapene rundt skogslokalitetene. \\label{fig:perim-shann-plot}", out.width = "0.7\\textwidth"----
locality_sampling_biomass_big_join_forest %>% 
  ggplot(.) +
  geom_point(aes(x = lsm_l_shdi,
                 y = tot_perimeter_length),
                 color = ninaColors()[1],
             size = 3) +
  ylab("Total lengde linjære element") +
  xlab("Artype Shannon indeks") 


## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
global_model_big_forest <- lmer(sum_wet_weight/no_trap_days ~ scale(avg_values_MX2301A_temperature) +
                       # scale(avg_values_MX2301A_rh) +
                        scale(avg_values_MX2202_light) +
                        scale(stddev_values_MX2301A_temperature) +
                        scale(stddev_values_MX2301A_rh) +
                        #scale(stddev_values_MX2202_light) +
                        scale(ar5_PC1) +
                        scale(ar5_PC2) +
                        scale(lsm_l_shdi) +
                        scale(tot_perimeter_length) +
                        scale(skog_PC1) +
                        scale(skog_PC2) +
                        scale(ano_mean_no_spec) +
                        scale(ano_mean_cover) +
                        scale(avg_age) +
                        dom_tree  + 
                        Year +
                        (1|locality) +
                        (1 | Year:Ukenummer)
                       ,
                      data = locality_sampling_biomass_big_join_forest,
                      na.action = "na.fail",
                      REML = F,
                      control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))) #Suppress warnings due to near zero estimate for locality random effect

# global_model_big_forest <- glm(sum_wet_weight/no_trap_days ~ scale(avg_values_MX2301A_temperature) + 
#                        # scale(avg_values_MX2301A_rh) + 
#                         scale(avg_values_MX2202_light) +
#                         scale(stddev_values_MX2301A_temperature) +
#                         scale(stddev_values_MX2301A_rh) + 
#                         scale(stddev_values_MX2202_light) +
#                         scale(ar5_PC1) +
#                         scale(ar5_PC2) +
#                         scale(lsm_l_shdi) +
#                         scale(tot_perimeter_length) +
#                         scale(skog_PC1) +
#                         scale(skog_PC2) +
#                         scale(ano_mean_no_spec) +
#                         scale(ano_mean_cover) +
#                         scale(avg_age) +
#                         dom_tree +
#                         Year, 
#                         data = locality_sampling_biomass_big_join_forest,
#                         na.action = "na.fail")




cand_models_big_forest <- dredge(global_model_big_forest,
                      fixed = c("scale(avg_values_MX2301A_temperature)",
                                #"scale(avg_values_MX2301A_rh)",
                                "scale(avg_values_MX2202_light)",
                                "scale(stddev_values_MX2301A_temperature)",
                                "scale(stddev_values_MX2301A_rh)"
                                #,"scale(stddev_values_MX2202_light)"
                                )
                      )



## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
model_averages_big_forest <- model.avg(cand_models_big_forest, subset = delta <= 2)
#summary(model_averages_big_forest)


## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
top_models_big_forest <- get.models(cand_models_big_forest, subset = delta <= 2)

biomass_temp_AIC_big_forest <- AICcmodavg::aictab(top_models_big_forest) %>% 
  as_tibble() %>% 
#  mutate(Modell = c("Full modell", 
 #                   "-AR5_PC2, -AR5_Skog_PC1, -AR5_Skog_PC2")
  #       ) %>% 
  select(Modnames,
         `Antall par.` = K,
         AICc,
         `Delta AICc` = Delta_AICc,
         `AIC-vekt` = AICcWt,
         `Kum.AIC-vekt` = Cum.Wt)



## ----temp-AICtab-big-forest, results = "asis"---------------------------------------------------------------------------------------------------

print(xtable(biomass_temp_AIC_big_forest, caption = "AICc-vekting av de beste modellene av biomasse fra skogslokaliteter. Kun de modeller med max 4 delta-AIC fra den beste modellen er tatt med. \\label{tab:temp-AICtab-big-forest}" ))



## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
par_weights_forest <- MuMIn::sw(get.models(cand_models_big_forest, subset = delta <= 4)) 

par_weights_forest <- tibble("Parameter" = attr(par_weights_forest, "names"),
                      "Total vekt" = as.numeric(par_weights_forest),
                      "Med i antall modeller" = round(attr(par_weights_forest, "n.models"))
)

par_vals_full <- model_averages_big_forest$coefficients[1,] %>% 
  data.frame() 

par_vals_full$Parameter <- rownames(par_vals_full) 
colnames(par_vals_full) <- c("Parameterestimat", "Parameter")



par_weights_forest <- par_weights_forest %>% 
  left_join(par_vals_full,
            by = c("Parameter" = "Parameter")) %>% 
  filter(Parameter != 'habitat_type',
         Parameter != 'offset(log(no_trap_days))',
         Parameter != 'scale(avg_values_MX2202_light)',
         Parameter != 'scale(avg_values_MX2301A_temperature)',
         Parameter != 'scale(stddev_values_MX2202_light)',
         Parameter != 'scale(stddev_values_MX2301A_rh)',
         Parameter != 'scale(stddev_values_MX2301A_temperature)',
         ) %>% 
  mutate(Parameter = ifelse(Parameter == "dom_tree", "Landsskog: dominerende treslag", Parameter),
         Parameter = ifelse(Parameter == "scale(avg_age)", "Landsskog: skogalder", Parameter),
         Parameter = ifelse(Parameter == "scale(ano_mean_cover)", "ANO: dekningsgrad karpl.", Parameter),
         Parameter = ifelse(Parameter == "scale(ano_mean_no_spec)", "ANO: artsantall karpl.", Parameter),
         Parameter = ifelse(Parameter == "scale(ar5_PC1)", "AR5-arealtype: PCA-akse 1", Parameter),
         Parameter = ifelse(Parameter == "scale(ar5_PC2)", "AR5-arealtype: PCA-akse 2", Parameter),
         Parameter = ifelse(Parameter == "scale(skog_PC1)", "AR5-skog: PCA-akse 1", Parameter),
         Parameter = ifelse(Parameter == "scale(skog_PC2)", "AR5-skog: PCA-akse 2", Parameter),
         Parameter = ifelse(Parameter == "scale(lsm_l_shdi)", "AR5-arealtype: Shannon indeks", Parameter),
         Parameter = ifelse(Parameter == "scale(tot_perimeter_length)", "AR5-arealtype: linjære elementer", Parameter),
         ) %>% 
  arrange(desc(`Total vekt`), Parameter)
  





## ----par_weights_forest, results = "asis"-------------------------------------------------------------------------------------------------------
par_weights_forest <- par_weights_forest %>% 
  mutate(Parameter = ifelse(Parameter == "Year", "År", Parameter),
         Parameter = ifelse(Parameter == "AR5-skog: PCA-akse 2", "AR5-skog: PCA-akse 2 (Bonitet)", Parameter),
         Parameter = ifelse(Parameter == "AR5-skog: PCA-akse 1", "AR5-skog: PCA-akse 1 (Løv/Barr)", Parameter),
         Parameter = ifelse(Parameter == "AR5-arealtype: PCA-akse 1", "AR5-arealtype: PCA-akse 1 (ikke skog/skog)", Parameter),
         Parameter = ifelse(Parameter == "AR5-arealtype: PCA-akse 2", "AR5-arealtype: PCA-akse 2 (åpen mark/fulldyrka mark)", Parameter))

print(xtable(par_weights_forest, caption = "Vekter av miljøparametere i de beste modellene av innsamlet biomasse i skog (delta-AICc <= 4). \\label{tab:par_weights_forest}."), size = "footnotesize")


## -----------------------------------------------------------------------------------------------------------------------------------------------
r2_verdi_biomasse_big_forest <- r.squaredGLMM(get.models(cand_models_big_forest, subset = delta < 0.01)[[1]])


## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
global_model_big_semi_nat <- lmer(sum_wet_weight/no_trap_days ~ scale(avg_values_MX2301A_temperature) +
                        scale(avg_values_MX2301A_rh) +
                        scale(avg_values_MX2202_light) +
                        scale(stddev_values_MX2301A_temperature) +
                        scale(stddev_values_MX2301A_rh) +
                        scale(stddev_values_MX2202_light) +
                        scale(ar5_PC1) +
                        scale(ar5_PC2) +
                        scale(lsm_l_shdi) +
                        scale(tot_perimeter_length) +
                        scale(skog_PC1) +
                        scale(skog_PC2) +
                        scale(ano_mean_no_spec) +
                        scale(ano_mean_cover) +
                         Year +
                        (1|locality) +
                         (1 | Year:Ukenummer),
                      data = locality_sampling_biomass_big_join_semi_nat,
                      na.action = "na.fail",
                      REML = F,
                      control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))) #Suppress warnings due to near zero estimate for locality random effect

# global_model_big_semi_nat <- glm(sum_wet_weight/no_trap_days ~ scale(avg_values_MX2301A_temperature) + 
#                        # scale(avg_values_MX2301A_rh) + 
#                         scale(avg_values_MX2202_light) +
#                         scale(stddev_values_MX2301A_temperature) +
#                         scale(stddev_values_MX2301A_rh) + 
#                         scale(stddev_values_MX2202_light) +
#                         scale(ar5_PC1) +
#                         scale(ar5_PC2) +
#                         scale(lsm_l_shdi) +
#                         scale(tot_perimeter_length) +
#                         scale(skog_PC1) +
#                         scale(skog_PC2) +
#                         scale(ano_mean_no_spec) +
#                         scale(ano_mean_cover), 
#                       data = locality_sampling_biomass_big_join_semi_nat,
#                       na.action = "na.fail"
# )
                     


cand_models_big_semi_nat <- dredge(global_model_big_semi_nat,
                      fixed = c("scale(avg_values_MX2301A_temperature)",
                                "scale(avg_values_MX2301A_rh)",
                                "scale(avg_values_MX2202_light)",
                                "scale(stddev_values_MX2301A_temperature)",
                                "scale(stddev_values_MX2301A_rh)",
                                "scale(stddev_values_MX2202_light)"
                                )
                      )



## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
model_averages_big_semi_nat <- model.avg(cand_models_big_semi_nat, subset = delta <= 2)
#summary(model_averages_big_semi_nat)


## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
top_models_big_semi_nat <- get.models(cand_models_big_semi_nat, subset = delta <= 2)

biomass_temp_AIC_big_semi_nat <- AICcmodavg::aictab(top_models_big_semi_nat) %>% 
  as_tibble() %>% 
  # mutate(Modell = c("Full modell", 
  #                   "-AR5_PC2, -AR5_Skog_PC1, -AR5_Skog_PC2")
  #        ) %>% 
  select(Modnames,
         `Antall par.` = K,
         AICc,
         `Delta AICc` = Delta_AICc,
         `AIC-vekt` = AICcWt,
         `Kum.AIC-vekt` = Cum.Wt)



## ----temp-AICtab-big-semi-nat, results = "asis"-------------------------------------------------------------------------------------------------

print(xtable(biomass_temp_AIC_big_semi_nat, caption = "AICc-vekting av de beste modellene av biomasse fra semi-naturelle lokaliteter. Kun de modeller med max 2 delta-AIC fra den beste modellen er tatt med. \\label{tab:temp-AICtab-big-semi-nat}" ))



## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
par_weights_semi_nat <- MuMIn::sw(get.models(cand_models_big_semi_nat, subset = delta <= 4)) 

par_weights_semi_nat <- tibble("Parameter" = attr(par_weights_semi_nat, "names"),
                      "Total vekt" = as.numeric(par_weights_semi_nat),
                      "Med i antall modeller" = attr(par_weights_semi_nat, "n.models")
)

par_vals_full <- model_averages_big_semi_nat$coefficients[1,] %>% 
  data.frame() 

par_vals_full$Parameter <- rownames(par_vals_full) 
colnames(par_vals_full) <- c("Parameterestimat", "Parameter")



par_weights_semi_nat <- par_weights_semi_nat %>% 
  left_join(par_vals_full,
            by = c("Parameter" = "Parameter")) %>% 
  filter(Parameter != 'habitat_type',
         Parameter != 'offset(log(no_trap_days))',
         Parameter != 'scale(avg_values_MX2202_light)',
         Parameter != 'scale(avg_values_MX2301A_temperature)',
         Parameter != 'scale(stddev_values_MX2202_light)',
         Parameter != 'scale(stddev_values_MX2301A_rh)',
         Parameter != 'scale(stddev_values_MX2301A_temperature)',
         Parameter != 'scale(avg_values_MX2301A_rh)'
         ) %>% 
  mutate(Parameter = ifelse(Parameter == "dom_tree", "Landsskog: dominerende treslag", Parameter),
         Parameter = ifelse(Parameter == "scale(avg_age)", "Landsskog: skogalder", Parameter),
         Parameter = ifelse(Parameter == "scale(ano_mean_cover)", "ANO: dekningsgrad kartpl.", Parameter),
         Parameter = ifelse(Parameter == "scale(ano_mean_no_spec)", "ANO: artsantall kartpl.", Parameter),
         Parameter = ifelse(Parameter == "scale(ar5_PC1)", "AR5-arealtype: PCA-akse 1", Parameter),
         Parameter = ifelse(Parameter == "scale(ar5_PC2)", "AR5-arealtype: PCA-akse 2", Parameter),
         Parameter = ifelse(Parameter == "scale(skog_PC1)", "AR5-skog: PCA-akse 1", Parameter),
         Parameter = ifelse(Parameter == "scale(skog_PC2)", "AR5-skog: PCA-akse 2", Parameter),
         Parameter = ifelse(Parameter == "scale(lsm_l_shdi)", "AR5-arealtype: Shannon indeks", Parameter),
         Parameter = ifelse(Parameter == "scale(tot_perimeter_length)", "AR5-arealtype: linjære elementer", Parameter),
         ) %>% 
  mutate(Parameter = ifelse(Parameter == "Year", "År", Parameter),
         Parameter = ifelse(Parameter == "AR5-skog: PCA-akse 2", "AR5-skog: PCA-akse 2 (Bonitet)", Parameter),
         Parameter = ifelse(Parameter == "AR5-skog: PCA-akse 1", "AR5-skog: PCA-akse 1 (Løv/Barr)", Parameter),
         Parameter = ifelse(Parameter == "AR5-arealtype: PCA-akse 1", "AR5-arealtype: PCA-akse 1 (ikke skog/skog)", Parameter),
         Parameter = ifelse(Parameter == "AR5-arealtype: PCA-akse 2", "AR5-arealtype: PCA-akse 2 (åpen mark/fulldyrka mark)", Parameter)) %>%   arrange(desc(`Total vekt`), Parameter)
  




## ----results = "asis"---------------------------------------------------------------------------------------------------------------------------
print(xtable(par_weights_semi_nat, caption = "Vekter av miljøparametere i de beste modellene av innsamlet biomasse i semi-naturelle lokaliteter (delta-AICc <= 4).\\label{tab:par-weights-semi-nat}"), size = "footnotesize")


## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
r2_verdi_biomasse_big_semi_nat <- r.squaredGLMM(get.models(cand_models_big_semi_nat, subset = delta < 0.1)[[1]])


## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
#Here summing all the no_reads per species_latin_fixed and sample

otu_loc_q <- "
SELECT l.locality, l.region_name, l.habitat_type, ls.sampling_name, date_part('week', ls.end_date) as end_week, st.sample_name, o.species_latin_fixed as otu, sum(o.no_reads)::integer no_reads
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
AND yl.year = '2022'
AND yl.project_short_name = 'NasIns'
AND traps.trap_short_name = 'MF1'
AND (o.taxonomic_level = 'species'
OR o.taxonomic_level = 'genus'
OR o.taxonomic_level = 'family'
OR o.taxonomic_level = 'order'
)
GROUP BY l.locality, l.region_name, l.habitat_type, ls.end_date, ls.sampling_name, st.sample_name, o.species_latin_fixed

" 



otu_loc_raw <- dbGetQuery(con,
                      otu_loc_q) %>% 
  as_tibble %>% 
  mutate(habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type))

otu_loc_wider <- otu_loc_raw %>% 
  pivot_wider(names_from = otu,
              values_from = no_reads,
              values_fill = 0,) %>% 
  as.data.frame()

rownames(otu_loc_wider) <- otu_loc_wider$sampling_name
otu_loc <- otu_loc_wider[-c(1:6)]
#dim(otu_loc)



## ---- message = FALSE, results = "hide", cache = cache_var--------------------------------------------------------------------------------------
ord_malaise <- vegan::metaMDS(otu_loc, 
                              distance = "bray",
                              trymax = 1000)
  
ord_malaise_aug <- cbind(otu_loc_wider[c(1:6)], vegan::scores(ord_malaise)[1])



## ---- cache = cache_var-------------------------------------------------------------------------------------------------------------------------
fig_habitats <- 
  ggplot(ord_malaise_aug,
                       aes(x = sites.NMDS1,
                           y = sites.NMDS2,
                           shape = habitat_type,
                           color = end_week)) +
  geom_point() +
  stat_ellipse(aes(x = sites.NMDS1,
                   y = sites.NMDS2,
                   group = habitat_type),
                   show.legend = F,
               data = ord_malaise_aug,
               inherit.aes = F) +
  scale_color_nina(name = "Tømmeuke",
  discrete = F, palette = "blue-orange") +
    scale_shape_discrete(name = "Habitattype") +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black")) 
  



## ---- message = FALSE, results = "hide", cache = cache_var--------------------------------------------------------------------------------------

otu_loc_wider_semi_nat <- otu_loc_raw %>% 
  filter(habitat_type == "Semi-nat") %>% 
  pivot_wider(names_from = otu,
              values_from = no_reads,
              values_fill = 0,) %>% 
  as.data.frame()

rownames(otu_loc_wider_semi_nat) <- otu_loc_wider_semi_nat$sampling_name
otu_loc_semi_nat <- otu_loc_wider_semi_nat[-c(1:6)]


ord_semi_nat <- vegan::metaMDS(otu_loc_semi_nat, 
                              distance = "bray",
                              trymax = 100)
  
ord_semi_nat_aug <- cbind(otu_loc_wider_semi_nat[c(1:6)], vegan::scores(ord_semi_nat)[1])
  
#points <- data.frame(ord_semi_nat$points)

points <- ord_semi_nat_aug %>% 
  left_join(locality_sampling_biomass_big_join_semi_nat,
             by = c("sampling_name" = "sampling_name",
                    "locality" = "locality",
                    "region_name" = "region_name",
                    "habitat_type" = "habitat_type"))  %>% 
  mutate(region_habitat = paste0(region_name, " - ", habitat_type))


fit <- vegan::envfit(ord_semi_nat ~ avg_values_MX2301A_temperature + avg_values_MX2202_light + ar5_PC1 + ar5_PC2,
                     data = points,
                     na.rm = TRUE)

fit_scores = data.frame(scores(fit, display = "vectors"))
fit_scores = cbind(fit_scores, vars = rownames(fit_scores))



## ---- results = "hide", warning = F, message = F, cache = cache_var-----------------------------------------------------------------------------
fig_seminat <- ggplot(points) +
  geom_point(data = points,
             aes(x = sites.NMDS1, y = sites.NMDS2,
                 color = region_name)) +
  # scale_color_manual(name = "Lokalitet",
  #                    values=paste0(Polychrome::palette36.colors())) +
  scale_color_nina(name = "Region") +
  coord_fixed() +
  theme(legend.position = "bottom",
       legend.key.size = unit(0.5, 'cm'),
       legend.text = element_text(size=7),
       legend.key.width = unit(0.5, 'cm')) +
  xlab("NMDS1") +
  ylab("NMDS2")

fig_seminat <- fig_seminat + 
  geom_segment(data = fit_scores[1:3,],
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               show.legend = F,
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_segment(data = fit_scores[4,],
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2 ),
               linetype=2,
               show.legend = F,
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_text(data = fit_scores,
            aes(x = NMDS1-0.05, y = NMDS2-0.05, label = c("temp","lys","ar5_PC1","ar5_PC2")),
            size = 3) + 
   ggtitle("Semi-Nat")


## ---- message = FALSE, results = "hide", cache = cache_var--------------------------------------------------------------------------------------

otu_loc_wider_skog <- otu_loc_raw %>% 
  filter(habitat_type == "Forest") %>% 
  pivot_wider(names_from = otu,
              values_from = no_reads,
              values_fill = 0,) %>% 
  as.data.frame()

rownames(otu_loc_wider_skog) <- otu_loc_wider_skog$sampling_name
otu_loc_skog <- otu_loc_wider_skog[-c(1:6)]


ord_skog <- vegan::metaMDS(otu_loc_skog, 
                              distance = "bray",
                              trymax = 100)
  
ord_skog_aug <- cbind(otu_loc_wider_skog[c(1:6)], vegan::scores(ord_skog)[1])
  
#points <- data.frame(ord_skog$points)

points <- ord_skog_aug %>% 
  left_join(locality_sampling_biomass_big_join_forest,
             by = c("sampling_name" = "sampling_name",
                    "locality" = "locality",
                    "region_name" = "region_name",
                    "habitat_type" = "habitat_type"))  %>% 
  mutate(region_habitat = paste0(region_name, " - ", habitat_type))


fit <- vegan::envfit(ord_skog ~ avg_values_MX2301A_temperature + avg_values_MX2202_light + ar5_PC1 + ar5_PC2,
                     data = points,
                     na.rm = TRUE)

fit_scores = data.frame(scores(fit, display = "vectors"))
fit_scores = cbind(fit_scores, vars = rownames(fit_scores))



## ---- results = "hide", warning = F, message = F------------------------------------------------------------------------------------------------
fig_skog <- ggplot(points) +
  geom_point(data = points,
             aes(x = sites.NMDS1, y = sites.NMDS2,
                 color = locality)) +
  # scale_color_manual(name = "Lokalitet",
  #                    values=paste0(Polychrome::palette36.colors())) +
  scale_color_nina(name = "Lokalitet") +
  coord_fixed() +
  theme(legend.position = "bottom",
       legend.key.size = unit(0.5, 'cm'),
       legend.text = element_text(size=7),
       legend.key.width = unit(0.5, 'cm')) +
  xlab("NMDS1") +
  ylab("NMDS2")

fig_skog <- fig_skog + 
  geom_segment(data = fit_scores[1:3,],
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               show.legend = F,
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_segment(data = fit_scores[4,],
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2 ),
               linetype=2,
               show.legend = F,
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_text(data = fit_scores,
            aes(x = NMDS1-0.05, y = NMDS2-0.05, label = c("temp","lys","ar5_PC1","ar5_PC2")),
            size = 3) + 
   ggtitle("Skog")


## ----beta-div-patterns-overall, fig.pos = "!h", fig.cap = "Ordinasjon som viser hovedmønstrene i betadiversitet i alle insektprøvene innsamlet i 2022. Punktene er farget etter innsamlingsdato og ellipsene viser et 95% konfidens-intervall for hvert habitat.\\label{fig:beta-div-patterns-overall}"----
fig_habitats


## ----betadiv-pattern-semi-nat, fig.pos = "!hb", fig.cap = "Ordinasjon som viser sammenheng mellom forklaringsvariabler og mønstrer i samfunnskomposisjon i semi-naturlig mark. Solide linjer viser signifikante sammenhenger med forklaringsvariabler.", out.width = "0.8\\textwidth"----
fig_seminat


## ----betadiv-pattern-skog, fig.pos = "!hb", fig.cap = "Ordinasjon som viser sammenheng mellom forklaringsvariabler og mønstrer i samfunnskomposisjon i skog. Solide linjer viser signifikante sammenhenger med forklaringsvariabler.", out.width = "0.8\\textwidth"----
fig_skog


## -----------------------------------------------------------------------------------------------------------------------------------------------
yl_catch_tab <- year_locality_biodiv %>% 
  left_join(year_locality_biomass) %>% 
  arrange(year, region_name, habitat_type, locality) %>% 
  select(År = year,
         Lokalitet = locality,
         Økosystemtype = habitat_type,
         Region = region_name,
         `Antall arter` = no_species,
         `Shannon div.` = shannon_div,
         `Våtvekt (middelv.)` = avg_wet_weight 
         )



## ----ylcatchtab, results = "asis"---------------------------------------------------------------------------------------------------------------
print(xtable(yl_catch_tab), 
      size = "footnotesize",
      floating = F,
      floating.environment = "table*",
      tabular.environment = "longtable")


## -----------------------------------------------------------------------------------------------------------------------------------------------
obs <- tbl(con,
           Id(schema = "occurrences",
              table = "observations")) %>% 
  filter(identification_confidence == 'HIGH')



## -----------------------------------------------------------------------------------------------------------------------------------------------
pot_alien_spec <- tbl(con,
                  Id(schema = "lookup",
                     table = "pot_alien_species2")) %>% 
  mutate(saul = as.character(saul),
         daisie = as.character(daisie)) %>% 
  mutate(saul = ifelse(saul == "false", "", "X"),
         daisie = ifelse(daisie == "false", "", "X")) %>% 
  select(Artsnavn = species_latin_fixed,
         "Listet i Saul" = saul,
         "Listet i Daisie" = daisie) %>% 
    arrange(Artsnavn) %>% 
  collect()


pot_alien_spec <- pot_alien_spec %>% 
  inner_join(obs,
            by = c("Artsnavn" = "species_latin_fixed"),
            copy = T) %>% 
  select(Order = "id_order",
         Artsnavn,
         "Listet i Saul",
         "Listet i Daisie") %>% 
  distinct()

fennoscand_spec <- tbl(con,
                  Id(schema = "lookup",
                     table = "fennoscand_species2"))  %>% 
    mutate(saul = as.character(saul),
         daisie = as.character(daisie)) %>% 
  mutate(saul = ifelse(saul == "false", "", "X"),
         daisie = ifelse(daisie == "false", "", "X")) %>% 
  select(Artsnavn = species_latin_fixed,
         "Listet i Saul" = saul,
         "Listet i Daisie" = daisie) %>% 
  arrange(Artsnavn) %>% 
  collect()

fennoscand_spec <- fennoscand_spec  %>% 
  inner_join(obs,
            by = c("Artsnavn" = "species_latin_fixed"),
            copy = T) %>% 
  select(Order = "id_order",
         Artsnavn,
         "Listet i Saul",
         "Listet i Daisie") %>% 
  distinct()



## -----------------------------------------------------------------------------------------------------------------------------------------------

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

alien_obs <- dbGetQuery(con,
                        alien_obs_q) %>% 
  mutate(habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type)) %>% 
  select(Orden = id_order,
         Artsnavn = species_latin_fixed,
         Risikovurdering = riskCategory,
         Habitatstype = habitat_type) %>% 
  arrange(Artsnavn)



## ---- results = "asis"--------------------------------------------------------------------------------------------------------------------------
print(xtable(alien_obs))

