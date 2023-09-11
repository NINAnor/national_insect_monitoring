##Temporary file to store R-functions that could go into a separate r-package for the insect monitoring

checkCon <- function() {if(!exists("con")){ stop("No connection!")} else{
  if(class(con)!= "PqConnection"){ stop("\"con\" is not of class \"PqConnection\". Have you connected to the database?")}
  if(!DBI::dbIsValid(con)) { stop("No connection")}
}
}


calc_shannon <- function(community) {
  p <- table(community)/length(community) # Find proportions
  p <- p[p > 0] # Get rid of zero proportions (log zero is undefined)
  -sum(p * log(p)) # Calculate index
}




obs_from_db <- function(id_type = c("metabarcoding"),
                        subset_orders = NULL,
                        subset_families = NULL,
                        subset_genus = NULL,
                        subset_species = NULL,
                        subset_year = NULL,
                        subset_region = c(NULL, "Østlandet", "Trøndelag"),
                        trap_type = c("All", "MF", "VF"), 
                        limit = NULL,
                        dataset = c("NasIns"),
                        agg_level = c("year_locality",
                                      "locality_sampling",
                                      "region_habitat",
                                      "region_habitat_year",
                                      "total"),
                        as_tibble = F){
  
  checkCon()

  
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
              by = c("identification_name", "identification_name"),
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
    mutate(weeks_sampled = ifelse(grepl("2020", year) & (grepl("1", trap_short_name) | grepl("3", trap_short_name)), 2, 4)) %>% 
    mutate(weeks_sampled = ifelse(grepl("2020", year), weeks_sampled, 2))
  
  joined <- joined %>% 
    filter(weeks_sampled == 2)
  
  
  if(id_type == "metabarcoding"){
    joined <- joined %>% 
      filter(identification_type == "metabarcoding")
  }
  
  #Filter on region name
  if(!is.null(subset_region)){
    subset_region <- c("", subset_region)
    joined <- joined %>% 
      filter(region_name %IN% subset_region)
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
      filter(species_latin_fixed %in% subset_species)
  }
  
  if(!is.null(subset_year)){
    subset_year <- c("", subset_year)
    joined <- joined %>% 
      filter(year %IN% subset_year)
  }
  
  if(!is.null(subset_genus)){
    subset_genus <- c("", subset_genus)
    joined <- joined %>% 
      filter(id_genus %IN% subset_genus)
  }
  
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
  
  
  ##Aggregate data to choosen level
  ##Add more choices?
  
  res <- joined
  

  ##This is slow because we have to collect the data before we calculate Shannon index.
  ##Best would be to do the Shannon calc on the database side. Seems harder than I first thought.
  if(agg_level == "year_locality"){
  
      res <- res %>% 
        collect() %>% 
      group_by(year_locality_id, locality_id) %>% 
      summarise(no_species = n_distinct(species_latin_fixed),
                shannon_div = calc_shannon(species_latin_fixed)) %>% 
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
               shannon_div) %>% 
        arrange(year,
                region_name,
                habitat_type,
                locality)
    
  }
  
  
  if(agg_level == "locality_sampling"){
    
    res <- res %>% 
      collect() %>% 
      group_by(sampling_name, year_locality_id, locality_id) %>% 
      summarise(no_trap_days = mean(as.numeric(end_date_obs - start_date_obs)), ##to get the mean trap days from all traps within the sampling event (should be the same for all traps)
                no_species = n_distinct(species_latin_fixed),
                shannon_div = calc_shannon(species_latin_fixed)) %>% 
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
             shannon_div) %>% 
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
               habitat_type) %>% 
      summarise(no_species = n_distinct(species_latin_fixed),
                shannon_div = calc_shannon(species_latin_fixed)) %>% 
      ungroup() %>% 
      select(habitat_type,
             region_name,
             no_species,
             shannon_div) %>% 
      arrange(habitat_type,
              region_name)

  }
  
  
  if(agg_level == "region_habitat_year"){
    
    res <- res %>% 
      collect() %>% 
      group_by(region_name,
               habitat_type,
               year) %>% 
      summarise(no_species = n_distinct(species_latin_fixed),
                shannon_div = calc_shannon(species_latin_fixed)) %>% 
      ungroup() %>% 
      select(year,
             habitat_type,
             region_name,
             no_species,
             shannon_div) %>% 
      arrange(year,
              habitat_type,
              region_name
              )
    
  }
  
  if(agg_level == "total"){
    
    res <- res %>% 
      collect() %>% 
      summarise(no_species = n_distinct(species_latin_fixed),
                shannon_div = calc_shannon(species_latin_fixed)) %>% 
      ungroup() %>% 
      select(no_species,
             shannon_div) 
    
    
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



weights_from_db <- function(limit = NULL,
                            trap_type = c("ALL", "MF", "VF"),
                            subset_region = c(NULL, "Østlandet", "Trøndelag"),
                            dataset = c("NasIns"),
                            agg_level = c("year_locality",
                                          "locality_sampling",
                                          "total"),
                            as_tibble = F){
    
    checkCon()
    
    

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
                           choices = c("ALL", "MF", "VF", NULL))
    
    
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
    
    if(as_tibble){
      res <- res %>% 
        as_tibble()
    }
    
    
    
    return(res)  
    
  }




##Currently the temperature data is in the long format, should be faster to rearrange it in the database and skip the pivot_wider here.
##Note that this splits up the temperature data from the different logger types (since they are not placed similarly in the field)
loggerdata_from_db <- function(limit = NULL,
                        dataset = c("NasIns"),
                        agg_level = c("year_locality",
                                      "locality_sampling"),
                        as_tibble = F){
  
  checkCon()
  
  
  dataset <- match.arg(dataset, choices = c("NasIns", 
                                            "OkoTrond",
                                            "TidVar",
                                            "Nerlandsøya"))
  agg_level <- match.arg(agg_level, choices = c("year_locality",
                                                "locality_sampling"))
  
  
  ##Set up table sources
  ##Probably needs updating after new batch of data. Also need to test filtering of different identification types 
  #locality_sampling <- dplyr::tbl(con, dbplyr::in_schema("events", "locality_sampling"))
  #year_locality <- dplyr::tbl(con, dbplyr::in_schema("events", "year_locality")) 
  localities <- dplyr::tbl(con, dbplyr::in_schema("locations", "localities"))

  #logger_deployments <- dplyr::tbl(con, dbplyr::in_schema("loggers", "logger_deployments"))
  #logger_data <- dplyr::tbl(con, dbplyr::in_schema("loggers", "logger_data"))
  
  ##Do the pivot_wider on this side, temporarily, before I set this up in the database
  ##Not finished. Better to wait on a wider format database logger format and write a custom sql for the time interval join
  ##We also have the fuzzyjoin package that could do the trick, but it is logical to do the join and the aggregations 
  ## on the database side anyway.
  
  
  ##Aggregate data to choosen level
  ##Add more choices!
  
  if(agg_level == "year_locality"){
 
    yl_query <- paste0("
    SELECT yl.year, l.id as locality_id,ld.logger_type, ld.data_type, 
    avg(ld.value) as avg_values, 
    max(ld.value) as max_values,
    min(ld.value) as min_values,
    stddev_pop(ld.value) as stddev_values
    FROM loggers.logger_data ld,
    loggers.logger_deployments ldep,
    events.year_locality yl,
    locations.localities l
    WHERE ld.logger_id = ldep.logger_id
    AND ldep.year_locality_id = yl.id
    AND yl.start_date::timestamp <= ld.logger_time
    AND yl.end_date::timestamp >= ld.logger_time
    AND yl.locality_id = l.id
    AND yl.project_short_name = '", 
                       dataset,
                       "'
    GROUP BY yl.id, l.id, ld.logger_type, ld.data_type
     ")
    
    yl_data_raw <- dbGetQuery(con,
                              yl_query)
    
    yl_data <- yl_data_raw %>% 
      pivot_wider(id_cols = c("locality_id", "year"),
                  values_from = c("avg_values", "max_values", "min_values", "stddev_values"),
                  names_from = c(logger_type, data_type),
                  names_prefix = "") 
    
    loc_info <- localities %>% 
      select(id, 
             locality,
             habitat_type,
             region_name)
    
    res <- yl_data %>% 
      left_join(loc_info,
                by = c("locality_id" = "id"),
                copy = T) %>% 
      select(-"locality_id") %>% 
      select(locality,
             habitat_type,
             region_name,
             everything()) %>% 
    arrange(year,
            region_name,
            habitat_type,
            locality)

   }
  
  
  ##Works but could be faster? check the sql join and indices
  if(agg_level == "locality_sampling"){
    
   ls_query <- paste0("
    SELECT yl.year, l.id as locality_id, ls.sampling_name, ld.logger_type, ld.data_type, 
    avg(ld.value) as avg_values, 
    max(ld.value) as max_values,
    min(ld.value) as min_values,
    stddev_pop(ld.value) as stddev_values
    FROM loggers.logger_data ld,
    loggers.logger_deployments ldep,
    events.year_locality yl,
    locations.localities l,
	events.locality_sampling ls
    WHERE ld.logger_id = ldep.logger_id
    AND ldep.year_locality_id = yl.id
	AND ls.year_locality_id = yl.id
	AND yl.locality_id = l.id
    AND ls.start_date::timestamp <= ld.logger_time
    AND ls.end_date::timestamp >= ld.logger_time
    AND yl.project_short_name = '", 
                       dataset,
                       "'
    GROUP BY ls.id, yl.year, l.id, ld.logger_type, ld.data_type
	ORDER BY sampling_name
    ")
    
    ls_data_raw <- dbGetQuery(con,
                              ls_query) %>% 
      as_tibble()
    
    ls_data <- ls_data_raw %>% 
      pivot_wider(id_cols = c("year", "locality_id", "sampling_name"),
                  values_from = c("avg_values", "max_values", "min_values", "stddev_values"),
                  names_from = c(logger_type, data_type),
                  names_prefix = "") 
    
    loc_info <- localities %>% 
      select(id, 
             locality,
             habitat_type,
             region_name)
    
    res <- ls_data %>% 
      left_join(loc_info,
                by = c("locality_id" = "id"),
                copy = T) %>% 
      select(-"locality_id") %>% 
      select(locality,
             habitat_type,
             region_name,
             everything()) %>% 
      arrange(year,
              region_name,
              habitat_type,
              locality)
    
  
    
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

#Skogtype agg from db

skogtype_from_db <- function(limit = NULL,
                             dataset = c("NasIns"),
                             stat_grouping = c("treslag", "bonitet"), 
                             as_tibble = F) {
  
  
  dataset <- match.arg(dataset, choices = c("NasIns", 
                                            "OkoTrond",
                                            "TidVar",
                                            "Nerlandsøya"))
  
  proj_filter <- dplyr::enquo(dataset)
  
  
  yl <- dplyr::tbl(con,
                   Id(schema = "events",
                      table = "year_locality")
  )
  l <- dplyr::tbl(con,
                  Id(schema = "locations",
                     table = "localities")
  )
  
  year_locality <- l %>% 
    left_join(yl,
              by = c("id" = "locality_id")) %>% 
    filter(project_short_name == !!proj_filter) %>% 
    select(locality) %>% 
    pull()
  
  
  skog <- read_sf(con,
                  Id(schema = "backgrounds",
                     table = "locality_1000m_buffer_skogtype_agg")
  ) %>%
    st_drop_geometry() %>% 
    filter(locality %in% year_locality)
  
  
  group_enq <- dplyr::enquo(stat_grouping)
  
  if(length(stat_grouping) == 1){
    skog <- skog %>% 
      group_by(across(all_of(c("locality", stat_grouping)))) %>% 
      summarise(area = sum(area))
    }
  
  
  out <-  skog %>% 
    pivot_wider(id_cols = locality,
                names_from = !!group_enq,
                values_from = area,
                names_prefix = "area_") 
  
  
  out <- out %>% 
    mutate_if(is.numeric,
              replace_na,
              replace = 0)
  
  if(as_tibble){
    out <- out %>% 
      as_tibble()
  }
  
  
  return(out)
}



#loggerdata_from_db(agg_level = "locality_sampling")



longerHobo2301 <- function(inputFile,
                           guess_max = 10000,
                           ...){
  
  rawDat <- read_csv(inputFile,
                     col_types = cols(.default = "c"),
                     guess_max = guess_max,
                     ...)
  
  dat <- rawDat %>%  
    select(-"Line#") %>% 
    mutate(date = as.POSIXct(Date, format = "%m/%d/%y %H:%M:%S")) %>% 
    mutate_if(is_character, as.double) %>% 
    select(-Date)
  
  
  temp <- dat %>% 
    pivot_longer(cols = starts_with("Temperature"),
                 names_to = "logger_id",
                 values_to = "temperature") %>% 
    select(date,
           logger_id,
           temperature) %>% 
    filter(!is.na(temperature))
  
  rh <- dat %>% 
    pivot_longer(cols = starts_with("RH"),
                 names_to = "logger_id",
                 values_to = "rh") %>% 
    select(date,
           logger_id,
           rh)%>% 
    filter(!is.na(rh))
  
  dew_point  <- dat %>% 
    pivot_longer(cols = starts_with("Dew"),
                 names_to = "logger_id",
                 values_to = "dew_point") %>% 
    select(date,
           logger_id,
           dew_point) %>% 
    filter(!is.na(dew_point))
  
  
  #Fix to allow for two deployments of same logger. gets duplicate column names from hobo-export
  temp <- temp %>% 
    mutate(logger_id = str_extract(logger_id,
                                   "[^, ]+$")) %>% 
    mutate(logger_id = str_extract(logger_id,
                                   "(^[0-9]*)"))
    
  rh <- rh %>% 
    mutate(logger_id = str_extract(logger_id,
                                   "[^, ]+$")) %>% 
    mutate(logger_id = str_extract(logger_id,
                                   "(^[0-9]*)"))
  dew_point <- dew_point %>% 
    mutate(logger_id = str_extract(logger_id,
                                   "[^, ]+$")) %>% 
    mutate(logger_id = str_extract(logger_id,
                                   "(^[0-9]*)"))
  
  if(!all(all(temp$date == rh$date),
          all(rh$date == dew_point$date))) stop("Tables datetimes doesn't match")
  
  combDat <- temp %>% 
    full_join(rh,
              by = c("date" = "date",
                     "logger_id" = "logger_id")) %>% 
    full_join(dew_point,
              by = c("date" = "date",
                     "logger_id" = "logger_id")) %>% 
    arrange(logger_id,
            date) %>% 
    mutate(logger_type = "MX2301A") %>% 
    select(date, 
           logger_type,
           logger_id,
           temperature,
           rh,
           dew_point)
  
  return(combDat)
}






longerHobo2202 <- function(inputFile,
                           guess_max = 10000){
  rawDat <- read_csv(inputFile,
                     guess_max = guess_max,
                     col_types = cols())
  
  dat <- rawDat %>%  
    select(-"Line#") %>% 
    mutate(date = as.POSIXct(Date, format = "%m/%d/%y %H:%M:%S")) %>% 
    mutate_if(is_character, as.double) %>% 
    select(-Date)
  
  
  
  temp <- dat %>% 
    pivot_longer(cols = starts_with("Temperature"),
                 names_to = "logger_id",
                 values_to = "temperature") %>% 
    select(date,
           logger_id,
           temperature) %>% 
    filter(!is.na(temperature))
  
  light <- dat %>% 
    pivot_longer(cols = starts_with("Light"),
                 names_to = "logger_id",
                 values_to = "light") %>% 
    select(date,
           logger_id,
           light)%>% 
    filter(!is.na(light))
  
  
  
  temp <- temp %>% 
    mutate(logger_id = str_extract(logger_id,
                                   "[^, ]+$"))
  light <- light %>% 
    mutate(logger_id = str_extract(logger_id,
                                   "[^, ]+$"))
  
  if(!all(temp$date == light$date)) stop("Tables datetimes doesn't match")
  
  combDat <- temp %>% 
    full_join(light,
              by = c("date" = "date",
                     "logger_id" = "logger_id")) %>% 
    arrange(logger_id,
            date) %>% 
    mutate(logger_type = "MX2202") %>% 
    select(date, 
           logger_type,
           logger_id,
           temperature,
           light)
  
  return(combDat)
}




#Function to make list of rasters out of many localities
#To be used with fasterize

rast_list <- function(input){
  
  landscapes <- input %>% 
    st_drop_geometry() %>% 
    select(locality) %>% 
    distinct() %>% 
    arrange(locality) %>% 
    pull()
  
  out <- list()
  
  for(i in landscapes){
    
    index <- dplyr::enquo(i)
    
    r <- input %>% 
      filter(locality == !!i) %>% 
      fasterize::raster(., res = 10)
    
    vect <- input %>% 
      filter(locality == i)
    
    out[[i]] <- fasterize::fasterize(vect,
                          r,
                          field = "arealressursArealtype"
    )
    
  } 
  
  return(out)
}





#######AR5 colors
ar5_colors <- tibble(arealtype = c("Fulldyrka jord",
                                   "Overflatedyrka jord",
                                   "Innmarksbeite",
                                   "Skog",
                                   "Myr",
                                   "Åpen fastmark",
                                   "Ferskvann",
                                   "Hav",
                                   "Snøisbre",
                                   "Bebygd",
                                   "Samferdsel",
                                   "Ikke kartlagt"),
                     red = c(255,
                             255,
                             255,
                             158,
                             209,
                             217,
                             145,
                             204,
                             230,
                             252,
                             179,
                             255),
                     green = c(209,
                               255,
                               255,
                               204,
                               209,
                               217,
                               231,
                               254,
                               255,
                               219,
                               120,
                               255),
                     blue = c(110,
                              76,
                              173,
                              115,
                              255,
                              217,
                              255,
                              254,
                              255,
                              214,
                              76,
                              255)) %>% 
  mutate(color = rgb(red, green, blue,
                     maxColorValue = 255))

ar5_colors_vect <- ar5_colors %>% 
  select(color) %>% 
  pull()
names(ar5_colors_vect) <- ar5_colors$arealtype


scale_fill_ar5 <- function(...) {
  scale_fill_manual(values = ar5_colors_vect,
                    ...)
}


#Treslag scale
treslag_colors <- tibble(treslag = c("Barskog",
                                   "Lauvskog",
                                   "Blandingsskog",
                                   "Ikke tresatt",
                                   "Ikke relevant",
                                   "Ikke registrert"),
                     red = c(125,
                             128,
                             158,
                             207,
                             255,
                             255
                             ),
                     green = c(191,
                               255,
                               204,
                               204,
                               255,
                               255),
                     blue = c(110,
                              8,
                              115,
                              145,
                              255,
                              255)) %>% 
  mutate(color = rgb(red, green, blue,
                     maxColorValue = 255))

treslag_colors_vect <- treslag_colors %>% 
  select(color) %>% 
  pull()
names(treslag_colors_vect) <- treslag_colors$treslag


scale_fill_treslag <- function(...) {
  scale_fill_manual(values = treslag_colors_vect,
                    ...)
}


#Bonitet scale
bonitet_colors <- tibble(bonitet = c("Særs høg",
                                     "Høg",
                                     "Middels",
                                     "Lav",
                                     "Impediment",
                                     "Ikke relevant",
                                     "Ikke registrert"),
                         red = c(0,
                                 125,
                                 158,
                                 209,
                                 235,
                                 255,
                                 255
                                 ),
                         green = c(173,
                                   191,
                                   204,
                                   232,
                                   245,
                                   255,
                                   255),
                         blue = c(59,
                                  110,
                                  115,
                                  181,
                                  209,
                                  255,
                                  255
                                  )
                         ) %>% 
  mutate(color = rgb(red, green, blue,
                     maxColorValue = 255))

bonitet_colors_vect <- bonitet_colors %>% 
  select(color) %>% 
  pull()
names(bonitet_colors_vect) <- bonitet_colors$bonitet


scale_fill_bonitet <- function(...) {
  scale_fill_manual(values = bonitet_colors_vect,
                    ...)
}




locality_ar5_plot <- function(locality,
                              ano_circles = T,
                              title = T, 
                              legend = T,
                              strip = F,
                              ...){
  
  loc_filter <- dplyr::enquo(locality)
  
  nas_ins_lok_1000m_agg_q <- paste0("
  SELECT ar5.*, ar5_kode.egenskapsverdi as arealtype
  FROM backgrounds.locality_1000m_buffer_arealtype_agg ar5,
  lookup.ar5_arealtype ar5_kode,
  events.year_locality yl
  WHERE ar5.locality_id = yl.locality_id
  AND yl.project_short_name = 'NasIns'
  AND ar5.\"arealressursArealtype\" = ar5_kode.kodeverdi",
                                    "\nAND ar5.locality = '", 
                                    locality,
                                    "';")
  
  
  lok_1000m_agg <- read_sf(con,
                           query = nas_ins_lok_1000m_agg_q)
  
  
  lok <- st_read(con, 
                 layer = Id(schema = "locations",
                            table = "localities"),
                 geometry_column = "geom") %>% 
    filter(locality == !!loc_filter) 
  
  
  traps <- st_read(con,
                   layer = Id(schema = "locations",
                              table = "traps"),
                   geometry_column = "geom") %>% 
    filter(locality == !!loc_filter) 
  
  
  
  p <- lok_1000m_agg %>% 
    ggplot() +
    geom_sf(aes(fill = arealtype)) +
    scale_fill_ar5(name = "Arealtype AR5") +
    geom_sf(col = "red",
            data = traps) +
    geom_sf(fill = NA,
            col = "blue",
            lwd = 1,
            data = lok) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  if(ano_circles){
    
    ano_geoms_q <- "
  SELECT l.locality,
  sp.ano_flate_id,
  sp.ano_punkt_id,
  ano_b.geom
  FROM ano.survey_points sp,
  ano.ano_18_250m_buffers ano_b,
  locations.localities l,
  events.year_locality yl
  WHERE sp.ano_flate_id = ano_b.ano_id
  AND sp.ano_punkt_id = ano_b.ano_point_id
  AND yl.ano_flate_id = sp.ano_flate_id
  AND yl.locality_id = l.id
  " 
    
    ano_geoms <- read_sf(con,
                         query = ano_geoms_q) %>% 
      filter(locality == !!loc_filter)
    
    
    p <- p +
      geom_sf(data = ano_geoms,
              fill = NA,
              col = nina_colors["purple"])
      
  }
  
  if(title){
    p <- p +
      ggtitle(locality)  
    
  }
  
  if(!legend){
    p <- p +
      guides(fill = "none")
  }
  
  if(strip){
    p <- p +
      theme(axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            legend.position="none",
            panel.background=element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank()) +
      theme(plot.margin = margin(0, 0, 0, 0))
  }
  
  return(p)
  
  
} 




chart.Correlation2 <- function (R, 
                                histogram = TRUE, 
                                method = c("pearson", 
                                           "kendall", 
                                           "spearman"), 
                                color = "red",
                                ...) 
{
  x = checkData(R, method = "matrix")
  if (missing(method)) 
    method = method[1]
  cormeth <- method
  panel.cor <- function(x, y, digits = 2, prefix = "", use = "pairwise.complete.obs", 
                        method = cormeth, cex.cor, ...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y, use = use, method = method)
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste(prefix, txt, sep = "")
    if (missing(cex.cor)) 
      cex <- 0.8/strwidth(txt)
    test <- cor.test(as.numeric(x), as.numeric(y), method = method)
    Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                     cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", 
                                                                              "**", "*", ".", " "))
    text(0.5, 0.5, txt, cex = cex * (abs(r) + 0.3)/1.3)
    text(0.8, 0.8, Signif, cex = cex, ..., col = color)
  }
  f <- function(t) {
    dnorm(t, mean = mean(x), sd = sd.xts(x))
  }
 dotargs <- list(...)
dotargs$method <- NULL
 rm(method)
  hist.panel = function(x, ... = NULL) {
    par(new = TRUE)
    hist(x, col = "light gray", probability = TRUE, axes = FALSE, 
         main = "", breaks = "FD")
    lines(density(x, na.rm = TRUE), lwd = 1, col = color)
    rug(x)
  }
  if (histogram) 
    suppressWarnings(pairs(x, gap = 0, lower.panel = panel.smooth, upper.panel = panel.cor, 
          diag.panel = hist.panel, col.smooth = color))
  else pairs(x, gap = 0, lower.panel = panel.smooth, upper.panel = panel.cor)
}






community_matrix_from_db <- function(limit = NULL,
                                     id_type = c("metabarcoding"),
                                     trap_type = c("MF", "VF", "All", NULL),
                                     dataset = c("NasIns"),
                                     subset_orders = NULL,
                                     subset_families = NULL,
                                     subset_species = NULL,
                                     subset_habitat = NULL,
                                     subset_region = c(NULL, "Østlandet", "Trøndelag"),
                                     exclude_singletons = F,
                                     transposed_matrix = F,
                                     as_tibble = F){
  
  checkCon()
  
  
  dataset <- match.arg(dataset, 
                       choices = c("NasIns", 
                                   "OkoTrond",
                                   "TidVar",
                                   "Nerlandsøya"))
  
  trap_type <- match.arg(trap_type,
                         choices = c("MF", "VF", "All", NULL))
  
  
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
              by = c("identification_name", "identification_name"),
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


toiNEXT <- function(input){
  
  all_names <- lapply(input, rownames) %>% unlist() %>% unique()
  
  missing_names <- lapply(input, function(x) all_names[!(all_names %in% rownames(x))])
  
  to_add <- list()
  
  for(i in 1:length(input)){
    to_add[[i]] <- matrix(0, 
                          nrow = length(missing_names[[i]]), 
                          ncol = ncol(input[[i]]))
    
    rownames(to_add[[i]]) = missing_names[[i]]
    
  }
  
  
  added = list()
  for(i in 1:length(input)){
    added[[i]] <- rbind(input[[i]], to_add[[i]])
    
  }
  
  out <- lapply(added, function(x) x[sort(row.names(x)), ])
  
  
  return(out)
}


plot_beta_part <- function(input){
  
  ylims <- c(0, max(input$beta_sim)*1.1)
  
  p_sor <- ggplot(input) +
    geom_point(aes(x = distance,
                   y = beta_sor)) +
    geom_smooth(aes(x = distance,
                    y = beta_sor),
                method = "lm",
                color = ninaPalette()[2]) +
    ylim(ylims) +
    ylab(unname(TeX("Total forskjell $\\beta_{sor}$"))) +
    xlab("Avstand mellom lokaliteter (km)") 
  
  
  
  p_sim <- ggplot(input) +
    geom_point(aes(x = distance,
                   y = beta_sim)) +
    geom_smooth(aes(x = distance,
                    y = beta_sim),
                method = "lm",
                color = ninaPalette()[2]) +
    ylim(ylims) +
    ylab(unname(TeX("Turnover $\\beta_{sim}$"))) +
    xlab("Avstand mellom lokaliteter (km)") 
  
  p_nes <- ggplot(input) +
    geom_point(aes(x = distance,
                   y = beta_sne)) +
    geom_smooth(aes(x = distance,
                    y = beta_sne),
                method = "lm",
                color = ninaPalette()[2]) +
    ylim(ylims) +
    ylab(unname(TeX("Nestedness $\\beta_{nes}$"))) +
    xlab("Avstand mellom lokaliteter (km)") 
  
  marrangeGrob(list(p_sor,
                    p_sim,
                    p_nes
  ),
  top = "",
  nrow = 3,
  ncol = 1)
  
} 

combine_dist_to_comm_mat <- function(comm_mat,
                                     region_name,
                                     habitat_type){
  
  #It ain't pretty but it saves some code in the document
  
  dist_q <- paste0("
  SELECT a.locality as loc_a, b.locality loc_b, 
  ST_Distance(ST_Centroid(a.geom), ST_Centroid(b.geom)) as distance
  FROM 
  (SELECT l.*
  FROM locations.localities l,
  events.year_locality yl
  WHERE yl.locality_id = l.id\
  AND region_name IN ", region_name," 
  AND habitat_type = '", habitat_type, "' 
  --AND yl.year = 2020
  AND yl.project_short_name = 'NasIns') a 
  CROSS JOIN (SELECT l.*
  FROM locations.localities l,
  events.year_locality yl
  WHERE yl.locality_id = l.id
  AND region_name IN ", region_name," 
  AND habitat_type = '", habitat_type, "' 

  AND yl.project_short_name = 'NasIns') as b 
  ORDER BY loc_a, loc_b
  ")
  
  dist <- dbGetQuery(con,
                     dist_q)
  
  
  ##Now arrange these together and plot
  #something like this?
  
  beta_pair <- beta.pair(betapart.core(comm_mat))
  
  beta_sim <- as.matrix(beta_pair$beta.sim) %>% 
    matrix(.,
           ncol = 1)
  beta_sne <- as.matrix(beta_pair$beta.sne) %>% 
    matrix(.,
           ncol = 1)
  beta_sor <- as.matrix(beta_pair$beta.sor) %>% 
    matrix(.,
           ncol = 1)

  dist_beta <- dist %>% 
    cbind(beta_sim) %>% 
    cbind(beta_sne) %>% 
    cbind(beta_sor) %>% 
    filter(distance != 0) %>% 
    as_tibble() %>% 
    mutate(distance = distance / 1000)
  
  return(dist_beta)
  
}
