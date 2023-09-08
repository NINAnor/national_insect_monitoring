GBIF data export from the Norwegian insect monitoring program (NorIns)
================
Jens Åström
9/8/23

``` r
suppressPackageStartupMessages({
require(tidyverse)
require(NinaR) #optional
require(DBI)
require(RPostgres)
require(sf)
require(maps)
require(tidyjson)
require(xml2)
})
```

# Introduction

This script downloads observation records from the Norwegian insect
monitoring program (NorIns), and unpacks the data for further use. The
data is stored within GBIF, and the first section shows how to download
the data from there. The rest of the script deals with rearranging the
data from GBIFs eventCore standard, where several hierarchical levels is
compressed into a single event table. Portions of the script can be used
as an example for projects with other data formatting preferences. The
data is freely available, given that you site the source appropriately.
See “Citation and meta-data” below.

## Data structure

The data collection of the program has a hierarchical nature, where the
data from each respective level is originally stored as separate tables
in a database (normalised data). For the purpose of this GBIF dataset,
at the topmost level 4), we have a whole sampling season (april-october)
of each locality, where data associated with the entire locality-season
is stored, such as habitat qualities. At level 3), we have individual
sampling periods (of typically 2 weeks) of one or multiple traps within
the locality. This data includes e.g. time period and associated
weather. Each season can have up to 14 such sampling periods. At level
2), we have data from a single trap within a sampling period, such as
biomass and ethanol concentration. Lastly, at level 1, we have insect
identifications of each trap sample. Most samples only have 1
identification “event”, although multiple identification methods can be
employed some samples.

The hierarchical structure is handled within GBIF through the Darwin
Event Core-standard. Here, all the actual observation data of the
species are located in an “occurrence”-table, which is linked to an
“event”-table containing the information linked to the observation. The
event-table has to handle all the hierarchical levels mentioned above,
and thus combines data from many original tables. The hierarcical nature
of these tables are maintained through “parent-event”-keys, that relate
lower levels to their respective higher “parent” levels. By this design
of GBIF, you can have an arbitrary number of parent events (hierarchical
levels) within the same table, instead of normalizing the data across
several tables. The downside is that you need to keep track of the
eventID’s and the parentEventID’s to recreate the original data
structure.

# Finding the raw data and how to cite it

The dataset can be found at http://gbif.org by searching for “National
insect monitoring in Norway”. This will take you to the webpage of the
dataset:
https://www.gbif.org/dataset/19fe96b0-0cf3-4a2e-90a5-7c1c19ac94ee Here
you can access the meta-data of the dataset. This data is freely
available, but please cite the source. The citation is visible at the
dataset’s webpage at GBIF. The citation can be fetched programmatically
this way.

``` r
dataset_id <- "19fe96b0-0cf3-4a2e-90a5-7c1c19ac94ee" ##From the webpage URL
# Suggested citation: Take the citation as from downloaded from GBIF website, replace "via GBIF.org" by endpoint url. 
tmp <- tempfile()
download.file(paste0("http://api.gbif.org/v1/dataset/",dataset_id,"/document"),tmp) # get medatadata from gbif api
meta <- read_xml(tmp) %>% as_list() # create list from xml schema
gbif_citation <- meta$eml$additionalMetadata$metadata$gbif$citation[[1]] # extract citation
```

The citation to use is then:

- Åström J (2023). National insect monitoring in Norway. Version 1.9.
  Norwegian Institute for Nature Research. Sampling event dataset
  https://doi.org/10.15468/2jwnc6 accessed via GBIF.org on 2023-09-08.

# Fetching the raw-data from GBIF

We need an “endpoint_url” to retrieve the actual data set. This is found
at the bottom of the webpage for the dataset. This can also be found
programatically as shown here.

``` r
dataset_url <-  paste0("http://api.gbif.org/v1/dataset/",dataset_id,"/endpoint")
dataset <- RJSONIO::fromJSON(dataset_url)
endpoint_url <- dataset[[1]]$url # extracting URL from API call result
```

The endpoint for our dataset looks like this: (Also visible on the
webpage for the dataset)

- https://ipt.nina.no/archive.do?r=national_insect_monitoring

With the endpoint, we can download the data in a zipped format.

``` r
# Download from dwc-a from IPT  putting into data folder
system("mkdir -p GBIF_data")
download.file(endpoint_url, destfile = "GBIF_data/gbif_download.zip", mode="wb")
```

We unzip it in the same folder.

``` r
unzip("GBIF_data/gbif_download.zip",
      exdir = "GBIF_data")
```

## Read in the raw GBIF files

The data is a simple tab delimited text file.

``` r
event_raw <- read_delim("GBIF_data/event.txt", 
                       delim = "\t",
                       locale = locale(encoding = "UTF-8"),
                       progress = FALSE,
                       show_col_types = FALSE)
```

The occurrence data follows the same procedure.

``` r
occurrence_raw <- read_delim("GBIF_data/occurrence.txt", 
                       delim = "\t",
                       locale = locale(encoding = "UTF-8"),
                       progress = FALSE,
                       show_col_types = FALSE
                       )
```

# Extracting the data from the GBIF raw format

We next split out the separate levels from the combined event table.

``` r
identifications_raw <- event_raw %>% 
  filter(samplingProtocol == "Level_1")

sampling_trap_raw <- event_raw %>% 
  filter(samplingProtocol == "Level_2")

locality_sampling_raw <- event_raw %>% 
  filter(samplingProtocol == "Level_3")

year_locality_raw <- event_raw %>% 
  filter(samplingProtocol == "Level_4")
```

## Extract the data that didn’t fit into the GBIF standard

In all levels, we have information that don’t fit into the GBIF data
format. These are bundled together in to a JSON-string in a GBIF column
named “dynamicProperties”. We need to spread this content into separate
columns.

``` r
identifications_raw %>% 
  select(dynamicProperties) %>% 
  head(1)
```

    # A tibble: 1 × 1
      dynamicProperties                                                             
      <chr>                                                                         
    1 "{\"id\":\"cc5c4506-d23c-4db4-be1c-cc1318c63ec5\",\"identification_name\":\"m…

Here, we split out the json string (dynamicProperties) into a separate
table, and join it back to the original data afterwards. I tried various
json-packages earlier with limited success, but this works and should be
stable.

``` r
identifications <- identifications_raw %>% 
  #collect() %>% 
  select(-dynamicProperties) %>% 
  mutate(`document.id` = 1:nrow(.))

tempDynamic <- identifications_raw %>% 
  #head() %>% 
  select(dyn = dynamicProperties) %>% 
  #collect() %>% 
  unlist() %>% 
  spread_all()  %>% 
  as_tibble()


identifications <- identifications %>% 
  left_join(tempDynamic,
            by = c("document.id" = "document.id",
                   "id" = "id"))
```

The new columns contain an identification method name, optional
identification comment, read_abundance (number of DNA-copies),
extraction date, extraction-kit, and optional extraction comment. A
sample of the extracted data:

``` r
identifications %>% 
  select(identification_name,
         identification_comment,
         read_abundance,
         ekstraksjonsdato,
         ekstraksjonskit,
         ekstraksjonskommentar) %>% 
  head()
```

    # A tibble: 6 × 6
      identification_name   identification_comment read_abundance ekstraksjonsdato
      <chr>                 <lgl>                           <dbl> <lgl>           
    1 metabarcoding_novaseq NA                             556742 NA              
    2 metabarcoding_novaseq NA                             390776 NA              
    3 metabarcoding_novaseq NA                             249041 NA              
    4 metabarcoding_myseq   NA                              31057 NA              
    5 metabarcoding_myseq   NA                              83008 NA              
    6 metabarcoding_novaseq NA                            3557357 NA              
    # ℹ 2 more variables: ekstraksjonskit <lgl>, ekstraksjonskommentar <lgl>

## Continue with other levels…