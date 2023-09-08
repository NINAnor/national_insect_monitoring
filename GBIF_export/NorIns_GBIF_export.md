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
datasetID <- "19fe96b0-0cf3-4a2e-90a5-7c1c19ac94ee" ##From the webpage URL
# Suggested citation: Take the citation as from downloaded from GBIF website, replace "via GBIF.org" by endpoint url. 
tmp <- tempfile()
download.file(paste0("http://api.gbif.org/v1/dataset/",datasetID,"/document"),tmp) # get medatadata from gbif api
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
datasetURL <-  paste0("http://api.gbif.org/v1/dataset/",datasetID,"/endpoint")
dataset <- RJSONIO::fromJSON(datasetURL)
endpoint_url <- dataset[[1]]$url # extracting URL from API call result
```

The endpoint for our dataset looks like this: (Also visible on the
webpage for the dataset)

- https://ipt.nina.no/archive.do?r=national_insect_monitoring

With the endpoint, we can download the data in a zipped format.

``` r
# Download from dwc-a from IPT  putting into data folder
system("mkdir -p GBIF_data")
download.file(endpoint_url, destfile="GBIF_data/gbif_download.zip", mode="wb")
```

We unzip it in the same folder.

``` r
unzip("GBIF_data/gbif_download.zip",
      exdir = "GBIF_data")
```

## Read in the raw GBIF files

The data is a simple tab delimited text file.

``` r
eventRaw <- read_delim("GBIF_data/event.txt", 
                       delim = "\t",
                       locale = locale(encoding = "UTF-8"))
```

    Rows: 4287 Columns: 26
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr  (21): id, type, ownerInstitutionCode, dynamicProperties, eventID, paren...
    dbl   (4): sampleSizeValue, decimalLatitude, decimalLongitude, coordinateUnc...
    dttm  (1): modified

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

The occurrence data follows the same procedure.

``` r
occurrenceRaw <- read_delim("GBIF_data/occurrence.txt", 
                       delim = "\t",
                       locale = locale(encoding = "UTF-8")
                       )
```

    Rows: 327417 Columns: 20
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: "\t"
    chr  (16): id, basisOfRecord, occurrenceID, organismQuantityType, occurrence...
    dbl   (1): organismQuantity
    lgl   (2): sex, lifeStage
    dttm  (1): modified

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
