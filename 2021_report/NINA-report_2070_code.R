## ---- include = F---------------------------------------------------------------------------------------------------------
# Some common packages, loading rmarkdown doesn't like the messages from tidyverse, so we don't include this in the document'
require(tidyverse)
require(DBI)
require(RPostgres)
require(ggplot2)
require(xtable)
require(NinaR)
require(sf)
require(ggplot2)
require(ggforce)
#require(bookdown)
#require(knitr)
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


## ---- include = F, eval = T-----------------------------------------------------------------------------------------------
#To log into something with your own user, it's good practice not to store passwords in the script.
##I store the username and passwords in a file in my home folder, that only the owner can read.
#Not the prettiest solution, since I read in the info in plain text, but we delete the objects in the same code chunk so it works. 
##Here I also specify that big integers should be converted to integers in R. 
##There are some functions that don't have support for bigints in R at the moment.

source("~/.rpgpass")

postgreSQLConnect(
  host = "ninradardata01.nina.no",
  dbname = "insect_monitoring",
  username = username,
  password = password,
  bigint  = "integer"
)

rm(list = c("username", "password"))


## ----setup, echo=F, message=F, warning=F----------------------------------------------------------------------------------
require(NinaR)
require(xtable)
require(knitr)
opts_chunk$set(tidy=T, 
               tidy.opts=list(width.cutoff=50),
               echo = F,
               message = F,
               warning = F,
               dpi = 600,
               dev = c("pdf", "png"),
               fig.path = "figure/")
options(width = 60,
        "xtable.comment" = F,
        "xtable.include.rownames" = F)

library(ggthemes)
theme_set(theme_few())



## -------------------------------------------------------------------------------------------------------------------------
source("../../functions.R")


## -------------------------------------------------------------------------------------------------------------------------
richn_2020 <- obs_from_db(subset_year = 2020,
                          agg_level = "total") %>% 
  select(no_species) %>% 
  pull()

richn_2021 <- obs_from_db(subset_year = 2021,
                          agg_level = "total") %>% 
  select(no_species) %>% 
  pull()

richn_tot <- obs_from_db(agg_level = "total") %>% 
  select(no_species) %>% 
  pull()



## ---- eval = F------------------------------------------------------------------------------------------------------------
## weights_from_db(agg_level = "total")


## -------------------------------------------------------------------------------------------------------------------------
#Note, includes all data not only NasIns
to_species_q <-"
SELECT distinct on(species_latin_fixed) obs.taxonomic_level, obs.species_latin_fixed
FROM occurrences.observations obs LEFT JOIN
events.identifications ids ON obs.identification_id = ids.id LEFT JOIN
lookup.identification_techniques id_tech ON ids.identification_name = id_tech.identification_name LEFT JOIN
events.sampling_trap st ON sampling_trap_id = st.id LEFT JOIN
events.locality_sampling ls ON locality_sampling_id = ls.id LEFT JOIN
events.year_locality yl ON year_locality_id = yl.id
WHERE yl.year = 2021
AND yl.project_short_name = 'NasIns'
"

to_species <- dbGetQuery(con,
                         to_species_q) %>% 
  group_by(taxonomic_level) %>% 
  summarise(no_records = n()) %>% 
  ungroup() %>% 
  mutate(perc = no_records / sum(no_records))



## -------------------------------------------------------------------------------------------------------------------------
load(file = "out/ost_iNEXT.rda")


no_locality_samplings_q <-"
SELECT count(distinct ls.id)
FROM events.locality_sampling ls,
events.year_locality yl
WHERE yl.id = ls.year_locality_id
AND yl.project_short_name = 'NasIns'

" 

no_locality_samplings <- dbGetQuery(con,
                                    no_locality_samplings_q)


## -------------------------------------------------------------------------------------------------------------------------
localities <- tbl(con,
                  Id(schema = "locations",
                     table = "localities"))
year_locality <- tbl(con,
                     Id(schema = "events",
                        table = "year_locality")) 

locality_sampling <- tbl(con,
                         Id(schema = "events",
                            table = "locality_sampling"))


## -------------------------------------------------------------------------------------------------------------------------
year_locality_biodiv <- obs_from_db(agg_level = "year_locality",
                                 dataset = "NasIns",
                                 as_tibble = T)


## -------------------------------------------------------------------------------------------------------------------------
locality_sampling_biodiv <- obs_from_db(agg_level = "locality_sampling",
                                 dataset = "NasIns",
                                 as_tibble = T)
#locality_sampling_biodiv


## -------------------------------------------------------------------------------------------------------------------------
biomass_mf_locality_sampling <- weights_from_db(agg_level = "locality_sampling",
                                                dataset = "NasIns",
                                                trap_type = "MF")

#biomass_mf_locality_sampling

## -------------------------------------------------------------------------------------------------------------------------
year_locality_biomass <- weights_from_db(dataset = "NasIns",
                                             trap_type = "MF",
                                             agg_level = "year_locality")


## -------------------------------------------------------------------------------------------------------------------------
nas_ins_lok <- read_sf(con,
                       Id(schema = "locations",
                       table = "localities")) %>% 
  filter(grepl("Semi-", locality) |
         grepl("Skog_", locality)
         ) %>% 
  mutate(habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type))



## -------------------------------------------------------------------------------------------------------------------------
norge <- read_sf(con,
                 Id(schema = "backgrounds",
                    table = "norway_terrestrial")
                 ) 


## ----lokaler_2020_2021, fig.cap = "Overvåkingslokaliteter i perioden 2020-2021. \\label{fig:lokaler_2020_2021}", out.width = "\\textwidth"----
ggplot(norge) +
  geom_sf() +
  geom_sf(data = nas_ins_lok,
          aes(fill = habitat_type,
              color = habitat_type),
          lwd = 2) +
  scale_color_nina(name = "Habitat",
                   palette = "orange-green") +
  scale_fill_nina(name = "Habitat",
                   palette = "orange-green")


## -------------------------------------------------------------------------------------------------------------------------
rand_500_ssb <- read_sf(con,
                        query = "
                        SELECT ano.*
                        FROM ano.ano_10000 ano
                        order by random() limit 500
                        "
                        )



## -------------------------------------------------------------------------------------------------------------------------
rand_500_ssb <- rand_500_ssb %>% 
  mutate(habitat_type = rep(c("Forest", "Semi-nat"), each = 250))



## -------------------------------------------------------------------------------------------------------------------------
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


## ----calc_ar5_shannon_index-----------------------------------------------------------------------------------------------
nas_ins_lok_1000m_agg_brick <- rast_list(nas_ins_lok_1000m_agg)
nas_ins_lok_1000m_arealtype_shannon_index <-lsm_l_shdi(nas_ins_lok_1000m_agg_brick)
nas_ins_lok_1000m_arealtype_shannon_index$layer <- names(nas_ins_lok_1000m_agg_brick)

locality_ar5_shdi <- nas_ins_lok_1000m_arealtype_shannon_index %>% 
  select(locality = layer,
         lsm_l_shdi = value) 



## ----calc_ar5_perimeter_length--------------------------------------------------------------------------------------------
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


## ----calc_ar5_area_sum----------------------------------------------------------------------------------------------------
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
  mutate_if(is.numeric , replace_na, replace = 0) %>% 
  mutate_if(is.numeric , units::drop_units)



## -------------------------------------------------------------------------------------------------------------------------
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
            "Hav"))

locality_ar5_area_group <- locality_ar5_area %>% 
  left_join(loc_types,
            by = c("locality" = "locality"),
            copy = T) %>% 
  ungroup() %>% 
  select(c("locality", "habitat_type")) %>% 
  mutate(habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type))




## ---- message= F, results = "hide"----------------------------------------------------------------------------------------
example_NMDS=metaMDS(locality_ar5_area_data, # Our community-by-species matrix
                     k=2)


## ----ar5_nmds, fig.pos = "!h", fig.cap = "NMDS over variasjonen i landskapskomposisjon i en sirkel med 1km radius rundt lokalitetene (arealtype i AR5). Siffrer i grønt viser semi-naturlige lokaliteter og siffrer i lilla viser lokaliteter i skog. \\label{fig:ar5_nmds}", out.width = "0.8\\textwidth"----
ordiplot(example_NMDS,type="n")

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


## ----semi_nat_ar5_plot, fig.cap = "Eksempel på variasjonen i landskapskomposition (arealtype i AR5) i semi-naturlig mark innen en 1 km buffer, SSB 500x500m grid i blått, felleplassering i rødt og de kartlagte ANO-sirklene i lilla. \\label{fig:semi_nat_ar5_plot}", fig.height = 9, fig.width = 5.5----
 
p1 <- locality_ar5_plot("Semi-nat_07",
                        legend = F,
                        strip = T)
p2 <- locality_ar5_plot("Semi-nat_22",
                        legend = F,
                        strip = T) 

p3 <- locality_ar5_plot("Semi-nat_20",
                        legend = F,
                        strip = T) 

p4 <- locality_ar5_plot("Semi-nat_03",
                        legend = F,
                        strip = T) 
 

leg <- locality_ar5_plot("Semi-nat_03",
                         legend = T,
                    strip = T) +
  theme(legend.position = "bottom",
       # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
       legend.key.size = unit(0.5, 'cm'),
       legend.text = element_text(size=7),
       legend.key.width = unit(0.5, 'cm'))

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




## ----skog_ar5_plot, fig.cap = "Eksempel på variasjonen i landskapskomposition (arealtype i AR5) i skogsmark innen en 1 km buffer, SSB 500x500m grid i blått, felleplassering i rødt og de kartlagte ANO-sirklene i lilla. \\label{fig:skog_ar5_plot}", fig.height = 9, fig.width = 5.5----

p1 <- locality_ar5_plot("Skog_20",
                        legend = F,
                        strip = T) 

p2 <- locality_ar5_plot("Skog_12",
                        legend = F,
                        strip = T) 
  

p3 <- locality_ar5_plot("Skog_14",
                        legend = F,
                        strip = T)


p4 <- locality_ar5_plot("Skog_05",
                        legend = F,
                        strip = T) 

leg <- locality_ar5_plot("Skog_05",
                         legend = T,
                    strip = T) +
  
  
  theme(legend.position = "bottom",
       # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
       legend.key.size = unit(0.5, 'cm'),
       legend.text = element_text(size=7),
       legend.key.width = unit(0.5, 'cm'))

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




## -------------------------------------------------------------------------------------------------------------------------
no_samplings_per_lok_q <- "

SELECT l.locality ,  
count(distinct ls.sampling_name),
min(ls.start_date) as  start_date,
max(ls.end_date) as end_date 
FROM events.locality_sampling ls,
events.year_locality yl,
locations.localities l,
events.sampling_trap st
WHERE ls.year_locality_id = yl.id
AND yl.locality_id = l.id
AND yl.year = 2021
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



## ----no_samplings_per_lok, results = "asis"-------------------------------------------------------------------------------
print(xtable(no_samplings_per_lok, caption = "Antall felletømminger per lokalitet i 2021. \\label{tab:no_samplings_per_lok}"))


## -------------------------------------------------------------------------------------------------------------------------
no_samples_q <-  "

SELECT trap_model as \"Felletype\", 
count(distinct st.sample_name) as \"Antall prosesserte felleprøver\"
FROM events.locality_sampling ls,
events.year_locality yl,
locations.traps,
events.sampling_trap st
WHERE ls.year_locality_id = yl.id
AND st.trap_id = traps.id
AND yl.year = 2021
AND yl.project_short_name = 'NasIns'
AND st.locality_sampling_id = ls.id
AND st.received_at_lab IS TRUE
GROUP BY trap_model
"

no_samples <- dbGetQuery(con,
                           no_samples_q) %>% 
  mutate(Felletype = ifelse(Felletype == "Malaise_Watkins_Black", "Malaisefelle", "Vindusfelle"))


## ---- antall_prover, results = "asis"-------------------------------------------------------------------------------------
print(xtable(no_samples, caption = "Antall prosesserte prøver i 2021. De 4 vindusfelleprøvene i hver tømmerunde på en plass ble slått sammmen og metastrekkodet sammen. \\label{tab:antall_prover}"))


## ---- fig.cap = "Flytskjema for identifisering av målgruppe sekvensene og kvalitetssikring av vår referansedatabase for kontaminasjon og referansesekvensmangel. \\label{fig:dna_flow_1}", out.width = "\\textwidth"----
include_graphics("figures/flyt_1.pdf")

## ---- fig.cap = "Flytskjema for klassifisering og kvalitetssikring av arter basert på DNA-metastrekkoding. Arter med lav arts-konfidens er ikke tatt med i den videre arbeidsflyten som avhenger av artsidentitet. \\label{fig:dna_flow_2}", out.width = "\\textwidth"----
include_graphics("figures/flyt_2.pdf")


## ---- fig.cap = "Flytskjema for horizon scanning for fremmede arter og nye arter for Norge. *I tillegg til å kryssjekke arter mot den norske fremmedartslista har vi også kryssjekket artsfunn i dette prosjektet mot databasene Daisie, Saul Pathways, Roy og List of Union Concern. Arter som er listet i en av disse databasene, er antatt å være fremmede arter i et eller flere land i Europa. Vi fant ingen arter lister i Roy et al. eller List of Union Concern i dette prosjektet. \\label{fig:dna_flow_3}", out.width = "0.9\\textwidth", fig.pos = "H"----
include_graphics("figures/flyt_3.pdf")


## -------------------------------------------------------------------------------------------------------------------------
tax_class_melt_detailed=readRDS("figures/overall_taxonomic_determinations.rds")
tax_class <- readRDS("figures/tax_class.rds")
fig_dat=readRDS("figures/figdata_species_level_determinations.rds")
custom_breaks=data.frame(values=c(0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.99),labels=c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
custom_breaks$transformed=car::logit(custom_breaks$values)
fig_dat$lower[fig_dat$lower<0]=0



## ----artsnivå_bionf, fig.pos = '!b', fig.cap = "Modellering av sannsynligheten for klassifisering til artsnivå (y-aksen) i forhold til hvor mange DNA-kopier vi finner fra arten (x-aksen). Til venstre på figurene er arter med veldig mange DNA-kopier (vanlige arter), mens til høyre er arter med få DNA-kopier (sjeldne arter). A) Totalt antall arter påvist i 2021, der x-aksen viser artsantall og B) fem utvalgte ordener som viser forskjellen mellom ulike grupper av insekter.\\label{fig:artsniva}"----

Afig <- ggplot(tax_class_melt_detailed,aes(x=rank_abundance,y=value,fill=variable)) + 
  geom_bar(stat="identity",position="fill") + 
  ylab("Andel arter") + 
  xlab("Arter i rekkefølge av abundans") + 
  scale_fill_nina() +
  #scale_fill_manual(values=rev(c("#008C95","#004F71","#E57200","#93328E","#7A9A01"))) + 
  theme(panel.background = element_blank(),
        axis.line=element_line(colour="black"),
        legend.title = element_blank(),
        legend.position="bottom") +
  guides(fill = guide_legend(override.aes = list(size=5), nrow = 2, byrow = T)) + 
   geom_vline(xintercept=min(tax_class$rank_abundance[tax_class$pct_all>0.899999999]),linetype="dashed") + 
  ggtitle("A")

Bfig <- ggplot(fig_dat[!fig_dat$set%in%"alle arter",],aes(x=scaled_x,y=y,colour=set)) + 
  geom_path(size=2.5) + 
  scale_color_nina() + 
  ylim(c(min(custom_breaks$transformed),max(custom_breaks$transformed))) + 
  scale_y_continuous(breaks=custom_breaks$transformed,labels=custom_breaks$labels) + 
  xlab("Arter i rekkefølge av abundans -->") + 
  theme(axis.text.x = element_blank(),
        axis.line=element_line(colour="black"),
        panel.background=element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        legend.position="bottom") + 
   guides(color = guide_legend(override.aes = list(size=5), nrow = 2, byrow = T)) + 
  ylab("Sannsynlighet for taksonomisk bestemmelse til artsnivå") + 
  ggtitle("B")

grid.arrange(Afig,Bfig,ncol=2)


## -------------------------------------------------------------------------------------------------------------------------
physeq.species=readRDS("figures/ord_physeq.rds")



## ----sekvenseringsdybde, fig.pos = '!h', fig.cap = "Sammenhengen mellom antall DNA-kopier (sekvenseringsdybde) per prøve på x-aksen og antall arter påvist i den samme prøven på y-aksen. \\label{fig:sekvenseringsdybde}", eval = T----

out.estimates.species=readRDS("figures/inext_results.rds")
tmp=left_join(data.frame(Site=physeq.species@sam_data$sample_name,seq_depth=sample_sums(physeq.species)),select(out.estimates.species$AsyEst[out.estimates.species$AsyEst$Diversity=="Species richness",],Site,Observed))
tmp$habitat[grepl("Semi-nat",tmp$Site)]="Semi-nat"
tmp$habitat[grepl("Skog",tmp$Site)]="Skog"
tmp=tmp[!is.na(tmp$habitat),]
tmp=tmp[grepl("MF",tmp$Site),]

ggplot(tmp,aes(x=seq_depth,y=Observed,col=habitat))+geom_point()+theme(panel.background = element_blank(),axis.line = element_line(color="black"))+xlab("Sekvenseringsdybde")+ylab("Antall arter")+scale_color_manual(values=c("#004F71","#008C95"),name="Habitat")




## -------------------------------------------------------------------------------------------------------------------------
col_summary <- readRDS("figures/col_morph_table.rds")  %>% 
  as_tibble()
lep_summary <- readRDS("figures/lep_morph_table.rds") %>% 
  as_tibble()

col_summary <- col_summary %>% 
  mutate(prøve = as.integer(gsub("(Prøve_)(.)", "\\2", prøve))) %>% 
  mutate(kategori = as.character(kategori)) %>% 
  mutate(kategori = ifelse(kategori == "mulig feil bestemt", "mulig feilbestemt", kategori)) %>% 
  mutate(kategori = str_to_sentence(kategori)) %>% 
  mutate(kategori = factor(kategori, levels = c("Ikke oppdaget", "Mulig feilbestemt", "Bestemt til slekt", "Bestemt til art")))

lep_summary <- lep_summary %>% 
  mutate(prøve = as.integer(gsub("(Prøve_)(.)", "\\2", prøve))) %>% 
  mutate(kategori = as.character(kategori)) %>% 
  mutate(kategori = ifelse(kategori == "mulig feil bestemt", "mulig feilbestemt", kategori)) %>% 
  mutate(kategori = str_to_sentence(kategori)) %>% 
  mutate(kategori = factor(kategori, levels = c("Ikke oppdaget", "Mulig feilbestemt", "Bestemt til slekt", "Bestemt til art")))

col = ggplot(col_summary,
             aes(x=method,
                 y=value,
                 fill=kategori)) +
  geom_bar(stat="identity") +
  facet_wrap(~prøve + felletype,
             ncol=15)+
  scale_fill_manual(values=c("grey85","#E57200","#008C95","#004F71"))+
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle=90),
        panel.background = element_blank(),
        axis.line = element_line(colour="black"),
        axis.title.x=element_blank(),
        strip.text=element_text(size=5.5))+
  scale_x_discrete(labels=c("DNA","Morf"))+
  ylab("Antall arter")+
  ggtitle("Coleoptera")

lep=ggplot(lep_summary,aes(x=method,y=value,fill=kategori))+
  geom_bar(stat="identity")+
  facet_wrap(~prøve,ncol=10)+
 scale_fill_manual(values=c("grey85","#E57200","#008C95","#004F71"))+
 theme(legend.title = element_blank(),
       axis.text.x = element_text(angle=90),
       panel.background = element_blank(),
       axis.line = element_line(colour="black"),
       axis.title.x=element_blank())+
  scale_x_discrete(labels=c("DNA","Morf"))+ylab("Antall arter")+
  ggtitle("Lepidoptera")


## ----dna_comp,fig.pos = "!hb", fig.cap = "Sammenligning av antall arter påvist for DNA-basert og morfologisk artsbestemmelse for biller (Coleoptera) og sommerfugler (Lepidoptera) fangt av malaise feller (M) og vindusfeller (V). \\label{fig:dna_comp}", fig.height = 5.5----
grid.arrange(col,
             lep,
             ncol = 1)


## -------------------------------------------------------------------------------------------------------------------------
df=readRDS("figures/df_table_for_morphology_DNA_curves.rds")

df.point <- df[which(df$method=="observed"),]
df.point$shape[grepl("Morph",df.point$site)]="Morfologi"
df.point$shape[grepl("DNA",df.point$site)]="DNA"
df.point$color=gsub("_Morph|_DNA","",df.point$site)
df.point$color=gsub("_DNA","",df.point$color)

df.line <- df[which(df$method!="observed"),]
df.line$method <- factor(df.line$method, 
                         c("interpolated", "extrapolated"),
                         c("interpolation", "extrapolation"))
df.line$color=gsub("_Morph|_DNA","",df.line$site)
df.line$shape[grepl("Morph",df.line$site)]="Morfologi"
df.line$shape[grepl("DNA",df.line$site)]="DNA"

df$shape[grepl("Morph",df$site)]="Morfologi"
df$shape[grepl("DNA",df$site)]="DNA"
df$color=gsub("_Morph|_DNA","",df$site)




## ----col_lepi_curves, fig.cap = "Artsakkumulasjonskurver for biller og sommerfugler fra malaise og vindusfeller ved bruk av morfologiske og DNA bestemmelser metoder. Notere at hver av vindusfellene egentlig er 4 feller slått sammen, slik at de tilsammen er 20 st. \\label{fig:col_lepi_curves}"----
ggplot(df, aes(x=x, y=y, colour=color,shape=shape,fill=color)) + 
  geom_point(aes(shape=shape,colour=color), size=5,data=df.point) +
  geom_line(aes(linetype=method), lwd=1.5, data=df.line) +
  geom_ribbon(aes(ymin=y.lwr, ymax=y.upr,
                  fill=color, colour=NULL), alpha=0.2,data=df) +
  labs(x="Antall feller", y="Antall arter") +
  theme(legend.position = "right", 
        legend.title=element_blank(),
        text=element_text(size=18),
        legend.box = "vertical") +
  scale_fill_manual(values=c("#004F71","#2DCCD3","#E57200"),
                    labels=c("Coleoptera malaisefeller","Coleoptera vindusfeller", "Lepidoptera malaisefeller")) +
  scale_color_manual(values=c("#004F71","#2DCCD3","#E57200"),
                     labels=c("Coleoptera malaisefeller","Coleoptera vindusfeller","Lepidoptera malaisefeller"))


## ----col_lepi_dna, fig.cap = "Korrelasjon mellom andel sommerfugler i forhold til antall biller på x-aksen og hvor stor andel av biller DNA-basert artsbestemmelse påviste i forhold til morfologisk artsbestemmelse \\label{fig:col_lepi_dna}", out.width = "0.8\\textwidth"----
include_graphics("figures/col_lepi_dna.png")


## -------------------------------------------------------------------------------------------------------------------------
ano_herb_agg <- tbl(con,
                    Id(schema = "views", 
                       table = "ano_herb_agg")) %>% 
  filter(project_short_name == "NasIns")



## ----ano_herb_spec_nor, eval = T, echo = F, fig.pos = "!b", fig.cap = "Dekningsgrad og artantall av karplanter i ANO-kartlegging i overvåkingsklokalitetene. For skogshabitat der en full ANO-kartlegging gjennomførtes viser vi gjennomsnittsverdier fra de 18 rutene som ble målt. Dekningsgrader på over 100 \\% er mulig fordi man registrerer dekning på flere høydelag i vegetasjonen.\\label{fig:ano_herb_spec_nor}", out.width = "0.7\\textwidth"----
ano_herb_agg %>% 
  ggplot(.) +
  geom_bar(aes(x = locality, y = ano_mean_no_spec, fill = ano_mean_cover), stat = "identity") +
  scale_fill_nina(name = "Middelverdi\ndekningsgrad\nkarpl.",
                  discrete = F,
                  palette = "darkblue-orange") +
  ylab("Middels antall arter karplanter per ANO sirkel") +
  xlab("Lokalitet") +
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust = 1,
                                   size = 6))


## -------------------------------------------------------------------------------------------------------------------------
year_locality_biodiv_ano <- year_locality_biodiv %>%
  left_join(ano_herb_agg,
            by = c("locality" = "locality"),
            copy = T)
  


## ----tot_insect_richn_herb_richn, echo = F, fig.pos = "!b", fig.cap = "Relasjon mellom artsantall insekter og artsantall karplanter i ANO-sirkel, fordelt på habitatstype. Insektdata kommer fra 1 malaisefelle i semi-naturlig mark og 1 malaisefelle og 4 vindusfeller i skog. \\label{fig:tot_insect_richn_herb_richn}", out.width = "0.7\\textwidth"----
year_locality_biodiv_ano %>% 
  group_by("habitat_type") %>% 
  ggplot(., 
         aes(x = ano_mean_no_spec,
             y = no_species)) +
  geom_point(aes(color = habitat_type)) +
  scale_color_nina(name = "Habitattype",
                   palette = "darkblue-orange") +
  ylab("Totalt artsantall av identifiserte insekter") +
  xlab("Middels antall arter karplanter per ANO sirkel") +
   geom_smooth(aes(color = habitat_type),
              method = "lm",
              se = T,
              lwd = 2)



## -------------------------------------------------------------------------------------------------------------------------
tree_agg <- tbl(con,
                Id(schema = "views",
                   table = "landsskog_agg"))


## -------------------------------------------------------------------------------------------------------------------------
year_locality_biomass_ano_tree <- tree_agg %>% 
  left_join(year_locality_biomass,
            by = c("locality" = "locality"),
            copy = T)


## ----skog_biomasse, fig.pos = "!b", fig.cap = "Sammenheng mellom målt alder på skogsbestand og gjennomsnittlig biomasse innsamlet insekter per felletømming (våtvekt i g). Data fra biomasse kommer kun fra 1 malaisefelle per tømmetilfelle per lokalitet.\\label{fig:skog_biomass}", out.width = "\\textwidth"----
year_locality_biomass_ano_tree %>% 
  group_by("habitat_type") %>% 
  ggplot(., 
         aes(x = avg_age,
             y = avg_wet_weight)) +
  geom_point(aes(color = dom_tree)) +
  scale_color_nina(name = "Dominerende treslag",
                   palette = "darkblue-orange") +
  ylab("Total biomasse (i malaisefelleer)") +
  xlab("Alder på største lokale trær") +
   geom_smooth(aes(color = dom_tree),
              method = "lm",
              se = T,
              lwd = 2)
  
# geom_smooth(aes(color = dom_tree),
#               method = "nls",
#               formula = y ~ a * x + b,
#               se = F,
#               method.args = list(start = list(a = 0.1,
#                                               b = 0.1)
#                                  ),
#               lwd = 2
#               )



## -------------------------------------------------------------------------------------------------------------------------
year_locality_biodiv_ano_tree <- tree_agg %>% 
  left_join(year_locality_biodiv_ano,
            by = c("locality" = "locality"),
            copy = T)


## ----skog_biodiv, fig.pos = "!t", fig.cap = "Sammenheng mellom målt alder på skogsbestand og totalt antall identifiserte arter. Insektdata kommer fra 1 malaisefelle og 4 vindusfeller i hver skogslokalitet.  \\label{fig:skog_biodiv}", out.width = "\\textwidth"----

year_locality_biodiv_ano_tree %>% 
  group_by(habitat_type)  %>% 
  ggplot(.,
         aes(x = avg_age,
                y = no_species)) +
  geom_point(aes( color = dom_tree),
             cex = 2) + 
     geom_smooth(aes(color = dom_tree),
              method = "lm",
              se = T,
              lwd = 2) +
    scale_color_nina(name = "Dominerende treslag",
                   palette = "darkblue-orange") +
  ylab("Totalt artsantall av identifiserte insekter") +
  xlab("Alder på største lokale trær")


## ---- figures-side, fig.show="hold", out.width="40%", fig.cap = "Temperatur og luftfuktighetsmåler av typen Hobo MX2301A (venstre) og temperatur og lysmåler av typen Hobo MX2202 (høyre). MX2301A er plassert i skygge og begge tømmes med blåtann til en skyløsning. \\label{fig:templogger}"----
include_graphics("figures/templogger.jpg")
include_graphics("figures/lyslogger_small.png")


## -------------------------------------------------------------------------------------------------------------------------
locality_sampling_loggerdata <- loggerdata_from_db(dataset = "NasIns",
                                                   agg_level = "locality_sampling") 



## ----logger_corr, fig.cap = "Korrelasjonsplot over loggerdata. Relativ luftfuktighet (%) måles av loggermodell MX2301A og lys (lux) av loggermodell MX2202. MX2301A er montert unner en lysskjerm og reagerer ikke like hyppig på sollys som en uskjermete loggeren MX2202 som også måler lys. \\label{fig:logger_corr}"----

locality_sampling_loggerdata %>% 
  select(Temp_MX2301A = avg_values_MX2301A_temperature,
         Temp_MX2202 = avg_values_MX2202_temperature,
         RL = avg_values_MX2301A_rh,
         Lys = avg_values_MX2202_light
         ) %>% 
  filter(Lys >= 100) %>% 
  chart.Correlation2(., 
                  histogram = TRUE,
                  method = "pearson",
                  color = ninaColors()[1])



## -------------------------------------------------------------------------------------------------------------------------
locality_sampling_weight_loggerdata <- biomass_mf_locality_sampling %>% 
  left_join(locality_sampling_loggerdata,
            by = c("sampling_name" = "sampling_name",
                   "year" = "year",
                   "locality" = "locality",
                   "habitat_type" = "habitat_type")) %>% 
  mutate(habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type))




## ----locality_sampling_biomass_loggerdata, fig.cap = "Sammenheng mellom temperatur og innsamlet biomasse per fangstperiode i malaisefelleer. \\label{fig:locality_sampling_biomass_loggerdata}", out.width = "0.8\\textwidth"----
locality_sampling_weight_loggerdata %>% 
  group_by(habitat_type) %>% 
  ggplot(aes(y = avg_wet_weight/no_trap_days, 
             x = avg_values_MX2301A_temperature)) +
  geom_point(aes(color = habitat_type)) +
  scale_color_nina(name = "Habitatstype",
                   palette = "darkblue-orange") +
     geom_smooth(aes(color = habitat_type),
              method = "lm",
              se = T,
              lwd = 2) +
  ylab("Avrunnen våtvekt (g.) per felledag") +
  xlab("Middeltemperatur ved felle unner fangstperioden")


## -------------------------------------------------------------------------------------------------------------------------
locality_sampling_biodiv_loggerdata <- locality_sampling_biodiv %>% 
  left_join(locality_sampling_loggerdata,
            by = c("sampling_name" = "sampling_name",
                   "year" = "year",
                   "locality" = "locality",
                   "habitat_type" = "habitat_type")) %>% 
  mutate(habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type))



## ----locality_sampling_biodiv_loggerdata, fig.cap = "Sammenheng mellom temperatur og mengde arter innsamlet per fangstperiode i malaisefelleer. \\label{fig:locality_sampling_biodiv_loggerdata}", out.width = "0.8\\textwidth"----
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


## -------------------------------------------------------------------------------------------------------------------------
#Sum the percentage number of reads within samples to individual species_latin names 
rel_no_reads_q <- "  
SELECT species_latin_fixed, 
id_order,
id_family,
sum(perc_reads) sum_perc_reads
FROM
(SELECT foo.*, 
tot_reads.tot_reads,
no_reads/tot_reads.tot_reads as perc_reads
FROM 
(SELECT ids.sampling_trap_id, 
species_latin_fixed,
id_order,
id_family,
no_reads
FROM occurrences.observations obs,
events.identifications ids,
events.sampling_trap st,
locations.traps 
WHERE obs.identification_id = ids.id
AND st.id = ids.sampling_trap_id
AND st.trap_id = traps.id
AND trap_short_name LIKE 'MF%' 
) foo LEFT JOIN
(SELECT ids.sampling_trap_id, 
sum(no_reads) tot_reads
FROM occurrences.observations obs,
events.identifications ids,
events.sampling_trap st,
locations.traps 
WHERE obs.identification_id = ids.id
AND st.id = ids.sampling_trap_id
AND st.trap_id = traps.id
AND trap_short_name LIKE 'MF%'
GROUP BY ids.sampling_trap_id) tot_reads
ON foo.sampling_trap_id = tot_reads.sampling_trap_id
 ) bar
GROUP BY species_latin_fixed, id_order, id_family
ORDER BY sum_perc_reads DESC

"
rel_no_reads_mf <- dbGetQuery(con,
                           rel_no_reads_q ) %>% 
  mutate(tot_sum_perc_reads = round(sum_perc_reads / sum(sum_perc_reads) * 100, 4))


## -------------------------------------------------------------------------------------------------------------------------
times_most_reads_q <- "


SELECT id_order,
id_family,
species_latin_fixed, count(*) times_abund
FROM 
(SELECT DISTINCT ON (ids.sampling_trap_id)
ids.sampling_trap_id,
species_latin_fixed,
id_order, 
id_family,
no_reads
FROM occurrences.observations obs,
events.identifications ids,
events.sampling_trap st,
locations.traps 
WHERE obs.identification_id = ids.id
AND st.id = ids.sampling_trap_id
AND st.trap_id = traps.id
AND trap_short_name LIKE 'MF%'
ORDER BY sampling_trap_id, no_reads DESC) foo
GROUP BY species_latin_fixed, id_order, id_family
ORDER BY times_abund DESC



"

times_most_reads_mf <- dbGetQuery(con,
                               times_most_reads_q)


## ----most_often_dom, fig.cap = "De arter som oftest dominerte Malaise-felleprøvene i mengde DNA-reads. Kun de arter som dominerte prøvene 4 eller fler ganger vises. NB at Sarcophaga lehmannni er kjent fra Norge tidligere (Knut Rognes), men ikke registrert i Artskart. Arter merket med sp.x er arter som ikke har match i referansebiblioteket men der man kan knytte dem til en høyere taksonomisk nivå. \\label{fig:most_often_dom}", fig.width = 7, fig.height = 4----
times_most_reads_mf %>% 
   filter(times_abund >= 4) %>% 
  ggplot(., aes(x = reorder(species_latin_fixed, -times_abund),
                y = times_abund)) +
  geom_bar(stat = "identity",
           aes(fill = id_order)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Art-id") +
  ylab("Antall felleprøver arten\nhade flest DNA-reads") +
  scale_fill_nina(name = "Orden")


## -------------------------------------------------------------------------------------------------------------------------
#Sum the percentage number of reads within samples to individual species_latin names 
rel_no_reads_vf_q <- "  
SELECT species_latin_fixed, 
id_order,
id_family,
sum(perc_reads) sum_perc_reads
FROM
(SELECT foo.*, 
tot_reads.tot_reads,
no_reads/tot_reads.tot_reads as perc_reads
FROM 
(SELECT ids.sampling_trap_id, 
species_latin_fixed,
id_order,
id_family,
no_reads
FROM occurrences.observations obs,
events.identifications ids,
events.sampling_trap st,
locations.traps 
WHERE obs.identification_id = ids.id
AND st.id = ids.sampling_trap_id
AND st.trap_id = traps.id
AND trap_short_name LIKE 'VF%' 
) foo LEFT JOIN
(SELECT ids.sampling_trap_id, 
sum(no_reads) tot_reads
FROM occurrences.observations obs,
events.identifications ids,
events.sampling_trap st,
locations.traps 
WHERE obs.identification_id = ids.id
AND st.id = ids.sampling_trap_id
AND st.trap_id = traps.id
AND trap_short_name LIKE 'VF%'
GROUP BY ids.sampling_trap_id) tot_reads
ON foo.sampling_trap_id = tot_reads.sampling_trap_id
 ) bar
GROUP BY species_latin_fixed, id_order, id_family
ORDER BY sum_perc_reads DESC

"
rel_no_reads_vf <- dbGetQuery(con,
                           rel_no_reads_vf_q ) %>% 
  mutate(tot_sum_perc_reads = round(sum_perc_reads / sum(sum_perc_reads)*100, 4))


## -------------------------------------------------------------------------------------------------------------------------
times_most_reads_vf_q <- "


SELECT id_order,
id_family,
species_latin_fixed, count(*) times_abund
FROM 
(SELECT DISTINCT ON (ids.sampling_trap_id)
ids.sampling_trap_id,
species_latin_fixed,
id_order, 
id_family,
no_reads
FROM occurrences.observations obs,
events.identifications ids,
events.sampling_trap st,
locations.traps 
WHERE obs.identification_id = ids.id
AND st.id = ids.sampling_trap_id
AND st.trap_id = traps.id
AND trap_short_name LIKE 'VF%'
ORDER BY sampling_trap_id, no_reads DESC) foo
GROUP BY species_latin_fixed, id_order, id_family
ORDER BY times_abund DESC

"

times_most_reads_vf <- dbGetQuery(con,
                               times_most_reads_vf_q)


## ----most_often_dom_vf, fig.cap = "De arter som oftest dominerte vindusfelleprøvene i mengde DNA-reads. Kun de arter som dominerte prøvene 2 eller fler ganger vises. Notere at Selatosomus impressus er synonym for Paraphotistus impressus, noe som tyder på at denne art er registert unner begge navn i referansebibliotekene. \\label{fig:most_often_dom_vf}", fig.width = 7, fig.height = 4----
times_most_reads_vf %>%
   filter(times_abund >= 2) %>% 
  ggplot(., aes(x = reorder(species_latin_fixed, -times_abund),
                y = times_abund)) +
  geom_bar(stat = "identity",
           aes(fill = id_order)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Art-id") +
  ylab("Antall felleprøver arten\nhade flest DNA-reads") +
  scale_fill_nina(name = "Orden")


## ----region_hab_div-------------------------------------------------------------------------------------------------------
region_habitat_div <- obs_from_db(agg_level = "region_habitat_year") %>% 
  mutate(habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type)) %>% 
  select(`År` = year,
         `Habitat-type` = habitat_type,
         `Region` = region_name,
         `Antall arter` = no_species)

total_div <- obs_from_db(agg_level = "total",
                         dataset = c("NasIns")) %>% 
  mutate(År = "Alle",
          `Habitat-type` = "Alle",
          Region = "Alle" ) %>% 
  select(`År`,
         `Habitat-type`,
         `Region`,
         `Antall arter` = no_species)

region_habitat_div <- region_habitat_div %>% 
  rbind(total_div)

save(region_habitat_div,
     file = "out/region_habitat_div.rda")


## ----div_tab, results = "asis"--------------------------------------------------------------------------------------------
xtable(region_habitat_div, caption = "Totalt antall identifiserte taksa fra overvåkingen i 2020-2021. \\label{tab:div_tab}")


## ----richness_locality, fig.cap = "Totalt observert artsantall per lokalitet. Notere at metastrekkodingsrutinen skiller seg mellom år. \\label{fig:richness_locality}", out.width = "0.9\\textwidth"----
year_locality_biodiv %>% 
  mutate(habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type),
         year = as.factor(year)) %>% 
  ggplot(.) +
  geom_bar(aes(x = locality,
               y = no_species,
               fill = year),
           stat = "identity"
           ) +
  scale_fill_nina(name = "År",
                  palette = "darkblue-orange") +
  ylab("Totalt antall observerte arter") +
    theme(axis.text.x = element_text(angle = 90, 
                                     vjust = 0.5, 
                                     hjust = 1,
                                     size = 8))


## -------------------------------------------------------------------------------------------------------------------------
semi_nat_ost <- community_matrix_from_db(trap_type = "MF",
                                         subset_habitat = "Semi-nat",
                                         subset_region = "Østlandet",
                                         transposed_matrix = T) 

skog_ost <- community_matrix_from_db(trap_type = "MF",
                                         subset_habitat = "Forest",
                                         subset_region = "Østlandet",
                                         transposed_matrix = T) 

semi_nat_trond <- community_matrix_from_db(trap_type = "MF",
                                         subset_habitat = "Semi-nat",
                                         subset_region = "Trøndelag",
                                         transposed_matrix = T) 




ost_iNEXT_data <- toiNEXT(input = list(semi_nat_ost,
                                       skog_ost, 
                                       semi_nat_trond))

names(ost_iNEXT_data) <- c("Semi-nat_Ost",
                           "Skog_Ost",
                           "Semi-nat_Trond")


#ost_iNEXT <- iNEXT(ost_iNEXT_data,
#                  datatype = "incidence_raw",
#                  endpoint = 50)

#save(ost_iNEXT,
#    file = "out/ost_iNEXT.rda")

load(file = "out/ost_iNEXT.rda")



## ----div_est, fig.cap = "Oppmålt og estimert total artsantall (Chao) i malaisefelleer i de undersøkte regionene. Notere  at extrapoleringen antar en perfekt deteksjonsevne.\\label{fig:div_est}", out.width = "0.9\\textwidth"----
ggiNEXT(ost_iNEXT) +
  scale_color_nina() +
  ylab("Artsantall") +
  xlab("Antall lokaliteter") +
   theme(legend.text = element_text(size = 10))



## -------------------------------------------------------------------------------------------------------------------------
all_comm_mat <- community_matrix_from_db(trap_type = "All",
                                         dataset = "NasIns",
                                         as_tibble = T)



## -------------------------------------------------------------------------------------------------------------------------
all_comm_mat_mat <- all_comm_mat %>% 
  select(-c(year, locality))

levels <- tibble(l1 = seq(1:nrow(all_comm_mat_mat)),
                 l2 = rep(1, nrow(all_comm_mat_mat))
)


## ---- eval = F------------------------------------------------------------------------------------------------------------
## div_part <- multipart(y = all_comm_mat_mat,
##           x = levels,
##           global = T,
##           nsimul = 999,
##           scales = 0)
## 
## save(div_part,
##      file = "out/div_part.rds")


## -------------------------------------------------------------------------------------------------------------------------
load(file = "out/div_part.rds")


## ----div_part_tab, results = 'asis'---------------------------------------------------------------------------------------
div_part_tab <- tibble(Nivå = c("$\\alpha$", 
                                "$\\beta$",
                                "$\\gamma$"),
                       "Målt verdi" = div_part$oecosimu$statistic[c(1, 3, 2)],
                       "Stand. Effekt" = div_part$oecosimu$z[c(1, 3, 2)],
                       "Forventet verdi" = div_part$oecosimu$means[c(1, 3, 2)],
                       "Prob (sim)" = div_part$oecosimu$pval[c(1, 3, 2)])

print(xtable(div_part_tab, caption = "Randomiseringstest av alpha, beta og gamma-komponenter av diversitet (N = 999).  \\label{tab:div_part_tab}",
             digits = 3),
       sanitize.text.function = function(x) as.character(x),
       sanitise.rownames.function = function(x) as.character(x),
       include.rownames = F)


## -------------------------------------------------------------------------------------------------------------------------
div_betacore <- betapart.core(all_comm_mat_mat)


## -------------------------------------------------------------------------------------------------------------------------
all_beta_part_freq <- beta.multi(div_betacore) %>% 
  unlist()

all_beta_part_freq <- as.data.frame(t(round(all_beta_part_freq / all_beta_part_freq["beta.SOR"] * 100, 2)))
names(all_beta_part_freq) <- c("$\\beta_{SIM}\\%$", "$\\beta_{NES}\\%$", "$\\beta_{SOR}\\%$")

  

## ----beta_tab, results = "asis"-------------------------------------------------------------------------------------------
print(xtable(all_beta_part_freq, caption = "Komponenter i oppmålt beta-diversitet mellom lokalitetetene i overvåkingen. $\\beta_{SIM}\\%$ står for 'species turnover', $\\beta_{NES}\\%$ for 'nestedness' og $\\beta_{SOR}\\%$ for den totale $\\beta-$diversiteten.\\label{tab:beta_tab}"), sanitize.colnames.function = function(x) as.character(x))


## -------------------------------------------------------------------------------------------------------------------------
skog_ost_comm_mat <- community_matrix_from_db(trap_type = "MF",
                                              dataset = "NasIns",
                                              subset_habitat = "Forest",
                                              subset_region = "Østlandet",
                                              as_tibble = T,
                                              transposed_matrix = F
                                              ) %>% 
  select(-c(year, locality))

semi_nat_comm_mat <- community_matrix_from_db(trap_type = "MF",
                                              dataset = "NasIns",
                                              subset_habitat = "Semi-nat",
                                              as_tibble = T,
                                              transposed_matrix = F
                                              ) %>% 
  select(-c(year, locality))





## -------------------------------------------------------------------------------------------------------------------------
skog_ost_dist_beta <- combine_dist_to_comm_mat(comm_mat = skog_ost_comm_mat,
                                               region_name = "('Østlandet')",
                                               habitat_type = "Forest")


## ---- eval = F------------------------------------------------------------------------------------------------------------
## beta_sne_lm <- glm(beta_sne~distance, data = skog_ost_dist_beta)
## summary(beta_sne_lm)
## beta_sim_lm <- glm(beta_sim~distance, data = skog_ost_dist_beta)
## summary(beta_sim_lm)
## beta_sor_lm <- glm(beta_sor~distance, data = skog_ost_dist_beta)
## summary(beta_sor_lm)
## 
## 
## with(summary(beta_sor_lm), 1 - deviance/null.deviance)


## ----skog_ost_dist_beta, fig.cap = "Oppdeling av $\\beta-$diversitet i turnover og nestedness, og dess sammenheng med avstand mellom lokaliteter. Skogslokaliteter på Østlandet. \\label{fig:skog_ost_dist_beta}", fig.height = 9----
plot_beta_part(skog_ost_dist_beta)


## -------------------------------------------------------------------------------------------------------------------------
semi_nat_dist_beta <- combine_dist_to_comm_mat(comm_mat = semi_nat_comm_mat,
                                               region_name = "('Østlandet', 'Trøndelag')",
                                               habitat_type = "Semi-nat")


## ----semi_nat_ost_dist_beta, fig.cap = " Oppdeling av $\\beta-$diversitet i turnover og nestedness, og dess sammenheng med avstand mellom lokaliteter. Semi-naturlige lokaliteter på Østlandet og Trøndelag. \\label{fig:semi_nat_ost_dist_beta}", fig.height = 9----
plot_beta_part(semi_nat_dist_beta)


## ---- eval = F------------------------------------------------------------------------------------------------------------
## 
## skog_ost_comm_mat_lepi <- community_matrix_from_db(trap_type = "MF",
##                                               dataset = "NasIns",
##                                               subset_habitat = "Forest",
##                                               subset_region = "Østlandet",
##                                               subset_orders = "Lepidoptera",
##                                               as_tibble = T,
##                                               transposed_matrix = F,
##                                               exclude_singletons = T
##                                               ) %>%
##   select(-c(year, locality))
## 
## semi_nat_ost_comm_mat_lepi <- community_matrix_from_db(trap_type = "MF",
##                                               dataset = "NasIns",
##                                               subset_habitat = "Semi-nat",
##                                               subset_region = "Østlandet",
##                                               subset_orders = "Lepidoptera",
##                                               as_tibble = T,
##                                               transposed_matrix = F,
##                                               exclude_singletons = T
##                                               ) %>%
##   select(-c(year, locality))
## 
## 


## ---- eval = F------------------------------------------------------------------------------------------------------------
## skog_ost_dist_beta_lepi <- combine_dist_to_comm_mat(comm_mat = skog_ost_comm_mat_lepi,
##                                                region_name = "Østlandet",
##                                                habitat_type = "Forest")
## 
## semi_nat_ost_dist_beta_lepi <- combine_dist_to_comm_mat(comm_mat = semi_nat_ost_comm_mat_lepi,
##                                                region_name = "Østlandet",
##                                                habitat_type = "Semi-nat")


## ----skog_ost_dist_beta_lepi, fig.cap = "Skogslokaliteter på Østlandet. Oppdeling av $\\beta-$diversitet i turnover og nestedness i Lepidoptera, og dess sammenheng med avstand mellom lokaliteter. Alle arter som er blitt funnet kun én gang i noen lokalitet exkludert. \\label{fig:skog_ost_dist_beta_lepi}", fig.height = 9, eval = F----
## plot_beta_part(skog_ost_dist_beta_lepi)


## ----semi_nat_ost_dist_beta_lepi, fig.cap = "Semi-naturlige lokaliteter på Østlandet. Oppdeling av $\\beta-$diversitet i turnover og nestedness hos Lepidoptera, og dess sammenheng med avstand mellom lokaliteter. Alle arter som er blitt funnet kun én gang i noen lokalitet exkludert. \\label{fig:semi_nat_ost_dist_beta_lepi}", fig.height = 9, eval = F----
## plot_beta_part(semi_nat_ost_dist_beta_lepi)


## ----haplotype_kart, fig.cap = "Frekvenser av haplotyper (varianter av ASV/DNA-kopier) for ulike lokaliteter som visere genetisk variasjon innen hver enkeltart. Her er de tre artene med flest DNA-kopier innenfor hver av ordenene sommerfugler (Lepidoptera), biller (Coleoptera), tovinger (Diptera) og vepser (Hymenoptera) valgt ut som eksempler. \\label{fig:haplotype_kart}", eval = T, out.width = "0.9\\textwidth"----
include_graphics("figures/composite_map_figure.pdf")


## ----raw_asv_haplotype, fig.cap = "Boksplot som viser avvik i haplotypediversitet (ASV-antall) i forhold til forventet antall haplotyper gitt antallet arter påvist i en lokalitet. Alle innsamlinger på ulike datoer på samme lokalitet er inkludert og står for variasjonen innen hver boks. Den vertikale linjen på x-aksen angir forventet haplotypediversitet på tvers av alle lokalitetene. Bokser til venstre for denne linjen viser altså lokaliteter som har lavere haplotypediversitet enn forventet ut ifra antall arter.\\label{fig:raw_asv_haplotype}", out.width = "0.7\\textwidth"----

sp.asv.rich = readRDS("figures/asvs_given_species_data.rds")

sp.asv.rich <- sp.asv.rich %>% 
  mutate(sig_colors = ifelse(sig_colors == "significant", "Signifikant", "Ikke sign."))


ggplot(sp.asv.rich, 
       aes(x=lokalitet,
           y=Estimator_haplo-predicted_Estimator_haplo,
           fill=sig_colors)
       ) + 
  geom_hline(yintercept=0)+
  geom_boxplot()+
  coord_flip()+
  scale_fill_manual(values=c("#A2AAAD","#008C95"))+
  theme(panel.background=element_blank(),
        axis.line=element_line(colour="black"),
        axis.title.y = element_blank(),
        legend.title=element_blank())+
  ylab("Residualer")


## ----shan_haplotype, fig.cap = "Boksplot som viser avvik i Shannon indeks vektet på haplotypediversitet (ASV-antall) i forhold til forventet Shannon indeks gitt antallet arter påvist i en lokalitet. Alle innsamlinger på ulike datoer på samme lokalitet er inkludert og står for variasjonen innen hver boks. Den vertikale linjen på x-aksen angir forventet Shannon diversitet på tvers av alle lokalitetene. Bokser til venstre for denne linjen viser altså lokaliteter som har lavere Shannon indeks enn forventet ut ifra antall arter.  \\label{fig:shan_haplotype}", out.width = "0.7\\textwidth"----

sp.gen.shan=readRDS("figures/shan_gen_given_shan_species_data.rds")

sp.gen.shan <- sp.gen.shan %>% 
   mutate(sig_colors = ifelse(sig_colors == "significant", "Signifikant", "Ikke sign."))

ggplot(sp.gen.shan, aes(x=lokalitet,
                        y=Estimator_gen-predicted_Estimator_gen,
                        fill=sig_colors))+ 
  geom_hline(yintercept=0)+
  geom_boxplot()+
  coord_flip()+
  scale_fill_manual(values=c("#A2AAAD","#008C95"))+
  theme(panel.background=element_blank(),axis.line=element_line(colour="black"),axis.title.y = element_blank(),legend.title=element_blank())+
  ylab("Residuals")


## -------------------------------------------------------------------------------------------------------------------------
biomass_mf_locality_sampling_time <- biomass_mf_locality_sampling %>% 
  left_join(locality_sampling,
            by = c("sampling_name" = "sampling_name"),
            copy = T) %>% 
  mutate(julian_day = lubridate::yday(end_date),
         habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type))  %>% 
  filter(sum_wet_weight > 0)
  




## ----biomass_mf_julian_day_nor, fig.pos = "!h" , fig.cap = "Tidstrend over innsamlet biomasse (i malaisefeller). \\label{fig:biomass_mf_julian_day_nor}", out.width = "\\textwidth"----
biomass_mf_locality_sampling_time %>% 
ggplot(.,
       aes(y = sum_wet_weight/no_trap_days,
           x = as.Date(julian_day, origin = as.Date("2010-01-01")),
           group = interaction(habitat_type, 
                               region_name,
                               sep = "_"),
           color = interaction(habitat_type, 
                               region_name,
                               sep = "_"))) +
  geom_smooth(formula = y ~ x,
             cex = 2) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_nina(name = "Tidsserie") +
  facet_wrap(~year,
             ncol = 1) +
  ylab("Biomasse per felledøgn (gram)") + 
  xlab("Tømmedato") +
  scale_x_date(date_labels = "%d.%b")
  


## -------------------------------------------------------------------------------------------------------------------------
biomass_mf_year_locality <- weights_from_db(trap_type = "MF",
                                            agg_level = "year_locality",
                                            subset_region = "Østlandet")



## ----biomass_year_hab, fig.pos = "!h", fig.cap = "Middelverdi av (middelverdi) biomasse insekter i malaisefeller per lokalitet og år i skog og semi-naturlig mark. Kun data fra Østlandet er tatt med. Notere  at data fra de to årene er innsamlet unner ulike tidsperioder og derfor ikke er helt sammenlignbare.\\label{fig:biomass_year_hab}", out.width = "0.7\\textwidth"----
year_mean_biomass <- biomass_mf_year_locality %>% 
  mutate(habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type)) %>% 
  group_by(year,
           habitat_type) %>% 
  summarise(mean_biomass = mean(avg_wet_weight),
            sd_biomass = sd(avg_wet_weight)) %>% 
  mutate(year = as.factor(year),
         habitat_type = as.factor(habitat_type)) %>% 
  ungroup()

  ggplot(year_mean_biomass,
         aes(x = year,
             y = mean_biomass,
             group = habitat_type)) +
    geom_point(aes(x = year,
               y = mean_biomass,
               color = habitat_type),
               cex = 3) +
    geom_errorbar(aes(x = year,
                      ymin = mean_biomass - sd_biomass,
                      ymax = mean_biomass + sd_biomass),
                   width=.1,
                  position=position_dodge(0.05)) +
    geom_line(aes(x = year,
                  y = mean_biomass,
                  group = habitat_type)) +
    scale_color_nina(name = "Habitatstype",
                     palette = "darkblue-orange") +
    ylab("Middelverdi biomasse per sampletilfelle og lokalitet") +
    xlab("År")



## ----skog_ost_biomass, fig.cap = "Variasjon i biomasse over tid for skogslokaliteter på Østlandet. \\label{fig:skog_ost_biomass}", fig.height = 8----

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
  ylab("Biomasse per\nfelledøgn (gram))") + 
  scale_x_date(date_labels = "%d.%m.%y",
               limits = as.Date(xlims, origin = as.Date("2020-01-01"))) +
  scale_y_continuous(limits = ylims) +
  xlab("Tømmedato") +
  scale_color_nina(name = "Lokalitet") +
theme(legend.position = "bottom",
       # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
       legend.key.size = unit(0.5, 'cm'),
       legend.text = element_text(size=7),
       legend.key.width = unit(0.5, 'cm')) 



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
       legend.key.width = unit(0.5, 'cm'))

marrangeGrob(list(p1,
             p2),
             top = "",
             ncol = 1,
             nrow = 2)



## ----semi_nat_ost_biomass,fig.pos = "!b",  fig.cap = "Variasjon i biomasse over tid for semi-naturlig mark på Østlandet. \\label{fig:semi_nat_ost_biomass}", fig.height = 8----

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
       legend.key.width = unit(0.5, 'cm')) 



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
       legend.key.width = unit(0.5, 'cm'))

marrangeGrob(list(p1,
             p2),
             top = "",
             ncol = 1,
             nrow = 2)



## ----semi_nat_trond_biomass, fig.cap = "Variasjon i biomasse over tid for semi-naturlig mark i Trøndelag. \\label{fig:semi_nat_trond_biomass}", fig.height = 4----
biomass_mf_locality_sampling_time %>% 
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
       legend.key.width = unit(0.5, 'cm')) 



## -------------------------------------------------------------------------------------------------------------------------
locality_sampling_biodiv_sampling_time <- locality_sampling_biodiv %>% 
    left_join(locality_sampling,
            by = c("sampling_name" = "sampling_name"),
            copy = T) %>% 
  mutate(julian_day = lubridate::yday(end_date),
         habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type)) 
  


## ----diversity_mf_julian_day_nor, fig.pos = "!h", fig.cap = "Tidstrend over artsantall (i malaisefeller). \\label{fig:diversity_mf_julian_day_nor}", out.width = "\\textwidth"----
locality_sampling_biodiv_sampling_time %>% 
ggplot(.,
       aes(y = no_species/no_trap_days,
           x = as.Date(julian_day, origin = as.Date("2010-01-01")),
           group = interaction(habitat_type, 
                               region_name,
                               sep = "_"),
           color = interaction(habitat_type, 
                               region_name,
                               sep = "_"))) +
  geom_smooth(formula = y ~ x,
             cex = 2) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_nina(name = "Tidsserie") +
  facet_wrap(~year,
             ncol = 1) +
  ylab("Artsantall per felledøgn") + 
  xlab("Tømmedato") +
  scale_x_date(date_labels = "%d.%b")
  


## -------------------------------------------------------------------------------------------------------------------------
year_locality_div <- obs_from_db(subset_region = "Østlandet",
                                 agg_level = "year_locality",
                                 trap_type = "MF") %>% 
  mutate(habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type))



## ----div_year_hab, fig.pos = "!h", fig.cap = "Middelverdi av antall arter fanget per lokalitet og år i skog og semi-naturlig mark. Kun malaisefeller. Oppgangen i 2021 forklares med en overgang av metastrekkoding til NovaSeq-maskin \\label{fig:div_year_hab}", out.width = "0.7\\textwidth"----
year_mean_div <- year_locality_div %>% 
  group_by(year,
           habitat_type) %>% 
  summarise(mean_richness = mean(no_species),
            sd_richness = sd(no_species)) %>% 
  mutate(year = as.factor(year),
         habitat_type = as.factor(habitat_type)) %>% 
  ungroup()

  ggplot(year_mean_div,
         aes(x = year,
             y = mean_richness,
             group = habitat_type)) +
    geom_point(aes(x = year,
               y = mean_richness,
               color = habitat_type),
               cex = 3) +
    geom_errorbar(aes(x = year,
                      ymin = mean_richness - sd_richness,
                      ymax = mean_richness + sd_richness),
                   width=.1,
                  position=position_dodge(0.05)) +
    geom_line(aes(x = year,
                  y = mean_richness,
                  group = habitat_type)) +
    scale_color_nina(name = "Habitatstype",
                     palette = "darkblue-orange") +
    ylab("Antall arter per lokalitet (middelv.)") +
    xlab("År")



## ----skog_ost_div, fig.pos = "!h", fig.cap = "Variasjon i artsantall over tid for skogslokaliteter på Østlandet. \\label{fig:skog_ost_div}", fig.height = 8----

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
       legend.key.width = unit(0.5, 'cm')) 



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
       legend.key.width = unit(0.5, 'cm'))

marrangeGrob(list(p1,
             p2),
             top = "",
             ncol = 1,
             nrow = 2)



## ----skog_semi_nat_div, fig.cap = "Variasjon i artsantall over tid for semi-naturlig mark på Østlandet. \\label{fig:skog_semi_nat_div}", fig.height = 8----
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
       legend.key.width = unit(0.5, 'cm')) 



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
       legend.key.width = unit(0.5, 'cm'))

marrangeGrob(list(p1,
             p2),
             top = "",
             ncol = 1,
             nrow = 2)


## ----semi_nat_trond_div, fig.cap = "Variasjon i artsantall over tid for semi-naturlig mark i Trøndelag. \\label{fig:semi_nat_trond_div}", fig.height = 4----
locality_sampling_biodiv_sampling_time %>% 
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
       legend.key.width = unit(0.5, 'cm'))


## -------------------------------------------------------------------------------------------------------------------------
pollinator_table <- tibble("Norsk navn" = c("Gravebier",
                                    "Langtunbebier",
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



## ----pollinator_table, results = "asis"-----------------------------------------------------------------------------------
print(xtable(pollinator_table, caption = "Familier som brukes for å ta ut pollinatorer i datasettet. Notere at vi bare tar ut en familie av tovinger, tross att mange flere arter er tenkbare pollinatorer, men at der deres funksjon ikke er homogen innen enkelte familier. Riodinidae har idag ikke noen kjent forekomst i Norge. \\label{tab:pollinator_table}"))



## -------------------------------------------------------------------------------------------------------------------------
locality_sampling_pollinators <- obs_from_db(subset_families = pollinators,
                                             agg_level = "locality_sampling",
                                             trap_type = "MF")


## -------------------------------------------------------------------------------------------------------------------------
locality_sampling_pollinators_sampling_time <- locality_sampling_pollinators %>% 
    left_join(locality_sampling,
            by = c("sampling_name" = "sampling_name"),
            copy = T) %>% 
  mutate(julian_day = lubridate::yday(end_date),
         habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type)) 
  


## ----pollinators_mf_julian_day_nor, fig.cap = "Tidstrend over artsantall pollinatorer (i malaisefeller). \\label{fig:pollinators_mf_julian_day_nor}", out.width = "0.8\\textwidth"----
locality_sampling_pollinators_sampling_time %>% 
ggplot(.,
       aes(y = no_species/no_trap_days,
           x = as.Date(julian_day, origin = as.Date("2010-01-01")),
           group = interaction(habitat_type, 
                               region_name,
                               sep = "_"),
           color = interaction(habitat_type, 
                               region_name,
                               sep = "_"))) +
  geom_smooth(formula = y ~ x,
             cex = 2) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_nina(name = "Tidsserie") +
  facet_wrap(~year,
             ncol = 1) +
  ylab("Artsantall per felledøgn") + 
  xlab("Tømmedato") +
  scale_x_date(date_labels = "%d.%b")
  


## -------------------------------------------------------------------------------------------------------------------------
year_locality_pollinators <- obs_from_db(subset_families = pollinators,
                                         subset_region = "Østlandet",
                                         agg_level = "year_locality",
                                         trap_type = "MF") %>% 
  mutate(habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type))



## ----poll_year_hab, fig.cap = "Middelverdi av antall pollinatorer fanget per lokalitet og år i skog og semi-naturlig mark. Kun malaisefeller. \\label{fig:poll_year_hab}", out.width = "0.7\\textwidth"----
year_mean_pollinators <- year_locality_pollinators %>% 
  group_by(year,
           habitat_type) %>% 
  summarise(mean_richness = mean(no_species),
            sd_richness = sd(no_species)) %>% 
  mutate(year = as.factor(year),
         habitat_type = as.factor(habitat_type)) %>% 
  ungroup()

  ggplot(year_mean_pollinators,
         aes(x = year,
             y = mean_richness,
             group = habitat_type)) +
    geom_point(aes(x = year,
               y = mean_richness,
               color = habitat_type),
               cex = 3) +
    geom_errorbar(aes(x = year,
                      ymin = mean_richness - sd_richness,
                      ymax = mean_richness + sd_richness),
                   width=.1,
                  position=position_dodge(0.05)) +
    geom_line(aes(x = year,
                  y = mean_richness,
                  group = habitat_type)) +
    scale_color_nina(name = "Habitatstype",
                     palette = "darkblue-orange") +
    ylab("Antall arter pollinatorer per lokalitet (middelv.)") +
    xlab("År")



## -------------------------------------------------------------------------------------------------------------------------
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
  ungroup()



## -------------------------------------------------------------------------------------------------------------------------
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
         expl_3_main = ifelse(expl_3_main == "Ukjent", NA, expl_3_main)) 
 



## -------------------------------------------------------------------------------------------------------------------------
no_detected_redlist_species <- redlisted_obs_2021 %>% 
  st_drop_geometry() %>% 
  summarise(no_spec = n_distinct(vitenskapelig_navn_id)) %>% 
  pull()


## -------------------------------------------------------------------------------------------------------------------------
norge <- read_sf(con,
                 Id(schema = "backgrounds",
                    table = "norway_terrestrial")
                 ) 


## -------------------------------------------------------------------------------------------------------------------------
redlisted_obs_2021_agg <- redlisted_obs_2021 %>% 
  group_by(locality,
           kategori_2021) %>% 
  summarise(no_spec = n_distinct(species_latin_fixed))  

kat_order <- tibble(kategori_2021 = c("DD",
                                      "NT",
                                      "VU",
                                      "EN",
                                      "CR"),
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


## -------------------------------------------------------------------------------------------------------------------------
no_locales_with_redlisted_species <- redlisted_obs_2021_agg %>% 
  st_drop_geometry() %>% 
  ungroup() %>% 
  summarise(n_distinct(locality)) %>% 
  pull()


## -------------------------------------------------------------------------------------------------------------------------
#redlist_colors <- c("#5a5b5d","#990013", "#d80024" , "#ff4c29", "#f36c21", "#f9a61a")
#names(redlist_colors) <- c("RE", "CR", "EN", "VU", "NT", "DD")

redlist_colors <- c("#f9a61a", "#f36c21", "#ff4c29", "#d80024" )
names(redlist_colors) <- levels(redlisted_obs_2021_agg$kategori_2021)


## ----redlist_map, fig.cap = "Fordeling over funnen av rødlistete arter. Sirkelstørrelsen viser antall funne arter på hver lokalitet, og fargen viser rødlistekategori. Posisjonene er spredt ut noe på kartet for at de ikke skal skjule hverandre. \\label{fig:redlist_map}", fig.height = 8----
set.seed(12345)
norge %>% 
  filter(navn == "Innlandet" |
         navn =="Vestfold og Telemark" |
         navn == "Oslo" |
         navn == "Viken" | 
        navn == "Trøndelag" ) %>% 
ggplot(.) +
  geom_sf() +
  geom_sf(aes(size = no_spec,
              color = kategori_2021),
              alpha = 0.8,
          data = redlisted_obs_2021_agg) +
  scale_color_manual(name = "Rødlistekategori",
                     values = redlist_colors) +
  scale_size(breaks = seq(from = 1,
                          to = max(redlisted_obs_2021_agg$no_spec), 
                          by = 2),
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
  
  


## -------------------------------------------------------------------------------------------------------------------------

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


## ----redlist_map_expl, fig.cap = "Fordeling over funnen av rødlistete arter etter påvirkningsfaktor. Sirkelstørrelsen viser antall funne arter på hver lokalitet, og fargen viser påvirkningsfaktor (flere nivåer) angitt i rødlisten for 2021.  Posisjonene er spredt ut noe på kartet for at de ikke skal skjule hverandre. \\label{fig:redlist_map_expl}", fig.height = 8----
set.seed(12345)
norge %>% 
  filter(navn == "Innlandet" |
         navn =="Vestfold og Telemark" |
         navn == "Oslo" |
         navn == "Viken" | 
        navn == "Trøndelag" ) %>% 
ggplot(.) +
  geom_sf() +
  geom_sf(aes(size = no_spec_loc_expl,
              color = expl_3_main),
              alpha = 0.7,
          data = redlisted_obs_2021_agg_expl_big) +
  scale_color_nina(name = "Påvirkning") +
  scale_size(breaks = seq(from = 1,
                          to = max(redlisted_obs_2021_agg_expl_big$no_spec_loc_expl), 
                          by = 2),
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
  
  


## ----redlisted_obs_agg----------------------------------------------------------------------------------------------------
kat_order <- tibble(kategori_2021 = c("DD",
                                      "NT",
                                      "VU",
                                      "EN",
                                      "CR"),
                    kat_order = 1:5)


redlisted_obs_2021_less_agg <-  redlisted_obs_2021 %>% 
  st_drop_geometry() %>% 
  select(species_latin = species_latin_fixed,
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




## ----redlist_detection, fig.cap = "Deteksjon av rødlistede arter i hver lokalitet. Hver sirkel viser funn av en enkelt art i en lokalitet. Størrelsen på sirklene viser antall ganger arten ble detektert på lokalet (enkelte felleobservasjoner).\\label{fig:redlist_detection}", fig.height = 7----
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



## -------------------------------------------------------------------------------------------------------------------------
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
  select(`Rødlistekategori 2021` = kategori_2021,
         `Orden` = id_order,
         `Latinsk navn` = species_latin,
         `Norsk navn` = popular_name,
         `Funnet i antall lokaler` = no_lok,
         `Antall ganger observert` = no_detected)

spec_no_loc <- spec_no_loc %>% 
  mutate(`Norsk navn` = str_to_title(`Norsk navn`))
  

spec_no_loc1 <- spec_no_loc %>% 
  filter(`Rødlistekategori 2021` == "DD" |
         `Rødlistekategori 2021` == "NT")

spec_no_loc2 <- spec_no_loc %>% 
  filter(`Rødlistekategori 2021` == "VU" |
         `Rødlistekategori 2021` == "EN")

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




## ----redlist_table1, results = "asis"-------------------------------------------------------------------------------------
print(xtable(spec_no_loc1, 
             caption = "Rødlistede arter i kategori DD og NT  i prosjektet så langt. \\label{tab:redlist_table1}",
             align = rep("l", 7)), 
      size = "footnotesize",
      sanitize.colnames.function=function(x){x},
      table.placement = "p")


## ----redlist_table2, results = "asis"-------------------------------------------------------------------------------------
print(xtable(spec_no_loc2, 
             caption = "Rødlistede arter funne i prosjektet så langt. NB at man ikke kan skille Ampedus sanguinolentus fra Ampedus pomonae med de primers vi bruker. \\label{tab:redlist_table2}",
             align = rep("l", 7)), 
      size = "footnotesize",
      table.placement = "p",
      sanitize.colnames.function=function(x){x},)



## ---- fig.cap = "Dørstokkarten Dohrniphora cornuta, som er spritt i store deler av verden gjennom mennesker og som ble funnet i en malaisefelle i 2020, 500 meter fra et hagesenter i Notodden. CC BY-SA 4.0 Aleksandrs Balodis \\label{fig:dohrniphora}", out.width = "0.6\\textwidth"----
include_graphics("figures/dohrniphora_cornuta.png")


## -------------------------------------------------------------------------------------------------------------------------

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
GROUP BY id_phylum, id_class, id_order
ORDER BY id_phylum, id_class, id_order 
"  

no_alien <- dbGetQuery(con,
                       no_alien_q) %>% 
  mutate(count = as.integer(count)) %>%
  select(everything(),
         no_alien = count)



## -------------------------------------------------------------------------------------------------------------------------

no_fennoscand_q <- "
SELECT obs.id_phylum, obs.id_class, obs.id_order, count(distinct obs.species_latin_fixed)
FROM occurrences.observations obs,
lookup.fennoscand_species alien,
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
GROUP BY id_phylum, id_class, id_order
ORDER BY id_phylum, id_class, id_order 
"  

no_fennoscand <- dbGetQuery(con,
                       no_fennoscand_q) %>%
  mutate(count = as.integer(count)) %>%
  select(everything(),
         no_fennoscand_alien = count)



## -------------------------------------------------------------------------------------------------------------------------

no_pot_alien_q <- "
SELECT obs.id_phylum, obs.id_class, obs.id_order, count(distinct obs.species_latin_fixed)
FROM occurrences.observations obs,
lookup.pot_alien_species alien,
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
GROUP BY id_phylum, id_class, id_order
ORDER BY id_phylum, id_class, id_order 
"  

no_pot_alien <- dbGetQuery(con,
                       no_pot_alien_q) %>%
  mutate(count = as.integer(count)) %>%
  select(everything(),
         no_pot_alien = count)



## -------------------------------------------------------------------------------------------------------------------------

no_all_species_q <- "
SELECT obs.id_phylum, obs.id_class, obs.id_order, count(distinct obs.species_latin_fixed)
FROM occurrences.observations obs,
events.identifications,
events.year_locality yl,
events.locality_sampling ls,
locations.localities l,
events.sampling_trap st
WHERE obs.identification_id = identifications.id
AND identifications.sampling_trap_id = st.id
AND st.locality_sampling_id = ls.id
AND ls.year_locality_id = yl.id
AND yl.locality_id = l.id
AND yl.project_short_name = 'NasIns'
GROUP BY id_phylum, id_class, id_order
ORDER BY id_phylum, id_class, id_order 
"  

no_all_species <- dbGetQuery(con,
                       no_all_species_q) %>% 
  mutate(count = as.integer(count)) %>% 
  select(everything(),
         no_all_species = count)



## -------------------------------------------------------------------------------------------------------------------------

no_all_species_to_species_q <- "
SELECT obs.id_phylum, obs.id_class, obs.id_order, count(distinct obs.species_latin_fixed)
FROM occurrences.observations obs,
events.identifications,
events.year_locality yl,
events.locality_sampling ls,
locations.localities l,
events.sampling_trap st
WHERE obs.identification_id = identifications.id
AND identifications.sampling_trap_id = st.id
AND st.locality_sampling_id = ls.id
AND ls.year_locality_id = yl.id
AND yl.locality_id = l.id
AND yl.project_short_name = 'NasIns'
AND obs.taxonomic_level = 'species'
GROUP BY id_phylum, id_class, id_order
ORDER BY id_phylum, id_class, id_order 
"  

no_all_species_to_species <- dbGetQuery(con,
                       no_all_species_to_species_q) %>% 
  mutate(count = as.integer(count)) %>% 
  select(everything(),
         no_all_species_to_species = count)



## -------------------------------------------------------------------------------------------------------------------------

no_all_native_species_to_species_q <- "
SELECT obs.id_phylum, obs.id_class, obs.id_order, count(distinct obs.species_latin_fixed)
FROM occurrences.observations obs,
events.identifications,
events.year_locality yl,
events.locality_sampling ls,
locations.localities l,
events.sampling_trap st
WHERE obs.identification_id = identifications.id
AND identifications.sampling_trap_id = st.id
AND st.locality_sampling_id = ls.id
AND ls.year_locality_id = yl.id
AND yl.locality_id = l.id
AND yl.project_short_name = 'NasIns'
AND obs.taxonomic_level = 'species'
AND obs.native IS TRUE
GROUP BY id_phylum, id_class, id_order
ORDER BY id_phylum, id_class, id_order 
"  

no_all_native_species_to_species <- dbGetQuery(con,
                       no_all_native_species_to_species_q) %>% 
  mutate(count = as.integer(count)) %>% 
  select(everything(),
         no_all_native_species_to_species = count)



## -------------------------------------------------------------------------------------------------------------------------
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




## ---- results = "asis"----------------------------------------------------------------------------------------------------
xtab <- xtable(no_all_aliens, caption = "Antall taksa og arter påvist med DNA-metastrekkoding av alle felleprøver i 2020 og 2021 sortert på orden, antall kjent fra Norge, kjente fremmede arter for Norge, tilstedeværelse i Fennoskandia, og potentielt fremmede arter. Resultatene inneholder også arter observert i uttesting av lengre tømmingsperioder. \\label{tab:no_species_order}")

align(xtab) <- "p{1.5cm}p{1.2cm}p{1.2cm}p{2.2cm}p{1cm}p{1cm}p{1cm}p{1cm}p{1.7cm}p{1cm}"

print.xtable(xtab, sanitize.text.function=function(x){x}, size = "scriptsize")



## -------------------------------------------------------------------------------------------------------------------------
pot_alien_spec <- tbl(con,
                  Id(schema = "lookup",
                     table = "pot_alien_species")) %>% 
  select(species_latin_fixed) %>% 
  pull()

fennoscand_spec <- tbl(con,
                  Id(schema = "lookup",
                     table = "fennoscand_species")) %>% 
  select(species_latin_fixed) %>% 
  pull()

alien_spec <- tbl(con,
                  Id(schema = "lookup",
                     table = "fremmedartslista_2018_artsdatabanken")) %>% 
  filter(expertGroup == "Ikke-marine invertebrater") %>% 
  select(scientificName) %>%
  pull()



## -------------------------------------------------------------------------------------------------------------------------

fennoskand_obs_q <- "
SELECT 
obs.species_latin_fixed,
yl.year,
ls.sampling_name,
l.locality,
ST_Centroid(l.geom) as geom
FROM occurrences.observations obs,
lookup.fennoscand_species fennoscand,
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
"

fennoskand_obs <- read_sf(con,
                          query = fennoskand_obs_q) %>% 
  mutate(alien_kat = "Fennoskandisk forek.")



## -------------------------------------------------------------------------------------------------------------------------

pot_alien_obs_q <- "
SELECT 
obs.species_latin_fixed,
yl.year,
ls.sampling_name,
l.locality,
ST_Centroid(l.geom) as geom
FROM occurrences.observations obs,
lookup.pot_alien_species pot_alien,
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
"

pot_alien_obs <- read_sf(con,
                          query = pot_alien_obs_q) %>% 
  mutate(alien_kat = "Ikke fennoskandisk forek.")



## -------------------------------------------------------------------------------------------------------------------------

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
AND alien.\"riskCategory\" IN ('SE', 'HI', 'PH', 'LO')
"  

alien_obs <- read_sf(con,
                     query = alien_obs_q) %>% 
  mutate(alien_kat = "På fremmedartslista")



## -------------------------------------------------------------------------------------------------------------------------
all_alien_obs <- fennoskand_obs %>% 
  rbind(pot_alien_obs)  %>% 
  rbind(alien_obs)


## -------------------------------------------------------------------------------------------------------------------------
kat_order <- tibble(alien_kat = c("Fennoskandisk forek.",
                                  "Ikke fennoskandisk forek.",
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


# kat_order <- tibble(expl_3_main = c("Ukjent/Ikke angitt",
#                                          "Buskap/dyrehold",
#                                          "Forurensing",
#                                          "Fremmede arter",
#                                          "Jordbruk",
#                                          "Opphørt/redusert drift",
#                                          "Påvirker habitatet (beite tråkk mm.)",
#                                          "Skogbruk (kommersielt)",
#                                          "Utbygging/utvinning"),
#                     kat_order = 1:9)
# 
# redlisted_obs_2021_agg_expl <- redlisted_obs_2021_agg_expl %>% 
#   left_join(kat_order,
#             by = c("expl_3_main" = "expl_3_main")) %>% 
#   group_by(expl_3_main) %>% 
#   mutate(no_spec_per_expl = n()) %>% 
#   ungroup() %>% 
#   mutate(expl_3_main = factor(paste0(expl_3_main, " (", no_spec_per_expl, " st.)"))) %>% 
#   mutate(expl_3_main = fct_reorder(expl_3_main, kat_order))
# 
# #levels(redlisted_obs_2021_agg_expl$expl_3_main)

all_alien_obs_agg <- all_alien_obs_agg %>% 
  st_jitter(all_alien_obs_agg, amount = 7000) 


## ----non_native_species_kat, fig.cap = "Observasjoner av insektarter med ikke tidligere kjent forekomst i Norge. Observasjonene er basert på metastrekkoding. Arter vurdert til NK (Ingen kjent risiko) eller NR (Ikke risikovurdert) på fremmedartslista er ikke tatt med. Posisjonene er spredt ut noe på kartet for at de ikke skal skjule hverandre.  \\label{fig:non_native_species_kat}", fig.height = 8----
set.seed(12345)
norge %>% 
  filter(navn == "Innlandet" |
         navn =="Vestfold og Telemark" |
         navn == "Oslo" |
         navn == "Viken" | 
        navn == "Trøndelag" ) %>% 
ggplot(.) +
  geom_sf() +
  geom_sf(aes(size = no_spec,
              color = alien_kat_append),
              alpha = 0.7,
          data = all_alien_obs_agg) +
  scale_color_nina(name = "Kategori") +
  scale_size(breaks = c(seq(from = min(all_alien_obs_agg$no_spec),
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



## ----alien_obs_agg--------------------------------------------------------------------------------------------------------
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



## ----alien_detection, fig.cap = "Deteksjon av arter som ikke tidligere har kjent forekomst i Norge. Hver sirkel viser funn av en enkelt art i en lokalitet. Arter vurdert til NK (Ingen kjent risiko) eller NR (Ikke risikovurdert) på fremmedartslista er ikke tatt med. Størrelsen på sirklene viser antall ganger arten ble detektert på lokalet (enkelte felleobservasjoner).\\label{fig:alien_detection}", fig.height = 7----
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



## ----richn_shannon_year_loc, fig.cap = "Sammenheng mellom Shannon-indeks basert på den genetiske diversiteten hos hver art og artsantall, summert over lokaliteter og år. \\label{fig:richn_shannon_year_loc}", out.width = "0.7\\textwidth"----
year_locality_biodiv %>%  
  mutate(habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type)) %>% 
  ggplot(aes(x = no_species, 
                 y = shannon_div)) +
  geom_point(aes(color = habitat_type),
             cex = 2) +
  geom_smooth(method = "lm",
              se = F,
              color = ninaPalette()[2]) +
  ylab("Shannon indeks (ASV-antall)") +
  xlab("Antall detekterte insektarter per lokal og år") +
  scale_color_nina(name = "Habitat type",
                   palette = "darkblue-orange") 


## ----richn_shannon_loc_sampl, fig.cap = "Sammenheng mellom Shannon-indeks basert på den genetiske diversiteten hos hver art og artsantall, summert over lokalitet og tømmeperiode. \\label{fig:richn_shannon_loc_sampling}", out.width = "0.7\\textwidth"----
locality_sampling_biodiv  %>% 
  mutate(habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type)) %>% 
  ggplot(aes(x = no_species, 
                 y = shannon_div)) +
  geom_point(aes(color = habitat_type),
             cex = 2) +
  geom_smooth(method = "lm",
              se = F,
              color = ninaPalette()[2]) +
  ylab("Shannon-indeks (ASV)") +
  xlab("Antall detekterte insektarter per lokal og sampletilfelle") +
  scale_color_nina(name = "Habitat type",
                   palette = "darkblue-orange") 


## -------------------------------------------------------------------------------------------------------------------------
# locality_sampling_weight_biodiv <- locality_sampling_biodiv %>% 
#   left_join(biomass_mf_locality_sampling,
#             by = c("year" = "year",
#                    "locality" = "locality",
#                    "sampling_name" = "sampling_name",
#                    "habitat_type" = "habitat_type",
#                    "region_name" = "region_name")) %>% 
#   mutate(habitat_type == ifelse(habitat_type == "Forest", "Skog", habitat_type))


#Joining ser ut å fungere, også uten å spesifisere. m.a.o. no_trap_days är lika mellan dataseten.
# locality_sampling_weight_biodiv %>% 
#   mutate(diff_trap_days = no_trap_days.x - no_trap_days.y) %>% 
#   select(diff_trap_days) %>% 
#   pull()
  

locality_sampling_weight_biodiv <- locality_sampling_biodiv %>% 
  left_join(biomass_mf_locality_sampling) %>% 
  mutate(habitat_type = ifelse(habitat_type == "Forest", "Skog", habitat_type))



## ----locality_sampling_weight_biodiv, fig.cap = "Sammenheng mellom biomasse og antall detekterte taksa i felleprøver. A estimert linjært sammenheng, B estimert polynomisk sammenheng. \\label{fig:locality_sampling_weight_biodiv}"----
p_lin <- locality_sampling_weight_biodiv %>% 
  ggplot(aes(x = avg_wet_weight,
             y = no_species)) +
  geom_point(aes(color = habitat_type)) +
  scale_color_nina(name = "Habitat type",
                   palette = "darkblue-orange") +
  ylab("Antall detekterte insektarter per lokal og sampletilfelle") +
  xlab("Avrunnen våtvekt (g.)") +
      geom_smooth(aes(color = habitat_type),
                  method = "glm",
                  formula = y ~ x,
              # method = "nls",
              # formula = y ~ a * x + b,
              # se = F,
              # method.args = list(start = list(a = 0.1,
              #                                 b = 0.1)
              #                    ),
              lwd = 2
              ) +
   theme(legend.position = "bottom",
       # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
       legend.key.size = unit(0.5, 'cm'),
       legend.text = element_text(size=7),
       legend.key.width = unit(0.5, 'cm')) +
  scale_y_continuous(breaks = c(0, 250, 500, 750, 1000), 
                   name = "Antall detekterte insektarter per lokal og sampletilfelle",
                   limits = c(0, 1000)) +
  ggtitle("A")
  



p_poly <- locality_sampling_weight_biodiv %>% 
  ggplot(aes(x = avg_wet_weight,
             y = no_species)) +
  geom_point(aes(color = habitat_type)) +
  scale_color_nina(name = "Habitat type",
                   palette = "darkblue-orange") +
  #ylab("Antall detekterte insektarter per lokal og sampletilfelle") +
  xlab("Avrunnen våtvekt (g.)") +
      geom_smooth(aes(color = habitat_type),
                  method = "glm",
                  formula = y ~ poly(x, 2),
              # method = "nls",
              # formula = y ~ a * x + b,
              # se = F,
              # method.args = list(start = list(a = 0.1,
              #                                 b = 0.1)
              #                    ),
              lwd = 2
              )  +
   theme(legend.position = "bottom",
       # legend.box.margin = unit(x = c(0,40,0,40), units = "mm"),
       legend.key.size = unit(0.5, 'cm'),
       legend.text = element_text(size=7),
       legend.key.width = unit(0.5, 'cm')) +
  scale_y_continuous(breaks = c(0, 250, 500, 750, 1000), 
                   name = "",
                   limits = c(0, 1000)) +
  ggtitle("B")
  



lm <-  matrix(c(1:4, 5, 5),
              ncol = 2,
              byrow = T)


marrangeGrob(list(p_lin,
                  p_poly
                  ),
             top = "",
             nrow = 1,
             ncol = 2)



## -------------------------------------------------------------------------------------------------------------------------
locality_sampling_biomass_temp <- biomass_mf_locality_sampling %>% 
  left_join(locality_sampling_loggerdata,
            by = c("sampling_name" = "sampling_name",
                   "year" = "year",
                   "locality" = "locality",
                   "habitat_type" = "habitat_type",
                   "region_name" = "region_name"))

locality_sampling_biomass_temp_compl_light <- locality_sampling_biomass_temp %>% 
  filter(!is.na(sum_wet_weight),
         !is.na(avg_values_MX2301A_temperature),
         !is.na(avg_values_MX2301A_rh),
         !is.na(avg_values_MX2202_light),
         !is.na(stddev_values_MX2301A_temperature),
         !is.na(stddev_values_MX2301A_rh),
         !is.na(stddev_values_MX2202_light))




## -------------------------------------------------------------------------------------------------------------------------
global_model <- lmer(sum_wet_weight/no_trap_days ~ 
                       scale(avg_values_MX2301A_temperature) +
                       scale(avg_values_MX2301A_rh) + 
                       scale(avg_values_MX2202_light) +
                       scale(stddev_values_MX2301A_temperature) +
                       scale(stddev_values_MX2301A_rh) + 
                       scale(stddev_values_MX2202_light) +
                       habitat_type + 
                       (1|locality), 
                     data = locality_sampling_biomass_temp_compl_light,
                     na.action = "na.fail",
                     REML = F)

cand_models <- dredge(global_model,
                      fixed = c("habitat_type"))


## -------------------------------------------------------------------------------------------------------------------------
model_averages <- model.avg(cand_models, subset = delta <= 4)
#summary(model_averages)


## -------------------------------------------------------------------------------------------------------------------------
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



## ----temp_AICtab, results = "asis"----------------------------------------------------------------------------------------
print(xtable(biomass_temp_AIC, caption = "AICc-vekting av de beste modellene av biomasse fra loggerdata. Kun de modeller med max 4 delta-AIC fra den beste modellen er tatt med. \\label{tab:temp_AICtab}" ))



## -------------------------------------------------------------------------------------------------------------------------
par_weights <- MuMIn::sw(get.models(cand_models, subset = delta <= 4)) 

par_weights <- tibble("Parameter" = attr(par_weights, "names"),
                      "Total vekt" = as.numeric(par_weights),
                      "Med i antall modeller" = attr(par_weights, "n.models")
  )

par_vals_full <- round(model_averages$coefficients[1,], 3) %>% 
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
  



## ----logger_model_par, results = "asis"-----------------------------------------------------------------------------------
xtable(par_weights, digits = 3, caption = "Vekter av parametere i de beste modellene (delta-AICc < 4) av innsamlet biomasse og deres vektede parameterestimat.\\label{tab:logger_model_par}")


## -------------------------------------------------------------------------------------------------------------------------
#sw(model_averages)

r2_verdi_biomasse_temp_glob <- r.squaredGLMM(global_model)[1]


## -------------------------------------------------------------------------------------------------------------------------
nas_ins_lok_1000m_skog_q <- "
SELECT skog.*
FROM backgrounds.locality_1000m_buffer_skogtype_agg skog,
events.year_locality yl
WHERE skog.locality_id = yl.locality_id
AND yl.project_short_name = 'NasIns'

"


nas_ins_lok_1000m_skog <- read_sf(con,
                                query = nas_ins_lok_1000m_skog_q)



## -------------------------------------------------------------------------------------------------------------------------

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
  mutate_if(is.numeric , replace_na, replace = 0) %>% 
  mutate_if(is.numeric , units::drop_units)


## -------------------------------------------------------------------------------------------------------------------------
locality_skog_area_data <- locality_skog_area %>% 
  ungroup() %>% 
  select(-"locality")

skog_rda <- decorana(locality_skog_area_data)

#summary(skog_rda) #hm, too large axis length?


## -------------------------------------------------------------------------------------------------------------------------
skog_PCA <- rda(locality_skog_area_data, scale = F)

#summary(skog_PCA)


## -------------------------------------------------------------------------------------------------------------------------

skog_scores <- scores(skog_PCA)$sites %>% 
  as_tibble() %>% 
  cbind(locality_skog_area["locality"]) %>% 
  select(locality,
         skog_PC1 = PC1,
         skog_PC2 = PC2)


locality_skog_area <- locality_skog_area %>% 
  left_join(skog_scores,
            by = c("locality" = "locality"))


## -------------------------------------------------------------------------------------------------------------------------
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


## -------------------------------------------------------------------------------------------------------------------------
hm <- autoplot(ar5_PCA) +
  scale_color_nina()


## ----ar5_PCA, fig.pos = "!bh",  fig.cap = "PCA-analyse over landskapskomposisjon i arealtype etter AR5. Verdiene er ikke normalisert, da alle er fra samme datakilde og samme skala (km^2). \\label{fig:ar5_PCA}", fig.dim = c(5, 4)----
hm



## -------------------------------------------------------------------------------------------------------------------------
tt <- autoplot(skog_PCA) +
  scale_color_nina()


## ----skog_PCA, fig.pos = "!hb", fig.cap = "PCA-analyse over landskapskomposisjon i skogtype og bonitet etter AR5 (basert på kartdata fra SR16). Verdiene er ikke normalisert, da alle er fra samme datakilde og samme skala (km^2). \\label{fig:skog_PCA}", fig.dim = c(5, 4)----
tt


## ----ar5_shannon_plot, fig.pos = "!h", fig.cap = "Fordeling av Shannon index av AR5:arealtype innen en buffer med radius 1000m kring lokalitetene. \\label{fig:ar5_shannon_plot}", out.width = "\\textwidth"----

locality_ar5_shdi %>% 
  mutate(habitat_type = ifelse(grepl("Semi-nat", locality), "Semi-nat", "Skog")) %>% 
  ggplot(.) +
  geom_point(aes(x = locality,
                 y = lsm_l_shdi,
                 color = habitat_type),
               lwd = 2)  + 
  scale_color_nina(name = "Habitattype")+
  xlab("Lokalitet") +
  ylab("Shannon indeks") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4)) 
  


## ----ar5_perimeter_plot, fig.pos = "!h", fig.cap = "Fordeling av total lengde av grenser mellom arealtyper, innen et buffertområde på  1 km kring lokalitetene. \\label{fig:ar5_perimeter_plot}", out.width = "\\textwidth"----

locality_ar5_tot_perimeter %>% 
  mutate(habitat_type = ifelse(grepl("Semi-nat", locality), "Semi-nat", "Skog")) %>% 
  ggplot(.) +
  geom_point(aes(x = locality,
                 y = tot_perimeter_length,
                 color = habitat_type),
               lwd = 2)  + 
  scale_color_nina(name = "Habitattype")+
  xlab("Lokalitet") +
  ylab("Tot. lengde linjære element") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4)) 
  


## -------------------------------------------------------------------------------------------------------------------------
locality_ar5_area <- locality_ar5_area %>% 
  left_join(ar5_scores,
            by = c("locality" = "locality"))



## -------------------------------------------------------------------------------------------------------------------------
locality_sampling_biomass_big_join <- locality_sampling_biomass_temp_compl_light %>% 
  left_join(locality_ar5_area,
            by = c("locality" = "locality")) %>% 
  left_join(locality_ar5_shdi,
            by = c("locality" = "locality")) %>% 
  left_join(locality_ar5_tot_perimeter,
            by = c("locality" = "locality")) %>% 
  left_join(locality_skog_area,
            by = c("locality" = "locality")) %>% 
  left_join(ano_herb_agg,
            by = c("locality" = "locality"),
            copy = T) %>% 
  left_join(tree_agg,
            by = c("locality" = "locality"),
            copy = T)

write_csv(locality_sampling_biomass_big_join,
          file = "out/locality_sampling_biomass_big_join.csv")



## -------------------------------------------------------------------------------------------------------------------------
locality_sampling_biomass_big_join_forest <- locality_sampling_biomass_big_join %>% 
  filter(habitat_type == "Forest")

locality_sampling_biomass_big_join_semi_nat <- locality_sampling_biomass_big_join %>% 
  filter(habitat_type == "Semi-nat")



## ----perim_shann_plot, fig.pos = "!h", fig.cap = "Sammenheng mellom diversiteten av arealtyper (Shannon index) og mengden linjære element i landskapene rundt skogslokalitetene. \\label{fig:perim_shann_plot}", out.width = "\\textwidth"----
locality_sampling_biomass_big_join_forest %>% 
  ggplot(.) +
  geom_point(aes(x = lsm_l_shdi,
                 y = tot_perimeter_length),
                 color = ninaColors()[1]) +
  ylab("Total lengde linjære element") +
  xlab("Artype Shannon indeks") 


## -------------------------------------------------------------------------------------------------------------------------
# global_model_big_forest <- lmer(sum_wet_weight/no_trap_days ~ scale(avg_values_MX2301A_temperature) + 
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
#                         dom_tree  +
#                         (1|locality), 
#                       data = locality_sampling_biomass_big_join_forest,
#                       na.action = "na.fail",
#                       REML = F,
#                       control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))) #Suppress warnings due to near zero estimate for locality random effect

global_model_big_forest <- glm(sum_wet_weight/no_trap_days ~ scale(avg_values_MX2301A_temperature) + 
                       # scale(avg_values_MX2301A_rh) + 
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
                        scale(avg_age) +
                        dom_tree, 
                        data = locality_sampling_biomass_big_join_forest,
                        na.action = "na.fail")




cand_models_big_forest <- dredge(global_model_big_forest,
                      fixed = c("scale(avg_values_MX2301A_temperature)",
                                #"scale(avg_values_MX2301A_rh)",
                                "scale(avg_values_MX2202_light)",
                                "scale(stddev_values_MX2301A_temperature)",
                                "scale(stddev_values_MX2301A_rh)",
                                "scale(stddev_values_MX2202_light)"
                                )
                      )



## -------------------------------------------------------------------------------------------------------------------------
model_averages_big_forest <- model.avg(cand_models_big_forest, subset = delta <= 4)
#summary(model_averages_big_forest)


## -------------------------------------------------------------------------------------------------------------------------
top_models_big_forest <- get.models(cand_models_big_forest, subset = delta <= 4)

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



## ----temp_AICtab_big_forest, results = "asis"-----------------------------------------------------------------------------

print(xtable(biomass_temp_AIC_big_forest, caption = "AICc-vekting av de beste modellene av biomasse fra skogslokaliteter. Kun de modeller med max 5 delta-AIC fra den beste modellen er tatt med. \\label{tab:temp_AICtab_big_forest}" ))



## -------------------------------------------------------------------------------------------------------------------------
par_weights_forest <- MuMIn::sw(get.models(cand_models_big_forest, subset = delta <= 4)) 

par_weights_forest <- tibble("Parameter" = attr(par_weights_forest, "names"),
                      "Total vekt" = as.numeric(par_weights_forest),
                      "Med i antall modeller" = attr(par_weights_forest, "n.models")
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
         Parameter = ifelse(Parameter == "scale(ano_mean_no_spec)", "ANO: artantall karpl.", Parameter),
         Parameter = ifelse(Parameter == "scale(ar5_PC1)", "AR5-arealtype: PCA-aksel 1", Parameter),
         Parameter = ifelse(Parameter == "scale(ar5_PC2)", "AR5-arealtype: PCA-aksel 2", Parameter),
         Parameter = ifelse(Parameter == "scale(skog_PC1)", "AR5-skog: PCA-aksel 1", Parameter),
         Parameter = ifelse(Parameter == "scale(skog_PC2)", "AR5-skog: PCA-aksel 2", Parameter),
         Parameter = ifelse(Parameter == "scale(lsm_l_shdi)", "AR5-arealtype: Shannon indeks", Parameter),
         Parameter = ifelse(Parameter == "scale(tot_perimeter_length)", "AR5-arealtype: linjære elementer", Parameter),
         ) %>% 
  arrange(desc(`Total vekt`), Parameter)
  





## ----par_weights_forest, results = "asis"---------------------------------------------------------------------------------
xtable(par_weights_forest, caption = "Vekter av miljøparametere i de beste modellene av innsamlet biomasse i skog (delta-AICc <= 4). \\label{tab:par_weights_forest}.")


## -------------------------------------------------------------------------------------------------------------------------
r2_verdi_biomasse_big_forest <- r.squaredGLMM(get.models(cand_models_big_forest, subset = delta < 0.01)[[1]])


## -------------------------------------------------------------------------------------------------------------------------
global_model_big_semi_nat <- lmer(sum_wet_weight/no_trap_days ~ scale(avg_values_MX2301A_temperature) +
                       # scale(avg_values_MX2301A_rh) +
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
                        (1|locality),
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
                                #"scale(avg_values_MX2301A_rh)",
                                "scale(avg_values_MX2202_light)",
                                "scale(stddev_values_MX2301A_temperature)",
                                "scale(stddev_values_MX2301A_rh)",
                                "scale(stddev_values_MX2202_light)"
                                )
                      )



## -------------------------------------------------------------------------------------------------------------------------
model_averages_big_semi_nat <- model.avg(cand_models_big_semi_nat, subset = delta <= 4)
#summary(model_averages_big_semi_nat)


## -------------------------------------------------------------------------------------------------------------------------
top_models_big_semi_nat <- get.models(cand_models_big_semi_nat, subset = delta <= 4)

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



## ----temp_AICtab_big_semi_nat, results = "asis"---------------------------------------------------------------------------

print(xtable(biomass_temp_AIC_big_semi_nat, caption = "AICc-vekting av de beste modellene av biomasse fra semi-naturelle lokaliteter. Kun de modeller med max 4 delta-AIC fra den beste modellen er tatt med. \\label{tab:temp_AICtab_big_semi_nat}" ))



## -------------------------------------------------------------------------------------------------------------------------
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
         ) %>% 
  mutate(Parameter = ifelse(Parameter == "dom_tree", "Landsskog: dominerende treslag", Parameter),
         Parameter = ifelse(Parameter == "scale(avg_age)", "Landsskog: skogalder", Parameter),
         Parameter = ifelse(Parameter == "scale(ano_mean_cover)", "ANO: dekningsgrad kartpl.", Parameter),
         Parameter = ifelse(Parameter == "scale(ano_mean_no_spec)", "ANO: artantall kartpl.", Parameter),
         Parameter = ifelse(Parameter == "scale(ar5_PC1)", "AR5-arealtype: PCA-aksel 1", Parameter),
         Parameter = ifelse(Parameter == "scale(ar5_PC2)", "AR5-arealtype: PCA-aksel 2", Parameter),
         Parameter = ifelse(Parameter == "scale(skog_PC1)", "AR5-skog: PCA-aksel 1", Parameter),
         Parameter = ifelse(Parameter == "scale(skog_PC2)", "AR5-skog: PCA-aksel 2", Parameter),
         Parameter = ifelse(Parameter == "scale(lsm_l_shdi)", "AR5-arealtype: Shannon indeks", Parameter),
         Parameter = ifelse(Parameter == "scale(tot_perimeter_length)", "AR5-arealtype: linjære elementer", Parameter),
         ) %>% 
  arrange(desc(`Total vekt`), Parameter)
  




## ----results = "asis"-----------------------------------------------------------------------------------------------------
xtable(par_weights_semi_nat, caption = "Vekter av miljøparametere i de beste modellene av innsamlet biomasse i semi-naturelle lokaliteter (delta-AICc <= 4).\\label{tab:par_weights_semi_nat}")


## -------------------------------------------------------------------------------------------------------------------------
r2_verdi_biomasse_big_semi_nat <- r.squaredGLMM(get.models(cand_models_big_semi_nat, subset = delta < 0.1)[[1]])


## ---- results = "hide", warning = F, message = F--------------------------------------------------------------------------
ord.malaise=scores(metaMDS(subset_samples(physeq.species,grepl("Trond|Ost",physeq.species@sam_data$fylke)&felletype=="malaise")@otu_table@.Data,distance="bray",trymax = 1000))
ord.malaise=data.frame(ord.malaise,data.frame(subset_samples(physeq.species,grepl("Trond|Ost",physeq.species@sam_data$fylke)&felletype=="malaise")@sam_data))


ord.seminat=(metaMDS(subset_samples(physeq.species,habitat=="Semi-nat"&felletype=="malaise")@otu_table@.Data,distance="bray",trymax = 1000))

points=data.frame(ord.seminat$points)
points=right_join(data.frame(physeq.species@sam_data),data.frame(sample=rownames(points),points))
fit=envfit(ord.seminat~avg_values_MX2301A_temperature+avg_values_MX2202_light+ar5_PC1+ar5_PC2,data=data.frame(subset_samples(physeq.species,habitat=="Semi-nat"&felletype=="malaise")@sam_data),na.rm=T)

fit.scores=data.frame(scores(fit, display = "vectors"))
fit.scores=cbind(fit.scores, vars = rownames(fit.scores))


## -------------------------------------------------------------------------------------------------------------------------
ord.malaise <- ord.malaise %>% 
  select(Habitattype = habitat,
         "Tømmetilfelle" = proevetaking,
         NMDS1,
         NMDS2)

fig.habitats <- 
  ggplot(ord.malaise,
                       aes(x=NMDS1,
                           y=NMDS2,
                           shape=Habitattype,
                           color=Tømmetilfelle)) +
  geom_point() +
  stat_ellipse(aes(x=NMDS1,
                   y=NMDS2,
                   group=Habitattype),
                   show.legend = F,
               data = ord.malaise,
               inherit.aes = F) +
  scale_color_manual(values=colorRampPalette(c("white","#93328E"))(12)[4:12]) + 
  theme(panel.background=element_blank(),axis.line = element_line(colour="black"))
  


## ---- results = "hide", warning = F, message = F--------------------------------------------------------------------------
fig.seminat <- ggplot(points) +
  geom_point(data=points,
             aes(x=MDS1,y=MDS2,
                 color=lokalitet)) +
  scale_color_manual(name = "Lokalitet",
                     values=paste0(Polychrome::palette36.colors())) +
  coord_fixed() +
  theme(legend.position = "bottom",
       legend.key.size = unit(0.5, 'cm'),
       legend.text = element_text(size=7),
       legend.key.width = unit(0.5, 'cm'))

fig.seminat <- fig.seminat + 
  geom_segment(data = fit.scores[1:3,],
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               show.legend = F,
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_segment(data = fit.scores[4,],
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2 ),
               linetype=2,
               show.legend = F,
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_text(data = fit.scores,
            aes(x = NMDS1-0.05, y = NMDS2-0.05, label = c("temp","lys","ar5_PC1","ar5_PC2")),
            size = 3) + 
   ggtitle("Semi-Nat")


## ---- results = "hide", warning = F, message = F--------------------------------------------------------------------------
ord.skog=(vegan::metaMDS(phyloseq::subset_samples(physeq.species,habitat=="Skog"&felletype=="malaise")@otu_table@.Data,distance="bray",trymax = 1000))

points=data.frame(ord.skog$points)
points=right_join(data.frame(physeq.species@sam_data),data.frame(sample=rownames(points),points))
fit=envfit(ord.skog~avg_values_MX2301A_temperature+avg_values_MX2202_light+skog_PC1+avg_age,data=data.frame(subset_samples(physeq.species,habitat=="Skog"&felletype=="malaise")@sam_data),na.rm=T)

fit.scores=data.frame(scores(fit, display = "vectors"))
fit.scores=cbind(fit.scores, vars = rownames(fit.scores))


## ---- results = "hide", warning = F, message = F--------------------------------------------------------------------------
fig.skog=ggplot(points) +
  geom_point(data=points,
             aes(x=MDS1,y=MDS2,color=lokalitet)) +
  scale_color_manual(name = "Lokalitet",
                     values=paste0(Polychrome::palette36.colors())) +
  coord_fixed()  +
  theme(legend.position = "bottom")

fig.skog=fig.skog + 
  geom_segment(data = fit.scores[-2,],
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               show.legend = F,
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_segment(data = fit.scores[2,],
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               show.legend =F,
               linetype=2,
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_text(data = fit.scores, aes(x = NMDS1-0.05, y = NMDS2-0.05, label = c("temp","lys","skog_PC1","age")),
            size = 3)+ggtitle("Skog")




## ---- fig.pos = "!h", fig.cap = "Ordinasjon som viser hovedmønstrene i betadiversitet i alle insektprøvene innsammlet i 2021. Punktene er farget etter innsamlingsdato og ellipsene viser den 95% konfidens interval for hver habitat.  \\label{fig:beta_div_patterns_overall}"----
fig.habitats


## ----betadiv_pattern_semi_nat, fig.pos = "!hb", fig.cap = "Ordinasjon som viser sammenheng mellom forklaringsvariabler og mønstrer i samfunnskomposisjon i semi-naturlig mark. Solide linier viser signifikant sammenhenger.  \\label{fig:betadiv_pattern_semi_nat}", out.width = "0.8\\textwidth"----
fig.seminat


## ----betadiv_pattern_skog, fig.pos = "!hb", fig.cap = "Ordinasjon som viser sammenheng mellom forklaringsvariabler og mønstrer i samfunnskomposisjon i skog. Solide linier viser signifikant sammenhenger. \\label{fig:betadiv_pattern_skog}", out.width = "0.8\\textwidth"----
fig.skog


## -------------------------------------------------------------------------------------------------------------------------
coleopt_species <- obs_from_db(subset_order = "Coleoptera",
                               agg_level = "year_locality")


## -------------------------------------------------------------------------------------------------------------------------
saprox_species <- tbl(con,
                      Id(schema = "lookup",
                         table = "nmbu_saproxylic_beetles"))  %>% 
  filter(saproxylic) %>% 
  collect() %>% 
  pull(species_latin)



## -------------------------------------------------------------------------------------------------------------------------
saprox_richness <- obs_from_db(subset_species = saprox_species,
                               agg_level = "year_locality")
#saprox_richness


## -------------------------------------------------------------------------------------------------------------------------
saprox_tree_agg <- tree_agg  %>% 
  left_join(saprox_richness,
            by = c("locality" = "locality",
                   "year" = "year"),
            copy = T)

#saprox_tree_agg


## -------------------------------------------------------------------------------------------------------------------------
redlisted_species_q <- "
SELECT DISTINCT ON(obs.species_latin_fixed) obs.species_latin_fixed
FROM occurrences.observations obs LEFT JOIN
lookup.redlist_2015 redlist
ON obs.species_latin_fixed = redlist.species_latin
WHERE redlist.kategori IN ('DD', 'NT', 'VU', 'EN', 'CR', 'RE')
GROUP BY obs.species_latin_fixed
"

redlisted_species <- dbGetQuery(con,
                                redlisted_species_q)



## -------------------------------------------------------------------------------------------------------------------------

saprox_redlisted_species <- redlisted_species %>% 
  filter(species_latin_fixed %in% saprox_species) %>% 
  pull()
 
#saprox_redlisted_species



## -------------------------------------------------------------------------------------------------------------------------
redlisted_saprox_richness <- obs_from_db(subset_species = saprox_redlisted_species,
                               agg_level = "year_locality")

#redlisted_saprox_richness


## -------------------------------------------------------------------------------------------------------------------------
redlist_saprox_tree_agg <- tree_agg  %>% 
  left_join(redlisted_saprox_richness,
            by = c("locality" = "locality",
                   "year" = "year"),
            copy = T)

#redlist_saprox_tree_agg


## -------------------------------------------------------------------------------------------------------------------------

redlisted_saprox_richness <- redlisted_saprox_richness %>% 
  select(no_species_redlist = no_species,
         everything())

coleopt_tree_agg <- tree_agg  %>% 
  left_join(coleopt_species,
            by = c("locality" = "locality",
                   "year" = "year"),
            copy = T) %>% 
  left_join(saprox_richness,
            by = c("year" = "year",
                   "locality" = "locality",
                   "habitat_type" = "habitat_type",
                   "region_name" = "region_name"),
            suffix = c("_coleopt", "_saprox"),
            copy = T
            ) %>% 
    left_join(redlisted_saprox_richness,
            by = c("year" = "year",
                   "locality" = "locality",
                   "habitat_type" = "habitat_type",
                   "region_name" = "region_name"),
            suffix = c("", "_redlist"),
            copy = T)

coleopt_richn_tree_agg <- coleopt_tree_agg %>% 
  select(year,
         locality,
         skogkarakter,
         dom_tree,
         avg_age,
         region_name,
         no_species_coleopt,
         no_species_saprox,
         no_species_redlist) %>% 
  pivot_longer(cols = c("no_species_coleopt", "no_species_saprox", "no_species_redlist"),
               names_to = "category",
               names_prefix = "no_species_",
               values_to = "no_spec") %>%
    collect() %>% 
  mutate(category = ifelse(category == "coleopt", "Alle", category),
         category = ifelse(category == "redlist", "Vedlevende på rødliste", category),
         category = ifelse(category == "saprox", "Vedlevende", category)) %>% 
  mutate(no_spec = ifelse(is.na(no_spec), 0, no_spec),
         category = factor(category, levels = c("Alle", "Vedlevende", "Vedlevende på rødliste")))



## ----coleopt_art_bar,fig.pos = "!h", fig.cap = "Antall identifiserte billearter per lokalitet. Fangstene kommer fra 1 malaisefelle og 4 vindusfeller per lokalitet.\\label{fig:coleopt_art_bar}", out.width = "\\textwidth"----

coleopt_richn_tree_agg %>% 
  #filter(category == "Alle") %>% 
ggplot(., aes(x = locality,
               y = no_spec,
               fill = category),
       width = 5) +
  geom_bar(position = "stack",
           stat = "identity") +
  scale_fill_nina(name = "Artskategori") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Antall identifiserte billearter") +
  xlab("Lokalitet")



## ----saprox_richness_forest_age, fig.cap = "Antall vedlevende arter biller funnet i skogslokalitetene i 2020-2021. Fangstene kommer fra 1 malaisefelle og 4 vindusfeller per lokalitet. \\label{fig:saprox_richness_forest_age}", out.width = "\\textwidth", fig.pos = "!h"----
coleopt_richn_tree_agg %>% 
  #filter(category == "Alle") %>% 
ggplot(., aes(x = avg_age,
               y = no_spec,
               color = category)) +
  geom_point(size = 2) +
  scale_fill_nina(name = "Artskategori") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Antall identifiserte billearter") +
  xlab("Alder på største trær") +
  scale_color_nina(name = "Artskategori") +
  geom_smooth(aes(x = avg_age,
                  y = no_spec,
                  color = category),
              method = "lm"
              )




## -------------------------------------------------------------------------------------------------------------------------
bees <- c("Andrenidae", 
          "Apidae",
          "Colletidae",
          "Halictidae",
          "Megachilidae",
          "Melittidae")

syrph <- c("Syrphidae")  

lepi <- c("Hesperiidae",
          "Lycaenidae",
          "Nymphalidae",
          "Papilionidae",
          "Pieridae",
          "Rionidae")

year_locality_bees <- obs_from_db(subset_families = bees,
                                             agg_level = "year_locality",
                                             trap_type = "MF") %>% 
  mutate(group = "Bier")
year_locality_lepi <- obs_from_db(subset_families = syrph,
                                             agg_level = "year_locality",
                                             trap_type = "MF") %>% 
  mutate(group = "Sommerfugler")

year_locality_syrph <- obs_from_db(subset_families = lepi,
                                             agg_level = "year_locality",
                                             trap_type = "MF") %>% 
  mutate(group = "Blomsterfluer")

year_locality_poll <- year_locality_bees %>% 
  bind_rows(year_locality_lepi) %>% 
  bind_rows(year_locality_syrph)



## ----poll_art_bar, fig.cap = "Antall identifiserte pollinatorer per lokal. Fangstene kommer fra 1 malaisefelle per lokalitet.\\label{fig:poll_art_bar}", out.width = "\\textwidth", fig.pos = "!h"----
year_locality_poll %>% 
 ggplot(., aes(x = locality,
               y = no_species,
               fill = group),
       width = 5) +
  geom_bar(position = "stack",
           stat = "identity") +
  scale_fill_nina(name = "Artskategori") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Antall identifiserte pollinatorarter") +
  xlab("Lokalitet")


## -------------------------------------------------------------------------------------------------------------------------
obs <- tbl(con,
           Id(schema = "occurrences",
              table = "observations"))

dung_species <- obs %>%  filter(id_family == "Geotrupidae" | id_family == "Scarabaeidae") %>% select(species_latin_fixed) %>% distinct() %>% pull()

year_locality_dungbeetles <- obs_from_db(subset_species = dung_species,
                                             agg_level = "year_locality") %>% 
  mutate(year = as.factor(year))


## ----dung_art_bar, fig.cap = "Antall identifiserte gjødselbiller per lokal. Fangstene kommer fra 1 malaisefelle per lokalitet.\\label{fig:dung_art_bar}", out.width = "\\textwidth"----
year_locality_dungbeetles %>% 
 ggplot(., aes(x = locality,
               y = no_species,
               fill = year),
       width = 5) +
  geom_bar(position = "stack",
           stat = "identity") +
  scale_fill_nina(name = "År",
                  palette = "darkblue-orange") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Antall identifiserte gjødselbiller") +
  xlab("Lokalitet")


## -------------------------------------------------------------------------------------------------------------------------
beetle_guilds_raw <- tbl(con,
                      Id(schema = "lookup",
                         table = "nmbu_beetle_traits"))

beetle_guilds <- beetle_guilds_raw %>% 
  select(species_latin, 
         feeding_type,
         mean_body_size) %>% 
  collect()

beetle_guilds_subset <- beetle_guilds %>% 
  pull(species_latin)


## -------------------------------------------------------------------------------------------------------------------------

beetle__richness <- obs_from_db(subset_orders = "Coleoptera",
                               subset_species = beetle_guilds_subset,
                               agg_level = "none")

#beetle__richness 


## -------------------------------------------------------------------------------------------------------------------------
beetle_guild_agg <- beetle__richness  %>% 
  left_join(beetle_guilds,
            by = c("species_latin_fixed" = "species_latin"),
            copy = T) %>% 
        collect() %>% 
      group_by(year_locality_id, locality_id, feeding_type) %>% 
      summarise(no_species = n_distinct(species_latin),
                shannon_div = calc_shannon(species_latin)) %>% 
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
               feeding_type,
               no_species,
               shannon_div) %>% 
        arrange(year,
                region_name,
                habitat_type,
                locality) %>% 
  mutate(feeding_type = recode_factor(feeding_type, 
                                      carc = "Åtsel",
                                      detr = "Detrivor",
                                      dung = "Møkk",
                                      herb = "Plantespiser",
                                      mixe = "Blandet",
                                      pred = "Predator"))

#beetle_guild_agg


## ----beetle_feeding_guild_skog, fig.cap = "Fordeling av funne biller i skog fordelt på kjent fødotype. Fangstene kommer fra 1 malaisefelle og 4 vindusfeller per lokalitet. \\label{fig:beetle_feeding_guild_skog}.", out.width = "\\textwidth", fig.pos = "!h"----
beetle_guild_agg %>% 
  filter(habitat_type == "Forest") %>% 
ggplot(.,
       aes(x = locality,
           y = no_species)) +
  geom_bar(aes(fill = feeding_type),
           stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_nina(name = "Fødotype") +
  xlab("Lokalitet") +
  ylab("Antall identifiserte billearter")



## ----beetle_feeding_guild_semi_nat, fig.cap = "Fordeling av funne biller i semi-naturlig mark fordelt på kjent fødotype. Fangstene kommer fra 1 malaisefelle per lokalitet. \\label{fig:beetle_feeding_guild_semi_nat}", out.width = "\\textwidth", fig.pos = "!h"----
beetle_guild_agg %>% 
  filter(habitat_type == "Semi-nat") %>% 
ggplot(.,
       aes(x = locality,
           y = no_species)) +
  geom_bar(aes(fill = feeding_type),
           stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_nina(name = "Fødotype") +
  xlab("Lokalitet") +
  ylab("Antall identifiserte billearter")



## -------------------------------------------------------------------------------------------------------------------------
beetle_size <- beetle__richness  %>% 
  left_join(beetle_guilds,
            by = c("species_latin_fixed" = "species_latin"),
            copy = T) %>% 
  distinct(year_locality_id,
           locality, 
           species_latin_fixed,
           mean_body_size) %>% 
  mutate(habitat_type = ifelse(grepl("Semi-nat", locality), "Semi-nat", "Skog")) %>% 
        collect() 
       


## ----beetle_size_eng, fig.cap = "Fordeling av kroppsstørrelse hos biller, fordelt på lokalitet. \\label{fig:beetle_size_eng}", fig.pos = "!h"----
beetle_size %>% 
 # filter(locality == "Skog_01") %>% 
ggplot(.) +
  geom_density(aes(x = mean_body_size,
           group = locality,
           color = locality)) +
  facet_wrap(~habitat_type,
             ncol = 1) +
  xlab("Kroppsstørrelse (mm)") +
  ylab("Densitet") +
  scale_color_nina(name = "Lokalitet") +
  theme(legend.position = "none")


## -------------------------------------------------------------------------------------------------------------------------
ind_tab <- tibble("Indikator" = c("Totalt artsantall",
                                  "Total biomasse",
                                  "Artsantall av rødlistede arter",
                                  "Artsantall av pollinatorer",
                                  "Artsantall hos sommerfugler (muligen videre delt opp)",
                                  "Artsantall hos sommerfugler spesialiert ift. vertsplante",
                                  "Artsantall av møkkbiller",
                                  "Artsantall av vedlevende biller",
                                  "Forhold i diversitet fluer og parasittveps",
                                  "Forhold i diversitet mellom feeding guilds hos biller"),
                  "Tenkt påvirkning" = c("Generell menneskelig og klimarelatert påvirkning",
                                         "Generell menneskelig og klimarelatert påvirkning",
                                         "Generell menneskelig og klimarelatert påvirkning",
                                         "Intensivt jordbruk, tap av tradisjonelle driftsformer, fragmentering.",
                                         "Intensivt jordbruk, tap av tradisjonelle driftsformer, fragmentering",
                                         "Intensivt jordbruk, tap av tradisjonelle driftsformer, fragmentering",
                                         "Intensivt jordbruk, tap av tradisjonelle driftsformer, fragmentering",
                                         "Intensivt skogsbruk",
                                         "Minsket mengde beite og forandret landskapsstruktur",
                                         "Størning av ekosystem og endret arealbruk påvirker trofisk struktur"),
                  "Status" = c("Vi mangler data for referansetilstand",
                               "Vi mangler data for referansetilstand",
                               "Vi mangler data for referansetilstand",
                               "Vi mangler data for referansetilstand",
                               "Vi mangler data for referansetilstand",
                               "Vi mangler data for referansetilstand",
                               "Lite data fra overvåkingen. Vi mangler data for referansetilstand",
                               "Vi mangler data for referansetilstand",
                               "Vi mangler data for referansetilstand",
                               "Vi mangler data for referansetilstand")
                  )




## ----ind_tab, results = "asis"--------------------------------------------------------------------------------------------
xtab <- xtable(ind_tab, 
             caption = "Tentativ liste over tenkbare indikatorer som overvåkingsprosjektet kan gi data for. \\label{tab:ind_tab}")


align(xtab) <- "p{4cm}p{4cm}p{4cm}p{4cm}"

print.xtable(xtab, sanitize.text.function=function(x){x})



## -------------------------------------------------------------------------------------------------------------------------
oppd_budj <- tibble("Aktivitet" = c("Forarbeid",
                                    "Feltarbeid",
                                    "Etterarbeid",
                                    "Rapportering",
                                    "Driftsmidler",
                                    "Sum kostnad"),
                    "Kostnader" = c(750,
                                    1500,
                                    1250,
                                    550,
                                    850,
                                    4900)
) %>% 
  mutate(Kostnader = as.integer(Kostnader))




## ----oppd_budj_tab, results = "asis"--------------------------------------------------------------------------------------

xtab <- xtable(oppd_budj, 
             caption = "Oppdatert budsjett for overvåking av semi-naturlig mark og skogmark i Østlandet. Priser i KNOK. \\label{tab:oppd_budj_tab}",
             digits = 0)


align(xtab) <- "p{2cm}p{3cm}p{2cm}"

print.xtable(xtab, sanitize.text.function=function(x){x})




## ----tenkbar_fremtid, fig.cap = "Skissert utvikling med fullt omløp i skog og semi-naturlig mark i hele landet. \\label{fig:tenkbar_fremtid}", out.width = "0.8\\textwidth"----
ggplot(norge) +
  geom_sf() +
  geom_sf(data = rand_500_ssb,
          aes(fill = habitat_type,
              color = habitat_type),
          lwd = 2) +
  scale_color_nina(name = "Habitat",
                   palette = "orange-green") +
  scale_fill_nina(name = "Habitat",
                   palette = "orange-green")
  


## -------------------------------------------------------------------------------------------------------------------------
utvid_tab = tibble("Aktivitet" = c("Forarbeid",
                                    "Feltarbeid",
                                    "Etterarbeid",
                                    "Rapportering",
                                    "Driftsmidler",
                                    "Sum kostnad"),
                    "Kostnad per del" = c(600,
                                    1500,
                                    900,
                                    300,
                                    900,
                                    4200),
                   "Antall deler" = rep(5, 6),
                   "Total kostnad for utvidelse" = c(2400,
                                                       6000,
                                                       3600,
                                                       1200,
                                                       3600,
                                                       16800)
) 



## ----utvid_budj_tab, results = "asis"-------------------------------------------------------------------------------------

xtab <- xtable(utvid_tab, 
             caption = "Estimert budsjett for en geografisk utvidelse av overvåking i semi-naturlig mark og skog i hele landet (kostnad utover Østlandet). Tallen er basert på erfaringene i Trøndelag. Feltkostnader for områdene «Vest» og «Nord» vil være høyere og vil avhenge av strategiske valg. Priser i KNOK. \\label{tab:utvid_budj_tab}",
             digits = 0)


align(xtab) <- "p{2cm}p{3cm}p{2cm}p{2cm}p{2cm}"

print.xtable(xtab, sanitize.text.function=function(x){x})




## ----prio-----------------------------------------------------------------------------------------------------------------
prio_tab <- tibble(Prioritering = c(1, 2, 3, 4, 4),
                   Del = c("Inkludere økosystemet skog i Trøndelag",
                           "Etablering av referanseområder i de overvåkede områdene",
                           "Utvidelse til region Sør (Agder, Rogaland)",
                           "Utvidelse til region Vest (Vestland, Møre og Romsdal)",
                           "Utvidelse til region Nord (Nordland, Troms og Finnmark)"),
                   Merknader = c("Praktisk og kostnadsmessig mest gunstig, og gir balansert design.",
                                 "Utvidelse til referanseområder bør kunne følge nåværende rutiner og kostnadsanslag, men vil kreve noe utvikling av utvalgskriterier",
                                 "Utvidelse i dette området bør kunne følge nåværende rutiner og kostnadsanslag.",
                                 "Utvidelse til denne regionen kan gi høyere feltkostnader på grunn av avstand og veinett. Kan være nødvendig med et spleiselag på feltarbeidet med andre aktører. Utvidelse til disse områdene bør kommuniseres i god tid da det krever mer planlegging og utøkning av kapasitet.", 
                                 "Utvidelse til denne regionen kan gi høyere feltkostnader på grunn av avstand og veinett. Kan være nødvendig med et spleiselag på feltarbeidet med andre aktører. Utvidelse til disse områdene bør kommuniseres i god tid da det krever mer planlegging og utøkning av kapasitet."))


## ----prio_tab, results = "asis"-------------------------------------------------------------------------------------------
xtab <- xtable(prio_tab, 
             caption = "Grov skisse av videre utvikling i prosjektet. Prioriteringen av regioner for utvidelse bør vurderes sammen med oppdragsgiver. \\label{tab:prio_tab}",
             digits = 0)


align(xtab) <- "p{2cm}p{2cm}p{4cm}p{6cm}"

print.xtable(xtab, sanitize.text.function=function(x){x})




## -------------------------------------------------------------------------------------------------------------------------
pot_alien_spec <- tbl(con,
                  Id(schema = "lookup",
                     table = "pot_alien_species")) %>% 
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
  left_join(obs,
            by = c("Artsnavn" = "species_latin_fixed"),
            copy = T) %>% 
  select(Order = "id_order",
         Artsnavn,
         "Listet i Saul",
         "Listet i Daisie") %>% 
  distinct()

fennoscand_spec <- tbl(con,
                  Id(schema = "lookup",
                     table = "fennoscand_species"))  %>% 
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
  left_join(obs,
            by = c("Artsnavn" = "species_latin_fixed"),
            copy = T) %>% 
  select(Order = "id_order",
         Artsnavn,
         "Listet i Saul",
         "Listet i Daisie") %>% 
  distinct()



## -------------------------------------------------------------------------------------------------------------------------

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



## ---- results = "asis"----------------------------------------------------------------------------------------------------
print(xtable(alien_obs))


## ---- results = "asis"----------------------------------------------------------------------------------------------------
print(xtable(fennoscand_spec),
      size = "footnotesize",
      floating = F,
      floating.environment = "table*",
      tabular.environment = "longtable")


## ---- results = "asis"----------------------------------------------------------------------------------------------------
print(xtable(pot_alien_spec), 
      size = "footnotesize",
      floating = F,
      floating.environment = "table*",
      tabular.environment = "longtable")

