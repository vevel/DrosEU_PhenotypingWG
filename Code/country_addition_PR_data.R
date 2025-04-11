
rm(list = ls())

##### set wd

setwd("~/Work/UNIFR/GitHub/DrosEU_PhenotypingWG/")


##### load packages and source functions

library(tidyverse)

##### load data


pops <- readRDS("InfoTables/Droseu_Populations.rds")$by_lat %>%
    select(Population, Country, Latitude, Longitude, Altitude)

pr <- read.csv("./Data/MasterSheets_Apr23_git/PR_MasterSheet_Apr23.csv")

pr_up <- inner_join(pr, pops, by = "Population")

write.csv(pr_up, "./Data/MasterSheets_Mar25_git/PR_MasterSheet_Mar25.csv", row.names = FALSE)
