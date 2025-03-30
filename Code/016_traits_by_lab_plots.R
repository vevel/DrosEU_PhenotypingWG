

# Generate plots for each trait
# per population split by labs

# This is to replace the code that used to be in the .rmd file when those plots
# were produced while knitting

# This code only deals with the traits that ED and EK analysed

##### clean workspace
rm(list = ls())

##### libraries
library(tidyverse)

##### set working directory
setwd("~/Work/UNIFR/GitHub/DrosEU_PhenotypingWG/")

##### source functions
source("Code/functions.R")

##### load latest data
droseu <- readRDS("Data/droseu_master_list_2025-03-24.rds")
pops <- readRDS("InfoTables/DrosEU_Populations.rds")


##### remove traits for which we do not need plots since we did not analyse
# them ourselves
# also filter values for PR_HR that are above 1

traits_to_rm <- c("cets", "lsm")
droseu <- droseu[!(names(droseu) %in% traits_to_rm)]
droseu$pr <- filter(droseu$pr, HostResistance <= 1 & HostResistance >= 0)


##### define trait and common variables

trait_var <- c(
  "CCRT_seconds", "CSM_PropDead_ED", "Prop_Max_Stage9", "DT_EggAdult", "DT_EggPupa",
  "DW_micrograms", "NumberOfAdultsEclosed", "TimeDeath_min", "Period",
  "CircPhase", "AbsPhase", "ND", "Activity", "LSL_AgeAtDeath_days",
  "LSP_AgeAtDeath_days", "PercT4", "PercT5", "PercT6",
  "TotalPerc", "AgeAtDeath_hours", "TL_micrometers",
  "ProportionEggtoAdultSurvival", "CentroidSizeLeft_micrometers",
  "CentroidSizeRight_micrometers", "HostResistance"
)

common_var <- c("Supervisor.PI", "Batch", "Population", "Line", "Sex")


##### add extra columns for easier batch processing and split data by single lab, trait and sex

droseu_long <- droseu
for (i in seq_along(droseu_long)){
  d <- droseu_long[[i]]
  if (!"Line" %in% colnames(d)) d$Line <- NA
  if (!"Batch" %in% colnames(d)) d$Batch <- NA
  if (!"Sex" %in% colnames(d)) d$Sex <- NA
  d <- d[, colnames(d) %in% c(common_var, trait_var)]
  d <- pivot_longer(d, cols = -all_of(common_var), names_to = "Trait_name", values_to = "Value")
  d$Trait <- names(droseu_long)[i]
  d <- split(d, d$Trait_name)
  droseu_long[[i]] <- d
}

droseu_long <- bind_rows(unlist(droseu_long, recursive = FALSE)) %>%
  rename(Trait_handle = Trait, Lab = Supervisor.PI)

trait_names <- read.csv("InfoTables/trait_names.csv") %>%
  dplyr::select(Trait_handle, Trait, Trait_name_raw, Legend, Title, Directory) %>%
  rename(Trait_name = Trait_name_raw)

droseu_long <- inner_join(
  droseu_long,
  trait_names
)

droseu_long <- inner_join(
  droseu_long,
  pops$by_lat
) %>%
  mutate(Population = factor(Population, levels = pops$by_lat$Population))

droseu_trait_sex <- group_by(droseu_long, Trait, Trait_name, Sex) %>%
  group_split()



##### populations as main groups, with labs as subgroups

for (i in seq_along(droseu_trait_sex)) {

  n <- names(droseu_trait_sex)[i]
  d <- droseu_trait_sex[[i]]

  title_text <- paste(d$Title[1], d$Sex[1], sep = " - ")
  title_text <- sub(" - F", " - Females", title_text)
  title_text <- sub(" - M", " - Males", title_text)
  title_text <- sub(" - NA", "", title_text)

  out_file_png <- paste(d$Directory[1], "/",
    "p_", d$Trait[1], "_", d$Sex[1], "_pop_lab", ".png",
    sep = ""
  )
  out_file_png <- sub("_NA", "", out_file_png)

  p <- ggplot(d, aes(x = Population, y = Value, fill = Lab)) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(0.9)) +
  labs(
      x = "Population",
      y = d$Legend[1],
      title = title_text
  ) +
  scale_fill_grey(start = 1, end = 0) +
  theme_classic(15)

  ggsave(
    p,
    filename = file.path(out_file_png),
    width = 7, height = 4.5
  )
}
