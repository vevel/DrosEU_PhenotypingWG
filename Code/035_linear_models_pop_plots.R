

# Generate plots for each lab and trait
# per population with pairwise compararisons (tukey) and
# per population split by lab phenotyping batches

# This is to replace the code that used to be in the .rmd file when those plots
# were produced while knitting

# This code only deals with tha trait and lab combinations correponding to
# the 97 linear models used in the manuscript


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

##### load tukey stats
tuk_lmers <- readRDS("LinearModelsPop/all_lmers_pop_tukey_list.rds")
tuk_glmers <- readRDS("LinearModelsPop/all_glmers_pop_tukey_list.rds")

##### create output directory
batch_pops_plots_dir <- "LinearModelsPopPlots"
dir.create(batch_pops_plots_dir, showWarnings = FALSE)




##### remove traits for which we do not need plots since we did not run models ourselves

traits_to_rm <- c("cets", "lsm")
droseu <- droseu[!(names(droseu) %in% traits_to_rm)]


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
  dplyr::select(Trait_handle, Trait, Trait_name_raw, Legend, Title) %>%
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

droseu_lab_trait_sex <- group_by(droseu_long, Trait, Trait_name, Lab, Sex) %>%
  group_split()


##### get all lab, sex, trait combinations to retrieve tukey statistics

lab_trait_sex <- list()
for (i in seq_along(droseu_lab_trait_sex)) {
  lab_trait_sex[[i]] <- paste(
    droseu_lab_trait_sex[[i]]$Trait,
    droseu_lab_trait_sex[[i]]$Sex,
    droseu_lab_trait_sex[[i]]$Lab,
    sep = "_"
  )[1]
}
lab_trait_sex <- unlist(lab_trait_sex)
lab_trait_sex <- gsub("_NA", "", lab_trait_sex)


names(droseu_lab_trait_sex) <- lab_trait_sex



##### prepare tukey statistics, exclude models that are not used
# not very elegant but I am bit rusty

tuk_all_models <- c(tuk_lmers, tuk_glmers)

models_to_rm <- c(
  "Dia_Flatt_lm_pop",
  "Dia_Bergland_lmer_pop",
  "Dia_Schlotterer_lmer_pop"
)

tuk_models <- tuk_all_models[!(names(tuk_all_models) %in% models_to_rm)]

tuk_letters <- lapply(tuk_models, function(x) x$letters)

trait_sex_lab <- str_match(names(tuk_letters), "(.*)_(.*?)_(.*?)$")[, -1][, 1]
trait_sex_lab <- sub("Pgm_T4_", "Pgm_T4_F_", trait_sex_lab)
trait_sex_lab <- sub("Pgm_T5_", "Pgm_T5_F_", trait_sex_lab)
trait_sex_lab <- sub("Pgm_T6_", "Pgm_T6_F_", trait_sex_lab)
trait_sex_lab <- sub("Pgm_Total_", "Pgm_Total_F_", trait_sex_lab)
trait_sex_lab <- sub("Dia_", "Dia_F_", trait_sex_lab)
trait_sex_lab <- sub("Fec_", "Fec_F_", trait_sex_lab)
trait_sex_lab <- sub("LS_F_Flatt", "LSP_F_Flatt", trait_sex_lab)
trait_sex_lab <- sub("LS_M_Flatt", "LSP_M_Flatt", trait_sex_lab)
trait_sex_lab <- sub("LS_F_Parsch", "LSL_F_Parsch", trait_sex_lab)
trait_sex_lab <- sub("LS_M_Parsch", "LSL_M_Parsch", trait_sex_lab)
trait_sex_lab <- sub("LS_F_Pasyukova", "LSL_F_Pasyukova", trait_sex_lab)
trait_sex_lab <- sub("LS_M_Pasyukova", "LSL_M_Pasyukova", trait_sex_lab)
trait_sex_lab <- sub("_Tauber", "_M_Tauber", trait_sex_lab)
trait_sex_lab <- sub("NDlog2", "ND", trait_sex_lab)

names(tuk_letters) <- trait_sex_lab




##### populations as main groups, with tukey results

for (i in seq_along(droseu_lab_trait_sex)) {

  n <- names(droseu_lab_trait_sex)[i]
  d <- droseu_lab_trait_sex[[i]]
  l <- bind_rows(tuk_letters[names(tuk_letters) == n])

  min_y <- min(d$Value, na.rm = TRUE)
  max_y <- max(d$Value, na.rm = TRUE)
  range <- max_y - min_y
  tuk_y <- max_y + range * 0.07
  lim_y <- max_y + range * 0.1


  title_text <- paste(d$Title[1], d$Sex[1], sep = " - ")
  title_text <- sub(" - F", " - Females", title_text)
  title_text <- sub(" - M", " - Males", title_text)
  title_text <- sub(" - NA", "", title_text)

  out_file_pdf <- paste(
    "p_", d$Trait[1], "_", d$Sex[1], "_pop_", d$Lab[1], ".pdf",
    sep = ""
  )
  out_file_pdf <- sub("_NA", "", out_file_pdf)

  p <- ggplot(d, aes(x = Population, y = Value, fill = Population)) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(0.9)) +
  labs(
      x = "Population",
      y = d$Legend[1],
      title = title_text
  ) +
  scale_y_continuous(limits = c(NA, lim_y)) +
  droseu_fill_scale_pop +
  geom_text(
    data = l, aes(x = Population, label = cld),
    y = tuk_y, size = 4
  ) +
  theme_classic(15)

  ggsave(
    p,
    filename = file.path(batch_pops_plots_dir, out_file_pdf),
    width = 7, height = 4.5
  )
}

##### populations as main groups, with batches as subgroups


for (i in seq_along(droseu_lab_trait_sex)) {

  n <- names(droseu_lab_trait_sex)[i]
  d <- droseu_lab_trait_sex[[i]]

  min_y <- min(d$Value, na.rm = TRUE)
  max_y <- max(d$Value, na.rm = TRUE)
  range <- max_y - min_y
  lim_y <- max_y + range * 0.1


  title_text <- paste(d$Title[1], d$Sex[1], sep = " - ")
  title_text <- sub(" - F", " - Females", title_text)
  title_text <- sub(" - M", " - Males", title_text)
  title_text <- sub(" - NA", "", title_text)

  out_file_pdf <- paste(
    "p_", d$Trait[1], "_", d$Sex[1], "_pop_batch_", d$Lab[1], ".pdf",
    sep = ""
  )
  out_file_pdf <- sub("_NA", "", out_file_pdf)

  p <- ggplot(d, aes(x = Population, y = Value, fill = Batch)) +
    geom_boxplot(outlier.shape = NA, position = position_dodge(0.9)) +
  labs(
      x = "Population",
      y = d$Legend[1],
      title = title_text
  ) +
  scale_y_continuous(limits = c(NA, lim_y)) +
  scale_fill_grey(start = 1, end = 0) +
  theme_classic(15)

  ggsave(
    p,
    filename = file.path(batch_pops_plots_dir, out_file_pdf),
    width = 7, height = 4.5
  )
}