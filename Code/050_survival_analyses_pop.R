
###############################################################
###################### SURVIVAL ANALYSES ######################
###############################################################



# Running survival analyses takes too long and clogs the computer / knitting. 
# Thus, all were done externally and output files are shown in .Rmd file. 




##### clean workspace
rm(list = ls())

##### libraries
library(tidyverse)
library(ggpubr)
library(coxme)
library(survival)
library(survminer)




##### set working directory
setwd("~/Work/UNIFR/GitHub/DrosEU_PhenotypingWG/")


##### load data
droseu <- readRDS("Data/droseu_master_list_2025-03-24.rds")

##### create output directory
surv_dir <- "SurvivalAnalyses"
dir.create(surv_dir, showWarnings = F) 

##### source functions
source("Code/functions.R")






############# STARVATION RESISTANCE ############# 

# create output directory
out_dir <- "Starvation"
dir.create(file.path(surv_dir, out_dir), showWarnings = F) 

droseu$sr <- droseu$sr %>% mutate(Censor = 1)

# initialize output list
SR_coxme_pop <- list()

# Gonzalez
SR_coxme_pop$SR_F_Gonzalez_coxme_pop <- coxme(Surv(AgeAtDeath_hours, Censor) ~ Population + (1|Batch) + (1|Population/Line), data = filter(droseu$sr, Supervisor.PI == "Gonzalez", Sex == "F"))

SR_coxme_pop$SR_M_Gonzalez_coxme_pop <- coxme(Surv(AgeAtDeath_hours, Censor) ~ Population + (1|Batch) + (1|Population/Line), data = filter(droseu$sr, Supervisor.PI == "Gonzalez", Sex == "M"))

# Onder
SR_coxme_pop$SR_F_Onder_coxme_pop <- coxme(Surv(AgeAtDeath_hours, Censor) ~ Population + (1|Batch) + (1|Population/Line), data = filter(droseu$sr, Supervisor.PI == "Onder", Sex == "F"))

SR_coxme_pop$SR_M_Onder_coxme_pop <- coxme(Surv(AgeAtDeath_hours, Censor) ~ Population + (1|Batch) + (1|Population/Line), data = filter(droseu$sr, Supervisor.PI == "Onder", Sex == "M"))

# Onder
SR_coxme_pop$SR_F_Pasyukova_coxme_pop <- coxme(Surv(AgeAtDeath_hours, Censor) ~ Population + (1|Batch) + (1|Population/Line), data = filter(droseu$sr, Supervisor.PI == "Pasyukova", Sex == "F"))

SR_coxme_pop$SR_M_Pasyukova_coxme_pop <- coxme(Surv(AgeAtDeath_hours, Censor) ~ Population + (1|Batch) + (1|Population/Line), data = filter(droseu$sr, Supervisor.PI == "Pasyukova", Sex == "M"))

# save output list
saveRDS(SR_coxme_pop, file = file.path(surv_dir, out_dir, "SR_coxmes_pop.rds"))









############# LIFESPAN ############# 

# create output directory
out_dir <- "Lifespan"
dir.create(file.path(surv_dir, out_dir), showWarnings = F) 

droseu$lsl <- droseu$lsl %>% mutate(Censor = ifelse(Censor == 0, 1, 0))
droseu$lsp <- droseu$lsp %>% mutate(Censor = ifelse(Censor == 0, 1, 0))

# initialize output list
LS_coxme_pop <- list()

# Flatt
LS_coxme_pop$LS_P_F_Flatt_coxme_pop <- coxme(Surv(LSP_AgeAtDeath_days, Censor) ~ Population + (1|Population/ReplicateCage), data = filter(droseu$lsp, Supervisor.PI == "Flatt", Sex == "F"))

LS_coxme_pop$LS_P_M_Flatt_coxme_pop <- coxme(Surv(LSP_AgeAtDeath_days, Censor) ~ Population + (1|Population/ReplicateCage), data = filter(droseu$lsp, Supervisor.PI == "Flatt", Sex == "M"))

# Parsch
LS_coxme_pop$LS_L_F_Parsch_coxme_pop <- coxme(Surv(LSL_AgeAtDeath_days, Censor) ~ Population + (1|Batch) + (1|Population/Line), data = filter(droseu$lsl, Supervisor.PI == "Parsch", Sex == "F"))

LS_coxme_pop$LS_L_M_Parsch_coxme_pop <- coxme(Surv(LSL_AgeAtDeath_days, Censor) ~ Population + (1|Batch) + (1|Population/Line), data = filter(droseu$lsl, Supervisor.PI == "Parsch", Sex == "M"))

# Pasyukova
LS_coxme_pop$LS_L_F_Pasyukova_coxme_pop <- coxme(Surv(LSL_AgeAtDeath_days, Censor) ~ Population + (1|Batch) + (1|Population/Line), data = filter(droseu$lsl, Supervisor.PI == "Pasyukova", Sex == "F"))

LS_coxme_pop$LS_L_M_Pasyukova_coxme_pop <- coxme(Surv(LSL_AgeAtDeath_days, Censor) ~ Population + (1|Batch) + (1|Population/Line), data = filter(droseu$lsl, Supervisor.PI == "Pasyukova", Sex == "M"))

# save output list
saveRDS(LS_coxme_pop, file = file.path(surv_dir, out_dir, "LS_coxmes_pop.rds"))





############# HEATSHOCK ############# 

# create output directory
out_dir <- "Heatshock"
dir.create(file.path(surv_dir, out_dir), showWarnings = F) 

droseu$hsm <- droseu$hsm %>% mutate(Censor = ifelse(Censor == 0, 1, 0))

# initialize output list
HSM_coxme_pop <- list()

# Parsch
HSM_coxme_pop$HSM_F_Parsch_coxme_pop <- coxme(Surv(TimeDeath_min, Censor) ~ Population + (1|Batch) + (1|Population/Line) , data = filter(droseu$hsm, Supervisor.PI == "Parsch", Sex == "F"))

HSM_coxme_pop$HSM_M_Parsch_coxme_pop <- coxme(Surv(TimeDeath_min, Censor) ~ Population + (1|Batch) + (1|Population/Line) , data = filter(droseu$hsm, Supervisor.PI == "Parsch", Sex == "M"))

# Vieira
HSM_coxme_pop$HSM_F_Vieira_coxme_pop <- coxme(Surv(TimeDeath_min, Censor) ~ Population + (1|Batch) + (1|Population/Line) , data = filter(droseu$hsm, Supervisor.PI == "Vieira", Sex == "F"))

HSM_coxme_pop$HSM_M_Vieira_coxme_pop <- coxme(Surv(TimeDeath_min, Censor) ~ Population + (1|Batch) + (1|Population/Line) , data = filter(droseu$hsm, Supervisor.PI == "Vieira", Sex == "M"))

# save output list
saveRDS(HSM_coxme_pop, file = file.path(surv_dir, out_dir, "HSM_coxmes_pop.rds"))





############# CCRT ############# 

# create output directory
out_dir <- "ChillComa"
dir.create(file.path(surv_dir, out_dir), showWarnings = F) 

droseu$ccrt <- droseu$ccrt %>% mutate(Censor = ifelse(Censor == 0, 1, 0))

# initialize output list

CCRT_coxme_pop <- list()

# Vieira
CCRT_coxme_pop$CCRT_F_Vieira_coxme_pop <- coxme(Surv(CCRT_seconds, Censor) ~ Population + (1|Batch) + (1|Population/Line) , data = filter(droseu$ccrt, Supervisor.PI == "Vieira", Sex == "F"))

CCRT_coxme_pop$CCRT_M_Vieira_coxme_pop <- coxme(Surv(CCRT_seconds, Censor) ~ Population + (1|Batch) + (1|Population/Line) , data = filter(droseu$ccrt, Supervisor.PI == "Vieira", Sex == "M"))

# Mensch
CCRT_coxme_pop$CCRT_F_Mensch_coxme_pop <- coxme(Surv(CCRT_seconds, Censor) ~ Population + (1|Batch) + (1|Population/Line) , data = filter(droseu$ccrt, Supervisor.PI == "Mensch", Sex == "F"))

CCRT_coxme_pop$CCRT_M_Mensch_coxme_pop <- coxme(Surv(CCRT_seconds, Censor) ~ Population + (1|Batch) + (1|Population/Line) , data = filter(droseu$ccrt, Supervisor.PI == "Mensch", Sex == "M"))


# save output list
saveRDS(CCRT_coxme_pop, file = file.path(surv_dir, out_dir, "CCRT_coxmes_pop.rds"))






######### combine all lmers into one list ###########


all_coxmes_pop <- rdsBatchReaderToList(path = surv_dir, recursive = T, full.names = T, pattern = "coxme_pop.rds")

all_coxmes_pop <- unlist(all_coxmes_pop, recursive=FALSE)
names(all_coxmes_pop) <- str_split(names(all_coxmes_pop), "\\.", simplify = T)[,2]

saveRDS(all_coxmes_pop, file = file.path(surv_dir, "all_coxmes_pop_list.rds"))




############# output all lmers summaries and anovas as global lists ############# 

all_coxmes_pop_summary <- lapply(all_coxmes_pop, summary)
saveRDS(all_coxmes_pop_summary, file = file.path(surv_dir, "all_coxmes_pop_summary_list.rds"))


all_coxmes_pop_anova <- lapply(all_coxmes_pop, anova)
saveRDS(all_coxmes_pop_anova, file = file.path(surv_dir, "all_coxmes_pop_anova_list.rds"))





######### plot linear models residuals per trait

# list models outputs to keep the directory structure
models <- list.files(path = surv_dir, recursive = T, full.names = T, pattern = "_coxmes_pop.rds")

######### output all models summaries, anovas and tukeys by trait

# loop over models to get statistics
for (m in 1:length(models)){
  outputModelsStats(models[m])
}



######### output all models summaries, anovas and tukeys by trait and lab

# loop over models to get statistics
for (m in 1:length(models)){
  outputModelsStatsLab(models[m])
}





######### output all models P values

all_coxmes_pop_anova <- readRDS(file.path(surv_dir, "all_coxmes_pop_anova_list.rds"))

pvalues <- combinePValues(all_coxmes_pop_anova)

saveRDS(pvalues, file = file.path(surv_dir, "all_coxmes_pop_pvalues.rds"))
write.csv(pvalues, file = file.path(surv_dir, "all_coxmes_pop_pvalues.csv"), row.names = F)







