


#########################################################################
######################  RUN ALL LMMs FOR Altitude ######################
#########################################################################


### NOT RUN FOR PR_HR AS THE DATA CAME AFTER wE DECIDED NOT TO USE DO THE BELOW ANALZSES


##### clean workspace
rm(list = ls())

##### libraries
library(tidyverse)
library(lme4)
library(lsmeans)
library(afex)
library(multcomp)


##### set working directory
setwd("~/Work/UNIFR/GitHub/DrosEU_PhenotypingWG/")


##### source functions
source("Code/functions.R")

##### load data
droseu <- readRDS("Data/droseu_master_list_2022-05-02.rds")

##### create output directory
lmer_dir <- "LinearModelsAlt"
dir.create(lmer_dir, showWarnings = F) 





############# VIABILITY #############

# create output directory
out_dir <- "Viability"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 

# initialize output list
Via_lmers_alt <- list()

#### Gibert Lab
Via_lmers_alt$Via_Gibert_lmer_alt <- lmer(ProportionEggtoAdultSurvival ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$via, Supervisor.PI == "Gibert"))

#### Grath Lab
Via_lmers_alt$Via_Grath_lmer_alt <- lmer(ProportionEggtoAdultSurvival ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$via, Supervisor.PI == "Grath"))

#### Hoedjes Lab, not converging, removed Population
Via_lmers_alt$Via_Hoedjes_lmer_alt <- lmer(ProportionEggtoAdultSurvival ~ Altitude + (1|Line:Population), data = filter(droseu$via, Supervisor.PI == "Hoedjes"))

#### Schmidt Lab
Via_lmers_alt$Via_Schmidt_lmer_alt <- lmer(ProportionEggtoAdultSurvival ~ Altitude + (1|Population), data = filter(droseu$via, Supervisor.PI == "Schmidt"))

#### Stamenkovic-Radak Lab
Via_lmers_alt$Via_StamenkovicRadak_lmer_alt <- lmer(ProportionEggtoAdultSurvival ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$via, Supervisor.PI == "StamenkovicRadak"))

#### Zwaan Lab
Via_lmers_alt$Via_Zwaan_lmer_alt <- lmer(ProportionEggtoAdultSurvival ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$via, Supervisor.PI == "Zwaan"))

# save output list
saveRDS(Via_lmers_alt, file = file.path(lmer_dir, out_dir, "Via_lmers_alt.rds"))








############# DEVELOPMENT TIME #############

# create output directory
out_dir <- "DevelopmentTime"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 


### Egg-to-Pupae

# initialize output list
DT_lmers_alt <- list()

##### Schmidt Lab
DT_lmers_alt$DT_P_Schmidt_lmer_alt <- lmer(DT_EggPupa ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$dtp, Supervisor.PI == "Schmidt"))

### Egg-to-adult developmental time

### Females

##### Gibert Lab
DT_lmers_alt$DT_A_F_Gibert_lmer_alt <- lmer(DT_EggAdult ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$dta, Supervisor.PI == "Gibert" & Sex == "F"))

##### Grath Lab
DT_lmers_alt$DT_A_F_Grath_lmer_alt <- lmer(DT_EggAdult ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$dta, Supervisor.PI == "Grath" & Sex == "F"))

##### Hoedjes Lab
DT_lmers_alt$DT_A_F_Hoedjes_lmer_alt <- lmer(DT_EggAdult ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$dta, Supervisor.PI == "Hoedjes" & Sex == "F"))

##### Schmidt Lab
DT_lmers_alt$DT_A_F_Schmidt_lmer_alt <- lmer(DT_EggAdult ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$dta, Supervisor.PI == "Schmidt" & Sex == "F"))

##### Stamenkovic-Radak Lab
DT_lmers_alt$DT_A_F_StamenkovicRadak_lmer_alt <- lmer(DT_EggAdult ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$dta, Supervisor.PI == "StamenkovicRadak" & Sex == "F"))

##### Zwaan Lab
DT_lmers_alt$DT_A_F_Zwaan_lmer_alt <- lmer(DT_EggAdult ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$dta, Supervisor.PI == "Zwaan" & Sex == "F"))


### Males

##### Gibert Lab
DT_lmers_alt$DT_A_M_Gibert_lmer_alt <- lmer(DT_EggAdult ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$dta, Supervisor.PI == "Gibert" & Sex == "M"))

##### Grath Lab
DT_lmers_alt$DT_A_M_Grath_lmer_alt <- lmer(DT_EggAdult ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$dta, Supervisor.PI == "Grath" & Sex == "M"))

##### Hoedjes Lab
DT_lmers_alt$DT_A_M_Hoedjes_lmer_alt <- lmer(DT_EggAdult ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$dta, Supervisor.PI == "Hoedjes" & Sex == "M"))

##### Schmidt Lab
DT_lmers_alt$DT_A_M_Schmidt_lmer_alt <- lmer(DT_EggAdult ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$dta, Supervisor.PI == "Schmidt" & Sex == "M"))

##### Stamenkovic-Radak Lab
DT_lmers_alt$DT_A_M_StamenkovicRadak_lmer_alt <- lmer(DT_EggAdult ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$dta, Supervisor.PI == "StamenkovicRadak" & Sex == "M"))

##### Zwaan Lab
DT_lmers_alt$DT_A_M_Zwaan_lmer_alt <- lmer(DT_EggAdult ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$dta, Supervisor.PI == "Zwaan" & Sex == "M"))

# save output list
saveRDS(DT_lmers_alt, file = file.path(lmer_dir, out_dir, "DT_lmers_alt.rds"))









############# DRY WEIGHT #############

# create output directory
out_dir <- "DryWeight"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 

# initialize output list
DW_lmers_alt <- list()

## Females

##### Colinet Lab
DW_lmers_alt$DW_F_Colinet_lmer_alt <- lmer(DW_micrograms ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$dw, Supervisor.PI == "Colinet" & Sex == "F"))

##### Hoedjes Lab
DW_lmers_alt$DW_F_Hoedjes_lmer_alt <- lmer(DW_micrograms ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$dw, Supervisor.PI == "Hoedjes" & Sex == "F"))

##### Onder Lab
DW_lmers_alt$DW_F_Onder_lmer_alt <- lmer(DW_micrograms ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$dw, Supervisor.PI == "Onder" & Sex == "F"))


## Males

##### Colinet Lab
DW_lmers_alt$DW_M_Colinet_lmer_alt <- lmer(DW_micrograms ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$dw, Supervisor.PI == "Colinet" & Sex == "M"))

##### Hoedjes Lab, singular fit removed Population
DW_lmers_alt$DW_M_Hoedjes_lmer_alt <- lmer(DW_micrograms ~ Altitude + (1|Line:Population), data = filter(droseu$dw, Supervisor.PI == "Hoedjes" & Sex == "M"))

##### Onder Lab
DW_lmers_alt$DW_M_Onder_lmer_alt <- lmer(DW_micrograms ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$dw, Supervisor.PI == "Onder" & Sex == "M"))

# save output list
saveRDS(DW_lmers_alt, file = file.path(lmer_dir, out_dir, "DW_lmers_alt.rds"))











############# THORAX LENGTH #############

# create output directory
out_dir <- "ThoraxLength"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 

# initialize output list
TL_lmers_alt <- list()

## Females

##### Kozeretska Lab
TL_lmers_alt$TL_F_Kozeretska_lmer_alt <- lmer(TL_micrometers ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$tl, Supervisor.PI == 'Kozeretska' & Sex == "F"))

##### Posnien Lab
TL_lmers_alt$TL_F_Posnien_lmer_alt <- lmer(TL_micrometers ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$tl, Supervisor.PI == 'Posnien' & Sex == "F"))

##### Ritchie Lab, singular fit, removed Population
TL_lmers_alt$TL_F_Ritchie_lmer_alt <- lmer(TL_micrometers ~ Altitude + (1|Line:Population), data = filter(droseu$tl, Supervisor.PI == 'Ritchie' & Sex == "F"))

##### Schmidt Lab
TL_lmers_alt$TL_F_Schmidt_lmer_alt <- lmer(TL_micrometers ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$tl, Supervisor.PI == 'Schmidt' & Sex == "F"))


## Males

##### Kozeretska Lab
TL_lmers_alt$TL_M_Kozeretska_lmer_alt <- lmer(TL_micrometers ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$tl, Supervisor.PI == 'Kozeretska' & Sex == "M"))

##### Posnien Lab
TL_lmers_alt$TL_M_Posnien_lmer_alt <- lmer(TL_micrometers ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$tl, Supervisor.PI == 'Posnien' & Sex == "M"))

##### Ritchie Lab, singular fit, removed Population
TL_lmers_alt$TL_M_Ritchie_lmer_alt <- lmer(TL_micrometers ~ Altitude + (1|Line:Population), data = filter(droseu$tl, Supervisor.PI == 'Ritchie' & Sex == "M"))

# save output list
saveRDS(TL_lmers_alt, file = file.path(lmer_dir, out_dir, "TL_lmers_alt.rds"))







############# WING AREA ############# 

# create output directory
out_dir <- "WingArea"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 

# initialize output list
WA_lmers_alt <- list()


## Females left

##### Onder Lab
WA_lmers_alt$WA_L_F_Onder_lmer_alt <- lmer(CentroidSizeLeft_micrometers ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "Onder" & Sex == "F"))

##### Posnien Lab
WA_lmers_alt$WA_L_F_Posnien_lmer_alt <- lmer(CentroidSizeLeft_micrometers ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "Posnien" & Sex == "F"))

##### Ritchie Lab
WA_lmers_alt$WA_L_F_Ritchie_lmer_alt <- lmer(CentroidSizeLeft_micrometers ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "Ritchie" & Sex == "F"))

##### Onder Lab
WA_lmers_alt$WA_L_F_StamenkovicRadak_lmer_alt <- lmer(CentroidSizeLeft_micrometers ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "StamenkovicRadak" & Sex == "F"))


## Males left

##### Onder Lab
WA_lmers_alt$WA_L_M_Onder_lmer_alt <- lmer(CentroidSizeLeft_micrometers ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "Onder" & Sex == "M"))

##### Posnien Lab, singular fit, removed Population
WA_lmers_alt$WA_L_M_Posnien_lmer_alt <- lmer(CentroidSizeLeft_micrometers ~ Altitude + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "Posnien" & Sex == "M"))

##### Ritchie Lab
WA_lmers_alt$WA_L_M_Ritchie_lmer_alt <- lmer(CentroidSizeLeft_micrometers ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "Ritchie" & Sex == "M"))

##### Onder Lab
WA_lmers_alt$WA_L_M_StamenkovicRadak_lmer_alt <- lmer(CentroidSizeLeft_micrometers ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "StamenkovicRadak" & Sex == "M"))



## Females right

##### Onder Lab
WA_lmers_alt$WA_R_F_Onder_lmer_alt <- lmer(CentroidSizeRight_micrometers ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "Onder" & Sex == "F"))

##### Posnien Lab
WA_lmers_alt$WA_R_F_Posnien_lmer_alt <- lmer(CentroidSizeRight_micrometers ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "Posnien" & Sex == "F"))

##### Ritchie Lab
WA_lmers_alt$WA_R_F_Ritchie_lmer_alt <- lmer(CentroidSizeRight_micrometers ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "Ritchie" & Sex == "F"))

##### Onder Lab
WA_lmers_alt$WA_R_F_StamenkovicRadak_lmer_alt <- lmer(CentroidSizeRight_micrometers ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "StamenkovicRadak" & Sex == "F"))


## Mmales right

##### Onder Lab
WA_lmers_alt$WA_R_M_Onder_lmer_alt <- lmer(CentroidSizeRight_micrometers ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "Onder" & Sex == "M"))

##### Posnien Lab, singular fit, removed Population
WA_lmers_alt$WA_R_M_Posnien_lmer_alt <- lmer(CentroidSizeRight_micrometers ~ Altitude + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "Posnien" & Sex == "M"))

##### Ritchie Lab
WA_lmers_alt$WA_R_M_Ritchie_lmer_alt <- lmer(CentroidSizeRight_micrometers ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "Ritchie" & Sex == "M"))

##### Onder Lab
WA_lmers_alt$WA_R_M_StamenkovicRadak_lmer_alt <- lmer(CentroidSizeRight_micrometers ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$wa, Supervisor.PI == "StamenkovicRadak" & Sex == "M"))

# save output list
saveRDS(WA_lmers_alt, file = file.path(lmer_dir, out_dir, "WA_lmers_alt.rds"))






############# FECUNDITY ############# 

# create output directory
out_dir <- "Fecundity"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 

# initialize output list
Fec_lmers_alt <- list()

##### Billeter Lab
Fec_lmers_alt$Fec_Billeter_lmer_alt <- lmer(NumberOfAdultsEclosed ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$fec, Supervisor.PI == "Billeter"))

##### Fricke Lab, singular fit, removed Population
Fec_lmers_alt$Fec_Fricke_lmer_alt <- lmer(NumberOfAdultsEclosed ~ Altitude + (1|Line:Population), data = filter(droseu$fec, Supervisor.PI == "Fricke"))

# save output list
saveRDS(Fec_lmers_alt, file = file.path(lmer_dir, out_dir, "Fec_lmers_alt.rds"))







############# LIFESPAN ############# 

# create output directory
out_dir <- "Lifespan"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 

# initialize output list
LS_lmers_alt <- list()


## Females

#### Flatt Lab
LS_lmers_alt$LS_F_Flatt_lmer_alt <- lmer(LSP_AgeAtDeath_days ~ Altitude + (1|Population), data = filter(droseu$lsp, Censor == "0" & Supervisor.PI == "Flatt" & Sex == "F"))

#### Parsch Lab
LS_lmers_alt$LS_F_Parsch_lmer_alt <- lmer(LSL_AgeAtDeath_days ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$lsl, Censor == "0" & Supervisor.PI == "Parsch" & Sex == "F"))

#### Pasyukova Lab
LS_lmers_alt$LS_F_Pasyukova_lmer_alt <- lmer(LSL_AgeAtDeath_days ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$lsl, Censor == "0" & Supervisor.PI == "Pasyukova" & Sex == "F"))


## Males

#### Flatt Lab
LS_lmers_alt$LS_M_Flatt_lmer_alt <- lmer(LSP_AgeAtDeath_days ~ Altitude + (1|Population), data = filter(droseu$lsp, Censor == "0" & Supervisor.PI == "Flatt" & Sex == "M"))

#### Parsch Lab
LS_lmers_alt$LS_M_Parsch_lmer_alt <- lmer(LSL_AgeAtDeath_days ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$lsl, Censor == "0" & Supervisor.PI == "Parsch" & Sex == "M"))

#### Pasyukova Lab
LS_lmers_alt$LS_M_Pasyukova_lmer_alt <- lmer(LSL_AgeAtDeath_days ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$lsl, Censor == "0" & Supervisor.PI == "Pasyukova" & Sex == "M"))


# save output list
saveRDS(LS_lmers_alt, file = file.path(lmer_dir, out_dir, "LS_lmers_alt.rds"))







############# COLD-SHOCK MORTALITY ############# 

# create output directory
out_dir <- "ColdShock"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 

# initialize output list
CSM_lmers_alt <- list()


## Females

#### Gonzalez Lab
CSM_lmers_alt$CSM_F_Gonzalez_lmer_alt <- lmer(CSM_PropDead_ED_asin ~ Altitude + (1|Line:Population) + (1|Population), data = filter(droseu$csm, Supervisor.PI == "Gonzalez" & Sex == "F"))

#### Kozeretska Lab, singular fit, removed population
CSM_lmers_alt$CSM_F_Kozeretska_lmer_alt <- lmer(CSM_PropDead_ED_asin ~ Altitude + (1|Line:Population), data = filter(droseu$csm, Supervisor.PI == "Kozeretska" & Sex == "F"))

#### Vieira Lab
CSM_lmers_alt$CSM_F_Vieira_lmer_alt <- lmer(CSM_PropDead_ED_asin ~ Altitude + (1|Line:Population) + (1|Population), data = filter(droseu$csm, Supervisor.PI == "Vieira" & Sex == "F"))


## Males

#### Gonzalez Lab
CSM_lmers_alt$CSM_M_Gonzalez_lmer_alt <- lmer(CSM_PropDead_ED_asin ~ Altitude + (1|Line:Population) + (1|Population), data = filter(droseu$csm, Supervisor.PI == "Gonzalez" & Sex == "M"))

#### Kozeretska Lab, singular fit, removed population
CSM_lmers_alt$CSM_M_Kozeretska_lmer_alt <- lmer(CSM_PropDead_ED_asin ~ Altitude + (1|Line:Population), data = filter(droseu$csm, Supervisor.PI == "Kozeretska" & Sex == "M"))

#### Vieira Lab
CSM_lmers_alt$CSM_M_Vieira_lmer_alt <- lmer(CSM_PropDead_ED_asin ~ Altitude + (1|Line:Population) + (1|Population), data = filter(droseu$csm, Supervisor.PI == "Vieira" & Sex == "M"))

# save output list
saveRDS(CSM_lmers_alt, file = file.path(lmer_dir, out_dir, "CSM_lmers_alt.rds"))






############# CHILL-COMA RECOVERY TIME ############# 

# create output directory
out_dir <- "ChillComa"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 

# initialize output list
CCRT_lmers_alt <- list()


## Females

#### Vieira Lab
CCRT_lmers_alt$CCRT_F_Vieira_lmer_alt <- lmer(CCRT_seconds ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$ccrt, Censor == "0" & Supervisor.PI == "Vieira" & Sex == "F"))

#### Mensch Lab
CCRT_lmers_alt$CCRT_F_Mensch_lmer_alt <- lmer(CCRT_seconds ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$ccrt, Censor == "0" & Supervisor.PI == "Mensch" & Sex == "F"))


## Males

#### Vieira Lab
CCRT_lmers_alt$CCRT_M_Vieira_lmer_alt <- lmer(CCRT_seconds ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$ccrt, Censor == "0" & Supervisor.PI == "Vieira" & Sex == "M"))

#### Mensch Lab
CCRT_lmers_alt$CCRT_M_Mensch_lmer_alt <- lmer(CCRT_seconds ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$ccrt, Censor == "0" & Supervisor.PI == "Mensch" & Sex == "M"))

# save output list
saveRDS(CCRT_lmers_alt, file = file.path(lmer_dir, out_dir, "CCRT_lmers_alt.rds"))






############# HEAT-SHOCK MORTALITY ############# 

# create output directory
out_dir <- "HeatShock"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 

# initialize output list
HSM_lmers_alt <- list()

## Females

#### Parsch Lab, not converging, removed Population
HSM_lmers_alt$HSM_F_Parsch_lmer_alt <- lmer(TimeDeath_min ~ Altitude + (1|Line:Population), data = filter(droseu$hsm, Censor == "0" & Supervisor.PI == "Parsch" & Sex == "F"))

#### Vieira Lab
HSM_lmers_alt$HSM_F_Vieira_lmer_alt <- lmer(TimeDeath_min ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$hsm, Censor == "0" & Supervisor.PI == "Vieira" & Sex == "F"))


## Males

#### Parsch Lab
HSM_lmers_alt$HSM_M_Parsch_lmer_alt <- lmer(TimeDeath_min ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$hsm, Censor == "0" & Supervisor.PI == "Parsch" & Sex == "M"))

#### Vieira Lab
HSM_lmers_alt$HSM_M_Vieira_lmer_alt <- lmer(TimeDeath_min ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$hsm, Censor == "0" & Supervisor.PI == "Vieira" & Sex == "M"))

# save output list
saveRDS(HSM_lmers_alt, file = file.path(lmer_dir, out_dir, "HSM_lmers_alt.rds"))






############# DIAPAUSE ############# 

# create output directory
out_dir <- "Diapause"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 

# initialize output list
Dia_lmers_alt <- list()

#### Bergland Lab
Dia_lmers_alt$Dia_Bergland_lmer_alt <- lmer(Prop_Max_Stage9_asin ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$dia, Supervisor.PI == "Bergland"))

#### Flatt Lab, removed Population because of singular fit
Dia_lmers_alt$Dia_Flatt_lm_alt <- lm(Prop_Max_Stage9_asin ~ Altitude, data = filter(droseu$dia, Supervisor.PI == "Flatt"))

#### Schlotterer Lab
Dia_lmers_alt$Dia_Schlotterer_lmer_alt <- lmer(Prop_Max_Stage9_asin ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$dia, Supervisor.PI == "Schlotterer"))

# save output list
saveRDS(Dia_lmers_alt, file = file.path(lmer_dir, out_dir, "Dia_lmers_alt.rds"))




## GLMERs

dia <- read.csv("Data/MasterSheets_May22_git/DIA_MasterSheet_Feb22.csv")

dia <- group_by(dia, Supervisor.PI, Diet, Batch, Population, Line, Altitude) %>%
  summarise(Prop_Max_Stage9 = mean(MostAdvancedStage <= 9 & NumberOfEggs == 0), 
            Prop_Max_Stage9_asin = asin(sqrt(Prop_Max_Stage9)),          
            n = n())

Dia_glmers_alt <- foreach(pi = unique(dia$Supervisor.PI)) %do% {
  m1 <- glmer(Prop_Max_Stage9 ~ Altitude + (1|Population) + (1|Line:Population), weights = n, family = binomial(), data = filter(dia, Supervisor.PI == pi)) 
  m2 <- glmer(Prop_Max_Stage9 ~ + (1|Population) + (1|Line:Population), weights = n, family = binomial(), data = filter(dia, Supervisor.PI == pi))
  a <- anova(m1, m2)
  l <- list(m1, a)
  names(l) <- c(paste0("Dia_", pi, "_glmer_alt"), paste0("Dia_", pi, "_glmer_alt_anova")) 
  l }

Dia_glmers_alt <- unlist(Dia_glmers_alt, recursive=FALSE) 
Dia_glmers_alt_anova <- Dia_glmers_alt[grep("anova", names(Dia_glmers_alt))]
names(Dia_glmers_alt_anova) <- sub("_anova", "", names(Dia_glmers_alt_anova))
Dia_glmers_alt <- Dia_glmers_alt[grep("anova", names(Dia_glmers_alt), invert = T)]

saveRDS(Dia_glmers_alt, file = file.path(lmer_dir, out_dir, "Dia_glmers_alt.rds"))
saveRDS(Dia_glmers_alt_anova, file = file.path(lmer_dir, out_dir, "Dia_glmers_alt_anova.rds"))
capture.output(Dia_glmers_alt_anova, file = file.path(lmer_dir, out_dir, "Dia_glmers_alt_anova.txt"))








############# CIRCADIAN ECLOSION TIMING ############# 

# create output directory
out_dir <- "CircadianEclosion"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 

# initialize output list
CET_lmers_alt <- list()




############# LOCOMOTOR ACTIVITY ############# 

# create output directory
out_dir <- "Locomotor"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 

# initialize output list
LA_lmers_alt <- list()

#### Tauber Lab, removed "-Inf"
LA_lmers_alt$LA_NDlog2_Tauber_lmer_alt <- lmer(ND_log2 ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$la, ND_log2 != -Inf)) # log2 transformation creates -Inf values (2 cases), should they be removed?

#### Tauber Lab, untransformed ND
#LA_lmers_alt$LA_ND_Tauber_lmer_alt <- lmer(ND ~ Altitude + (1|Population) + (1|Line:Population), data = droseu$la)


LA_lmers_alt$LA_Period_Tauber_lmer_alt <- lmer(Period ~ Altitude + (1|Population) + (1|Line:Population), data = droseu$la)

LA_lmers_alt$LA_CircPhase_Tauber_lmer_alt <- lmer(CircPhase ~ Altitude + (1|Population) + (1|Line:Population), data = droseu$la)

# singular fit, removed Line
LA_lmers_alt$LA_AbsPhase_Tauber_lmer_alt <- lmer(AbsPhase ~ Altitude + (1|Population), data = droseu$la)

LA_lmers_alt$LA_Activity_Tauber_lmer_alt <- lmer(Activity ~ Altitude + (1|Population) + (1|Line:Population), data = droseu$la)

# save output list
saveRDS(LA_lmers_alt, file = file.path(lmer_dir, out_dir, "LA_lmers_alt.rds"))







############# STARVATION RESISTANCE ############# 

# create output directory
out_dir <- "Starvation"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 

# initialize output list
SR_lmers_alt <- list()

## Females

#### Gonzalez Lab
SR_lmers_alt$SR_F_Gonzalez_lmer_alt <- lmer(AgeAtDeath_hours ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$sr, Supervisor.PI == "Gonzalez" & Sex == "F"))

#### Onder Lab
SR_lmers_alt$SR_F_Onder_lmer_alt <- lmer(AgeAtDeath_hours ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$sr, Supervisor.PI == "Onder" & Sex == "F"))

#### Pasyukova Lab
SR_lmers_alt$SR_F_Pasyukova_lmer_alt <- lmer(AgeAtDeath_hours ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$sr, Supervisor.PI == "Pasyukova" & Sex == "F"))

## Males

#### Gonzalez Lab
SR_lmers_alt$SR_M_Gonzalez_lmer_alt <- lmer(AgeAtDeath_hours ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$sr, Supervisor.PI == "Gonzalez" & Sex == "M"))

#### Onder Lab
SR_lmers_alt$SR_M_Onder_lmer_alt <- lmer(AgeAtDeath_hours ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$sr, Supervisor.PI == "Onder" & Sex == "M"))

#### Pasyukova Lab
SR_lmers_alt$SR_M_Pasyukova_lmer_alt <- lmer(AgeAtDeath_hours ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$sr, Supervisor.PI == "Pasyukova" & Sex == "M"))

# save output list
saveRDS(SR_lmers_alt, file = file.path(lmer_dir, out_dir, "SR_lmers_alt.rds"))






############# PIGMENTATION ############# 

# create output directory
out_dir <- "Pigmentation"
dir.create(file.path(lmer_dir, out_dir), showWarnings = F) 

# initialize output list
Pgm_lmers_alt <- list()

#### Abbott Lab
Pgm_lmers_alt$Pgm_T4_Abbott_lmer_alt <- lmer(PercT4_asin ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$pgm, Supervisor.PI == "Abbott"))

Pgm_lmers_alt$Pgm_T5_Abbott_lmer_alt <- lmer(PercT5_asin ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$pgm, Supervisor.PI == "Abbott"))

Pgm_lmers_alt$Pgm_T6_Abbott_lmer_alt <- lmer(PercT6_asin ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$pgm, Supervisor.PI == "Abbott"))

Pgm_lmers_alt$Pgm_Total_Abbott_lmer_alt <- lmer(TotalPerc_asin ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$pgm, Supervisor.PI == "Abbott"))

#### Gibert Lab
Pgm_lmers_alt$Pgm_T4_Gibert_lmer_alt <- lmer(PercT4_asin ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$pgm, Supervisor.PI == "Gibert"))

Pgm_lmers_alt$Pgm_T5_Gibert_lmer_alt <- lmer(PercT5_asin ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$pgm, Supervisor.PI == "Gibert"))

Pgm_lmers_alt$Pgm_T6_Gibert_lmer_alt <- lmer(PercT6_asin ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$pgm, Supervisor.PI == "Gibert"))

# singular fit, removed population
Pgm_lmers_alt$Pgm_Total_Gibert_lmer_alt <- lmer(TotalPerc_asin ~ Altitude + (1|Line:Population), data = filter(droseu$pgm, Supervisor.PI == "Gibert"))

#### Schmidt Lab
Pgm_lmers_alt$Pgm_T4_Schmidt_lmer_alt <- lmer(PercT4_asin ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$pgm, Supervisor.PI == "Schmidt"))

Pgm_lmers_alt$Pgm_T5_Schmidt_lmer_alt <- lmer(PercT5_asin ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$pgm, Supervisor.PI == "Schmidt"))

Pgm_lmers_alt$Pgm_T6_Schmidt_lmer_alt <- lmer(PercT6_asin ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$pgm, Supervisor.PI == "Schmidt"))

Pgm_lmers_alt$Pgm_Total_Schmidt_lmer_alt <- lmer(TotalPerc_asin ~ Altitude + (1|Population) + (1|Line:Population), data = filter(droseu$pgm, Supervisor.PI == "Schmidt"))


# save output list
saveRDS(Pgm_lmers_alt, file = file.path(lmer_dir, out_dir, "Pgm_lmers_alt.rds"))










###############################################################################
###############################################################################
###############################################################################

######### combine all linear models into global lists, one for all LMERs and one for all GLMERs

### LMERs

# list with 2 levels where each element is a trait and sub element is a lab
all_lmers_alt <- rdsBatchReaderToList(path = lmer_dir, recursive = T, full.names = T, pattern = "_lmers_alt.rds")
# flatten the list and rename elements
all_lmers_alt <- unlist(all_lmers_alt, recursive=FALSE)
names(all_lmers_alt) <- str_split(names(all_lmers_alt), "\\.", simplify = T)[,2]
# output the list
saveRDS(all_lmers_alt, file = file.path(lmer_dir, "all_lmers_alt_list.rds"))

### GLMERs - currently for diapause only
all_glmers_alt <- rdsBatchReaderToList(path = lmer_dir, recursive = T, full.names = T, pattern = "_glmers_alt.rds")
all_glmers_alt <- unlist(all_glmers_alt, recursive=FALSE)
names(all_glmers_alt) <- str_split(names(all_glmers_alt), "\\.", simplify = T)[,2]
saveRDS(all_glmers_alt, file = file.path(lmer_dir, "all_glmers_alt_list.rds"))



######### output all model summaries, anovas and tukeys as global lists
######### anovas only for LMERs as they are already performed in the trait sections for GLMERs 
######### tukeys only when applicable (eg for Population)


### LMERs

all_lmers_alt <- readRDS(file.path(lmer_dir, "all_lmers_alt_list.rds"))
# summary
all_lmers_alt_summary <- lapply(all_lmers_alt, summary)
saveRDS(all_lmers_alt_summary, file = file.path(lmer_dir, "all_lmers_alt_summary_list.rds"))
# anova
all_lmers_alt_anova <- lapply(all_lmers_alt, anova)
saveRDS(all_lmers_alt_anova, file.path(lmer_dir, "all_lmers_alt_anova_list.rds"))


### GLMERs

all_glmers_alt <- readRDS(file.path(lmer_dir, "all_glmers_alt_list.rds"))
# summary
all_glmers_alt_summary <- lapply(all_glmers_alt, summary)
saveRDS(all_glmers_alt_summary, file = file.path(lmer_dir, "all_glmers_alt_summary_list.rds"))
# anova
all_glmers_alt_anova <- rdsBatchReaderToList(path = lmer_dir, recursive = T, full.names = T, pattern = "_glmers_alt_anova.rds")
all_glmers_alt_anova <- unlist(all_glmers_alt_anova, recursive=FALSE)
names(all_glmers_alt_anova) <- str_split(names(all_glmers_alt_anova), "\\.", simplify = T)[,2]
saveRDS(all_glmers_alt_anova, file = file.path(lmer_dir, "all_glmers_alt_anova_list.rds"))





######### plot linear models residuals per trait

# list models outputs to keep the directory structure
lmers <- list.files(path = lmer_dir, recursive = T, full.names = T, pattern = "_lmers_alt.rds")
glmers <- list.files(path = lmer_dir, recursive = T, full.names = T, pattern = "_glmers_alt.rds")
models <- c(lmers, glmers)

# loop over models to get residuals
for (m in 1:length(models)){
  plotResiduals(models[m])
}


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

all_lmers_alt_anova <- readRDS(file.path(lmer_dir, "all_lmers_alt_anova_list.rds"))
all_glmers_alt_anova <- readRDS(file.path(lmer_dir, "all_glmers_alt_anova_list.rds"))
all_alt_anova <- c(all_lmers_alt_anova, all_glmers_alt_anova)

pvalues <- combinePValues(all_alt_anova)

saveRDS(pvalues, file = file.path(lmer_dir, "all_models_alt_pvalues.rds"))
write.csv(pvalues, file = file.path(lmer_dir, "all_models_alt_pvalues.csv"), row.names = F)






