#rm(list=ls())
library(RColorBrewer)
library(immunotation) #freq allele
library(genepop)
library(tidyverse)
library(janitor)
library(ggplot2)
library(ggpubr) #qqplot
library(patchwork)
library(cowplot)
library(gridExtra) #tablegrob
library(grid)
library(gtsummary)
library(flextable)
library(rstatix) #pvale
library(stringr)
library(HIBAG)
library(correlation) #correlation analyses
library(corrplot)
library(xlsx)
library(emmeans) #lm
library(data.table) #lm
library(car) #glm
library(MASS) #glm.nb
library(ggplotify)
library(DHARMa)
library(rmdformats)

#sessioninfo::package_info()

#libraries <- c("RColorBrewer")
#lapply(x, require, character.only = TRUE)
# COLORS -----------------------------------------------------------------------

#https://coolors.co/palette/fdba31-370075-620097-e2dada-000000
# https://coolors.co/gradient-palette/fdba31-95722e?number=7
gradient_grey <- RColorBrewer::brewer.pal(n = 8, name = "Greys")
gradient_orange <- RColorBrewer::brewer.pal(n = 8, name = "Oranges")
gradient_purple <- RColorBrewer::brewer.pal(n = 8, name = "Purples")
gradiente_diverg <- RColorBrewer::brewer.pal(n = 8, name = "PuOr")


# Control - Allele Frequencies Brazil REDOME Parana (Bone Marrow Donors
control_freq_redome <-
  immunotation::query_allele_frequencies(hla_locus = "B", hla_population = 3181)


# DATASET
data0 <- read.csv("/home/schagas/hlafozpublic/dataset_final.csv", sep = ",")

data_final <- data0 %>%
  filter(HospitalPeriod_days < 65)
data0 %>%
  filter(SamplesGroupYear == 2020) %>%
  nrow()

data0 %>%
  filter(SamplesGroupYear == 2021) %>%
  nrow()
