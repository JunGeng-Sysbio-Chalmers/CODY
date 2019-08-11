library(ggplot2)

library(grid)
library(gridExtra)
library(ggrepel)
library(ggsci)
library(cowplot)
library(ggthemes)

library(dplyr)
library(tidyr)
library(MASS)

library(extrafont)

library(ggpmisc)

####################
# loading functions for plotting
####################
source("speciesPlot_12M_Jun.R")
source("metsPlot_12M_Jun.R")

####################
# loading data
####################

load("./Rdata_12M/12M_lumen_T1.Rdata")
load("./Rdata_12M/12M_lumen_T2.Rdata")
load("./Rdata_12M/12M_lumen_T3.Rdata")
load("./Rdata_12M/12M_lumen_T4.Rdata")

load("./Rdata_12M/12M_mucosa_T1.Rdata")
load("./Rdata_12M/12M_mucosa_T2.Rdata")
load("./Rdata_12M/12M_mucosa_T3.Rdata")
load("./Rdata_12M/12M_mucosa_T4.Rdata")

load("./Rdata_12M/12M_blood.Rdata")
load("./Rdata_12M/12M_feces.Rdata")
load("./Rdata_12M/12M_rectum.Rdata")

#############
# preprocessing data
#############

# only select 10% of samples from the dataset, and only keep 17 columns
sample_index <- seq(1, nrow(lumen_T1), by=10)

##############################
# lumen
#
##########
# lumen Tank 1 
lumen_T1_s <- lumen_T1[sample_index, 1:19]
lumen_T1_bac <- lumen_T1_s[,1:9]
colnames(lumen_T1_bac) <- c("t", "Btheta", "Bfragilis", "Blongum", "Bbreve", "Bado", "Ehalli", "Fpransnitzii", "Rintestinalls")
lumen_T1_mets <- lumen_T1_s[, c(1, 10:19)]
colnames(lumen_T1_mets) <- c("t", "macs", "hexose", "succinate", "acetate", "propionate", "lactate", "formate", "ethanol", "butyrate", "h2")

speciesPlot_12M_Jun(lumen_T1_bac, "./fig_12M_Jun/12M_lumen_T1_species", "Ascending Lumen")
metsPlot_12M_Jun(lumen_T1_mets, "./fig_12M_Jun/12M_lumen_T1_mets", "Ascending Lumen")

##########
# lumen Tank 2
lumen_T2_s <- lumen_T2[sample_index, 1:19]
lumen_T2_bac <- lumen_T2_s[,1:9]
colnames(lumen_T2_bac) <- c("t", "Btheta", "Bfragilis", "Blongum", "Bbreve", "Bado", "Ehalli", "Fpransnitzii", "Rintestinalls")
lumen_T2_mets <- lumen_T2_s[, c(1, 10:19)]
colnames(lumen_T2_mets) <- c("t", "macs", "hexose", "succinate", "acetate", "propionate", "lactate", "formate", "ethanol", "butyrate", "h2")

speciesPlot_12M_Jun(lumen_T2_bac, "./fig_12M_Jun/12M_lumen_T2_species", "Transverse Lumen")
metsPlot_12M_Jun(lumen_T2_mets, "./fig_12M_Jun/12M_lumen_T2_mets", "Transverse Lumen")

##########
# lumen Tank 3
lumen_T3_s <- lumen_T3[sample_index, 1:19]
lumen_T3_bac <- lumen_T3_s[,1:9]
colnames(lumen_T3_bac) <- c("t", "Btheta", "Bfragilis", "Blongum", "Bbreve", "Bado", "Ehalli", "Fpransnitzii", "Rintestinalls")
lumen_T3_mets <- lumen_T3_s[, c(1, 10:19)]
colnames(lumen_T3_mets) <- c("t", "macs", "hexose", "succinate", "acetate", "propionate", "lactate", "formate", "ethanol", "butyrate", "h2")

speciesPlot_12M_Jun(lumen_T3_bac, "./fig_12M_Jun/12M_lumen_T3_species", "Descending Lumen")
metsPlot_12M_Jun(lumen_T3_mets, "./fig_12M_Jun/12M_lumen_T3_mets", "Descending Lumen")

##########
# lumen Tank 4
lumen_T4_s <- lumen_T4[sample_index, 1:19]
lumen_T4_bac <- lumen_T4_s[,1:9]
colnames(lumen_T4_bac) <- c("t", "Btheta", "Bfragilis", "Blongum", "Bbreve", "Bado", "Ehalli", "Fpransnitzii", "Rintestinalls")
lumen_T4_mets <- lumen_T4_s[, c(1, 10:19)]
colnames(lumen_T4_mets) <- c("t", "macs", "hexose", "succinate", "acetate", "propionate", "lactate", "formate", "ethanol", "butyrate", "h2")

speciesPlot_12M_Jun(lumen_T4_bac, "./fig_12M_Jun/12M_lumen_T4_species", "Sigmoid Lumen")
metsPlot_12M_Jun(lumen_T4_mets, "./fig_12M_Jun/12M_lumen_T4_mets", "Sigmoid Lumen")


##############################
# mucosa
#
##########
# mucosa Tank 1 
mucosa_T1_s <- mucosa_T1[sample_index, 1:19]
mucosa_T1_bac <- mucosa_T1_s[,1:9]
colnames(mucosa_T1_bac) <- c("t", "Btheta", "Bfragilis", "Blongum", "Bbreve", "Bado", "Ehalli", "Fpransnitzii", "Rintestinalls")
mucosa_T1_mets <- mucosa_T1_s[, c(1, 10:19)]
colnames(mucosa_T1_mets) <- c("t", "macs", "hexose", "succinate", "acetate", "propionate", "lactate", "formate", "ethanol", "butyrate", "h2")

speciesPlot_12M_Jun(mucosa_T1_bac, "./fig_12M_Jun/12M_mucosa_T1_species", "Ascending Mucus")
metsPlot_12M_Jun(mucosa_T1_mets, "./fig_12M_Jun/12M_mucosa_T1_mets", "Ascending Mucus")

##########
# mucosa Tank 2
mucosa_T2_s <- mucosa_T2[sample_index, 1:19]
mucosa_T2_bac <- mucosa_T2_s[,1:9]
colnames(mucosa_T2_bac) <- c("t", "Btheta", "Bfragilis", "Blongum", "Bbreve", "Bado", "Ehalli", "Fpransnitzii", "Rintestinalls")
mucosa_T2_mets <- mucosa_T2_s[, c(1, 10:19)]
colnames(mucosa_T2_mets) <- c("t", "macs", "hexose", "succinate", "acetate", "propionate", "lactate", "formate", "ethanol", "butyrate", "h2")

speciesPlot_12M_Jun(mucosa_T2_bac, "./fig_12M_Jun/12M_mucosa_T2_species", "Transverse Mucus")
metsPlot_12M_Jun(mucosa_T2_mets, "./fig_12M_Jun/12M_mucosa_T2_mets", "Transverse Mucus")

##########
# mucosa Tank 3
mucosa_T3_s <- mucosa_T3[sample_index, 1:19]
mucosa_T3_bac <- mucosa_T3_s[,1:9]
colnames(mucosa_T3_bac) <- c("t", "Btheta", "Bfragilis", "Blongum", "Bbreve", "Bado", "Ehalli", "Fpransnitzii", "Rintestinalls")
mucosa_T3_mets <- mucosa_T3_s[, c(1, 10:19)]
colnames(mucosa_T3_mets) <- c("t", "macs", "hexose", "succinate", "acetate", "propionate", "lactate", "formate", "ethanol", "butyrate", "h2")

speciesPlot_12M_Jun(mucosa_T3_bac, "./fig_12M_Jun/12M_mucosa_T3_species", "Descending Mucus")
metsPlot_12M_Jun(mucosa_T3_mets, "./fig_12M_Jun/12M_mucosa_T3_mets", "Descending Mucus")

##########
# mucosa Tank 4
mucosa_T4_s <- mucosa_T4[sample_index, 1:19]
mucosa_T4_bac <- mucosa_T4_s[,1:9]
colnames(mucosa_T4_bac) <- c("t", "Btheta", "Bfragilis", "Blongum", "Bbreve", "Bado", "Ehalli", "Fpransnitzii", "Rintestinalls")
mucosa_T4_mets <- mucosa_T4_s[, c(1, 10:19)]
colnames(mucosa_T4_mets) <- c("t", "macs", "hexose", "succinate", "acetate", "propionate", "lactate", "formate", "ethanol", "butyrate", "h2")

speciesPlot_12M_Jun(mucosa_T4_bac, "./fig_12M_Jun/12M_mucosa_T4_species", "Sigmoid Mucus")
metsPlot_12M_Jun(mucosa_T4_mets, "./fig_12M_Jun/12M_mucosa_T4_mets", "Sigmoid Mucus")

##############################
# blood
#
blood_s <- blood[sample_index, 1:11]
blood_mets <- blood_s
colnames(blood_mets) <- c("t", "macs", "hexose", "succinate", "acetate", "propionate", "lactate", "formate", "ethanol", "butyrate", "h2")

metsPlot_12M_Jun(blood_mets, "./fig_12M_Jun/12M_blood_mets", "Blood")

##############################
# rectum
#
rectum_s <- rectum[sample_index, 1:19]
rectum_bac <- rectum_s[,1:9]
colnames(rectum_bac) <- c("t", "Btheta", "Bfragilis", "Blongum", "Bbreve", "Bado", "Ehalli", "Fpransnitzii", "Rintestinalls")
rectum_mets <- rectum_s[, c(1, 10:19)]
colnames(rectum_mets) <- c("t", "macs", "hexose", "succinate", "acetate", "propionate", "lactate", "formate", "ethanol", "butyrate", "h2")

speciesPlot_12M_Jun(rectum_bac, "./fig_12M_Jun/12M_rectum_species", "Rectum")
metsPlot_12M_Jun(rectum_mets, "./fig_12M_Jun/12M_rectum_mets", "Rectum")


##############################
# feces
#
feces_s <- feces[sample_index, 1:19]
feces_bac <- feces_s[,1:9]
colnames(feces_bac) <- c("t", "Btheta", "Bfragilis", "Blongum", "Bbreve", "Bado", "Ehalli", "Fpransnitzii", "Rintestinalls")
feces_mets <- feces_s[, c(1, 10:19)]
colnames(feces_mets) <- c("t", "macs", "hexose", "succinate", "acetate", "propionate", "lactate", "formate", "ethanol", "butyrate", "h2")

speciesPlot_12M_Jun(feces_bac, "./fig_12M_Jun/12M_feces_species", "Feces")
metsPlot_12M_Jun(feces_mets, "./fig_12M_Jun/12M_feces_mets", "Feces")

