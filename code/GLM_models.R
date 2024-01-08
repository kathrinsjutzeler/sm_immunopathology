# GLM models for all paramters
# Author: Kathrin Jutzeler
# Date: May 16, 2023
# Last updated: September 7, 2023
# R version 4.2.0, tidyverse version 1.3.2, ggplot2 3.3.6      ✔ purrr   0.3.4 
#✔ tibble  3.1.8      ✔ dplyr   1.0.10
#✔ tidyr   1.2.0      ✔ stringr 1.4.1 
#✔ readr   2.1.2      ✔ forcats 0.5.2 

# Load packages
library(tidyverse)
library(readxl)   # Import Excel data
library(gtsummary) # To export pretty table
library(ggpubr)


# Define labels
labels <-  c(pop ~ 'Parasite', host ~ 'Host', total_eggs ~ 'Total eggs')
CBC_labels <- c(pop ~ 'Parasite', host ~ 'Host', total_eggs ~ 'Total eggs')

# Import the data ####

all_df <- read_csv("all_data.csv")

all_df$host <- as.factor(all_df$host)
all_df$pop <- factor(all_df$pop, levels = c('Control', 'BRE', 'EG', 'LE', 'OR'), 
                     labels = c('Control', 'BRE', 'EG', 'LE', 'OR'))

all_df <- all_df %>%
  filter(pop != 'Control') %>%
  droplevels()

CBC <- read_csv("CBC_data.csv")
CBC$time <- factor(CBC$time, levels = c('Baseline', 'Week 2', 'Week 4', 'Week 6', 'Week 8', 'Week 10'))

# Manipulate CBC data
cells <- c('EO#(K/uL)', 'NEUT#(K/uL)', 'MONO#(K/uL)', 'LYMPH#(K/uL)', 'RET#(K/uL)', 'BASO#(K/uL)') 

cell_df <- list()

for (i in cells) {
  cell_df[[i]] <- CBC %>%
    dplyr::select(sample, host, pop, time, i) %>% 
    spread(time, i) 
}

CBC_df <- cell_df$`EO#(K/uL)` %>%
  left_join(cell_df$`NEUT#(K/uL)`, by = 'sample', suffix = c('', '_NEUT')) %>%
  left_join(cell_df$`MONO#(K/uL)`, by = 'sample', suffix = c('', '_MONO')) %>%
  left_join(cell_df$`LYMPH#(K/uL)`, by = 'sample', suffix = c('', '_LYMPH')) %>%
  left_join(cell_df$`BASO#(K/uL)`, by = 'sample', suffix = c('', '_BASO')) %>%
  left_join(cell_df$`RET#(K/uL)`, by = 'sample', suffix = c('_EO', '_RET')) 

CBC_df <- CBC_df %>% rename(c(host = host_EO, pop = pop_EO)) %>%
  dplyr::select(-host_NEUT, -host_LYMPH, -host_MONO, -pop_NEUT, -pop_LYMPH, -pop_MONO,
                -host_RET, -pop_RET,-host_BASO, -pop_BASO)


CBC_df$host <- as.factor(CBC_df$host)
CBC_df$pop <- factor(CBC_df$pop, levels = c('Control', 'BRE', 'EG', 'LE', 'OR'), 
                     labels = c('Control', 'BRE', 'EG', 'LE', 'OR'))

CBC_df <- CBC_df %>%
  filter(pop != 'Control') %>%
  left_join(dplyr::select(all_df,sample, total_eggs), by = 'sample') %>%
  droplevels()

# Test relationship between intestine eggs and liver eggs
ggplot(all_df, aes(intestine_eggs, liver_eggs)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  stat_cor(method = "spearman", label.x = 5000, label.y = 12) 


# Set up functions ============================================================

# Generate all GLMs
f_glm <- function(trait){
  f <- paste0(trait, "~ pop * host + total_eggs")
  
  trait_glm <- glm(f,  data = subset(all_df))
}

# Use log transformation for all cytokines
f_glm_cyto  <- function(trait){
  f <- paste0('log(', trait, ") ~ pop * host + total_eggs")
  
  trait_glm <- glm(f,  data = subset(all_df))
}

# Function to generate all cytokine tables
f_tbl_cytokine <- function(input){
  tbl_regression(input, label = labels) %>% 
    bold_p(t = 0.05) %>%
    bold_labels() %>%
    italicize_levels() %>%
    add_global_p(include = c(pop, host, `pop:host`))
}

# List of parameters

traits <- c("liver_wt_norm", "spleen_wt_norm", "intestine_length", "fibrosis", "granuloma")
traits_list <- as.list(traits)
names(traits_list) <- traits

# Run all GLMs and then assess fit individually
glms <- lapply(traits_list, f_glm)

# Do the same with cytokines

cytokines <- c("IFNy", "TNFa", "IL2", "IL4", "IL5", "IL6", "IL10", "IL13")
cytokine_list <- as.list(cytokines)
names(cytokine_list) <- cytokines

cytokine_glms <- lapply(cytokine_list, f_glm_cyto)


# LIVER WEIGHT ============================================================

# Models
#liver_glm <- glm(liver_wt_norm ~ pop * host + total_eggs, data = all_df)
liver_glm1 <- glm(sqrt(liver_wt_norm) ~ pop * host + total_eggs, data = all_df, family = gaussian)
#liver_glm2 <- glm(log(liver_wt_norm) ~ pop * host * total_eggs, data = all_df, family = gaussian)

summary(liver_glm1)

plot(liver_glm1, 1) # evaluate plots 1:3

# Use transformed model, better fit and not much of a different outcome
tbl_liver <- 
tbl_regression(liver_glm1, label = labels) %>% 
  modify_caption("**Table 2. Generalized linear model output for all traits**") %>%
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  italicize_levels() %>%
  add_global_p()

#glm(sqrt(liver_wt_norm) ~ pop * host * total_eggs, data = all_df, family = gaussian) %>%
#  tbl_regression() %>%
#  add_global_p(type = "II")

  
#SPLEEN WEIGHT ============================================================

# Models
#spleen_glm <- glm(log(spleen_wt_norm) ~ pop * host +  total_eggs , data = all_df)
spleen_glm2 <- glm(sqrt(spleen_wt_norm) ~ pop * host + total_eggs , data = all_df)

summary(spleen_glm2)

plot(spleen_glm,3)
plot(spleen_glm2,3)

# Sqrt transformed model, better fit 

tbl_spleen <- tbl_regression(spleen_glm2, label = labels) %>%
  modify_caption("**Table 2. Generalized linear model output for all traits**") %>%
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  italicize_levels() %>%
  add_global_p(include = c(pop, host, `pop:host`))

#em <- emmeans(spleen_glm, "pop")
#contrast(em, "pairwise", adjust = "BH")


# INTESTINE LENGTH ============================================================

# Models
intestine_glm <- glm(intestine_length ~ pop * host + total_eggs,  data = all_df)

summary(intestine_glm)

# Use regular model - outcome doesn't change with transformed data
tbl_intestine <- tbl_regression(intestine_glm, label = labels) %>% 
  modify_caption("**Table 2. Generalized linear model output for all traits**") %>%
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  italicize_levels() %>%
  add_global_p(include = c(pop, host, `pop:host`))


#FIBROTIC AREA ============================================================

# Models
fibrosis_glm <- glm(fibrosis ~ pop * host + total_eggs,  data = all_df)

summary(fibrosis_glm)

# Use regular model - outcome doesn't change with transformed data
tbl_fibrosis <- tbl_regression(glms$fibrosis, label = labels) %>% 
  #modify_caption("**Table 2. Generalized linear model output for spleen weight**") %>%
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  italicize_levels() %>%
  add_global_p()


# GRANULOMA AREA =================================================

# Models

granuloma_glm <- glm(granuloma ~ pop * host + total_eggs,  data = all_df)

summary(granuloma_glm)

# Use regular model - outcome doesn't change with transformed data
tbl_granuloma <- tbl_regression(glms$granuloma, label = labels) %>% 
  #modify_caption("**Table 2. Generalized linear model output for spleen weight**") %>%
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  italicize_levels() %>%
  add_global_p()


# IFN-y ============================================================

## Models

#IFN_glm <- glm(log(IFNy) ~ pop * host + total_eggs,  data = all_df)

#summary(IFN_glm)

# Use predefined model (with log transformation) for all cytokines
tbl_IFNy <- tbl_regression(cytokine_glms$IFNy, labels = labels) %>% 
#  modify_caption("**Table 2. Generalized linear model output for spleen weight**") %>%
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  italicize_levels() %>%
  add_global_p(include = c(pop, host, `pop:host`))


# TNF-a ============================================================

## Models 

#TNF_glm <- glm(log(TNFa) ~ pop * host + total_eggs,  data = all_df)

#summary(TNF_glm)

tbl_TNFa <- tbl_regression(cytokine_glms$TNFa, label = labels) %>% 
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  italicize_levels() %>%
  add_global_p(include = c(pop, host, `pop:host`))

# Eosinophils ===========================================================
#eo_glm <- glm(`Week 10_EO` ~ pop * host + total_eggs + Baseline_EO,  data = CBC_df)
eo_glm2 <- glm(sqrt(`Week 10_EO`) ~ pop * host + total_eggs + Baseline_EO,  data = CBC_df)

summary(eo_glm2)

# Chose transformed model due to better fit

tbl_EO <- tbl_regression(eo_glm2, label = c(CBC_labels, Baseline_EO ~ 'Baseline')) %>% 
  modify_caption("**Table 3. Generalized linear model output for CBC**") %>%
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  italicize_levels() %>%
  add_global_p(include = c(pop, host, `pop:host`))

#Neutrophils ===========================================================
neut_glm <- glm(sqrt(`Week 10_NEUT`) ~ pop * host + total_eggs + Baseline_NEUT,  data = CBC_df)

summary(neut_glm)

plot(neut_glm, 1)

# Chose transformed model due to better fit

tbl_NEUT <- tbl_regression(neut_glm, label = c(CBC_labels, Baseline_NEUT ~ 'Baseline')) %>% 
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  italicize_levels() %>%
  add_global_p(include = c(pop, host, `pop:host`))


# Monocytes ===========================================================
mono_glm <- glm(sqrt(`Week 10_MONO`) ~ pop * host + total_eggs + Baseline_MONO,  data = CBC_df)

summary(mono_glm)


# Chose transformed model due to better fit

tbl_MONO <- tbl_regression(mono_glm, label = c(CBC_labels, Baseline_MONO ~ 'Baseline')) %>% 
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  italicize_levels() %>%
  add_global_p(include = c(pop, host, `pop:host`))


# Lymphocytes ===========================================================

lymph_glm <- glm(log(`Week 10_LYMPH`) ~ pop * host + total_eggs + Baseline_LYMPH,  data = CBC_df)

summary(lymph_glm)

# Chose transformed model due to better fit

tbl_LYMPH <- tbl_regression(lymph_glm, label = c(CBC_labels, Baseline_LYMPH ~ 'Baseline')) %>% 
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  italicize_levels() %>%
  add_global_p(include = c(pop, host, `pop:host`))

#Basophils ===========================================================
baso_glm <- glm(sqrt(`Week 10_BASO`) ~ pop * host + total_eggs + Baseline_BASO,  data = CBC_df)
BIC(baso_glm)
summary(baso_glm)

# Chose transformed model due to better fit

tbl_baso <- tbl_regression(baso_glm, label = c(CBC_labels, Baseline_BASO ~ 'Baseline')) %>% 
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  italicize_levels() %>%
  add_global_p(include = c(pop, host, `pop:host`))


# Reticulocytes ===========================================================
ret_glm <- glm(log(`Week 10_RET`) ~ pop * host + total_eggs + Baseline_RET,  data = CBC_df)

summary(ret_glm)



# Chose transformed model due to better fit

tbl_ret <- tbl_regression(ret_glm, label = c(CBC_labels, Baseline_RET ~ 'Baseline')) %>% 
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  italicize_levels() %>%
  add_global_p(include = c(pop, host, `pop:host`))



tbl_regression(model)

#Cytokine tables ====================================================

tbl_cytokines <- lapply(cytokine_glms, f_tbl_cytokine)


# Merge tables ============================================================
## Table 2 ####
tables <- list(tbl_liver, tbl_spleen, tbl_intestine, tbl_fibrosis, tbl_granuloma)

tbl_merge(tbls = tables, 
          tab_spanner = c("**Liver weight**", "**Spleen weight**", '**Intestine length**',
                          '**Fibrotic area**', '**Granuloma area**'  ) )
                      
## Table 3 ####
CBC_tables <- list(tbl_EO, tbl_NEUT, tbl_MONO, tbl_LYMPH, tbl_baso, tbl_ret)

tbl_merge(tbls = CBC_tables, 
          tab_spanner = c("**Eosinophils**", "**Neutrophils**", '**Monocytes**',
                          '**Lymphocytes**', '**Basophils**', '**Reticulocytes**'))

## Table 4 ####
tbl_merge(tbls = tbl_cytokines, 
          tab_spanner = c('**IFNy**', '**TNFa**', '**IL-2**', '**IL4**', '**IL5**', '**IL6**',
                          '**IL10**', '**IL13**') )




# To export and keep formatting, save as web page and then select everything and paste into Word or Excel
