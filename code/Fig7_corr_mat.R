# Figure 7 - Correleation matrix 
# Author: Kathrin Jutzeler
# Date: May 26, 2023
# Last updated: October 17, 2023

# Load packages
library(tidyverse) # for everything
library(patchwork) # to export plots
library(Hmisc) # to add p value to correlation matrix 
library(corrplot) # for correlation plot
library(ggplotify) # to convert plot to ggplot object

# Define theme, colors, labels etc. ####
my_theme <- function(){theme(panel.grid = element_blank(), legend.position = "none", strip.placement = 'outside',
                             axis.line = element_line(1),
                             panel.spacing = unit(0,'lines'), 
                             text = element_text(size = 14),
                             axis.text = element_text(size=12, color="black"),
                             legend.text = element_text(size=12, color="black"),
                             legend.title = element_text(size=13, color="black"),
                             plot.title = element_text(size=14, face='bold')) }

# Labels
labs <- c("Intestine eggs", "Liver eggs", 'Total worms', 'Liver weight', 'Spleen weight', 'Intestine length', 
          'Fecundity', 'Fibrotic area', 'Granuloma area', 
          'IFNy', 'TNFa', 'IL2', 'IL4', 'IL5', 'IL6', 'IL10', 'IL13')

# Import data #### 
all_df <- read.csv('all_data.csv') 

all_df$host <- as.factor(all_df$host)
all_df$pop <- as.factor(all_df$pop)

all_df <- all_df %>% filter(pop != 'Control') %>% droplevels()

#+++++++++++++++
# Formulas ####
#+++++++++++++++

# This one selects the data we want to show and normalizes non-parameteric data
f_cormatrix <- function(mouse){
  em_full <- all_df %>%
    filter(host == mouse ) %>%
      dplyr::select(intestine_eggs, liver_eggs, total_worms, liver_wt_norm, spleen_wt_norm, intestine_length,
                  fecundity, fibrosis, granuloma,
                  IFNy, TNFa, IL2, IL4, IL5, IL6, IL10, IL13)

  # Perform normality test
  test <- em_full %>%
      summarize_all(list(shapiro = ~shapiro.test(.)$p.value) )

  test <- gather(test, parameter, value)

  transform <- filter(test, value <= 0.05)
  transform <- str_remove(transform$parameter, "_shapiro")

  # Log transform non-parametric data
  em_full <- em_full %>%
    mutate_at(transform, list( ~ log(.)))

  em_full <- em_full %>% mutate_if(is.numeric, function(x) ifelse(is.infinite(x), 0, x))

  em_full <- as.matrix(em_full)
  return(em_full)
}

# This one generates the matrix and plot

f_corplot <- function(matrix){
  
  matrix <- lapply(matrix, function(x) `rownames<-` (x, labs))
  matrix <-lapply(matrix, function(x) `colnames<-` (x, labs))
  
  matrix$P[is.na(matrix$P)] <- 0
  
  corrplot(matrix$r, p.mat = matrix$P, sig.level = 0.05, #order = 'hclust', 
           type = 'upper', insig = 'blank', tl.srt =45, tl.col = 'black',
           tl.cex = 0.75)
}

#++++++++++++++++++++++++++
# Correlation matrix  ####
#++++++++++++++++++++++++++

cormat_black6 <- f_cormatrix('C57BL/6')
cormat_balbc <- f_cormatrix('BALB/c')

res_black6 <- rcorr(cormat_black6, type = "pearson")
res_balbc <- rcorr(cormat_balbc, type = "pearson")

#++++++++++++++++++++++
#Export the plots #### 
#++++++++++++++++++++++

png('black_6.png', width = 5, height =4, units = 'in', res = 150)
f_corplot(res_black6)

dev.off()

png('balb_c.png', width = 5, height =4, units = 'in', res = 150)

f_corplot(res_balbc)

dev.off()

