# Figure 5 - Cytokines 
# Author: Kathrin Jutzeler
# Date: May 23, 2023
# Last updated: March 7, 2024

# Load packages ####
library(tidyverse)
library(readxl)   # Import Excel data
library(rstatix)  # Pipe frindly stats 
library(wesanderson)  # Color palettes
library(ggpubr) # ggplot add-on 
library(ggprism) # add_pvalue to plots
library(emmeans) # Mean adjustment for linear models 
library(multcomp) # Perform mutliple comparisons 
library(multcompView) # Add lettes to stats
library(plotrix) # for standard error
library(patchwork)

# Define theme and colors ####
my_theme <- function(){theme_minimal() + theme(panel.grid = element_blank(), legend.position = "none", strip.placement = 'outside',
                             axis.title.x = element_blank(), axis.line = element_line(1),
                             panel.spacing = unit(0,'lines'), 
                             text = element_text(size = 14)) }

my_color = c("Control" = wes_palette("Darjeeling2")[1],
             "LE" = wes_palette("Darjeeling2")[2],
             "EG" = wes_palette("Darjeeling2")[3],
             "OR" = "#710193",
             "BRE" = wes_palette("Darjeeling2")[5])

# Import the data ####

all_df <- read_csv('all_data.csv')

all_df$pop <- factor(all_df$pop, levels = c('Control', 'BRE', 'EG', 'LE', 'OR'),
                          labels = c('Ctrl', 'BRE', 'EG', 'LE', 'OR'))

cytokine_col <- c("IFNy", "IL5",  "TNFa", "IL2",  "IL6",  "IL4",  "IL10", "IL13")

# Normalize by penetration rate
all_df <- all_df %>%
  mutate_at(vars(cytokine_col), ~ ifelse(pop != 'Ctrl', ./penrate, .))

# Perform normality tests 
norms <- all_df %>%
  group_by(host, pop) %>%
  summarize_at(vars(cytokine_col,-IL13), list(normality = ~shapiro.test(.)$p) )

shapiro.test(all_df$IL13)

# Non-parametric groups -> use Kruskal and Wilcox for all cytokines

#+++++++++++++++++
# Functions ####
#+++++++++++++++++

# To plot 

f_plot <- function(cytokine) {
  
  placeholder <- deparse(substitute(cytokine))
  
  # Get the maximum value for the vertical line
  yend <- round(max(all_df[placeholder], na.rm = T))
  
  # Get  y-axis label
  ylab <- paste0(placeholder, ' (pg/ml)')

  ggplot(all_df, aes(pop, {{cytokine}})) +
    geom_boxplot(aes(fill = pop)) +
    geom_vline(data=all_df[all_df$host=="BALB/c",], aes(xintercept=5.5, linetype="solid")) +
    geom_vline(aes(xintercept=1.5), linetype="dashed") +
    facet_wrap(~host, strip.position = 'bottom') + 
    scale_fill_manual(values = my_color) +
    ylab(ylab) +
    my_theme() 
}

# To convert p table to letters 
f_letters <- function(p_table) {
  
  empty <- NULL
  # Subset the data to separate hosts
  for (i in c('BALB/c', 'C57BL/6')) {
  
  output <- subset(p_table, host == i)
  
  # Convert the tibble to a data frame with row names 
  res <- data.frame(Group1=output$group1, Group2=output$group2, PValue=output$p.adj) %>% 
    pivot_wider(values_from="PValue", names_from=c("Group1"))
  
  df <- as.data.frame(res)
  rownames(df) <- df$Group2
  df$Group2 <- NULL
  
  # Make a table suitable for multcompletters
  ft <- fullPTable(df)
  # Get letters
  res2 <- multcompLetters(ft)
  df <- as.data.frame(res2$Letters)
  colnames(df) <- 'letter'
  df$pop <- rownames(df)
  df$host <- i
  
  empty <- bind_rows(empty, df)
  }
  return(empty)
}

#+++++++++++++++++
# Plot IFN-y ####
#+++++++++++++++++
# Perform stats
stat_IFN <- all_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(host) %>%
  kruskal_test(IFNy ~ pop) %>%
  mutate(pop = 'EG', method = 'K-W') # Arbitrary value for plotting

IFN_host <- all_df %>%
  ungroup() %>%
  filter(pop != 'Ctrl') %>%
  group_by(pop) %>%
  wilcox_test(IFNy ~ host) %>%
  adjust_pvalue(method = 'BH') %>%
  add_significance(symbols = c('####', '###', '##', '#', 'ns') ) %>%
  mutate(host = "C57BL/6")


# Plot
plot_IFN <- f_plot(IFNy) +
  ylab("IFN\u03B3 (pg/ml)") +
  geom_text(data = stat_IFN, aes(label = paste0(method, ', p = ',format(p, nsmall =3))), y = 47000) +
  geom_text(data = IFN_host, aes(label = p.adj.signif), vjust = -0.5, y = 7000)

#+++++++++++++++++
# Plot IL-2 ####
#+++++++++++++++++
# Perform stats
stat_IL2 <- all_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(host) %>%
  kruskal_test(IL2 ~ pop) %>%
  mutate(pop = 'EG', method = 'K-W')

IL2_dunn <- all_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(host) %>%
  dunn_test(IL2 ~ pop) %>%
  mutate(pop = rep(c('BRE', 'EG', 'LE', 'OR'),3))

IL2_host <- all_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(pop) %>%
  wilcox_test(IL2 ~ host) %>%
  adjust_pvalue(method = 'BH') %>%
  add_significance(symbols = c('####', '###', '##', '#', 'ns') ) %>%
  filter(p.adj <= 0.05) %>%
  mutate(host = "C57BL/6")

# Plot
plot_IL2 <- f_plot(IL2) +
  geom_text(data = stat_IL2, aes(label = paste0(method, ', p = ',round(p, 3))), y = 370) +
  geom_text(data = subset(IL2_dunn, host == 'BALB/c'), aes(label = p.adj.signif), vjust = -0.5, y = 300)

#+++++++++++++++++
# Plot IL-4 ####
#+++++++++++++++++
# Perform stats
stat_IL4 <- all_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(host) %>% 
  kruskal_test(IL4 ~ pop) %>%
  mutate(pop = 'EG', method = 'K-W')

IL4_host <- all_df %>%
  group_by(pop) %>% 
  #ungroup() %>%
  filter(pop != 'Ctrl') %>%
  wilcox_test(IL4 ~ host) %>%
  adjust_pvalue(method = 'BH') %>%
  add_significance(symbols = c('####', '###', '##', '#', 'ns') ) %>%
  filter(p.adj <= 0.05) %>%
  mutate(host = "C57BL/6")

# Plot
plot_IL4 <- f_plot(IL4) +
  geom_text(data = stat_IL4, aes(label = paste0(method, ', p = ',format(p,nsmall =3))), y = 20000) +
  geom_text(data = IL4_host, aes(label = p.adj.signif), vjust = -0.5, y = 7000)

#+++++++++++++++
# Plot IL-5 ####
#+++++++++++++++

# Perform stats
stat_IL5 <- all_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(host) %>%
  kruskal_test(IL5 ~ pop) %>%
  mutate(pop = 'EG', method = 'K-W')

IL5_host <- all_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(pop) %>%
  wilcox_test(IL5 ~ host) %>%
  adjust_pvalue(method = 'BH') %>%
  add_significance(symbols = c('####', '###', '##', '#', 'ns') ) %>%
  filter(p.adj <= 0.05) %>%
  mutate(host = "C57BL/6") %>%
  add_xy_position()

# Plot
plot_IL5 <- f_plot(IL5) +
  ylim(0,3000) +
  geom_text(data = stat_IL5, aes(label = paste0(method, ', p = ',round(p,3))), y = 3000) +
  geom_text(data = subset(IL5_host, pop != 'LE'), aes(label = p.adj.signif), vjust = -0.5, y = 1400) +
  geom_text(data = subset(IL5_host, pop == 'LE'), aes(label = p.adj), vjust = -0.5, y = 1400) 

#++++++++++++++
# Plot IL6 ####
#++++++++++++++
# Perform stats
stat_IL6 <- all_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(host) %>%
  kruskal_test(IL6 ~ pop) %>%
  mutate(pop = 'EG', method = 'K-W')

IL6_host <- all_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(pop) %>%
  wilcox_test(IL6 ~ host) %>%
  adjust_pvalue(method = 'BH') %>%
  add_significance(symbols = c('####', '###', '##', '#', 'ns') ) %>%
  filter(p.adj <= 0.05) %>%
  mutate(host = "C57BL/6") %>%
  add_y_position()

# Plot
plot_IL6 <- f_plot(IL6) +
  geom_text(data = stat_IL6, aes(label = paste0(method, ', p = ',format(p, 1))), y = 120000) +
  geom_text(data = subset(IL6_host, pop %in% c('EG', 'OR')), aes(label = p.adj.signif), vjust = -0.5, y = 20000) +
  geom_text(data = subset(IL6_host, !pop %in% c('EG', 'OR')), aes(label = p.adj.signif), vjust = -0.5, y = 20000)

#+++++++++++++++++
# Plot IL10 ####
#+++++++++++++++++

# Perform stats
stat_IL10 <- all_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(host) %>%
  kruskal_test(IL10 ~ pop) %>%
  mutate(pop = 'EG', method = 'K-W')

IL10_host <- all_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(pop) %>%
  wilcox_test(IL10 ~ host) %>%
  adjust_pvalue(method = 'BH') %>%
  add_significance(symbols = c('####', '###', '##', '#', 'ns') ) %>%
  filter(p <= 0.05) %>%
  mutate(host = "C57BL/6") %>%
  add_y_position()

plot_IL10 <- f_plot(IL10) +
  geom_text(data = stat_IL10, aes(label = paste0(method, ', p = ',p)), y = 270)  +
  geom_text(data = IL10_host, aes(label = p.adj.signif), vjust = -0.5, y = 100)

#+++++++++++++++
# Plot IL13 ####
#+++++++++++++++

# Perform stats
stat_IL13 <- all_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(host) %>%
  kruskal_test(IL13 ~ pop) %>%
  mutate(pop = 'EG', method = 'K-W')

IL13_host <- all_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(pop) %>%
  wilcox_test(IL13 ~ host) %>%
  adjust_pvalue(method = 'BH') %>%
  add_significance(symbols = c('####', '###', '##', '#', 'ns') ) %>%
  filter(p <= 0.05) %>%
  mutate(host = "C57BL/6") %>%
  add_y_position()

# Plot 
plot_IL13 <- f_plot(IL13) +
  geom_text(data = stat_IL13, aes(label = paste0(method, ', p = ',p)), y = 105) +
  geom_text(data = IL13_host, aes(label = p.adj.signif), vjust = -0.5, y = 60)

#+++++++++++++++++
# Plot TNFa ####
#+++++++++++++++++

# Perform stats
stat_TNFa <- all_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(host) %>%
  kruskal_test(TNFa ~ pop) %>%
  mutate(pop = 'EG', method = 'K-W')

TNFa_host <- all_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(pop) %>%
  wilcox_test(TNFa ~ host) %>%
  adjust_pvalue(method = 'BH') %>%
  add_significance(symbols = c('####', '###', '##', '#', 'ns') ) %>%
  filter(p <= 0.05) %>%
  mutate(host = "C57BL/6") %>%
  add_y_position()

# Plot
plot_TNFa <- f_plot(TNFa) +
  ylab("TNF\u03B1 (pg/ml)") +
  geom_text(data = stat_TNFa, aes(label = paste0(method, ', p = ',format(p,1))), y = 20500) +
  geom_text(data = TNFa_host, aes(label = p.adj.signif), vjust = -0.5, y = 7000) 

#++++++++++++++++++++++
# EXPORT PLOTS ####
#++++++++++++++++++++++

jpeg("Figure5.jpg", width = 10, height = 12, unit = 'in', res =300)

(plot_IFN + plot_TNFa) / (plot_IL4 | plot_IL5) / (plot_IL6 | plot_spacer())  +
 plot_annotation(tag_levels = c("A")) &  theme(plot.tag = element_text(face = 'bold'), axis.title.x = element_blank())

dev.off()


jpeg("Additional file FigureS4.jpg", width = 10, height = 8, unit = 'in', res =300)

(plot_IL2 + plot_IL10) / (plot_IL13 | plot_spacer()) +
  plot_annotation(tag_levels = c("A")) &  theme(plot.tag = element_text(face = 'bold'), axis.title.x = element_blank())

dev.off()
