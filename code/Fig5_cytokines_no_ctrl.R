# Figure 5 - Cytokines 
# Author: Kathrin Jutzeler
# Date: May 23, 2023
# Last updated

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
    geom_segment(data=all_df[all_df$host=="BALB/c",], aes(x=5.5, y=0, xend=5.5, yend=yend), linetype="solid",
                 color = "black", linewidth=0.5) +
    geom_segment(data=all_df[all_df$host=="BALB/c",], aes(x=1.5, y=0, xend=1.5, yend=yend), linetype="dashed",
                 color = "black", linewidth=0.5) +
    geom_segment(data=all_df[all_df$host=="C57BL/6",], aes(x=1.5, y=0, xend=1.5, yend=yend), linetype="dashed",
                 color = "black", linewidth=0.5) +
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
  geom_text(data = stat_IFN, aes(label = paste0(method, ', p = ',p)), y = 42000) +
  geom_text(data = IFN_host, aes(label = p.adj.signif), vjust = -0.5, y = 7000)

#+++++++++++++++++
# Plot IL-2 ####
#+++++++++++++++++
# Perform stats
stat_IL2 <- cytokine_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(host) %>%
  kruskal_test(IL2 ~ pop) %>%
  mutate(pop = 'EG', method = 'K-W')

IL2_host <- cytokine_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(pop) %>%
  wilcox_test(IL2 ~ host) %>%
  adjust_pvalue(method = 'BH') %>%
  add_significance(symbols = c('####', '###', '##', '#', 'ns') ) %>%
  mutate(host = "C57BL/6")

# Plot
plot_IL2 <- f_plot(IL2) +
  geom_text(data = stat_IL2, aes(label = paste0(method, ', p = ',p)), y = 380) +
  geom_text(data = IL2_host, aes(label = p.adj.signif), vjust = -0.5, y = 300)

#+++++++++++++++++
# Plot IL-4 ####
#+++++++++++++++++
# Perform stats
stat_IL4 <- cytokine_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(host) %>% 
  kruskal_test(IL4 ~ pop) %>%
  mutate(pop = 'EG', method = 'K-W')

IL4_host <- cytokine_df %>%
  group_by(pop) %>% 
  #ungroup() %>%
  filter(pop != 'Ctrl') %>%
  wilcox_test(IL4 ~ host) %>%
  adjust_pvalue(method = 'BH') %>%
  add_significance(symbols = c('####', '###', '##', '#', 'ns') ) %>%
  mutate(host = "C57BL/6")

# Plot
plot_IL4 <- f_plot(IL4) +
  geom_text(data = stat_IL4, aes(label = paste0(method, ', p = ',p)), y = 17000) +
  geom_text(data = IL4_host, aes(label = p.adj.signif), vjust = -0.5, y = 7000)


# Plot IL-5 ===================================================================

# Perform stats
stat_IL5 <- cytokine_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(host) %>%
  kruskal_test(IL5 ~ pop) %>%
  mutate(pop = 'EG', method = 'K-W')

IL5_host <- cytokine_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(pop) %>%
  wilcox_test(IL5 ~ host) %>%
  adjust_pvalue(method = 'BH') %>%
  add_significance(symbols = c('####', '###', '##', '#', 'ns') ) %>%
  mutate(host = "C57BL/6") %>%
  add_xy_position()

# Plot
plot_IL5 <- f_plot(IL5) +
  geom_text(data = stat_IL5, aes(label = paste0(method, ', p = ',p)), y = 2100) +
  geom_text(data = IL5_host, aes(label = p.adj.signif), vjust = -0.5, y = 1400) 

#+++++++++++++++++
# Plot IL6 =================================================================
#+++++++++++++++++
# Perform stats
stat_IL6 <- cytokine_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(host) %>%
  kruskal_test(IL6 ~ pop) %>%
  mutate(pop = 'EG', method = 'K-W')

IL6_host <- cytokine_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(pop) %>%
  wilcox_test(IL6 ~ host) %>%
  adjust_pvalue(method = 'BH') %>%
  add_significance(symbols = c('####', '###', '##', '#', 'ns') ) %>%
  mutate(host = "C57BL/6") %>%
  add_y_position()

# Plot
plot_IL6 <- f_plot(IL6) +
  geom_text(data = stat_IL6, aes(label = paste0(method, ', p = ',p)), y = 120000) +
  geom_text(data = IL6_host, aes(label = p.adj.signif), vjust = -0.5, y = 20000) 

#+++++++++++++++++
# Plot IL10 ####
#+++++++++++++++++

# Perform stats
stat_IL10 <- cytokine_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(host) %>%
  kruskal_test(IL10 ~ pop) %>%
  mutate(pop = 'EG', method = 'K-W')

IL10_host <- cytokine_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(pop) %>%
  wilcox_test(IL10 ~ host) %>%
  adjust_pvalue(method = 'BH') %>%
  add_significance(symbols = c('####', '###', '##', '#', 'ns') ) %>%
  mutate(host = "C57BL/6") %>%
  add_y_position()

plot_IL10 <- f_plot(IL10) +
  geom_text(data = stat_IL10, aes(label = paste0(method, ', p = ',p)), y = 280)  +
  geom_text(data = IL10_host, aes(label = p.adj.signif), vjust = -0.5, y = 100)


# Plot IL13 ==================================================================

# Perform stats
stat_IL13 <- cytokine_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(host) %>%
  kruskal_test(IL13 ~ pop) %>%
  mutate(pop = 'EG', method = 'K-W')

IL13_host <- cytokine_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(pop) %>%
  wilcox_test(IL13 ~ host) %>%
  adjust_pvalue(method = 'BH') %>%
  add_significance(symbols = c('####', '###', '##', '#', 'ns') ) %>%
  mutate(host = "C57BL/6") %>%
  add_y_position()

# Plot 
plot_IL13 <- f_plot(IL13) +
  geom_text(data = stat_IL13, aes(label = paste0(method, ', p = ',p)), y = 115) +
  geom_text(data = IL13_host, aes(label = p.adj.signif), vjust = -0.5, y = 60)

#+++++++++++++++++
# Plot TNFa ####
#+++++++++++++++++

# Perform stats
stat_TNFa <- cytokine_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(host) %>%
  kruskal_test(TNFa ~ pop) %>%
  mutate(pop = 'EG', method = 'K-W')

TNFa_host <- cytokine_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(pop) %>%
  wilcox_test(TNFa ~ host) %>%
  adjust_pvalue(method = 'BH') %>%
  add_significance(symbols = c('####', '###', '##', '#', 'ns') ) %>%
  mutate(host = "C57BL/6") %>%
  add_y_position()

# Plot
plot_TNFa <- f_plot(TNFa) +
  geom_text(data = stat_TNFa, aes(label = paste0(method, ', p = ',p)), y = 19000) +
  geom_text(data = TNFa_host, aes(label = p.adj.signif), vjust = -0.5, y = 7000) 


#++++++++++++++++++++++
# EXPORT PLOTS ####
#++++++++++++++++++++++

pdf("Figure5_ctrl_excluded.pdf", width = 10, height = 4)

(plot_IFN + plot_TNFa)  +
 plot_annotation(tag_levels = c("A")) &  theme(plot.tag = element_text(face = 'bold'), axis.title.x = element_blank())

dev.off()


pdf("FigureS4.pdf", width = 10, height = 12)

(plot_IL2 + plot_IL4) / (plot_IL5 | plot_IL6) / (plot_IL10 | plot_IL13) +
  plot_annotation(tag_levels = c("A")) &  theme(plot.tag = element_text(face = 'bold'), axis.title.x = element_blank())

dev.off()
