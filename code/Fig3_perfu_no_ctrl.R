# Figure 3 - Tissue and histology data
# Author: Kathrin Jutzeler
# Date: May 16, 2023
# Last updated: October 19, 2023
# R version 4.2.0, tidyverse version 1.3.2, ggplot2 3.3.6      ✔ purrr   0.3.4 
#✔ tibble  3.1.8      ✔ dplyr   1.0.10
#✔ tidyr   1.2.0      ✔ stringr 1.4.1 
#✔ readr   2.1.2      ✔ forcats 0.5.2 

# Load packages ####
library(tidyverse)
#library(readxl)   # Import Excel data
library(rstatix)  # Pipe frindly stats 
library(wesanderson)  # Color palettes
library(ggpubr) # ggplot add-on 
library(ggprism) # add_pvalue to plots
#library(magick) # Import images and convert to ggplot
library(emmeans) # Mean adjustment for linear models 
library(rcompanion) # For letter display
library(multcomp) # Perform mutliple comparisons 
library(multcompView)
library(ggrepel) # To add labels at the end line plots
library(plotrix) # for standard error
library(PMCMRplus) # Post-hoc test for Friedman
library(patchwork)


# Define theme and colors ####
my_theme <- function(){theme_minimal() +theme(panel.grid = element_blank(), legend.position = "none", strip.placement = 'outside',
        axis.title.x = element_blank(), axis.line = element_line(1),
        panel.spacing = unit(0,'lines'), 
        text = element_text(size = 14)) }

my_color = c("Control" = wes_palette("Darjeeling2")[1],
             "Ctrl"= wes_palette("Darjeeling2")[1],
             "LE" = wes_palette("Darjeeling2")[2],
             "EG" = wes_palette("Darjeeling2")[3],
             "OR" = "#710193",
             "BRE" = wes_palette("Darjeeling2")[5])

#++++++++++++++++++++++++
# Import the data ####
#++++++++++++++++++++++++
# Weight data
weight <- read_csv('body_weight_data.csv')
weight$host <- as.factor(weight$host)
weight$pop <- factor(weight$pop,  levels = c('Control', 'BRE', 'EG', 'LE', 'OR'))
weight$time <- as.factor(weight$time)

# Images
#liverimage <- magick::image_read('livers.png', TRUE)
#spleenimage <- magick::image_read('spleens.png', TRUE)

# Granuloma data 
df_granulomas <- read_csv("granuloma_data.csv")

  # Change column name 
  colnames(df_granulomas)[4] = 'area'

df_granulomas$pop <- factor(df_granulomas$pop, levels = c("BRE", "EG", "LE", "OR"))
                            
# Remaining data
all_df <- read_csv('all_data.csv')

all_df$pop <- factor(all_df$pop, levels = c('Control', 'BRE', 'EG', 'LE', 'OR'), 
                    labels = c('Ctrl', 'BRE', 'EG', 'LE', 'OR'))

#+++++++++++++++++++++++++++
# Function for letters ####
#+++++++++++++++++++++++++++
f_letters <- function(output, mouse) {
  # Subset the data to separate hosts
  output <- subset(output, host == mouse)
  
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
  df$host <- mouse
  return(df)
}


#++++++++++++++++++++++
# BODY WEIGHT  ####
#++++++++++++++++++++++
# Perform normality tests
weight %>%
  filter(time != 0) %>%
  group_by(host, pop,time) %>%
  shapiro_test(gain) %>%
  filter(p <= 0.05)

# Calculate mean and se per group
weight_summary <- weight %>%
  na.omit() %>%
  group_by(host, pop, time) %>%
  summarize(mean_weight = mean(weight), se = std.error(weight), 
            mean_gain = mean(gain), se_gain = std.error(gain))

##### This is the final plot #####
# Repeated measures 
fgain <- weight_summary %>%
  filter(pop != 'Control') %>%
  group_by(host) %>%
  friedman_test(mean_gain ~ pop | time) %>%
  p_round(p, 3)

#weight_summary %>%
#  group_by(pop) %>%
#  friedman_test(mean_gain ~ host | time) %>%
#  p_round(p, 3)

# Artificial values for plotting
fgain$pop <- c(1, 2)
fgain$time <- c(2,2)

#black <- filter(weight_summary, host == "C57BL/6", pop != 'Control') %>% droplevels()
#res.b <- PMCMRplus::frdAllPairsConoverTest(y = black$mean_gain, groups = black$pop, 
                           #      blocks = black$time, p.adjust.method = "BH")

white <- filter(weight_summary, host == "BALB/c", pop != 'Control') %>% droplevels()
res.w <- PMCMRplus::frdAllPairsConoverTest(y = white$mean_gain, groups = white$pop, 
                                         blocks = white$time, p.adjust.method = "BH")

resw <- rcompanion::fullPTable(res.w$p.value)
res <- multcompLetters(resw)
w <- as.data.frame(res$Letters)

#resb <- fullPTable(res.b$p.value)
#res <- multcompLetters(resb)
#b <- as.data.frame(res$Letters)

df <-  data.frame(letter= c(' ', w$`res$Letters`)) #,' ', b$`res$Letters`))

stat <- weight_summary %>%
  filter(host == 'BALB/c', time == 12) %>%
  bind_cols( df)

# Weight gain
plot_weight <- 
  ggplot(weight_summary, aes(time, mean_gain, group = interaction(pop))) +
  geom_point(size =1.5, position = position_dodge(width = 0.5), aes(color = pop)) +
  geom_line(data = subset(weight_summary, pop != 'Control'), position = position_dodge(width = 0.5), aes(color = pop)) +
  geom_line(data=subset(weight_summary,pop == 'Control'), position = position_dodge(width = 0.5), aes(color = pop),
            linetype = 'dashed') +
  geom_errorbar(aes(ymin=mean_gain-se_gain, ymax=mean_gain+se_gain, 
                    color = pop), width=.2,
                position=position_dodge(0.5)) +
  facet_wrap(host~.,  strip.position = 'bottom') +
  ylab("Mean weight gain (%)") +
  xlab("Weeks post infection") +
  labs(color = "Population") + 
  scale_color_manual(values = my_color) +
  my_theme() +
  geom_text(data = fgain, aes(label = paste0('Friedman p = ', format(round(p, 3),3)), x = 4, y =25)) +
  geom_text_repel(
    aes(label = letter, color = pop), data = stat, position = position_dodge(width = 0.5), hjust = 0.5,
    show.legend = F) +
  geom_vline(data=weight_summary[weight_summary$host=="BALB/c",], aes(xintercept=13.6), 
             linetype = 'solid')+
  theme(legend.position = 'top', axis.title.x = element_text()) 

#++++++++++++++++++++++
# LIVER WEIGHT ####
#++++++++++++++++++++++
all_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(host, pop) %>%
  shapiro_test(liver_wt_norm)

bartlett.test(data = subset(all_df, host == "BALB/c" & pop != 'Ctrl'), liver_wt_norm ~ pop)
#p-value = 0.1839
# Anova for BALB/c and Kruskal for C57BL/6

# Perform stats tests
all_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(host) %>%
  anova_test(liver_wt_norm ~ pop)

all_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(host) %>%
  kruskal_test(liver_wt_norm ~ pop)

liver_res <- data.frame(host = c('BALB/c', 'C57BL/6'), test = c('ANOVA', 'K-W'), p = c(0.003, 0.021),
                        pop = 'EG')

p_liver <- all_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(host) %>%
  tukey_hsd(liver_wt_norm ~ pop, p.adjust.method = 'BH') %>%
  add_xy_position(x = 'pop')

p_liver2 <- all_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(host) %>%
  dunn_test(liver_wt_norm ~ pop, p.adjust.method = 'BH') %>%
  add_xy_position(x = 'pop')

# Add letters for BALB/c
letters_w <- f_letters(p_liver, 'BALB/c')

# Add letters for C57BL/6
letters_b <- f_letters(p_liver2, 'C57BL/6')

# Combine letter df
liver_p <- bind_rows(letters_w, letters_b)

p_liver_h <- all_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(pop) %>%
  t_test(liver_wt_norm ~ host) %>%
  add_significance(symbols = c('####', '###', '##', '#', 'ns') ) %>%
  adjust_pvalue(method = "BH") %>%
  mutate(host = "C57BL/6") %>%
  filter(p.adj <= 0.05)

# Plot normalized liver weight
y_B <- all_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(host, pop) %>%
  summarize(max = max(liver_wt_norm)+0.5)

plot_liver <- 
  ggplot(all_df, aes(pop, liver_wt_norm)) +
  geom_boxplot(aes(fill = pop)) +
  facet_wrap(~host, strip.position = 'bottom') + 
  geom_vline(data=all_df[all_df$host=="BALB/c",], aes(xintercept=5.5), linetype="solid") +
  geom_vline(xintercept = 1.5, linetype = 'dashed') +
  scale_fill_manual(values = my_color) +
  ylim(0, 15) +
  geom_text(data = liver_p, aes(label = str_trim(letter), y = y_B$max)) +
  geom_text(data = p_liver_h, aes(label = p.signif), vjust = -0.5, y = 8) +
  geom_text(data = liver_res, aes(label = paste0(test, ', p = ',p)), y = 15) +
  ylab("Liver weight (% BW)") +
  my_theme() 
  

# Plot liver image 
#img_A <- image_ggplot(liverimage)

# This ensures that the image leaves some space at the edges
#theme(plot.margin = margin(t=1, l=1, r=1, b=1, unit = "cm"))

#++++++++++++++++++++++
# SPLEEN WEIGHT ####
#++++++++++++++++++++++
# Perform normality tests for worm counts
all_df %>%
  group_by(host, pop) %>%
  shapiro_test(spleen_wt_norm)

bartlett.test(data = subset(all_df, host == "BALB/c" & pop != 'Ctrl' ), spleen_wt_norm ~ pop)
#p-value = 0.04872
bartlett.test(data = subset(all_df, host == "C57BL/6" & pop != 'Ctrl'), spleen_wt_norm ~ pop)
#p-value = 0.04588

# Use ANOVA for both
# Perform stats test
spleen_res <- all_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(host) %>%
  anova_test(spleen_wt_norm ~ pop)

spleen_res$pop <- 'EG'

p_spleen <- all_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(host) %>%
  tukey_hsd(spleen_wt_norm ~ pop, p.adjust.method = 'BH')

spleen_p1 <- f_letters(p_spleen, 'BALB/c')
spleen_p2 <- f_letters(p_spleen, 'C57BL/6')

spleen_p <- bind_rows(spleen_p1, spleen_p2)

p_spleen_h <- all_df %>%
  group_by(pop) %>%
  t_test(spleen_wt_norm ~ host) %>%
  adjust_pvalue(method = 'BH') %>%
  add_significance(symbols = c('####', '###', '##', '#', 'ns')) %>%
  filter(p.adj <= 0.05)

# Plot spleen weight normalized
y_C <- all_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(host, pop) %>%
  summarize(max = max(spleen_wt_norm)+0.5)

plot_spleen <- 
  ggplot(all_df, aes(pop, spleen_wt_norm)) +
  geom_boxplot(aes(fill = pop)) +
  scale_fill_manual(values = my_color) +
  ylab("Spleen weight (% BW)") +
  ylim(0, 6.5) +
  facet_wrap(~host, strip.position = 'bottom') + 
  geom_vline(data=all_df[all_df$host=="BALB/c",], aes(xintercept =5.5)) +
  geom_vline(xintercept = 1.5, linetype="dashed") +
  geom_text(data = spleen_res, aes(label = paste0('ANOVA', ', p = ',p)), y = 6.4) +
  geom_text(data = spleen_p, aes(label = str_trim(letter), y = y_C$max)) +
  my_theme() 

# Plot spleen image 
#img_B <- image_ggplot(spleenimage) 

# This ensures that the image leaves some space at the edges
#theme(plot.margin = margin(t=1, l=1, r=1, b=1, unit = "cm"))

#++++++++++++++++++++++
#### Intestine length ####
#++++++++++++++++++++++
# Perform normality tests 
all_df %>%
  group_by(host, pop) %>%
  shapiro_test(intestine_length)

bartlett.test(data = subset(all_df, host == "BALB/c" & pop != 'Ctrl'), intestine_length ~ pop)
#p-value = 0.3103
bartlett.test(data = subset(all_df, host == "C57BL/6" & pop != 'Ctrl'), intestine_length ~ pop)
#p-value = 0.06166

# ANOVA for both

# Perform stats test
intestine_res <- all_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(host) %>%
  anova_test(intestine_length ~ pop)

# Artificial value for plotting
intestine_res$pop <- 'EG'
intestine_res$sig<- 'ns'

#p_intestine <- all_df %>%
#  filter(pop != 'Ctrl') %>%
#  group_by(host) %>%
#  tukey_hsd(intestine_length ~ pop)

#intestine_p1 <- f_letters(p_intestine, 'BALB/c')
#intestine_p2 <- f_letters(p_intestine, 'C57BL/6')

#intestine_p <- bind_rows(intestine_p1, intestine_p2)

p_intestine_h <- all_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(pop) %>%
  t_test(intestine_length ~ host) %>%
  adjust_pvalue(method = 'BH') %>%
  add_significance(symbols = c('####', '###', '##', '#', 'ns')) %>%
  filter(p.adj <= 0.05) %>%
  mutate(host = "C57BL/6")


# Plot intestine length normalized
#y_X <- all_df %>%
#  filter(host == "C57BL/6", pop != 'Ctrl') %>%
#  group_by(pop) %>%
#  summarize(max = max(intestine_length, na.rm = T)+10)

plot_S1 <- 
  ggplot(all_df, aes(pop, intestine_length)) +
  geom_boxplot(aes(fill = pop)) +
  scale_fill_manual(values = my_color) +
  ylab("Intestine length (mm)") +
  facet_wrap(~host, strip.position = 'bottom') + 
  geom_vline(data=all_df[all_df$host=="BALB/c",], aes(xintercept=5.5), linetype="solid") +
  geom_vline(aes(xintercept = 1.5), linetype="dashed") +
  geom_text(data = intestine_res, aes(label = paste0('ANOVA = ', format(p, 1)), y= 550)) +
  geom_text(data = p_intestine_h, aes(label = p.adj.signif), vjust = -0.5, y = 500) +
  my_theme() 

#++++++++++++++++++
# FIBROSIS #####
#++++++++++++++++++
# Perform normality tests
all_df %>%
  group_by(host, pop) %>%
  shapiro_test(fibrosis)

bartlett.test(data = subset(all_df, host == 'BALB/c' & pop != 'Ctrl'), fibrosis ~ pop)
#0.3517

# ANOVA for BALB/c, KW for C57BL/6

# Perform stats test
all_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(host) %>%
  anova_test(fibrosis ~ pop) 

anova_fibrosis <- data.frame(method = 'ANOVA', host = 'BALB/c', p = 0.073, pop = 'EG')

kw_fibrosis <- all_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(host) %>%
  kruskal_test(fibrosis ~ pop) %>%
  filter(host == 'C57BL/6') %>%
  mutate(pop = 'EG', method = 'K-W')

# Perform host comparisons

p_fibrosis_h <- all_df %>%
  filter(pop != 'Ctrl') %>%
  group_by(pop) %>%
  t_test(fibrosis ~ host) %>%
  add_significance(symbols = c('####', '###', '##', '#', 'ns')) %>%
  adjust_pvalue(method = 'BH') %>%
  filter(p <= 0.05) %>%
  mutate(host = "C57BL/6")

plot_fibrosis <- 
  ggplot(all_df, aes(pop,fibrosis)) +
  geom_boxplot(aes(fill = pop)) +
  facet_wrap(~host, strip.position = 'bottom') +
  ylab("Fibrotic area (%)") +
  ylim(0, 60)+
  geom_vline(data=stain_1[stain_1$host=="BALB/c",], aes(xintercept =5.5)) +
  geom_vline(xintercept =1.5, linetype="dashed") +
  geom_text(data = kw_fibrosis, aes(label = paste0(method, ', p = ', format(round(p,3)),1), y = 60)) +
  geom_text(data = anova_fibrosis, aes(label = paste0(method, ', p = ',p), y = 60)) +
  geom_text(data = p_fibrosis_h, aes(label = p.signif), vjust = -0.5, y = 4) +
  scale_fill_manual(values = my_color) +
  my_theme()
  
#ggsave("Plots/fibrosis.png", width = 8, height = 6)


#++++++++++++++++++
# GRANULOMAS ####
#++++++++++++++++++

# Overview of granulomas
df_granulomas %>%
  group_by(host,pop) %>%
  count()

# Randomly select x number of granulomas 
#set.seed(1231)

#df_rand <- df_granulomas %>%
#  group_by(host,pop) %>%
#  sample_n(80) %>%
#  ungroup()

# Perform normality tests
df_granulomas %>%
  group_by(host, pop) %>%
  shapiro_test(area)

# Perform stats
#df_rand %>%
kw_granuloma <- df_granulomas %>% 
  group_by(host) %>%
  kruskal_test(area ~ pop) %>%
  mutate(pop = 'EG',  method = 'K-W')

# A tibble: 2 × 7
#host    .y.       n statistic    df       p method        
#* <chr>   <chr> <int>     <dbl> <int>   <dbl> <chr>         
# 1 BALB/c  area    513      4.88     3 0.181   Kruskal-Wallis
#2 C57BL/6 area    523     16.2      3 0.00105 Kruskal-Wallis

p_gran <- df_granulomas %>%
  group_by(host) %>%
  dunn_test(area ~ pop, p.adjust.method = 'BH') %>%
  add_xy_position()

gran_p <- f_letters(p_gran, 'C57BL/6')

p_gran_h <- df_granulomas %>%
  group_by(pop) %>%
  wilcox_test(area ~ host) %>%
  adjust_pvalue(method = 'BH') %>%
  add_significance(symbols = c('####', '###', '##', '#', 'ns')) %>%
  filter(p.adj <= 0.05) %>%
  mutate(host = "C57BL/6")

# Plot granulomas
plot_granuloma <- 
  ggplot(df_granulomas, aes(pop,area)) +
  geom_boxplot(aes(fill = pop)) +
  facet_wrap(~host, strip.position = 'bottom') +
  ylab("Granuloma area (\u03BCm2)") +
  geom_vline(data=df_granulomas[df_granulomas$host=="BALB/c",], aes(xintercept=4.5)) +
  geom_text(data = p_gran_h, aes(label = p.adj.signif), vjust = -0.5, y = 4) +
  geom_text(data = kw_granuloma, aes(label = paste0(method, ", p = ", round(p,3)), y = 2.6E05) )+
  scale_fill_manual(values = my_color) +
  my_theme() +
  #theme(text = element_text(size =18)) +
  geom_text(data = gran_p, aes(label = letter), vjust = -0.5, y = 2.17E+05)  # +
  #add_pvalue(p_gran, tip.length = 0.01, label.size = 5)


ggsave('Plots/granuloma.png', width = 8, height =6)
#+++++++++++++++++++++++++++++++++++++
#### CORRELEATION EGGS VS FIBROSIS ####
#+++++++++++++++++++++++++++++++++++++

temp <- stain_1 %>%
  left_join(perfu, by = c("Sample" = "animal_ID"))

ggplot(subset(temp, pop.x != "Ctrl"), aes(liver_eggs, Value)) +
  geom_point(aes(color = pop.x)) +
  geom_smooth(method = "lm", se = FALSE, color = "grey") +
  facet_wrap(~host.x) +
  scale_color_manual(values = my_color, limits = force) +
  ylab("Fibrotic area") +
  xlab("Liver eggs") + 
  labs(color = 'Population') +
  stat_cor(method = "spearman", label.x = 5000, label.y = 12) +
  theme_minimal()

#+++++++++++++++++++++++
# EXPORT PLOTS ####
#+++++++++++++++++++++++

png("Figure3.png", width = 10, height = 12, unit = 'in', res =300)

(plot_weight + theme(legend.position = 'bottom')) / (plot_liver | plot_spacer()) / (plot_spleen | plot_spacer()) /
(plot_fibrosis | plot_granuloma) +
 plot_annotation(tag_levels = c("A")) &  theme(plot.tag = element_text(face = 'bold'))

dev.off()

#pdf("Figure3.pdf", width = 10, height = 12)

#ggarrange(ggarrange(plot_A, ncol = 1, labels = c("A.")),
          #ggarrange(plot_B, img_A, ncol = 2, labels = 'B.'),
          #ggarrange(plot_C, img_B, ncol = 2, labels = 'C.'),
          #ggarrange(plot_D, plot_E, ncol = 2, labels = c("D.", "E.")), nrow = 4)
#dev.off()


pdf("Additional file FigureS2.pdf", width = 5, height = 4)

(plot_S1) + #| plot_spacer()) / (plot_spacer()) / (plot_spacer()) + 
  plot_annotation(tag_levels = c('A')) & theme(plot.tag = element_text(face = 'bold'))
dev.off()