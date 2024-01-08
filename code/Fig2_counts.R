# Figure 2 - Parasite data
# Author: Kathrin Jutzeler
# Date: May 16, 2023
# Last updated
# R version 4.2.0, tidyverse version 1.3.2, ggplot2 3.3.6      ✔ purrr   0.3.4 
#✔ tibble  3.1.8      ✔ dplyr   1.0.10
#✔ tidyr   1.2.0      ✔ stringr 1.4.1 
#✔ readr   2.1.2      ✔ forcats 0.5.2 


# Load packages ####
library(tidyverse) # for everything
#library(readxl) # maybe not nee
library(rstatix) # stats
library(wesanderson) # colors
library(ggpubr) # stats/plots
library(ggprism) # add p-value to plot 
library(plotrix) # for standard error
library(patchwork) # for output layout

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
infctn <- read_csv("infection_data.csv")

infctn$pop <- as.factor(infctn$pop)
infctn$host <- as.factor(infctn$host)


#all_df <- read_csv("all_dfsion_data.csv")
all_df <- read_csv("all_data.csv")


all_df$pop <- as.factor(all_df$pop)
all_df$host <- as.factor(all_df$host)

# Functions ####

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

#++++++++++++++++++++++++
# INFECTION RATE ####
#++++++++++++++++++++++++
# Perform stats tests for infection rate
infctn %>%
  filter(pop != "Control") %>%
  group_by(host, pop) %>%
  shapiro_test(penrate)

p_infctn <- infctn %>%
  group_by(host) %>%
  kruskal_test(penrate ~ pop) %>%
  mutate(pop = 'EG', method = 'K-W')

# Between hosts 
infctn %>%
  filter(pop != 'Control') %>%
  group_by(pop) %>%
  wilcox_test(penrate ~ host) %>%
  adjust_pvalue(method = 'BH')

# Plot infection rate 
plot_infection <- 
  ggplot(subset(infctn, pop != "Control"), aes(pop, penrate)) +
  geom_boxplot(aes(fill = pop)) +
  geom_vline(data=infctn[infctn$host=="BALB/c",], aes(xintercept = 4.5), linetype="dashed") +
  scale_fill_manual(values = my_color, limits = force) +
  ylab("Penetration rate (%)") +
  facet_wrap(~host, strip.position = 'bottom') + 
  my_theme() +
  geom_text(data = p_infctn, aes(label = paste0(method, ", p = ", p), y =1.05)) +
  labs(fill = "Population")

#++++++++++++++++++++++++
# WORM COUNTS  ####
#++++++++++++++++++++++++

# Perform normality tests for worm counts
all_df %>%
  filter(pop != 'Control') %>%
  group_by(host, pop) %>%
  shapiro_test(total_worms)

bartlett.test(data = subset(all_df, host == "BALB/c" & pop != "Control"), total_worms ~ pop)
bartlett.test(data = subset(all_df, host == "C57BL/6" & pop != "Control"), total_worms ~ pop)

# Perform stats test
p_worm <- all_df %>%
  filter(pop != "Control") %>%
  group_by(host) %>%
  anova_test(total_worms ~ pop) 

p_worm$method <- 'ANOVA'
p_worm$pop <- 'EG'

# A tibble: 2 x 8
#host    Effect   DFn   DFd     F     p `p<.05`   ges
#* <fct>   <chr>  <dbl> <dbl> <dbl> <dbl> <chr>   <dbl>
#1 BALB/c  pop        3    16  5.40 0.009 "*"     0.503
#2 C57BL/6 pop        3    16  2.60 0.088 ""      0.328

p_worms <- all_df %>%
  filter(pop != "Control") %>%
  group_by(host) %>%
  tukey_hsd(total_worms ~ pop, p.adjust.methods = 'BH') %>%
  #filter(p.adj <= 0.05, host == "BALB/c") %>%
  add_x_position(x = 'pop')

# Add letters for BALB/c
letters_w <- f_letters(p_worms, 'BALB/c')

# Add letters for C57BL/6
letters_b <- f_letters(p_worms, 'C57BL/6')

# Combine letter df
worms_p <- bind_rows(letters_w, letters_b)

# Between hosts 
all_df %>%
  filter(pop != 'Control') %>%
  group_by(pop) %>%
  t_test(total_worms ~ host) %>%
  adjust_pvalue(method = 'BH')

# Plot total worms
plot_worms <- 
  ggplot(subset(all_df, pop != "Control"), aes(pop, total_worms)) +
  geom_boxplot(aes(fill = pop)) +
  scale_fill_manual(values = my_color, limits = force) +
  ylab("Total worms") +
  facet_wrap(~host, strip.position = 'bottom') + 
  geom_vline(data=all_df[all_df$host=="BALB/c",], aes(xintercept=4.5), linetype="dashed") +
  my_theme() +
  geom_text(data = worms_p, aes(label = str_trim(letter), y = 25)) +
  geom_text(data = p_worm, aes(label = paste0(method, ", p = ", p), y =30)) +
  labs(fill = "Population")

#++++++++++++++++++++++++
# EGG COUNTS #####
#++++++++++++++++++++++++
# Perform normality tests for egg counts
all_df %>%
  filter(pop != 'Control') %>%
  group_by(host, pop) %>%
  shapiro_test(total_eggs)

bartlett.test(data = subset(all_df, host == "BALB/c" & pop != "Control"), total_eggs ~ pop)
bartlett.test(data = subset(all_df, host == "C57BL/6" & pop != "Control"), total_eggs ~ pop)

# Perform stats test
eggs_res <- all_df %>%
  filter(pop != "Control") %>%
  group_by(host) %>%
  anova_test(total_eggs ~ pop)

eggs_res$method <- 'ANOVA'
eggs_res$pop <- 'EG'

# A tibble: 2 x 8
#host    Effect   DFn   DFd     F     p `p<.05`   ges
#* <fct>   <chr>  <dbl> <dbl> <dbl> <dbl> <chr>   <dbl>
#1 BALB/c  pop        3    14  3.42 0.047 "*"     0.423
#2 C57BL/6 pop        3    12  2.72 0.091 ""      0.405

p_eggs <- all_df %>%
  filter(pop != "Control") %>%
  group_by(host) %>%
  tukey_hsd(total_eggs ~ pop, p.adjust.methods = 'BH') %>%
  filter(p.adj <= 0.06, host == "BALB/c") 

# Between hosts 
all_df %>%
  filter(pop != 'Control') %>%
  group_by(pop) %>%
  t_test(total_eggs ~ host) %>%
  adjust_pvalue(method = 'BH')
  
  
# Plot total eggs
plot_C <- 
  ggplot(subset(all_df, pop != "Control"), aes(pop, total_eggs)) +
  geom_boxplot(aes(fill = pop)) +
  scale_fill_manual(values = my_color, limits = force) +
  ylab("Eggs per g of liver and intestine") +
  facet_wrap(~host, strip.position = 'bottom') + 
  geom_segment(data=all_df[all_df$host=="BALB/c",], aes(x=4.5, y=0, xend=4.5, yend=50000), linetype="dashed",
               color = "black", size=0.5) +
  my_theme() +
  add_pvalue(p_eggs, tip.length = 0.01, label = {"p.adj"}, label.size = 4,
                        step.group.by = 'host', y.position = 48000) +
  geom_text(data =eggs_res, aes(label = paste0(method, ", p = ", p), y =52000)) 

# Perform normality tests for liver eggs only
all_df %>%
  filter(pop != 'Control') %>%
  group_by(host, pop) %>%
  shapiro_test(liver_eggs)

bartlett.test(data = subset(all_df, host == "BALB/c" & pop != "Control"), total_eggs ~ pop)
bartlett.test(data = subset(all_df, host == "C57BL/6" & pop != "Control"), total_eggs ~ pop)

# Perform stats test
all_df %>%
  filter(pop != "Control") %>%
  group_by(host) %>%
  anova_test(liver_eggs ~ pop)

# A tibble: 2 × 8
#host    Effect   DFn   DFd     F     p `p<.05`   ges
#* <chr>   <chr>  <dbl> <dbl> <dbl> <dbl> <chr>   <dbl>
#1 BALB/c  pop        3    15  6.44 0.005 *       0.563
#2 C57BL/6 pop        3    12  4.99 0.018 *       0.555

p_liver_eggs <- all_df %>%
  filter(pop != 'Control') %>%
  group_by(host) %>%
  tukey_hsd(liver_eggs ~ pop, p.adjust.methods = 'BH') %>%
  filter(p.adj <= 0.06)

# Between hosts 
all_df %>%
  filter(pop != 'Control') %>%
  group_by(pop) %>%
  t_test(total_eggs ~ host)


# Plot liver eggs

ggplot(subset(all_df, pop != "Control"), aes(pop, liver_eggs)) +
  geom_boxplot(aes(fill = pop)) +
  scale_fill_manual(values = my_color, limits = force) +
  ylab("Eggs per g of liver") +
  facet_wrap(~host, strip.position = 'bottom') + 
  geom_segment(data=all_df[all_df$host=="BALB/c",], aes(x=4.5, y=0, xend=4.5, yend=35000), linetype="dashed",
               color = "black", size=0.5) +
  theme_minimal() + 
  my_theme() +
  theme(text = element_text(size = 18)) +
  add_pvalue(p_liver_eggs, tip.length = 0.01, label = {"p.adj.signif"}, label.size = 5,
             step.group.by = 'host', y.position = 31000) 

ggsave('Plots/livereggs.png', width = 8, height =6)

#++++++++++++++++++++++++
# FECUNDITY #####
#++++++++++++++++++++++++
# Perform normality tests for egg counts
all_df %>%
  filter(pop != 'Control') %>%
  group_by(host, pop) %>%
  shapiro_test(fecundity)

# KW for BALB/c ANOVA for C57Bl/6

# Perform stats test
all_df %>%
  filter(pop != "Control") %>%
  group_by(host) %>%
  kruskal_test(fecundity ~ pop) 

all_df %>%
  filter(pop != "Control") %>%
  group_by(host) %>%
  anova_test(fecundity ~ pop) 

fec_res <- data.frame(host = c('BALB/c', 'C57BL/6'), pop = 'EG', method = c('K-W', 'ANOVA'),
                      p = c(0.0291, 0.102))

# A tibble: 2 x 7
#host    .y.           n statistic    df      p method        
#* <fct>   <chr>     <int>     <dbl> <int>  <dbl> <chr>         
#1 BALB/c  fecundity    20      9.02     3 0.0291 Kruskal-Wallis

p_fecundity <- all_df %>%
  filter(pop != "Control") %>%
  group_by(host) %>%
  dunn_test(fecundity ~ pop, p.adjust.method = 'BH') 
  #filter(p.adj <= 0.05)

# Add letters for BALB/c
letters_w <- f_letters(p_fecundity, 'BALB/c')

# Add letters for C57BL/6
letters_b <- f_letters(p_fecundity, 'C57BL/6')

# Combine letter df
fecundity_p <- bind_rows(letters_w, letters_b)


# Between hosts 
all_df %>%
  filter(pop != 'Control') %>%
  group_by(pop) %>%
  wilcox_test(fecundity ~ host) %>%
  adjust_pvalue(method = 'BH')
  
# Plot fecundity
plot_fecundity <- 
  ggplot(subset(all_df, pop != "Control"), aes(pop, fecundity/1000)) +
  geom_boxplot(aes(fill = pop)) +
  scale_fill_manual(values = my_color, limits = force) +
  ylab("Fecundity (x 1000)") +
  ylim(0, 25) +
  facet_wrap(~host, strip.position = 'bottom') + 
  geom_vline(data=fec[fec$host=="BALB/c",], aes(xintercept =4.5), linetype="dashed") +
  my_theme() +
  geom_text(data =fec_res, aes(label = paste0(method, ", p = ", round(p,3))), y = 25)  +
  geom_text(data = fecundity_p, aes(label = str_trim(letter), y = 13)) 
  #add_pvalue(p_fecundity, tip.length = 0.01, label = {"p.adj.signif"}, label.size = 5,
   #          step.group.by = 'host', y.position = 13) 

#++++++++++++++++++++++++
# EGG DISTRIBUTION ####
#++++++++++++++++++++++++
# Create new data frame
all_df_eggs <- all_df %>%
  filter(pop != 'Control') %>%
  gather("organ", "eggs", c(liver_eggs, intestine_eggs )) %>%
  droplevels()

# Perform normality tests
all_df_eggs %>%
  group_by(host, pop, organ) %>%
  shapiro_test(eggs)

bartlett.test(data = subset(all_df, host == "BALB/c" & pop != "Control"), liver_eggs ~ pop)
bartlett.test(data = subset(all_df, host == "C57BL/6" & pop != "Control"), liver_eggs ~ pop)

# Perform stats test
all_df_eggs %>%
  group_by(pop, host) %>%
  wilcox_test(eggs ~ organ) %>%
  adjust_pvalue(method = "BH") 

#p_organ$group1 = c('LE','OR')
#p_organ$group2 = c('LE', 'OR')

# Between hosts 
all_df_eggs %>%
  group_by(pop, organ) %>%
  wilcox_test(eggs ~ host) %>%
  adjust_pvalue(method = "BH") 

# Between populations (only liver eggs significant)
all_df_eggs %>%
  group_by(host, organ) %>%
  anova_test(eggs ~ pop) 
  
anova_res <- data.frame(host = c('BALB/c', 'C57BL/6'), 
           method = c('ANOVA', 'ANOVA'),
           result = c(0.005, 0.018),
           pop = c('EG', 'EG'))

p_pop <- all_df_eggs %>%
  group_by(host, organ) %>%
  tukey_hsd(eggs ~ pop,  p.adjust.methods = 'BH') %>%
  filter(organ == 'liver_eggs') %>%
  add_xy_position(x = "pop")

# Add letters for BALB/c
letters_w <- f_letters(p_pop, 'BALB/c')

# Add letters for C57BL/6
letters_b <- f_letters(p_pop, 'C57BL/6')

# Combine letter df
pop_p <- bind_rows(letters_w, letters_b)

    
## Plot egg distribution
# Pop on the X axis
plot_eggs <- 
  ggplot(all_df_eggs, aes(pop, eggs)) +
  geom_boxplot(aes(fill = organ)) +
  scale_fill_manual(values = c("wheat","indianred"), label = c('Intestine', 'Liver')) +
  facet_wrap(.~host, strip.position = 'bottom') + 
  geom_vline(data=all_df[all_df$host=="BALB/c",], aes(xintercept  = 4.5), linetype = 'dashed') +
  geom_text(data = anova_res, aes(label = paste0(method,", p = ",result)), 
            y = 33500, color = 'indianred') +
  ylim(0, 34000) +
  my_theme() +
  geom_text(data = pop_p, aes(label = str_trim(letter), y = 28000),color = 'indianred') +
  theme(legend.position = 'right') +
  labs(fill = 'Organ',  y = "Eggs per g of tissue")


# Plot egg ratio ####
# Pop on the X axis
all_df <- all_df %>%
  mutate(ratio = liver_eggs/intestine_eggs)

# Perform normality tests
all_df %>%
  filter(pop != 'Control') %>%
  group_by(host, pop) %>%
  shapiro_test(ratio)

# ANOVA for C57BL/6 and KW for BALB/c
bartlett.test(data = subset(all_df, host == "C57BL/6" & pop != "Control"), liver_eggs ~ pop)

all_df %>%  group_by(host) %>%
  anova_test(ratio ~ pop) 

all_df %>%  group_by(host) %>%
  kruskal_test(ratio ~ pop) 

all_df %>% 
  #group_by(pop) %>%
  wilcox_test(ratio ~ host)

p_ratio <- data.frame(pop = 'EG', p = c(0.235, 0.306 ), method = c('K-W', 'ANOVA'), host =
                        c('BALB/c', 'C57BL/6')) 

plot_ratio <- 
  ggplot(subset(all_df, pop != 'Control'), aes(pop, ratio)) +
  geom_boxplot(aes(fill = pop)) +
  scale_fill_manual(values = my_color) +
  facet_wrap(.~host, strip.position = 'bottom') + 
  geom_vline(data=all_df[all_df$host=="BALB/c",], aes(xintercept=4.5), linetype="solid") +
  geom_text(data = p_ratio, aes(label = paste0(method, ", p = ", p), y = 10) )+
  my_theme() +
  labs(y = 'Liver eggs:Intestine eggs')

# organ on the x-axis 
plot_F <- ggplot(all_df_eggs, aes(organ, eggs)) +
  geom_boxplot(aes(fill = pop)) +
  facet_wrap(.~host, strip.position = 'bottom') + 
  scale_fill_manual(values = my_color) +
  scale_color_manual(values = my_color, guide = "none") +
  scale_x_discrete(labels = c("intestine_eggs"= "Intestine", "liver_eggs" = "Liver")) +
  geom_segment(data=all_df[all_df$host=="BALB/c",], aes(x=2.5, y=0, xend=2.5, yend=32000), linetype="dashed",
               color = "black", size=0.5) +
  my_theme() +
  stat_pvalue_manual(p_organ2, y.position = c(31000, 32000), label.size =  5, tip.length = 0, xmin = 'xmin', xmax = 'xmax', 
             color = "pop") +
  add_pvalue(p_pop2, y.position = c(29000, 30100), label.size =  4, tip.length = 0.01, xmin = 'xmin', 
             xmax = 'xmax') +
  labs(fill = 'Population', y = "Eggs per g of tissue")

#+++++++++++++++++++++
# Export plots ####
#+++++++++++++++++++++

pdf("Figure2.pdf", width = 10, height = 12)

(plot_infection | plot_worms) / plot_eggs / (plot_fecundity  | plot_spacer()) + plot_annotation(tag_levels = c("A")) &
#((plot_A | plot_B)  + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')) / 
  #plot_F / (plot_D | plot_spacer()  )  + plot_annotation(tag_levels = c("A"))  & theme(plot.tag = element_text(face = 'bold'))
#(plot_F | plot_D)  + plot_annotation(tag_levels = c("A"))
   theme(plot.tag = element_text(face = 'bold'))

dev.off()
  

pdf("Additional file FigureS1.pdf", width = 5, height = 4)

(plot_ratio) +# | plot_spacer()) / (plot_spacer()) / (plot_spacer())  
  plot_annotation(tag_levels = c("A")) &
  theme(plot.tag = element_text(face = 'bold'))

dev.off()