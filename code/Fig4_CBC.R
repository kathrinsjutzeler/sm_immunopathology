# Figure 4 - CBC 
# Author: Kathrin Jutzeler
# Date: May 23, 2023
# Last updated: September 7, 2023

# Load packages ####
library(tidyverse)
library(rstatix)  # Pipe frindly stats 
library(wesanderson)  # Color palettes
library(ggpubr) # ggplot add-on 
library(ggprism) # add_pvalue to plots
library(emmeans) # Mean adjustment for linear models 
library(multcomp) # Perform mutliple comparisons 
library(multcompView)
library(plotrix) # for standard error
library(patchwork) # to arrange plots
library(PMCMRplus) # for Conover test
library(rcompanion) # full ptable
library(ggrepel)

# Define theme and colors ####
my_theme <- function(){theme_minimal() + theme(panel.grid = element_blank(), legend.position = "none", strip.placement = 'outside',
                             axis.line = element_line(1),
                             panel.spacing = unit(0,'lines'), 
                             text = element_text(size = 14)) }

my_color = c("Control" = wes_palette("Darjeeling2")[1],
             "LE" = wes_palette("Darjeeling2")[2],
             "EG" = wes_palette("Darjeeling2")[3],
             "OR" = "#710193",
             "BRE" = wes_palette("Darjeeling2")[5])


# Import the data ####
CBC_all <- read_csv("CBC_data.csv")

CBC_all$time <- factor(CBC_all$time, levels = c("Baseline", "Week 2", "Week 4", "Week 6", "Week 8", "Week 10"))
CBC_all$pop <- factor(CBC_all$pop, levels = c('Control', 'BRE', 'EG', 'LE', 'OR'))
CBC_all$sample <- factor(CBC_all$sample)

#CBC_all_gather <- CBC_all %>%
#  gather("cell_type", "concentration", 3:14 )

# Remove mice with missing values
CBC_refined <- CBC_all %>%
  filter(!sample %in% c('2D', '10E'))

# Define parameters
parameters <- c('LYMPH#(K/uL)','LYMPH%(%)','RET#(K/uL)', 'RET%(%)', 'NEUT#(K/uL)', 'NEUT%(%)',
                'MONO#(K/uL)','MONO%(%)', 'EO#(K/uL)', 'EO%(%)', 'HCT(%)', 'BASO#(K/uL)', 'BASO%(%)')

# Add summary stats for plotting
CBC_summary <- CBC_refined %>%
  group_by(host, pop, time) %>%
  summarize_at(vars(parameters), c(mean = mean, se = std.error)) %>%
  ungroup()

# Perform normality tests 
norms_CBC <- CBC_refined %>%
  filter(pop != 'Control') %>%
  group_by(host,pop, time) %>%
  summarize_at(parameters, list(shapiro = ~shapiro.test(.)$p.value) ) 

#++++++++++++++++++++++++++++++++++++++
# Function for Conover and letters ####
#++++++++++++++++++++++++++++++++++++++

f_conover <- function(parameter, mouse){
  df <- filter(CBC_summary, host == mouse, pop != 'Control') %>% droplevels()
  
  z <- df[[parameter]]
  
  res <- PMCMRplus::frdAllPairsConoverTest(y = z, groups = df$pop, 
                                           blocks = df$time, p.adjust.method = "BH")
  
  res.p <- fullPTable(res$p.value)
  res.l <- multcompLetters(res.p)
  res.df <- as.data.frame(res.l$Letters)
  colnames(res.df) <- 'letter'
  
  #res.df$pop <- rownames(res.df)
  #res.df$host <- mouse
  
  return(res.df)
}

#+++++++++++++++++++++++++
# Plot lymphocytes in % ####
#+++++++++++++++++++++++++

# Perform stats
# Repeated measures 
f_LYMPH_per <- CBC_summary %>%
  group_by(host) %>%
  friedman_test(`LYMPH%(%)_mean` ~ pop | time) %>%
  p_round(p, 3)

# Artificial value for plotting
f_LYMPH_per$pop <- c(1, 2)
f_LYMPH_per$time <- c(2,2)

lymph_b <- f_conover("LYMPH%(%)_mean", "C57BL/6")
lymph_w <- f_conover("LYMPH%(%)_mean", "BALB/c")

df <-  bind_rows(lymph_b, lymph_w)

stat <- CBC_summary %>%
  filter(time == 'Week 10') %>%
  bind_cols( df)

# Plot lymphocytes
plot_A <- 
  ggplot(CBC_summary, aes(time, `LYMPH%(%)_mean`, group = interaction(pop))) +
  geom_point(size =1.5, position = position_dodge(width = 0.5), aes(color = pop)) +
  geom_line(position = position_dodge(width = 0.5), aes(color = pop)) +
  geom_errorbar(aes(ymin=`LYMPH%(%)_mean`-`LYMPH%(%)_se`, ymax=`LYMPH%(%)_mean`+`LYMPH%(%)_se`, 
                    color = pop), width=.2,
                position=position_dodge(0.5)) +
  facet_grid(host~., scale ='fixed') +
  geom_segment(data=CBC_all[CBC_all$host=="BALB/c",], aes(x=0, y=35, xend=6.5, yend=35), linetype="dashed",
               color = "black", linewidth=0.5) +
  ylab("Mean lymphocytes (%)") +
  labs(color = "Population") +
  scale_color_manual(values = my_color) +
  scale_x_discrete(labels=c("0","2","4","6", "8", "10")) +
  my_theme() +
  geom_text(data = f_LYMPH_per, aes(label = paste0('Friedman p = ', round(p, digits =10)), y =40)) +
  geom_text_repel(
    aes(label = letter), data = stat, position = position_dodge(width = 0.5), hjust = 0.5)

#++++++++++++++++++++++++++++++++
# Plot lymphocytes in units ####
#++++++++++++++++++++++++++++++++

# Repeated measures 
f_LYMPH_uni <- CBC_summary %>%
  filter(pop != 'Control') %>%
  group_by(host) %>%
  friedman_test(`LYMPH#(K/uL)_mean` ~ pop | time) %>%
  p_round(p, 3)

# Repeated measures for host
CBC_summary %>%
  filter(pop != 'Control') %>%
  group_by(pop) %>%
  friedman_test(`LYMPH#(K/uL)_mean` ~ host | time) %>%
  p_round(p, 3)

CBC_summary %>%
  filter(pop != 'Control') %>%
  group_by(time, pop) %>%
  wilcox_test(`LYMPH#(K/uL)_mean` ~ host) 

# Artificial value for plotting
f_LYMPH_uni$pop <- c(1, 2)
f_LYMPH_uni$time <- c(2,2)

lymph_b <- f_conover("LYMPH#(K/uL)_mean", "C57BL/6")
#lymph_w <- f_conover("LYMPH#(K/uL)_mean", "BALB/c")

#df <-  bind_rows(lymph_b, lymph_w)

stat <- CBC_summary %>%
  filter(time == 'Week 10', host == "C57BL/6", pop != 'Control') %>%
  bind_cols(lymph_b)

# Between hosts 
p_lymph_h1 <- CBC_refined %>%
  filter(pop != 'Control') %>%
  group_by( time) %>%
  wilcox_test(`LYMPH#(K/uL)` ~ host) %>%
  adjust_pvalue(method = 'BH') %>%
  add_significance(symbols = c('####', '###', '##', '#', 'ns')) %>%
  filter(p.adj <= 0.05) %>%
  mutate(host = "C57BL/6", pop = 'OR')

# Plot 
plot_lymph <- 
  ggplot(CBC_summary, aes(time, `LYMPH#(K/uL)_mean`, group = interaction(pop))) +
  geom_point(size =1.5, position = position_dodge(width = 0.5), aes(color = pop)) +
  geom_line(data=subset(CBC_summary, pop != 'Control'),position = position_dodge(width = 0.5), aes(color = pop)) +
  geom_line(data=subset(CBC_summary, pop == 'Control'),position = position_dodge(width = 0.5), aes(color = pop),
            linetype = 'dashed') +
  geom_errorbar(aes(ymin=`LYMPH#(K/uL)_mean`-`LYMPH#(K/uL)_se`, ymax=`LYMPH#(K/uL)_mean`+`LYMPH#(K/uL)_se`, 
                    color = pop), width=.2,
                position=position_dodge(0.5)) +
  facet_grid(host~., scale ='fixed') +
  geom_segment(data=CBC_all[CBC_all$host=="BALB/c",], aes(x=0, y=3, xend=6.5, yend=3), linetype="dashed",
               color = "black", size=0.5) +
  ylab("Mean lymphocytes (K/uL)") +
  xlab('Weeks post infection') +
  labs(color = "Population") +
  geom_text(data = f_LYMPH_uni, aes(label = paste0('Friedman p = ', round(p, digits =10)), y =16)) +
  geom_text(data = p_lymph_h1, aes(label = p.adj.signif), vjust = -0.5, y = 4) +
  geom_text_repel(
    aes(label = letter, color = pop), show.legend = F, data = stat, position = position_dodge(width = 0.5), hjust = 0.5) +
  scale_color_manual(values = my_color) +
  scale_x_discrete(labels=c("0","2","4","6", "8", "10")) +
  theme_minimal() +
  my_theme()


#++++++++++++++++++++++++++++
# Plot eosinophils in % ####
#++++++++++++++++++++++++++++

# Repeated measures 
f_EO_per <- CBC_summary %>%
  group_by(host) %>%
  friedman_test(`EO%(%)_mean` ~ pop | time) %>%
  p_round(p, 3)

# Artificial value for plotting
f_EO_per$pop <- c(1, 2)
f_EO_per$time <- c(2,2)

eo <- f_conover('EO%(%)_mean', 'BALB/c')

stat_eo <- CBC_summary %>%
  filter(time == 'Week 10', host == 'BALB/c') %>%
  bind_cols(eo)

plot_C <- 
  ggplot(CBC_summary, aes(time, `EO%(%)_mean`, group = interaction(pop))) +
  geom_point(size =1.5, position = position_dodge(width = 0.5), aes(color = pop)) +
  geom_line(data=subset(CBC_summary, pop != 'Control'),position = position_dodge(width = 0.5), aes(color = pop)) +
  geom_line(data=subset(CBC_summary, pop == 'Control'),position = position_dodge(width = 0.5), aes(color = pop),
            linetype = 'dashed') +
  geom_errorbar(aes(ymin=`EO%(%)_mean`-`EO%(%)_se`, ymax=`EO%(%)_mean`+`EO%(%)_se`, 
                    color = pop), width=.2,
                position=position_dodge(0.5)) +
  facet_grid(host~., scale ='fixed') +
  geom_segment(data=CBC_all[CBC_all$host=="BALB/c",], aes(x=0, y=0, xend=6.5, yend=0), linetype="dashed",
               color = "black", size=0.5) +
  geom_text(data = f_EO_per, aes(label = paste0('Friedman p = ', round(p, digits =10)), y =5)) +
  geom_text_repel(
    aes(label = letter), data = stat_eo, position = position_dodge(width = 0.5), hjust = 0.5) +
  ylab("Mean eosinophils (%)") +
  labs(color = "Population") +
  scale_color_manual(values = my_color) +
  scale_x_discrete(labels=c("0","2","4","6", "8", "10")) +
  theme_minimal() +
  my_theme()

#++++++++++++++++++++++++++++++++
# Plot eosinophils in units ####
#++++++++++++++++++++++++++++++++

# Repeated measures 
f_EO_uni <- CBC_summary %>%
  filter(pop != 'Control') %>%
  group_by(host) %>%
  friedman_test(`EO#(K/uL)_mean` ~ pop | time) %>%
  p_round(p, 3)

##### Trying stuff=================================================
# Repeated measures for host
temp <- CBC_refined %>%
  filter(pop != 'Control') %>%
  group_by(host, time) %>%
  summarize_at(vars(parameters), c(mean = mean, se = std.error)) %>%
  ungroup()

CBC_summary %>%
  filter(pop != 'Control') %>%
  group_by(pop) %>%
  friedman_test(`MONO#(K/uL)_mean` ~ host | time) %>%
  p_round(p, 3)

temp <- filter(CBC_summary, pop != 'Control', host == 'C57BL/6') %>% droplevels()

#PMCMRplus::frdAllPairsConoverTest(y = temp$`LYMPH#(K/uL)_mean`, groups = temp$pop, 
#                                  blocks = temp$time), p.adjust.method = "BH")

##### End=================================================
# Artificial value for plotting
f_EO_uni$pop <- c(1, 2)
f_EO_uni$time <- c(2,2)

eo_w <- f_conover("EO#(K/uL)_mean", "BALB/c")
eo_b <- f_conover("EO#(K/uL)_mean", "C57BL/6")

df <-  bind_rows(eo_w, eo_b)

stat <- CBC_summary %>%
  filter(pop != 'Control', time == 'Week 10') %>%
  bind_cols( df)

# Between hosts 
p_eo_h1 <- CBC_refined %>%
  filter(pop != 'Control') %>%
  group_by(time) %>%
  wilcox_test(`EO#(K/uL)` ~ host) %>%
  adjust_pvalue(method = 'BH') %>%
  add_significance(symbols = c('####', '###', '##', '#', 'ns')) %>%
  filter(p.adj <= 0.05) %>%
  mutate(host = "C57BL/6", pop = 'OR')

# Plot 
plot_eo <- 
  ggplot(CBC_summary, aes(time, `EO#(K/uL)_mean`, group = interaction(pop))) +
  geom_point(size =1.5, position = position_dodge(width = 0.5), aes(color = pop)) +
  geom_line(data=subset(CBC_summary, pop != 'Control'),position = position_dodge(width = 0.5), aes(color = pop)) +
  geom_line(data=subset(CBC_summary, pop == 'Control'),position = position_dodge(width = 0.5), aes(color = pop),
            linetype = 'dashed') +
  geom_errorbar(aes(ymin=`EO#(K/uL)_mean`-`EO#(K/uL)_se`, ymax=`EO#(K/uL)_mean`+`EO#(K/uL)_se`, 
                    color = pop), width=.2,
                position=position_dodge(0.5)) +
  facet_grid(host~., scale ='fixed') +
  geom_segment(data=CBC_all[CBC_all$host=="BALB/c",], aes(x=0, y=0, xend=6.5, yend=0), linetype="dashed",
               color = "black", size=0.5) +
  geom_text(data = f_EO_uni, aes(label = paste0('Friedman p = ', round(p, digits =10)), y =1.5)) +
  geom_text(data = p_eo_h1, aes(label = p.adj.signif), vjust = -0.5, y = 1.5) +
  ylab("Mean eosinophils (K/uL)") +
  xlab("Weeks post infection") +
  labs(color = "Population") +
  scale_color_manual(values = my_color) +
  scale_x_discrete(labels=c("0","2","4","6", "8", "10")) +
  my_theme()

  
#+++++++++++++++++++++++
# Plot monocytes in % 
#+++++++++++++++++++++++

# Repeated measures 
f_MONO_per <- CBC_summary %>%
  group_by(host) %>%
  friedman_test(`MONO%(%)_mean` ~ pop | time) %>%
  p_round(p, 3)

# Artificial value for plotting
f_MONO_per$pop <- c(1, 2)
f_MONO_per$time <- c(2,2)

mono_w <- f_conover('MONO%(%)_mean', 'BALB/c')
mono_b <- f_conover('MONO%(%)_mean', 'C57BL/6')

mono <- bind_rows(mono_w, mono_b)

stat_mono <- CBC_summary %>%
  filter(time == 'Week 10') %>%
  bind_cols(mono)

plot_E <- 
  ggplot(CBC_summary, aes(time, `MONO%(%)_mean`, group = interaction(pop))) +
  geom_point(size =1.5, position = position_dodge(width = 0.5), aes(color = pop)) +
  geom_line(position = position_dodge(width = 0.5), aes(color = pop)) +
  geom_errorbar(aes(ymin=`MONO%(%)_mean`-`MONO%(%)_se`, ymax=`MONO%(%)_mean`+`MONO%(%)_se`, 
                    color = pop), width=.2,
                position=position_dodge(0.5)) +
  facet_grid(host~., scale ='fixed') +
  geom_segment(data=CBC_all[CBC_all$host=="BALB/c",], aes(x=0, y=0, xend=6.5, yend=0), linetype="dashed",
               color = "black", size=0.5) +
  geom_text(data = f_MONO_per, aes(label = paste0('Friedman p = ', round(p, digits =10)), y =5)) +
  geom_text_repel(
    aes(label = letter), data = stat_mono, position = position_dodge(width = 0.5), hjust = 0.5) +
  ylab("Mean monocytes (%)") +
  labs(color = "Population") +
  scale_color_manual(values = my_color) +
  scale_x_discrete(labels=c("0","2","4","6", "8", "10")) +
  theme_minimal() +
  my_theme()


#+++++++++++++++++++++++++++++
# Plot monocytes in unit ####
#+++++++++++++++++++++++++++++

# Repeated measures 
f_MONO_uni <- CBC_summary %>%
  filter(pop != 'Control') %>%
  group_by(host) %>%
  friedman_test(`MONO#(K/uL)_mean` ~ pop | time) %>%
  p_round(p, 3)

# Repeated measures for host
CBC_summary %>%
  filter(pop != 'Control') %>%
  group_by(pop) %>%
  friedman_test(`NEUT#(K/uL)_mean` ~ host | time) %>%
  p_round(p, 3)

# Artificial value for plotting
f_MONO_uni$pop <- c(1, 2)
f_MONO_uni$time <- c(2,2)

mono_w <- f_conover('MONO%(%)_mean', 'BALB/c')
mono_b <- f_conover('MONO%(%)_mean', 'C57BL/6')

mono <- bind_rows(mono_w, mono_b)

stat_mono <- CBC_summary %>%
  filter(time == 'Week 10', pop != 'Control') %>%
  bind_cols(mono)

# Between hosts 
p_mono_h1 <- CBC_refined %>%
  filter(pop != 'Control') %>%
  group_by(time) %>%
  wilcox_test(`MONO#(K/uL)` ~ host) %>%
  adjust_pvalue(method ='BH') %>%
  add_significance(symbols = c('####', '###', '##', '#', 'ns')) %>%
  filter(p.adj <= 0.05) %>%
  mutate(host = "C57BL/6", pop = 'OR')

# Plot
plot_mono <- 
  ggplot(CBC_summary, aes(time, `MONO#(K/uL)_mean`, group = interaction(pop))) +
  geom_point(size =1.5, position = position_dodge(width = 0.5), aes(color = pop)) +
  geom_line(data=subset(CBC_summary, pop != 'Control'),position = position_dodge(width = 0.5), aes(color = pop)) +
  geom_line(data=subset(CBC_summary, pop == 'Control'),position = position_dodge(width = 0.5), aes(color = pop),
            linetype = 'dashed') +
  geom_errorbar(aes(ymin=`MONO#(K/uL)_mean`-`MONO#(K/uL)_se`, ymax=`MONO#(K/uL)_mean`+`MONO#(K/uL)_se`, 
                    color = pop), width=.2,
                position=position_dodge(0.5)) +
  facet_grid(host~., scale ='fixed') +
  geom_segment(data=CBC_all[CBC_all$host=="BALB/c",], aes(x=0, y=0, xend=6.5, yend=0), linetype="dashed",
               color = "black", size=0.5) +
  geom_text(data = f_MONO_uni, aes(label = paste0('Friedman p = ', round(p, digits =10)), y =2)) +
  geom_text(data = p_mono_h1, aes(label = p.adj.signif), vjust = -0.5, y = 1.5) +
  ylab("Mean monocytes (K/uL)") +
  xlab('Weeks post infection') +
  labs(color = "Population") +
  scale_color_manual(values = my_color) +
  scale_x_discrete(labels=c("0","2","4","6", "8", "10")) +
  my_theme()

#++++++++++++++++++++++++++++
# Plot reticulocytes in % ####
#+++++++++++++++++++++++++++

# Perform stats
CBC_all %>%
  group_by(host, time) %>%
  kruskal_test(`RET%(%)` ~ pop)

CBC_all %>%
  group_by(host, time) %>%
  dunn_test(`RET%(%)` ~ pop) %>%
  filter(p.adj <= 0.05)

plot_G <- 
  ggplot(CBC_summary, aes(time, `RET%(%)_mean`, group = interaction(pop))) +
  geom_point(size =1.5, position = position_dodge(width = 0.5), aes(color = pop)) +
  geom_line(position = position_dodge(width = 0.5), aes(color = pop)) +
  geom_errorbar(aes(ymin=`RET%(%)_mean`-`RET%(%)_se`, ymax=`RET%(%)_mean`+`RET%(%)_se`, 
                    color = pop), width=.2,
                position=position_dodge(0.5)) +
  facet_grid(host~., scale ='fixed') +
  geom_segment(data=CBC_all[CBC_all$host=="BALB/c",], aes(x=0, y=0, xend=6.5, yend=0), linetype="dashed",
               color = "black", size=0.5) +
  stat_compare_means(label = 'p.signif', vjust = 5, data = CBC_all, aes(time,`RET%(%)`)) +
  ylab("Mean reticulocytes (%)") +
  xlab("Weeks post infection") +
  labs(color = "Population") +
  scale_color_manual(values = my_color) +
  scale_x_discrete(labels=c("0","2","4","6", "8", "10")) +
  theme_minimal() +
  my_theme()

#++++++++++++++++++++++++++++
# Plot reticulocytes in unit ####
#++++++++++++++++++++++++++++

# Repeated measures 
f_RET_uni <- CBC_summary %>%
  filter(pop != 'Control') %>%
  group_by(host) %>%
  friedman_test(`RET#(K/uL)_mean` ~ pop | time) %>%
  p_round(p, 3)

# Artificial value for plotting
f_RET_uni$pop <- c(1, 2)
f_RET_uni$time <- c(2,2)

ret_w <- f_conover('RET#(K/uL)_mean', 'BALB/c')
ret_b <- f_conover('RET#(K/uL)_mean', 'C57BL/6')

ret <- bind_rows(ret_w, ret_b)

stat_ret<- CBC_summary %>%
  filter(time == 'Week 10') %>%
  bind_cols(ret)


# Between hosts 
p_ret_h1 <- CBC_refined %>%
  filter(pop != 'Control') %>%
  group_by(time) %>%
  wilcox_test(`RET#(K/uL)` ~ host) %>%
  add_significance(symbols = c('####', '###', '##', '#', 'ns')) %>%
  filter(p <= 0.05) %>%
  mutate(host = "C57BL/6", pop = 'OR')

# Plot
plot_ret <- 
  ggplot(CBC_summary, aes(time, `RET#(K/uL)_mean`, group = interaction(pop))) +
  geom_point(size =1.5, position = position_dodge(width = 0.5), aes(color = pop)) +
  geom_line(data=subset(CBC_summary, pop != 'Control'),position = position_dodge(width = 0.5), aes(color = pop)) +
  geom_line(data=subset(CBC_summary, pop == 'Control'),position = position_dodge(width = 0.5), aes(color = pop),
            linetype = 'dashed') +
  geom_errorbar(aes(ymin=`RET#(K/uL)_mean`-`RET#(K/uL)_se`, ymax=`RET#(K/uL)_mean`+`RET#(K/uL)_se`, 
                    color = pop), width=.2,
                position=position_dodge(0.5)) +
  facet_grid(host~., scale ='fixed') +
  geom_segment(data=CBC_all[CBC_all$host=="BALB/c",], aes(x=0, y=0, xend=6.5, yend=0), linetype="dashed",
               color = "black", size=0.5) +
  geom_text(data = f_RET_uni, aes(label = paste0('Friedman p = ', round(p, digits =10)), y =1500)) +
  ylab("Mean reticulocytes (K/uL)") +
  xlab('Weeks post infection') +
  labs(color = "Population") +
  scale_color_manual(values = my_color) +
  scale_x_discrete(labels=c("0","2","4","6", "8", "10")) +
  my_theme()      

#+++++++++++++++
# Plot HCT ####
#+++++++++++++++
# Repeated measures 
f_HCT <- CBC_summary %>%
  filter(pop != 'Control') %>%
  group_by(host) %>%
  friedman_test(`HCT(%)_mean` ~ pop | time) %>%
  p_round(p, 3)

# Artificial value for plotting
f_HCT$pop <- c(1, 2)
f_HCT$time <- c(2,2)

# Between hosts 
p_HCT_h1 <- CBC_refined %>%
  filter(pop != 'Control') %>%
  group_by(time) %>%
  wilcox_test(`HCT(%)` ~ host) %>%
  adjust_pvalue(method = 'BH') %>%
  add_significance(symbols = c('####', '###', '##', '#', 'ns')) %>%
  filter(p.adj <= 0.05) %>%
  mutate(host = "C57BL/6", pop = 'OR')

# Plot HCT
plot_HCT <- 
  ggplot(CBC_summary, aes(time, `HCT(%)_mean`, group = interaction(pop))) +
  geom_point(size =1.5, position = position_dodge(width = 0.5), aes(color = pop)) +
  geom_line(position = position_dodge(width = 0.5), aes(color = pop)) +
  geom_errorbar(aes(ymin=`HCT(%)_mean`-`HCT(%)_se`, ymax=`HCT(%)_mean`+`HCT(%)_se`, 
                    color = pop), width=.2,
                position=position_dodge(0.5)) +
  facet_grid(host~., scale ='fixed') +
  geom_segment(data=CBC_all[CBC_all$host=="BALB/c",], aes(x=0, y=0, xend=6.5, yend=0), linetype="dashed",
               color = "black", linewidth=0.5) +
  geom_text(data = f_HCT, aes(label = paste0('Friedman p = ', round(p, digits =10)), y =10)) +
  geom_text(data = p_HCT_h1, aes(label = p.adj.signif), vjust = -0.5, y = 20) +
  ylab("Mean HCT (%)") +
  xlab('Weeks post infection') +
  labs(color = "Population") +
  scale_color_manual(values = my_color) +
  scale_x_discrete(labels=c("0","2","4","6", "8", "10")) +
  my_theme()
 
#+++++++++++++++++++++++++
# Plot neutrophils (units) ####
#+++++++++++++++++++++++++

# Repeated measures 
f_NEUT_uni <- CBC_summary %>%
  filter(pop != 'Control') %>%
  group_by(host) %>%
  friedman_test(`NEUT#(K/uL)_mean` ~ pop | time) %>%
  p_round(p, 3)

# Artificial value for plotting
f_NEUT_uni$pop <- c(1, 2)
f_NEUT_uni$time <- c(2,2)


# Between hosts 
p_neut_h1 <- CBC_refined %>%
  filter(pop != 'Control') %>%
  group_by(time) %>%
  wilcox_test(`NEUT#(K/uL)` ~ host) %>%
  adjust_pvalue(method = 'BH') %>%
  add_significance(symbols = c('####', '###', '##', '#', 'ns')) %>%
  filter(p.adj <= 0.05) %>%
  mutate(host = "C57BL/6", pop = 'OR')

# Plot
plot_neut <- ggplot(CBC_summary, aes(time, `NEUT#(K/uL)_mean`, group = interaction(pop))) +
  geom_point(size =1.5, position = position_dodge(width = 0.5), aes(color = pop)) +
  geom_line(data=subset(CBC_summary, pop != 'Control'),position = position_dodge(width = 0.5), aes(color = pop)) +
  geom_line(data=subset(CBC_summary, pop == 'Control'),position = position_dodge(width = 0.5), aes(color = pop),
            linetype = 'dashed') +
  geom_errorbar(aes(ymin=`NEUT#(K/uL)_mean`-`NEUT#(K/uL)_se`, ymax=`NEUT#(K/uL)_mean`+`NEUT#(K/uL)_se`, 
                    color = pop), width=.2,
                position=position_dodge(0.5)) +
  facet_grid(host~., scale ='fixed') +
  geom_segment(data=CBC_all[CBC_all$host=="BALB/c",], aes(x=0, y=0, xend=6.5, yend=0), linetype="dashed",
               color = "black", size=0.5) +
  geom_text(data = f_NEUT_uni, aes(label = paste0('Friedman p = ', round(p, digits =10)), y = 0.5)) +
  geom_text(data = p_neut_h1, aes(label = p.adj.signif), vjust = -0.5, y = 4.5) +
  ylab("Mean neutrophils (K/uL)") +
  xlab('Weeks post infection') +
  labs(color = "Population") +
  scale_color_manual(values = my_color) +
  scale_x_discrete(labels=c("0","2","4","6", "8", "10")) +
  my_theme()

#++++++++++++++++++++++++
# Plot basophils (units) ####
#+++++++++++++++++++++++++

# Repeated measures 
f_BASO_uni <- CBC_summary %>%
  filter(pop != 'Control') %>%
  group_by(host) %>%
  friedman_test(`BASO#(K/uL)_mean` ~ pop | time) %>%
  p_round(p, 3)


# Artificial value for plotting
f_BASO_uni$pop <- c(1, 2)
f_BASO_uni$time <- c(2,2)

temp <- CBC_summary %>% filter(host == 'C57BL/6')

baso_b <- f_conover("BASO#(K/uL)_mean", "C57BL/6")

baso_b

stat <- CBC_summary %>%
  filter(time == 'Week 10', pop != 'Control', host == 'C57BL/6') %>%
  bind_cols( baso_b)

# Between hosts 
p_baso_h1 <- CBC_refined %>%
  filter(pop != 'Control') %>%
  group_by(time) %>%
  wilcox_test(`BASO#(K/uL)` ~ host) %>%
  adjust_pvalue(method ='BH') %>%
  add_significance(symbols = c('####', '###', '##', '#', 'ns')) %>%
  filter(p.adj <= 0.05) %>%
  mutate(host = "C57BL/6", pop = 'OR')

# Plot
plot_baso <-
  ggplot(CBC_summary, aes(time, `BASO#(K/uL)_mean`, group = interaction(pop))) +
  geom_point(size =1.5, position = position_dodge(width = 0.5), aes(color = pop)) +
  geom_line(data=subset(CBC_summary, pop != 'Control'),position = position_dodge(width = 0.5), aes(color = pop)) +
  geom_line(data=subset(CBC_summary, pop == 'Control'),position = position_dodge(width = 0.5), aes(color = pop),
            linetype = 'dashed') +
  geom_errorbar(aes(ymin=`BASO#(K/uL)_mean`-`BASO#(K/uL)_se`, ymax=`BASO#(K/uL)_mean`+`BASO#(K/uL)_se`, 
                    color = pop), width=.2,
                position=position_dodge(0.5)) +
  facet_grid(host~., scale ='fixed') +
  geom_segment(data=CBC_all[CBC_all$host=="BALB/c",], aes(x=0, y=0, xend=6.5, yend=0), linetype="dashed",
               color = "black", size=0.5) +
  geom_text(data = f_BASO_uni, aes(label = paste0('Friedman p = ', round(p, digits =10)), y = 0.03)) +
  geom_text(data = p_baso_h1, aes(label = p.adj.signif), vjust = -0.5, y = 0.025) +
  geom_text_repel(
    aes(label = letter, color = pop), show.legend = F, data = stat, position = position_dodge(width = 0.5), hjust = 0.5) +
  ylab("Mean basophils (K/uL)") +
  xlab('Weeks post infection') +
  labs(color = "Population") +
  scale_color_manual(values = my_color) +
  scale_x_discrete(labels=c("0","2","4","6", "8", "10")) +
  my_theme()

#++++++++++++++++++
# EXPORT PLOTS ####
#++++++++++++++++++
# 4 panels
pdf("Figure4.pdf", width = 10, height = 8)

((plot_eo | plot_baso) & theme(axis.title.x = element_blank())) / 
  ((plot_lymph | plot_mono)  + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')) + plot_annotation(tag_levels = c("A")) &
  theme(plot.tag = element_text(face = 'bold'))

dev.off()

# Supplemental
pdf("Additional file FigureS3.pdf", width = 10, height = 8)

((plot_neut | plot_HCT)) / 
  ((plot_ret + plot_spacer())  + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')) + plot_annotation(tag_levels = c("A")) &
  theme(plot.tag = element_text(face = 'bold'))

dev.off()

# 6 panels
#((plot_eo | plot_baso) & theme(axis.title.x = element_blank())) / ((plot_lymph | plot_mono) & theme(axis.title.x = element_blank())) /
#  ((plot_neut | plot_ret)  + plot_layout(guides = 'collect') & theme(legend.position = 'bottom')) + plot_annotation(tag_levels = c("A")) &
#  theme(plot.tag = element_text(face = 'bold'))
