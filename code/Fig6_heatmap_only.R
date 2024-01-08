# Figure 6 - Effect size
# Author: Kathrin Jutzeler
# Date: May 26, 2023
# Last updated: September 6, 2023

# Load packages
library(tidyverse) # for everything
library(MDMR) # To calculate effect size
library(patchwork) # to export plots
#library(PoiClaClu)
library(RColorBrewer) # To define colors
library(ggplotify) # To convert to ggplot
library(pheatmap) # To plot heatmap

# Define theme ####

my_theme <- function(){theme(panel.grid = element_blank(), legend.position = "none", strip.placement = 'outside',
                             axis.line = element_line(1),
                             panel.spacing = unit(0,'lines'), 
                             text = element_text(size = 14),
                             axis.text = element_text(size=12, color="black"),
                             legend.text = element_text(size=12, color="black"),
                             legend.title = element_text(size=13, color="black"),
                             plot.title = element_text(size=14, face='bold')) }

my_palette <- colorRampPalette(brewer.pal(7, "YlOrRd"))(255)


# Import data #### 
all_df <- read.csv('all_data.csv') 

all_df$host <- as.factor(all_df$host)
all_df$pop <- as.factor(all_df$pop)

all_df <- all_df %>% filter(pop != 'Control') %>% droplevels()

# Select data
em_full <- all_df %>%
  dplyr::select(total_worms, intestine_eggs, liver_eggs, fecundity, liver_wt_norm, spleen_wt_norm, intestine_length,
                gain, fibrosis, granuloma,
                4:8, # CBC data
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

#+++++++++++++
# Heatmap ####
#+++++++++++++

# Labels for plot
labs <- c("Total worms", "Intestine eggs", "Liver eggs", 'Fecundity', 'Liver weight', 'Spleen weight', 'Intestine length', 
          'Body weight gain', 'Fibrotic area', 'Granuloma area', 
          'Wk 10 EO', 'Wk 10 NEUT', 'Wk 10 MONO', 'Wk 10 LYMPH', 'Wk 10 BASO',
          'IFNy', 'TNFa', 'IL2', 'IL4', 'IL5', 'IL6', 'IL10', 'IL13')
          
# Calculate the model
model <- aov(em_full ~ all_df$pop * all_df$host)

aovall <- summary.aov(model)

# Get p values
pvalues <- 
  lapply(aovall, function(x) x[1:3,5])

pvalues <- bind_rows(pvalues)
pvalues <- gather(pvalues, 'response', 'pvalue')

# Calcualate effect size
output <- 
  effectsize::eta_squared(model, partial = T)

output <- bind_cols(output, pvalues)

output <- output %>%
  mutate(significant = ifelse(pvalue <= 0.05 & pvalue > 0.01, paste0(round(Eta2_partial,3),"*"),
                              ifelse(pvalue <= 0.01 & pvalue > 0.001, paste0(round(Eta2_partial,3), "**"),
                                     ifelse(pvalue <= 0.001, paste0(round(Eta2_partial, 3),"***"), 
                                            round(Eta2_partial, 3)))))
#mtest <- car::Anova(model)
#effectsize::eta_squared(mtest)

# Transform output and convert to matrix
output <- as.data.frame(output)

# This is for the labels
output_signif <- output %>% 
  dplyr::select(1:2, 9) %>%
  pivot_wider(
    names_from = Response, 
    values_from = significant)

mat_label <- as.matrix(output_signif)
rownames(mat_label) <- mat_label[,1]
mat_label <- mat_label[,-1]

# This is for the acutal heatmap
output_eta <- output %>% 
  dplyr::select(1:3) %>%
  pivot_wider(
    names_from = Response, 
    values_from = Eta2_partial)

mat <- as.matrix(output_eta)
rownames(mat) <- mat[,1]
mat <- mat[,-1]

mat %>% data.frame() %>% mutate(across(where(is.character), as.numeric)) %>% as.matrix() -> mat.num
mat.num <- round(mat.num,3)

# Plot 
p_bypheno_int2 <- as.ggplot(pheatmap(t(mat.num), color = my_palette, #scale = 'column' , 
         cluster_rows = F, cluster_cols = F, 
   labels_col =  c("Parasite", "Host", "Interaction"),  
  display_numbers = t(mat_label), labels_row = labs,
 number_color = 'black', fontsize = 16, main = "Effect size", legend = F, cellwidth = 50,
cellheight = 30, angle_col = 45))

#+++++++++++++++++++
# Export plots  ####
#+++++++++++++++++++
# As jpg
#ggsave('Figure6.jpg', p, width =10, height = 12)

fig_6_int <- as.ggplot(((p_bypheno_int) ) + plot_layout(widths = c(4, 0.5, 4)) +
                         plot_annotation(tag_levels = c("A")) &  theme(plot.tag = element_text(face = 'bold', size = 16)))

ggsave('Figure6_int.jpg', fig_6_int, width =10, height = 12)

png('Fig6.png', width = 5, height = 11, unit = 'in', res =300)
p_bypheno_int2
dev.off()
