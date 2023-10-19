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
  dplyr::select(intestine_eggs, liver_eggs, total_worms, liver_wt_norm, spleen_wt_norm, intestine_length,
                fecundity, fibrosis, granuloma,
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
labs <- c("Intestine eggs", "Liver eggs", 'Total worms', 'Liver weight', 'Spleen weight', 'Intestine length', 
          'Fecundity', 'Fibrotic area', 'Granuloma area', 
          'Wk 10 EO', 'Wk 10 NEUT', 'Wk 10 MONO', 'Wk 10 LYMPH', 'Wk 10 BASO',
          'IFNy', 'TNFa', 'IL2', 'IL4', 'IL5', 'IL6', 'IL10', 'IL13')
          
# Calculate the model
model <- aov(em_full ~ all_df$pop * all_df$host)

aovall <- summary.aov(model)

# Calcualate effect size
output <- 
  effectsize::eta_squared(model, partial = T)

#mtest <- car::Anova(model)
#effectsize::eta_squared(mtest)

# Transform output and convert to matrix
output <- as.data.frame(output)

output <- output %>% 
  dplyr::select(1:3) %>%
  pivot_wider(
    names_from = Response, 
    values_from = Eta2_partial)

mat <- as.matrix(output)
rownames(mat) <- mat[,1]
mat <- mat[,-1]

mat %>% data.frame() %>% mutate(across(where(is.character), as.numeric)) %>% as.matrix() -> mat.num
mat.num <- round(mat.num,3)

# Plot 
p_bypheno_int2 <- as.ggplot(pheatmap(t(mat.num), color = my_palette, #scale = 'column' , 
         cluster_rows = F, cluster_cols = F, 
   labels_col =  c("Parasite", "Host", "Interaction"),  
  display_numbers = t(mat.num),labels_row = labs,
 number_color = 'black', fontsize = 12, main = "Effect size", legend = F, cellwidth = 50,
cellheight = 30, angle_col = 45))

#+++++++++++++++++++
# Export plots  ####
#+++++++++++++++++++
# As jpg
#ggsave('Figure6.jpg', p, width =10, height = 12)

fig_6_int <- as.ggplot(((p_bypheno_int) ) + plot_layout(widths = c(4, 0.5, 4)) +
                         plot_annotation(tag_levels = c("A")) &  theme(plot.tag = element_text(face = 'bold', size = 16)))

ggsave('Figure6_int.jpg', fig_6_int, width =10, height = 12)

png('Fig6.png', width = 4, height = 11, unit = 'in', res =300)
p_bypheno_int
dev.off()
