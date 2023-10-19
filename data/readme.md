# Datasets for the generation of figures/models for the manuscript

## all_data.csv
This file includes data that was used for trait-by-trait and comprehensive analyses.    
  
**Description of columns**  
sample: unique identifier of mouse. Number indicates mouse in consecutive order of the infection and letter indicates parasite line.   
host: indicates mouse strain  
pop:  indicates control or parasite population used for the infection  
Week.10_EO: eosinophil counts in K/ul at week 10 post-infection	Week  
Week.10_NEUT: neutrophil counts in K/ul at week 10 post-infection  
Week.10_MONO: monocyte counts in K/ul at week 10 post-infection  
Week.10_BASO: basophil counts in K/ul at week 10 post-infection  
columns IFNy, IL5, TNFa, IL2, IL6, IL4, IL10, IL13: respective cytokine levels in pg/ml  
body_wt: body weight in gram on day of euthanasia, 12 weeks post-infection  
female_worms: number of female worms counted after the perfusion  
male_worms: number of male worms counted after the perfusion  
total_worms: sum of female_worms and male_worms  
liver_wt: liver weight measured in grams
liver_wt_norm: liver weight as percentage of total body weight  
spleen_wt: spleen weight measured in grams  
spleen_wt_norm: spleen weight as percentage of total body weight  
intestine_wt: weight of the intestine (to calculate egg burden in the intestine)  
liver_eggs: number of eggs per gram of liver  
intestine_eggs: number of eggs per gram of intestine  
intestine_length: length of intestine measured in mm  
total_eggs: sum of liver_eggs and intestine_eggs  
fecundity: number of eggs divided by number of female_worms  
fibrosis: area of blue collagen stain as percentage of total stained area  
granuloma: average granuloma area per mouse 

## CBC_data.csv
This file includes data that was used for Figure 4 - Complete blood count with differential. The output is as generated by the IDEXX ProCyte hematology analyzer.  

## body_weight.csv  
This file includes all body weight measurements and was used to generate figure 3A.  

**Description of columns**  
animal_ID: unique identifier of mouse. Number indicates mouse in consecutive order of the infection and letter indicates parasite line.   
host: indicates mouse strain  
pop:  indicates control or parasite population used for the infection  
time: week number post-infection when weight was measured  
weight: mouse weight in grams  
gain: calculated weight gain in percent    

## infection_data.csv  
This file includes all body weight measurements and was used to generate figure 3A.  

**Description of columns**  
animal_ID: unique identifier of mouse. Number indicates mouse in consecutive order of the infection and letter indicates parasite line.   
host: indicates mouse strain  
pop:  indicates control or parasite population used for the infection  
no_cercariae: number of cercariae used for the infection  
no_cercariae_remain: number of cercariae counted in the water after the infection    
no_heads_remain: number of cercariae heads in the water after the infection  
days_infected: number of days the mouse was infected for  
status: mouse lived (0) or died (1) before the infection period was over  
cercariae	penrate: calculated infection rate