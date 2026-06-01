
library(ggplot2)
library(ggpubr) 
library(tidyr)
library(forcats) 
library(glue) 


# Load the dataset


allele_data <- read_excel("/Users/marcosdacruz/Downloads/daCruz_manuscript_data.xlsx", sheet = "fig4")


groups_to_omit2 <- c("Western isothermality", "Carpathian isothermality")

allele_data <- subset(allele_data, !Climatic_Variable %in% groups_to_omit2) 

allele_data$Climatic_Variable <- factor(allele_data$Climatic_Variable, levels = c("Carpathian mean temperature","Carpathian max temperature",
                                                                                  "Carpathian precipitation", "Western mean temperature", 
                                                                                  "Western max temperature", "Western precipitation"))
                                                                         
             
                                                                   
# Change the order in which the legend appears 

allele_data <- allele_data %>%
  mutate(
    category = interaction(category), # Create combined factor
    category = fct_relevel(category,
                                   "Carpathian outlier",
                                   "Carpathian neutral", 
                                   "Western outlier",
                                   "Western neutral") # Reorder levels
  )

# Create the scatter plot
ggplot(allele_data, aes(x = Climatic_Value, y = allele_frequency, color = category)) +
geom_point(size = 0.6, alpha = 1, stat = "summary", fun = "mean") + # Points for each Refugia 
  geom_errorbar(stat = "summary",
                fun.data = "mean_se", 
                size = 0.6) + 
  geom_smooth(aes(), method = "lm", se = TRUE, linetype = "dashed", size = 0.3) + # Trend lines. 
  stat_cor(aes(label = paste("R =", ..r..), fontface = "bold.italic"),
           method = "pearson",  
           show.legend = FALSE,
           label.x.npc = 0.45,
           label.y.npc = 1, 
           data = subset(allele_data,  category == "Carpathian outlier" & Climatic_Variable == "Carpathian max temperature"),
           parse = FALSE,  
           r.digits = 1, 
           output.type = "text", 
           size = 2.8
  ) +    
  stat_cor(aes(label = paste("  p = ", ..p..), fontface = "bold.italic"),
           method = "pearson",  
           show.legend = FALSE,
           label.x.npc = 0.69,
           label.y.npc = 1, 
           data = subset(allele_data,  category == "Carpathian outlier" & Climatic_Variable == "Carpathian max temperature"),
           parse = FALSE,  
           r.digits = 1, 
           output.type = "text", 
           size = 2.8
  )  + 
  stat_cor(aes(label = paste("R =", ..r..), fontface = "bold.italic"),
           method = "pearson",  
           show.legend = FALSE,
           label.x.npc = 0.45,
           label.y.npc = 0.9, 
           data = subset(allele_data,  category == "Carpathian neutral" & Climatic_Variable == "Carpathian max temperature"),
           parse = FALSE,  
           r.digits = 1, 
           output.type = "text", 
           size = 2.8
  ) +    
  stat_cor(aes(label = paste(
    "  p = ", sprintf("%.2f", after_stat(p))), fontface = "bold.italic"),
           method = "pearson",  
           show.legend = FALSE,
           label.x.npc = 0.69,
           label.y.npc = 0.9, 
           data = subset(allele_data,  category == "Carpathian neutral" & Climatic_Variable == "Carpathian max temperature"),
           parse = FALSE,  
           r.digits = 1, 
           output.type = "text", 
           size = 2.8
  ) + 
  stat_cor(aes(label = paste("R =", ..r..), fontface = "bold.italic"),
           method = "pearson",  
           show.legend = FALSE,
           label.x.npc = 0.45,
           label.y.npc = 1, 
           data = subset(allele_data,  category == "Western outlier" & Climatic_Variable == "Western max temperature"),
           parse = FALSE,  
           r.digits = 1, 
           output.type = "text",
           size = 2.8
  ) +    
  stat_cor(aes(label = paste(
    "  p = ", sprintf("%.2f", after_stat(p))), fontface = "bold.italic"),
           method = "pearson",  
           show.legend = FALSE,
           label.x.npc = 0.70,
           label.y.npc = 1, 
           data = subset(allele_data,  category == "Western outlier" & Climatic_Variable == "Western max temperature"),
           parse = FALSE,  
           r.digits = 1, 
           output.type = "text", 
           size = 2.8
  )  + 
  stat_cor(aes(label = paste(
    "R = ", sprintf("%.2f", after_stat(r))), fontface = "bold.italic"),
           method = "pearson",  
           show.legend = FALSE,
           label.x.npc = 0.45,
           label.y.npc = 0.9, 
           data = subset(allele_data,  category == "Western neutral" & Climatic_Variable == "Western max temperature"),
           parse = FALSE,  
           r.digits = 1, 
           output.type = "text", 
           size = 2.8
  ) +    
  stat_cor(aes(label = paste(
    "  p = ", sprintf("%.2f", after_stat(p))), fontface = "bold.italic"),
    method = "pearson",  
    show.legend = FALSE,
    label.x.npc = 0.69,
    label.y.npc = 0.9, 
    data = subset(allele_data,  category == "Western neutral" & Climatic_Variable == "Western max temperature"),
    parse = FALSE,  
    r.digits = 1, 
    output.type = "text", 
    size = 2.8
  ) + 
  stat_cor(aes(label = paste(
    "R =  ", sprintf("%.2f", after_stat(r))), fontface = "bold.italic"),
           method = "pearson",  
           show.legend = FALSE,
           label.x.npc = 0.51,
           label.y.npc = 1, 
           data = subset(allele_data,  category == "Carpathian outlier" & Climatic_Variable == "Carpathian mean temperature"),
           parse = FALSE,  
           r.digits = 1, 
           output.type = "text", 
           size = 2.8
  ) +    
  stat_cor(aes(label = paste(
    "  p = ", sprintf("%.3f", after_stat(p))), fontface = "bold.italic"),
           method = "pearson",  
           show.legend = FALSE,
           label.x.npc = 0.75,
           label.y.npc = 1, 
           data = subset(allele_data,  category == "Carpathian outlier" & Climatic_Variable == "Carpathian mean temperature"),
           parse = FALSE,  
           r.digits = 1, 
           output.type = "text", 
           size = 2.8
  )  + 
  stat_cor(aes(label = paste("R =", ..r..), fontface = "bold.italic"),
           method = "pearson",  
           show.legend = FALSE,
           label.x.npc = 0.51,
           label.y.npc = 0.9, 
           data = subset(allele_data,  category == "Carpathian neutral" & Climatic_Variable == "Carpathian mean temperature"),
           parse = FALSE,  
           r.digits = 1, 
           output.type = "text", 
           size = 2.8
  ) +    
  stat_cor(aes(label = paste(
    "  p = ", sprintf("%.3f", after_stat(p))), fontface = "bold.italic"),
    method = "pearson",  
    show.legend = FALSE,
    label.x.npc = 0.75,
    label.y.npc = 0.9, 
    data = subset(allele_data,  category == "Carpathian neutral" & Climatic_Variable == "Carpathian mean temperature"),
    parse = FALSE,  
    r.digits = 1, 
    output.type = "text", 
    size = 2.8
  ) + 
  stat_cor(aes(label = paste("R =", ..r..), fontface = "bold.italic"),
           method = "pearson",  
           show.legend = FALSE,
           label.x.npc = 0.51,
           label.y.npc = 1, 
           data = subset(allele_data,  category == "Western outlier" & Climatic_Variable == "Western mean temperature"),
           parse = FALSE,  
           r.digits = 1, 
           output.type = "text",
           size = 2.8
  ) +    
  stat_cor(aes(label = paste(
    "  p = ", sprintf("%.3f", after_stat(p))), fontface = "bold.italic"),
           method = "pearson",  
           show.legend = FALSE,
           label.x.npc = 0.75,
           label.y.npc = 1, 
           data = subset(allele_data,  category == "Western outlier" & Climatic_Variable == "Western mean temperature"),
           parse = FALSE,  
           r.digits = 1, 
           output.type = "text", 
           size = 2.8
  )  + 
  stat_cor(aes(label = paste("R =", ..r..), fontface = "bold.italic"),
           method = "pearson",  
           show.legend = FALSE,
           label.x.npc = 0.51,
           label.y.npc = 0.9, 
           data = subset(allele_data,  category == "Western neutral" & Climatic_Variable == "Western mean temperature"),
           parse = FALSE,  
           r.digits = 1, 
           output.type = "text", 
           size = 2.8
  ) +    
  stat_cor(aes(label = paste(
    "  p = ", sprintf("%.3f", after_stat(p))), fontface = "bold.italic"),
    method = "pearson",  
    show.legend = FALSE,
    label.x.npc = 0.75,
    label.y.npc = 0.9, 
    data = subset(allele_data,  category == "Western neutral" & Climatic_Variable == "Western mean temperature"),
    parse = FALSE,  
    r.digits = 1, 
    output.type = "text",
    size = 2.8
  ) + 
  stat_cor(aes(label = paste(
    "R =  ", sprintf("%.3f", after_stat(r))), fontface = "bold.italic"),
    method = "pearson",  
    show.legend = FALSE,
    label.x.npc = 0.52,
    label.y.npc = 1, 
    data = subset(allele_data,  category == "Carpathian outlier" & Climatic_Variable == "Carpathian precipitation"),
    parse = FALSE,  
    r.digits = 1, 
    output.type = "text", 
    size = 2.8
  ) +    
  stat_cor(aes(label = paste(
    "  p = ", sprintf("%.2f", after_stat(p))), fontface = "bold.italic"),
           method = "pearson",  
           show.legend = FALSE,
           label.x.npc = 0.76,
           label.y.npc = 1, 
           data = subset(allele_data,  category == "Carpathian outlier" & Climatic_Variable == "Carpathian precipitation"),
           parse = FALSE,  
           r.digits = 1, 
           output.type = "text", 
           size = 2.8
  )  + 
  stat_cor(aes(label = paste(
    "R =", sprintf("%.3f", after_stat(r))), fontface = "bold.italic"),
           method = "pearson",  
           show.legend = FALSE,
           label.x.npc = 0.52,
           label.y.npc = 0.9, 
           data = subset(allele_data,  category == "Carpathian neutral" & Climatic_Variable == "Carpathian precipitation"),
           parse = FALSE,  
           r.digits = 1, 
           output.type = "text", 
           size = 2.8
  ) +    
  stat_cor(aes(label = paste(
    "  p = ", sprintf("%.2f", after_stat(p))), fontface = "bold.italic"),
    method = "pearson",  
    show.legend = FALSE,
    label.x.npc = 0.76,
    label.y.npc = 0.9, 
    data = subset(allele_data,  category == "Carpathian neutral" & Climatic_Variable == "Carpathian precipitation"),
    parse = FALSE,  
    r.digits = 1, 
    output.type = "text", 
    size = 2.8
  ) + 
  stat_cor(aes(label = paste(
    "R =  ", sprintf("%.3f", after_stat(r))), fontface = "bold.italic"),
    method = "pearson",  
    show.legend = FALSE,
    label.x.npc = 0.52,
    label.y.npc = 1, 
    data = subset(allele_data,  category == "Western outlier" & Climatic_Variable == "Western precipitation"),
    parse = FALSE,  
    r.digits = 1, 
    output.type = "text", 
    size = 2.8
  ) +    
  stat_cor(aes(label = paste(
    "  p = ", sprintf("%.2f", after_stat(p))), fontface = "bold.italic"),
    method = "pearson",  
    show.legend = FALSE,
    label.x.npc = 0.76,
    label.y.npc = 1, 
    data = subset(allele_data,  category == "Western outlier" & Climatic_Variable == "Western precipitation"),
    parse = FALSE,  
    r.digits = 1, 
    output.type = "text", 
    size = 2.8
  )  + 
  stat_cor(aes(label = paste(
    "R =", sprintf("%.3f", after_stat(r))), fontface = "bold.italic"),
           method = "pearson",  
           show.legend = FALSE,
           label.x.npc = 0.52,
           label.y.npc = 0.9, 
           data = subset(allele_data,  category == "Western neutral" & Climatic_Variable == "Western precipitation"),
           parse = FALSE,  
           r.digits = 1, 
           output.type = "text", 
           size = 2.8
  ) +    
  stat_cor(aes(label = paste(
    "  p = ", sprintf("%.2f", after_stat(p))), fontface = "bold.italic"),
    method = "pearson",  
    show.legend = FALSE,
    label.x.npc = 0.76,
    label.y.npc = 0.9, 
    data = subset(allele_data,  category == "Western neutral" & Climatic_Variable == "Western precipitation"),
    parse = FALSE,  
    r.digits = 1, 
    output.type = "text", 
    size = 2.8
  ) +
  # Correlation values for each facet + # Correlation values for each facet
  facet_wrap(~ Climatic_Variable, scales = "free_x", nrow = 2) + # Facet by climatic variable  
  theme(text = element_text(family = "Helvetica")) +
  labs(
    x = "Climatic value",
    y = "Average allele frequency",
    color = "Allele origin and type",
    linetype = "Refugia"
  ) + 
  #theme_minimal(base_size = 5) + 
  scale_color_manual(values = c("Western neutral" = "#f4b8a7", "Western outlier" = "#b85042", "Carpathian neutral" = "#92c7dd",  "Carpathian outlier" = "#1978a5")) +
  theme(  
    strip.background = element_rect(fill = "white", color = "white"), # Set background for facet labels
    panel.background = element_rect(fill = "white", color = "white"), 
    strip.text = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"), 
    axis.line = element_line(color = "black", size = 0.2), 
    legend.text = element_text(size = 11, face = "bold"), 
    legend.title = element_text(size = 11, face = "bold"),
    panel.grid.major = element_line(color = "gray", size = 0.25, linetype = "dotted"),
    panel.grid.minor = element_line(color = "gray", size = 0.25, linetype = "dotted"), 
    legend.background = element_rect(fill = "white", color = "white"), # Set white background for the legend
    legend.box.background = element_rect(fill = "white", color = "white"), # White background around the legend box
    plot.background = element_rect(fill = "white", color = "white") # Ensures the entire plot area has a white background
  )

ggsave("/Users/marcosdacruz/Dropbox/PhD_project/bank_vole_paper/figures_new/neutral_outlier_mean_standard_error_allele_freq_climate_corr_corrected_mean_v8_fixed.png", device = "png", width = 26, height = 16, units = "cm", dpi = 1200) 

ggsave("/Users/marcosdacruz/Dropbox/PhD_project/bank_vole_paper/figures_new/neutral_outlier_mean_standard_error_allele_freq_climate_corr_corrected_mean_v8_fixed.pdf", device = "pdf", width = 26, height = 16, units = "cm", dpi = 1200) 

ggsave("/Users/marcosdacruz/Dropbox/PhD_project/bank_vole_paper/figures_new/neutral_outlier_mean_standard_error_allele_freq_climate_corr_corrected_mean_v8_fixed.svg", device = "svg", width = 26, height = 16, units = "cm", dpi = 1200)

