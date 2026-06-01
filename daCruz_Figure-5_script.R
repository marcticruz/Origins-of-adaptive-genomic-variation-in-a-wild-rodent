# Load necessary libraries 


library(ggplot2)
library(ggridges)  
library(dplyr) 
library(svglite)


#with shades of gray in the background   

# Load the dataset
data <- read_excel("/Users/marcosdacruz/Downloads/daCruz_manuscript_data.xlsx", sheet = "fig5")



#data$category

  # Specify the desired order for the categories
  desired_order <- c( "Western precipitation", "Carpathian precipitation", 
                      "Western max temperature", "Carpathian max temperature", 
                      "Western mean temperature", "Carpathian mean temperature")
  
  # Convert the category column to a factor with the specified order
  data$category <- factor(data$category, levels = desired_order)
  
  # Assign shades of gray for each climatic variable group

  
  # Ridgeline Plot
  ggplot(data, aes(x = correlation, y = category, fill = refugia)) +
    # Add shaded rectangles in the background for each climatic variable group  
    #geom_rect(
      #data = background_data,
      #aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
      #inherit.aes = FALSE,
      #fill = background_data$shade,
      #alpha = 0.3
    #)
    # Add ridgeline plot with increased spacing
    geom_density_ridges(alpha = 0.7, scale = 0.9, rel_min_height = 0.01) +
    # Specify colors for the ridgeline distributions
    scale_fill_manual(
      values = c(
        "Carpathian neutral" = "#c2e4f5",     # Navy
        "Carpathian outlier" = "#1978a5", # Powder Blue
        "Western neutral" = "#f4b8a7",       # Dark Red
        "Western outlier" = "#b85042"     # Light Salmon
      )
    ) +    
    scale_y_discrete(expand = c(0, 0)) +  
    expand_limits(y = 7) +
    theme_minimal() + 
    theme(text = element_text(family = "Helvetica")) +
    labs(
      x = "Correlation of allele frequency with climatic variable",
      y = "Density",
      fill = "Allele origin and type"
    ) +
    theme(
      plot.title = element_text(size = 8, face = "bold", hjust = 0.5), # Larger title
      axis.title = element_text(size = 8, face = "bold"),                            # Larger axis titles
      axis.text = element_text(size = 8, face = "bold"),                             # Larger axis text
      legend.title = element_text(size = 8, face = "bold"),                          # Larger legend title
      legend.text = element_text(size = 8, face = "bold"),                           # Larger legend text
      legend.position = "bottom",  
      legend.box = "horizontal", 
      legend.title.position = "top",
      panel.grid.major = element_line(color = "black", size = 0.2, linetype = "dotted"),    # Darker major gridlines
      panel.grid.minor = element_line(color = "black", size = 0.2, linetype = "dotted"),     # Darker minor gridlines 
      panel.border = element_blank(), # Removes the full border
      axis.line = element_line(color = "black", size = 0.8),  
      plot.background = element_rect(fill = "white", color = "white"), # Adds axis lines  
      axis.text.y = element_blank()
      
    ) +    
    guides(fill = guide_legend(ncol = 2)) +
    geom_hline(yintercept = 3, linetype = "solid", color = 	"black", size = 0.8) + 
    geom_hline(yintercept = 5, linetype = "solid", color = 	"black", size = 0.8) +
    annotate(geom = "text", x = 1.27, y = 6.9, label = "p = 0.00", color = "black", size = 2, fontface = "bold.italic") + 
    annotate(geom = "text", x = 1.27, y = 5.9, label = "p = 0.053", color = "black", size = 2, fontface = "italic") +  
    annotate(geom = "text", x = 1.27, y = 4.9, label = "p = 0.18", color = "black", size = 2, fontface = "italic") + 
    annotate(geom = "text", x = 1.27, y = 3.9, label = "p = 0.24", color = "black", size = 2, fontface = "italic") +
    annotate(geom = "text", x = 1.27, y = 2.9, label = "p = 0.03", color = "black", size = 2, fontface = "bold.italic") +  
    annotate(geom = "text", x = 1.27, y = 1.9, label = "p = 0.001", color = "black", size = 2, fontface = "bold.italic") + 
    annotate(geom = "text", x = -0.93, y = 6.8, label = "Mean temperature", color = "black", size = 4) + 
    annotate(geom = "text", x = -0.97, y = 4.8, label = "Max temperature", color = "black", size = 4) + 
    annotate(geom = "text", x = -1.11, y = 2.8, label = "Precipitation", color = "black", size = 4)
    
  
ggsave("/Users/marcosdacruz/Dropbox/PhD_project/bank_vole_paper/figures_new/ridgeline_plot_null_outlier_correlations_with_climate_with_gridlines_v10_withpvalues_fixed_test.pdf", device = "pdf", width = 10, height = 15, units = "cm", dpi = 1200)

ggsave("/Users/marcosdacruz/Dropbox/PhD_project/bank_vole_paper/figures_new/ridgeline_plot_null_outlier_correlations_with_climate_with_gridlines_v10_withpvalues_fixed.png", device = "png", width = 10, height = 15, units = "cm", dpi = 1200) 

ggsave("/Users/marcosdacruz/Dropbox/PhD_project/bank_vole_paper/figures_new/ridgeline_plot_null_outlier_correlations_with_climate_with_gridlines_v10_withpvalues_fixed.svg", device = "svg", width = 10, height = 15, units = "cm", dpi = 1200)
