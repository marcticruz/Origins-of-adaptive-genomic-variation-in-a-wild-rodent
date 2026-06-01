 

#install.packages("readxl")

library("readxl")
library("ggplot2")
library("RColorBrewer")
library("ggtext")
library("ggpubr")

# Import data  


runtc_dat <- read_excel("/Users/marcosdacruz/Downloads/daCruz_manuscript_data.xlsx", sheet = "fig3")


runtc_dat$Locality <- factor(runtc_dat$Locality, levels = c("Britain outlier", "Britain neutral", "Western outlier", "Western neutral", "Carpathian outlier", 
                                                            "Carpathian neutral", "Widespread outlier", "Widespread neutral"))



# Add column for colors  


runtc_dat$coloring <- vector("character", length = nrow(runtc_dat)) 


runtc_dat[runtc_dat$Locality == "Britain outlier",]$coloring <- "#a7beae"   


runtc_dat[runtc_dat$Locality == "Britain neutral",]$coloring <- "#dbe8df"   

runtc_dat[runtc_dat$Locality == "Western outlier",]$coloring <- "#b85042"     

runtc_dat[runtc_dat$Locality == "Western neutral",]$coloring <- "#f4b8a7" 


runtc_dat[runtc_dat$Locality == "Carpathian outlier",]$coloring <- "#1978a5"   

runtc_dat[runtc_dat$Locality == "Carpathian neutral",]$coloring <- "#c2e4f5"


runtc_dat[runtc_dat$Locality == "Widespread outlier",]$coloring <-  "#967bb6"     

runtc_dat[runtc_dat$Locality == "Widespread neutral",]$coloring <-  "#d2c7e3" 

# Add shapes column 

runtc_dat$shapes <- vector("numeric", length = nrow(runtc_dat)) 


runtc_dat[runtc_dat$Locality == "Britain outlier",]$shapes <- 21 

runtc_dat[runtc_dat$Locality == "Britain neutral",]$shapes <- 21 

runtc_dat[runtc_dat$Locality == "Western outlier",]$shapes <- 24 

runtc_dat[runtc_dat$Locality == "Western neutral",]$shapes <- 24  

runtc_dat[runtc_dat$Locality == "Carpathian outlier",]$shapes <- 22  

runtc_dat[runtc_dat$Locality == "Carpathian neutral",]$shapes <- 22 


runtc_dat[runtc_dat$Locality == "Widespread outlier",]$shapes <-  23     

runtc_dat[runtc_dat$Locality == "Widespread neutral",]$shapes <-  23 

# Add x-adjusted position column 

runtc_dat$x_adjusted <- vector("numeric", length = nrow(runtc_dat)) 


runtc_dat[runtc_dat$Locality == "Britain outlier",]$x_adjusted <- 1 

runtc_dat[runtc_dat$Locality == "Britain neutral",]$x_adjusted <- 1.69 

runtc_dat[runtc_dat$Locality == "Western outlier",]$x_adjusted <- 3 

runtc_dat[runtc_dat$Locality == "Western neutral",]$x_adjusted <- 3.69 

runtc_dat[runtc_dat$Locality == "Carpathian outlier",]$x_adjusted <- 5  

runtc_dat[runtc_dat$Locality == "Carpathian neutral",]$x_adjusted <- 5.69 


runtc_dat[runtc_dat$Locality == "Widespread outlier",]$x_adjusted <- 7     

runtc_dat[runtc_dat$Locality == "Widespread neutral",]$x_adjusted <- 7.69     



custom_labels <- c("outlier", "outlier", " ", "neutral", "neutral")     

sublabels <- c(" ", " ", "Britain", " ", " ",  " ", " ", " ", "Western", " ", " ", " ", " ", " ", "Carpathian", " ", " ",
               " ", " ", " ", "Widespread", " ", " ")



options(scipen = 999)

g1 <- ggplot(runtc_dat, aes(y = Age, x = x_adjusted, fill = Locality, colour = Locality)) + 
  stat_boxplot(geom = "errorbar", color = "black", lwd = 0.2, width = 0.3) + 
  geom_jitter(colour = "black", shape = runtc_dat$shapes, size = 0.7, alpha = 1, width = 0.3, stroke = 0) +
  geom_boxplot(colour = "black", outlier.shape=NA,lwd = 0.2, fatten = 1.5, alpha=0.7, width = 0.6) + 
  scale_x_continuous(
    breaks = c(1.0, 1.35, 1.69, 3, 3.35, 3.69, 5, 5.35, 5.69, 7, 7.35, 7.69), # Center positions for each group
    labels = function(x) {
      # Combine custom labels and sublabels into multi-line tick labels
      paste0(custom_labels[seq(1, length(custom_labels), 2)], "<br>", 
             "<b>", sublabels[seq(1, length(sublabels), 2)], "</b>")
    }
  ) + 
  scale_y_log10() +   
  #annotate("text", x = 0.55, y = 10000000, size = 3, label = "(b)") +
  theme(panel.background = element_blank(), 
        plot.title = element_text(face = "bold", size = 12),
        strip.text = element_text(size = 12, face = "bold"),
        axis.line = element_line(color = "black", linewidth = 0.5), 
        axis.title.x = element_text(size = 12, face = "bold"), 
        axis.title.y = element_text(size = 12, face = "bold"), 
        axis.text.x = element_markdown(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, face = "bold"), 
        legend.text = element_text(size = 12, face = "bold"),  
        #plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm") 
  ) +
  scale_fill_manual(values = c("#a7beae", "#dbe8df", "#b85042", "#f4b8a7" , "#1978a5", "#c2e4f5",  "#967bb6", "#d2c7e3")) +  
  scale_color_manual(values = c("#a7beae", "#dbe8df", "#b85042", "#f4b8a7" , "#1978a5", "#c2e4f5", "#967bb6", "#d2c7e3")) +
  labs(title = "(A)", x = "Allele origin and type", y = 'Allele age (Years before present)', size = 5)  +   
  theme(text = element_text(family = "Helvetica")) +
  geom_hline(yintercept = 8000, linetype = "dashed", color = 	"#FF7F50", size = 0.2) + 
  geom_hline(yintercept = 14700, linetype = "dashed", color = 	"black", size = 0.2) +
  guides(fill="none") + 
  annotate("text", x = c(1.5, 3.5, 5.5, 7.5), y = -10, label = c("Britain", "Western", "Carpathian", "Widespread"),
           size = 3, fontface = "bold", hjust = 0)

#Importing data for heatmap

data <-  read_excel("/Users/marcosdacruz/Downloads/daCruz_manuscript_data.xlsx", sheet = "fig3_heatmap")

# Check if the dataset has at least 3 columns
if (ncol(data) < 3) {
  stop("Your dataset must have at least three columns for x, y, and value.")
}

# Automatically assign the first three columns
x_col <- names(data)[1]
y_col <- names(data)[2]
value_col <- names(data)[3]

message("Mapping columns: x = ", x_col, ", y = ", y_col, ", value = ", value_col)

# Convert x and y to factors
data[[x_col]] <- factor(data[[x_col]], levels = unique(data[[x_col]]))
data[[y_col]] <- factor(data[[y_col]], levels = unique(data[[y_col]]))

# --- NEW STEP: Create a category for coloring ---
# This creates a TRUE/FALSE column: TRUE if significant, FALSE if not
data$is_significant <- data[[value_col]] <= 0.05

# Define the scale
# We use scale_fill_manual to assign exact colors to TRUE and FALSE
if (is.numeric(data[[value_col]])) {
  fill_scale <- scale_fill_manual(
    values = c("FALSE" = "white", "TRUE" = "#08306b"), # White for >0.05, Dark Blue for <=0.05
    labels = c("FALSE" = "> 0.05", "TRUE" = "\u2264 0.05"), # Custom labels for legend
    name = "p-value"
  )
  title_addition <- ""
} else {
  fill_scale <- scale_fill_identity()
  title_addition <- "(Using Identity Scale)"
}

# Create the heatmap
heatmap_plot <- ggplot(data, aes_string(x = x_col, y = y_col, fill = "is_significant")) +
  geom_tile(color = "white") +             
  geom_text(aes(label = round(P.adj, 3)),
            # Text color logic: White text on blue tiles, Black text on white tiles
            color = I(ifelse(data$is_significant, "white", "black")),
            size = 2.5, fontface = "bold") +
  fill_scale +                             
  labs(title = paste("(B)", title_addition),
       x = " ", 
       y = "Allele origin and type"
  ) + 
  theme(text = element_text(family = "Helvetica")) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12), 
        plot.title = element_text(face = "bold", size = 12),
        axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12), 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        # Optional: Add a border to the legend keys so the white one is visible
        legend.key = element_rect(color = "black", size = 0.2)) 




ggarrange(g1, heatmap_plot, nrow = 2, align = "hv") 


ggsave("/Users/marcosdacruz/Dropbox/PhD_project/bank_vole_paper/figures_new/logscale_neutral_outlier_boxjitter_heatmap_runtc_results_fixed_v12.png", device = "png", width = 27, height = 19.5, units = "cm", dpi = 1200) 

ggsave("/Users/marcosdacruz/Dropbox/PhD_project/bank_vole_paper/figures_new/logscale_neutral_outlier_boxjitter_heatmap_runtc_results_fixed_v12.pdf", device = "pdf", width = 27, height = 19.5, units = "cm", dpi = 1200)

ggsave("/Users/marcosdacruz/Dropbox/PhD_project/bank_vole_paper/figures_new/logscale_neutral_outlier_boxjitter_heatmap_runtc_results_fixed_v12.svg", device = "svg", width = 27, height = 19.5, units = "cm", dpi = 1200)



