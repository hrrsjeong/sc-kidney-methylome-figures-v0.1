library(tidyverse)
library(ggrepel)

INPUT_CSV <- "all_genes_comprehensive_summary.csv" 
OUTPUT_PLOT <- "age_vs_deg_scatter.pdf"

# DEG Significance Threshold
PADJ_CUTOFF <- 0.01

main <- function() {
  message("--- Age vs Disease Correlation Analysis (R) ---")
  
  if (!file.exists(INPUT_CSV)) {
    stop(paste("Error:", INPUT_CSV, "not found."))
  }
  
  df <- read_csv(INPUT_CSV, show_col_types = FALSE)
  
  # 1. Filter Data
  # Must have valid Age and Disease stats
  plot_df <- df %>%
    filter(!is.na(Spearman_Rho), !is.na(log2FoldChange), !is.na(DEG_padj))
  
  message(paste("Plotting", nrow(plot_df), "genes..."))
  
  # 2. Process Data
  # Invert LFC: Original is Healthy > Altered. We want Altered > Healthy.
  plot_df <- plot_df %>%
    mutate(
      LFC_Altered = -log2FoldChange,
      
      Category = case_when(
        DEG_padj < PADJ_CUTOFF & LFC_Altered > 0 ~ "DEG Up (Altered)",
        DEG_padj < PADJ_CUTOFF & LFC_Altered < 0 ~ "DEG Down (Altered)",
        TRUE ~ "Not Significant"
      ),
      
      # Calculate distance from origin for labeling top hits
      dist = sqrt(Spearman_Rho^2 + LFC_Altered^2)
    )
  
  # 3. Calculate Correlation
  cor_res <- cor.test(plot_df$Spearman_Rho, plot_df$LFC_Altered, method = "pearson")
  r_val <- cor_res$estimate
  p_val <- cor_res$p.value
  
  message(paste0("Correlation: R=", sprintf("%.3f", r_val), ", P=", sprintf("%.2e", p_val)))
  
  # 4. Select Top Genes for Labeling (Top 10 most extreme)
  top_hits <- plot_df %>%
    arrange(desc(dist)) %>%
    slice_head(n = 10)
  
  # 5. Plotting
  # Set Factor levels for legend order
  plot_df$Category <- factor(plot_df$Category, 
                             levels = c("DEG Up (Altered)", "DEG Down (Altered)", "Not Significant"))
  
  # Colors matching the Python script
  # Red (#d73027) for Up, Blue (#4575b4) for Down, Grey for NS
  my_colors <- c("DEG Up (Altered)" = "#d73027", 
                 "DEG Down (Altered)" = "#4575b4", 
                 "Not Significant" = "lightgrey")
  
  # Alpha (transparency) mapping
  my_alpha <- c("DEG Up (Altered)" = 0.6, 
                "DEG Down (Altered)" = 0.6, 
                "Not Significant" = 0.3)
  
  p <- ggplot(plot_df, aes(x = Spearman_Rho, y = LFC_Altered)) +
    # Add Regression Line
    
    # Scatter Points
    # Note: We map alpha manually or split layers, but ggplot handles alpha mapping too.
    # For simplicity, we can just plot NS first then Sig on top by arranging.
    geom_point(data = plot_df %>% arrange(desc(Category == "Not Significant")),
               aes(color = Category, fill = Category),
               shape = 21, color = "transparent", size = 2, alpha = 0.5) +
    
    scale_fill_manual(values = my_colors) +
    geom_smooth(method = "lm", color = "black", size = 0.5) +
    
    # Quadrant Lines
    geom_hline(yintercept = 0, color = "black", size = 0.3,linetype = "dashed") +
    geom_vline(xintercept = 0, color = "black", size = 0.3,linetype = "dashed") +
    xlim(c(-0.5,0.6))+
    ylim(c(-3,5))+
    # Labels with ggrepel
    #geom_text_repel(data = top_hits, aes(label = gene_name),
    ##                size = 3, fontface = "bold", color = "black",
    #                box.padding = 0.5, max.overlaps = Inf) +
    
    # Annotations (Quadrants)
    #annotate("text", x = 0.95, y = max(plot_df$LFC_Altered), label = "Discordant\n(Up/Down)", 
    #         hjust = 1, vjust = 1, color = "green", fontface = "bold") +
    #annotate("text", x = -0.95, y = max(plot_df$LFC_Altered), label = "Discordant\n(Down/Up)", 
    #         hjust = 0, vjust = 1, color = "green", fontface = "bold") +
    #annotate("text", x = 0.95, y = min(plot_df$LFC_Altered), label = "Concordant\n(Up/Up)", 
    #         hjust = 1, vjust = 0, color = "purple", fontface = "bold") +
    #annotate("text", x = -0.95, y = min(plot_df$LFC_Altered), label = "Concordant\n(Down/Down)", 
    #         hjust = 0, vjust = 0, color = "purple", fontface = "bold") +

    # Theme and Labels
    theme_classic()+
    #scale_fill_manual(values=c("black","firebrick2"))+
    #scale_color_manual(values=c("black","firebrick2"))+ 
    #scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")+
    #scale_color_aaas()+
    #scale_fill_aaas()+
    
    theme_bw() + theme(
      plot.title = element_text(size=16),
      #legend.position='in',
      legend.position='right', 
      #legend.position = c(0.80, 0.6),
      #legend.justification='left',
      #legend.direction='vertical',
      plot.background = element_blank(),
      axis.text.x = element_text(color = "black",size=12),
      #axis.text.x=element_text(color = "black", size = 16),
      axis.text.y=element_text(color = "black", size = 12),
      strip.text = element_text(size = 12),
      text = element_text(size = 12),
      axis.title.x = element_text(size = 12),
      #axis.title.x=element_blank(),
      axis.title.y = element_text(size = 12),
      legend.title=element_text(size=12),
      legend.text=element_text(size=12),
      #legend.box.background = element_rect(colour = "black"),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    labs(
      #title = "Convergence of Aging and Disease in PT Cells",
      subtitle = paste0("R=", sprintf("%.3f", r_val)),
      x = "Age Association (Spearman Rho)",
      y = "Disease Association (log2FC Altered/Healthy)",
      fill = "DEG Status"
    )
  p
  # Save
  ggsave(OUTPUT_PLOT, plot = p, width = 10, height = 8, dpi = 300, bg = "white")
  message(paste("Plot saved to:", OUTPUT_PLOT))
}

# Run
main()
