Plot alpha diversity indices in a boxplot
================

After checking the rarefaction and species richness, we can also check other alpha diversity measurements. Here, we use the Simpson index, which we calculate with the `vegan` package and visualise the diversity in a boxplot.

Load Data
---------

``` r
rm(list = ls())

library(ggplot2)
library(vegan)
library(RColorBrewer)
library(ggpubr)
library(viridis)
library(agricolae)

#setwd("04_Alpha_Diversity/")

OTU_Table = as.data.frame(read.csv("../00_Data/Oomycota/05_OwnSamples_OTU_Table_min-freq-9588_transposed_withMetadata.tsv", 
                     header = T, 
                     sep = "\t", 
                     stringsAsFactors = T))

SampleMetadata = OTU_Table[,1:5]
Microhabitat = SampleMetadata$Microhabitat
TreeSpecies = SampleMetadata$TreeSpecies
OTU_Table = OTU_Table[,6:ncol(OTU_Table)]
```

Calculate Alpha Diversity
-------------------------

To run the diversity analyses, simply load the table and specify the `index` - in this case: The Simpson index. Then convert it into a dataframe and add the metadata and group.

``` r
simpson = diversity(OTU_Table, index = "simpson")
simpson = as.data.frame(simpson)
simpson$Microhabitat = Microhabitat
rownames(simpson) = SampleMetadata$SampleID
simpson$Group = "Oomycota"
#simpson$TreeSpecies = TreeSpecies

# perform Tukey's Honest Differences Test (HSD) 
# this defines significant groups which we can plot on top of the boxplots

TukeyLetters = HSD.test(aov(simpson$simpson ~ simpson$Microhabitat),
                        "simpson$Microhabitat")$group
TukeyLetters$Microhabitat = rownames(TukeyLetters)
```

Plot the Figure
---------------

Now we put the diverity measurements into a habitat specific context. It can be easiest visualised in a boxplot:

``` r
g = ggplot(simpson, aes(x = Microhabitat, y = simpson, fill = Microhabitat)) + 
  stat_boxplot(geom = "errorbar", width = 0.2, show.legend = F, lwd = 0.25) +
  geom_boxplot(show.legend = F, lwd = 0.25) + 
  geom_text(data = TukeyLetters, aes(label = groups), y = 1.025) +
  scale_fill_manual(values = c(viridis(7, direction = -1), 
                               "#8e8878", "#524640"), 
                    limits = c("Arboreal Soil", "Bark", "Deadwood", 
                                "Fresh Leaves", "Hypnum", "Lichen", 
                                "Orthotrichum", "Leaf Litter", "Soil")) + 
  scale_x_discrete(limits = c("Arboreal Soil", "Bark", "Deadwood", 
                                "Fresh Leaves", "Hypnum", "Lichen", 
                                "Orthotrichum", "Leaf Litter", "Soil")) + 
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() + 
  labs(#title = "Alpha Diversity of Samples", 
       #subtitle = "Oomycota"
       y = "Simpson Index (1-D)", 
       x = "Microhabitat") +
  #geom_dotplot(aes(x = Microhabitat, y = simpson, fill = TreeSpecies), 
  #             binaxis = "y", stackdir = "center", binwidth = 1, 
  #             binpositions = "all", dotsize = 0.01, 
  #             position = "stack") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14, face = "bold"), 
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
        plot.subtitle = element_text(size = 14, hjust = 0.5), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"))

g
```

![](AlphaBoxplot_files/figure-markdown_github/OomycotaAlphaBoxPlot-1.png)

All samples show a very high diversity. Only some of the leaf litter samples have a quite low alpha diversity. This could maybe be due to the lifestyle of some oomycetes: Saprotrophic oomycetes could experience a strong competition in this habitat.

Cercozoa
--------

Let's check if the Cercozoa show a similar pattern:

``` r
Cerco_OTU_Table = as.data.frame(read.csv("../00_Data/Cercozoa/05_Cercozoa_OwnSamples_OTU_Table_min-freq-15684_transposed_withMetadata.tsv", 
                     header = T, 
                     sep = "\t", 
                     stringsAsFactors = T))

Cerco_SampleMetadata = Cerco_OTU_Table[,1:5]
Cerco_Microhabitat = Cerco_SampleMetadata$Microhabitat
Cerco_OTU_Table = Cerco_OTU_Table[,6:ncol(Cerco_OTU_Table)]

Cerco_simpson = diversity(Cerco_OTU_Table, index = "simpson")
Cerco_simpson = as.data.frame(Cerco_simpson)
rownames(Cerco_simpson) = Cerco_SampleMetadata$SampleID
colnames(Cerco_simpson) = "simpson"
Cerco_simpson$Microhabitat = Cerco_Microhabitat
Cerco_simpson$Group = "Cercozoa"

TukeyLetters_cerco = HSD.test(aov(Cerco_simpson$simpson ~ 
                                  Cerco_simpson$Microhabitat),
                        "Cerco_simpson$Microhabitat")$group
TukeyLetters_cerco$Microhabitat = rownames(TukeyLetters_cerco)
```

Plot Cercozoa Figure
--------------------

``` r
g_cerco = ggplot(Cerco_simpson, aes(x = Microhabitat, y = simpson, fill = Microhabitat)) + 
  stat_boxplot(geom = "errorbar", width = 0.2, show.legend = F, lwd = 0.25) +
  geom_boxplot(show.legend = F, lwd = 0.25) + 
  geom_text(data = TukeyLetters_cerco, aes(label = groups), y = 1.025) +
  scale_fill_manual(values = c(viridis(7, direction = -1), 
                               "#8e8878", "#524640"), 
                    limits = c("Arboreal Soil", "Bark", "Deadwood", 
                                "Fresh Leaves", "Hypnum", "Lichen", 
                                "Orthotrichum", "Leaf Litter", "Soil")) + 
  scale_x_discrete(limits = c("Arboreal Soil", "Bark", "Deadwood", 
                                "Fresh Leaves", "Hypnum", "Lichen", 
                                "Orthotrichum", "Leaf Litter", "Soil")) + 
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() + 
  labs(#title = "Alpha Diversity of Samples",
       #subtitle = "Cercozoa",
       y = "Simpson Index (1-D)", 
       x = "Microhabitat") +
  theme(axis.text=element_text(size=12, face = "bold"), 
        axis.title=element_text(size=14, face = "bold"), 
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
        plot.subtitle = element_text(size = 14, hjust = 0.5), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"))

g_cerco
```

![](AlphaBoxplot_files/figure-markdown_github/CercoAlphaBoxPlot-1.png)

For the cercozoa, all samples show an even higher diversity. That's interesting. We could also plot both groups next to each other for a direct comparison, like this:

``` r
simpson_both = rbind(simpson, Cerco_simpson)

g_both = ggplot(simpson_both, aes(x = Microhabitat, y = simpson, fill = Microhabitat, color = Group)) + 
  stat_boxplot(geom = "errorbar", width = 0.5) +
  geom_boxplot(width = 0.5) + 
  scale_fill_manual(values = c(viridis(7, direction = -1), 
                               "#8e8878", "#524640"), 
                    limits = c("Arboreal Soil", "Bark", "Deadwood", 
                                "Fresh Leaves", "Hypnum", "Lichen", 
                                "Orthotrichum", "Leaf Litter", "Soil")) + 
  scale_x_discrete(limits = c("Arboreal Soil", "Bark", "Deadwood", 
                                "Fresh Leaves", "Hypnum", "Lichen", 
                                "Orthotrichum", "Leaf Litter", "Soil"), 
                   expand = expand_scale()) +
  scale_color_manual(values = c("indianred4", "darkslategray"), 
                     limits = c("Cercozoa", "Oomycota")) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() + 
  labs(title = "Alpha Diversity of Samples", 
       y = "Simpson Index (1-D)", 
       x = NULL, 
       subtitle = "Cercozoa & Oomycota") +
  theme(axis.text.x = element_text(size=12, face = "bold"), 
        axis.title = element_text(size=14, face = "bold"), 
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5), 
        plot.subtitle = element_text(size = 14, hjust = 0.5), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14, face = "bold"))
```

    ## Warning: `expand_scale()` is deprecated; use `expansion()` instead.

``` r
g_both
```

![](AlphaBoxplot_files/figure-markdown_github/CombinedAlphaBoxplot-1.png)

For the final visualisation, we combine the two plots not within the same plot like above, but next to each other:

``` r
g$theme$axis.text.x$angle = 45
g$theme$axis.text.x$hjust = 1
g$theme$axis.text.x$vjust = 1
g_cerco$theme$axis.text.x$angle = 45
g_cerco$theme$axis.text.x$hjust = 1
g_cerco$theme$axis.text.x$vjust = 1

g$theme$plot.subtitle$hjust = 0
g_cerco$theme$plot.subtitle$hjust = 0

combi = ggarrange(g_cerco, g, 
                  labels = c("A", "B"), 
                  ncol = 2, nrow = 1, 
                  align = "h") #%>%
  #annotate_figure(top = text_grob("Alpha diversity of samples", 
  #                                face = "bold", size = 20), 
  #                fig.lab = "Figure X", fig.lab.face = "bold", 
  #                fig.lab.size = 18)

#ggsave("AlphaBoxplotCombined.tif", plot = combi, 
#       device = "tiff", dpi = 600, width = 16, height = 9, 
#       units = "cm")
ggsave("AlphaBoxplotCombined.png", plot = combi, 
       device = "png", dpi = 300, width = 16, height = 9, 
       units = "cm")
ggsave("AlphaBoxplotCombined.jpeg", plot = combi, 
       device = "jpeg", dpi = 300, width = 16, height = 9, 
       units = "cm")
ggsave("AlphaBoxplotCombined.pdf", plot = combi, 
       device = "pdf", dpi = 300, width = 16, height = 9, 
       units = "cm")

combi
```

![](AlphaBoxplot_files/figure-markdown_github/unnamed-chunk-1-1.png)
