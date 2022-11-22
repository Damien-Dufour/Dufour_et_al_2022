Correlation DEG and PCA, Senp2cKO mice
================
Damien, Dufour

# Description

Here I described the methods I used to produce some of the figures of
Dufour et al, 2022

This is based on data from RNA-seq performed on male and female
4-week-old WT and *Senp2cKO* adrenals.

# Data availability

The raw files and table containing normalised counts are available on
the GEO platform under the accession code GSE193480

# Required packages

``` r
library(ggplot2)
library(tidyverse)
library(dplyr) 

library(pheatmap) # to plot nice heatmaps
library(viridis)  # with nice colours 

library(ggrepel) # so that text can be displayed on the scatter plot 

library(Vennerable) # to plot the Venn/euler diagram
library(broom)
```

# 1. Load and clean the data

``` r
#load the dataframe available on the GEO
raw <- read.csv("S21045_alldata.tsv", 
                sep ="\t",
                dec = ",")


# create a bool column T for duplicated F for unique
raw <- data.frame(raw,
                   duplicated(raw$Gene.name))   


# select only the non duplicated protein coding genes 
raw <- subset(raw, 
               duplicated.raw.Gene.name. == "FALSE") %>%
   subset(Gene.biotype == "protein_coding") 

rpkm <-raw[,c(34:49)] # extract the read counts

FC <- raw[,c(60:62,69:71)] # extract the comparisons


names(rpkm) <- c("WTM1",
                "WTM2",
                "WTM3",
                "WTM4",
                "KOM1",
                "KOM2",                  #rename cols
                "KOM3",
                "KOM4",
                "WTF1",
                "WTF2",
                "WTF3",
                "WTF4",
                "KOF1",
                "KOF2",                                                         
                "KOF3",
                "KOF4")
row.names(rpkm) <- raw$Gene.name  

names(FC) <- c("log2FoldChangeM",
                "pvalM",
                "padjM",
                "log2FoldChangeF",
                "pvalF",
                "padjF")

df <- data.frame(rpkm,
                   FC)              #merge rpkm and analysis

df$Gene.name <- row.names(df)
```

# 2. Plot the data

``` r
# create a column with specific name based on FC and P-value
# this will come in handy to easily display color based on these parameter

df$col <- "other"
df$col[df$padjF <= 0.05 & abs(df$log2FoldChangeF) >= 0.58] <- "Regulated genes in female"
df$col[df$padjM <= 0.05 & abs(df$log2FoldChangeM) >= 0.58] <- "Regulated genes in male"
df$col[df$padjM <= 0.05 & df$padjF <= 0.05 & abs(df$log2FoldChangeM) >= 0.58 & abs(df$log2FoldChangeF) >= 0.58] <- "Regulated genes in both conditions"

# create an empty column for gene name that will be displayed 

df$nameview <- NA

# list of genes to be displayed 

genelist <- c("Akr1d1", "Hsd3b1", "Star", "Cdkn1a", "Cyp11b2", "Cyp11a1", "Cyp11b1") 

# loop to assign the gene names in the nameview column 

for(i in genelist){
   #if(df[i,24] == "Regulated genes in both conditions")  #uncomment to select only the genes that are regulated in both conditions 
      df[i,25] <- i 
}

# create a column which takes as value the mean expression of the genes in all condition
# this is usefull to change the color based on this parameter and display expression on top of FC

df$mean <- rowMeans(df[,c(1:16)])


# create the scatter plot

ggplot(df)+
   aes(x = -log2FoldChangeM,
       y = -log2FoldChangeF,
       fill = col, #could be nameview or mean with a continuous fill scale
       #label = nameview, #uncomment to display names from nameview col
       col = col)+
   geom_point(shape = 21,
              size = .5)+
   scale_fill_manual(values = c("#cad6e2",
                                "dark red", 
                                "red",
                                "black"))+
   scale_colour_manual(values = c("#cad6e2",
                                  "dark red", 
                                  "red",
                                  "black"))+
   theme_classic(base_size = 8)+
   theme(legend.position = c(0.32,1),
         legend.title = element_blank(),
         legend.text = element_text(size = 4),
         legend.key.height = unit(0.2,"cm"),
         legend.background = element_blank(),
         legend.spacing.y = unit(1, 'cm'),
         legend.spacing.x = unit(0, 'cm')
   )+
   geom_hline(yintercept = 0,
              linetype = 2,
              size = 0.05)+
   geom_vline(xintercept = 0,
              linetype = 2,
              size = 0.05)+
#geom_text_repel(size=.75,max.overlaps = 10,segment.colour = NA,force = 0.05)+ #uncomment to display names from nameview col
   geom_smooth(inherit.aes = F,
               aes(x = -log2FoldChangeM,
                   y = -log2FoldChangeF),
               colour = "black",
               method = "lm",
               size = 0.5)
```

# 2. Test correlation

## 2.1 Test normality of the two parameters

``` r
# qqnorm allows visualisation of normality and can be considered more reliable for so many samples
# Kolmogoroff-Smrinoff test is used since Shapito & Wilk cannot handle sample higher than 5000


qqnorm(df$log2FoldChangeM)
ks.test(df$log2FoldChangeM,
        pnorm)

qqnorm(df$log2FoldChangeF)
ks.test(df$log2FoldChangeF,
             pnorm)
```

## 2.2 Test correlation between Log2(FC) in males and females

``` r
# We use Kendall method because of non-normality of the samples and presence of ties (prevents use of Spearmann)

cor.test(df$log2FoldChangeM,
         df$log2FoldChangeF,
         method = "kendall")
```

# 3. Principal component analysis

``` r
#select only expressed genes and remove NA value

PCAdf <- subset(rpkm,
             rowMeans(rpkm) != 0)%>%
   t() %>% 
   data.frame() %>%
  drop_na()

# create the PCA object 

pca_fit <- PCAdf %>%
  select(where(is.numeric)) %>%
  scale() %>%
  prcomp()

# get the implication of every component into observed variance

summary(pca_fit)

# create a dataframe with all information that will be plottable 

x <- pca_fit %>%
  augment(PCAdf) %>%
  rename_at(vars(starts_with(".fitted")),
            list(~str_replace(.,".fitted","")))

# remove all the columns corresponding to read counts

x <- x[,c(18432:18447)]

# create a column for genotypes

x$Genotype <- c(rep("WTM",4),
                rep("KOM",4),
                rep("WTF",4),
                rep("KOF",4))


# plot PC1 vs PC2 

ggplot(x)+
  aes(x = PC1,
       y = PC2,
      fill = Genotype,
      color = Genotype)+
  geom_point(shape = 21)+
  theme_classic(base_size = 6)+
  scale_fill_manual(values = c("#a00000",
                               "black",
                               "white",
                               "white"))+
  scale_color_manual(values = c("black",
                               "black",
                               "#a00000",
                               "black"))+
  stat_ellipse()
```
