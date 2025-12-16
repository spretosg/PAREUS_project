############################################################################
# Pareus preliminary analysis
# Dato: 25.11.2025
############################################################################

# Bytter til norsk locale (kun hvis støttet på systemet)
Sys.setlocale("LC_ALL", "no_NO.UTF-8")

# ---------------------------------------------------------
# Pakker
# ---------------------------------------------------------
library(tidyverse)
library(readxl)
library(summarytools)
library(writexl)
library(magrittr)
library(cluster)
library(factoextra)
library(reshape)
library(ggplot2)
library(ggrepel)
library(ggforce)

# 1. Read all sheets
# path <- "C:/Users/trond.simensen/NINA/312204 - PAREUS - Documents/WP3/TASK_3_2/POLICY_COHERENCE_NOR/Interviews with practitioners/results_clean 2025.xlsx"
path <- "C:/Users/roel.may/OneDrive - NINA/Werk/Projects/PAREUS/PCA_slovakia.xlsx"
path <- "C:/Users/roel.may/OneDrive - NINA/Werk/Projects/PAREUS/results_clean 2025.xlsx"
sheet_names <- excel_sheets(path)
tools_list <- map(sheet_names, ~ read_excel(path, sheet = .x))

# 2. Convert to numeric BEFORE merging
tools_list <- map(
  tools_list,
  ~ .x %>%
    mutate(across(starts_with("es_"), ~ as.numeric(.x)))
)

# 3. Merge
merged <- reduce(
  tools_list,
  ~ left_join(.x, .y, by = c("tool_id", "tool"))
)

# 4. Convert ALL es_-columns to numeric AGAIN after merge
merged <- merged %>%
  mutate(across(matches("^es_"), ~ as.numeric(.x)))

# 5. Identify es-variables
es_vars <- grep("^es_", names(tools_list[[1]]), value = TRUE)

# 6. Compute row-wise means (NA removed automatically)
tools_mean <- map_dfc(
  es_vars,
  function(v) {
    cols <- merged %>% select(matches(paste0("^", v)))
    tibble(!!v := rowMeans(cols, na.rm = TRUE))
  }
)

# 7. Add metadata
tools_mean <- bind_cols(
  merged %>% select(tool_id, tool) %>% distinct(),
  tools_mean
)

tools_mean <- as.data.frame(tools_mean)



# PCA ---------------------------------------------------------------------

# Crawley, the R-book, s. 809

# --- 1. Prepare data: remove empty rows and set tool as rownames ---
tools_pca <- tools_mean %>%
  filter(if_any(starts_with("es_"), ~ !is.na(.))) %>%  # remove rows with all NA
  select(tool, starts_with("es_"))                     # keep relevant cols

# set rownames for PCA
tool_labels <- tools_pca$tool
es_data <- tools_pca %>% select(-tool)
rownames(es_data) <- tool_labels

# Center and scale the data
es_data <- t(es_data)
es_data <- scale(es_data,center=T,scale=T)
es_data <- t(es_data)

# --- 2. Run PCA ---
pca <- prcomp(es_data)
# pca <- prcomp(es_data, scale = TRUE)

# First look
summary(pca)
plot(pca, main = "Pareus scree-plot", col = "steelblue")
biplot(pca)

# --- 3. Extract results ---
# Variance explained
pca_summary <- summary(pca)

# Loadings (variables → PC axes)
pca_loadings <- as_tibble(pca$rotation, rownames = "variable")

# Scores (tools → PC axes), includes the tool name
pca_scores <- as_tibble(pca$x, rownames = "tool")

# Data frame for scree plot
scree_df <- tibble(
  PC = paste0("PC", seq_along(pca$sdev)),
  eigenvalue = pca$sdev^2,
  prop_var = (pca$sdev^2) / sum(pca$sdev^2),
  cum_var = cumsum((pca$sdev^2) / sum(pca$sdev^2))
)

# --- 4. Also provide biplot-ready data ---
biplot_vars  <- pca_loadings      # arrows (ecosystem services)
biplot_tools <- pca_scores        # points (tools)

biplot_vars
biplot_tools

ggplot(biplot_tools, aes(PC1, PC2, label = tool)) + geom_text()


# Hierarchical K-Means Cluster Analysis
# Both for plans and ES
# These are shown in the plot by their colours (plans) or linetype (ES)
# Also, ellipses are drawn around them in a similar fashion
wss <- (nrow(es_data)-1)*sum(apply(es_data,2,var))
nn <- nrow(es_data)
for (i in 2:nn){
  if(i<nrow(es_data)){wss[i] <- sum(hkmeans(es_data,k=i)$withinss)}
  if(i==nrow(es_data)){wss[i] <- 0}
}
dd<-(1-wss/max(wss))*nn
group <- min(which(diff(dd)<1))
km.pl <- hkmeans(es_data, group)
cs.pl <- round(km.pl$betweenss/km.pl$tot.withinss,3)

wss <- (nrow(t(es_data))-1)*sum(apply(t(es_data),2,var))
nn <- nrow(t(es_data))
for (i in 2:nn){
  if(i<nrow(t(es_data))){wss[i] <- sum(hkmeans(t(es_data),k=i)$withinss)}
  if(i==nrow(t(es_data))){wss[i] <- 0}
}
dd<-(1-wss/max(wss))*nn
group <- min(which(diff(dd)<1))
km.es <- hkmeans(t(es_data), group)
cs.es <- round(km.es$betweenss/km.es$tot.withinss,3)

# Calculate mean scores to indicate how much influence plans have
# or how much ES are being influenced
# This is indicated by the font size of the text in the plot
es_weights <- colMeans(abs(tools_mean[,-c(1:2)]))
pl_weights <- rowMeans(abs(tools_mean[,-c(1:2)]))

# --- 1. Prepare scaling factor so arrows and points fit same plot ---
arrow_scale <- 
  max(abs(pca_scores$PC1), abs(pca_scores$PC2)) /
  max(abs(biplot_vars$PC1), abs(biplot_vars$PC2))

biplot_vars_scaled <- biplot_vars %>%
  mutate(
    PC1 = PC1 * arrow_scale,
    PC2 = PC2 * arrow_scale
  )

biplot_vars_scaled$variable <- unlist(lapply(biplot_vars_scaled$variable,function(x)capitalize(tail(unlist(strsplit(x,"_")),n=1))))

# --- 2. Biplot with text labels and arrows ---
ggplot() +
  # tools (points)
  # geom_point(
  #   data = pca_scores,
  #   aes(x = PC1, y = PC2),
  #   size = 2
  # ) +
  geom_text_repel(
    data = pca_scores,
    aes(x = PC1, y = PC2, label = tool, colour=factor(km.pl$cluster)),
    size = sqrt(pl_weights)*4
#    nudge_x = sign(pca_scores$PC1)*0.1, nudge_y = sign(pca_scores$PC2)*0.1
  ) +
 
  geom_mark_ellipse(
    data = pca_scores,
    aes(x = PC1, y = PC2, colour = factor(km.pl$cluster)),
    linewidth=0.5
  ) +
  
  # loadings (arrows)
  geom_segment(
    data = biplot_vars_scaled,
    aes(x = 0, y = 0, xend = PC1, yend = PC2),
    arrow = arrow(type="closed", length = unit(0.25, "cm")),
    linewidth = 1, linetype=km.es$cluster,
    color = "darkgrey"
  ) +
  geom_text(
    data = biplot_vars_scaled,
    aes(x = PC1, y = PC2, label = variable),
    color = "darkgrey", size = sqrt(es_weights)*4,
    nudge_x = sign(biplot_vars_scaled$PC1)*0.1, nudge_y = sign(biplot_vars_scaled$PC2)*0.1
  ) +
  
  geom_mark_ellipse(
    data = biplot_vars_scaled,
    aes(x = PC1, y = PC2, group=km.es$cluster), 
    colour = "darkgrey", linewidth=0.5, linetype=km.es$cluster
  ) +
  
  # axis zero lines (optional but tidy)
  geom_hline(yintercept = 0, linetype = "solid") +
  geom_vline(xintercept = 0, linetype = "solid") +

  labs(color='Cluster') +

  labs(
    title = "PCA Biplot: Policy Tools and Ecosystem Service Loadings",
    subtitle = paste("Clustering strength of Plans (",cs.pl,") and ES (",cs.es,")",sep=""),
    x = "PC1",
    y = "PC2"
  ) +
  theme_minimal()

es_data


# Correlation plot with plan or ES clusters:
type <- "es" # Change this trigger to either 'pl' or 'es' to adjust the plot
dat <- es_data
if(type=="pl"){dat <- t(es_data)}
clus <- get(paste0("km.",type))$cluster
names(clus) <- colnames(dat)

cormat <- round(cor(dat),2)
ord.clus <- clus[order(clus)]
cormat <- cormat[names(ord.clus),names(ord.clus)]
diag(cormat) <- 1.1
melted_cormat <- reshape::melt(cormat)
colnames(melted_cormat) <- c("Var1","Var2","value")
melted_cormat$text <- melted_cormat$value
melted_cormat$text <- format(melted_cormat$text,nsmall=2)
melted_cormat$text[melted_cormat$value=="1.1"] <- ""
melted_cormat$row <- 1:nrow(melted_cormat)
if(type=="es"){
  melted_cormat$Var1 <- unlist(lapply(melted_cormat$Var1,function(x)capitalize(tail(unlist(strsplit(as.character(x),"_")),n=1))))
  melted_cormat$Var2 <- unlist(lapply(melted_cormat$Var2,function(x)capitalize(tail(unlist(strsplit(as.character(x),"_")),n=1))))
}

# Heatmap
ggheatmap <- ggplot(melted_cormat, aes(reorder(Var2,row), reorder(Var1,row), fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  geom_text(aes(reorder(Var2,row), reorder(Var1,row), label = text), color = "black", size = 4) +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
        axis.text.y = element_text(size=12),
        plot.title = element_text(size=24),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        #    legend.justification = c(1, 0),
        #    legend.position = c(0.6, 0.7),
        #    legend.direction = "horizontal"
  )+
  ggtitle(paste0(ifelse(type=="pl","Policy","ES")," bundle synergies and trade-offs"))+
  #  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
  #         title.position = "top", title.hjust = 0.5)) +
  coord_fixed()

nr.clus <- unique(clus)
loc.clus <- 1:length(clus)
plot.cor <- ggheatmap
for(c in 1:length(nr.clus)){
  plot.cor <- plot.cor + 
    annotate("rect",xmin=min(loc.clus[which(ord.clus==c)])-0.45,xmax=max(loc.clus[which(ord.clus==c)])+0.45,
             ymin=min(loc.clus[which(ord.clus==c)])-0.45,ymax=max(loc.clus[which(ord.clus==c)])+0.45,
             alpha=0,colour="black",linetype=2,size=2)
}
plot.cor

####################################################
# compare values of overall versus mean across ES...
####################################################

# Get overall matrix loaded:
path <- "C:/Users/roel.may/OneDrive - NINA/Werk/Projects/PAREUS/Results_conflictevaluation_clean.xlsx"
sheet_names <- excel_sheets(path)
tools_list <- map(sheet_names, ~ read_excel(path, sheet = .x))

# 2. Convert to numeric BEFORE merging
# tools_list <- map(
#   tools_list,
#   ~ .x %>%
#     mutate(across(starts_with("es_"), ~ as.numeric(.x)))
# )

# 3. Merge
merged <- reduce(
  tools_list,
  ~ left_join(.x, .y, by = c("tool_id", "tool"))
)

# 4. Convert ALL es_-columns to numeric AGAIN after merge
# merged <- merged %>%
#   mutate(across(matches("^es_"), ~ as.numeric(.x)))

# 5. Identify es-variables
es_vars <- colnames(tools_list[[1]])[-c(1:2)]

# 6. Compute row-wise means (NA removed automatically)
tools_all <- map_dfc(
  es_vars,
  function(v) {
    cols <- merged %>% select(matches(paste0("^", v)))
    tibble(!!v := rowMeans(cols, na.rm = TRUE))
  }
)

tools_all <- tools_all %>%
  mutate_all(~ifelse(is.nan(.), 0, .))

tools_all <- tools_all + t(tools_all)

# 7. Add metadata
# tools_all <- bind_cols(
#   merged %>% select(tool_id, tool) %>% distinct(),
#   tools_all
# )
# 
# tools_all <- as.data.frame(tools_all)

# Plot the overall plan loadings:
fviz_pca_var(prcomp(tools_all, scale=T),repel = TRUE,col.var="blue",
             arrowsize = 1,title="Overall policy plan coherence")

# Compare overall scoring to correlation of ESs:
tools_compare <- cor(t(tools_mean[,-c(1:2)]))
rownames(tools_compare) <- colnames(tools_compare) <- tools_mean$tool
incl <- which(tools_mean$tool %in% colnames(tools_all))
tools_compare <- tools_compare[incl,incl]
diag(tools_compare) <- 0

cor.test(as.matrix(tools_all),tools_compare)
