############################################################################
# Pareus WP3
# Date: 06.08.2026
# Author: Roel May (Norwegian Institute for Nature Research)
############################################################################

library(tidyverse)
library(Hmisc)
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
library(ggridges)
library(ahpsurvey)

################################
## Prepare the data for analysis
################################

# 1. Read all sheets
country <- "France" # Norway, Slovakia, France

esnms <- read_excel("data/WP3/PCA_info.xlsx", sheet = "esnames")
plannms <- read_excel("data/WP3/PCA_info.xlsx", sheet = "plannames")
plannms <- subset(plannms,Country==country)
plannms <- plannms[plannms$Include==1,]

path <- paste0("data/WP3/PCA_",country,".xlsx")
sheet_names <- excel_sheets(path)
tools_list <- map(sheet_names, ~ read_excel(path, sheet = .x))

tools_list <- lapply(tools_list, function(x) subset(x,tool %in% plannms$Column))

# 2. Convert to numeric BEFORE merging
tools_list <- map(
  tools_list,
  ~ .x %>%
    mutate(across(starts_with("es_"), ~ suppressWarnings(as.numeric(.x))))
)

# 3. Merge lans
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

################################
## Produce data figures
################################

tools_df <- tools_mean
tools_df <- pivot_longer(tools_df, 
                        cols = 3:ncol(tools_df), 
                        names_to = "ES", 
                        values_to = "Score")
tools_df$Score <- as.numeric(tools_df$Score)
tools_df <- na.omit(tools_df)
tools_df$tool <- plannms$Code[match(tools_df$tool,plannms$Column)]
tools_df$tool <- factor(tools_df$tool,levels=rev(plannms$Code))
tools_df$ES <- esnms$Code[match(tools_df$ES,esnms$Column)]
tools_df$ES <- factor(tools_df$ES,levels=rev(esnms$Code))



p3 <- ggplot(data=tools_df, aes(y=ES, x=Score)) +
  geom_vline(xintercept = 0,linetype = "dashed") +
  geom_density_ridges_gradient(
    aes(fill = ..x..), scale = 1) +
  scale_fill_gradientn(
    colours = c("#D55E00", "#F0E442", "#009E73"),
    name = "Score",
    breaks = c(-3,0,3),
    limits=c(-3,3),
    labels=c("Negative","Neutral","Positive")
  )+
  scale_x_continuous(breaks = seq(-3, 3, by = 1),limits=c(-3,3)) +
  labs(title = "Ecosystem services", y=country) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size=12),
        plot.title = element_text(size=18),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=24),
        legend.title = element_text(size=14,face = "bold"),
        legend.text = element_text(size=12)
  )

p4 <- ggplot(data=tools_df, aes(y=tool, x=Score)) +
  geom_vline(xintercept = 0,linetype = "dashed") +
  geom_density_ridges_gradient(
    aes(fill = ..x..), scale = 1) +
  scale_fill_gradientn(
    colours = c("#D55E00", "#F0E442", "#009E73"),
    name = "Score",
    breaks = c(-3,0,3),
    limits=c(-3,3),
    labels=c("Negative","Neutral","Positive"),
    guide = guide_legend(override.aes = list(alpha=0))
  )+
  guides(fill=guide_legend(override.aes = list(fill="transparent",color="transparent"))) +
  scale_x_continuous(breaks = seq(-3, 3, by = 1),limits=c(-3,3)) +
  scale_y_discrete(position = "right") +
  labs(title = "Key planning tools", y=NULL) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size=12),
        plot.title = element_text(size=18),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_text(color = "transparent"),
        legend.title = element_text(color = "transparent"),
        legend.key = element_rect(fill = "white", color = "white")
  )

tiff(file=paste("outputs/WP3/Figure_Data_",country,".tiff",sep=""),width=4000,height=1200,units="px",res=400,compression="lzw")
gridExtra::grid.arrange(p3,p4,ncol=2)
dev.off()

################################
## Produce the PCA figures
################################

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
plot.pca <- ggplot() +
  # tools (points)
  # geom_point(
  #   data = pca_scores,
  #   aes(x = PC1, y = PC2),
  #   size = 2
  # ) +
  geom_text_repel(
    data = pca_scores,
    aes(x = PC1, y = PC2, label = plannms$Code, colour=factor(km.pl$cluster)),
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
    aes(x = PC1, y = PC2, label = esnms$Code),
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
  #  title = "PCA Biplot: Policy Tools and Ecosystem Service Loadings",
    title = country,
    subtitle = paste("Clustering strength of Plans (",cs.pl,") and ES (",cs.es,")",sep=""),
    x = "PC1",
    y = "PC2"
  ) +
  
  theme_minimal() +
  theme(legend.position="none")

assign(paste0("plot.pca.",country),plot.pca)

tiff(file=paste0("outputs/WP3/Figure_PCA_",country,".tiff"),width=3000,height=3000,units="px",res=400,compression="lzw")
plot.pca
dev.off()

################################
## Produce the correlation plots
################################

# Correlation plot with plan or ES clusters:
tps <- c("es","pl") # Change this trigger to either 'pl' or 'es' to adjust the plot
for(type in tps){
  dat <- es_data
  if(type=="pl"){dat <- t(es_data)}
  clus <- get(paste0("km.",type))$cluster
  names(clus) <- colnames(dat)
  
  cormat <- round(cor(dat),2)
  ord.clus <- clus[order(clus)]
  cormat <- cormat[names(ord.clus),names(ord.clus)]
  diag(cormat) <- 1.1
  
  if(type=="pl"){rownames(cormat) <- colnames(cormat) <- plannms$Code}
  if(type=="es"){rownames(cormat) <- colnames(cormat) <- esnms$Code}
  
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
          plot.title = element_text(size=18),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=24),
          panel.grid.major = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank(),
          #    legend.justification = c(1, 0),
          #    legend.position = c(0.6, 0.7),
          #    legend.direction = "horizontal"
    )+
    ggtitle(paste0(ifelse(type=="pl","Policy","ES")," bundle synergies and trade-offs"))+
    ylab(ifelse(type=="es",country,"")) +
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
  if(type=="pl"){
    plot.cor <- plot.cor +
      scale_y_discrete(position = "right") +
   #   guide = guide_legend(override.aes = list(alpha=0)))+
  guides(fill=guide_legend(override.aes = list(fill="transparent",color="transparent"))) +
  theme(legend.text = element_text(color = "transparent"),
        legend.title = element_text(color = "transparent"),
        legend.key = element_rect(fill = "white", color = "white"))
        
  }
#  plot.cor
  
  assign(paste0("plot.cor.",type),plot.cor)
}

tiff(file=paste("outputs/WP3/Figure_Cor_",country,".tiff",sep=""),width=6000,height=3000,units="px",res=400,compression="lzw")
gridExtra::grid.arrange(plot.cor.es,plot.cor.pl,ncol=2)
dev.off()

################################
## Produce the piechart plots
################################
type <- c("plan","es")

# Calculate the overall consensus across respondents
# based on weight matrix; columns are criteria; rows are respondents
ahp.consensus <- function(w){
  M <- 9
  n <- ncol(w)
  k <- nrow(w)
  
  Hamin <- -(M/(n+M-1))*log((M/(n+M-1))) - ((n-1)/(n+M-1))*log((1/(n+M-1)))
  cor <- exp(-Hamin)/n
  
  Ha <- sum(rowSums(-w*log(w)))/k
  Hg <- sum(-colMeans(w)*log(colMeans(w)))
  S <- exp(-(Hg-Ha))
  
  round(Scons <- (S-cor)/(1-cor),3)
}

for(xtype in type){
  dat <- tools_mean[,-c(1:2)]
  rownames(dat) <- tools_mean[,2]
  if(xtype=="es"){dat <- as.data.frame(t(dat))}
  R <- nrow(dat)
  kpt <- rownames(dat)
  kpt <- gsub("_","-",kpt)
  res <- as.data.frame(matrix(NA,nrow=ncol(dat), ncol=0))
  rownames(res) <- colnames(dat)
  colnames(dat) <- gsub("_","-",colnames(dat))
  
  for (r1 in 1:(R-1)){
    for (r2 in (r1+1):R){
      tmp <- t(dat[r1,] - dat[r2,])
      colnames(tmp) <- paste0(kpt[r1], "_", kpt[r2])
      res <- cbind(res,tmp)
    }
  }
  
  res <- (round((abs(res) / 6 * 8), 0) + 1)^sign(res)
  
  atts <- get(paste0(xtype,"nms"))$Code
  
  out <- ahp(df=res, atts=kpt, negconvert = FALSE, agg = TRUE)
  out$aggpref <- cbind.data.frame(out$aggpref,kpt=atts)
  
  out$aggpref
  Cons <- round(ahp.consensus(out$indpref[,1:7]),2) # should be above 0.7
  CR <- round(mean(out$indpref$CR),2) # should be below 0.1
  
  out$aggpref$kpt <- factor(out$aggpref$kpt,levels=atts)
  
  out$aggpref <- out$aggpref %>%
    mutate(
      prop = AggPref/sum(AggPref)*100,
      ypos = cumsum(prop)- 0.5*prop,
      angle = 90 - (ypos / 100 * 360) 
    )
  
  hj <- rep(1,nrow(out$aggpref))
  hj[which(out$aggpref$angle < -90)] <- 0
  out$aggpref$angle[which(out$aggpref$angle < -90)] <- out$aggpref$angle[which(out$aggpref$angle < -90)]+180
  
  p.pie <- ggplot(out$aggpref, aes(x="", y=prop, fill=rev(kpt))) +
    geom_col(width=1) +
    geom_text(aes(x=1.5, y = ypos, label = paste0("  ",kpt, "  ", "\n","  ", round(prop,0), "%  "), angle = angle),
              size = 5, color = "black",hjust=hj
    ) +
    coord_polar("y", start=0) +
    ggtitle(ifelse(xtype=="plan","Key planning tools","Ecosystem services"),subtitle=paste0("Consensus: ",Cons,"\nConsistency:",CR)) +
    #  annotate("text", x = 0, y = 0,label = country, 
    #    angle = 90, hjust = 0.5,vjust = -8,size = 10,fontface = "bold") +
    theme_void() +
    theme(#text = element_text(size=24,face = "bold"),
      plot.title=element_text(size=18,face = "bold",hjust=0.5),
      plot.subtitle = element_text(size=12,hjust=0.5,face = "italic",margin = margin(t = 0, r = 0, b = -25, l = 0))) +
    theme(legend.position="none")
  
  assign(paste0("plot.pie.",xtype),p.pie)
}

pt <- ggplot(data.frame(x=1,y=1,my_text=country), aes(x=x,y=y,label=my_text)) +
             geom_text(size=10,fontface = "bold",angle=90) +
             theme_void()

tiff(file=paste("outputs/WP3/Figure_Pies_",country,".tiff",sep=""),width=4000,height=2000,units="px",res=400,compression="lzw")
gridExtra::grid.arrange(pt,plot.pie.es,plot.pie.plan,ncol=3,widths=c(1,4,4))
dev.off()
