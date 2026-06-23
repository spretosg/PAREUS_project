source("WP2/2_ecosystem_service/1_code/utils.R")
library(tidyr)
library(tidyverse)
library(dplyr)


stud_id<-"FRA_BAR2"
main_path<-paste0("P:/312204_pareus/WP2/T2.2/PGIS_ES_mapping/",stud_id,"/raw_data_backup")

es_pair<-read.csv(paste0(main_path,"/es_pair.csv"))%>%filter(siteID == stud_id)

#make analysis groups
att_main <- c("Cultural","Provisioning","Regulation")
att_gr1 <-  c("erosion","flood","habitat")
att_gr2 <- c("aest","sense","recr")
att_gr3<-c("farm","wild_hunt","wild_col","mat")
att_list<-list(att_gr1,att_gr2,att_gr3,att_main)

grp_names_vec<-c("regulating","cultural","provisioning","main")

#vector with individual decision makers (userIDs)
dec_makers<-es_pair%>%distinct(userID)

pref_list<-list() #list to store tmp pref
agree_list<-list() #list to store tmp consistencies
#the same order as atts_list
group_vec<-es_pair%>%distinct(ahp_section)



weights_list <- list()
cr_vec <- c()


for(g in 1:nrow(group_vec)){
  
  num_grp <- as.numeric(group_vec[g,])
  tmp_attr <- att_list[[num_grp]]   # <-- MUST come before matrix build
  
  tmp_dat <- es_pair %>%
    dplyr::filter(ahp_section %in% num_grp) %>%
    dplyr::select(userID, selection_val, ES_left, ES_right)
  
  tmp_dat$pair <- paste0(tmp_dat$ES_left, "_", tmp_dat$ES_right)
  
  tDat <- tmp_dat %>%
    dplyr::select(selection_val, pair, userID) %>%
    dplyr::group_by(pair, userID) %>%
    summarise(selection_val = mean(selection_val, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(id_cols = userID,
                names_from = pair,
                values_from = selection_val)
  
  tDat <- tDat[, -1]  # remove userID
  
  # --- AHP replacement ---
  weights_list <- list()
  cr_vec <- c()
  
  for(i in 1:nrow(tDat)) {
    
    mat <- build_ahp_matrix(tDat[i, ], tmp_attr, negconvert = TRUE)
    res <- ahp_eigen(mat)
    cr  <- ahp_cr(res$lambda, length(tmp_attr))
    
    weights_list[[i]] <- res$weights
    cr_vec[i] <- cr
  }
  
  # Optional consistency filtering (recommended)
  # keep <- cr_vec < 0.1
  # weights_list <- weights_list[keep]
  
  # geometric mean aggregation (matches ahpsurvey)
  pref_tmp <- aggregate_geom_mean(weights_list)
  
  # match original output structure
  a <- as.data.frame(pref_tmp)
  rownames(a) <- tmp_attr
  colnames(a) <- "pref_tmp"
  
  pref_list[[g]] <- a
  
  # consistency output (same idea as ahp.cr output)
  df_cr <- data.frame(CR = cr_vec)
  
  # --- consensus (unchanged) ---
  icc_res <- irr::icc(as.data.frame(tDat),
                      model = "twoway",
                      type = "agreement",
                      unit = "average")
  
  agree_list[[g]] <- c(icc_res$value,
                       grp_names_vec[num_grp],
                       stud_id)
}


# OPTIONAL: filter inconsistent respondents (same logic users often apply)
# weights_list <- weights_list[cr_vec < 0.1]

# geometric mean aggregation (matches ahpsurvey)
pref_tmp <- aggregate_geom_mean(weights_list)

df_cr <- data.frame(CR = cr_vec)




ind_prov<-3
ind_cult<-4
ind_all<-1
ind_reg<-2


main<-as.data.frame(t(pref_list[[ind_all]]))

cul<-as.data.frame(pref_list[[ind_cult]])
cul$global<-main$Cultural
cul$class<-rep("cul",nrow(cul))
prov<-as.data.frame(pref_list[[ind_prov]])
prov$global<-main$Provisioning
prov$class<-rep("prov",nrow(prov))

reg<-as.data.frame(t(pref_list[[ind_reg]]))
reg<-as.data.frame(pref_list[[ind_reg]])
reg$global<-main$Regulation
reg$class<-rep("reg",nrow(reg))

## final ranking of es per respondent

all<-as.data.frame(rbind(reg,prov,cul))
n_es_grp<-table(all$class)
sum_glob<-sum(all$global)

reg$adj_glob<-n_es_grp[3]*(reg$global/sum_glob)
cul$adj_glob<-n_es_grp[1]*(cul$global/sum_glob)
prov$adj_glob<-n_es_grp[2]*(prov$global/sum_glob)
all_fin<-rbind(reg,cul,prov)
all_fin$pref_adj<-all_fin$pref_tmp*all_fin$adj_glob
# all<-t(reg*main$Regulation)



all_fin<-all_fin%>%rownames_to_column(var = "esID")
all_fin$siteID<-rep(stud_id,nrow(all_fin))
sum(all_fin$pref_adj)

write.csv(all_fin,paste0(main_path,"/ahp_weights.csv"))
