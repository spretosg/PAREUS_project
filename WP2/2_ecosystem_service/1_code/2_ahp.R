library(ahpsurvey)
library(tidyr)
library(tidyverse)
library(irr)
stud_id<-"FRL04"
main_path<-paste0("P:/312204_pareus/WP2/T2.2/PGIS_ES_mapping/",stud_id,"/raw_data_backup")

es_pair<-read.csv(paste0(main_path,"/es_pair.csv"))

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

  ## store preference and agreement

for(g in 1:nrow(group_vec)){
  num_grp<-as.numeric(group_vec[g,])
  tmp_dat<-es_pair%>%filter(ahp_section %in% num_grp)%>%select(userID,selection_val, ES_left, ES_right)
  tmp_dat$pair<-paste0(tmp_dat$ES_left,"_",tmp_dat$ES_right)
  tDat<-tmp_dat%>%select(selection_val, pair, userID)%>%
    pivot_wider(id_cols = userID, id_expand = F,  names_from = pair, values_from = selection_val)
  tDat<-tDat[, -1]

  ## preferences
  tmp_attr<-att_list[[num_grp]]

  ahp_tmp <- ahp.mat(df = tDat, atts = tmp_attr, negconvert = TRUE)
  pref_tmp<-ahp.aggpref(ahp_tmp, tmp_attr, method = "eigen", aggmethod = "eigen")

  # store these:
  df_cr<-ahp.cr(ahp_tmp, tmp_attr, ri = NULL)

  a<-as.data.frame(pref_tmp)
  pref_list[[g]]<-a
  # cons_list[group_vec[g,]]<-ahp.cr(ahp_all, atts, ri = NULL)

  ## consensus
  icc_res<-icc(as.data.frame(tDat), model = "twoway", type = "agreement", unit = "average")
  agree_list[[g]]<-c(icc_res$value,grp_names_vec[num_grp],stud_id)

}
ind_reg<-1
ind_cult<-2
ind_prov<-3
ind_all<-4

main<-as.data.frame(t(pref_list[[ind_all]]))
reg<-as.data.frame(pref_list[[ind_reg]])
reg$global<-main$Regulation
reg$class<-rep("reg",nrow(reg))
cul<-as.data.frame(pref_list[[ind_cult]])
cul$global<-main$Cultural
cul$class<-rep("cul",nrow(cul))
prov<-as.data.frame(pref_list[[ind_prov]])
prov$global<-main$Provisioning
prov$class<-rep("prov",nrow(prov))

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

write.csv(all_fin,paste0(main_path,"/ahp_weights.csv"))

# upload agreement among participants
# Convert list to data frame
data_df <- as.data.frame(do.call(rbind, agree_list))

# Rename columns
colnames(data_df) <- c("icc_val", "es_Category")

# Convert Value column to numeric
data_df$icc_val <- as.numeric(data_df$icc_val)
write.csv(data_df,paste0(main_path,"/icc_weights.csv"))


