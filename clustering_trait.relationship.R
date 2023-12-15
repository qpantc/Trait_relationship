
rm(list = ls())

if(!"tidyverse" %in% rownames(installed.packages())){install.packages("tidyverse")}
if(!"foreach" %in% rownames(installed.packages())){install.packages("foreach")}
if(!"doParallel" %in% rownames(installed.packages())){install.packages("doParallel")}
if(!"smatr" %in% rownames(installed.packages())){install.packages("smatr")}

# 指定要清空的文件夹路径
folder_path <- "./Results/clustering/"

# 获取文件夹中的所有文件
files_to_remove <- list.files(folder_path, full.names = TRUE)

# 删除文件
if (length(files_to_remove) > 0) {
  file.remove(files_to_remove)
  print("Files in the folder have been removed.")
} else {
  print("Folder is already empty.")
}

library(tidyverse)
library(foreach)
library(doParallel)
source("./function_cluster.R")

registerDoParallel(cores = 32)
sample_times = 10
sample_increasement = 0.2

used_traits <- c("LSLA", "LNC.m", "LL", # "LSC", "LPC.m",
                 "SSD", "SCdia","SBthick",
                 "RNC",  "RSRL", "RTD",
                 "Max.H", "RD","SD","LA","SdM")#Rdia, 

# DATA LOAD
# species 
# analysised_species <- head(read.csv('./Data/species_list_analysed.csv'),10)
analysised_species <- read.csv('./Data/species_list_analysed.csv')

analysised_species$cluster <- as.numeric(rownames(analysised_species))

# trait data
gapfilled_species_trait <- read.csv('./Data/woody_gapfilled_species_trait_controled_&_loged.csv') %>% 
  filter(species %in% analysised_species$species)
gapfilled_species_trait <- left_join(gapfilled_species_trait,analysised_species)

# TRAIT PAIRS
TRAIT_PAIRS <- combn(used_traits,2)


# 计算距离矩阵
cluster_list <- unique(analysised_species$cluster)
n <- length(cluster_list)
distance_matrix <- matrix(Inf, n, n)


for (i in 1:(n-1)) {
  cluster_i <- cluster_list[i]
  for (j in (i+1):n) {
    
    cluster_j <-  cluster_list[j]
    two_cluster_species_df <- subset(gapfilled_species_trait,cluster %in% c(cluster_i,cluster_j))
    
    # trait_pair_distence <- list()
    # for (pari_i in c(1:length(TRAIT_PAIRS[1,]))) {
    #   one_trait_pair <- TRAIT_PAIRS[,pari_i]
    #   sma_data <- two_cluster_species_df[c('cluster',one_trait_pair)] %>% na.omit()
    #   distance_pair <- sma_test(sma_data, sample_times = sample_times)
    #   trait_pair_distence[paste(one_trait_pair,collapse = '_')] = distance_pair
    #   print(paste(c(pari_i,one_trait_pair,distance_pair)))
    # }
    trait_pair_distence <- foreach(params = 1:length(TRAIT_PAIRS[1,]), .combine = c) %dopar% {
      piar_parallel(pari_i = params, TRAIT_PAIRS = TRAIT_PAIRS,
                    two_cluster_species_df=two_cluster_species_df,sample_times=sample_times)}
    
    print(paste(Sys.time(),": ",cluster_i,cluster_j,mean(trait_pair_distence)))
    # distance_matrix[i, j] <- round(rowMeans(as.data.frame(trait_pair_distence)),4)
    distance_matrix[i, j] <- round(mean(trait_pair_distence),4)
    distance_matrix[j, i] <- distance_matrix[i, j]
  }
}

save(distance_matrix, cluster_list, sample_times, file = "./Results/clustering_data.RData")

load("./Results/clustering_data.RData")

while (length(cluster_list)>1){
  # 找到最小距离的簇
  # merge_clusters <- c(0, 0)
  
  min_dist <- min(distance_matrix)
  
  merge_list <- list()
  merged_species <- c()
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      if (distance_matrix[i, j] <= min_dist) {
        merged_cluster <- cluster_list[c(i,j)]
        break
        # merge_clusters <- c(i, j)
      }
    }
  }
  
  # update the cluster list
  new_cluster <- max(as.numeric(cluster_list)) + 1
  # for the trait dataset 
  gapfilled_species_trait$cluster[which(gapfilled_species_trait$cluster %in% merged_cluster)] <- new_cluster
  
  analysised_species$old_cluster <- analysised_species$cluster
  analysised_species$cluster[which(analysised_species$cluster %in% merged_cluster)] <- new_cluster
  analysised_species$distance = min_dist
  write.csv(analysised_species,paste0('./Results/clustering/',new_cluster,'_',min_dist,'.csv'),row.names = FALSE)
  
  # update the new distance matrix
  new_cluster_list <- unique(analysised_species$cluster)
  n <- length(new_cluster_list)
  new_distance_matrix <- matrix(Inf, n, n)
  
  if (n <=2) {break}
  # if they are including the new cluster, recalculate the  distance
  for (i in 1:(n-1)) {
    cluster_i <- new_cluster_list[i]
    # first check if the cluster_i is new cluster?
    if (cluster_i == new_cluster){
      # yes, recalculate
      for (j in (i+1):n) {
        cluster_j <-  new_cluster_list[j]
        
        trait_pair_distence <- list()
        two_cluster_species_df <- subset(gapfilled_species_trait,cluster %in% c(cluster_i,cluster_j))
        
        trait_pair_distence <- foreach(params = 1:length(TRAIT_PAIRS[1,]), .combine = c) %dopar% {
          piar_parallel(pari_i = params, TRAIT_PAIRS = TRAIT_PAIRS,
                        two_cluster_species_df=two_cluster_species_df,sample_times=sample_times)}
        
        print(paste(Sys.time(),": ",cluster_i,cluster_j,mean(trait_pair_distence)))
        
        new_distance_matrix[i, j] <- round(mean(trait_pair_distence),4)
        new_distance_matrix[j, i] <- new_distance_matrix[i, j]
      }
      
    }
    else{
      # no?  then check if the cluster_j is new cluster?
      for (j in (i+1):n) {
        cluster_j <-  new_cluster_list[j]
        if (cluster_j == new_cluster){
          # if yes:
          trait_pair_distence <- list()
          two_cluster_species_df <- subset(gapfilled_species_trait,cluster %in% c(cluster_i,cluster_j))
          
          trait_pair_distence <- foreach(params = 1:length(TRAIT_PAIRS[1,]), .combine = c) %dopar% {
            piar_parallel(pari_i = params, TRAIT_PAIRS = TRAIT_PAIRS,
                          two_cluster_species_df=two_cluster_species_df,sample_times=sample_times)}
          
          print(paste(Sys.time(),": ",cluster_i,cluster_j,mean(trait_pair_distence)))
          
          new_distance_matrix[i, j] <- round(mean(trait_pair_distence),4)
          new_distance_matrix[j, i] <- new_distance_matrix[i, j]
        }
        else{
          # no? now we can quire the distance from the former matrix where the cluster is the same
          new_distance_matrix[i, j] <- distance_matrix[which(cluster_list == new_cluster_list[i]), 
                                                       which(cluster_list == new_cluster_list[j])]
          new_distance_matrix[j, i] <- new_distance_matrix[i, j]
        }
      }
    }
  }
  distance_matrix <- new_distance_matrix
  cluster_list <- new_cluster_list
  sample_times <- ifelse(sample_times <= 30, sample_times + sample_increasement, sample_times)
  save(distance_matrix, cluster_list, sample_times, file = "./Results/clustering_data.RData")
}

registerDoParallel()
