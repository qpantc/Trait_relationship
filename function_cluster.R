
library(tidyverse)
library(smatr)

# function for sample integer from a range with minimal value and maximal value.
get_random_integer <- function(maximal = 100,minimal = 1,size = 10){
  integer_list <- c()
  while (length(integer_list) < size) {
    one_integer <- as.integer(runif(1, min = minimal, max = maximal))
    if (!(one_integer %in% integer_list) & !is.na(one_integer)){
      integer_list <- c(integer_list,one_integer)
    }
  }
  return(integer_list)
}
# get_random_integer(maximal =1000,size = 20)
# get_random_integer(maximal = 10000000000000,size = 30)

# call calculate the distance of a pair of trait relationship with SMA regression
sma_test <- function(sma_data,sample_times = 10){
  colnames(sma_data) <- c('cluster','x','y')
  cluster_list <- unique(sma_data$cluster)
  if (any(table(sma_data['cluster']) < 3) | length(cluster_list) < 2){
    # if the some trait records are less than 3, there will be no SMA
    pair_distence = 0.5
  }else{
    # do SMA to calculate the distance based on their regression similarity
    cluster_1_data <- subset(sma_data,cluster == cluster_list[1])
    cluster_2_data <- subset(sma_data,cluster == cluster_list[2])
    
    # preparing sampling for the different size of trait records
    sample_size <- min(nrow(cluster_1_data),nrow(cluster_2_data))
    
    max_cluster_1 <- tryCatch(expr = {
      max_cluster_1 <- c(1:choose(nrow(cluster_1_data),sample_size))
      },
     error = function(e) {
       cat("An error occurred when set the maximal samples:", conditionMessage(e), "\n")
       return(c(1:1000000)) })
    
    max_cluster_2 <- tryCatch(expr = {
      max_cluster_2 <- c(1:choose(nrow(cluster_2_data),sample_size))
      },
    error = function(e) {
      cat("An error occurred when set the maximal samples:", conditionMessage(e), "\n")
      return(c(1:1000000)) })
    
    # start sample to get the average results
    distance_sample_list <- c()
    for (i_ in max_cluster_1) {
      if (length(max_cluster_1) > 10000){
        cluster_1_data_sampled <- cluster_1_data[get_random_integer(maximal = nrow(cluster_1_data),size = sample_size),]
      }else{
        cluster_1_data_sampled <- cluster_1_data[sample(c(1:nrow(cluster_1_data)),sample_size),]
      }
      for (j_ in max_cluster_2) {
        if (length(max_cluster_2) > 10000) {
          cluster_2_data_sampled <- cluster_2_data[get_random_integer(maximal = nrow(cluster_2_data),size = sample_size),]
        }else{
          cluster_2_data_sampled <- cluster_2_data[sample(c(1:nrow(cluster_2_data)),sample_size),]
        }
        
        # merge all data
        all_sampled_data <- rbind(cluster_1_data_sampled,cluster_2_data_sampled)
        all_sampled_data$cluster <- "0"
        
        # if there still some reason that the SMA can't be ferformed, we skip the sample
        tryCatch(
          expr = {
            one_distance <- 0
            
            # check the significance for each cluster
            # cluster 1 
            cluster1_check <- rbind(all_sampled_data,cluster_1_data_sampled)
            # slope
            slope_test <- sma(y~x*cluster,data=cluster1_check)
            slope_significance <- slope_test$commoncoef$p
            # elevation
            elevation_test <- sma(y~x+cluster,type = "elevation", data=cluster1_check)
            intercept_significance <- elevation_test$gtr$p[[1]]
            
            one_distance = ifelse(slope_significance < 0.05,one_distance+0.25,one_distance)
            one_distance = ifelse(intercept_significance < 0.05,one_distance+0.25,one_distance)
            
            # cluster 2
            cluster2_check <- rbind(all_sampled_data,cluster_2_data_sampled)
            # slope
            slope_test <- sma(y~x*cluster,data=cluster2_check)
            slope_significance <- slope_test$commoncoef$p
            # elevation
            elevation_test <- sma(y~x+cluster,type = "elevation", data=cluster2_check)
            intercept_significance <- elevation_test$gtr$p[[1]]
            
            one_distance = ifelse(slope_significance < 0.05,one_distance+0.25,one_distance)
            one_distance = ifelse(intercept_significance < 0.05,one_distance+0.25,one_distance)
            
            distance_sample_list <- c(distance_sample_list,one_distance)
          },
          error = function(e) {
            # Code to execute when an error occurs
            cat("An error occurred when perfrom SMA:", conditionMessage(e), "\n")
            return(NULL)  # You can return a default value or do other error handling here
          }
        )
        
        # when we have enough samples, we stop
        if (length(distance_sample_list) > sample_size/10 & length(distance_sample_list) > sample_times){
          break
        }
      }
      # check again when we have enough samples, we stop
      if (length(distance_sample_list) > sample_size/10 & length(distance_sample_list) > sample_times){
        break
      }
    }
    # if the sma is not successed, we let the distance as 0.5
    pair_distence = ifelse(length(distance_sample_list)>0,mean(distance_sample_list),0.5)
  }
  return(pair_distence)
}

# 定义一个需要并行计算的函数，接受多个参数
piar_parallel <- function(pari_i,TRAIT_PAIRS, two_cluster_species_df,sample_times) {
  one_trait_pair <- TRAIT_PAIRS[,pari_i]
  sma_data <- two_cluster_species_df[c('cluster',one_trait_pair)] %>% na.omit()
  distance_pair <- sma_test(sma_data, sample_times = sample_times)
  return(distance_pair)
}