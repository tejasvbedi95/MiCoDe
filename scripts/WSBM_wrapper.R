##################################################################################
# WSBM Wrapper function

WSBM_wrapper <- function(data, K = "auto", cor = "spearman", transform = T,
                         K_max = 20, alpha = 0.05){
  
  require(mcclust)
  require(Rcpp)
  
  # Load functions
  Rcpp::sourceCpp("scripts/SBM_cpp_v3.4.cpp") # CPP function for WSBM (auto)
  Rcpp::sourceCpp("scripts/SBM_cpp_v2.4.cpp") # CPP function for WSBM (fixed)
  source("scripts/functions.R")
  
  # data = n by p taxonomic abundance count table
  # K = "auto" for automatic inference, numeric values between 2 - 10 for fixed communities
  # cor = spearman / pearson correlation method
  # transform = TRUE for CLR transformation on count table
  # K_max = max value of K if K_max = "auto"
  # alpha = DP concentration parameter if K_max = "auto"
  # OUTPUT: a list of correlation matrix and community label vector
  
  if(transform){
    if(cor == "spearman"){
      cor_data <- count_to_cor(data, method = "spearman")$CLR
    }else{
      cor_data <- count_to_cor(data, method = "pearson")$CLR
    }
  }else{
    if(cor == "spearman"){
      cor_data <- count_to_cor(data, method = "spearman")$comp
    }else{
      cor_data <- count_to_cor(data, method = "pearson")$comp
    }
  }
  
  if(K == "auto"){
    res <- auto_WSBM(cor_data, K_max, alpha, T)
    diag(res$ppm_store) <- 1000
    clust_res <- minbinder(res$ppm_store/1000, method = "comp")$cl
  }else{
    if(K < 2 | K > 10){
      stop("Enter the value of K b/w 2 and 10 or choose 'auto' ")
    }else{
      res <- WSBM(cor_data, K, T)
      clust_res <- res$z+1
    }
  }
  
  # Relabel clusters so that the smallest label has the largest no. of observations
  
  # browser()
  
  # names(clust_res) <- 1:length(clust_res)
  # clust_res_temp <- clust_res
  # clust_true_factor <- as.factor(clust_res_temp)
  # recode_levels <- names(sort(summary(as.factor(clust_res_temp)), T))
  # levels(clust_true_factor) <- recode_levels
  # clust_res <- as.numeric(as.character(clust_true_factor))
  
  return(list(cor_mat = cor_data, cluster_labels = clust_res))
  
}



