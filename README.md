# Microbiome Network Analysis
 
## Tutorial

The following script is used to fit a Weighted Stochastic Infinite Block Model (WSIBM) to analyze real valued weighted adjacency matrices. We consider data in the form of a correlation matrix obtained by performing central log ratio transformation (CLR) on microbiome count data (taxonomic abundance table). WSIBM performs model-based clustering on the correlation matrix in order to detect the latent community structure as well as estimate the distribution of number of communities via Dirichlet Process stick-breaking. We further evaluate the posterior summaries via Posterior Probability Matrix (PPM) which enables us to solve the label switching problem of clustering. 

## Required Packages

```r
library(Rcpp)
library(tidyverse)
library(mcclust)
library(MASS)
library(ggpubr)
```

## Required Functions

```r
sourceCpp("scripts/SBM_cpp_v3.4.cpp") # CPP function for WSBM
source("scripts/functions.R")
```

## Simulation Study

The following are the data components and model parameters:

-   `n`: Number of nodes (taxa) in the network 
-   `W`: Weighted adjacency matrix (correlation matrix with diagonals set to 0)
-   `K`: Number of communities
-   `z`: Latent community membership vector
-   `mu`: Matrix of weighted block means
-   `var`: Matrix of weighted block variance
-   `Kmax`: Upper bound over K
-   `eta0`: Dirichlet process concentration parameter 


``` r
# Simulating the data 
n <- 100
K <- 4
mu_true <- diag(c(-0.9, 0.7, 0.4, -0.5))
mu_true[1, 2] <- 0.3
mu_true[2, 3] <- -0.3
mu_true[1, 3] <- mu_true[1, 4] <- mu_true[2, 4] <- mu_true[3, 4] <- 0.001

cor.sim <- WSBM_sim(n = 100, K = 4, mu_true = mu_true, var_true = matrix(0.1, K, K),
                    seed = 1) # True correlation matrix
cor.mat <- cor.sim$cor.mat

# Visualizing the true correlation matrix

cor.mat.melt <- reshape2::melt(cor.mat)
colnames(cor.mat.melt) <- c("Row", "Col", "Correlation")

g1 <- ggplot(data = cor.mat.melt, aes(x = Row, y = Col)) +
  geom_tile(aes(fill = Correlation), size = 0.25) +
  theme_light() +
  scale_fill_gradient2(low = "pink", mid = "white", high = "green", limits = c(-1, 1)) +
  labs(title = "Real Data")

# Permutating the matrix to mimic real data

sample_id <- sample.int(n = n, size = n, replace = F)
cor.mat.temp <- cor.mat[sample_id, sample_id]
rownames(cor.mat.temp) <- colnames(cor.mat.temp) <- 1:n

cor.mat.melt <- reshape2::melt(cor.mat.temp)
colnames(cor.mat.melt) <- c("Row", "Col", "Correlation")

g2 <- ggplot(data = cor.mat.melt, aes(x = Row, y = Col)) +
  geom_tile(aes(fill = Correlation), size = 0.25) +
  theme_light() +
  scale_fill_gradient2(low = "pink", mid = "white", high = "green", limits = c(-1, 1)) +
  labs(title = "Permutated Data")

ggarrange(g1, g2)
```
![alt text](https://github.com/tejasvbedi95/MicrobiomeNetworkAnalysis/blob/b86e8f316be064c700355a4ff974a93973bb3533/Results/sim_plot_1.png)

``` r
# Fit the permutated data to WSBM function

res <- auto_WSBM(cor.mat.temp, K_max = 20, eta0 = 1, store = T) 

# Obtaining z_ppm (clustering result)

diag(res$ppm_store) <- 500
clust_res <- minbinder(res$ppm_store/500, method = "comp")$cl

names(clust_res) <- sample_id
clust_res_arranged <- clust_res[order(as.numeric(names(clust_res)))]

W_data_temp <- cbind(order = clust_res, cor.mat.temp)
W_data_order <- W_data_temp[, -1][order(W_data_temp[, 1]), order(W_data_temp[, 1])]
rownames(W_data_order) <- colnames(W_data_order) <- 1:n

cor.res.melt <- reshape2::melt(W_data_order)
colnames(cor.res.melt) <- c("Row", "Col", "Correlation")

g3 <- ggplot(data = cor.res.melt, aes(x = Row, y = Col)) +
  geom_tile(aes(fill = Correlation), size = 0.25) +
  theme_light() +
  scale_fill_gradient2(low = "pink", mid = "white", high = "green", limits = c(-1, 1)) +
  labs(title = "Clustering Result")

g3 # Plot of clustering result
```
![alt text](https://github.com/tejasvbedi95/MicrobiomeNetworkAnalysis/blob/b86e8f316be064c700355a4ff974a93973bb3533/Results/sim_plot_2.png)

``` r
# External cluster validation
z_true <- cor.sim$z_true

ARI(z_true, clust_res_arranged) # = 1, perfect agreement
NMI(z_true, clust_res_arranged) # = 1, perfect agreement
NVI(z_true, clust_res_arranged) # = 0, perfect agreement

```

## Real Data Analysis
``` r
# Loading Count Table
data <- read.csv("Data/metaphlan_qc.csv", header = T)[, -c(1:2)]

# Calculating spearman correlation on CLR transformed data
cor_data <- count_to_cor(data, method = "spearman")$CLR 

## write.csv(cor_data, file = "Results/CLR_corr.csv", row.names = T)

# Visualizing the data

cor.mat.melt <- reshape2::melt(cor_data)
colnames(cor.mat.melt) <- c("Row", "Col", "Correlation")

g4 <- ggplot(data = cor.mat.melt, aes(x = Row, y = Col, fill = Correlation)) +
  geom_tile(color = "white", size = 0.25) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, size = 5),
        axis.text.y = element_text(size = 5)) +
  labs(x="", y="") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "green", limits = c(-1, 1)) +
  labs(title = "Raw Correlation Matrix")
  

# Fitting WSBM model

res <- auto_WSBM(cor_data, 20, 1, T)

diag(res$ppm_store) <- 500
clust_res <- minbinder(res$ppm_store/500, method = "comp")$cl

ppm_mat <- res$ppm_store/500
colnames(ppm_mat) <- rownames(ppm_mat) <- colnames(cor_data)
## write.csv(ppm_mat, file = "Results/CLR_ppm.csv")

W_data_temp <- cbind(order = clust_res, cor_data)
W_data_order <- W_data_temp[, -1][order(W_data_temp[, 1]), order(W_data_temp[, 1])]

sort_vec <- cumsum(tapply(sort(clust_res), as.factor(sort(clust_res)), length)) + 0.5

cor.res.melt <- reshape2::melt(W_data_order)
colnames(cor.res.melt) <- c("Row", "Col", "Correlation")

g5 <- ggplot(data = cor.res.melt, aes(x = Row, y = Col, fill = Correlation)) +
  geom_tile(color = "white", size = 0.25) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, size = 5),
        axis.text.y = element_text(size = 4.5)) +
  labs(x="", y="") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "green", limits = c(-1, 1)) +
  geom_hline(yintercept = sort_vec, size = 0.7, linetype = "dashed", col = "red", alpha = 0.4) +
  geom_vline(xintercept = sort_vec, size = 0.7, linetype = "dashed", col = "red", alpha = 0.4) +
  labs(title = "Clustered Correlation Matrix")

ggarrange(g4, g5)

```
![alt text](https://github.com/tejasvbedi95/MicrobiomeNetworkAnalysis/blob/41c32bb3b7c1deca1489dc75137cbe648c73b835/Results/res_plot_1.png)

###Using WSBM wrapper function

#### This function enables user to input the raw counts data and obtain the above results

``` r
data <- read.csv("Data/metaphlan_qc.csv", header = T)[, -c(1:2)]

res <- WSBM_wrapper(data, K = "auto", cor = "spearman", transform = T)
clust_res <- res$cluster_labels

```









