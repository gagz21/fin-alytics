# comments
library(ggplot2)
clt_sim <- function(n, source=NULL, param1=NULL, param2=NULL) {
  r <- 10000  # Number of replications/samples - DO NOT ADJUST
  # This produces a matrix of observations with  
  # n columns and r rows. Each row is one sample:
  my.samples <- switch(source,
                       "E" = matrix(rexp(n*r,param1),r),
                       "N" = matrix(rnorm(n*r,param1,param2),r),
                       "U" = matrix(runif(n*r,param1,param2),r),
                       "P" = matrix(rpois(n*r,param1),r),
                       "C" = matrix(rcauchy(n*r,param1,param2),r),
                       "B" = matrix(rbinom(n*r,param1,param2),r),
                       "G" = matrix(rgamma(n*r,param1,param2),r),
                       "X" = matrix(rchisq(n*r,param1),r),
                       "T" = matrix(rt(n*r,param1),r))
  all.sample.sums <- apply(my.samples,1,sum)
  all.sample.means <- apply(my.samples,1,mean)   
  all.sample.vars <- apply(my.samples,1,var) 
  par(mfrow=c(2,2))
  hist(my.samples[1,],col="gray",main="Distribution of One Sample")
  hist(all.sample.sums,col="gray",main="Sampling Distribution \n of
       the Sum")
  hist(all.sample.means,col="gray",main="Sampling Distribution \n of the Mean")
  hist(all.sample.vars,col="gray",main="Sampling Distribution \n of
       the Variance")
}
clt_sim(100,source="E",param1=1)


require(pdfetch)
require(xts)
require(zoo)
EIUIR <- pdfetch_BLS(c("EIUIR", "EIUIR100"), 2000, 2016) # start and end years
head(EIUIR)
xmprice <- na.omit(EIUIR) # to clean up any missing data
xmprice.r <- as.zoo(na.omit((diff(log(xmprice[, 1])))))
head(xmprice.r)

write.csv(EIUIR, "data/EIUIR-test.csv")

EIUIR_test <- read.csv("data/EIUIR-test.csv")
head(EIUIR_test)

str(EIUIR_test)

library(ggplot2)
library(dplyr)
library(moments)
n <- 100
r <- 100000
sample_sim <- matrix(rexp(n*r, 1), ncol = n)
cor_rw <- function (variable){
  var_lagged <- na.omit(lag(variable))
  variable <- variable[-1]
  return(cor(variable, var_lagged))
}
x <- rexp(n*r, 1)
y <- matrix(1 + 2 * x, nrow = r) 
mean_sim <- apply(y, 1, mean)
sd_sim <- apply(y, 1, sd)
skew_sim <- apply(y, 1, skewness)
kurt_sim <- apply(y, 1, kurtosis)
IQR_sim <- apply(y, 1, IQR)
cor_sim <- apply(y, 1, cor_rw)
cor_sim <- sapply(1:nrow(y), cor_rw)
sim_df <- data.frame(y = y[,1], mean = mean_sim, sd = sd_sim, skewness = skew_sim, kurtosis = kurt_sim, IQR = IQR_sim, acf_1 = cor_sim)
data_moments <- function(data){
  library(moments)
  library(matrixStats)
  mean.r <- colMeans(data)
  median.r <- colMedians(data)
  sd.r <- colSds(data)
  IQR.r <- colIQRs(data)
  skewness.r <- skewness(data)
  kurtosis.r <- kurtosis(data)
  result <- data.frame(mean = mean.r, median = median.r, std_dev = sd.r, IQR = IQR.r, skewness = skewness.r, kurtosis = kurtosis.r)
  return(result)
}
# Run data_moments()
answer <- data_moments(as.matrix(sim_df))
# Build pretty table
answer <- round(answer, 4)
knitr::kable(answer)
knitr::kable(summary(sim_df))
ggplot(sim_df , aes(x = y)) + geom_histogram()
ggplot(sim_df , aes(x = mean)) + geom_histogram()
ggplot(sim_df , aes(x = sd)) + geom_histogram()
ggplot(sim_df , aes(x = skewness)) + geom_histogram()
ggplot(sim_df , aes(x = kurtosis)) + geom_histogram()
ggplot(sim_df , aes(x = IQR)) + geom_histogram()
ggplot(sim_df , aes(x = acf_1)) + geom_histogram()

summary(cor_sim)
skewness(cor_sim)
summary(sim_df)
quantile(sim_df$kurtosis, 0.025)
quantile(sim_df$kurtosis, 0.9725)
mean(sim_df$kurtosis <= 3)
ecdf(sim_df$kurtosis)(3)

set.seed(1016)
quantile_95_sim <- replicate(10000, quantile(sample(exrates.r[, 1], size = 100, replace = FALSE), 0.95))
quantile_99_sim <- replicate(10000, quantile(sample(exrates.r[, 1], size = 100, replace = FALSE), 0.99))
quantile_df <- data.frame(quantile_95_sim = quantile_95_sim, quantile_99_sim = quantile_99_sim)
data_moments(as.matrix(quantile_df))

