rm(list = ls())
#setwd("/media/ornela/caa36dee-ee54-44e7-8fa8-872f75e66246/ornela/Dropbox/MSC Statistics/4th semester/RobustStatistics/project/code")

library(mrfDepthLight)
library(ggplot2)
library(class)
source('DO_code_clean.R')
source('covariance.r')

n_obs = 500
groups = 2

# generate data
# group 1: bivariate normal distribution, mean = 5, sd = 1, corr  = -0.5
# group 2: bivariate normal distribution, mean = 8, sd = 1, corr  = 0.7
# group 3: bivariate exponential distribution, lambda = 1, corr = 0

g1 = mvrnorm(n = n_obs, mu = rep(4, 2), Sigma = create_covariance(p = 2, rho = -0.5)) 
# g2 = mvrnorm(n = n_obs, mu = rep(8, 2), Sigma = create_covariance(p = 2, rho = 0.7)) 
g3_1 = rexp(n = n_obs, rate = 1)
g3_2 = rexp(n = n_obs, rate = 1)
# g3 = cbind(g3_1, g3_2)
g2 = cbind(g3_1, g3_2)

my_data = data.frame(x1 = g1[,1], x2 = g1[,2], class = 1)
my_data2 = data.frame(x1 = g2[,1], x2 = g2[,2], class = 2)
# my_data3 = data.frame(x1 = g3[,1], x2 = g3[,2], class = 3)

my_data = rbind(my_data, my_data2)#, my_data3)
# remove(my_data2, my_data3)

p <- ggplot(data = my_data, aes(x = x1, y = x2)) + geom_point()
p <- p + aes(colour = factor(class))
p

class1_ix = which(my_data$class == 1)
class2_ix = which(my_data$class == 2)

train_ix1 = sample(class1_ix, 200)
train_ix2 = sample(class2_ix, 200)
train_ix = c(train_ix1, train_ix2)

#################################################################################################
####################  SDO #######################################################################

sdo1 <- outlyingness(my_data[class1_ix, c(1, 2)], my_data[, c(1, 2)])
sdo2 <- outlyingness(my_data[class2_ix, c(1, 2)], my_data[, c(1, 2)])
# sdo3 <- outlyingness(my_data[which(my_data$class == 3), c(1, 2)], my_data[, c(1, 2)])

dist_sdo = data.frame(sdoG1 = sdo1$outlyingnessZ, sdoG2 = sdo2$outlyingnessZ)#, sdoG3 = sdo3$outlyingnessZ)

dist_sdo$class = c(rep(1, 500), rep(2, 500))

g = ggplot(data = dist_sdo, aes(x = sdoG1, y = sdoG2)) + geom_point()
g = g + aes(colour = factor(class))
g

sdo_class = knn(dist_sdo[train_ix, c(1, 2)], dist_sdo[-train_ix, c(1, 2)], dist_sdo[train_ix, 3],
                k = 5)
summary(sdo_class)

sdo_results = table(dist_sdo[-train_ix, 3], sdo_class)

# miss classifications
missclassif = sdo_results[1, 2] + sdo_results[2, 1]
missclassif
length(missclassif)/(n_obs*groups - length(train_ix))

#################################################################################################
####################  AO  #######################################################################

ao1 <- adjOutlyingness(my_data[class1_ix, c(1, 2)], my_data[, c(1, 2)])
ao2 <- adjOutlyingness(my_data[class2_ix, c(1, 2)], my_data[, c(1, 2)])

dist_ao = data.frame(aoG1 = ao1$outlyingnessZ, aoG2 = ao2$outlyingnessZ)#, sdoG3 = sdo3$outlyingnessZ)

dist_ao$class = c(rep(1, 500), rep(2, 500))

g = ggplot(data = dist_ao, aes(x = aoG1, y = aoG2)) + geom_point()
g = g + aes(colour = factor(class))
g

ao_class = knn(dist_ao[train_ix, c(1, 2)], dist_ao[-train_ix, c(1, 2)], dist_ao[train_ix, 3],
                k = 5)
summary(ao_class)

ao_results = table(dist_ao[-train_ix, 3], ao_class)

# miss classifications
missclassif_ao = ao_results[1, 2] + ao_results[2, 1]
missclassif_ao
length(missclassif_ao)/(n_obs*groups - length(train_ix))

accuracy <- rep(0, 10)
k <- 1:10

for(x in k){
    prediction <- knn(dist_ao[train_ix, c(1, 2)], dist_ao[-train_ix, c(1, 2)],
                      dist_ao[train_ix, 3], k = x)
    accuracy[x] <- mean(prediction == dist_ao[-train_ix, 3])
}
plot(k, accuracy, type = 'b')

which.min(accuracy)

#################################################################################################
####################  bag distance  #######################################################################

bd1 <- bagdistance(my_data[class1_ix, c(1, 2)], my_data[, c(1, 2)])
bd2 <- bagdistance(my_data[class2_ix, c(1, 2)], my_data[, c(1, 2)])

dist_bd = data.frame(bdG1 = bd1$bagdistance, bdG2 = bd2$bagdistance)#, sdoG3 = sdo3$outlyingnessZ)

dist_bd$class = c(rep(1, 500), rep(2, 500))

g = ggplot(data = dist_bd, aes(x = bdG1, y = bdG2)) + geom_point()
g = g + aes(colour = factor(class))
g

bd_class = knn(dist_bd[train_ix, c(1, 2)], dist_bd[-train_ix, c(1, 2)], dist_bd[train_ix, 3],
               k = 5)
summary(bd_class)

bd_results = table(dist_bd[-train_ix, 3], bd_class)

# miss classifications
missclassif_bd = bd_results[1, 2] + bd_results[2, 1]
missclassif_bd
length(missclassif_bd)/(n_obs*groups - length(train_ix))

#################################################################################################
####################  DO #######################################################################
