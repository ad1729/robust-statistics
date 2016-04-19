library(dplyr)
library(ggplot2)
library(readr)
library(magrittr)
library(faraway)
data("motorins")
data("pima")

library(mrfDepthLight)
library(class)
library(MASS)

## Tumor Classification data
m = read_csv("/home/ad/Desktop/KUL Course Material/Robust Statistics/Project/Data/TutorialRCC366/testData_366.csv") %>% dplyr::select(-`NA`) %>% mutate(ClearCell = as.factor(ClearCell), event = as.factor(event))
glimpse(m)
pairs(formula = ~., data = m[,1:8], col = m$ClearCell)
table(m$ClearCell) # G = 2, p = 8 (8 markers)

## Motorins data from package faraway
str(motorins)
pairs(formula = ~., data = motorins[,5:7], col = as.factor(motorins$Bonus))
table(motorins$Bonus)
pairs(formula = ~., data = motorins[,5:7], col = motorins$Kilometres)

## PIMA Indians dataset

pima_dat = pima
summary(pima)
pima_dat$test = as.factor(pima_dat$test)
pairs(formula = ~., data = pima, col = as.factor(pima$test))

## Hubert 2010 dataset (belgian national household 2005 survey)

vars = c("clothing", "alcoholic_drinks", "durable_consumer_goods", "energy", "food", "health", "leisure", "nonalcoholic_drinks", "transport", "housing", "income")

employed = read_csv("/home/ad/Desktop/KUL Course Material/Robust Statistics/Project/Data/Hubert2010_emp.csv", col_names = FALSE)

unemployed = read_csv("/home/ad/Desktop/KUL Course Material/Robust Statistics/Project/Data/Hubert2010_unemp.csv", col_names = FALSE)

names(employed) = vars
names(unemployed) = vars

employed$employment = 1

unemployed$employment = 0

hubert_data = rbind(employed, unemployed) %>% mutate(employment = as.factor(employment))

glimpse(hubert_data)

pairs(formula = ~., data = hubert_data[,-12], col = hubert_data$employment)

#####-------------------------------------------------------------------
## Code taken from the simulation file v4simulationRobustStats.R
#####-------------------------------------------------------------------

n_obs = dim(hubert_data)[1]
n_sim = 2

data_results = data.frame(
  missclassSdo = rep(NA, n_sim), missclassPercSdo = rep(NA, n_sim), Sdo_k = rep(NA, n_sim),
  missclassBd = rep(NA, n_sim), missclassPercBd = rep(NA, n_sim), Bd_k = rep(NA, n_sim),
  missclassAo = rep(NA, n_sim), missclassPercAo = rep(NA, n_sim), Ao_k = rep(NA, n_sim),
  missclassDo = rep(NA, n_sim), missclassPercDo = rep(NA, n_sim), Do_k = rep(NA, n_sim),
  missclassKnn = rep(NA, n_sim), missclassPercKnn = rep(NA, n_sim), Knn_k = rep(NA, n_sim)
)

for (i in 1:n_sim) {  
  
  print(paste('Simulation ', i))
  
  # divide into train (40%) and test set (60%)
  # class1_ix = which(pima_dat$test == 0)
  # class2_ix = which(pima_dat$test == 1)
  # 
  # train_ix1 = base::sample(class1_ix, 0.5 * length(class1_ix))
  # train_ix2 = base::sample(class2_ix, 0.5 * length(class2_ix))
  # train_ix = c(train_ix1, train_ix2)
  
  # training and test sets 50%
  
  train_ix = sample(length(hubert_data$employment), 0.8*length(hubert_data$employment))
  train = hubert_data[train_ix,]
  test = hubert_data[-train_ix,]
  
  grp0 = which(train$employment == 0)
  grp1 = which(train$employment == 1)
  
  #################################################################################################
  ####################  SDO #######################################################################
  
  # distance from unemployed to unemployed
  sdo00 <- data.frame(SDO = outlyingness(hubert_data[grp0, -12], hubert_data[grp0, -12])$outlyingnessZ, 
                      class = 0)
  
  # distance from employed to unemployed
  sdo01 <- data.frame(SDO = outlyingness(hubert_data[grp0, -12], hubert_data[grp1, -12])$outlyingnessZ, 
                      class = 1)
  
  # distance from unemployed to employed
  sdo10 <- data.frame(SDO = outlyingness(hubert_data[grp1, -12], hubert_data[grp0, -12])$outlyingnessZ, 
                      class = 0)
  
  # distance from employed to employed
  sdo11 <- data.frame(SDO = outlyingness(hubert_data[grp1, -12], hubert_data[grp1, -12])$outlyingnessZ, 
                      class = 1)
  
  # sdo12 <- outlyingness(hubert_data[g1, -12], hubert_data[g2, -12])
  # sdo21 <- outlyingness(hubert_data[g2, -12], hubert_data[g1, -12])
  # sdo22 <- outlyingness(hubert_data[g2, -12], hubert_data[g2, -12])
  
  # # SDO from unemployed to unemployed
  # sdo00 <- data.frame(distance = sdo12$outlyingnessX, employment = 0)
  # 
  # # SDO from employed to unemployed
  # sdo10 <- data.frame(distance = sdo12$outlyingnessZ, employment = 1)
  # 
  # # SDO from employed to unemployed
  # sdo01 <- data.frame(distance = sdo21$outlyingnessZ, employment = 1)
  
  sdo0 = rbind(sdo00, sdo01)
  sdo1 = rbind(sdo10, sdo11)
  
  dist_sdo = data.frame(sdoG0 = sdo0, sdoG1 = sdo1) %>% 
    dplyr::select(-sdoG1.class) %>% 
    rename(class = sdoG0.class, sdoG0 = sdoG0.SDO, sdoG1 = sdoG1.SDO) %>%
    mutate(class = as.factor(class))
  
  #dist_sdo$class = c(rep(1, 500), rep(2, 500))
  
  #knn 
  
  accuracy <- rep(0, 10)
  k <- 1:10
  
  for (x in k) {
    prediction <- knn(dist_sdo[train_ix, c(1, 2)], 
                      dist_sdo[-train_ix, c(1, 2)],
                      dist_sdo[train_ix, 3], k = x)
    accuracy[x] <- mean(prediction == dist_sdo[-train_ix, 3])
  }
  # plot(k, accuracy, type = 'b')
  
  k_best = which.max(accuracy)
  
  sdo_class = knn(dist_sdo[train_ix, c(1, 2)], 
                  dist_sdo[-train_ix, c(1, 2)], 
                  dist_sdo[train_ix, 3],
                  k = k_best)
  # summary(sdo_class)
  
  sdo_results = table(dist_sdo[-train_ix, 3], sdo_class)
  
  # miss classifications
  missclassif = sdo_results[1, 2] + sdo_results[2, 1]
  # missclassif
  
  sim_results[i, 'missclassSdo'] = missclassif
  sim_results[i, 'missclassPercSdo'] = (missclassif)/(n_obs*groups - length(train_ix))    
  sim_results[i, 'Sdo_k'] = k_best
  
  #################################################################################################
  ####################  AO  #######################################################################
  
  ao1 <- adjOutlyingness(hubert_data[class1_ix, c(1, 2)], hubert_data[, c(1, 2)])
  ao2 <- adjOutlyingness(hubert_data[class2_ix, c(1, 2)], hubert_data[, c(1, 2)])
  
  dist_ao = data.frame(aoG1 = ao1$outlyingnessZ, aoG2 = ao2$outlyingnessZ)#, sdoG3 = sdo3$outlyingnessZ)
  
  dist_ao$class = c(rep(1, 500), rep(2, 500))
  
  accuracy <- rep(0, 10)
  k <- 1:10
  
  for (x in k) {
    
    prediction <- knn(dist_ao[train_ix, c(1, 2)], 
                      dist_ao[-train_ix, c(1, 2)],
                      dist_ao[train_ix, 3], k = x)
    
    accuracy[x] <- mean(prediction == dist_ao[-train_ix, 3])
  }
  
  k_best = which.max(accuracy)   
  
  ao_class = knn(dist_ao[train_ix, c(1, 2)], 
                 dist_ao[-train_ix, c(1, 2)], 
                 dist_ao[train_ix, 3],
                 k = k_best)
  
  ao_results = table(dist_ao[-train_ix, 3], ao_class)
  
  # miss classifications
  missclassif_ao = ao_results[1, 2] + ao_results[2, 1]
  sim_results[i, 'missclassAo'] = missclassif_ao
  sim_results[i, 'missclassPercAo'] = (missclassif_ao)/(n_obs*groups - length(train_ix))    
  sim_results[i, 'Ao_k'] = k_best
  
  #################################################################################################
  ####################  bag distance  #######################################################################
  
  bd1 <- bagdistance(hubert_data[class1_ix, c(1, 2)], hubert_data[, c(1, 2)])
  bd2 <- bagdistance(hubert_data[class2_ix, c(1, 2)], hubert_data[, c(1, 2)])
  
  dist_bd = data.frame(bdG1 = bd1$bagdistance, bdG2 = bd2$bagdistance)#, sdoG3 = sdo3$outlyingnessZ)
  
  dist_bd$class = c(rep(1, 500), rep(2, 500))
  
  accuracy <- rep(0, 10)
  k <- 1:10
  
  for(x in k){
    prediction <- knn(dist_bd[train_ix, c(1, 2)], dist_bd[-train_ix, c(1, 2)],
                      dist_bd[train_ix, 3], k = x)
    accuracy[x] <- mean(prediction == dist_bd[-train_ix, 3])
  }
  
  k_best = which.max(accuracy)   
  
  bd_class = knn(dist_bd[train_ix, c(1, 2)], dist_bd[-train_ix, c(1, 2)], dist_bd[train_ix, 3],
                 k = k_best)
  
  bd_results = table(dist_bd[-train_ix, 3], bd_class)
  
  # miss classifications
  missclassif_bd = bd_results[1, 2] + bd_results[2, 1]
  sim_results[i, 'missclassBd'] = missclassif_bd
  sim_results[i, 'missclassPercBd'] = (missclassif_bd)/(n_obs*groups - length(train_ix))    
  sim_results[i, 'Bd_k'] = k_best
  
  #################################################################################################
  ####################  DO   #######################################################################
  do1 <- DOclassif(data.matrix(hubert_data[class1_ix, c(1, 2)]), data.matrix(hubert_data[, c(1, 2)]))
  do2 <- DOclassif(data.matrix(hubert_data[class2_ix, c(1, 2)]), data.matrix(hubert_data[, c(1, 2)]))
  
  dist_do = data.frame(doG1 = do1, doG2 = do2)#, sdoG3 = sdo3$outlyingnessZ)
  
  dist_do$class = c(rep(1, 500), rep(2, 500))
  
  accuracy <- rep(0, 10)
  k <- 1:10
  
  for(x in k){
    prediction <- knn(dist_do[train_ix, c(1, 2)], dist_do[-train_ix, c(1, 2)],
                      dist_do[train_ix, 3], k = x)
    accuracy[x] <- mean(prediction == dist_do[-train_ix, 3])
  }
  
  k_best = which.max(accuracy)   
  
  do_class = knn(dist_do[train_ix, c(1, 2)], dist_do[-train_ix, c(1, 2)], dist_do[train_ix, 3],
                 k = k_best)
  
  do_results = table(dist_do[-train_ix, 3], do_class)
  
  # miss classifications
  missclassif_do = do_results[1, 2] + do_results[2, 1]
  sim_results[i, 'missclassDo'] = missclassif_do
  sim_results[i, 'missclassPercDo'] = (missclassif_do)/(n_obs*groups - length(train_ix))    
  sim_results[i, 'Do_k'] = k_best
  
  #################################################################################################
  ###################### KNN  #####################################################################
  accuracy <- rep(0, 10)
  k <- 1:10
  
  for(x in k){
    prediction <- knn(train[,-12], test[,-12], train[,12]$employment, k = x) #employment status is column 12
    accuracy[x] <- mean(prediction == test[,12]$employment)
  }
  
  k_best = which.max(accuracy)
  
  knn_class = knn(train[,-12], test[,-12], train[,12]$employment, k = k_best) #employment status is column 12
  
  knn_results = table(test[, 12]$employment, knn_class)
  
  (knn_results[1, 2] + knn_results[2, 1])/length(test$employment) #missclassification %
  
  # miss classifications
  missclassif_knn = knn_results[1, 2] + knn_results[2, 1]
  sim_results[i, 'missclassKnn'] = missclassif_knn
  sim_results[i, 'missclassPercKnn'] = (missclassif_knn)/(n_obs*groups - length(train_ix))    
  sim_results[i, 'Knn_k'] = k_best
  
}


#################################################################################################


par(mfrow = c(1, 1))
boxplot(sim_results$missclassPercSdo, main = 'SDO DistSpace', ylim = c(0.10, 0.25),
        ylab = '% Missclassification')
# adjbox(sim_results$missclassPercSdo)
boxplot(sim_results$missclassPercAo, main = 'AO DistSpace', ylim = c(0.10, 0.25),
        ylab = '% Missclassification')
boxplot(sim_results$missclassPercBd, main = 'BD DistSpace', ylim = c(0.10, 0.25),
        ylab = '% Missclassification')
boxplot(sim_results$missclassPercDo, main = 'DO DistSpace', ylim = c(0.10, 0.25),
        ylab = '% Missclassification')
boxplot(sim_results$missclassPercKnn, main = 'Knn', ylim = c(0.10, 0.25),
        ylab = '% Missclassification')