library(dplyr)
library(ggplot2)
library(readr)
library(magrittr)
library(mrfDepthLight)
library(class)
source('/home/ad/Desktop/KUL Course Material/Robust Statistics/Project/DO_code_clean.R')
source('/home/ad/Desktop/KUL Course Material/Robust Statistics/Project/classificationcode.R')

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

## specifying the stuff for the simulation
n_obs = dim(hubert_data)[1]
n_sim = 10

data_results = data.frame(
  missclassSdo = rep(NA, n_sim), missclassPercSdo = rep(NA, n_sim), Sdo_k = rep(NA, n_sim),
  missclassBd = rep(NA, n_sim), missclassPercBd = rep(NA, n_sim), Bd_k = rep(NA, n_sim),
  missclassAo = rep(NA, n_sim), missclassPercAo = rep(NA, n_sim), Ao_k = rep(NA, n_sim),
  missclassDo = rep(NA, n_sim), missclassPercDo = rep(NA, n_sim), Do_k = rep(NA, n_sim),
  missclassKnn = rep(NA, n_sim), missclassPercKnn = rep(NA, n_sim), Knn_k = rep(NA, n_sim)
)

# function which transforms the data
dist_transform = function(transform.method = outlyingness) {
  
  # converting the string to name object
  # dist_method = as.name(method)
  
  # distance from unemployed to unemployed
  dist00 <- data.frame(
    distance = transform.method(hubert_data[grp0, -12], hubert_data[grp0, -12])$outlyingnessZ,
                       class = rep(0, nrow(hubert_data[grp0,]))
    )
  
  # distance from employed to unemployed
  dist01 <- data.frame(
    distance = transform.method(hubert_data[grp0, -12], hubert_data[grp1, -12])$outlyingnessZ,
                       class = rep(1, nrow(hubert_data[grp1,])))
  
  # distance from unemployed to employed
  dist10 <- data.frame(
    distance = transform.method(hubert_data[grp1, -12], hubert_data[grp0, -12])$outlyingnessZ,
                       class = rep(0, nrow(hubert_data[grp0,])))
  
  # distance from employed to employed
  dist11 <- data.frame(
    distance = transform.method(hubert_data[grp1, -12], hubert_data[grp1, -12])$outlyingnessZ,
                       class = rep(1, nrow(hubert_data[grp1,])))
  
  dist0 = rbind(dist00, dist01)
  dist1 = rbind(dist10, dist11)
  
  transformed = data.frame(distG0 = dist0, distG1 = dist1) %>% 
    dplyr::select(-distG1.class) %>% 
    rename(class = distG0.class, distG0 = distG0.distance, distG1 = distG1.distance) %>%
    mutate(class = as.factor(class)) 
  
  return(transformed)
}

# getting the indices for the groups
grp0 = which(hubert_data$employment == 0)
grp1 = which(hubert_data$employment == 1)

for (i in 1:n_sim) {  
  
  print(paste('Simulation ', i))
  
  # training and test sets 50%
  #train_ix = sample(length(hubert_data$employment), 0.8*length(hubert_data$employment))
  
  # grp0 = which(hubert_data$employment == 0)
  # grp1 = which(hubert_data$employment == 1)
  
  #################################################################################################
  ####################  SDO #######################################################################
  
  dist_sdo = dist_transform(transform.method = outlyingness)
  
  train_ix = sample(length(hubert_data$employment), 0.8*length(hubert_data$employment))
  
  train = dist_sdo[train_ix,]
  test = dist_sdo[-train_ix,]
  
  #knn 
  
  accuracy <- rep(0, 10)
  k <- 1:10
  
  for (x in k) {
    prediction <- knn(train[,c(1,3)], test[,c(1,3)], train[,2], k = x)
    accuracy[x] <- mean(prediction == test[, 2])
  }
  # plot(k, accuracy, type = 'b')
  
  k_best = which.max(accuracy)
  
  sdo_class = knn(train[,c(1,3)], test[,c(1,3)], train[,2], k = k_best)
  # summary(sdo_class)
  
  sdo_results = table(test[,2], sdo_class)
  
  # missclassifications
  missclassif_sdo = sdo_results[1, 2] + sdo_results[2, 1]
  #missclassif_sdo/length(test$class)
  
  data_results[i, 'missclassSdo'] = missclassif_sdo
  data_results[i, 'missclassPercSdo'] = missclassif_sdo/length(test$class)    
  data_results[i, 'Sdo_k'] = k_best
  
  #################################################################################################
  ####################  AO  #######################################################################
  
  dist_ao = dist_transform(transform.method = adjOutlyingness)
  
  train = dist_ao[train_ix,]
  test = dist_ao[-train_ix,]
  
  accuracy <- rep(0, 10)
  k <- 1:10
  
  for (x in k) {
    prediction <- knn(train[,c(1,3)], test[,c(1,3)], train[,2], k = x)
    accuracy[x] <- mean(prediction == test[, 2])
  }
  # plot(k, accuracy, type = 'b')
  
  k_best = which.max(accuracy)
  
  ao_class = knn(train[,c(1,3)], test[,c(1,3)], train[,2], k = k_best)
  # summary(sdo_class)
  
  ao_results = table(test[,2], ao_class)
  
  # missclassifications
  missclassif_ao = ao_results[1, 2] + ao_results[2, 1]
  #missclassif_ao/length(test$class)
  
  data_results[i, 'missclassAo'] = missclassif_ao
  data_results[i, 'missclassPercAo'] = (missclassif_ao)/length(test$class)    
  data_results[i, 'Ao_k'] = k_best
  
  #################################################################################################
  ####################  bag distance  #############################################################
  
  #dist_bd = dist_transform(transform.method = bagdistance, element = bagdistance)
  # this failed so doing it manually
  
  bd00 <- data.frame(
    distance = bagdistance(hubert_data[grp0, -12], hubert_data[grp0, -12])$bagdistance, class = 0)
  
  bd01 <- data.frame(
    distance = bagdistance(hubert_data[grp0, -12], hubert_data[grp1, -12])$bagdistance, class = 1)
  
  bd10 <- data.frame(
    distance = bagdistance(hubert_data[grp1, -12], hubert_data[grp0, -12])$bagdistance, class = 0)
  
  bd11 <- data.frame(
    distance = bagdistance(hubert_data[grp1, -12], hubert_data[grp1, -12])$bagdistance, class = 1)
  
  bd0 = rbind(bd00, bd01)
  bd1 = rbind(bd10, bd11)
  
  dist_bd = data.frame(distG0 = bd0, distG1 = bd1) %>% 
    dplyr::select(-distG1.class) %>% 
    rename(class = distG0.class, distG0 = distG0.distance, distG1 = distG1.distance) %>%
    mutate(class = as.factor(class)) 
  
  train = dist_bd[train_ix,]
  test = dist_bd[-train_ix,]
  
  accuracy <- rep(0, 10)
  k <- 1:10
  
  for (x in k) {
    prediction <- knn(train[,c(1,3)], test[,c(1,3)], train[,2], k = x)
    accuracy[x] <- mean(prediction == test[, 2])
  }
  
  # plot(k, accuracy, type = 'b')
  
  k_best = which.max(accuracy)
  
  bd_class = knn(train[,c(1,3)], test[,c(1,3)], train[,2], k = k_best)
  # summary(sdo_class)
  
  bd_results = table(test[,2], bd_class)
  
  # missclassifications
  missclassif_bd = bd_results[1, 2] + bd_results[2, 1]
  #missclassif_bd/length(test$class)
  
  data_results[i, 'missclassBd'] = missclassif_bd
  data_results[i, 'missclassPercBd'] = (missclassif_bd)/length(test$class)    
  data_results[i, 'Bd_k'] = k_best
  
  #################################################################################################
  ####################  DO   #######################################################################
  
  #dist_do = dist_transform(transform.method = DOclassif)
  
  do00 <- data.frame(
    distance = DOclassif(data.matrix(hubert_data[grp0, -12]), data.matrix(hubert_data[grp0, -12])), 
    class = 0)
  
  do01 <- data.frame(
    distance = DOclassif(data.matrix(hubert_data[grp0, -12]), data.matrix(hubert_data[grp1, -12])), 
    class = 1)
  
  do10 <- data.frame(
    distance = DOclassif(data.matrix(hubert_data[grp1, -12]), data.matrix(hubert_data[grp0, -12])), 
    class = 0)
  
  do11 <- data.frame(
    distance = DOclassif(data.matrix(hubert_data[grp1, -12]), data.matrix(hubert_data[grp1, -12])), 
    class = 1)
  
  do0 = rbind(do00, do01)
  do1 = rbind(do10, do11)
  
  dist_do = data.frame(distG0 = do0, distG1 = do1) %>% 
    dplyr::select(-distG1.class) %>% 
    rename(class = distG0.class, distG0 = distG0.distance, distG1 = distG1.distance) %>%
    mutate(class = as.factor(class)) 
  
  train = dist_do[train_ix,]
  test = dist_do[-train_ix,]
  
  accuracy <- rep(0, 10)
  k <- 1:10
  
  for (x in k) {
    prediction <- knn(train[,c(1,3)], test[,c(1,3)], train[,2], k = x)
    accuracy[x] <- mean(prediction == test[, 2])
  }
  
  # plot(k, accuracy, type = 'b')
  
  k_best = which.max(accuracy)
  
  do_class = knn(train[,c(1,3)], test[,c(1,3)], train[,2], k = k_best)
  # summary(sdo_class)
  
  do_results = table(test[,2], do_class)
  
  # missclassifications
  missclassif_do = do_results[1, 2] + do_results[2, 1]
  #missclassif_do/length(test$class)
  
  data_results[i, 'missclassDo'] = missclassif_do
  data_results[i, 'missclassPercDo'] = (missclassif_do)/length(test$class)    
  data_results[i, 'Do_k'] = k_best
  
  #################################################################################################
  ###################### KNN  #####################################################################
  
  train = hubert_data[train_ix,]
  test = hubert_data[-train_ix,]
  
  accuracy <- rep(0, 10)
  k <- 1:10
  
  for (x in k) {
    prediction <- knn(train[,-12], test[,-12], train[,12]$employment, k = x)
    accuracy[x] <- mean(prediction == test[,12]$employment)
  }
  
  # plot(k, accuracy, type = 'b')
  
  k_best = which.max(accuracy)
  
  knn_class = knn(train[,-12], test[,-12], train[,12]$employment, k = k_best)
  # summary(sdo_class)
  
  knn_results = table(test[,12]$employment, knn_class)
  
  # missclassifications
  missclassif_knn = knn_results[1, 2] + knn_results[2, 1]
  #missclassif_knn/length(test$employment)
  
  data_results[i, 'missclassKnn'] = missclassif_knn
  data_results[i, 'missclassPercKnn'] = (missclassif_knn)/length(test$employment)   
  data_results[i, 'Knn_k'] = k_best
  
}

beepr::beep(3)
#################################################################################################

results = data_results %>% dplyr::select(contains("perc"))
label_vec = c("SDO DistSpace", "BD DistSpace", "AO DistSpace", "DO DistSpace", "kNN")

par(mfrow = c(1,5))

for (i in 1:5) {
  boxplot(results[[i]])#, xlab = label_vec[i], ylab = "% Misclassification")
  mtext(label_vec[i], side = 1, line = 1) #xlab
  mtext("% Misclassification", side = 2, line = 2.4)
}

par(mfrow = c(1,1))

title(main = "Misclassification % Using DistSpace Method", line = 2.3)

# boxplot(
#   data.frame(data_results$missclassPercSdo, data_results$missclassPercAo, data_results$missclassPercBd, data_results$missclassPercDo, data_results$missclassPercKnn),
#  , main="% Missclassification")