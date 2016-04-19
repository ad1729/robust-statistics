library(dplyr)
library(ggplot2)
library(readr)
library(magrittr)
library(mrfDepthLight)
library(class)

## Function for simulation

dist_space = function(data = hubert_data, method = "outlyingess") {
  
  # creating the training and test sets
  train_ix = sample(length(hubert_data$employment), 0.8*length(hubert_data$employment))
  train = hubert_data[train_ix,]
  test = hubert_data[-train_ix,]
  
  grp0 = which(train$employment == 0)
  grp1 = which(train$employment == 1)
  
  # converting the string to name object
  dist_method = as.name(method)
  
  # distance from unemployed to unemployed
  dist00 <- data.frame(distance = dist_method(hubert_data[grp0, -12], hubert_data[grp0, -12])$outlyingnessZ,
                      class = 0)
  
  # distance from employed to unemployed
  dist01 <- data.frame(distance = dist_method(hubert_data[grp0, -12], hubert_data[grp1, -12])$outlyingnessZ,
                       class = 1)
  
  # distance from unemployed to employed
  dist10 <- data.frame(distance = dist_method(hubert_data[grp1, -12], hubert_data[grp0, -12])$outlyingnessZ,
                      class = 0)
  
  # distance from employed to employed
  dist11 <- data.frame(distance = dist_method(hubert_data[grp1, -12], hubert_data[grp1, -12])$outlyingnessZ,
                      class = 1)
  
  dist0 = rbind(dist00, dist01)
  dist1 = rbind(dist10, dist11)
  
  dist_main = data.frame(distG0 = dist0, distG1 = dist1) %>% 
    dplyr::select(-distG1.class) %>% 
    rename(class = distG0.class, distG0 = distG0.distance, distG1 = distG1.distance) %>%
    mutate(class = as.factor(class))
  
  ## knn on the computed distances
  accuracy <- rep(0, 10)
  k <- 1:10
  
  # for(x in k){
  #   prediction <- knn(train[,-12], test[,-12], train[,12]$employment, k = x) #employment status is column 12
  #   accuracy[x] <- mean(prediction == test[,12]$employment)
  # }
  # 
  # k_best = which.max(accuracy)
  # 
  # knn_class = knn(train[,-12], test[,-12], train[,12]$employment, k = k_best) #employment status is column 12
  # 
  # knn_results = table(test[, 12]$employment, knn_class)
  # 
  # (knn_results[1, 2] + knn_results[2, 1])/length(test$employment) #missclassification %
  # 
  # # miss classifications
  # missclassif_knn = knn_results[1, 2] + knn_results[2, 1]
  # sim_results[i, 'missclassKnn'] = missclassif_knn
  # sim_results[i, 'missclassPercKnn'] = (missclassif_knn)/(n_obs*groups - length(train_ix))    
  # sim_results[i, 'Knn_k'] = k_best
  
  
  
  
  
  for (x in k) {
    prediction <- knn(dist_main[train_ix, -12], 
                      dist_main[-train_ix, -12],
                      dist_main[train_ix, 3], k = x)
    accuracy[x] <- mean(prediction == dist_main[-train_ix, 3])
  }
  # plot(k, accuracy, type = 'b')
  
  # picking the neighbourhood value of k with best accuracy
  k_best = which.max(accuracy)
  
  # rerunning the knn classifier with k = k_best
  sdo_class = knn(dist_main[train_ix, c(1, 2)], 
                  dist_main[-train_ix, c(1, 2)], 
                  dist_main[train_ix, 3],
                  k = k_best)
  # summary(sdo_class)
  
  sdo_results = table(dist_main[-train_ix, 3], sdo_class)
  
  # miss classifications
  missclassif = sdo_results[1, 2] + sdo_results[2, 1]
  # missclassif
  
  sim_results[i, 'missclassSdo'] = missclassif
  sim_results[i, 'missclassPercSdo'] = (missclassif)/(n_obs*groups - length(train_ix))    
  sim_results[i, 'Sdo_k'] = k_best
}

## Main code
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