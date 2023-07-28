install.packages("corrr")
library('corrr')
library(tidyverse)
install.packages('glmnet')

install.packages("ggcorrplot")
library(ggcorrplot)

install.packages("FactoMineR")
library("FactoMineR")


# Library required for fviz_cluster function
install.packages("factoextra")
library(factoextra)
library(dplyr)

paths = as.data.frame(pathway.scores)
pca = as.data.frame(pc_scores)
survive = as.data.frame(Y)


library(glmnet)

#perform k-fold cross-validation to find optimal lambda value
cv_model <- cv.glmnet(pc_scores, Y, alpha = 0.9)

#find optimal lambda value that minimizes test MSE
best_lambda <- cv_model$lambda.min
best_lambda

plot(cv_model) 

best_model <- glmnet(pc_scores, Y, alpha = 0.9, lambda = best_lambda)
coef(best_model)


a <- matrix(rep(-0.5, 143), nrow=1, ncol=143)


predict(best_model, s = best_lambda, newx = a)

y_predicted <- predict(best_model, s = best_lambda, newx = pc_scores)


#find SST and SSE
sst <- sum((Y - mean(Y))^2)
sse <- sum((y_predicted - Y)^2)

#find R-Squared
rsq <- 1 - sse/sst
rsq










# 
# path_means = paths %>% 
#   rowMeans()
# data1 <- data.frame(Y, path_means)
# mod = lm(path_means ~ Y)
# summary(mod)
# data1 %>% 
#   ggplot(mapping = aes(x = path_means, y = Y)) + geom_point() 
# 
# km <- kmeans(paths, centers = 3, nstart = 25)
# 
# # Visualize the clusters
# fviz_cluster(km, data = paths)
# 
# data2 <- pca                                           # Duplicate example data
# data2 <- tibble::rownames_to_column(data2, "row_names") # Apply rownames_to_column
# data2                                             
# 
# path_scores = data2 %>% 
#   pivot_longer(cols=colnames(select(data2, -row_names)),
#                names_to='MRI',
#                values_to='value')
# 
# path_scores$MRI = gsub("\\.\\d", "", path_scores$MRI)
# 
# what = path_scores %>% 
#   group_by(row_names, MRI) %>% 
#   summarise(avg = mean(value)) %>% View()
# 
# df = data.frame(what, Y) %>% View()
#   ggplot(mapping = aes(x = avg, y = Y)) + geom_point()
# 
# patient_means = pca %>% 
#   rowMeans()
# 
# data<-data.frame(Y,patient_means)
# 
# model = glm(patient_means ~ Y)
# summary(model)
# 
# data %>% 
#   ggplot(mapping = aes(x = patient_means, y = Y)) + geom_point() + geom_smooth(method = 'lm')
# 
# km2 <- kmeans(pca, centers = 3, nstart = 25)
# 
# # Visualize the clusters
# fviz_cluster(km2, data = pca)





