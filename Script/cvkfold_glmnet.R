cvKfold_glmnet <- function(k, kfold_sel, PCvar, glmnet_alpha, glmnet_lambda) {
  
  #### Prep data ####
  train <- analysis(kfold_sel$splits[[k]])
  test <- assessment(kfold_sel$splits[[k]])
  
  # Prep train data
  pca <- prcomp(train[, -(1:3)])
  prVar <- pca$sdev^2/sum(pca$sdev^2)
  if(PCvar == 1) {
    pcax <- pca$x
  } else {
    cutoff <- which(cumsum(prVar)>PCvar)[1]
    pcax <- pca$x[, seq_len(cutoff)]
  }
  
  y <- train[, AgeYrs]
  n <- length(y)
  
  # Prep test data
  testx <- predict(object = pca, test[, -(1:3)])
  testx <- testx[, seq_len(cutoff)]
  testy <- test[, AgeYrs]
  
  #### Fit glmnet ####
  set.seed(21)
  fit1=glmnet(as.matrix(train[, -(1:3)]), y, alpha=glmnet_alpha)
  cv.fit1 <- cv.glmnet(as.matrix(train[, -(1:3)]), y, alpha=glmnet_alpha)
  lam <- if(glmnet_lambda == "min") cv.fit1$lambda.min else cv.fit1$lambda.1se
  pred_y_raw <- predict(fit1, newx = as.matrix(test[, -(1:3)]), s=cv.fit1$lambda.1se)  # make predictions
  
  
  #### Fit PC glmnet ####
  set.seed(20)
  fitPca=glmnet(pcax, y, alpha=0.5)
  cv.fitPca <- cv.glmnet(pcax, y, alpha=0.5)
  pred_y_pca <- predict.glmnet(fitPca, newx = testx, s=cv.fitPca$lambda.min)  # make predictions
  
  res <- test[, 1:2]
  res[, AgePred:=pred_y_raw]
  res[, PC:=FALSE]
  tmp <- test[, 1:2]
  tmp[, AgePred:=pred_y_pca]
  tmp[, PC:=TRUE]
  res <- rbind(res, tmp)
  res[, ':='(Model="EN", Data="Test", UncertAge=FALSE)]
  return(res)
}
