# functions for amb algorithms



preparte_layout <- function(tag_ls,sampleWeight_ls,predictor_ls,setID_ls){
  layout <- list()
  layout$tag_ls <- tag_ls
  layout$sampleWeight_ls <- sampleWeight_ls
  layout$predictor_ls <- predictor_ls
  layout$setID_ls <- setID_ls
  
  return(layout)
}




prepare_abm_data <- function(data,taget_ls,sampleWeight_ls,predictor_ls,setID_ls){
  tag_ls <- copy(taget_ls)
  sampleWeight <- copy(sampleWeight_ls)
  num_var_ls <- copy(predictor_ls)

  setID <- copy(setID_ls)
  
  
  train <- which(data[[setID]]==1)
  
  tag_ddt <- cbind(data[,..tag_ls],data[,..sampleWeight])  
  num_ddt <- data[,..num_var_ls]
  
  num_ddt_train <- cbind(tag_ddt[train,],num_ddt[train,])
  num_ddt_validaton <- cbind(tag_ddt[-train,],num_ddt[-train,])

  
  abm_data <- list()
  abm_data$train_data <- num_ddt_train
  abm_data$validation_data <- num_ddt_validaton
  
  return(abm_data)
}


prepare_coarse_bin_abm_data <- function(abm_data,layout,nbins){
  train_data <- abm_data$train_data
  validation_data <- abm_data$validation_data
  
  
  ls1 <- layout$predictor_ls
  ls2 <- layout$tag_ls
  ls3 <- layout$sampleWeight_ls
  train_X_dt <- (abm_data$train_data)[,..ls1]
  train_Y_dt <- (abm_data$train_data)[,..ls2]
  train_W_dt <- (abm_data$train_data)[,..ls3]
  
  
  validation_X_dt <- (abm_data$validation_data)[,..ls1]
  validation_Y_dt <- (abm_data$validation_data)[,..ls2]
  validation_W_dt <- (abm_data$validation_data)[,..ls3]
  
  
  bin_split_pts <- get_split_pts(origin_bin_library(train_X_dt,nbins))
  
  
  # apply bining for train and test dataset 
  train_out <- do_cut_with_split_pts(train_X_dt,bin_split_pts)
  validation_out <- do_cut_with_split_pts(validation_X_dt,bin_split_pts)
  
  
  
  abm_data_bined <- list()
  abm_data_bined$train_X_dt_bined <- train_out$bine_dt
  abm_data_bined$train_Y_dt <- train_Y_dt
  abm_data_bined$train_W_dt <- train_W_dt
  
  abm_data_bined$validation_X_dt_bined <- validation_out$bine_dt
  abm_data_bined$validation_Y_dt <- validation_Y_dt
  abm_data_bined$validation_W_dt <- validation_W_dt
  
  abm_data_bined$group <- do_cut_with_split_pts(train_X_dt,bin_split_pts)$bined_group
  abm_data_bined$vnames <- colnames(abm_data_bined$train_X_dt_bined)
  return(abm_data_bined)
}



prepare_params <- function(lambda1,lambda2){
  params <- list()
  params$lambda1 <- lambda1
  params$lambda2 <- lambda2
  
  return(params)
}



train_fg_lasso_model_v1 <- function(abm_data_bined,params){
  
  lambda1 <- params$lambda1
  lambda2 <- params$lambda2
  
  X <- as.matrix(abm_data_bined$train_X_dt_bined)
  Y <- as.matrix(abm_data_bined$train_Y_dt)
  W <- as.matrix(abm_data_bined$train_W_dt)
  
  group <- abm_data_bined$group
  vnames <- abm_data_bined$vnames
  
  
  coef <- round(c(get_wt_group_lasso_and_group_fused_lasso_logistic_beta_with_cvxr(X,Y,W,group,lambda1,lambda2)),3)
  names(coef) <- vnames
  
  
  
  # model performance part 
  train_X <- as.matrix(abm_data_bined$train_X_dt_bined)
  train_Y <- c(as.matrix(abm_data_bined$train_Y_dt))
  train_W <- c(as.matrix(abm_data_bined$train_W_dt))

  idx1 <- c()
  for(i in c(1:length(train_W))){
    wgt <- round(train_W[i])
    idx1 <- c(idx1,rep(i,wgt))
  }

  train_ot <- test_logistic_model(coef,train_X[idx1,],train_Y[idx1])
  
  
  test_X <- as.matrix(abm_data_bined$validation_X_dt_bined)
  test_Y <- c(as.matrix(abm_data_bined$validation_Y_dt))
  test_W <- c(as.matrix(abm_data_bined$validation_W_dt))
  
  idx2 <- c()
  for(i in c(1:length(test_W))){
    wgt <- round(test_W[i])
    idx2 <- c(idx2,rep(i,wgt))
  }

  test_ot <- test_logistic_model(coef,test_X[idx2,],test_Y[idx2])
  
  
  
  cat('train\t',train_ot$prt_line)
  cat('test\t',test_ot$prt_line)
  
  
  
  train_report <- list()
  train_report$acc <- train_ot$acc
  train_report$auc <- train_ot$auc
  train_report$ks <- train_ot$ks
  train_report$loss <- train_ot$logit_loss
  
  test_report <- list()
  test_report$acc <- test_ot$acc
  test_report$auc <- test_ot$auc
  test_report$ks <- test_ot$ks
  test_report$loss <- test_ot$logit_loss
  
  model <- list()
  model$coef <- coef
  model$vnames <- vnames
  model$group <- group
  model$train_report <- train_report
  model$test_report <- test_report
  
  return(model)
}



train_fg_lasso_model_v2 <- function(abm_data_bined,params){
  
  lambda1 <- params$lambda1
  lambda2 <- params$lambda2
  
  X <- as.matrix(abm_data_bined$train_X_dt_bined)
  Y <- as.matrix(abm_data_bined$train_Y_dt)
  W <- as.matrix(abm_data_bined$train_W_dt)
  
  group <- abm_data_bined$group
  vnames <- abm_data_bined$vnames
  
  
  coef <- round(c(get_wt_group_lasso_and_group_fused_lasso_linear_regression_beta_with_cvxr(X,Y,W,group,lambda1,lambda2)),3)
  names(coef) <- vnames
  
  
  
  # model performance part 
  train_X <- as.matrix(abm_data_bined$train_X_dt_bined)
  train_Y <- c(as.matrix(abm_data_bined$train_Y_dt))
  train_W <- c(as.matrix(abm_data_bined$train_W_dt))
  
  idx1 <- c()
  for(i in c(1:length(train_W))){
    wgt <- round(train_W[i])
    idx1 <- c(idx1,rep(i,wgt))
  }
  
  train_ot <- test_logistic_model(coef,train_X[idx1,],train_Y[idx1])
  
  
  test_X <- as.matrix(abm_data_bined$validation_X_dt_bined)
  test_Y <- c(as.matrix(abm_data_bined$validation_Y_dt))
  test_W <- c(as.matrix(abm_data_bined$validation_W_dt))
  
  idx2 <- c()
  for(i in c(1:length(test_W))){
    wgt <- round(test_W[i])
    idx2 <- c(idx2,rep(i,wgt))
  }
  
  test_ot <- test_logistic_model(coef,test_X[idx2,],test_Y[idx2])
  
  
  
  cat('train\t',train_ot$prt_line)
  cat('test\t',test_ot$prt_line)
  
  
  
  train_report <- list()
  train_report$acc <- train_ot$acc
  train_report$auc <- train_ot$auc
  train_report$ks <- train_ot$ks
  train_report$loss <- train_ot$logit_loss
  
  test_report <- list()
  test_report$acc <- test_ot$acc
  test_report$auc <- test_ot$auc
  test_report$ks <- test_ot$ks
  test_report$loss <- test_ot$logit_loss
  
  model <- list()
  model$coef <- coef
  model$vnames <- vnames
  model$group <- group
  model$train_report <- train_report
  model$test_report <- test_report
  
  return(model)
}




show_bin <- function(model){
  coef <- model$coef
  grp <- levels(factor(model$group))
  drop_list <- get_drop_list(model$coef,model$group)
  
  idx <- which(model$group %in% drop_list == FALSE)
  
  ccoef <- coef[idx]
  
  print(ccoef)
  return(data.table(bin_nm<-names(ccoef),coef<-ccoef))
}

show_bin_plot <- function(model){
  coef <- model$coef
  plot(c(1:length(coef)),coef)
}

show_drop_list <- function(model){
  drop_list <- get_drop_list(model$coef,model$group)
  print(drop_list)
}

show_keep_list <- function(model){
  drop_list <- get_drop_list(model$coef,model$group)
  grp <- levels(factor(model$group))
  keep_list <- grp[grp %in% drop_list == FALSE]
  print(keep_list)
}





