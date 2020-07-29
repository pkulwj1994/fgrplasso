##############  step0: load data
library(data.table)

data <- data.table(read_autoloan_dataset())

data[dealLoanToVal< -10000,dealLoanToVal:=NA]
data[cbMosAvg< -10000,cbMosAvg:=NA]
data[cbMosDlq< -10000,cbMosDlq:=NA]
data[cbMosInq< -10000,cbMosInq:=NA]
data[cbUtilizn< -10000,cbUtilizn:=NA]
data[cbPctGood< -10000,cbPctGood:=NA]
data[cbInq5Mos< -10000,cbInq5Mos:=NA]
data[cb90Ever< -10000,cb90Ever:=NA]
data[cbTimeFile< -10000,cbTimeFile:=NA]
data[appIncome< -10000,appIncome:=NA]
data[appTimeAddress< -10000,appTimeAddress:=NA]
data[appAge< -10000,appAge:=NA]




##############  step1: prepare layout 
tag_ls <- c('target')
sampleWeight_ls <- c('sampwt')
predictor_ls <- c('appAge','appTimeAddress','appIncome','dealLoanToVal','cbFICO','cbTimeFile','cbMosAvg','cbUtilizn','cb90Ever','cbPctGood','cbMosDlq','cbInq5Mos','cbMosInq')
#predictor_ls <-c('cbFICO')
#predictor_ls <-c('cbTimeFile')
#predictor_ls <-c('dealLoanToVal')
#predictor_ls <-c('cbMosAvg')
setID_ls <- c('trainFlg')


layout <- preparte_layout(tag_ls,sampleWeight_ls,predictor_ls,setID_ls)
##############  step2: prepare abm_data 

abm_data <- prepare_abm_data(data,layout$tag_ls,layout$sampleWeight_ls,layout$predictor_ls,layout$setID_ls)





############# step3: prepare coarse bined data
abm_data_bined <- prepare_coarse_bin_abm_data(abm_data,layout,100)



############ step4: run abm training 
lambda1 <- 0.0001
lambda2 <- 0.001

params <- prepare_params(lambda1,lambda2)






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


