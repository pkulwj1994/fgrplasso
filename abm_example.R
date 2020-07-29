
# abm_example


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


model1 <- train_fg_lasso_model_v1(abm_data_bined,params)
model2 <- train_fg_lasso_model_v2(abm_data_bined,params)


############ step5: model visualization

bin_lib1 <- show_bin(model1)
View(bin_lib1)

bin_lib2 <- show_bin(model2)
View(bin_lib2)

show_bin_plot(model)

show_drop_list(model)

show_keep_list(model)




# suggestions:
# 1*.smart rounding;
# 2*.discrete vars;
# 3.bin plot drops all remaining vars;
# 4.visualize as histogram ??;
# 5.visualize woes on plots; (before and after with whole vars);
# 6.visualize binning;
# 7.use -999-like number to substitute NA;
# 8.





# german credit example
##############  step0: load data
library(data.table)

path1 <- './data/german.data'
path2 <- './data/german.data-numeric'

data1 <- data.table(read.table(path1))

colnames(data1) <- c('Account_status',
                     'Duration',
                     'Purpose',
                     'Credit_history',
                     'Credit_amount',
                     'Savings',
                     'Employment_since',
                     'Installment_rate',
                     'Personal_status_and_sex',
                     'Other_debtors',
                     'Residence_since',
                     'Property',
                     'Age',
                     'Installment_plans',
                     'Housing',
                     'Existing_credits',
                     'Job',
                     'Number_of_liable_people',
                     'Telephone',
                     'Foreign_worker',
                     'target')
#target 2 bad 1 good

data1[,target:= ifelse(target %in% c(1),0,1)]

data1[,weight:=1.0]
data1[target %in% c(1),weight:=5.0]


numerical_vars <- c('Duration','Credit_amount','Installment_rate','Residence_since',
              'Age','Existing_credits','Number_of_liable_people')

discrete_vars <- c('Account_status',
                   'Purpose',
                   'Credit_history',
                   'Savings',
                   'Employment_since',
                   'Personal_status_and_sex',
                   'Other_debtors',
                   'Property',
                   'Installment_plans',
                   'Housing',
                   'Job',
                   'Telephone',
                   'Foreign_worker')


set.seed(2020)
train <- sample(nrow(data1),0.6*nrow(data1))

data1[,setID:=2]
data1[train,setID:=1]


fwrite(data1,'./data/german_processed.csv')


##############  step1: prepare layout 
tag_ls <- c('target')
sampleWeight_ls <- c('weight')
predictor_ls <- c('Duration','Credit_amount','Installment_rate','Residence_since','Age','Existing_credits','Number_of_liable_people')
setID_ls <- c('setID')


layout <- preparte_layout(tag_ls,sampleWeight_ls,predictor_ls,setID_ls)
##############  step2: prepare abm_data 
path <- './data/german_processed.csv'
data <- fread(path)[,c('target','weight','setID',
                       'Duration','Credit_amount','Installment_rate','Residence_since','Age','Existing_credits','Number_of_liable_people')]


abm_data <- prepare_abm_data(data,layout$tag_ls,layout$sampleWeight_ls,layout$predictor_ls,layout$setID_ls)





############# step3: prepare coarse bined data
abm_data_bined <- prepare_coarse_bin_abm_data(abm_data,layout,10)



############ step4: run abm training 
lambda1 <- 0.01
lambda2 <- 0.01

params <- prepare_params(lambda1,lambda2)


model <- train_fg_lasso_model(abm_data_bined,params)


############ step5: model visualization

bin_lib <- show_bin(model)

View(bin_lib)

show_bin_plot(model)

show_drop_list(model)











