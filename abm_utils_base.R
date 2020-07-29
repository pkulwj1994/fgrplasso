#####################

############ our algorithms for fused group lasso binning 
library(data.table)
library(CVXR)
library(ROCit)




read_autoloan_dataset <- function(){
  path <- './data/AutoLoan.csv'
  ddt <- fread(path)
  return(ddt)
}





discrete_encoding <- function(char_dt){
  char_ddt <- copy(char_dt)
  
  char_var_ls <- copy(colnames(char_ddt))
  
  char_var_values <- list()
  for(i in c(1:length(char_var_ls))){
    nm <- char_var_ls[i]
    tmp <- unique(char_ddt[[nm]])
    cmd_line <- paste0('char_var_values$',nm,' <- tmp')
    eval(parse(text = cmd_line))
  }
  
  
  
  for(i in c(1:length(char_var_ls))){
    
    vnm <- char_var_ls[i] 
    vvls <- char_var_values[[i]]
    
    for(j in c(1:(length(vvls)+1))){
      if(j != (length(vvls)+1)){
        prefix <- paste0(vnm,'_')
        nm <- paste0(prefix,vvls[j])
        vvl <- vvls[j]
        
        cmd_line <- paste0('char_ddt[,',nm,':=ifelse((is.na(',vnm,')==FALSE) & ',vnm,'==','\'',vvl,'\'',',1,0)]')
        eval(parse(text=cmd_line))
      }
      
      if(j == (length(vvls)+1)){
        prefix <- paste0(vnm,'_')
        nm <- paste0(prefix,'mv')
        
        cmd_line <- paste0('char_ddt[,',nm,':=ifelse(is.na(',vnm,'),1,0)]')
        eval(parse(text=cmd_line))
      }
      
      
    }
  }
  
  char_dddt <- char_ddt[,-..char_var_ls]
  
  
  return(char_dddt)
}










origin_bin_library <- function(ddt,nbins){
  dt <- copy(ddt)
  
  bin_libraries <- list()
  
  vnames <- colnames(dt)
  for(i in c(1:length(vnames))){
    vname <- vnames[i]
    a <- dt[[vname]]
    
    bin_library <- c(levels(cut(a,nbins,right=FALSE)))
    
    cmd_line <- paste0('bin_libraries$',vname,' <- bin_library')
    eval(parse(text = cmd_line))
  }
  
  return(bin_libraries)
}





get_split_pts <- function(bin_libraries){
  
  lis <- bin_libraries
  vnames <- names(lis)
  
  
  bin_split_pts <- list()
  for(i in c(1:length(vnames))){
    vname <- vnames[i]
    
    cmd_line <- paste0('bl <- lis$',vname)
    eval(parse(text = cmd_line))
    
    
    split_pts <- c(-Inf)
    for(i in c(1:length(bl))){
      bin <- bl[i]
      bds <- as.numeric(strsplit(gsub('\\[','',gsub(')','',bin)),split=',')[[1]])
      lft_b <- bds[1]
      split_pts <- c(split_pts,lft_b)
      
    }
    split_pts <- c(split_pts,Inf)
    
    cmd_line <- paste0('bin_split_pts$',vname,' <- split_pts')
    eval(parse(text = cmd_line))
  }
  return(bin_split_pts)
}





do_cut_with_split_pts <- function(ddt,bin_split_pts){
  dt <- copy(ddt)
  
  # initialize group and fuse pos
  group <- c()
  fuse_pos <- c()
  
  vnames <- names(bin_split_pts)
  for(i in c(1:length(vnames))){
    vname <- vnames[i]
    
    eval(parse(text = paste0('bsp <- bin_split_pts$',vname)))
    
    cmd_line <- paste0('dt[,',vname,':=cut(',vname,',bsp,right = FALSE)]')
    eval(parse(text = cmd_line))
  }
  
  for(i in c(1:length(vnames))){
    vname <- vnames[i]
    
    eval(parse(text = paste0('levs <- levels(dt$',vname,')')))
    
    prefix <- paste0(vname,'_')
    
    eval(parse(text = paste0('dt[,paste0(prefix,\'mv\'):=ifelse(is.na(',vname,'),1,0)]')))
    group <- c(group,vname)
    for(j in c(1:length(levs))){
      subvname <- paste0(prefix,levs[j])
      subvname <- gsub('\\[','',subvname)
      subvname <- gsub('-','..',subvname)
      subvname <- gsub(')','',subvname)
      subvname <- gsub(',','_',subvname)
      subvname <- gsub(' ','',subvname)
      subvname <- gsub('\\+','',subvname)
      
      eval(parse(text = paste0('dt[,',subvname,':=ifelse(is.na(',vname,')==F & (',vname,'==levs[j]),1,0)]')))
      group <- c(group,vname)
    }
  }
  dt[,Intercept:=1]
  group <- c(group,c('Intercept'))
  
  
  var_ls <- colnames(dt)
  var_ls <- var_ls[var_ls %in% vnames == F]
  
  
  out <- list()
  out$bine_dt <- dt[,..var_ls]
  out$bined_group <- group
  return(out)
}





##########################

############### on optimize solver routines
get_wt_group_lasso_and_group_fused_lasso_logistic_beta_with_cvxr <- function(X,Y,W,group,lambda1=0.0,lambda2=0.0){
  # logistic regression 
  suppressWarnings(library(CVXR, warn.conflicts=FALSE))
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  # c('age_mv','age_1','age_2'), group == c('age','age','age')
  #group <- factor(c('age','age','age','height','height','width','width','width','width','width'))
  lev <- levels(factor(group))
  
  
  beta <- Variable(p)
  gplasso <- function(beta,idx) {
    lasso <- p_norm(beta[idx],2)
    return(lasso)
  }
  
  fusedlasso <- function(beta,idx){
    if(length(idx)<=1){lasso <- 0}
    if(length(idx)==2){lasso <- 0}
    if(length(idx)>2){
      f_idx <- idx[c(2:length(idx))]
      lasso <- p_norm(diff(beta[f_idx]),1)}
    return(lasso)
  }
  
  constrs <- list()
  
  #k = 1
  obj <- -sum(c(W[Y == 0])*X[Y == 0,] %*% beta)/sum(W) - sum(c(W)*logistic(-X %*% beta))/sum(W)
  for(i in 1:length(lev)){
    idx <- which(group==lev[i]) 
    obj <- obj - lambda1*gplasso(beta,idx) - lambda2*fusedlasso(beta,idx)
    
    #if(length(idx)>1){
    #constrs[k] <- diff(beta[idx])>=0
    #k <- k+1
    #}
  }
  
  prob <- Problem(Maximize(obj))
  result <- solve(prob)
  return(result$getValue(beta))
}





get_wt_group_lasso_and_group_fused_lasso_linear_regression_beta_with_cvxr <- function(X,Y,W,group,lambda1=0.0,lambda2=0.0){
  # logistic regression
  suppressWarnings(library(CVXR, warn.conflicts=FALSE))
  n <- dim(X)[1]
  p <- dim(X)[2]

  # c('age_mv','age_1','age_2'), group == c('age','age','age')
  #group <- factor(c('age','age','age','height','height','width','width','width','width','width'))
  lev <- levels(factor(group))


  beta <- Variable(p)
  gplasso <- function(beta,idx) {
    lasso <- p_norm(beta[idx],2)
    return(lasso)
  }

  fusedlasso <- function(beta,idx){
    if(length(idx)<=1){lasso <- 0}
    if(length(idx)==2){lasso <- 0}
    if(length(idx)>2){
      f_idx <- idx[c(2:length(idx))]
      lasso <- p_norm(diff(beta[f_idx]),1)}
    return(lasso)
  }

  constrs <- list()

  #k = 1
  obj <- -sum(W*square(Y-X%*%beta))/sum(W)
  #obj <- -sum(c(W[Y == 0])*X[Y == 0,] %*% beta)/sum(W) - sum(c(W)*logistic(-X %*% beta))/sum(W)
  for(i in 1:length(lev)){
    idx <- which(group==lev[i])
    obj <- obj - lambda1*gplasso(beta,idx) - lambda2*fusedlasso(beta,idx)

    #if(length(idx)>1){
    #constrs[k] <- diff(beta[idx])>=0
    #k <- k+1
    #}
  }

  prob <- Problem(Maximize(obj))
  result <- solve(prob)
  return(result$getValue(beta))
}








get_group_lasso_and_group_fused_lasso_logistic_beta_with_cvxr <- function(X,Y,W,group,lambda1=0.0,lambda2=0.0){
  # logistic regression 
  suppressWarnings(library(CVXR, warn.conflicts=FALSE))
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  
  #group <- factor(c('age','age','age','height','height','width','width','width','width','width'))
  lev <- levels(factor(group))
  
  
  beta <- Variable(p)
  gplasso <- function(beta,idx) {
    lasso <- p_norm(beta[idx],2)
    return(lasso)
  }
  fusedlasso <- function(beta,idx){
    if(length(idx)<=1){lasso <- 0}
    if(length(idx)==2){lasso <- 0}
    if(length(idx)>2){
      f_idx <- idx[c(2:length(idx))]
      lasso <- p_norm(diff(beta[f_idx]),1)}
    return(lasso)
  }
  
  constrs <- list()
  
  #k = 1
  obj <- -sum(X[Y == 0,] %*% beta)/n - sum(logistic(-X %*% beta))/n
  for(i in 1:length(lev)){
    idx <- which(group==lev[i]) 
    obj <- obj - lambda1*gplasso(beta,idx) - lambda2*fusedlasso(beta,idx)
    
    #if(length(idx)>1){
    #constrs[k] <- diff(beta[idx])>=0
    #k <- k+1
    #}
  }
  
  prob <- Problem(Maximize(obj))
  result <- solve(prob)
  return(result$getValue(beta))
}




##########################

############### on model evaluation routines
sigmoid_v1 <- function(logit){
  return(exp(logit)/(exp(logit)+1))
}

softmax_v1 <- function(logit){
  return(exp(logit)/sum(exp(logit)))
}



test_logistic_model <- function(coef,test_X,test_Y){
  library(ROCit)
  
  test_logit <- c(test_X%*%coef)
  test_prob <- sigmoid_v1(test_logit)
  test_pred <- ifelse(test_prob>0.5,1,0)
  
  ROC <- rocit(c(test_prob),test_Y)
  plot(ROC)
  
  KS_Plot <- ksplot(ROC)
  
  
  
  acc <- sum(test_Y==test_pred)/sum(is.na(test_Y)==F)
  #auc <- compute_auc(test_prob,test_Y)
  aauc <- 100*ROC$AUC
  #ks <- compute_ks(test_Y,test_prob)
  ks <- 100*KS_Plot$`KS stat`
  #plot(kscurve)
  
  
  logit_loss <- (sum(test_X[test_Y == 0,] %*% as.matrix(coef))/length(c(test_Y)) + sum(sigmoid_v1(test_X %*% as.matrix(coef))))/length(c(test_Y))
  
  
  
  prt_line <- paste0('acc is ',round(100*acc,2),' ','auc is ',round(aauc,2),' ','ks is ',round(ks,2),' ','logit loss is ',round(logit_loss,4),'\n')
  
  
  out <- list()
  out$acc <- acc
  out$auc <- aauc
  out$ks <- ks
  out$logit_loss <- logit_loss
  out$prt_line <- prt_line
  
  return(out)
}


#################

########### model interpretation 

print_coef <- function(coef){
  nms <- names(coef)
  for(i in c(1:length(coef))){
    prt_line <- paste0(nms[i],'\t',coef[i],'\n')
    cat(prt_line)
  }
}


plot_coef <- function(coef){
  plot(c(1:length(coef)),coef)
}


get_drop_list <- function(coef,group){
  
  drop_list <- c()
  levs <- levels(factor(group))
  for(i in c(1:length(levs))){
    lev <- levs[i]
    idx <- which(group==lev)
    if(sum(coef[idx])==0){
      drop_list <- c(drop_list,lev)
    }
  }
  
  return(drop_list)
}



###############################
#### Purpose: Calculate KS ####
#### By: William
#### Date: 
###############################

creat_binaryTarget <- function(dt, target, bad, good, varName=NULL){
  if(is.null(varName)){varNm <- paste0(target, '_tmp')}
  else{varNm <- varName}
  set(dt, j=varNm, value=NA_integer_)
  set(dt, i=which(dt[[target]] %in% bad),j=varNm, value=0)
  set(dt, i=which(dt[[target]] %in% good),j=varNm, value=1)
}


performance_calc <- function(dt, scores, targets, byVars=NULL, weight=NULL, plot=F){
  dt_tmp <- dt[,.SD, .SDcols=c(scores, targets, byVars, weight)]
  if(is.null(weight)){dt_tmp[['weight']] <- 1}
  else{setnames(dt_tmp,weight, 'weight')}
  KS_summary_all <- list()
  KS_chart_all <- list()
  for(y in targets){
    for(x in scores){
      eval(parse(text=paste0('dt_tmp[,\'',x,'\':=round(',x, ',digits=0)]')))
      set(dt_tmp, i=which(dt_tmp[[x]] %in% c(-99000800,0)),j=x,value=NA)
      eval(parse(text=paste0('tmp <- dt_tmp[!(is.na(',x,')|is.na(',y,')),.(B_cnt=sum(',y,'==0),G_cnt=sum(',y,'==1),B_wt=sum(weight[',y,'==0]),G_wt=sum(weight[',y,'==1]),Tot_wt=sum(weight[',y,' %in% c(0,1)])), keyby=c(byVars,\'',x,'\')][,c(\'B_cumpct\',\'G_cumpct\', \'G_cumpct\', \'Tot_cumpct\'):=.(cumsum(B_wt)/sum(B_wt), cumsum(G_wt)/sum(G_wt), cumsum(Tot_wt)/sum(Tot_wt)), by=byVars]' )))
      eval(parse(text=paste0('KS_summary <- tmp[,.(target=\'',y,'\',score=\'',x,'\',KS=(B_cumpct-G_cumpct)[which.max(abs(B_cumpct-G_cumpct))],ROC=sum((B_cumpct+shift(B_cumpct, fill=0))*0.5*diff(c(0,G_cumpct))),KS_score=',x,'[which.max(abs(B_cumpct-G_cumpct))],KS_percentile=Tot_cumpct[which.max(abs(B_cumpct-G_cumpct))],bad_rawCnt=sum(B_Cnt),good_rawCnt=sum(G_cnt),bad_ftrCnt=sum(B_wt),good_ftrCnt=sum(G_wt)), by=byVars')))
      KS_chart <- tmp[,c(list(target_nm=y,score_nm=x),.SD), .SDcols=c(byVars, x, 'B_cumpct', 'G_cumpct', 'Tot_cumpct')]
      setnames(KS_chart,x,'value')
      if(length(KS_summary_all)==0){
        KS_summary_all <- KS_summary
      }else{
        KS_summary_all <- rbindlist(list(KS_summary_all,ks_summary), use.names = T)
      }
      if(length(KS_chart_all)==0){
        KS_chart_all <- KS_chart;
      }
      else{
        KS_chart_all <- rbindlist(list(KS_chart_all,KS_chart), use.names = T);
      }
    }
  }
  out <- list()
  out$KS_summary <- KS_summary_all
  out$KS_chart <- KS_chart_all
  
  
  if(plot==T){
    eval(parse(text=paste0('KS_chart_all[, \'Segment\':=paste(',paste0(byVars,collapse = ','),',score_nm,sep=\'+\')]')))
    for(y in targets){
      tmp <- ggplot(data=KS_chart_all[target_nm==y], aes(x=G_cumpct, y=B_cumpct, colour=Segment))+geom_line(size=1.2)+xlim(c(0,1))+ylim(c(0,1))+geom_abline()+ggtitle(y)+xlab('FPR')+ylab('TPR')
      print(tmp)
      tmp <- ggplot(data=KS_chart_all[target_nm==y], aes(x=value, y=Tot_cumpct, colour=Segment))+geom_line(size=1.2)+xlim(c(300,1000))+ylim(c(0,1)) + geom_abline() + ggtitle(y) + xlab('score')+ylab('Cum %')
      print(tmp)
    }
  }
  return(out)
}














