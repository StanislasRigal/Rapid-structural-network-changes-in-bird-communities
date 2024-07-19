# Functions for multispatial CCM and S-map to obtain significatn temporal associations

# Detrend data

detrend_data <- function(series_all){
  series <- series_all$abond
  data <- data.frame(x=1:length(series), y=series)
  fit <- lm(y~x, data=data)
  fit_summary <- summary(fit)
  p_value <- fit_summary$coefficient[, 4]['x']
  
  if(!is.na(p_value)){
    if (p_value < 0.05){
      series_all$abond <- as.numeric(series - fit$fitted.values)
    }
  }

  return(series_all)
}

# Multispatial CCM

multisp_CCM_intern <- function(x, niter, i){
  
  x_multi_ccm <- droplevels(na.omit(x[,c(1,2,to_test[,i])]))
  
  if(nrow(x_multi_ccm)>500){
    while(nrow(x_multi_ccm)>500){
      lev <- levels(as.factor(x_multi_ccm$code_square))
      x_multi_ccm <- droplevels(x_multi_ccm[as.factor(x_multi_ccm$code_square) != lev[length(lev)],])
    }
  }
  
  if(length(levels(as.factor(x_multi_ccm$code_square)))>1 & nrow(x_multi_ccm)>=50){
    
    x_multi_ccm <- ddply(x_multi_ccm, .(code_square), .fun=add_row, .parallel = F)
    Accm <- x_multi_ccm[-nrow(x_multi_ccm),3]
    Bccm <- x_multi_ccm[-nrow(x_multi_ccm),4]
    
    maxE <- 5
    
    Emat <- matrix(nrow=maxE-1, ncol=2)
    
    colnames(Emat) <- c("A", "B")
    
    for(E in 2:maxE) {
      Emat[E-1,"A"] <- SSR_pred_boot(A=Accm, E=E, predstep=1, tau=1)$rho 
      Emat[E-1,"B"] <- SSR_pred_boot(A=Bccm, E=E, predstep=1, tau=1)$rho
    }
    
    E_A <- which.max(na.omit(Emat[,1]))+1
    E_B <- which.max(na.omit(Emat[,2]))+1
    
    if(length(E_A)>0 & length(E_B)>0){
      signal_A_out <- SSR_check_signal(A=Accm, E=E_A, tau=1, predsteplist=1:10)
      signal_B_out <- SSR_check_signal(A=Bccm, E=E_B, tau=1, predsteplist=1:10)
      
      if(summary(lm(signal_A_out$predatout$rho~signal_A_out$predatout$predstep))$coeff[2,1]<0 &
         summary(lm(signal_B_out$predatout$rho~signal_B_out$predatout$predstep))$coeff[2,1]<0){
        
        CCM_boot_A <- CCM_boot(Accm, Bccm, E_A, tau=1, iterations=niter)
        CCM_boot_B <- CCM_boot(Bccm, Accm, E_B, tau=1, iterations=niter)
        CCM_significance_test <- ccmtest(CCM_boot_A,CCM_boot_B)
      }
      if(summary(lm(signal_A_out$predatout$rho~signal_A_out$predatout$predstep))$coeff[2,1]<0 &
         summary(lm(signal_B_out$predatout$rho~signal_B_out$predatout$predstep))$coeff[2,1]>=0){
        
        CCM_boot_A <- CCM_boot(Accm, Bccm, E_A, tau=1, iterations=niter)
        CCM_significance_test <- ccmtest(CCM_boot_A,CCM_boot_A)
        CCM_significance_test[2] <- 1
      }
      if(summary(lm(signal_A_out$predatout$rho~signal_A_out$predatout$predstep))$coeff[2,1]>=0 &
         summary(lm(signal_B_out$predatout$rho~signal_B_out$predatout$predstep))$coeff[2,1]<0){
        
        CCM_boot_B <- CCM_boot(Bccm, Accm, E_B, tau=1, iterations=niter)
        CCM_significance_test <- ccmtest(CCM_boot_B,CCM_boot_B)
        CCM_significance_test[1] <- 1
      }
      if(summary(lm(signal_A_out$predatout$rho~signal_A_out$predatout$predstep))$coeff[2,1]>=0 &
         summary(lm(signal_B_out$predatout$rho~signal_B_out$predatout$predstep))$coeff[2,1]>=0){
        CCM_significance_test <- c(1,1)
      }
    }
  }else{
    CCM_significance_test <- c(1,1)
  }
  
  result <- data.frame(a_cause_b=CCM_significance_test[1], # a = species 1, b = species 2
                       b_cause_a=CCM_significance_test[2])
  names(result) <- c(paste0(names(x_multi_ccm)[3],"_cause_",names(x_multi_ccm)[4]),
                     paste0(names(x_multi_ccm)[4],"_cause_",names(x_multi_ccm)[3]))
  return(result)
}

multisp_CCM_intern2 <- function(x, niter, i){
  tryCatch(multisp_CCM_intern(x, niter, i),
           error=function(e) NA)
}


multisp_CCM <- function(x, niter=100){
  
  x <- droplevels(x[,c("code_sp","code_square","year","abond")])
  x <- dcast(x, code_square+year ~ code_sp, value.var="abond")
  
  if(ncol(x)>3){
    to_test <- combn(3:ncol(x),2)
    
    library(parallel)
    
    no_cores <- detectCores() - 1
    
    cl <- makeCluster(no_cores, type = "FORK")
    clusterExport(cl, c("multisp_CCM_intern","to_test","x","niter"), envir=environment())
    
    result_tot <- parLapply(cl, c(1:ncol(to_test)),
                            fun=function(y){
                              multisp_CCM_intern2(x = x,
                                                 niter = niter,
                                                 i = y)
                            })
    
    stopCluster(cl)
    
    
    res_multiCCM <- cbind(unlist(result_tot))
    
    res_multiCCM <- data.frame(species = row.names(res_multiCCM), rho = res_multiCCM[,1])
    
  }else{
    res_multiCCM <- data.frame(species = NA, rho = NA)
  }
  
  
  return(res_multiCCM)
}


E_dim <- function(x, i){
  
  x_multi_ccm <- droplevels(na.omit(x[,c(1,2,to_test[,i])]))
  
  
  if(length(levels(as.factor(x_multi_ccm$code_square)))>1 & nrow(x_multi_ccm)>=50){
    
    x_multi_ccm <- ddply(x_multi_ccm, .(code_square), .fun=add_row, .parallel = F)
    Accm <- x_multi_ccm[-nrow(x_multi_ccm),3]
    Bccm <- x_multi_ccm[-nrow(x_multi_ccm),4]
    
    maxE <- 5
    
    Emat <- matrix(nrow=maxE-1, ncol=2)
    
    colnames(Emat) <- c("A", "B")
    
    for(E in 2:maxE) {
      Emat[E-1,"A"] <- SSR_pred_boot(A=Accm, E=E, predstep=1, tau=1)$rho 
      Emat[E-1,"B"] <- SSR_pred_boot(A=Bccm, E=E, predstep=1, tau=1)$rho
    }
    
    E_A <- which.max(na.omit(Emat[,1]))+1
    E_B <- which.max(na.omit(Emat[,2]))+1
  }else{
    E_A <- E_B <- NA
  }
  
  result <- data.frame(E_a_cause_b=E_A, # a = species 1, b = species 2
                       E_b_cause_a=E_B)
  names(result) <- c(paste0("E_",names(x_multi_ccm)[3],"_cause_",names(x_multi_ccm)[4]),
                     paste0("E_",names(x_multi_ccm)[4],"_cause_",names(x_multi_ccm)[3]))
  return(result)
}

E_dim2 <- function(x, i){
  tryCatch(E_dim(x, i),
           error=function(e) NA)
}

Edim_multisp_CCM <- function(x){
  
  x <- droplevels(x[,c("code_sp","code_square","year","abond")])
  x <- dcast(x, code_square+year ~ code_sp, value.var="abond")
  
  if(ncol(x)>3){
    to_test <- combn(3:ncol(x),2)
    
    library(parallel)
    
    no_cores <- detectCores() - 1
    
    cl <- makeCluster(no_cores, type = "FORK")
    clusterExport(cl, c("E_dim2","to_test","x"), envir=environment())
    
    result_tot <- parLapply(cl, c(1:ncol(to_test)),
                            fun=function(y){
                              E_dim2(x = x,i = y)
                            })
    
    stopCluster(cl)
    
    
    res_multiCCM <- cbind(unlist(result_tot))
    
    res_multiCCM <- data.frame(species = row.names(res_multiCCM), rho = res_multiCCM[,1])
    
  }else{
    res_multiCCM <- data.frame(species = NA, rho = NA)
  }
  
  
  return(res_multiCCM)
}


# S-map

smap_fun <- function(block){
  
  block <- block[,c("spA","spB")]
  
  if(sum(abs(diff(block$spA))) != 0 & sum(abs(diff(block$spB))) != 0){
    block_sub <- data.frame(apply(block,2,scale))
    
    # Determine the best theta
    theta.examined <- seq(0, 10, by = 0.1)
    th.test <-
      pforeach(
        i      = theta.examined,
        .c     = rbind,
        .cores = 7
      )({
        th.test0 <- block_lnlp(block_sub, method = "s-map", tp = 1,columns=c(1:dim(block_sub)[2]),
                               target_column = 1,theta  = i,
                               silent = T, num_neighbors = 0)
      })
    
    best.th <- th.test[th.test$mae == min(th.test$mae), 'theta']
    
    # Perform multivariate S-map to quantify interaction strengths
    smapc.res <- block_lnlp(block_sub, method = "s-map", tp = 1, columns=c(1:dim(block_sub)[2]),
                            target_column = 1, theta  = best.th,
                            num_neighbors = 0, silent = T,
                            save_smap_coefficients = T)
    smapc.tmp <- data.frame(smapc.res[[1]]$smap_coefficients)
    
    colnames(smapc.tmp) <- c(colnames(block_sub), "Constant")
    
    rho <- smapc.res[[1]]$stats$rho
    
    n_pred <- smapc.res[[1]]$stats$num_pred
    if(!is.na(rho)){
      if(abs(rho)<1){
        t <- rho*sqrt(n_pred-2)/sqrt(1-rho^2)
        if (t >= 0){
          pvalue <- 1 - pt(t, df = n_pred-2)    
        } else {
          pvalue <- pt(t, df = n_pred-2)
        }
      }else{pvalue <- 0}
    }else{pvalue <- 1}
    
    theta <- round(best.th, 2)
  }else{
    smapc.tmp <- data.frame(spA=NA,spB=NA,Constant=NA)
    rho <- pvalue <- theta <- NA
  }
  
  return(list(coefficients = smapc.tmp, rho = rho, pvalue = pvalue, theta = theta))
}

smap_fun_signif <- function(res_ccm,data_ccm){
  res_ccm <- droplevels(res_ccm)
  data_ccm <- droplevels(data_ccm[which(data_ccm$code_sp %in% c(res_ccm$spA,res_ccm$spB) &
                                          data_ccm$zonebio==res_ccm$zonebio[1] &
                                          data_ccm$habit==res_ccm$habit[1]),])
  
  data_ccm$code_sp <- as.character(data_ccm$code_sp)
  data_ccm$code_sp[data_ccm$code_sp==res_ccm$spA] <- "spA"
  data_ccm$code_sp[data_ccm$code_sp==res_ccm$spB] <- "spB"
  
  data_ccm2 <- droplevels(na.omit(dcast(data_ccm, code_square+year ~ code_sp, value.var="abond")))
  
  res_smap <- dlply(data_ccm2, .(code_square), .fun=smap_fun)
  
  smap_sp_res <- ldply(res_smap, .fun=function(x){
    
    if(!is.na(x$pvalue)){
      if(x$pvalue<0.05){
        spB_smap <- na.omit(x$coefficients$spB)
        spB_smap2 <- na.omit(x$coefficients$spB)
      }else{
        spB_smap <- NA
        spB_smap2 <- na.omit(x$coefficients$spB)
      }
    }else{
      spB_smap <- spB_smap2 <- NA
    }
    
  
  return(data.frame(spB_smap, spB_smap2))
  })
    
  return(smap_sp_res)
}


smap_fun_rand <- function(block){
  
  block <- block[,c("spA","spB")]
  
  if(sum(abs(diff(block$spA))) != 0 & sum(abs(diff(block$spB))) != 0){
    block_sub <- data.frame(apply(block,2,scale))
    
    # Determine the best theta
    theta.examined <- seq(0, 10, by = 0.1)
    th.test <-
      pforeach(
        i      = theta.examined,
        .c     = rbind,
        .cores = 7
      )({
        th.test0 <- block_lnlp(block_sub, method = "s-map")
      })
    th.test$mae <- unlist(th.test$mae)
    best.th <- th.test[th.test$mae == min(th.test$mae), 'theta']
    
    # Perform multivariate S-map to quantify interaction strengths
    smapc.res <- block_lnlp(block_sub, method = "s-map",  theta  = best.th[1],
                            num_neighbors = 0, silent = T,
                            save_smap_coefficients = T)
    smapc.tmp <- data.frame(smapc.res$smap_coefficients)
    
    colnames(smapc.tmp) <- c(colnames(block_sub), "Constant")
    
    rho <- smapc.res$rho[[1]]
    
    n_pred <- smapc.res$num_pred
    if(!is.na(rho)){
      if(abs(rho)<1){
        t <- rho*sqrt(n_pred-2)/sqrt(1-rho^2)
        if (t >= 0){
          pvalue <- 1 - pt(t, df = n_pred-2)    
        } else {
          pvalue <- pt(t, df = n_pred-2)
        }
      }else{pvalue <- 0}
    }else{pvalue <- 1}
    
    theta <- round(best.th[1], 2)
  }else{
    smapc.tmp <- data.frame(spA=NA,spB=NA,Constant=NA)
    rho <- pvalue <- theta <- NA
  }
  
  return(list(coefficients = smapc.tmp, rho = rho, pvalue = pvalue, theta = theta))
}

smap_fun_signif_rand <- function(res_ccm,data_ccm){
  res_ccm <- droplevels(res_ccm)
  data_ccm <- droplevels(data_ccm[which(data_ccm$code_sp %in% c(res_ccm$spA,res_ccm$spB)),])
  
  data_ccm$code_sp <- as.character(data_ccm$code_sp)
  data_ccm$code_sp[data_ccm$code_sp==res_ccm$spA] <- "spA"
  data_ccm$code_sp[data_ccm$code_sp==res_ccm$spB] <- "spB"
  
  data_ccm2 <- droplevels(na.omit(dcast(data_ccm, code_square+year ~ code_sp, value.var="abond")))
  
  res_smap <- dlply(data_ccm2, .(code_square), .fun=smap_fun_rand)
  
  smap_sp_res <- ldply(res_smap, .fun=function(x){
    
    if(!is.na(x$pvalue)){
      if(x$pvalue<0.05){
        spB_smap <- na.omit(x$coefficients$spB)
        spB_smap2 <- na.omit(x$coefficients$spB)
      }else{
        spB_smap <- NA
        spB_smap2 <- na.omit(x$coefficients$spB)
      }
    }else{
      spB_smap <- spB_smap2 <- NA
    }
    
    
    return(data.frame(spB_smap, spB_smap2))
  })
  
  return(smap_sp_res)
}

