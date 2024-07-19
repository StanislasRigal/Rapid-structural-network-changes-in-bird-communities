# Function to randomise communities
# Based on `make_netassoc_network` (Morueta-Holmes et al., 2016). The idea is that a species enters a community with a probability that is proportional to its abundance in the regional pool (the whole dataset).

randomizeM <- function(mat) { 
  # Compute total species frequencies
  total_freqs <- apply(mat, 2, sum)
  total_freqs <- total_freqs / sum(total_freqs)
  result <- apply(mat, 1, function(X) { 
    # We sample sum(X) individuals from the regional pool. With the probability 
    # of a species entering the community proportional to its abundance in the 
    # pool (the whole dataset)
    picked_species <- sample(seq.int(ncol(mat)), 
                             replace = TRUE, 
                             size = sum(X),
                             prob = total_freqs)
    tabulate(picked_species, nbins = ncol(mat))
  })
  t(result) 
}

# Function used to present results


tabularize <- function(coocmat, name = "assoc") { 
  # Zero-size matrix
  if(is.null(dim(coocmat))){return(NA)}
  # Add numbers if species names are absent in the original matrix
  if(is.null(colnames(coocmat))){ 
    spnames <- seq.int(ncol(coocmat))
  } else { 
    spnames <- colnames(coocmat)
  }
  # Create new table
  tab <- data.frame(expand.grid(spi = spnames, spj = spnames), 
                    as.vector(coocmat))
  names(tab) <- c('spi', 'spj', name)
  return(tab)
}

# Function to calculate species associations



main_fun_space <- function(mat, # observed data
            N.null, # number of null matrices to simulate
            occmin){ # minimum occurence number for a species to be selected
  
  mat[,c("code_point","year","habit","zonebio")] <- list(NULL)
  
  mat <- round(mat)
  
  # select species with more occurence than occmin
  matpa <- replace(mat, mat != 0, 1)
  mat <- mat[, colSums(matpa) > occmin]
  name <- names(mat)
  
  # log-transformation of abundance data as recommanded in Morueta-Holmes et al. 2016 (https://doi.org/10.1111/ecog.01892)
  mat <- log(mat+1e-6) - log(1e-6)
  
  # observed association matrix
  mat <- t(as.matrix(mat))
  D.obs <- partial_correlation(mat,"shrinkage")
  
  # random association matrices
  D.null <- array(rep(NA,nrow(mat)*nrow(mat)*N.null),
                  dim = c(nrow(mat),nrow(mat),N.null))
  
  for(null in 1:N.null) {
    D.null[ , ,null] <- partial_correlation(apply(replicate(1000, t(randomizeM(t(mat)))), c(1, 2), mean),
                                            "shrinkage")
  }
  
  D.null.mean <- apply(D.null, c(1,2), mean)
  D.null.sd <- apply(D.null, c(1,2), sd)
  
  SES <- (D.obs-D.null.mean)/D.null.sd
  
  # pvalue matrix
  M <- array(c(D.obs, D.null), dim=c(nrow(mat),nrow(mat),N.null+1))
  p.val <- apply(M, c(1,2), function(x) {
    if ( (rank(x)[1]/(N.null+1)) < 0.5 ) { 
      rank(x)[1]/(N.null+1)
    } else { 
      ( 1 - rank(x)[1]/(N.null+1) )
    } 
  }) 
  
  p.val <- 2*p.val
  diag(p.val) <- NA
  
  # adjusted pvalues
  p.val.adj <- matrix(p.adjust(p.val, method = "BH"), nrow=nrow(p.val), ncol=ncol(p.val))
  
  # select significant SES
  SES[p.val.adj>0.05] <- NA
  SES <- matrix(SES, nrow = nrow(p.val.adj), ncol = ncol(p.val.adj))
  
  # retrieve selected species names
  attr(SES,"dimnames") <- list(name)
  attr(SES,"dimnames")[[2]] <- attr(SES,"dimnames")[[1]]
  
  results <- data.frame(tabularize(D.obs, name = "obs"), 
                        spatial_asso = tabularize(SES)[ ,3], 
                        pval = tabularize(p.val.adj)[ ,3], 
                        nullmean = tabularize(D.null.mean)[ ,3], 
                        nullsd = tabularize(D.null.sd)[ ,3])
  
  return(results)
}
