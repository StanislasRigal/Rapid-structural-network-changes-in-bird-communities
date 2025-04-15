# Function to calculate metrics of species association networks

community_nb_link_init <- function(x, data_interaction_possible){
  
  # calculate indices for the observed data
  
  b <- levels(droplevels(x$code_sp))
  
  datasso <- data_interaction_possible[which(data_interaction_possible$spA %in% b & data_interaction_possible$spB %in% b),c("spA","spB","obs_temp_asso")]
  
  adj_mat <- datasso
  for(i in 1:length(b)){
    adj_mat <- rbind(adj_mat, data.frame(spA=b[i],spB=b[i],obs_temp_asso=0))
  }
  adj_mat <- reshape2::dcast(adj_mat, spA ~ spB, value.var="obs_temp_asso")
  row.names(adj_mat) <- adj_mat$spA
  adj_mat$spA <- NULL
  
  # degree ditribution
  
  adj_mat[lower.tri(adj_mat)] <- t(adj_mat)[lower.tri(adj_mat)]
  data_names <- data.frame(t(setNames(rep(0,length(unique(dataprp$code_sp))), sort(unique(dataprp$code_sp)))))
  deg_dist <- data.frame(t(apply(adj_mat,1,sum)))
  deg_dist <- rbind.fill(data_names,deg_dist)[2,]
  
  # modularity
  
  adj_mat[lower.tri(adj_mat)] <- 0
  modularity <- netcarto(as.matrix(adj_mat))
  
  # connectance
  
  connectance <- sum(datasso$obs_temp_asso)/nrow(datasso)
  
  result <- cbind(data.frame(nb_association = sum(datasso$obs_temp_asso), # total number of observed associations
                             nb_sp = length(b), # total number of species
                             nb_asso_pot = nrow(datasso), # total number of potential association
                             connectance = connectance,
                             modularity = modularity[[2]]),
                  deg_dist)
  return(result)
}

community_nb_link <- function(x,data_interaction_possible){tryCatch(community_nb_link_init(x,data_interaction_possible),
                                          error=function(e) cbind(data.frame(nb_association = NA,
                                                                             nb_sp = NA,
                                                                             nb_asso_pot = NA,
                                                                             connectance = NA,
                                                                             modularity = NA)),
                                          data.frame(t(setNames(rep(NA,length(unique(dataprp$code_sp))), sort(unique(dataprp$code_sp))))))}
