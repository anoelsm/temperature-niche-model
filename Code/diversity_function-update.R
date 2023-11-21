diversity_function <- function(x, y, z) {
  extinction <- 1E-4
  
  box_mat <- lapply(1:bx, function(r){
    
    box_sp <- which(unlist(lapply(1:num_species, function(d){
      cond_1 <- all(x[,d,r] >= 0 | diff(round(y[,d,r], 3)) >= 0)
      cond_2 <- all(y[,d,r] >= (max(y[,,r], na.rm = T) * extinction))
      return(all(cond_1, cond_2))
        })))
    
    val <- mean(apply(apply(
      as.matrix(z[, box_sp, r]), 2, function(l) {
        cond <- which(l >= (max(as.matrix(z[, box_sp, r]),
                                na.rm = T) * extinction))
        rows_all <- 1:length(l)
        same_vals <- rows_all %in% cond}), 1, sum))
    
    val_mat <- apply(as.matrix(z[, box_sp, r]), 2,
                     function(l) {
                       cond <- which(l >= (max(z[, box_sp, r], na.rm = T) * extinction))
                       rows_all <- 1:length(l)
                       same_vals <- rows_all %in% cond})
    
    z[,-box_sp,r] <- 0
    time_shannon <- diversity(z[,,r], MARGIN = 1, "shannon")
    avg_shannon <- mean(time_shannon)
    max_shannon <- max(time_shannon)
    
    #val_mat[val_mat<1] = 1

  return(list(
    alpha = val,
    total_alpha = val_mat,
    species = box_sp,
    gamma = length(box_sp),
    avg_shannon = avg_shannon,
    max_shannon = max_shannon))
  })
  
  return(list(
    species = lapply(box_mat,function(c){c$species}),
    alpha = unlist(sapply(box_mat,function(c){c$alpha})),
    gamma = unlist(sapply(box_mat,function(c){c$gamma})),
    total_alpha = lapply(box_mat,function(c){c$total_alpha}),
    tot_avg_shannon = unlist(sapply(box_mat,function(c){c$avg_shannon})),
    tot_max_shannon = unlist(sapply(box_mat,function(c){c$max_shannon}))))
}