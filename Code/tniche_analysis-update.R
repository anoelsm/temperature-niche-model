tniche_analysis <- function(file) {
  print("Loading....")
  nc <- nc_open(paste0(nc_path, file))
  ids <- rbind(stringr::str_extract(file, c("T\\d", "D\\d", "E\\d", "exp\\d")))
  colnames(ids) <- c("var", "disp", "mag", "m")
  print(ids)
  time_lastyr <- ((oneyr_iter*5)-oneyr_iter):(oneyr_iter*5)
  ### Convergence stuff for diversity purposes
  p_timeavg <- array(
    apply(ncvar_get(nc, "biomass")[2:14601,,], 3, function(x) {
      apply(x, 2, function(y) {
        as.vector(by(y, ceiling(five_years / oneyr_iter), function(z) {
          mean(z, na.rm = TRUE)
        }))})}), dim = c(5, num_species, bx))
  
  p_timediff <- array(apply(
    p_timeavg, 3, diff), dim = c(4, num_species, bx))
  
  div_list_here <- diversity_function(
    x = p_timediff,
    y = p_timeavg,
    z = ncvar_get(nc, "biomass")[2:14601,,])
  div_list_here$ids <- ids
  sp_here <- unique(unlist(div_list_here$species))
  which(unlist(lapply(div_list_here$species,length)) >3)

    converged_df <- do.call("rbind", lapply(1:bx, function(v){
      p1 <- data.frame(
        div_list_here$total_alpha[[v]][time_lastyr,] * ncvar_get(nc, "biomass")[time_lastyr, div_list_here$species[[v]], v])
      colnames(p1) <- paste0("Species.", div_list_here$species[[v]])
      p1$time_vec <- as.vector(time_vec[lastyear] / 365)
      p1$box <- v
      p1$temp <- ncvar_get(nc, "temp")[time_lastyr,v]
      p1 <- p1 %>%
        pivot_longer(cols = !c("box", "time_vec", "temp"), names_to = "species", values_to = "bio") %>%
        filter(bio > 0)
      return(p1)
    }))
    
  converged_df <- cbind(converged_df, ids)
  converged_df <- converged_df %>% 
    tidyr::drop_na(time_vec)
  
  
  converged_vals <- do.call("rbind", converged_df %>%
                              filter(bio > 0) %>%
                              group_by(species) %>%
                              group_map(~{
                                fit <- tryCatch(sm.regression(
                                  .x$temp,
                                  .x$bio,
                                  display = "none",
                                  eval.points = seq(
                                    (min(.x$temp)),
                                    (max(.x$temp)), by = 0.01), h = 1),
                                  error = function(e) "TryCatch")
                                if (fit[1] == "TryCatch") {
                                  df_outt <- data.frame(x = NA, y = NA)
                                  vals_out <- data.frame(
                                    tobs = .x$temp[1],
                                    wobs = 0,
                                    tdiff = diff(range(.x$temp)))
                                }else{
                                  reg_out <- data.frame(x = fit$eval.points, y = fit$estimate)
                                  tobs <- reg_out$x[which.max(reg_out$y)]
                                  wobs <- diff(quantile(reg_out$x, c(0.01, 0.99)))
                                  df_outt <- data.frame(reg_out)
                                  vals_out <- data.frame(
                                    tobs = tobs,
                                    wobs = wobs,
                                    tdiff = diff(range(.x$temp)))
                                }
                                return(vals_out)
                              }, .groups = "keep")) %>%
    mutate(name = group_keys(converged_df %>%
                               filter(bio > 0) %>%
                               group_by(species)) %>% pull(species))
  converged_vals$name <- as.numeric(gsub(
    "[^\\d]+", "",
    converged_vals$name,
    perl = TRUE))
  converged_vals <- cbind(converged_vals, ids)
  
  
  print("General Kenobi")
  
  if (all(c("T0","D0","E0","exp1") %in% c(ids))) {
    write.table(
      converged_vals,
      file = paste0(
      "./documents/metacomm_model/csv_output/converged_values_",
      subtitle,
      ".csv"),
      sep = ",", row.names = FALSE)

  }else{
    write.table(
      converged_vals,
      file = paste0(
        "./documents/metacomm_model/csv_output/converged_values_",
        subtitle,
        ".csv"),
      sep = ",",
      row.names = FALSE, append = TRUE, quote = FALSE, col.names = FALSE)
  }
  
  print("I have the high ground Anakin")
  
  nc_close(nc)
  
  return(list(div = div_list_here))
}