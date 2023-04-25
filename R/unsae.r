#' @importFrom plgp randomForest lme4
# '
#' inverse canonical link
#'
#' a simple function for inverse canonical link
#'
#' calculate the inverse canonical link \code{\link{unsae}}
#'
#' @param x a numeric vector
#'
#' @return calculated value
#'
#' @examples
#' hc(1)
#'
#' @export
hc <- function(x) as.numeric(1 /(1 + exp(-x))) # inverse canonical link


#' spatial likelihood
#'
#' a simple (log) spatial likelihood
#' used within multilevel_EM function of \code{\link{unsae}}
#'
#' @param par parameter values
#'
#' @param D distance matrix from plgp::distance
#'
#' @param Y "observed" response vector
#'
#' @param mu current estimate of mu
#'
#' @return negative log likelihood
#'
spatial_lik <- function(par, D, Y, mu) {
  ## from Gramacy (2020)
  # par: parameters (range and nugget)
  # D: coordinates
  # Y: responses
  # mu: current estimate of mu
  theta <- par[1]
  g <- par[2]
  n <- length(Y)
  K <- exp(-D/theta) + g*diag(n)
  K_inv <- solve(K)
  ldetK <- determinant(K, logarithm=TRUE)$modulus
  loglik <- -(n/2)*log(t(Y - mu) %*% K_inv %*% Y - mu) - (1/2)*ldetK
  return(-loglik)
}

#' multilevel EM
#'
#' a main function to fit multilevel em with gaussian spatial correlation. Some chunk of the code has been dedicated to handle the syntax. glmer function is called to get the initial value of the estimates.
#'
#' @param formula formula argument following the glmer syntax.
#'
#' @param tol tolerence level to control the em convergence. default value 1e-3.
#'
#' @param data name of the data set.
#'
#' @param coordinates specifying the name of coordinates within the data. These columns *should* exist within the data set provided.
#'
#' @return list of element that are needed for the prediction.
#'
#' @export


multilevel_EM <-
  function(formula, data, tol = 1e-3, coordinates){
    re_pattern <- "\\|(.*?)\\)"
    random_effect <- str_squish(str_match(formula, re_pattern)[,
                                                               2])
    if (!is.factor(data[[random_effect]])) {
      stop("area index must be a factor")
    }
    data$area_index <- droplevels(data[[random_effect]]) %>%
      as.numeric
    M <- data$area_index %>% unique %>% length
    mod_glmer <- lme4::glmer(formula = as.formula(formula), data = data,
                             family = "binomial")
    beta_1_t <- mod_glmer@beta
    sigma2.a_t <- as.numeric(VarCorr(mod_glmer)[1])
    a_i_vec <- ranef(mod_glmer)[[random_effect]][, 1]
    a_i.init_vec <- a_i_vec
    M <- length(a_i_vec)
    area_index_long <- data$area_index
    area_index_unique <- unique(area_index_long)
    a_i_ini_tbl <- tibble(a_i = a_i_vec, area_index = area_index_unique)
    a_i_long_tbl <- tibble(area_index = area_index_long)
    a_i_tbl <- left_join(a_i_long_tbl, a_i_ini_tbl, by = "area_index")
    re_form_pattern <- "\\((.*?)\\)"
    re_part_form <- str_match(formula, re_form_pattern)[1]
    fe_form <- str_remove_all(string = formula, pattern = coll(re_part_form)) %>%
      str_remove_all(pattern = "\\)") %>% str_remove_all(pattern = "\\(") %>%
      str_remove_all(pattern = "\\|") %>% str_squish()
    error_pattern <- "\\~(.*?)\\+"
    fixed_form <- ifelse(str_ends(string = fe_form, pattern = "\\+"),
                         str_sub(string = fe_form, start = 1L, end = (str_length(fe_form) -
                                                                        1)), str_replace(string = fe_form, pattern = error_pattern,
                                                                                         replace = "~")) %>% str_squish()
    response_var <- str_squish(str_split(fixed_form, "\\~",
                                         simplify = TRUE)[1])
    x <- model.matrix(as.formula(fixed_form), data = data)
    y <- data[[response_var]]
    if (length(unique(y)) != 2) {
      stop("too many values in response")
    }
    if (!is.numeric(y)) {
      y <- factor(y)
    }
    if (is.factor(y)) {
      y <- as.numeric(y == levels(factor(y))[2])
    }
    data$y <- y
    pred_part <- str_squish(str_split(fixed_form, "\\~",
                                      simplify = TRUE)[2])
    beta1.init <- rep(0, ncol(x))
    a_i.old_vec_long <- rep(0, length(y))
    change <- Inf
    beta1.old <- beta1.init
    a_i.old_vec <- a_i.init_vec
    outer_change <- Inf
    a_i.old_vec <- rep(0, M)
    cnt <- 0
    mu_hat <- 0
    tau2_hat <- 1
    params <- c(0.2, 2)
    while (outer_change > tol & cnt <= 20) {
      change <- Inf
      while (change > tol) {
        x_beta1 <- x %*% beta1.old %>% as.numeric
        eta <- x_beta1 + a_i.old_vec_long
        y.hat <- hc(eta)
        h.prime_eta <- ifelse(near(y.hat * (1 - y.hat), 0),
                              1e-04, y.hat * (1 - y.hat))
        z <- (y - y.hat)/h.prime_eta
        temp_tbl <- data
        temp_tbl$z <- z
        wls_form <- paste("z~", pred_part)
        beta1.new <- beta1.old + lm(formula = as.formula(wls_form),
                                    weights = h.prime_eta, data = data)$coef
        change <- sqrt(sum((beta1.new - beta1.old)^2))
        beta1.old <- beta1.new
      }
      beta_1.t <- beta1.old
      wk_tbl <- data %>% rename(a_index = area_index)
      wk_tbl$w_resid <- y - x %*% beta_1.t %>% as.numeric
      wk_tbl <- wk_tbl %>% mutate(pi_w = 1 - hc(w_resid))
      a_i.new_vec <- rep(NA, length(a_i.old_vec))
      v_i_vec <- rep(NA, length(a_i.old_vec))
      for (m in 1:M) {
        wk_tbl_i <- wk_tbl %>% filter(a_index == m)
        x_i <- model.matrix(as.formula(fixed_form), data = wk_tbl_i)
        y_i <- wk_tbl_i$y
        if (length(unique(y_i))>1){a_i.old <- a_i.old_vec[m]
        change <- Inf
        inner_cnt <- 0
        while (change > tol & inner_cnt <= 100) {
          x_beta1 <- x_i %*% beta1.old
          eta <- x_beta1 + a_i.old %>% as.numeric
          y.hat_i <- hc(eta)
          h.prime_eta_i <- ifelse(near(y.hat_i * (1 - y.hat_i),
                                       0), 1e-04, y.hat_i * (1 - y.hat_i))
          z_i <- (y_i - y.hat_i)/h.prime_eta_i
          a_i.new <- a_i.old + lm(z_i ~ 1, weights = h.prime_eta_i)$coef %>%
            as.numeric
          change <- sqrt(sum((a_i.new - a_i.old)^2))
          a_i.old <- a_i.new
          inner_cnt <- inner_cnt + 1
          a_i.new_vec[m] <- a_i.old
          v_i_vec[m] <- 1/sum(h.prime_eta_i)
        }
        }else{
          a_i.new_vec[m] <- a_i.old
          v_i_vec[m] <- 1e6
        }
      }
      V <- diag(v_i_vec)

      spat_coord <- data %>% dplyr::select(all_of(coordinates)) %>%
        unique
      D <- plgp::distance(spat_coord)
      K <- exp(-D/params[1]) + diag(params[2], nrow(spat_coord))
      K_inv <- solve(K)
      mu_hat <- drop(t(rep(1, nrow(K_inv))) %*% K_inv %*% a_i.new_vec/sum(K_inv))
      a_star <- drop(mu_hat + solve(tau2_hat * K + V) %*% (tau2_hat *
                                                             K) %*% (a_i.new_vec - mu_hat))
      V_star <- tau2_hat * K - tau2_hat * K %*% solve(tau2_hat *
                                                        K + V) %*% (tau2_hat * K)
      tau2_hat <- drop(sum(diag(solve(tau2_hat * K) %*% V_star)) +
                         t(a_star) %*% solve(tau2_hat * K) %*% a_star)/nrow(K)
      out <- optim(c(0.1, 0.01), spatial_lik, method = "L-BFGS-B",
                   lower = 1e-07, upper = c(10, 10), D = D, Y = a_star,
                   mu = mu_hat)
      params <- out$par
      a_i_ini_tbl_new <- tibble(a_i = a_star, area_index = unique(data$area_index))
      a_i_long_tbl_new <- tibble(area_index = data$area_index)
      a_i_tbl_new <- left_join(a_i_long_tbl_new, a_i_ini_tbl_new,
                               by = "area_index")
      a_i.old_vec_long <- a_i_tbl_new$a_i
      outer_change <- mean((a_star - a_i.old_vec)^2)
      # plot(a_i.old_vec, a_star)
      a_i.old_vec <- a_star
      cnt <- cnt + 1
      cat("...EM iteration:", cnt, "\n")
    }

    outcome <- list(cnt = cnt, area_tbl = a_i_tbl_new, area_re = a_i.old_vec,
                    params = c(out$par), mu_hat = mu_hat, spat_coord = spat_coord,
                    fomula = formula, V = v_i_vec, beta_hat = beta1.new,
                    tau2_hat = tau2_hat, D = D, coordinates = coordinates)

  return(outcome)
}

#' spatial prediction
#'
#' to calculate the spatial random effect for a new data set.
#'
#' @param em_output the outcome from the multilevel_EM function.
#'
#' @param test_set the test set (or new data set). The data set *should* have the predictor columns as well as the coordinate columns, all with the same names with the training data set.
#'
#' @return predicted spatial random effect at new data set

spatial_pred <- function(em_output, test_set){
  coordinates <- em_output$coordinates
  spat_coord <- em_output$spat_coord
  new_site <- test_set %>% dplyr::select(all_of(coordinates))
  if (is.null(nrow(new_site))){ # to make sure it works for a single site
    new_site <- as.data.frame(matrix(new_site, nrow = 1))
  }
  par <- em_output$params
  D_new <- plgp::distance(new_site)
  D <-  em_output$D
  mu_hat <- em_output$mu_hat
  a_star <- em_output$area_re
  K <- exp(- D/par[1]) + par[2]*diag(nrow(spat_coord))
  spat_coord <- em_output$spat_coord
  K_new <- exp(-D_new/par[1]) + par[2]*diag(ncol(D_new))
  D_new_train <- plgp::distance(new_site, spat_coord)
  K_new_train <- exp(-D_new_train/par[1])
  K_inv <- solve(K)
  a_hat <- mu_hat + K_new_train %*% K_inv %*% (a_star - mu_hat)
  a_hat <- drop(a_hat)
  return(a_hat)
}


#' covariate prediction
#'
#' to calculate the spatial random effect for a new data set.
#'
#' @param em_output the outcome from the multilevel_EM function.
#'
#' @param test_set the test set (or new data set). The data set *should* have the predictor columns as well as the coordinate columns, all with the same names with the training data set.
#'
#' @return predicted covariates parts Xb at new data set

covariate_pred <- function(em_output, test_set){
  beta_hat <- em_output$beta_hat
  formula <- em_output$fomula
  re_form_pattern <- "\\((.*?)\\)"
  re_part_form <- str_match(formula, re_form_pattern)[1]
  fe_form <-
    str_remove_all(string = formula, pattern = coll(re_part_form)) %>%
    # coll is to force the matching the "full expression", not as a regex
    str_remove_all(pattern = "\\)") %>%
    str_remove_all(pattern = "\\(") %>%
    str_remove_all(pattern = "\\|") %>%
    str_squish()

  error_pattern <- "\\~(.*?)\\+"
  fixed_form <-
    ifelse(str_ends(string = fe_form, pattern = "\\+"),
           str_sub(string = fe_form, start = 1L, end = (str_length(fe_form)-1)),
           str_replace(string = fe_form, pattern = error_pattern, replace = "~")) %>%
    str_squish()

  mm_form <- str_squish(str_split(fixed_form, "\\~", simplify = TRUE)[2])

  x <- model.matrix(as.formula(paste("~", mm_form)), data = test_set)
  out <- as.numeric(x %*% beta_hat)
  return(out)
}

#' main prediction
#'
#' a wrapper to calculate the spatial random effect for a new data set.
#'
#' @param em_output the outcome from the multilevel_EM function.
#'
#' @param test_set the test set (or new data set). The data set *should* have the predictor columns as well as the coordinate columns, all with the same names with the training data set.
#'
#' @return predicted outcome combining covariates parts Xb and spatial random effects at new data set
#'
#' @export

predict_em <- function(em_output, test_set){
  yhat <- covariate_pred(em_output, test_set)
  spat_hat <- spatial_pred(em_output, test_set)
  combined_hat <- yhat + spat_hat
  phat <- hc(combined_hat)
  return(phat)
}



#' simulate the spatial random effects for pixel within each area
#' internal wrapper
#'
#' @param em_output the outcome from the multilevel_EM function.
#'
#' @param test_set the test set (or new data set).
#' The data set *should* have the predictor columns as well
#' as the coordinate columns,
#' all with the same names with the training data set.
#'
#' @size number of simulated spatial random effect of size 100 (100 is default)
#' at new data set
#'
#'
#'
simulate_p <- function(em_output, test_set, size = 100){
  coordinates <- em_output$coordinates
  spat_coord <- em_output$spat_coord
  tau2_hat <- em_output$tau2_hat
  new_site <- test_set %>% dplyr::select(all_of(coordinates))
  if (is.null(nrow(new_site))){ # to make sure it works for a single site
    new_site <- as.data.frame(matrix(new_site, nrow = 1))
  }

  par <- em_output$params
  D_new <- plgp::distance(new_site)
  D <-  em_output$D
  mu_hat <- em_output$mu_hat
  a_star <- em_output$area_re
  K <- exp(- D/par[1]) + par[2]*diag(nrow(spat_coord))

  K_new <- exp(-D_new/par[1]) + par[2]*diag(ncol(D_new))
  D_new_train <- plgp::distance(new_site, spat_coord)
  K_new_train <- exp(-D_new_train/par[1])
  K_inv <- solve(K)
  a_hat <- mu_hat + K_new_train %*% K_inv %*% (a_star - mu_hat)
  a_hat <- drop(a_hat)

  Sigmap <- tau2_hat*(K_new - K_new_train %*% K_inv %*% t(K_new_train))

  D_new <- plgp::distance(new_site)
  K_new <- exp(-D_new/par[1]) + diag(par[2], nrow(D_new))

  Sigmap <- tau2_hat*(K_new - K_new_train %*% K_inv %*% t(K_new_train))

  Sigma.int <- tau2_hat*(exp(-D_new) + diag(par[2], nrow(D_new))  - (K_new_train %*% K_inv %*% t(K_new_train)))

  ## to ensure A = LL'; not an ideal fix but does the job..
  S <- svd(Sigma.int)
  Sigma.int <- S$u %*% diag(S$d) %*% t(S$u) # enforce the symmetry

  a_k <- rmvnorm(size, a_hat, Sigma.int)

  # routine to calculate the MSE part 1
  yhat <- covariate_pred(em_output, test_set)
  rep_yhat_mat <- matrix(rep(yhat, size), nrow = size, byrow = TRUE)
  combined_mat <- a_k + rep_yhat_mat

  simulated_p <- 1 / (1 + exp(-combined_mat))

  return(simulated_p)
}


#' get the small area estimates for each area
#'
#' @param em_output the outcome from the multilevel_EM function.
#'
#' @param level_2_geom the area shapefile-simple feature
#'#'
#' @param coord_info coordinates for each DHS survey -- must be indexed by column
#' named DHSCLUST
#'
#' @param frac how much cluster we want to choose.
#'
#' @param resampling resampling rate.
#'
#' @counter to track the computation progression -- logical; TRUE or FALSE
#'
#'
#'

produce_m_v_sae <- function(em_output, level_2_geom, coord_info,
                            frac = 0.1, resampling = 0.3, counter = FALSE){
  output <- level_2_geom

  cluster_coord <- coord_info %>% select(-DHSCLUST)

  output$m <- NA
  output$v <- NA

  grid_list <- list()
  for (admin_id in 1:nrow(level_2_geom)){

    gr_admin_id <- level_2_geom[admin_id,] %>%
      st_make_grid(n = c(7, 7), what = "centers") %>%     # get the box shape grid of 7 by 7
      st_intersection(level_2_geom[admin_id,]) # get the part intersects with the admin area

    # it stores in grid_list
    grid_list[[admin_id]] <- gr_admin_id
    grid_coord_admin_id <- st_coordinates(gr_admin_id)

    # coordinates should be converted to a matrix; em model only takes the matrix
    grid_coord_admin_id <- as_tibble(grid_coord_admin_id[,1:2])
    names(grid_coord_admin_id) <- names(em_output$spat_coord)

    # prediction for each grid point within the small area

    predicted_list <- list()
    simulated_out <- list()
    for (c_i in 1:nrow(grid_coord_admin_id)){
      current_point <- matrix(grid_coord_admin_id[c_i,], ncol = 2)
      names(current_point) <- names(em_output$spat_coord)
      cluster_dist_temp <- coord_info
      cluster_dist_temp$p_rank <-
        rdist(as.matrix(cluster_coord), current_point) %>% percent_rank()

      # we choose the points that belong to nearest 10%
      chosen_coord <- cluster_dist_temp %>% filter(p_rank <= frac)

      # now we choose the survey portion
      chosen_tbl <- valid_set_with_sp %>% filter(DHSCLUST %in% chosen_coord$DHSCLUST)

      # resampling
      temp_tbl <- chosen_tbl %>% sample_frac(resampling)
      # prediction
      predicted_list[[c_i]] <- predict_em(em_output, test_set = temp_tbl)
      simulated_out[[c_i]] <- simulate_p(em_output, test_set = temp_tbl, size = 100)
    }

    point_mean <- sapply(simulated_out, mean)
    v1 <- var(point_mean) / length(point_mean)
    v2 <- sapply(sapply(simulated_out, function(x) apply(x, 2, mean), simplify = FALSE), var)

    # mean estimates and MSE
    output$m[admin_id] <- do.call("c", predicted_list) %>% mean
    output$v[admin_id] <- v1 + mean(v2)

    if (counter & admin_id%%5==0) cat("---", admin_id, "th area /", nrow(level_2_geom), " -- \n")
  }

  return(output)
}
