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
  ll <- -(n/2)*log(t(Y - mu) %*% K_inv %*% Y - mu) - (1/2)*ldetK
  return(-ll)
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
  random_effect <- str_squish(str_match(formula, re_pattern)[,2])

  # create a col named "area_index" for convenience
  if (!is.factor(data[[random_effect]])){stop("area index must be a factor")}
  data$area_index <- droplevels(data[[random_effect]]) %>%
    as.numeric
  M <- data$area_index %>%
    unique %>% length


  # initial fit using glmer
  mod_glmer <- lme4::glmer(formula = as.formula(formula), data = data, family = "binomial")

  beta_1_t <- mod_glmer@beta
  sigma2.a_t <- as.numeric(VarCorr(mod_glmer)[1]) # estimated variance of area ranef
  a_i_vec <- ranef(mod_glmer)[[random_effect]][,1]
  a_i.init_vec <- a_i_vec

  M <- length(a_i_vec)

  area_index_long <- data$area_index
  area_index_unique <- unique(area_index_long)

  a_i_ini_tbl <- tibble(a_i =  a_i_vec, area_index = area_index_unique)
  a_i_long_tbl <- tibble(area_index = area_index_long)
  a_i_tbl <- left_join(a_i_long_tbl, a_i_ini_tbl, by = "area_index")

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

  response_var <- str_squish(str_split(fixed_form, "\\~", simplify = TRUE)[1])

  x <- model.matrix(as.formula(fixed_form), data = data)
  y <- data[[response_var]]

  if (length(unique(y))!= 2 ){stop("too many values in response")}
  if (!is.numeric(y)) {y <- factor(y) }
  if (is.factor(y)) {y <- as.numeric(y == levels(factor(y))[1]) }

  # this will be used in the m = 1, .. M loop
  data$y <- y

  pred_part <- str_squish(str_split(fixed_form, "\\~", simplify = TRUE)[2])

  # Setting up the inital params
  theta1.init <- rep(0, ncol(x))
  # here should be the current value of X beta_t + E(a_i| ...)
  a_i.old_vec_long <- rep(0, length(y))
  change <- Inf
  beta1.old <- theta1.init
  a_i.old_vec <- a_i.init_vec

  outer_change <- Inf
  a_i.old_vec <- rep(0, M)
  cnt <- 0

  # initialize
  mu_hat <- 0
  tau2_hat <- 1
  params <- c(0.2, 2)

  while (outer_change > tol & cnt <= 20){
    # M step in f1

    change <- Inf
    while(change > tol  ) { # M step for Q1 using IRWLS
      x_beta1 <- x %*% beta1.old %>% as.numeric
      # eta is a linear predictor with xb + a, a = area ranef, a ~ N(0, sigma2a_t)
      eta <- x_beta1 + a_i.old_vec_long # E(a_i | a_i_hat...) = a_i hat assuming normal
      y.hat <- hc(eta)
      h.prime_eta <- ifelse(near(y.hat * (1 - y.hat), 0), 1e-5, y.hat * (1 - y.hat))
      z <- (y - y.hat) / h.prime_eta
      temp_tbl <- data
      temp_tbl$z <- z
      wls_form <- paste("z~", pred_part)
      beta1.new <- beta1.old + lm(formula = as.formula(wls_form), weights = h.prime_eta, data = data)$coef  # WLS regression
      change <- sqrt(sum((beta1.new - beta1.old)^2))
      beta1.old <- beta1.new
      # cat(change, "\n")
      # M step associated with Q1 function
    } # at the end of this while loop, x_beta1[t+1] is obtained

    beta_1.t <- beta1.old
    wk_tbl <- data %>% rename(a_index = area_index)
    wk_tbl$w_resid <- y - x %*% beta_1.t %>% as.numeric # working residual
    wk_tbl <- wk_tbl %>% mutate(pi_w = 1 - hc(w_resid))


    # M step in f2
    a_i.new_vec <- rep(NA, length(a_i.old_vec))
    v_i_vec <- rep(NA, length(a_i.old_vec))

    for (m in 1:M){
      wk_tbl_i <- wk_tbl %>% filter(a_index == m) # for m-th area

      x_i <- model.matrix(as.formula(fixed_form), data = wk_tbl_i) # model matrix for m-th area
      y_i <- wk_tbl_i$y # response for m-th area

      a_i.old <- a_i.old_vec[m]
      change <- Inf
      while(change > tol) { # M step for Q1 using IRLS
        x_beta1 <- x_i %*% beta1.old
        # eta is a linear predictor with xb + a, a = area ranef, a ~ N(0, sigma2a_t)
        eta <- x_beta1 + a_i.old %>% as.numeric # E(a_i | a_i_hat...) = a_i hat assuming normal
        y.hat_i <- hc(eta)
        h.prime_eta_i <- ifelse(near(y.hat_i * (1 - y.hat_i), 0), 1e-5, y.hat_i * (1 - y.hat_i))
        z_i <- (y_i - y.hat_i) / h.prime_eta_i
        a_i.new <- a_i.old + lm(z_i ~ 1, weights = h.prime_eta_i)$coef %>% as.numeric # WLS regression


        change <- sqrt(sum((a_i.new - a_i.old)^2))
        a_i.old <- a_i.new
        # cat(change, "\n")
        # M step associated with Q1 function
      } # at the end of this while loop, a_i[t+1] is obtained
      a_i.new_vec[m] <- a_i.old
      v_i_vec[m] <- 1/sum(h.prime_eta_i)
    }

    V <- diag(v_i_vec)

    sigma_area_t <- mean(a_i.new_vec^2)

    spat_coord <- data %>%
      dplyr::select(all_of(coordinates)) %>%
      unique

    D <- plgp::distance(spat_coord)
    K <- exp(- D/params[1]) + diag(params[2], nrow(spat_coord))
    Ki <- solve(K)
    mu_hat <- drop(t(rep(1, nrow(Ki))) %*% Ki %*% a_i.new_vec / sum(Ki))

    a_star <- drop(mu_hat + solve(tau2_hat*K + V) %*% (tau2_hat*K) %*% (a_i.new_vec-mu_hat))
    V_star <- tau2_hat*K - tau2_hat*K %*% solve(tau2_hat*K + V) %*% (tau2_hat*K)

    tau2_hat <- drop(sum(diag( solve(tau2_hat*K) %*% V_star )) + t(a_star) %*% solve(tau2_hat*K) %*% a_star)/nrow(K)

    out <- optim(c(0.1, 0.01), spatial_lik, method="L-BFGS-B", lower=1e-7,
                 upper = c(10, 10), D = D, Y = a_star, mu = mu_hat)

    params <- out$par
    a_i_ini_tbl_new <- tibble(a_i =  a_star, area_index = unique(data$area_index))
    a_i_long_tbl_new <- tibble(area_index = data$area_index)
    a_i_tbl_new <- left_join(a_i_long_tbl_new, a_i_ini_tbl_new, by = "area_index")

    a_i.old_vec_long <- a_i_tbl_new$a_i
    outer_change <- mean((a_star-a_i.old_vec)^2)
    a_i.old_vec <- a_star


    cnt <- cnt + 1
    #cat(mean(abs(a_i.new_vec)), "\n")
    #cat(outer_change, "\n")
    #cat(cnt, "\n")

  }

  outcome <- list(cnt = cnt, area_tbl = a_i_tbl_new, area_re = a_i.old_vec,
                  params = c(out$par), mu_hat = mu_hat,
                  spat_coord = spat_coord,
                  fomula = formula,
                  V = v_i_vec, beta_hat = beta1.new,
                  tau2_hat = tau2_hat, D = D,
                  coordinates = coordinates)
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
  new_site <- test_set %>% dplyr::select(all_of(coordinates))
  if (is.null(nrow(new_site))){ # to make sure it works for a single site
    new_site <- as.data.frame(matrix(new_site, nrow = 1))
  }
  par <- em_output$params
  DXX <- plgp::distance(new_site)
  D <-  em_output$D
  mu_hat <- em_output$mu_hat
  a_star <- em_output$area_re
  K <- exp(- D/par[1]) + par[2]*diag(nrow(spat_coord))
  spat_coord <- em_output$spat_coord
  KXX <- exp(-DXX/par[1]) + par[2]*diag(ncol(DXX))
  DX <- plgp::distance(new_site, spat_coord)
  KX <- exp(-DX/par[1])
  Ki <- solve(K)
  a_hat <- mu_hat + KX %*% Ki %*% (a_star - mu_hat)
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
  return(combined_hat)
}




