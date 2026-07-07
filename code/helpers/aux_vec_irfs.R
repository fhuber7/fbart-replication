remove_outliers_draws <- function(X, k = 1.5, drop_rows = FALSE) {
  X <- as.matrix(X)
  
  X_clean <- apply(X, 2, function(col) {
    q1 <- quantile(col, 0.25, na.rm = TRUE)
    q3 <- quantile(col, 0.75, na.rm = TRUE)
    iqr <- q3 - q1
    lo <- q1 - k * iqr
    hi <- q3 + k * iqr
    
    col[col < lo | col > hi] <- NA
    col
  })
  
  if (drop_rows) {
    keep <- complete.cases(X_clean)
    X_clean <- X_clean[keep, , drop = FALSE]
  }
  
  X_clean
}

update_pi_power <- function(Lambda, psi,
                            kappa = 2,
                            a_pi = 3,
                            b_pi = 0.03,
                            tiny = 1e-12) {
  # Lambda : M x Q factor loading matrix
  # psi    : M x Q local variance matrix (tau_i * phi_ij)
  # kappa  : integer controlling deterministic column decay
  # c      : scale parameter of Inv-Gamma prior on pi
  #
  # Prior:     pi ~ Inv-Gamma(1/2, c)
  # Posterior: pi | . ~ Inv-Gamma((MQ + 1)/2, c + S/2)
  
  if (!is.matrix(Lambda) || !is.matrix(psi))
    stop("Lambda and psi must be matrices.")
  
  if (any(dim(Lambda) != dim(psi)))
    stop("Lambda and psi must have the same dimensions.")
  
  M <- nrow(Lambda)
  Q <- ncol(Lambda)
  
  ## Defensive cleaning
  psi_clean <- psi
  psi_clean[!is.finite(psi_clean) | psi_clean <= 0] <- tiny
  
  Lambda_clean <- Lambda
  Lambda_clean[!is.finite(Lambda_clean)] <- 0
  
  ## Compute sufficient statistic S
  S <- 0
  for (j in seq_len(Q)) {
    S <- S + (j^kappa) * sum(Lambda_clean[, j]^2 / psi_clean[, j])
  }
  if (!is.finite(S) || S <= 0) S <- tiny
  
  ## Posterior parameters
  shape_post <- (M * Q) / 2 + a_pi
  scale_post <- b_pi + 0.5 * S
  
  ## Sample pi from Inv-Gamma
  # If X ~ Gamma(shape, rate), then 1/X ~ Inv-Gamma(shape, scale = 1/rate)
  pi_draw <- 1 / rgamma(1, shape = shape_post, rate = scale_post)
  
  if (!is.finite(pi_draw) || pi_draw <= 0)
    pi_draw <- tiny
  
  ## Column-wise shrinkage factors
  pi_vec <- as.numeric(pi_draw / (seq_len(Q)^kappa))
  
  list(
    pi     = as.numeric(pi_draw),
    pi_vec = pi_vec,
    S      = S,
    shape  = shape_post,
    scale  = scale_post
  )
}

update_pi_power_G <- function(Lambda, psi,
                            kappa = 2,
                            c = 1,
                            tiny = 1e-12) {
  # Lambda : M x Q factor loading matrix
  # psi    : M x Q local variance matrix (horseshoe output)
  # kappa  : integer controlling deterministic column decay
  # c      : variance of half-Normal prior on sqrt(pi)
  #
  # Prior: pi ~ Gamma(1/2, 1/(2c))
  # Posterior: pi | . ~ GIG(1/2 - MQ/2, S, 1/c)

  if (!requireNamespace("GIGrvg", quietly = TRUE)) {
    stop("Package 'GIGrvg' is required.")
  }

  if (!is.matrix(Lambda) || !is.matrix(psi))
    stop("Lambda and psi must be matrices.")

  if (any(dim(Lambda) != dim(psi)))
    stop("Lambda and psi must have the same dimensions.")

  M <- nrow(Lambda)
  Q <- ncol(Lambda)

  # Defensive cleaning
  psi_clean <- psi
  psi_clean[!is.finite(psi_clean) | psi_clean <= 0] <- tiny

  Lambda_clean <- Lambda
  Lambda_clean[!is.finite(Lambda_clean)] <- 0

  # Compute sufficient statistic S
  S <- 0
  for (j in seq_len(Q)) {
    S <- S + (j^kappa) * sum(Lambda_clean[, j]^2 / psi_clean[, j])
  }
  if (!is.finite(S) || S <= 0) S <- tiny

  # GIG parameters
  lambda_gig <- 0.5 - 0.5 * M * Q
  chi_gig    <- S
  psi_gig    <- 1 / c

  # Sample pi
  pi_draw <- GIGrvg::rgig(
    n      = 1,
    lambda = lambda_gig,
    chi    = chi_gig,
    psi    = psi_gig
  )

  if (!is.finite(pi_draw) || pi_draw <= 0)
    pi_draw <- tiny

  # Column-wise shrinkage factors
  pi_vec <- as.numeric(pi_draw / (seq_len(Q)^kappa))

  list(
    pi        = as.numeric(pi_draw),
    pi_vec    = pi_vec,
    S         = S,
    lambda_gig = lambda_gig
  )
}

# update_pi_power <- function(Lambda, psi,
#                             kappa = 2,
#                             a_phi = 1,
#                             b_phi = 1,
#                             tiny = 1e-12) {
#   # Lambda : M x Q loading matrix
#   # psi    : M x Q local variance matrix
#   # kappa  : column shrinkage exponent
#   # a_phi  : Gamma shape for precision phi = pi^{-1}
#   # b_phi  : Gamma rate  for precision phi = pi^{-1}
#   #
#   # Returns:
#   #   pi      : global variance draw
#   #   pi_vec  : column-wise variances pi / j^kappa
#   
#   if (!is.matrix(Lambda) || !is.matrix(psi))
#     stop("Lambda and psi must be matrices.")
#   if (any(dim(Lambda) != dim(psi)))
#     stop("Lambda and psi must have the same dimensions.")
#   
#   M <- nrow(Lambda)
#   Q <- ncol(Lambda)
#   
#   # defensive cleaning
#   psi_clean <- psi
#   psi_clean[!is.finite(psi_clean) | psi_clean <= 0] <- tiny
#   
#   Lambda_clean <- Lambda
#   Lambda_clean[!is.finite(Lambda_clean)] <- 0
#   
#   # sufficient statistic
#   S <- 0
#   for (j in seq_len(Q)) {
#     S <- S + (j^kappa) * sum(Lambda_clean[, j]^2 / psi_clean[, j])
#   }
#   if (!is.finite(S) || S < 0) S <- 0
#   
#   # posterior for precision phi = pi^{-1}
#   shape_post <- a_phi + 0.5 * M * Q
#   rate_post  <- b_phi + 0.5 * S
#   
#   phi_draw <- rgamma(1, shape = shape_post, rate = rate_post)
#   if (!is.finite(phi_draw) || phi_draw <= 0)
#     phi_draw <- tiny
#   
#   # transform back to variance
#   pi_draw <- 1 / phi_draw
#   
#   # column-wise variances
#   pi_vec <- as.numeric(pi_draw / (seq_len(Q)^kappa))
#   
#   list(
#     pi     = pi_draw,
#     pi_vec = pi_vec
#   )
# }



update_column_hs <- function(Lambda,
                             pi_prev,
                             delta_prev,
                             a_delta = 1,
                             b_delta = 1) {
  # Lambda : M x Q matrix of standardized loadings (Lambda.mu / sqrt(psi))
  # pi_prev: length-Q vector (previous column scales)
  # delta_prev: length-Q vector (previous increments)
  #
  # Prior:
  #   delta_h ~ half-Cauchy(0,1)
  #   pi_h = prod_{l <= h} delta_l
  #
  # Returns:
  #   list(delta, pi)
  
  M <- nrow(Lambda)
  Q <- ncol(Lambda)
  
  if (length(delta_prev) != Q) delta_prev <- rep(1, Q)
  if (length(pi_prev) != Q) pi_prev <- cumprod(delta_prev)
  
  delta <- numeric(Q)
  pi <- numeric(Q)
  
  # sufficient statistics
  col_ss <- colSums(Lambda^2)
  
  for (h in 1:Q) {
    # effective variance contribution of delta_h
    denom <- prod(delta_prev[1:h], na.rm = TRUE)
    denom <- max(denom, 1e-8)
    
    shape <- (M / 2) + a_delta
    rate  <- (col_ss[h] / (2 * denom)) + b_delta
    
    # sample delta_h^{-1} (inverse-gamma representation of half-Cauchy)
    delta_inv <- rgamma(1, shape = shape, rate = rate)
    delta[h] <- 1 / sqrt(delta_inv)
    
    # numerical guard
    delta[h] <- min(max(delta[h], 1e-3), 1e3)
  }
  
  pi[1] <- delta[1]
  if (Q > 1) {
    for (h in 2:Q) pi[h] <- pi[h - 1] * delta[h]
  }
  
  # normalize for numerical stability
  pi <- pi / max(pi)
  
  list(delta = delta, pi = pi)
}

update_CUSP <- function(Lambda,
                               w_prev,
                               theta_inv_prev = NULL, # not used inside but kept for API compatibility
                               alpha = 1,
                               theta_inf = 1e-4,
                               a_theta = 2,
                               b_theta = 2) {
  # Lambda: p x H loading matrix
  # w_prev: length-H vector of previous stick-breaking weights (will be normalized)
  # returns: list(z, v, w, theta_inv)
  #
  # Robustified: uses log-likelihoods, handles NA/Inf, stable normalization.
  
  if (!is.matrix(Lambda)) stop("Lambda must be a p x H matrix.")
  p <- nrow(Lambda)
  H <- ncol(Lambda)
  
  # defensive w_prev
  if (length(w_prev) != H) {
    warning("w_prev length != ncol(Lambda). Resetting to uniform.")
    w_prev <- rep(1 / H, H)
  }
  # replace NA/NaN and negative entries in w_prev
  w_prev[is.na(w_prev) | is.nan(w_prev) | (w_prev < 0)] <- 0
  if (sum(w_prev) <= 0) w_prev <- rep(1 / H, H) else w_prev <- w_prev / sum(w_prev)
  log_w_prev <- log(w_prev)
  
  # Precompute log-likelihoods for each column under spike and slab
  # spike: Gaussian N(0, theta_inf)
  # slab: Student-t multivariate; using dmvt with log=TRUE (from mvtnorm)
  log_lhd_spike <- rep(-Inf, H)
  log_lhd_slab  <- rep(-Inf, H)
  
  for (h in seq_len(H)) {
    # spike: elementwise normals, sum log-densities
    ll_spike <- tryCatch({
      sum(dnorm(Lambda[, h], mean = 0, sd = sqrt(theta_inf), log = TRUE), na.rm = FALSE)
    }, error = function(e) { -Inf })
    if (!is.finite(ll_spike)) ll_spike <- -Inf
    log_lhd_spike[h] <- ll_spike
    
    # slab: multivariate t; dmvt from mvtnorm supports log=TRUE
    # If dmvt not present in namespace, try extra packages will be required by user.
    ll_slab <- tryCatch({
      dmvt(x = Lambda[, h],
           delta = rep(0, p),
           sigma = (b_theta / a_theta) * diag(p),
           df = 2 * a_theta,
           log = TRUE)
    }, error = function(e) {
      # fallback: approximate slab with multivariate normal log-density (less ideal)
      sum(dnorm(Lambda[, h], mean = 0, sd = sqrt((b_theta / a_theta)), log = TRUE), na.rm = FALSE)
    })
    if (!is.finite(ll_slab)) ll_slab <- -Inf
    log_lhd_slab[h] <- ll_slab
  }
  
  ## --- 4) sample z (robustly) ---
  z <- integer(H)
  for (h in seq_len(H)) {
    # construct log-prob vector of length H:
    # positions 1..h: use spike(log_lhd_spike[h])
    # positions (h+1)..H: use slab(log_lhd_slab[h])
    log_comp <- c(rep(log_lhd_spike[h], h),
                  rep(log_lhd_slab[h], max(0, H - h)))
    # make length exactly H (in case h==H or h==0)
    if (length(log_comp) < H) log_comp <- c(log_comp, rep(-Inf, H - length(log_comp)))
    
    # total log-prob = log w_prev + log_comp
    log_prob_raw <- log_w_prev + log_comp
    # replace non-finite with -Inf
    log_prob_raw[!is.finite(log_prob_raw)] <- -Inf
    
    # stable normalization via log-sum-exp
    m <- max(log_prob_raw)
    if (!is.finite(m)) {
      # everything is -Inf: fallback deterministic assignment (last position)
      prob_vec <- rep(0, H); prob_vec[H] <- 1
    } else {
      exps <- exp(log_prob_raw - m)
      sum_exps <- sum(exps)
      if (!is.finite(sum_exps) || sum_exps <= .Machine$double.eps) {
        # degenerate underflow; pick the max entry deterministically
        idx_max <- which.max(log_prob_raw)
        prob_vec <- rep(0, H); prob_vec[idx_max] <- 1
      } else {
        prob_vec <- exps / sum_exps
      }
    }
    # sample index robustly (if rounding errors make sum !=1, sample handles it)
    z[h] <- sample.int(H, size = 1, prob = prob_vec)
  }
  
  ## --- 5) sample v and update w (stick-breaking) ---
  v <- numeric(H)
  for (h in 1:(H - 1)) {
    v[h] <- rbeta(1,
                  shape1 = 1 + sum(z == h),
                  shape2 = alpha + sum(z > h))
  }
  v[H] <- 1
  w <- numeric(H)
  w[1] <- v[1]
  if (H > 1) {
    for (h in 2:H) {
      w[h] <- v[h] * prod(1 - v[1:(h - 1)])
    }
  }
  # defensive normalization
  if (sum(w) <= 0 || any(!is.finite(w))) {
    w <- rep(1 / H, H)
  } else {
    w <- w / sum(w)
  }
  
  ## --- 6) sample theta^{-1} ---
  theta_inv <- numeric(H)
  for (h in seq_len(H)) {
    if (z[h] <= h) {
      theta_inv[h] <- 1 / theta_inf
    } else {
      # robust scalar quadratic form
      quad <- as.numeric(crossprod(Lambda[, h]))
      rate_post <- b_theta + 0.5 * quad
      # guard against non-finite rate
      if (!is.finite(rate_post) || rate_post <= 0) rate_post <- b_theta + 0.5 * sum((Lambda[, h])^2, na.rm = TRUE) + 1e-8
      theta_inv[h] <- rgamma(1, shape = a_theta + 0.5 * p, rate = rate_post)
    }
  }
  
  list(z = z, v = v, w = w, theta_inv = theta_inv)
}


  
  ## ---
  
signident <- function(x, method = "maximin", implementation = 3) {
  
  # method <- "maximin"
  # implementation <- 3
  # x <- list("facload"=aperm(lambda.array,c(2,3,1)), "fac"=aperm(factors.array,c(3,2,1)))
  
  
  if (method != "diagonal" & method != "maximin")
    stop("Argument 'method' must either be 'diagonal' or 'maximin'.")
  r <- dim(x$facload)[2]
  ftpoints <- dim(x$fac)[2]
  
  if (r == 0) {
    x$identifier <- matrix(NA_real_, nrow = 0, ncol = 2)
  } else {
    identifier <- 1:r  # if method == "diagonal", overwritten otherwise
    distance <- rep(NA_real_, r)
    
    for (i in 1:r) {
      faccol <- matrix(x$facload[,i,,drop=FALSE], nrow = nrow(x$facload))
      if (method == "maximin") { # for each factor, look for the series where the
        # minimum absolute loadings are biggest
        identifier[i] <- which.max(apply(abs(faccol), 1, min))
      }
      
      distance[i] <- max(apply(abs(faccol), 1, min))
      mysig <- sign(faccol[identifier[i],])
      x$facload[,i,] <- t(t(faccol) * mysig)
      if (implementation == 1) {
        x$fac[i,,] <- x$fac[i,,] * rep(mysig, each = ftpoints)
      } else if (implementation == 2) {
        for (j in seq(along = mysig)) {
          x$fac[i,,j] <- x$fac[i,,j] * mysig[j]
        }
      } else if (implementation == 3) {
        for (j in 1:ftpoints) {
          x$fac[i,j,] <- x$fac[i,j,] * mysig
        }
      } else stop("Err0r.")
    }
    x$identifier <- matrix(c(identifier, distance), ncol = 2)
  }
  
  colnames(x$identifier) <- c("identifier", "distance")
  
  x
}

log_marginal_lik_many <- function(y, Xmat, sigma2, tau2,
                                  include_null = FALSE) {
  # y: T x 1 vector
  # Xmat: T x K matrix (each column a candidate regressor)
  # Returns: vector of log marginal likelihoods
  #           [null, X1, X2, ..., XK] if include_null = TRUE
  
  T <- length(y)
  K <- if (is.null(Xmat)) 0 else ncol(Xmat)
  
  yty <- sum(y^2)
  
  # ---- null model ----
  loglik_null <- -0.5 * (
    T * log(2 * pi) +
      T * log(sigma2) +
      yty / sigma2
  )
  
  if (K == 0) {
    return(loglik_null)
  }
  
  # ---- single-regressor models ----
  XtX <- colSums(Xmat^2)          # K x 1
  ytX <- as.numeric(crossprod(y, Xmat))  # K x 1
  
  c <- 1 + (tau2 / sigma2) * XtX
  
  logdet <- T * log(sigma2) + log(c)
  
  quad <- (yty / sigma2) -
    (tau2 / sigma2^2) * (ytX^2 / c)
  
  loglik_X <- -0.5 * (
    T * log(2 * pi) + logdet + quad
  )
  
  if (include_null) {
    return(c(loglik_null, loglik_X))
  } else {
    return(loglik_X)
  }
}
# get posteriors for the horseshoe prior (see Makalic & Schmidt, 2015)
# prior moments for Minnesota prior
prior_a_minn <- function(a_bar_1, a_bar_2, a_bar_3){
  # a_bar_1 <- 0.04
  # a_bar_2 <- 0.0016
  # a_bar_3 <- 10
  p.1 <- dgamma(a_bar_1, 4 , 100, log=TRUE) # Prior mean of 0.04; Prior variance .4/10^2
  p.2 <- dgamma(a_bar_2, 0.16, 100, log=TRUE) # Prior mean of 0.0016; 
  p.3 <- dgamma(a_bar_3, 0.01, 0.01, log=TRUE)
  
  return(p.1+p.2+p.3)
}
get_V <- function(a_bar_1,a_bar_2, a_bar_3 , ind1,sigma_sq1,p1,M1,K1){
  #a_bar_1;a_bar_2; ind1<-ind; sigma_sq1=sigma.sd^2; p; M1=M; K1=K.var
  V_i <- matrix(0,K1,M1)
  #this double loop fills the prior covariance matrix
  for (i in 1:M1){ #for each i equation
    for (j in 1:K1){ #for each variable on the rhs
      if (j > (M1*p)) {
        V_i[j,i] <- a_bar_3 * sigma_sq1[i,1] #variance on constant, trend, dummies and on ex variables         
      }    
      else if (any(j==ind1[i,])) {
        ll <- which(ind1[i,]==j)
        V_i[j,i] <- a_bar_1/(ll^2) #variance on own lags
      }else{
        ll <- which(ind1==j,arr.ind=TRUE)[2]
        kj <- ind1[which(ind1==j,arr.ind=TRUE)[1],1]
        V_i[j,i] <- (a_bar_2*sigma_sq1[i,1])/((ll^2)*sigma_sq1[kj,1])
      }
    }
  }  
  return(V_i)
}
get.hs <- function(bdraw,lambda.hs,nu.hs,tau.hs,zeta.hs){
  k <- length(bdraw)
  if (is.na(tau.hs)){
    tau.hs <- 1   
  }else{
    tau.hs <- invgamma::rinvgamma(1,shape=(k+1)/2,rate=1/zeta.hs+sum(bdraw^2/lambda.hs)/2) 
  }
  
  lambda.hs <- invgamma::rinvgamma(k,shape=1,rate=1/nu.hs+bdraw^2/(2*tau.hs))
  
  nu.hs <- invgamma::rinvgamma(k,shape=1,rate=1+1/lambda.hs)
  zeta.hs <- invgamma::rinvgamma(1,shape=1,rate=1+1/tau.hs)
  
  ret <- list("psi"=(lambda.hs*tau.hs),"lambda"=lambda.hs,"tau"=tau.hs,"nu"=nu.hs,"zeta"=zeta.hs)
  return(ret)
}

fred.trans <- function(x, trans){
  X <- matrix(NA, nrow(x), ncol(x), dimnames = list(NULL, colnames(x)))
  for (i in 1:ncol(x)){
    if (trans[[i]]==1){
      X[, i] <- x[, i]
    }else if (trans[[i]]==2){
      X[2:nrow(X), i] <- diff(x[,i])
    }else if (trans[[i]]==3){
      X[3:nrow(X), i] <- diff(diff(x[,i]))
    }else if (trans[[i]]==4){
      X[, i] <- log(x[, i])
    }else if(trans[[i]]==5){
      X[2:nrow(X), i] <- diff(log(x[, i]))
    }else if (trans[[i]] == 6){
      X[3:nrow(X), i] <- diff(diff(log(x[, i])))
    }else if (trans[[i]]==7){
      time_series <- x[,i]
      X[13:nrow(X), i] <- ((time_series[13:length(time_series)]) / time_series[1:(length(time_series) - 12)]) - 1
    }
  }
  #X <- X[-c(1:13), ]
  X <- X[!is.na(apply(X,1,sum)), ]
  return(X)
}



get.factors <- function(e,S,H,L,q,t){
  F_raw <- matrix(0,t,q)
  for (tt in 1:t){
    normalizer <- exp(-S[tt,]/2)
    Lt <- L*normalizer
    yt <- e[tt,]*normalizer
    
    if (q==1) fac.varinv <-  1/H[tt] else fac.varinv <- diag(q)/H[tt,]
    fac.Sigma <-  solve(crossprod(Lt)+fac.varinv)
    fac.mean <- fac.Sigma%*%crossprod(Lt,yt)
    
    F_temp <- try(fac.mean + t(chol(fac.Sigma)) %*% rnorm(q),silent=TRUE)
    if (is(F_temp,"try-error")) F_temp <- fac.mean + t(chol(fac.Sigma+diag(q)*1e-6)) %*% rnorm(q)
    F_raw[tt,] <- F_temp
  }
  return(F_raw)
}

get.fac.prop <- function(e,S,H,L,q,t){
  F_raw <- matrix(0,t,q)
  f.mom <- list()
  for (tt in 1:t){
    normalizer <- exp(-S[tt,]/2)
    Lt <- L*normalizer
    yt <- e[tt,]*normalizer
    
    if (q==1) fac.varinv <-  1/H[tt] else fac.varinv <- diag(q)/H[tt]
    fac.Sigma <-  solve(crossprod(Lt)+fac.varinv)
    fac.mean <- fac.Sigma%*%crossprod(Lt,yt)
    
    f.mom[[tt]] <- list("mean"=fac.mean, "sigma"=fac.Sigma)
  }
  return(f.mom)
}

# function to draw the factor loadings (basic linear regression)
get.facload <- function(yy,xx,l_sd){
  V_prinv <- diag(NCOL(xx))/l_sd
  V_lambda <- solve(crossprod(xx) + V_prinv)
  lambda_mean <- V_lambda %*% (crossprod(xx,yy))
  
  lambda_draw <- lambda_mean + t(chol(V_lambda)) %*% rnorm(NCOL(xx))
  return(lambda_draw)
}

# factor loadings draw
get.Lambda <- function(eps,fac,S,pr,m,q,id.fac){
  L <- matrix(0,m,q)
  if(id.fac){
    for(jj in 1:m){
      if (jj<=q){
        normalizer <- exp(0.5*S[,jj])
        yy0 <- (eps[,jj]-fac[,jj])/normalizer
        xx0 <- fac[,1:(jj-1),drop=FALSE]/normalizer
        if (jj>1){
          l_sd <- pr[jj,1:(jj-1)]
          lambda0 <- get.facload(yy0,xx0,l_sd=l_sd)
        }else{
          lambda0 <- 1
        }
        
        if (jj>1){
          L[jj,1:(jj-1)] <- lambda0
          L[jj,jj] <- 1
        }else if (jj==1){
          L[jj,jj] <- 1
        }
      }else{
        normalizer <- exp(0.5*S[,jj])
        yy0 <- (eps[,jj])/normalizer
        xx0 <- fac[,,drop=FALSE]/normalizer
        l_sd <- pr[jj,]
        lambda0 <- get.facload(yy0,xx0,l_sd=l_sd)
        L[jj,] <- lambda0
      }
    }
  }else{
    for(jj in 1:m){
      normalizer <- exp(0.5*S[,jj])
      yy0 <- (eps[,jj])/normalizer
      xx0 <- fac[,,drop=FALSE]/normalizer
      l_sd <- pr[jj,]
      lambda0 <- get.facload(yy0,xx0,l_sd=l_sd)
      L[jj,] <- lambda0
    }
  }
  return(L)
}

get.q.dens <- function(y, mean.y, cov.adapt, cov.fixed, alpha){
  (1-alpha)*mvnfast::dmvn(y, mean.y, cov.adapt, isChol=TRUE) + alpha * mvnfast::dmvn(y, mean.y, cov.fixed, isChol = TRUE)
}

mlag <- function(X,lag){
  p <- lag
  X <- as.matrix(X)
  Traw <- nrow(X)
  N <- ncol(X)
  Xlag <- matrix(NA,Traw,p*N)
  for (ii in 1:p){
    Xlag[(p+1):Traw,(N*(ii-1)+1):(N*ii)]=X[(p+1-ii):(Traw-ii),(1:N)]
  }
  return(Xlag)
}
get_companion <- function(Beta_,varndxv){
  nn <- varndxv[[1]]
  nd <- varndxv[[2]]
  nl <- varndxv[[3]]
  
  nkk <- nn*nl+nd
  
  Jm <- matrix(0,nkk,nn)
  Jm[1:nn,1:nn] <- diag(nn)
  
  MM <- rbind(t(Beta_),cbind(diag((nl-1)*nn), matrix(0, (nl-1)*nn, nn)))
  
  return(list(MM=MM,Jm=Jm))
}
compute_coverage_rate <- function(A, C) {
  # Check if the true outcome A is within the 16th and 84th percentiles
  coverage <- ((A >= C[, 1]) & (A <= C[, 2]))*1
  # Compute the mean coverage rate
  return(coverage)
  #mean(coverage)
}
compute_auc_draw <- function(scores, zero_idx) {
  
  labels <- rep(0, length(scores))
  labels[zero_idx] <- 1   # 1 = should be zero
  
  # ranks automatically handle ties
  r <- rank(scores)
  
  n1 <- sum(labels == 1)   # zeros
  n0 <- sum(labels == 0)   # nonzeros
  
  auc <- (sum(r[labels == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  
  return(auc)
}


compute_posterior_auc <- function(posterior_array, sparse_M) {
  
  # flip sign so larger = stronger shrinkage
  scores_mat <- -posterior_array
  
  auc_draws <- apply(scores_mat, 1, compute_auc_draw, 
                     zero_idx = sparse_M)
  
  list(
    mean_auc = mean(auc_draws),
    median_auc = median(auc_draws),
    ci_95 = quantile(auc_draws, c(0.025, 0.975)),
    draws = auc_draws
  )
}
