####  All functions ####

#' eta_function
#' 
#' This is an internal function used to generate \eqn{$\eta$} vector in the GEE.
#' @param x  a numeric matrix containing the design matrix in the missing data logistic regression model.
#' @param z_cont a vector of column names from x representing continuous covariates that are non-instrumental variables. If there is no continuous covariate, set it as NULL.
#' @param z_dis a vector of column names from x representing categorical covariates that are non-instrumental variables. If there is no categorical covariate, set it as NULL.
#' @param u_type a character string indicating the type of instrumental variable, set to "continuous" or "dicrete".
#' @param u a numeric vector representing the instrumental variable. 
#' @param V a numeric vector of missing data indicators, V = 1 if observed and V = 0 if missing.
#' @param futime a numeric vector of observed event time.
#' @param status a numeric vector of censoring indicator, status = 1 if event occurs and status = 0 if censored.
#' @returns A list, where each component is a numeric vector representing one component of the \eqn{$\eta$} vector.
#' @examples
#' \dontrun{
#' eta_function( x, c("age","bmi"), c("sex","race") ,"discrete",u, V, futime, status)
#'
#' } 
#' @noRd
eta_function <- function(x, z_cont, z_dis, u_type, u, V, futime, status){
  eta_list <- list()
  if (!is.null(z_cont)) {
    for (var in z_cont) {
      eta_list[[paste0("g_", var)]] <- x[, var]
    }
  }
  
  if (!is.null(z_dis)) {
    for (var in z_dis) {
      levels <- unique(x[, var])
      for (level in levels[-1]) {  # Exclude the reference category
        eta_list[[paste0("g_", var, "_", level)]] <- as.integer(x[, var] == level) 
      }
    }
  }
  # Continuous u variables
  if (u_type == "continuous") {
    u2 = u**2
    eta_list[["g_u"]] <- u 
    eta_list[["g_u2"]] <- u2 
  }else if (u_type == "discrete") {
    levels <- unique(u)
    for (level in levels) {
      eta_list[[paste0("g_u_", level)]] <- as.integer(u == level) 
    }
  } else {
    stop("u_type must be either 'continuous' or 'discrete'")
  }
  
  
  eta_list[["g_futime"]]= futime
  eta_list[["g_status"]] = status
  return(eta_list)
}



#' g_function
#' 
#' This is an internal function called by the est.pi.func function and weighting_matrix function.
#' @param phi   parameters for the missing data logistic regression model.
#' @param x  a numeric matrix containing the design matrix in the missing data logistic regression model.
#' @param z_cont a vector of column names from x representing continuous covariates that are non-instrumental variables. If there is no continuous covariate, set it as NULL.
#' @param z_dis a vector of column names from x representing categorical covariates that are non-instrumental variables. If there is no categorical covariate, set it as NULL.
#' @param u_type a character string indicating the type of instrumental variable, set to "continuous" or "dicrete".
#' @param u a numeric vector representing the instrumental variable. 
#' @param V a numeric vector of missing data indicators, V = 1 if observed and V = 0 if missing.
#' @param futime a numeric vector of observed event time.
#' @param status a numeric vector of censoring indicator, status = 1 if event occurs and status = 0 if censored.
#' @returns A numeric matrix with n rows and q columns, where n is the number of observations and q is the number of estimating equations in the GEE for the missing data logistic regression model.
#' @examples
#' \dontrun{
#' g_function(phi, x, c("age","bmi"), c("sex","race") ,"discrete",u, V, futime, status)
#'
#' } 
#' @noRd
g_function <- function(phi, x, z_cont, z_dis, u_type, u, V, futime, status) {
  pi <- as.vector(1 / (1 + exp(-x %*% phi)))
  temp = ifelse(V==1, V/pi - 1, -1) 
  
  temp_list = eta_function(x, z_cont, z_dis, u_type, u, V, futime, status)
  g_list = lapply(temp_list, function(g) g*temp)
  g_matrix = do.call(cbind, g_list)
  return(g_matrix) 
}

#' g_function_optim
#' 
#' This is an internal function that called by the objective function and DG_optim function.
#' @param phi   parameters for the missing data logistic regression model.
#' @param x  a numeric matrix containing the design matrix in the missing data logistic regression model.
#' @param z_cont a vector of column names from x representing continuous covariates that are non-instrumental variables. If there is no continuous covariate, set it as NULL.
#' @param z_dis a vector of column names from x representing categorical covariates that are non-instrumental variables. If there is no categorical covariate, set it as NULL.
#' @param u_type a character string indicating the type of instrumental variable, set to "continuous" or "dicrete".
#' @param u a numeric vector representing the instrumental variable. 
#' @param V a numeric vector of missing data indicators, V = 1 if observed and V = 0 if missing.
#' @param futime a numeric vector of observed event time.
#' @param status a numeric vector of censoring indicator, status = 1 if event occurs and status = 0 if censored.
#' @returns A numeric matrix with 1 row and q columns, where q is the number of estimating equations in the GEE for the missing data logistic regression model.
#' @examples
#' \dontrun{
#' g_function_optim(phi, x, c("age","bmi"), c("sex","race") ,"discrete",u, V, futime, status)
#'
#' } 
#' @noRd
g_function_optim <- function(phi, x, z_cont, z_dis, u_type, u, V, futime, status) {
  pi <- as.vector(1 / (1 + exp(-x %*% phi)))
  temp = ifelse(V==1, V/pi - 1, -1) 
  temp_list = eta_function(x, z_cont, z_dis, u_type, u, V, futime, status)
  g_list = lapply(temp_list, function(g) mean(g*temp))
  g_matrix = do.call(c, g_list)
  return(g_matrix)
}

#' objective
#' 
#' This is an internal function used to obtain the objective function to estimate parameters in the missing data logistic regression model.
#' @param phi   parameters for the missing data logistic regression model.
#' @param x  a numeric matrix containing the design matrix in the missing data logistic regression model.
#' @param z_cont a vector of column names from x representing continuous covariates that are non-instrumental variables. If there is no continuous covariate, set it as NULL.
#' @param z_dis a vector of column names from x representing categorical covariates that are non-instrumental variables. If there is no categorical covariate, set it as NULL.
#' @param u_type a character string indicating the type of instrumental variable, set to "continuous" or "dicrete".
#' @param u a numeric vector representing the instrumental variable. 
#' @param V a numeric vector of missing data indicators, V = 1 if observed and V = 0 if missing.
#' @param futime a numeric vector of observed event time.
#' @param status a numeric vector of censoring indicator, status = 1 if event occurs and status = 0 if censored.
#' @param W a numeric matrix of weight matrix.
#' @returns objective function to be used in optim() to estimate the parameters of the missing data logistic regression model.
#' @examples
#' \dontrun{
#' objective(phi, x, c("age","bmi"), c("sex","race") ,"discrete",u, V, futime, status,  diag(1, 8,8))
#'
#' } 
#' @noRd
objective <- function(phi, x, z_cont, z_dis, u_type, u, V, futime, status, W){
  moment = g_function_optim(phi, x, z_cont, z_dis, u_type, u, V, futime, status)
  return(as.numeric(t(moment)%*%W%*%moment ))
}

#' derivative_pi
#' 
#' This is an internal function used within the function Dg_gmm. 
#' @param phi   parameters for the missing data logistic regression model.
#' @param x  a numeric matrix containing the design matrix in the missing data logistic regression model.
#' @param V a numeric vector of missing data indicators, V = 1 if observed and V = 0 if missing.
#' @param eta a variable from the \eqn{$\eta$} vector in the GEE.
#' @returns  A numeric vector of size \(1 \times k\), where \(k\) is the number of parameters in `phi`. 
#' 
#' @noRd
derivative_pi <- function(phi, x, V, eta) {
  
  
  n <- nrow(x)
  pi <- numeric(n)
  derivative_pi <- matrix(nrow = n, ncol = length(phi))
  
  xtheta <- x %*% phi
  pi <- as.vector(1 / (1 + exp(- xtheta)))
  derivative.pi <- x * (pi * (1 - pi))
  
  
  temp = -(V * eta / pi^2) * derivative.pi
  temp[is.na(temp)] <- 0
  
  derivative_Y <- colMeans(temp)
  
  return(derivative_Y) 
}

#' Dg_gmm
#' 
#' This is an internal function called by Dg_optim.
#' @param phi   parameters for the missing data logistic regression model.
#' @param x  a numeric matrix containing the design matrix in the missing data logistic regression model.
#' @param z_cont a vector of column names from x representing continuous covariates that are non-instrumental variables. If there is no continuous covariate, set it as NULL.
#' @param z_dis a vector of column names from x representing categorical covariates that are non-instrumental variables. If there is no categorical covariate, set it as NULL.
#' @param u_type a character string indicating the type of instrumental variable, set to "continuous" or "dicrete".
#' @param u a numeric vector representing the instrumental variable. 
#' @param V a numeric vector of missing data indicators, V = 1 if observed and V = 0 if missing.
#' @param futime a numeric vector of observed event time.
#' @param status a numeric vector of censoring indicator, status = 1 if event occurs and status = 0 if censored.
#' @returns  derivative of the moment functions with respect to the parameters of the missing data logistic regression model.

#' @noRd
Dg_gmm <- function(phi, x, z_cont, z_dis, u_type, u, V, futime, status){
  
  temp_list = eta_function(x, z_cont, z_dis, u_type, u, V, futime, status)
  dg_list = lapply(temp_list, function(g) derivative_pi(phi, x, V, g))
  dg_matrix = do.call(rbind, dg_list)
  return(dg_matrix)  # q * k
}

#' DG_optim
#' 
#' This is an internal function used to obtain the derivative of the objective function to estimate parameters in the missing data logistic regression model.
#' @param phi   parameters for the missing data logistic regression model.
#' @param x  a numeric matrix containing the design matrix in the missing data logistic regression model.
#' @param z_cont a vector of column names from x representing continuous covariates that are non-instrumental variables. If there is no continuous covariate, set it as NULL.
#' @param z_dis a vector of column names from x representing categorical covariates that are non-instrumental variables. If there is no categorical covariate, set it as NULL.
#' @param u_type a character string indicating the type of instrumental variable, set to "continuous" or "dicrete".
#' @param u a numeric vector representing the instrumental variable. 
#' @param V a numeric vector of missing data indicators, V = 1 if observed and V = 0 if missing.
#' @param futime a numeric vector of observed event time.
#' @param status a numeric vector of censoring indicator, status = 1 if event occurs and status = 0 if censored.
#' @param W a numeric matrix of weight matrix
#' @returns objective function to be used in optim() to estimate the parameters of the missing data logistic regression model.
#' @examples
#' \dontrun{
#' DG_optim(phi, x, c("age","bmi"), c("sex","race") ,"discrete",u, V, futime, status,  diag(1, 8,8))
#'
#' } 
#' @noRd
DG_optim <- function(phi,x,z_cont, z_dis, u_type, u, V, futime, status, W){
  moment = g_function_optim(phi, x, z_cont, z_dis, u_type, u, V, futime, status) 
  g_derivative = Dg_gmm(phi, x, z_cont, z_dis, u_type, u, V, futime, status)
  return(as.vector(2*t(g_derivative)%*%W%*%moment))
}

#' weighting_matrix
#' 
#' This is an internal function used to obtain the weighting matrix for the second step in 2-step generalized method of moments.
#' @param phi   parameters for the missing data logistic regression model.
#' @param x  a numeric matrix containing the design matrix in the missing data logistic regression model.
#' @param z_cont a vector of column names from x representing continuous covariates that are non-instrumental variables. If there is no continuous covariate, set it as NULL.
#' @param z_dis a vector of column names from x representing categorical covariates that are non-instrumental variables. If there is no categorical covariate, set it as NULL.
#' @param u_type a character string indicating the type of instrumental variable, set to "continuous" or "dicrete".
#' @param u a numeric vector representing the instrumental variable. 
#' @param V a numeric vector of missing data indicators, V = 1 if observed and V = 0 if missing.
#' @param futime a numeric vector of observed event time.
#' @param status a numeric vector of censoring indicator, status = 1 if event occurs and status = 0 if censored.
#' @returns  weighting matrix for gmm step 2
#' @noRd
weighting_matrix <- function(phi, x, z_cont, z_dis, u_type, u, V, futime, status){
  moment = g_function(phi, x, z_cont, z_dis, u_type, u, V, futime, status)
  W <- t(moment) %*% moment
  return(solve(W)*nrow(x))
}

#' VCOV 
#' 
#' This is to estimate the variance-covariance matrix for the parameter estimates in the missing data logistic regression model
#' @param phi   parameters for the missing data logistic regression model.
#' @param x  a numeric matrix, including intercept and all covariates in the missing data logistic regression model.
#' @param z_cont a vector of column names from x representing continuous covariates that are non-instrumental variables. If there is no continuous covariate, set it as NULL.
#' @param z_dis a vector of column names from x representing categorical covariates that are non-instrumental variables. If there is no categorical covariate, set it as NULL.
#' @param u_type a character string indicating the type of instrumental variable, set to "continuous" or "dicrete".
#' @param u a numeric vector representing the instrumental variable. 
#' @param V a numeric vector of missing data indicators, V = 1 if observed and V = 0 if missing.
#' @param futime a numeric vector of observed event time.
#' @param status a numeric vector of censoring indicator, status = 1 if event occurs and status = 0 if censored.
#' @returns  variance-covariance matrix for the parameter estimates in the missing data logistic regression model
#' @examples
#' \dontrun{
#' VCOV(phi, x, c("age","bmi"), c("sex","race") ,"discrete",u, V, futime, status)
#'
#' } 
#' @noRd
VCOV <- function(phi, x, z_cont, z_dis, u_type, u, V, futime, status){
  
  W <- weighting_matrix(phi, x, z_cont, z_dis, u_type, u, V, futime, status)
  g_derivative = Dg_gmm(phi, x, z_cont, z_dis, u_type, u, V, futime, status)
  return(solve(t(g_derivative)%*%W%*%g_derivative)/nrow(x))
}

#' est.pi.func
#' 
#' This is to perform two-step GMM estimation for the missing data logistic regression model using optim().
#' @param x  a numeric matrix containing the design matrix in the missing data logistic regression model.
#' @param z_cont a vector of column names from x representing continuous covariates that are non-instrumental variables. If there is no continuous covariate, set it as NULL.
#' @param z_dis a vector of column names from x representing categorical covariates that are non-instrumental variables. If there is no categorical covariate, set it as NULL.
#' @param u_type a character string indicating the type of instrumental variable, set to "continuous" or "dicrete".
#' @param u a numeric vector representing the instrumental variable. 
#' @param V a numeric vector of missing data indicators, V = 1 if observed and V = 0 if missing.
#' @param futime a numeric vector of observed event time.
#' @param status a numeric vector of censoring indicator, status = 1 if event occurs and status = 0 if censored.
#' @param q a numeric value indicating the number of estimating equations.
#' @param initial a numeric vector of initial values for the parameters in the missing data logistic regression model.
#' @param method a character string indicating the optimization method, set to  one of the "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent".
#' @returns  a list, containing the following items:
#' \itemize{
#'   \item \code{phi.hat}: estimated regression coefficients from the missing data model.
#'   \item \code{convergence}: convergence code from the second step optimization.
#'   \item \code{convergence.message}: convergence message from the second step optimization.
#'   \item \code{objective}: value of the objective function from the second step optimization.
#'   \item \code{vcov.hat}: variance-covariance matrix for the parameter estimates in the missing data logistic regression model. 
#'   \item \code{pi.hat}: estimated probability of observing the biomarker
#'   \item \code{h.hat}: a numeric matrix with n rows and q columns, each row is \eqn{h(\hat{\phi},y_i)} 
#'   \item \code{composit}:this is \eqn{{(\Gamma^T{\Omega}^{-1}\Gamma)}^{-1}\Gamma^T{\Omega}^{-1}}.

#' }
#' @examples
#' \dontrun{
#' est.pi.func(x, c("age","bmi"), c("sex","race") ,"discrete",u, V, futime, status,8, c(0,0,0,0,0,0,0), "L-BFGS-B")
#'
#' } 
#' @export
est.pi.func <- function(x, z_cont, z_dis, u_type, u, V, futime, status, q,
                        initial, method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent")
                        
){
  # step1
  step1.result <- tryCatch({
    optim(initial, objective, DG_optim, x, z_cont, z_dis, u_type, u, V, futime, status,
          diag(1, q, q), method= method,
          control = list(maxit = 10000))
  }, error = function(e) {
    message(paste("step 1Error:", e$message))
    return(NULL)
  })
  
  
  # step2
  step2.result <- tryCatch({
    optim(step1.result$par, objective, DG_optim, x, z_cont, z_dis, u_type, u, V, futime, status,
          weighting_matrix(step1.result$par, x, z_cont, z_dis, u_type, u, V, futime, status),
          method= method,
          control = list(maxit = 10000))
  }, error = function(e) {
    message(paste("step 2 Error occurred:", e$message))
    return(NULL)
  })
  
  vcov =tryCatch({
    VCOV(step2.result$par, x, z_cont, z_dis, u_type, u, V, futime, status)}, error = function(e) {
      message(paste("Variance estimating error:", e$message))
      return(NULL)
    })
  
  
  g.function = g_function(step2.result$par,x, z_cont, z_dis, u_type, u, V, futime, status)
  
  xtheta <- x %*% step2.result$par
  pi <- as.vector(exp(xtheta) / (1 + exp(xtheta)))
  
  omega <- t(g.function) %*% g.function
  
  
  gamma = Dg_gmm(step2.result$par,x,z_cont, z_dis, u_type, u, V, futime, status)
  omega_inv = solve(omega)*nrow(x)
  
  composit = solve(t(gamma) %*% omega_inv %*% gamma)%*%t(gamma)%*%omega_inv
  
  outcome = list(phi.hat =  step2.result$par, convergence = step2.result$convergence,
                 convergence.message = step2.result$message,
                 objective = step2.result$value, vcov.hat = vcov,
                 pi.hat = pi,h.hat = g.function, composit = composit)
  return(outcome)
}



#' data_crossingdc
#' 
#' This is to reshape data for each cause of failure, the result reshaped data is ready for parameter estimation.
#' @import dplyr
#' @import rlang
#' @import tidyr
#' @param V  a numeric vector of missing data indicators, V = 1 if observed and V = 0 if missing
#' @param vt  a numeric vector of observed time
#' @param vc a numeric vector of censoring indicator
#' @param vm a numeric vector of biomarker values
#' @param vxs a covariate matrix used for the AUC logistic regression, do not include intercept or observed time
#' @param pi estimated probability of observing biomarker values
#' @param cont_vars a vector of column names from vxs representing continuous covariates
#' @param discrete_vars a vector of column names from vxs representing discrete covariates
#' @param c0 A constant value used for the bandwidth h, h = c0*n^(-1/3)
#' @returns A dataframe where each row corresponds to a case-control pair at each observed event time.
#' The dataframe has the following components:
#' \itemize{
#'   \item \code{i}: numeric index for cases.
#'   \item \code{idi}: patient id for cases.
#'   \item \code{idl}: patient id for controls.
#'   \item \code{yi}: observed event time for cases.
#'   \item \code{yl}: observed event time for controls in the risk set of i.
#'   \item \code{Iil}: indicator for comparison of baseline biomarker values (1 if case biomarker > control biomarker, 0 otherwise).
#'   \item \code{pii}: estimated probability of observing the biomarker for cases.
#'   \item \code{pil}: estimated probability of observing the biomarker for controls in the risk set of i.
#'   \item \code{event_indicator_i}: censoring indicator for cases, this is always 1.
#'   \item \code{event_indicator_l}: censoring indicator for controls in the risk set of i. 1 if event occurred, 0 if censored.
#'   \item \code{mi}: biomarker values for cases.
#'   \item \code{ml}: biomarker values for controls in the risk set of i.
#'   \item \code{Covariates with "i" suffix}: covariate columns, with names based on the original covariate name and an "i" suffix for case values.
#'   \item \code{rho01} and \code{rho02}: numeric vectors representing concordant and discordant event values, respectively, when continuous covariates are present.
#'   \item \code{rhoweight}: Kernel weight vector, calculated with continuous covariates.
#' }
#' @examples
#' \dontrun{
#' Y0<-x[,"time_year"]   
#' C0<-x[,"status"]
#' M0<-x[,"TAU_bl"]
#' VXS<-x[,c("age","sex", "race")];
#' V0 <- dat$V_tau
#' pi.hat <- est.pi.func(x, c("age","bmi"), c("sex","race") ,"discrete",u, V, futime, status,8, c(0,0,0,0,0,0,0), "L-BFGS-B")$pi.hat
#' data_crossingdc(V, Y0,C0,M0,VXS, pi.hat, , c("age"),c("sex", "race"), 1.5)
#'
#' } 
#' @export

data_crossingdc = function(V, vt,vc,vm,vxs,pi, cont_vars,discrete_vars, c0)
{
  
  nx<-ncol(vxs) 
  ni<-length(vc) 
  dat<-cbind(id = 1:ni,  vt,vm,vc, pi,vxs, V) # id: to match zeta and h in the V
  if (!is.data.frame(dat)) {
    dat <- as.data.frame(dat)
  }
  
  # remove data with missing vm
  dat <- dat%>%filter(V == 1)%>%select(-V)
  
  # new: rank dat by ascending order of vt
  dat<-dat[order(dat$vt),]
  ni<-sum(V) 
  dat <- cbind(l=1:ni, dat)
  dati<-dat%>%filter(vc==1)%>%rename(i=l)
  
  names(dati)[1:6]<-c("i", "idi", "yi","mi","event_indicator_i", "pii")
  names(dati)[7:(nx+6)]<-paste0(colnames(vxs),"i") 
  
  names(dat)[1:6]<-c("l","idl","yl","ml","event_indicator_l","pil") 
  names(dat)[7:(nx+6)]<-paste0(colnames(vxs),"l") 
  
  
  if (!is.null(cont_vars)) {
    h <- c0 * ni^(-1 / 3)
    
    squared_diffs <- vector("list", length(cont_vars))
    names(squared_diffs) <- cont_vars
    for (var in cont_vars) {
      case_var <- paste0(var, "i")
      control_var <- paste0(var, "l")
      squared_diffs[[var]] <- paste0("(", case_var, " - ", control_var, ")^2")
    }
    sum_of_squares_expr <- paste(squared_diffs, collapse = " + ")
    final_expr <- paste0("sqrt(", sum_of_squares_expr, ")")
  }
  
  if (!is.null(discrete_vars)) {
    conditions <- lapply(discrete_vars, function(var) {
      sym_case <- rlang::sym(paste0(var, "i"))
      sym_control <- rlang::sym(paste0(var, "l"))
      rlang::expr(!!sym_case == !!sym_control)
    })
  }
  
  
  
  tableij <- tidyr::crossing(dati, dat, .name_repair = "universal") %>%
    dplyr::filter(yi < yl) %>%
    dplyr::mutate(Iil = as.numeric(mi > ml))
  if (!is.null(discrete_vars)) {
    tableij <- tableij %>% dplyr::filter(!!!conditions)
  }
  tableij <- tableij %>% dplyr::group_by(i) %>%
    dplyr::arrange(i, ml)
  
  if (!is.null(cont_vars)) {
    tableij <- tableij %>%
      dplyr::mutate(
        dif = eval(parse(text = final_expr)),
        kh = 0.75 * (1 - (dif / h)^2) / h * as.numeric(abs(dif / h) < 1),
        rho01 = Iil * kh,
        rho02 = (1 - Iil) * kh,
        rhoweight = rho01 + rho02
      )
  }
  
  select_columns <- c("i", "idi", "idl", "yi","yl", "Iil", "pii", "pil","event_indicator_i","event_indicator_l", "mi","ml"
  )
  if (!is.null(cont_vars)) {
    select_columns <- c(select_columns, paste0(cont_vars, "i"),paste0(cont_vars, "l"),"rho01", "rho02", "rhoweight")
  }
  if (!is.null(discrete_vars)) {
    select_columns <- c(select_columns, paste0(discrete_vars, "i"), paste0(discrete_vars, "l"))
  }
  
  tableij <- tableij %>%
    dplyr::select(dplyr::any_of(select_columns))
  
  
  return(data.frame(tableij))
  
}




#' auc_pred
#' 
#' This is to estimate the time-dependent AUC with its 95% confidence interval.

#' @param beta.hat   estimated regression coefficients from the AUC regression model.
#' @param V the estimated asymptotic variance from the AUC regression model.
#' @param tt a vector of time.
#' @param coeffs a vector of covariate values, do not include intercept or observed time, ordered to correspond to the elements in \code{beta.hat} for each covariate
#' @param nf a numeric value of 3 or 7 for the number of polynomials. When nf = 7, the vector of polynomials = c(t^{-2},t^{-1},t^{-.5},log(t),t^{0.5},t,t^2). When nf = 3, the vector of polynomials = c(t^{0.5},t, t^2)
#' @returns a matrix of 95% confidence intervals for the time-dependent AUC, where the first row represents point estimate, the second row is the lower level, and the third row is the upper level.
#' @examples
#' \dontrun{
#' 
#' auc_pred(beta.hat,V.hat, seq(0.05,1,0.01),c((mean_age-min_age)/(max_age - min_age),1,1),3)
#'
#' } 
#' @export
#' 

auc_pred <- function(beta.hat,V,tt,coeffs, 
                     nf=3)
{
  se.store.1<- matrix(0,ncol=length(tt),nrow=3)
  for (i in 1:length(tt)){
    t <- tt[i]
    if (nf == 3) {
      polytt <- c(1,t^(.5),t,t^2,coeffs)
    }else if (nf == 7){
      polytt <- c(1,t^(-2),t^(-1),t^(-.5),log(t),t^(.5),t,t^2,coeffs)
    }
    temp <- c(beta.hat %*% polytt)
    lower <- temp- 1.96*sqrt(max(0,c(polytt %*% V %*% polytt)) )
    upper <- temp+ 1.96*sqrt(max(0,c(polytt %*% V %*% polytt)) )
    
    se.store.1[1,i] <- 1/(1+exp(-temp))
    se.store.1[2,i] <- 1/(1+exp(-lower))
    se.store.1[3,i] <-  1/(1+exp(-upper))
  }
  return(se.store.1)
}




#' Sum a Numeric Vector Using a C Function
#'
#' This wrapper provides an R interface to the `covariance_cal` C++ function.
#' @name covariance_cal_wrapper
#' @param a a numeric vector of regression coefficients from the AUC logistic regression model.
#' @param b a numeric vector of biomarker values.
#' @param c the \code{rho01} vector from the dataframe returned by the \code{data_crossingdc} function. When there are no continuous variables, use the \code{Iil} vector from the dataframe returned by the \code{data_crossingdc} function.
#' @param d the \code{rho02} vector from the dataframe returned by the \code{data_crossingdc} function. When there are no continuous variables, use the \code{1-Iil} vector from the dataframe returned by the \code{data_crossingdc} function.
#' @param e a numeric matrix with k rows and n columns. Each row represents a covariate used in the AUC logistic regression model, including the intercept as the first column. The number of columns is equal to the number of rows in the dataframe returned by the \code{data_crossingdc} function. 
#' @param f the \code{idi - 1} vector from the dataframe returned by the \code{data_crossingdc} function.
#' @param g the \code{idl - 1} vector from the dataframe returned by the \code{data_crossingdc} function.
#' @param h a numeric matrix with q rows and n columns. Each row represents a covariate used in the missing data logistic regression model for cases, including the intercept as the first column. The number of columns is equal to the number of rows in the dataframe returned by the `data_crossingdc` function 
#' @param i a numeric matrix with q rows and n columns. Each row represents a covariate used in the missing data logistic regression model for controls in the risk set, including the intercept as the first column. The number of columns is equal to the number of rows in the dataframe returned by the `data_crossingdc` function 
#' @param k a vector of estimated probability of observing biomarkers among case patients.
#' @param k a vector of estimated probability of observing biomarkers among control patients in the risk set of i.
#' @param l a vector of regression coefficiens from the missing data model.
#' @param m a numeric matrix, returned from the \code{est.pi.func} function, accessed through \code{composit}.
#' @param n a numeric matrix,, returned from the \code{est.pi.func} function, accessed through \code{h.hat}.
#' @return A list containing:
#'\itemize{
#'  \item \code{sigma1}: The matrix \eqn{\Sigma_1}, used in calculating the asymptotic variance.
#'  \item \code{sigma2}: The matrix \eqn{\Sigma_2}, also used in calculating the asymptotic variance.
#'  \item \code{V}: The asymptotic variance for the regression coefficients.
#' }
#' @export
covariance_cal_wrapper <- function(a, b, c, d, e, f, g,h,i,j,k,l,m,n) {
  # Call the C function
  covariance_cal(a, b, c, d, e, f, g,h,i,j,k,l,m,n)
}