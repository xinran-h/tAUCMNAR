#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


using namespace Rcpp;
using namespace arma;



//' covariance_cal
//'
//' This function computes three matrices, \eqn{\Sigma_1},\eqn{\Sigma_2}, and \eqn{V}.  \eqn{V} is the asymptotic variance, defined as \eqn{\Sigma_1^{-1} \Sigma_2 \Sigma_1^{-1}}.
//'
//' @param a a vector of regression coefficients from the AUC logistic regression model.
//' @param b a vector of biomarker values.
//' @param c the \code{rho01} vector from the dataframe returned by the \code{data_crossingdc} function. When there are no continuous variables, use the \code{Iil} vector from the dataframe returned by the \code{data_crossingdc} function.
//' @param d the \code{rho02} vector from the dataframe returned by the \code{data_crossingdc} function. When there are no continuous variables, use the \code{1-Iil} vector from the dataframe returned by the \code{data_crossingdc} function.
//' @param e A numeric matrix with k rows and n columns. Each row represents a covariate used in the AUC logistic regression model, including the intercept as the first column. The number of columns is equal to the number of rows in the dataframe returned by the \code{data_crossingdc} function. 
//' @param f the \code{idi - 1} vector from the dataframe returned by the \code{data_crossingdc} function.
//' @param g the \code{idl - 1} vector from the dataframe returned by the \code{data_crossingdc} function.
//' @param h a numeric matrix with q rows and n columns. Each row represents a covariate used in the missing data logistic regression model for cases, including the intercept as the first column. The number of columns is equal to the number of rows in the dataframe returned by the `data_crossingdc` function 
//' @param i a numeric matrix with q rows and n columns. Each row represents a covariate used in the missing data logistic regression model for controls in the risk set, including the intercept as the first column. The number of columns is equal to the number of rows in the dataframe returned by the `data_crossingdc` function 
//' @param j a vector of estimated probability of observing biomarkers among case patients.
//' @param k a vector of estimated probability of observing biomarkers among control patients.
//' @param l a vector of regression coefficiens from the missing data model.
//' @param m a numeric matrix, returned from the \code{pi.gamma.omega.est} function, accessed through \code{pi.gamma.omega.est$composit}.
//' @param n a numeric matrix,, returned from the \code{pi.gamma.omega.est} function, accessed through \code{pi.gamma.omega.est$h.hat}.
//' @return A list containing:
//' \itemize{
//'   \item \code{sigma1}: The matrix \eqn{\Sigma_1}, used in calculating the asymptotic variance.
//'   \item \code{sigma2}: The matrix \eqn{\Sigma_2}, also used in calculating the asymptotic variance.
//'   \item \code{V}: The asymptotic variance for the regression coefficents.
//' }
// [[Rcpp::export]]
Rcpp::List covariance_cal(arma::vec a, arma::vec b, arma::vec c,
                            arma::vec d, 
                            arma::mat e,
                            arma::ivec f, arma::ivec g,
                            arma::mat h, arma::mat i,
                            arma::vec j, arma::vec k,
                            arma::vec l, arma::mat m, arma::mat n) {
  // Define the inputs
  arma::vec beta = a;
  arma::vec M = b;
  arma::vec rho01= c;
  arma::vec rho02 = d;
  arma::mat XS = e;
  arma::ivec dati = f;  
  arma::ivec datj = g;  
  arma::mat XI = h;
  arma::mat XL = i;
  arma::vec piI = j;
  arma::vec piL = k;
  arma::vec phi = l;
  arma::mat composit = m;
  arma::mat h_function = n;


  int n0 = M.size();
  int nij = rho01.size();
  int nb = beta.size();
  int np = phi.size();

  arma::mat sigma1(nb, nb, arma::fill::zeros);
  arma::mat sigma2(nb, np, arma::fill::zeros);
  arma::cube f_store(n0, n0, nb, arma::fill::zeros);
  arma::vec poly(nb, arma::fill::zeros);
  arma::vec poly2(np, arma::fill::zeros);
  arma::vec poly3(np, arma::fill::zeros);
  arma::vec zeta(nb, arma::fill::zeros);
  arma::vec expectation_gij(nb, arma::fill::zeros);
  arma::mat V(nb, nb, arma::fill::zeros);
  arma::mat Vleft(nb, nb, arma::fill::zeros);
  arma::mat Vright(nb, nb, arma::fill::zeros);


  // First loop
  for (int i = 0; i < nij; ++i) {
    // sigma1
    poly = XS.col(i);
    double temp = exp(arma::dot(beta, poly));
    sigma1 += (poly * poly.t()) * (rho02[i] + rho01[i]) * temp / ((1.0 + temp) * (1.0 + temp));
    arma::vec f_storeij = poly / (1.0 + temp) * (rho01[i] - rho02[i] * temp);
    f_store.subcube(dati[i], datj[i], 0, dati[i], datj[i], nb - 1) = f_storeij;

    // sigma2
    poly2 = XI.col(i);
    double temp2 = piI[i];
    arma::vec delta_pii = poly2 * (1 - temp2);
  
    arma::vec poly3 = XL.col(i);
    double temp3 = piL[i];
    arma::vec delta_pil = poly3 * (1 - temp3);
  
    arma::vec delta_pi = delta_pii + delta_pil;
    
    sigma2 = sigma2 + f_storeij * delta_pi.t(); 
  }

  arma::mat sigma1_inv = arma::inv(sigma1);
  // V
  for (int i = 0; i < n0; ++i) {
    zeta.fill(0.0);
    for (int j = 0; j < n0; ++j) {
      if (j != i) {
        arma::vec gij = f_store.subcube(i, j, 0, i, j, nb - 1) + f_store.subcube(j, i, 0, j, i, nb - 1);
        zeta += gij;
      }
    }
    arma::vec h = h_function.row(i).t();
    arma::vec out = n0 * sigma1_inv * zeta + sigma1_inv * sigma2 * composit * h;
    V += out * out.t();
  }
  
return Rcpp::List::create(Named("sigma1") = sigma1, Named("sigma2") = sigma2, Named("V") = V/n0);
}

