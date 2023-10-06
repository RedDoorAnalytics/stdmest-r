#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(predictSurv1)]]
arma::vec predictSurv1(
    const double& t,
    const arma::mat& X,
    const arma::vec& betaX,
    const double& b,
    const double& ln_p) {
  // Survival
  arma::colvec S = exp(-exp(X * betaX + b) * pow(t, exp(ln_p)));
  // Return
  return(S);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(predictMeanSurv1)]]
arma::mat predictMeanSurv1(
    const arma::vec& t,
    const arma::mat& X,
    const arma::mat& betaX,
    const arma::vec& b,
    const arma::vec& ln_p) {
  // Check
  Rcpp::checkUserInterrupt();
  // Create matrix for output
  arma::mat output(t.n_rows, betaX.n_rows);
  // Loop over reps
  for (int j = 0; j < betaX.n_rows; j++) {
    // Check
    Rcpp::checkUserInterrupt();
    // Loop over times
    for (int i = 0; i < t.n_rows; i++) {
      output(i, j) = arma::mean(predictSurv1(t(i), X, betaX.row(j).t(), b(j), ln_p(j)));
    }
  }
  return(output);
}
