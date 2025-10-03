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
    // Loop over times
    for (int i = 0; i < t.n_rows; i++) {
      output(i, j) = arma::mean(predictSurv1(t(i), X, betaX.row(j).t(), b(j), ln_p(j)));
    }
  }
  return(output);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(predictIntSurv)]]
arma::vec predictIntSurv(
    const double& t,
    const arma::mat& X,
    const arma::vec& betaX,
    const double& b,
    const double& ln_p,
    const double& varmarg,
    const arma::vec& dnrm,
    const arma::vec& GHx,
    const arma::vec& GHw) {
  Rcpp::checkUserInterrupt();
  unsigned int GHxnr = GHx.n_rows;
  unsigned int Xnr = X.n_rows;
  arma::mat S1(Xnr, GHxnr);
  double tmp = sqrt(varmarg) * sqrt(2);
  for (unsigned int i = 0; i < GHxnr; i++) {
    double this_b = b + GHx(i) * tmp;
    S1.col(i) = predictSurv1(t, X, betaX, this_b, ln_p) * dnrm(i);
  }
  arma::vec Sint = tmp * (S1 * (GHw % exp(pow(GHx, 2.0))));
  return(Sint);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export(predictMeanIntSurv)]]
arma::mat predictMeanIntSurv(
    const arma::vec& t,
    const arma::mat& X,
    const arma::mat& betaX,
    const arma::vec& b,
    const arma::vec& ln_p,
    const arma::vec& varmarg,
    const arma::vec& GHx,
    const arma::vec& GHw
    ) {
  // Check
  Rcpp::checkUserInterrupt();
  // Create matrix for output
  unsigned int tnr = t.n_rows;
  unsigned int bXnr = betaX.n_rows;
  arma::mat output(tnr, bXnr);
  // dnorm
  arma::mat dnrm(varmarg.n_rows, GHx.n_rows);
  for (unsigned int i = 0; i < varmarg.n_rows; i++) {
    dnrm.row(i) = normpdf(GHx * sqrt(varmarg(i)) * sqrt(2), 0, sqrt(varmarg(i))).t();
  }
  // Loop over reps
  for (unsigned int j = 0; j < bXnr; j++) {
    // Loop over times
    for (unsigned int i = 0; i < tnr; i++) {
      output(i, j) = arma::mean(predictIntSurv(t(i), X, betaX.row(j).t(), b(j), ln_p(j), varmarg(j), dnrm.row(j).t(), GHx, GHw));
    }
  }
  return(output);
}
