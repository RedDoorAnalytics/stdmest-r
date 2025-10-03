#' @title `merlin`-Compatible Restricted Cubic Splines
#'
#' @param variable Numeric vector to calculate spline expansion for.
#'
#' @param knots Knots used by `merlin`; for baseline hazard functions (e.g., when
#'     fitting Royston-Parmar mixed-effects models) this can be extracted with
#'     the [read_e()] function.
#'
#' @param deriv Pass values 1, 2, or 3 if you need derivatives of the spline.
#'     Otherwise, defaults to 0.
#'
#' @param rmatrix Matrix used to orthogonalise the splines (if needed). Defaults
#'     to `NULL`, in which case no orthogonalisation is performed.
#'
#' @return A matrix with the spline expansion of `variable`.
#'
#' @export
#'
merlin_rcs = function(variable, knots, deriv = 0, rmatrix = NULL) {

}
# `RM' merlin_rcs(real colvector variable,	///
#                 real rowvector knots,| 		///
#                 real scalar deriv,			///
#                 real matrix rmatrix			///
#         )
# {
# 	real scalar Nobs, Nknots, kmin, kmax, interior, Nparams
# 	real matrix splines, knots2
#
# 	if (args()==2) deriv = 0
#
# 	//======================================================================================================================================//
# 	// Extract knot locations
#
# 		Nobs 	= rows(variable)
# 		Nknots 	= cols(knots)
# 		kmin 	= knots[1,1]
# 		kmax 	= knots[1,Nknots]
#
# 		if (Nknots==2) interior = 0
# 		else interior = Nknots - 2
# 		Nparams = interior + 1
#
# 		splines = J(Nobs,Nparams,.)
#
# 	//======================================================================================================================================//
# 	// Calculate splines
#
# 		if (Nparams>1) {
# 			lambda = J(Nobs,1,(kmax:-knots[,2..Nparams]):/(kmax:-kmin))
# 			knots2 = J(Nobs,1,knots[,2..Nparams])
# 		}
#
# 		if (deriv==0) {
# 			splines[,1] = variable
# 			if (Nparams>1) {
# 				splines[,2..Nparams] = (variable:-knots2):^3 :* (variable:>knots2) :- lambda:*((variable:-kmin):^3):*(variable:>kmin) :- (1:-lambda):*((variable:-kmax):^3):*(variable:>kmax)
# 			}
# 		}
# 		else if (deriv==1) {
# 			splines[,1] = J(Nobs,1,1)
# 			if (Nparams>1) {
# 				splines[,2..Nparams] = 3:*(variable:-knots2):^2 :* (variable:>knots2) :- lambda:*(3:*(variable:-kmin):^2):*(variable:>kmin) :- (1:-lambda):*(3:*(variable:-kmax):^2):*(variable:>kmax)
# 			}
# 		}
# 		else if (deriv==2) {
# 			splines[,1] = J(Nobs,1,0)
# 			if (Nparams>1) {
# 				splines[,2..Nparams] = 6:*(variable:-knots2) :* (variable:>knots2) :- lambda:*(6:*(variable:-kmin)):*(variable:>kmin) :- (1:-lambda):*(6:*(variable:-kmax)):*(variable:>kmax)
# 			}
# 		}
# 		else if (deriv==3) {
# 			splines[,1] = J(Nobs,1,0)
# 			if (Nparams>1) {
# 				splines[,2..Nparams] = 6:*(variable:>knots2) :- lambda:*6:*(variable:>kmin) :- (1:-lambda):*6:*(variable:>kmax)
# 			}
# 		}
#
# 		//orthog
# 		if (args()==4) {
# 			real matrix rmat
# 			rmat = luinv(rmatrix)
# 			if (deriv==0) splines = (splines,J(Nobs,1,1)) * rmat[,1..Nparams]
# 			else splines = splines * rmat[1..Nparams,1..Nparams]
# 		}
# 		return(splines)
# }
#
# `RM' merlin_orthog(`RM' x) 	//from rcsgen.ado SSC
# {
# 	meanx = mean(x)
# 	v = x :- meanx ,J(rows(x),1,1)
# 	q = J(rows(v),0,.)
# 	R = J(cols(v),cols(v),0)
# 	R[cols(v),] = (meanx,1)
# 	for (i=1;i<=cols(x);i++){
#                 r = norm(v[,i])/sqrt(rows(v))
#                 q = q, (v[,i]:/ r)
#                 R[i,i] = r
#                 for (j = i + 1; j<=cols(x); j++){
#                         r = (q[,i]' * v[,j])/rows(v)
# v[,j] = v[,j] - r*q[,i]
# R[i,j] = r
# }
# }
# return(R)
# }
