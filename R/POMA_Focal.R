# Main Function (ParetoR) - Pareto-Optimization via Normal Boundary Intersection
# Developer: Q. Chelsea Song
# Contact: qianqisong@gmail.com


#' ParetoR_MultiAI
#'
#' Command function to run Pareto-Optimal algorithm.
#' @param prop Proportion of minority applicants in full applicant pool
#' @param sr Selection ratio
#' @param d Subgroup difference
#' @param R Correlation matrix
#' @param graph If TRUE, plots will be generated for Pareto-optimal curve and predictor weights
#' @return out Pareto-Optimal solution with criterion values and predictor weights
#' @examples
#' \dontrun{
#' # Specify inputs
#' # (1) Proportion of minority applicants (prop) = (# of minority applicants)/(total # of applicants)
#' prop <- 1/4
#' # (2) Selection ratio (sr) = (# of selected applicants)/(total # of applicants)
#' sr <- 0.10
#' # (3) Subgroup differences (d): standardized mean differences between minority and majority subgroups (i.e., majority - minority), on each predictor (in applicant pool)
#' d <- c(1.00, 0.23, 0.09, 0.33)
#' # (4) Correlation matrix (R) = criterion & predictor inter-correlation matrix (in applicant pool)
#' # Format: Predictor_1, ..., Predictor_n, Criterion
#' R <- matrix(c(1, .24, .00, .19, .30,
#'               .24, 1, .12, .16, .30,
#'               .00, .12, 1, .51, .18,
#'               .19, .16, .51, 1, .28,
#'               .30, .30, .18, .28, 1),
#'             (length(d)+1),(length(d)+1))
#' # Fit Pareto-optimal model
#' out = ParetoR(prop, sr, d, R)
#' }
#' @export
POMA_Focal = function(sr, prop, Rx, Rxy, subd, Spac = Spac, graph = TRUE, display_solution = TRUE){

  sr_ParetoR <<- sr
  Rx_ParetoR <<- Rx
  Rxy_ParetoR <<- Rxy
  prop_ParetoR <<- prop
  subd_ParetoR <<- subd

  # Tolerance Level for Algorithm
  TolCon 	= 1.0000e-6 # tolerance of constraint
  TolF 	= 1.0000e-6 # tolerance of objective function
  TolX 	= 1.0000e-7 # tolerance of predictor variable

  # Do not change
  Fnum 	= 2
  Xnum = dim(Rx)[1]
  VLB = c(-(Xnum+1),rep(0,Xnum))
  VUB = c((Xnum+1),rep(1,Xnum))
  X0 = c(1,rep(1/Xnum,Xnum))

  ###### Find Pareto-Optimal Solution ######

  out = NBI(X0, Spac, Fnum, VLB, VUB, TolX, TolF, TolCon, graph=graph, display_solution=display_solution)
  return(out)

}
