# Pareto-Optimization via Normal Boundary Intersection
# Developer: Q. Chelsea Song
# Contact: qianqisong@gmail.com

####################################### NBI Main Function ####################################

#' NBI Main Function
#'
#' Main function for obtaining pareto-optimal solution via normal boundary intersection.
#' @param X0 Initial input for preditor weight vector
#' @param Spac Number of Pareto spaces (i.e., number of Pareto points minus one)
#' @param Fnum Number of criterions
#' @param VLB Lower boundary for weight vector estimation
#' @param VUB Upper boundary for weight vector estimation
#' @param TolX Tolerance index for estimating weight vector, default is 1e-4
#' @param TolF Tolerance index for estimating criterion, default is 1e-4
#' @param TolCon Tolerance index for constraint conditions, default is 1e-7
#' @param graph If TRUE, plots will be generated for Pareto-optimal curve and predictor weights
#' @param display_solution If TRUE, Pareto-optimal solution will be displayed
#' @import nloptr
#' @return Pareto-Optimal solutions
NBI = function(X0,Spac,Fnum,VLB=vector(),VUB=vector(),TolX=1e-4,TolF=1e-4,TolCon=1e-7,graph=TRUE,display_solution=TRUE){

cat('\n Estimating Pareto-Optimal Solution ... \n')

#------------------------------Initialize Options-------------------------------#

  X0 = assert_col_vec(X0)
  VLB = assert_col_vec(VLB)
  VUB = assert_col_vec(VUB)

  # Number of variables
  nvars = length(X0)

  suppressMessages(nloptr::nl.opts(optlist = list(
             maxeval = (nvars+1)*1000
             ,xtol_rel = TolX
             ,ftol_rel = TolF
             )))

  #Initialize PHI

  PHI = matrix(0,Fnum,Fnum)

  #------------------------------Shadow Point-------------------------------#

  # cat('\n ----Step 1: find shadow minimum---- \n')

  ShadowF = matrix(0,Fnum)
  ShadowX = matrix(0,nvars,Fnum)
  xstart  = X0
  out = WeightsFun(Fnum,Spac)
  Weight = out$Weights
  Near = out$Formers
  rm(out)
  Weight = Weight/Spac

  for(i in 1:Fnum){

    temp = c(1,dim(Weight)[2])
    j = temp[i]
    g_Weight <<- Weight[,j]
    fmin = 9999
    ntr = nvars-1
    fminv = matrix(0,ntr,1)
    fminx = matrix(0,nvars,ntr)

    for(k in 1:ntr){

      xstart = runif(length(X0))

      out = suppressMessages(nloptr::slsqp(x0 = X0, fn = myLinCom
                 ,lower = VLB, upper = VUB
                 ,hin = myCon_ineq
                 ,heq = myCon_eq
                  ))
      x = out$par
      f = out$value
      rm(out)

      fminv[k] = f
      fminx[,k] = x

      if(f <= fmin){

        fmin = f
        reps = k

      }

    }

    x = fminx[,reps]
    som = 0

    for(k in 2:nvars){
      som = som + x[k]
    }

    for(k in 2:nvars){
      x[k] = x[k]/som
    }
    # to make sum of x = 1

    ShadowX[,i] = x
    ShadowX = round(ShadowX,4)

    tempF = -myFM(x)$f
    ShadowF[i] = round(tempF[i],4)

  }

  #------------------------------Matrix PHI-------------------------------#

  # cat('\n ----Step 2: find PHI---- \n')

  for(i in 1:Fnum){

    PHI[,i] = myFM(ShadowX[,i])$f + ShadowF
    PHI[i,i] = 0

  }

  # Check to make sure that QPP is n-1 dimensional
  if(rcond(PHI) < 1e-8){stop(' Phi matrix singular, aborting.')}

  #------------------------------Quasi-Normal Direction-------------------------------#

  # cat('\n ----Step 3: find Quasi-Normal---- \n')

  g_Normal <<- -PHI%*%matrix(1,Fnum,1)

  #------------------------------weights-------------------------------#

  # cat('\n ----Step 4: create weights---- \n')

  out = WeightsFun(Fnum,Spac)
  Weight = out$Weight
  Near = out$Formers
  Weight = Weight/Spac
  num_w = dimFun(Weight)[2]

  #------------------------------NBI Subproblems-------------------------------#

  # cat('\n ----Step 5: solve NBI sub-problems---- \n')

  # Starting point for first NBI subproblem is the minimizer of f_1(x)
  xstart = c(ShadowX[,1],0)

  Pareto_Fmat = vector()       # Pareto Optima in F-space
  Pareto_Xmat = vector()       # Pareto Optima in X-space
  X_Near      = vector()

  # solve NBI subproblems
  for(k in 1:num_w){

    w  = Weight[,k]

    # Solve problem only if it is not minimizing one of the individual objectives
    indiv_fn_index = which(w == 1)
    # the boundary solution which has been solved

    if(length(indiv_fn_index) != 0){

      # w has a 1 in indiv_fn_index th component, zero in rest
      # Just read in solution from shadow data
      Pareto_Fmat = cbind(Pareto_Fmat, (-PHI[,indiv_fn_index] + ShadowF))
      Pareto_Xmat = cbind(Pareto_Xmat, ShadowX[,indiv_fn_index])
      X_Near = cbind(X_Near, c(ShadowX[,indiv_fn_index],0))

    }else{

      w = rev(w)

      if(Near[k] > 0){

        xstart = X_Near[,Near[k]]
        # start X is the previous weight-order's X

      }

      #start point in F-space
      g_StartF <<- PHI%*%w + ShadowF

      # SOLVE NBI SUBPROBLEM

      out = suppressMessages(nloptr::slsqp(x0 = xstart, fn = myT
                  ,lower = c(VLB,-Inf)
                  ,upper = c(VUB,Inf)
                  ,hin = myCon_ineq
                  ,heq = myTCon_eq))

      x_trial = out$par
      f = out$value
      rm(out)

      Pareto_Fmat = cbind(Pareto_Fmat, -myFM(x_trial[1:nvars])$f)  # Pareto optima in F-space
      som = 0

      for(k in 2:nvars){som = som + x_trial[k]}

      for(k in 2:nvars){x_trial[k] = x_trial[k]/som}

      Pareto_Xmat = cbind(Pareto_Xmat, x_trial[1:nvars])        # Pareto optima in X-space
      X_Near = cbind(X_Near,x_trial)

      }

    }

  #------------------------------Plot Solutions-------------------------------#

#   cat('\n ----Step 6: plot---- \n')

  if(graph==TRUE){plotPareto(Pareto_Fmat, Pareto_Xmat)}

  #------------------------------Output Solutions-------------------------------#

#   Output Solution

  Pareto_Fmat = t(Pareto_Fmat)
  Pareto_Xmat = t(Pareto_Xmat[2:nrow(Pareto_Xmat),])
  colnames(Pareto_Fmat) = c("Ry","AIR1")
  colnames(Pareto_Xmat) = c(paste0(rep("P",(nvars-1)),1:(nvars-1)))

  if(display_solution == TRUE){

    solution = round(cbind(Pareto_Fmat,Pareto_Xmat),3)
    colnames(solution) = c("Ry","AIR1", paste0(rep("P",(nvars-1)),1:(nvars-1)))
    cat("\n Pareto-Optimal Solution \n \n")
    print(solution)

  }else{
    cat("\n Done. \n \n")
  }


  return(list(Pareto_Fmat = round(Pareto_Fmat, 3),
              Pareto_Xmat = round(Pareto_Xmat, 3)))

}

########################### Supporting Functions (A) ########################

# User-Defined Input for NBI.r - Pareto-Optimization via Normal Boundary Intersection

# Input:
## 1) Population correlation matrix (R): criterion & predictor inter-correlation
## 2) Population subgroup difference (d): criterion & predictor mean difference
## between minority and majority subgroups
## 3) Proportion of minority applicants (prop):
## prop = (# of minority applicants)/(total # of applicants)
## 4) Selection ratio (sr): sr = (# of selected applicants)/(total # of applicants)

# Related functions:
# myFM
# myCon

###### combR()######

#' combR
#'
#' Support function to create predictor-criterion matrix
#' @param Rx Predictor inter-correlation matrix
#' @param Ry Predictor-criterion correlation (validity)
#' @return Rxy Predictor-criterion correlation matrix
combR <- function(Rx, Ry){
  cbind(rbind(Rx,c(Ry)),c(Ry,1))
}

###### myFM() ######

#' myFM
#'
#' Supporting function, defines criterion space
#' @param x Input predictor weight vector
#' @return f Criterion vector
myFM = function(x){

  # b = x
  b = x[-1] # Predictor weight
  p_c = x[1] # Cutoff score

  # Obtain within-package 'global' variables
  subd <- subd_ParetoR
  R <- combR(Rx_ParetoR, Rxy_ParetoR)
  R_u <- Rx_ParetoR

  # variance of minority and majority applicant weighted predictor
  # composite (P) distribution (DeCorte, 1999)
  sigma_p = sqrt(t(b)%*%R_u%*%b)

  # mean of Majority weighted predictor composite distribution (DeCorte, 1999)
  subd_sum = colSums(subd)
  p_a_bar0 = subd_sum%*%b/sigma_p

  # Majority group selection ratio
  SR0 = 1 - pnorm(p_c, p_a_bar0)

  p_i_bar = rep(NA,dim(subd)[1])
  SR_i = rep(NA,dim(subd)[1])

  for(i in 1:dim(subd)[1]){

    # mean of Minority_i weighted predictor composite distribution (DeCorte, 1999)
    p_i_bar[i] = (subd_sum-subd[i,])%*%b/sigma_p

    # Minority_i selection ratio (denoted as h_i in DeCorte et al., 1999)
    SR_i[i] = 1 - pnorm(p_c, p_i_bar[i])

  }

  a_g1 = SR_i[1]/SR0


  # Composite Validity R_xy
  Ry = t(c(t(b),0)%*%R%*%c(t(matrix(0,dimFun(R_u)[1],1)),1))/sqrt(t(b)%*%R_u%*%b) # DeCorte et al., 2007

  f = matrix(1,2,1)
  f[1,] = -Ry
  f[2,] = -a_g1

  return(list(f=f))

}

####### myCon_ineq() ######

# Nonlinear inequalities at x

#' myCon_ineq
#'
#' Support function, defines inequal constraint condition
#' @param x Input predictor weight vector
#' @return Inequal constraint condition for use in NBI()
myCon_ineq = function(x){return(vector())}

####### myCon_eq() ######

# Nonlinear equalities at x

#' myCon_eq
#'
#' Support function, defines equal constraint condition
#' @param x Input predictor weight vector
#' @return Equal constraint condition for use in NBI()
myCon_eq = function(x){

  # Obtain within-package 'global' variable
  sr <- sr_ParetoR
  R_u <- Rx_ParetoR

  prop <- prop_ParetoR
  subd <- subd_ParetoR

  b <- x[-1] # Predictor weight
  p_c <- x[1] # Cutoff score

  # variance of minority and majority applicant weighted predictor
  # composite (P) distribution (DeCorte, 1999)
  sigma_p = sqrt(t(b)%*%R_u%*%b)

  # mean of Majority weighted predictor composite distribution (DeCorte, 1999)
  subd_sum = colSums(subd)
  p_a_bar0 = subd_sum%*%b/sigma_p

  # Majority group selection ratio
  SR0 = 1 - pnorm(p_c, p_a_bar0)

  p_i_bar = rep(NA,dim(subd)[1])
  SR_i = rep(NA,dim(subd)[1])

  for(i in 1:dim(subd)[1]){

    # mean of Minority_i weighted predictor composite distribution (DeCorte, 1999)
    p_i_bar[i] = (subd_sum-subd[i,])%*%b/sigma_p

    # Minority_i selection ratio (denoted as h_i in DeCorte et al., 1999)
    SR_i[i] = 1 - pnorm(p_c, p_i_bar[i])

  }

  # Nonlinear equalities at x

  ceq = matrix(1,2,1)
  ceq[1,] = prop%*%SR_i + SR0*(1 - sum(prop)) - sr # DeCorte et al. (2007)
  ceq[2,] = (t(b)%*%R_u%*%b) - 1

  return(ceq)

}
