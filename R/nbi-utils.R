########################### Supporting Functions (B) ########################

# Supplementary Functions for NBI.r - Pareto-Optimization via Normal Boundary Intersection

# Function List
## assert_col_vec
## dimFun
## WeightsFun
## Weight_Generate
## myLinCom
## myT
## myTCon_eq
## plotPareto

###### assert_col_vec() ######

#' assert_col_vec
#'
#' Support function, refines intermediate variable for use in NBI()
#' @param v Intermediate variable v
#' @return Refined variable v
assert_col_vec = function(v){
  if(is.null(dimFun(v))){
    v=v
  }else if(dimFun(v)[1] < dimFun(v)[2]){v = t(t)}
  return(v)}

###### dimFun() ######

#' dimFun
#'
#' Support function, checks input predictor weight vector x
#' @param x Input predictor weight vector
#' @return x Checked and refined input predictor weight vector
dimFun = function(x){
  if(is.null(dim(x))){
    return(c(0,0))
  }else(return(dim(x)))
}

###### WeightsFun() ######

#' WeightsFun
#'
#' Support function, generates all possible weights for NBI subproblems
#' @param n Number of objects (i.e., number of predictor and criterion)
#' @param k Number of Pareto points
#' @return Weights All possible weights for NBI subproblem
WeightsFun = function(n, k){

  # global variables
  # weight, Weights, Formers, Layer, lastone, currentone
  #
  # Generates all possible weights for NBI subproblems given:
  # n, the number of objectives
  # 1/k, the uniform spacing between two w_i (k integral)
  # This is essentially all the possible integral partitions
  # of integer k into n parts.

  WeightSub <<- matrix(0,1,n)
  Weights <<- vector()
  # assign("Formers", vector(), envir = .GlobalEnv)
  Formers <<- vector()
  # assign("Layer", n, envir = .GlobalEnv)
  Layer <<- n
  # assign("lastone", vector(), envir = .GlobalEnv)
  lastone <<- -1
  # assign("currentone", -1, envir = .GlobalEnv)
  currentone <<- -1

  Weight_Generate(1, k)

  return(list(Weights = Weights, Formers = Formers))

}

###### Weight_Generate() ######

#' Weight_Generate
#'
#' Function intended to test the weight generation scheme for NBI for > 2 objectives
#' @param n Number of objects (i.e., number of predictor and criterion)
#' @param k Number of Pareto points
#' @return Weight_Generate
Weight_Generate = function(n, k){

  # global variables:
  # weight Weights Formers Layer lastone currentone

  # wtgener_test(n,k)
  #
  # Intended to test the weight generation scheme for NBI for > 2 objectives
  # n is the number of objectives
  # 1/k is the uniform spacing between two w_i (k integral)

  if(n == Layer){

    if(currentone >= 0){
      Formers <<- c(Formers,lastone)
      lastone <<- currentone
      currentone <<- -1
    }else{
      num = dimFun(Weights)[2]
      Formers <<- c(Formers,num)
    }

    WeightSub[(Layer - n + 1)] <<- k
    Weights <<- cbind(Weights,t(WeightSub))

  }else{

    for(i in 0:k){
      if(n == (Layer - 2)){
        num = dimFun(Weights)[2]
        currentone <<- num+1
      }

      WeightSub[(Layer - n + 1)] <<- i
      Weight_Generate(n+1, k-i)
    }

  }

}

###### myLinCom() ######

#' myLincom
#'
#' Support function
#' @param x Input predictor weight vector
#' @return f Criterion vector
myLinCom = function(x){

  # global variable: g_Weight
  f0 = myFM(x)$f
  f = t(g_Weight)%*%f0
  return(f)

}

###### myT() ######

#' myT
#'
#' Support function, define criterion space for intermediate step in NBI()
#' @param x_t Temporary input weight vector
#' @return f Temporary criterion space
myT = function(x_t){

  f = x_t[length(x_t)]
  return(f)

}

###### myTCon_eq() ######

#' myTCon_eq
#'
#' Support function, define constraint condition for intermediate step in NBI()
#' @param x_t Temporary input weight vector
#' @return ceq Temporary constraint condition
myTCon_eq = function(x_t){

  # global variables:
  # g_Normal g_StartF

  t = x_t[length(x_t)]
  x = x_t[1:(length(x_t)-1)]

  fe  = -myFM(x)$f - g_StartF - t * g_Normal

  # c = myCon_ineq(x)
  ceq1 = myCon_eq(x)
  ceq = c(ceq1,fe)

  return(ceq)

}

###### plotPareto() ######

#' plotPareto
#'
#' Function for plotting Pareto-optimal curve and predictor weights
#' @param CriterionOutput Pareto-Optimal criterion solution
#' @param ParetoWeights Pareto-Optimal predictor weight solution
#' @return Plot of Pareto-optimal curve and plot of predictor weights
plotPareto = function(Pareto_Fmat, Pareto_Xmat){

  par(mfrow=c(1,2))

  Ry = t(Pareto_Fmat[1,])
  AIR1 = t(Pareto_Fmat[2,])
  X = t(Pareto_Xmat[2:nrow(Pareto_Xmat),])

  # AI ratio - Composite Validity trade-off

  plot(AIR1, Ry,
       xlim = c(min(AIR1),max(AIR1)),
       main = "Composite Validity - AI ratio trade-off",
       xlab = "AIR1",
       ylab = "Composite Validity",
       type='c',col='blue')

  points(AIR1, Ry,
         pch=8,col='red')

  # Predictor weights

  plot(AIR1,X[,1],
       xlim=c(min(AIR1),max(AIR1)),ylim=c(0,1),
       main = "Predictor weights trade-off function",
       xlab = "AIR1",
       ylab = "Predictor weight",
       type='c',col='red')
  points(AIR1, X[,1], pch=8, col=rainbow(1))

  for(i in 2:ncol(X)){

    lines(AIR1, X[,i],type='c',
          col=rainbow(1, start=((1/ncol(X))*(i-1)), alpha=1))
    points(AIR1,X[,i],pch=8,
           col=rainbow(1, start=((1/ncol(X))*(i-1)), alpha=1))

  }

  legend('topleft',
         legend=c(paste0('Predictor ',1:ncol(X))),
         lty=c(rep(2,ncol(X))),lwd=c(rep(2,ncol(X))),
         col=rainbow(ncol(X)))

}
