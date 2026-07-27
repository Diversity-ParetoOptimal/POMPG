#' calc_biserial
#'
#' Calculate the biserial correlation corresponding to a standardized mean
#' difference (d), for use in data generation with subgroup mean differences.
#' @param d Standardized mean difference between subgroups
#' @param p Proportion of the focal (minority) subgroup
#' @return Biserial correlation corresponding to \code{d}
#' @export
calc_biserial <- function(d, p) {
    a <- p * (1 - p)
    lambda <- dnorm(qnorm(1 - p))

    r_pb <- d / sqrt(d^2 + (1 / a))
    r_bis <- r_pb * sqrt(a) / lambda
    return(r_bis)
}

#' calc_d
#'
#' Reverses \code{\link{calc_biserial}}, converting a biserial correlation
#' back to a standardized mean difference (d).
#' @param r Biserial correlation
#' @param p Proportion of the focal (minority) subgroup
#' @return Standardized mean difference corresponding to \code{r}
#' @export
calc_d <- function(r, p) {
    a <- p * (1 - p)
    lambda <- dnorm(qnorm(1 - p))
    d <- lambda * r / sqrt(a * (a - r^2 * lambda^2))
    return(d)
}

#' ai_ratio
#'
#' Convert a standardized mean subgroup difference into an adverse impact
#' ratio, given a selection ratio and subgroup proportion.
#' @param d Standardized mean difference between subgroups
#' @param sr Selection ratio
#' @param p Proportion of the focal (minority) subgroup
#' @return Adverse impact ratio
ai_ratio <- function(d, sr, p)
{
    xcut <- (qnorm(1 - sr) * sqrt(1 + d^2 * (p * (1 - p)))) - (d * p)
    air <- (1.64 * xcut + sqrt(.76 * xcut ^ 2 + 4)) /
        (sqrt(exp(d^2 + 2 * xcut * d)) *
             (1.64 * (xcut + d) + sqrt(.76 * (xcut + d)^2 + 4)))
    return(air)
}

#' quiet
#'
#' Suppress function outputs written via \code{cat}/\code{print}. Acquired
#' from: http://r.789695.n4.nabble.com/Suppressing-output-e-g-from-cat-td859876.html
#' @param x Expression whose printed/cat output should be suppressed
#' @return The (invisible) result of evaluating \code{x}
quiet <- function(x) {
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
}

#' get_sample
#'
#' Draw a stratified sample from subgroup-specific multivariate normal
#' distributions, given subgroup proportions, means, and a shared covariance
#' matrix.
#' @param n Total sample size
#' @param props Vector of subgroup proportions (first element is the reference/majority group)
#' @param mus Matrix of subgroup means (one row per subgroup)
#' @param num_groups Currently unused inside the function body; retained for interface compatibility with the ported simulation code
#' @param sig_mat Shared covariance matrix used for all subgroups
#' @importFrom MASS mvrnorm
#' @return A data frame of sampled observations with subgroup membership columns
#' @export
get_sample <- function(n, props, mus, num_groups, sig_mat) {
    # calculate subgroup sizes
    subn <- vector()
    subn[1] <- NA
    for(i in 2:length(props)){
        subn[i] <- round(n * props[i])
    }
    subn[1] <- n - sum(subn[-1])

    # generate subsamples
    subsamples <- list()
    groups <- vector()
    for(i in 1:length(subn)){
        tmp <- MASS::mvrnorm(n=subn[i], mu=mus[i, ], Sigma=sig_mat)
        groups <- rep(i, subn[i])
        subsamples[[i]] <- cbind(tmp, groups)
    }

    out <- as.data.frame(Reduce(rbind, subsamples))

    for(i in 1:length(subsamples)) {
        out[, paste0('group', i)] <- ifelse(out$groups==i, 1, 0)
    }

    return(out)
}

#' get_weights
#'
#' Obtain PO predictor weights, dispatching to \code{ParetoR::ParetoR()} for
#' a single subgroup or \code{\link{POMA_Focal}} for multiple subgroups.
#' @param sr Selection ratio
#' @param props Subgroup proportion(s)
#' @param ds Subgroup difference(s)
#' @param R Predictor/criterion correlation matrix
#' @param spac Number of Pareto solutions
#' @param graph If TRUE, plots will be generated for the Pareto-optimal curve and predictor weights
#' @param display_solution If TRUE, the Pareto-optimal solution will be displayed
#' @return Matrix of Pareto-optimal predictor weights
#' @export
get_weights <- function(sr, props, ds, R, spac,
                        graph=FALSE,
                        display_solution=FALSE) {

    n_preds <- ncol(R) - 1
    Rx <- R[1:n_preds, 1:n_preds]
    Rxy <- R[n_preds + 1, 1:n_preds]

    if(length(props)==1) {
        output <- quiet(ParetoR::ParetoR(prop=props,
                                Spac=spac,
                                sr=sr,
                                d=ds,
                                R=R,
                                graph=graph,
                                display_solution=display_solution))
    } else {
        output <- quiet(POMA_Focal(sr=sr,
                                   prop=props,
                                   Rx=Rx,
                                   Rxy=Rxy,
                                   subd=ds,
                                   Spac=spac,
                                   graph=graph,
                                   display_solution=display_solution))
    }
    weights <- as.matrix(output$Pareto_Xmat)
    return(weights)
}

#' hire
#'
#' Simulate hiring using a selection ratio and predictor weights, returning
#' composite validity and subgroup adverse impact ratios.
#' @param sr Selection ratio
#' @param b Predictor weight vector
#' @param samp Applicant sample data frame (predictors, then criterion, then a group column)
#' @param group_col Name of the subgroup membership column in \code{samp}
#' @return A list with \code{validity} (composite validity) and \code{airs} (adverse impact ratios relative to the reference group)
#' @importFrom stats cor
#' @export
hire <- function(sr, b, samp, group_col="groups") {
    groups <- unique(samp[, group_col])
    yhat <- as.matrix(samp[, 1:length(b)]) %*% b
    hired <- samp[order(yhat, decreasing=FALSE)[1:ceiling(nrow(samp)* sr)], ]

    sub_srs <- vector()
    for(i in groups) {
        sub_srs[i] <- sum(hired[, group_col]==i) / sum(samp[, group_col]==i)
    }

    validity <- cor(yhat, samp[, length(b) + 1])

    airs <- vector()
    for(i in 2:length(sub_srs)) {
        airs[i - 1] <- sub_srs[i] / sub_srs[1]
    }

    return(list(validity=validity, airs=airs))
}
