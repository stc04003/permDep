#' Permutation test for general dependent truncation
#'
#' Perform permutation test based on conditional or unconditional approach.
#'
#' @param trun is the left truncation time.
#' @param obs is the observed failure time.
#' @param permSize is the number of permutations.
#' @param cens is the status indicator; 0 = censored, 1 = event.
#' @param sampling a character string specifying the sampling method used in permutation.
#' The following are permitted:
#' \describe{
#' \item{\code{conditional}}{conditional permutation;}
#' \item{\code{uconditional}}{unconditional permutation;}
#' \item{\code{is.conditional}}{importance sampling version of conditional permutation;}
#' \item{\code{is.uconditional}}{importance sampling version of unconditional permutation;}
#' }
#' @param kendallOnly,minp1Only,minp2Only optional values indicating which test statistics to be used.
#' If all leave as \code{FALSE}, \code{permDep} will use all three test statistics in each permutation.
#' @param nc is the number of cores used in permutation.
#' When \code{nc > 1}, permutation is carried out with parallel computing.
#' @param seed an optional vector containing random seeds to be used to generate permutation samples.
#' Random seeds will be used when left unspecified.
#' @param minp.eps an optional value indicating the width of the intervals used in minp2 procedure.
#' The following input are allowed:
#' \describe{
#' \item{a single numerical number}{sets the common width for all intervals}
#' \item{a numerical vector}{the length of the vector needs to be a sub-multiple or
#' multiple of the number of uncensored events. This input is useful for user specified widths.}
#' \item{a character string}{This option can be either "in" or "out", which require the minimum number
#' of uncensored events inside or outside of each interval, respectively.}
#' \item{NULL}{automatic selection using the algorithm outlined in Chiou (2018) will be used. }
#' }
#' @param plot.int an optional logical value indicating whether an animated scaterplot will be produced
#' to how the minp intervals are choosen for the observed data. 
#' 
#' @return A list containing output with the following components:
#' \describe{
#' \item{obsKen}{the observed p-value using Kendall's tau test statistic.}
#' \item{obsP1}{the observed p-value using minp1 test statistic.}
#' \item{obsP2}{the observed p-value using minp2 test statistic.}
#' \item{obsTest1}{the observed minp1 test statistic.}
#' \item{obsTest2}{the observed minp2 test statistic.}
#' \item{permKen}{Kendall's tau test statistics from permutation samples.}
#' \item{permP1}{minp1 test statistics from permutation samples.}
#' \item{permP2}{minp2 test statistics from permutation samples.}
#' }
#'
#' @references Chiou, S.H., Qian, J., and Betensky, R.A. (2018).
#' Permutation Test for General Dependent Truncation.
#' \emph{Computational Statistics \& Data Analysis},128, p308--324.
#' 
#' @importFrom BB spg
#' @importFrom survival Surv survfit basehaz coxph
#' @importFrom stats model.matrix pchisq pnorm rexp runif var
#' @importFrom parallel detectCores parLapply makeCluster clusterExport stopCluster
#'
#' @export
#' @keywords permDep
#' 
#' @examples
#' simDat <- function(n) {
#'   k <- s <- 1
#'   tt <- xx <- yy <- cc <- delta <- rep(-1, n)
#'   while(k <= n){
#'     tt[k] <- runif(1, 0, 3.5)
#'     xx[k] <- 1.95 + 0.65 * (tt[k] - 1.25)^2 + rnorm(1, sd = 0.1)
#'     cc[k] <- runif(1, 0, 10)
#'     delta[k] <- (xx[k] <= cc[k])
#'     yy[k] <- pmin(xx[k], cc[k])
#'     s <- s + 1
#'     if(tt[k] <= yy[k]) k = k+1
#'   }
#'  data.frame(list(trun = tt, obs = yy, delta = delta))
#' }
#'
#' set.seed(123)
#' dat <- simDat(50)
#' B <- 20
#'
#' ## Perform conditional permutation with Kendall's tau, minp1 and minp2
#' set.seed(123)
#' system.time(fit <- with(dat, permDep(trun, obs, B, delta, nc = 1)))
#' fit
permDep <- function(trun, obs, permSize, cens,
                    sampling = c("conditional", "unconditional", "is.conditional", "is.unconditional"),
                    kendallOnly = FALSE, minp1Only = FALSE, minp2Only = FALSE,
                    nc = ceiling(detectCores() / 2), seed = NULL,
                    minp.eps = NULL, plot.int = FALSE) {
    sampling <- match.arg(sampling)
    n <- length(obs)
    ## if (!(sampling %in% c("cond", "ucond", "ucond1", "ucond2",
    ##                       "ucond3", "ucond4", "iscond", "isucond")))
    obsKen <- obsKenp <- obs1 <- obs2 <- obsP1 <- obsP2 <- obsTest1 <- obsTest2 <- p1 <- p2 <- NULL
    if (sum(kendallOnly, minp1Only, minp2Only) == 0) {
        tmp <- getKendall(trun, obs, cens)
        obsKen <- tmp[1]
        obsKenp <- 2 - 2 * pnorm(abs(tmp[1] / sqrt(tmp[2])))
        obs1 <- getMinP(trun, obs, cens, eps = minp.eps, plotInt = plot.int)
        obsP1 <- obs1$minP
        obsTest1 <- obs1$maxTest
        p1 <- obs1$p
        obs2 <- getMinP(trun, obs, cens, minp1 = FALSE, eps = minp.eps, plotInt = plot.int)
        obsP2 <- obs2$minP
        obsTest2 <- obs2$maxTest
        p2 <- obs2$p
    }
    if (sum(kendallOnly, minp1Only, minp2Only) > 0) {
        if (kendallOnly) {
            tmp <- getKendall(trun, obs, cens)
            obsKen <- tmp[1]
            obsKenp <- 2 - 2 * pnorm(abs(tmp[1] / sqrt(tmp[2])))
        }
        if (minp1Only) {        
            obs1 <- getMinP(trun, obs, cens, eps = minp.eps, plotInt = plot.int)
            obsP1 <- obs1$minP
            obsTest1 <- obs1$maxTest
            p1 <- obs1$p
        }
        if (minp2Only) {
            obs2 <- getMinP(trun, obs, cens, minp1 = FALSE, eps = minp.eps, plotInt = plot.int)
            obsP2 <- obs2$minP
            obsTest2 <- obs2$maxTest
            p2 <- obs2$p
        }
    }
    MinP1Time <- MinP2Time <- matrix(NA, permSize, length(trun))
    if (is.null(seed) || length(seed) != permSize) seed <- sample(.Machine$integer.max, size = permSize)
    if (nc == 1) {
        permKen <- permKenp <- permP1 <- permP2 <- permTest1 <- permTest2<- rep(NA, permSize)
        for(i in 1:permSize){
            perm.data <- getPerm(trun, obs, cens, sampling, seed[i])    
            if (sum(kendallOnly, minp1Only, minp2Only) > 0) {
                if (kendallOnly) {
                    tmp <- getKendall(perm.data$trun, perm.data$obs, perm.data$cens)
                    permKen[i] <- tmp[1]
                    permKenp[i] <- 2 - 2 * pnorm(abs(tmp[1] / sqrt(tmp[2])))
                }
                if (minp1Only) {
                    p1fit <- getMinP(perm.data$trun, perm.data$obs, perm.data$cens, eps = minp.eps)
                    permTest1[i] <- p1fit$maxTest
                }
                if (minp2Only) {
                    p2fit <- getMinP(perm.data$trun, perm.data$obs, perm.data$cens,
                                     minp1 = FALSE, eps = minp.eps)
                    permTest2[i] <- p2fit$maxTest
                }
            }
            if (sum(kendallOnly, minp1Only, minp2Only) == 0) {
                kendallOnly <- minp1Only <- minp2Only <- TRUE
                tmp <- getKendall(perm.data$trun, perm.data$obs, perm.data$cens)
                permKen[i] <- tmp[1]
                permKenp[i] <- 2 - 2 * pnorm(abs(tmp[1] / sqrt(tmp[2])))
                p1fit <- getMinP(perm.data$trun, perm.data$obs, perm.data$cens, eps = minp.eps)
                ## obsTest = ifelse(surv, NA, obsTest1))
                p2fit <- getMinP(perm.data$trun, perm.data$obs, perm.data$cens, minp1 = FALSE, eps = minp.eps)
                ## obsTest = ifelse(surv, NA, obsTest2), minp1 = FALSE)  
                permTest1[i] <- p1fit$maxTest
                permTest2[i] <- p2fit$maxTest
                ## MinP1Time[i,] <- p1fit$MinPTime
                ## MinP2Time[i,] <- p2fit$MinPTime
            }
        }
    }    
    if (nc > 1) {
        permKen <- permP1 <- permP2 <- permTest1 <- permTest2<- rep(NULL, permSize)
        cl <- makeCluster(nc)
        clusterExport(cl = cl, varlist=c("trun", "obs", "cens", "sampling", "seed", "obsTest1",
                                         "obsTest2", "minp.eps"),
                      envir = environment())
        clusterExport(cl = cl, varlist=c("getPerm", "getKendall", "getMinP"), envir = environment())
        out <- unlist(parLapply(cl, 1:permSize, function(x) {
            perm.data <- getPerm(trun, obs, cens, sampling, seed[x])
            if (sum(kendallOnly, minp1Only, minp2Only) > 0) {
                tmpK <- tmpP1 <- tmpP2 <- NULL
                if (kendallOnly)
                    tmpK <- getKendall(perm.data$trun, perm.data$obs, perm.data$cens)[1]
                if (minp1Only) 
                    tmpP1 <- getMinP(perm.data$trun, perm.data$obs, perm.data$cens, eps = minp.eps)
                if (minp2Only)                                     
                    tmpP2 <- getMinP(perm.data$trun, perm.data$obs, perm.data$cens,
                                     minp1 = FALSE, eps = minp.eps)
                return(c(tmpK, tmpP1$maxTest, tmpP2$maxTest))
            }
            if (sum(kendallOnly, minp1Only, minp2Only) == 0) {
                kendallOnly <- minp1Only <- minp2Only <- TRUE
                tmpK <- getKendall(perm.data$trun, perm.data$obs, perm.data$cens)[1]
                tmpP1 <- getMinP(perm.data$trun, perm.data$obs, perm.data$cens, eps = minp.eps)
                ## obsTest = ifelse(surv, 1e5, obsTest1))
                tmpP2 <- getMinP(perm.data$trun, perm.data$obs, perm.data$cens,
                                 minp1 = FALSE, eps = minp.eps)  
                ## obsTest = ifelse(surv, 1e5, obsTest2), minp1 = FALSE)  
                return(c(tmpK, tmpP1$maxTest, tmpP2$maxTest))
            }}))
        stopCluster(cl)
        out <- matrix(out, nrow = ifelse(sum(kendallOnly, minp1Only, minp2Only) == 0,
                                         3, sum(kendallOnly, minp1Only, minp2Only)))
        tmp <- matrix(NA, 3, permSize)
        if (sum(kendallOnly, minp1Only, minp2Only) == 0) {
            tmp <- out
            kendallOnly <- minp1Only <- minp2Only <- TRUE
        } else {
            tmp[(1:3)[c(kendallOnly, minp1Only, minp2Only)],] <- out
        }
        permKen <- tmp[1,]
        permTest1 <- tmp[2,]
        permTest2 <- tmp[3,]
    }
    if (!is.null(permTest1) && all(!is.na(permTest1))) {
        ## perm.data.P1 <- getPerm(trun, obs, cens, sampling, seed[which.min(permTest1)])
        permP1 <- 1 - pchisq(permTest1, df = 1)
    }
    if (!is.null(permTest2) && all(!is.na(permTest2))) {
        ## perm.data.P2 <- getPerm(trun, obs, cens, sampling, seed[which.min(permTest2)])
        permP2 <- 1 - pchisq(permTest2, df = 1)
    }
    out <- list(obsKen = obsKen, obsKenp = obsKenp, obsP1 = obsP1, obsP2 = obsP2,
                obsTest1 = obsTest1, obsTest2 = obsTest2,
                permKen = permKen, permP1 = permP1, permP2 = permP2,
                permTest1 = permTest1, permTest2 = permTest2,
                p.valueKen = ifelse(kendallOnly, (sum(abs(obsKen) < abs(permKen)) + 1)/(length(permKen) + 1), NA),
                p.valueMinp1 = ifelse(minp1Only, (sum(obsP1 > permP1) + 1) / (length(permKen) + 1), NA),
                p.valueMinp2 = ifelse(minp2Only, (sum(obsP2 > permP2) + 1) / (length(permKen) + 1), NA),
                kendallOnly = kendallOnly, minp1Only = minp1Only, minp2Only = minp2Only,
                sampling = sampling, permSize = permSize, random.seed = seed, p1 = p1, p2 = p2)
    class(out) <- "permDep"
    return(out)
}

mysample <- function(x, n) {
    res <- vector("integer", length(x))
    .C("mysampleC", as.integer(x), as.integer(n),
       as.integer(length(x)), out = as.integer(res), PACKAGE = "permDep")$out
}

getKendall <- function(trun, obs, cens = NULL) {
    if (is.null(cens)) cens <- rep(1, length(trun))
    res <- vector("double", 3)
    .C("kendallTrun", as.double(trun), as.double(obs), as.double(cens), as.integer(length(trun)), 
       out = as.double(res), PACKAGE = "permDep")$out
    
}

getKendallWgt <- function(trun, obs, cens = NULL, scX = NULL, scT = NULL) {
    if (is.null(cens)) cens <- rep(1, length(trun))
    if (is.null(scX)) scX <- rep(1, length(trun))
    if (is.null(scT)) scT <- rep(1, length(trun))
    res <- vector("double", 2)
    .C("kendallTrunWgt", as.double(trun), as.double(obs), as.double(cens), as.integer(length(trun)),
       as.double(scX), as.double(scT),
       out = as.double(res), PACKAGE = "permDep")$out
}

getScore <- function(x, y) {
    n <- nrow(y)
    nvar <- ncol(x)
    start <- y[, 1]
    stopp <- y[, 2]
    event <- y[, 3]
    sort.end <- order(-stopp, event) - 1L
    sort.start <- order(-start) - 1L
    newstrat <- n
    offset <- rep(0, n)
    weights <- rep(1, n)
    maxiter <- 20
    eps <- 1e-09
    toler.chol <- 1.8e-12
    storage.mode(y) <- storage.mode(x) <- "double"
    storage.mode(offset) <- storage.mode(weights) <- "double"
    storage.mode(newstrat) <- "integer"
    agfit <- .Call("agfit4", y, x, newstrat, weights, offset, 
                   as.double(rep(0, nvar)), sort.end, sort.start, as.integer(1),
                   as.integer(maxiter), as.double(eps), as.double(toler.chol), as.integer(1))
    agfit$sctest
}

## control <- list(eps = 1e-09, toler.chol = 1.818989e-12, iter.max = 20, toler.inf = 3.162278e-05, outer.max = 10)

#' @importFrom ggplot2 ggplot geom_point scale_color_manual labs xlab ylab aes
#' @importFrom gganimate transition_time animate
#' 
getMinP <- function(trun, obs, cens, obsTest = NA, minp1 = TRUE,
                    eps = NULL, plotInt = FALSE) {
    E <- min(10, round(0.2 * sum(cens)), round(0.2 * length(obs)))
    n <- length(obs)
    p <- NULL
    data0 <- data.frame(cbind(trun, obs, cens))[order(trun),]
    survObj <- with(data0, Surv(trun, obs, cens))
    log.rank.test <- log.rank.pval <- rep(NA, n)
    if (minp1) {
        for (j in E:(n-E)) {
            group <- rep(1, n)
            group[-(1:j)] <- 2
            if (sum(table(data0$cens, group)["1",] < E) == 0) {
                log.rank.test[j] <- getScore(model.matrix(~ factor(group) - 1), survObj)
                if (!is.na(obsTest) && !is.na(log.rank.test[j]) && log.rank.test[j] >= obsTest)
                    break
            }
        }
    }
    if (plotInt) groups <- NULL
    if (!minp1) {
        trun.tmp <- with(data0, trun[cens == 1])
        if (is.null(eps)) {
            eps1 <- with(data0, sapply(1:n, function(x) sort(abs(trun[x] - trun.tmp))[sum(cens) - E]))
            eps2 <- with(data0, sapply(1:n, function(x) sort(abs(trun[x] - trun.tmp))[E]))
            eps <- pmax(eps1, eps2)
            eps <- rep(min(eps), length(eps))
        }
        if ((eps == "in")[1])  ## interior has at least sum(cens) - E events
            eps <- with(data0, sapply(1:n, function(x) sort(abs(trun[x] - trun.tmp))[sum(cens) - E]))
        if ((eps == "out")[1])  ## interior has at least E events
            eps <- with(data0, sapply(1:n, function(x) sort(abs(trun[x] - trun.tmp))[E])) ## interior has E events
        if (length(eps) != n) {
            if (n %% length(eps) > 0)
                print("eps length is not a sub-multiple or multiple of the number of uncensored events.")
            eps <- rep(eps, n)[1:n]
        }
        for (j in 1:n) {
            group <- rep(1, n)
            group[abs(data0[,"trun"] - data0[j,"trun"]) <= eps[j]] <- 2
            if (plotInt) groups[[j]] <- group
            if (length(unique(group)) > 1 & sum(table(data0$cens, group)["1",] < E) == 0) {
                log.rank.test[j] <- getScore(model.matrix(~ factor(group) - 1), survObj)                
                if (!is.na(obsTest) && !is.na(log.rank.test[j]) && log.rank.test[j] >= obsTest)
                    break
            }
            if (length(unique(group)) == 1 || sum(table(data0$cens, group)["1",] < E) > 0) {
                log.rank.test[j] <- -999
                if (plotInt) groups[[j]] <- rep(NA, n)
            }
        }
    }
    if (max(log.rank.test, na.rm = TRUE) < -990)
        stop("minp.eps yields invalid insufficient intervals.")
    else log.rank.pval <- 1 - pchisq(log.rank.test, df = 1)
    if (plotInt) {        
        dat <- data.frame(trun = rep(data0$trun,n), 
                          obs = rep(data0$obs, n),
                          grp = unlist(groups),
                          stp = rep(1:n, each = n))
        dat <- dat[complete.cases(dat),]
        p <- ggplot(dat, aes(x = trun, y = obs, color = factor(grp))) +
            geom_point(show.legend = FALSE, alpha = 0.7, cex = 2) +
            scale_color_manual(breaks = c("1", "2"), values = c("black", "red")) +
            transition_time(stp) +
            labs(title = "Cutpoint: {frame_time}") +
            xlab("Truncation time") +
            ylab("Survival time")
        print(animate(p, start_pause = 10, end_pause = 10))
    }
    list(maxP = max(log.rank.pval, na.rm = TRUE),
         minP = min(log.rank.pval, na.rm = TRUE),
         maxTest = max(log.rank.test, na.rm = TRUE),
         minTest = min(log.rank.test, na.rm = TRUE),
         p = p)
}

getPerm <- function(trun, obs, cens, sampling, seed = NULL) {
    n <- length(obs)
    tt <- trun
    if (!is.null(seed))
        set.seed(seed)
    if (sampling == "conditional") {
        obs.sort <- sort(obs)
        cens.sort <- cens[order(obs)]
        trun.sort <- sort(trun)
        perm.trun.index <- rep(NULL, n)
        m <- colSums((matrix(rep(obs.sort, each=n), n,n) - trun.sort) > 0)
        for (j in 1:n){
            set.1 <- perm.trun.index[1:j] ## truncation time being selected for X[1:(j-1)]
            set.2 <- setdiff((1:m[j]), set.1)
            perm.trun.index[j] <- ifelse(length(set.2) == 1, set.2, sample(set.2, 1))
        }
        perm.trun <- trun.sort[perm.trun.index]
        out <- list(trun = perm.trun, obs = obs.sort, cens = cens.sort)
    }
    if (sampling == "unconditional") {
        perm.trun.index <- sample(1:n, n) 
        perm.trun.0 <- trun[perm.trun.index]
        out <- list(trun = perm.trun.0, obs = obs, cens = cens)
    }
    if (sampling == "unconditional1") {
        perm.trun.index <- sample(1:n, n) 
        perm.trun.0 <- trun[perm.trun.index]
        out <- list(trun = perm.trun.0, obs = obs, cens = cens)
        out <- subset(data.frame(out), trun < obs & tt < trun)
    }
    if (sampling == "unconditional2") {
        perm.trun.index <- sample(1:n, n) 
        perm.trun.0 <- trun[perm.trun.index]
        out <- list(trun = perm.trun.0, obs = obs, cens = cens)
        out$trun <- pmax(tt, out$trun)
    }
    if (sampling == "unconditional3") {
        perm.trun.index <- sample(1:n, n, replace = TRUE)
        perm.trun.0 <- trun[perm.trun.index]
        out <- list(trun = perm.trun.0, obs = obs, cens = cens)        
    }
    if (sampling == "unconditional4") {
        perm.trun.index <- unlist(lapply(sapply(1:n, function(x) which(trun[x] < obs)), function(y) ifelse(length(y) == 1, y, sample(y, 1))))
        perm.trun.0 <- trun[perm.trun.index]
        out <- list(trun = perm.trun.0, obs = obs, cens = cens)        
    }
    if (sampling == "is.conditional") {
        A <- 1 * sapply(obs, function(x) x > trun)
        trun <- trun[order(rowSums(A))]
        for (maxit in 1:50) {
            A2 <- A[order(rowSums(A)),]
            sp <- NULL
            for (i in 1:n) {
                if (sum(A2[i,] > 0) > 0) {
                    sz <- which(A2[i,] > 0)
                    prob <- exp((n - i + 1 - colSums(A2)) / (n - i))
                    if (length(sz) == 1) sp <- c(sp, sz)
                    else sp <- c(sp, sample(sz, 1, TRUE, prob[sz]))
                    } else break
                A2[,sp] <- 0
                A2[i,] <- 0
            }
            if (length(sp) == n)
                break
        }
        ## print(maxit)
        if (maxit == 50)
            stop("Permutation sample not found", call. = FALSE)
        out <- list(trun = trun, obs = obs[sp], cens = cens[sp])
    }
    if (sampling == "is.unconditional") {
        A <- 1 * sapply(obs, function(x) x > trun)
        trun <- trun[order(rowSums(A))]
        A <- A[order(rowSums(A)),]
        sp <- NULL
        for (i in 1:n) {
            prob <- ifelse(colSums(A) != 0, exp((n - i + 1 - colSums(A)) / (n - i)), 0)
            if (i == n) sp <- c(sp, which(!(1:n %in% sp)))
            else sp <- c(sp, sample(1:n, 1, TRUE, prob))
            A[,sp] <- 0
            A[i,] <- 0            
        }
        out <- list(trun = trun, obs = obs[sp], cens = cens[sp])        
    }
    out <- subset(data.frame(out), trun < obs )
    out[order(out[,1]),]
}
