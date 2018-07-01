mysample <- function(x, n) {
    res <- vector("integer", length(x))
    .C("mysample", as.integer(x), as.integer(n),
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
    ## x is a covariate matrix
    ## y is a survival object with left truncation, right censoring
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

getMinP <- function(trun, obs, cens, obsTest = NA, minp1 = TRUE, eps = NULL) {
    ## trun (n by 1) is the left truncation times
    ## obs (n by 1) is the observed failure time
    ## cens (n by 1) is the censoring indicator: 1-failure, 0-censored
    E <- min(10, round(0.2 * sum(cens)), round(0.2 * length(obs)))
    n <- length(obs)
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
    if (!minp1) {
        if (is.null(eps)) {
            trun.tmp <- with(data0, trun[cens == 1])
            ## eps <- with(data0, sapply(1:n, function(x) sort(abs(trun[x] - trun.tmp))[E])) ## interior has E events
            eps <- with(data0, sapply(1:n, function(x) sort(abs(trun[x] - trun.tmp))[sum(cens) - E])) ## interior has sum(cens) - E events
            ## trun.tmp <- trun[cens == 1]
            ## trun.tmp <- trun.tmp[order(trun.tmp)]
            ## eps <- max(sapply(E:(n-E), function(x) sort(abs(trun.tmp[x] - trun.tmp))[E]), na.rm = TRUE)
            ## tmp <- obs[cens == 1]
            ## eps <- max(diff(tmp[order(tmp)], E))
            ## IQR(data0[, "trun"]) / 2
        }
        for (j in 1:n) {
            group <- rep(1, n)
            group[abs(data0[,"trun"] - data0[j,"trun"]) <= eps[j]] <- 2
            if (length(unique(group)) > 1 & sum(table(data0$cens, group)["1",] < E) == 0) {
                log.rank.test[j] <- getScore(model.matrix(~ factor(group) - 1), survObj)
                if (!is.na(obsTest) && !is.na(log.rank.test[j]) && log.rank.test[j] >= obsTest)
                    break
            }            
        }
    }
    log.rank.pval <- 1 - pchisq(log.rank.test, df = 1)
    ## print(sum(!is.na(log.rank.test)))
    list(## MinPTime = MinPTime,
         maxP = max(log.rank.pval, na.rm = TRUE),
         minP = min(log.rank.pval, na.rm = TRUE),
         maxTest = max(log.rank.test, na.rm = TRUE),
         minTest = min(log.rank.test, na.rm = TRUE))
}

getPerm <- function(trun, obs, cens, sampling, seed = NULL) {
    n <- length(obs)
    tt <- trun
    if (!is.null(seed))
        set.seed(seed)
    if (sampling == "cond") {
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
    if (sampling == "ucond") {
        perm.trun.index <- sample(1:n, n) 
        perm.trun.0 <- trun[perm.trun.index]
        out <- list(trun = perm.trun.0, obs = obs, cens = cens)
    }
    if (sampling == "ucond1") {
        perm.trun.index <- sample(1:n, n) 
        perm.trun.0 <- trun[perm.trun.index]
        out <- list(trun = perm.trun.0, obs = obs, cens = cens)
        out <- subset(data.frame(out), trun < obs & tt < trun)
    }
    if (sampling == "ucond2") {
        perm.trun.index <- sample(1:n, n) 
        perm.trun.0 <- trun[perm.trun.index]
        out <- list(trun = perm.trun.0, obs = obs, cens = cens)
        out$trun <- pmax(tt, out$trun)
    }
    if (sampling == "ucond3") {
        perm.trun.index <- sample(1:n, n, replace = TRUE)
        perm.trun.0 <- trun[perm.trun.index]
        out <- list(trun = perm.trun.0, obs = obs, cens = cens)        
    }
    if (sampling == "ucond4") {
        perm.trun.index <- unlist(lapply(sapply(1:n, function(x) which(trun[x] < obs)), function(y) ifelse(length(y) == 1, y, sample(y, 1))))
        perm.trun.0 <- trun[perm.trun.index]
        out <- list(trun = perm.trun.0, obs = obs, cens = cens)        
    }
    if (sampling == "iscond") {
        A <- 1 * sapply(obs, function(x) x > trun)
        trun <- trun[order(rowSums(A))]
        for (maxit in 1:50) {
            A2 <- A[order(rowSums(A)),]
            sp <- NULL
            for (i in 1:n) {
                ## prob <- exp((n - colSums(A)) / (n - 1))
                ## if (max(0, prob[sz], na.rm = TRUE) > 0 & sum(A2[i,] > 0) > 0)
                if (sum(A2[i,] > 0) > 0) {
                    sz <- which(A2[i,] > 0)
                    prob <- exp((n - i + 1 - colSums(A2)) / (n - i))
                    if (length(sz) == 1) 
                        sp <- c(sp, sz)
                    else
                        sp <- c(sp, sample(sz, 1, TRUE, prob[sz]))
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
    if (sampling == "isucond") {
        A <- 1 * sapply(obs, function(x) x > trun)
        trun <- trun[order(rowSums(A))]
        A <- A[order(rowSums(A)),]
        sp <- NULL
        for (i in 1:n) {
            ## prob <- exp((n - colSums(A)) / (n - 1))
            prob <- ifelse(colSums(A) != 0, exp((n - i + 1 - colSums(A)) / (n - i)), 0)
            if (i == n)
                sp <- c(sp, which(!(1:n %in% sp)))
            else
                sp <- c(sp, sample(1:n, 1, TRUE, prob))
            A[,sp] <- 0
            A[i,] <- 0            
        }
        out <- list(trun = trun, obs = obs[sp], cens = cens[sp])        
    }
    out <- subset(data.frame(out), trun < obs )
    out[order(out[,1]),]
}

permDep <- function(trun, obs, permSize, cens, sampling = "cond", 
                     kendallOnly = FALSE, minp1Only = FALSE, minp2Only = FALSE,
                     nc = ceiling(detectCores() / 2), seed = NULL) {
    ## obs: the failure time, which is subject to left truncation by trun and right censroing
    ## observable:  obs >= trun
    ## trun: the truncation time
    ## permSize:  number of random samples from permutation distribution
    n <- length(obs)
    ## error check
    if (!(sampling %in% c("cond", "ucond", "ucond1", "ucond2", "ucond3", "ucond4", "iscond", "isucond")))
        stop("Invalid sampling method", call. = FALSE)
    obsKen <- obsKenp <- obs1 <- obs2 <- obsP1 <- obsP2 <- obsTest1 <- obsTest2 <- NULL
    if (sum(kendallOnly, minp1Only, minp2Only) == 0) {
        tmp <- getKendall(trun, obs, cens)
        obsKen <- tmp[1]
        obsKenp <- 2 - 2 * pnorm(abs(tmp[1] / sqrt(tmp[2])))
        obs1 <- getMinP(trun, obs, cens)
        obsP1 <- obs1$minP
        obsTest1 <- obs1$maxTest
        obs2 <- getMinP(trun, obs, cens, minp1 = FALSE)
        obsP2 <- obs2$minP
        obsTest2 <- obs2$maxTest
    }
    if (sum(kendallOnly, minp1Only, minp2Only) > 0) {
        if (kendallOnly) {
            tmp <- getKendall(trun, obs, cens)
            obsKen <- tmp[1]
            obsKenp <- 2 - 2 * pnorm(abs(tmp[1] / sqrt(tmp[2])))
        }
        if (minp1Only) {        
            obs1 <- getMinP(trun, obs, cens)
            obsP1 <- obs1$minP
            obsTest1 <- obs1$maxTest
        }
        if (minp2Only) {
            obs2 <- getMinP(trun, obs, cens, minp1 = FALSE)
            obsP2 <- obs2$minP
            obsTest2 <- obs2$maxTest
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
                    p1fit <- getMinP(perm.data$trun, perm.data$obs, perm.data$cens)
                    permTest1[i] <- p1fit$maxTest
                }
                if (minp2Only) {
                    p2fit <- getMinP(perm.data$trun, perm.data$obs, perm.data$cens, minp1 = FALSE)
                    permTest2[i] <- p2fit$maxTest
                }
            }
            if (sum(kendallOnly, minp1Only, minp2Only) == 0) {
                kendallOnly <- minp1Only <- minp2Only <- TRUE
                tmp <- getKendall(perm.data$trun, perm.data$obs, perm.data$cens)
                permKen[i] <- tmp[1]
                permKenp[i] <- 2 - 2 * pnorm(abs(tmp[1] / sqrt(tmp[2])))
                p1fit <- getMinP(perm.data$trun, perm.data$obs, perm.data$cens)
                                 ## obsTest = ifelse(surv, NA, obsTest1))
                p2fit <- getMinP(perm.data$trun, perm.data$obs, perm.data$cens, minp1 = FALSE)
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
        clusterExport(cl = cl, varlist=c("trun", "obs", "cens", "sampling", "seed", "obsTest1", "obsTest2"),
                      envir = environment())
        clusterExport(cl = cl, varlist=c("getPerm", "getKendall", "getMinP"), envir = environment())
        out <- unlist(parLapply(cl, 1:permSize, function(x) {
            perm.data <- getPerm(trun, obs, cens, sampling, seed[x])
            if (sum(kendallOnly, minp1Only, minp2Only) > 0) {
                tmpK <- tmpP1 <- tmpP2 <- NULL
                if (kendallOnly)
                    tmpK <- getKendall(perm.data$trun, perm.data$obs, perm.data$cens)[1]
                if (minp1Only) 
                    tmpP1 <- getMinP(perm.data$trun, perm.data$obs, perm.data$cens)
                if (minp2Only)                                     
                    tmpP2 <- getMinP(perm.data$trun, perm.data$obs, perm.data$cens, minp1 = FALSE)
                return(c(tmpK, tmpP1$maxTest, tmpP2$maxTest))
            }
            if (sum(kendallOnly, minp1Only, minp2Only) == 0) {
                kendallOnly <- minp1Only <- minp2Only <- TRUE
                tmpK <- getKendall(perm.data$trun, perm.data$obs, perm.data$cens)[1]
                tmpP1 <- getMinP(perm.data$trun, perm.data$obs, perm.data$cens)
                ## obsTest = ifelse(surv, 1e5, obsTest1))
                tmpP2 <- getMinP(perm.data$trun, perm.data$obs, perm.data$cens, minp1 = FALSE)  
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
    perm.data.Ken <- perm.data.P1 <- perm.data.P2 <- NULL
    if (!is.null(permKen) && all(!is.na(permKen))) 
        perm.data.Ken <- getPerm(trun, obs, cens, sampling, seed[which.min(abs(permKen))])
    if (!is.null(permTest1) && all(!is.na(permTest1))) {
        perm.data.P1 <- getPerm(trun, obs, cens, sampling, seed[which.min(permTest1)])
        permP1 <- 1 - pchisq(permTest1, df = 1)
    }
    if (!is.null(permTest2) && all(!is.na(permTest2))) {
        perm.data.P2 <- getPerm(trun, obs, cens, sampling, seed[which.min(permTest2)])
        permP2 <- 1 - pchisq(permTest2, df = 1)
    }
    out <- list(obsKen = obsKen, obsKenp = obsKenp, obsP1 = obsP1, obsP2 = obsP2,
                obsTest1 = obsTest1, obsTest2 = obsTest2,
                permKen = permKen, permP1 = permP1, permP2 = permP2, permTest1 = permTest1, permTest2 = permTest2,
                p.valueKen = ifelse(kendallOnly, (sum(abs(obsKen) < abs(permKen)) + 1)/(length(permKen) + 1), NA),
                p.valueMinp1 = ifelse(minp1Only, (sum(obsP1 > permP1) + 1) / (length(permKen) + 1), NA),
                p.valueMinp2 = ifelse(minp2Only, (sum(obsP2 > permP2) + 1) / (length(permKen) + 1), NA),
                perm.data.Ken = perm.data.Ken, perm.data.P1 = perm.data.P1, perm.data.P2 = perm.data.P2,
                kendallOnly = kendallOnly, minp1Only = minp1Only, minp2Only = minp2Only,
                sampling = sampling, permSize = permSize, random.seed = seed)
    class(out) <- "permDep"
    return(out)
}

## gpdPv <- function(permT, obsT) {
##     permT <- permT[order(permT, decreasing = TRUE)]
##     ## default starting value is 250
##     for (i in 0:24) {
##         tmp <- permT[1:(250 - 10 * i)]
##         ## Use bootstrap gof test, this depends on package goft
##         ## or we can use ks.test but this require to estimate the mle first
##         gofp <- gpd.test(tmp, 1000)$boot.test$p.value
##         if (gofp > 0.05) {
##             init <- c(sqrt(6 * var(tmp)) / pi, 0.1)
##             x <- optim(init, gpdlik, dat = tmp - mean(permT[250 - 10 * i], permT[251 - 10 * i]),
##                        hessian = TRUE, method = "Nelder-Mead")
##             break
##         }
##         ## if(i == 25) {
##         ##     tmp <- permT[1:300]
##         ##     init <- c(sqrt(6 * var(tmp)) / pi, 0.1)
##         ##     x <- optim(init, gpdlik, dat = tmp - mean(permT[300], permT[301]),
##         ##                hessian = TRUE, method = "Nelder-Mead")
##         ## }
##     }
##     (250 - 10 * i) * (1 - pgpd(obsT, -x$par[1], x$par[2])) / length(permT)
## }
