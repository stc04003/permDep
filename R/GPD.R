rgpd <- function(n, xi, beta) {
    if (xi == 0) {
        rexp(n, beta)
    } else {
        beta * (1 - runif(n) ^ xi) / xi
    }
}

pgpd <- function(q, xi, beta) {
    if (xi == 0) {
        cdf <- 1 - exp(-q / beta)
    } else {
        cdf <- ifelse(q < 0, 0, 1 - (1 - xi * q / beta) ^ xi ^-1)
    }
    cdf <- ifelse(is.na(cdf) | cdf > 1, 1, cdf)
    ifelse(cdf < 0, 0, cdf)
}

dgpd <- function(x, xi, beta) {
    if (xi == 0) {
        ifelse(x <= 0, 0, exp(-x / beta) / beta)
    } else {
        ifelse(x <= 0 | x >= max(beta / xi, -xi * Inf), 0, (1 - xi * x / beta) ^ (xi^-1 - 1) / beta)
    }
}

gpdlik <- function(parm, dat) {
    xi <- parm[1]
    beta <- parm[2]
    y <- 1 + xi * dat / beta
    if (beta <= 0 || min(y) <= 0)
        1e6
    else {
        length(dat) * log(beta) + sum(log(y) * (1/xi + 1))
    }
}

