#' @export
print.permDep <- function(x, ...) {
    if (!is.permDep(x)) stop("Must be a permDep object")
    if (all(x$kendallOnly, x$minp1Only, x$minp2Only) || all(!x$kendallOnly, !x$minp1Only, !x$minp2Only)) {
        x$kendallOnly <- TRUE
        x$minp1Only <- TRUE
        x$minp2Only <- TRUE
    }
    if (x$kendallOnly & x$sampling == "conditional") {
        cat("\n Hypothesis test of quasi-independence based on conditional Kendall's tau statistics\n")
        cat(paste(" Conditional permutation with size = ", x$permSize, "\n"))
        cat(paste(" p-value =", sprintf("%.4f", x$p.valueKen), "\n"))
    }
    if (x$kendallOnly & x$sampling == "unconditional") {
        cat("\n Hypothesis test of quasi-independence based on unconditional Kendall's tau statistics\n")
        cat(paste(" Unconditional permutation size = ", x$permSize, "\n"))
        cat(paste(" p-value =", sprintf("%.4f", x$p.valueKen), "\n"))
    }
    if (x$minp1Only & x$sampling == "conditional") {
        cat("\n Hypothesis test of quasi-independence based on minp1 statistics\n")
        cat(paste(" Conditional permutation size = ", x$permSize, "\n"))
        cat(paste(" p-value =", sprintf("%.4f", x$p.valueMinp1), "\n"))
    }
    if (x$minp1Only & x$sampling == "unconditional") {
        cat("\n Hypothesis test of quasi-independence based on minp1 statistics\n")
        cat(paste(" Unconditional permutation size = ", x$permSize, "\n"))
        cat(paste(" p-value =", sprintf("%.4f", x$p.valueMinp1), "\n"))
    }
    if (x$minp2Only & x$sampling == "conditional") {
        cat("\n Hypothesis test of quasi-independence based on minp2 statistics\n")
        cat(paste(" Conditional permutation size = ", x$permSize, "\n"))
        cat(paste(" p-value =", sprintf("%.4f", x$p.valueMinp2), "\n"))
    }
    if (x$minp2Only & x$sampling == "unconditional") {
        cat("\n Hypothesis test of quasi-independence based on minp2 statistics\n")
        cat(paste(" Unconditional permutation size = ", x$permSize, "\n"))
        cat(paste(" p-value =", sprintf("%.4f", x$p.valueMinp2), "\n"))
    }    
    cat("\n")
}

is.permDep <- function(x) inherits(x, "permDep")
