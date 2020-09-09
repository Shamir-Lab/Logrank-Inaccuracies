# This is an R implementation of the Heinze conditional log-rank test for at least two groups of samples.
# The Heinze algorithm is not assymptotic, unlike the log-rank implementation from the survival package.
# The assymptotic test tends to overestimate the difference in survival between groups. The Heinze
# method should therefore be preferred when the sample size is small, the group sizes are unbalanced,
# or the number of events in each group is small. For more details on the inaccuracy of the log-rank test
# even in large data sizes see:
# https://doi.org/10.15252/msb.20188754
# The algorithm is based on permuting the data and imputing event times. Full details on the algorithm
# appear in:
# 10.1111/j.0006-341x.2003.00132.x
# The implementation is an extension of the implementation from the R package permGS:
# https://github.com/cran/permGS
# The main difference is that this implementation supports more than two groups of samples.


#' @title Get p-value from survdiff
#' @name get.logrank.pvalue
#' @description Get log-rank's p-value from an object returned by survdiff.
#' @param survdiff.res return value from the function survidff from the library survival.
#' @return survdiff's log-rank p-value.
#' @export
get.logrank.pvalue <- function(survdiff.res) {
  1 - pchisq(survdiff.res$chisq, length(survdiff.res$n) - 1)
}

#' @title Asymptotic log-rank test (deprecated)
#' @name get.cond.asym.surv.given.data
#' @description Get the p-value for the asymptotic log-rank test. The test is performed using the survival
#' package survdiff function. This function is deprecated, and survdiff should be used for the asymptotic
#' log-rank test.
#' @param surv.data A data.frame containing three columns:
#' Survival - the event time.
#' Death - binary value. 1 is the event time is observed, 0 if censored.
#' cluster - variable indicating the group of each sample.
#' @return asymptotic log-rank p-value.
#' @examples
#' data('ColonSurvivalData')
#' get.cond.asym.surv.given.data(colon.survival.data)
#' @export
get.cond.asym.surv.given.data <- function(surv.data) {
  .Deprecated(new='', msg=paste0('This function is only maintained for backward compatibility. ',
                                 'Use survdiff to perform the asysmptotic log-rank test.'))
  .inner.get.cond.asym.surv.given.data(surv.data)
}

.inner.get.cond.asym.surv.given.data <- function(surv.data) {
  death.times = surv.data$Survival[surv.data$Death == 1]
  max.event.time = max(surv.data$Survival)
  if ((length(death.times) == 0) || (all(death.times == death.times[1]) & (death.times[1] == max.event.time))) {
    return(1)
  } else {
    surv.ret = survival::survdiff(survival::Surv(Survival, Death) ~ cluster, data=surv.data)
    return(get.logrank.pvalue(surv.ret))
  }
}


#' @title Heinze log-rank test
#' @name logrank.heinze
#' @description Get the p-value for the Heinze log-rank test.
#' @param formula A formula of the form time + event.type ~ group. time is the time of the event,
#' event.type is 1 if the event was observed and 0 if it was censored, and group is the group of each sample.
#' @param data A data.frame in which the formula is evaluated. If formula is NULL, the first column of data
#' is used as time, the second as event.type, and the third as group.
#' @param max.num.perms The maximal number of permutations to perform.
#' @param num.perms.between.continue.tests The number of permutations performed between tests of whether
#' more permutations should be performed. Defaults to 1000.
#' @param significance.threshold If this number is set, permutations continue as long as the confidence interval
#' for the p-value contain this number (or until max.num.perms permutations are performed).
#' @param confidence.relative.size If this number is set, permutations continue as long as the 95\%
#' confidence interval for the p-value is "relatively small". That is, its boundaries are within a
#' confidence.relative.size factor of the estimated p-value. For example, if the estimated p-value is 0.01
#' and confidence.relative.size=2, then the 95\% confidence interval must be contained in [0.005, 0.02].
#' @param should.continue.func A function used to decide whether more permutations should be performed.
#' If this is set, confidence.relative.size and should.continue.func cannot be set. The function receives three
#' arguments: the number of permutations performed so far, the number of permutations with more extreme
#' log-rank statistic than the original, and the original p-value.
#' The function should return a binary value of whether more permutations should be performed.
#' @param mc.cores Number of cores to use to parallelize permutations.
#' @param seed Initialization of R's random seed, or NULL to not set the seed.
#' @param verbose True <=> info on the number of permutations should be printed.
#' @return A list with the following entries:
#' pvalue - the estimation for log-rank's p-value.
#' conf.int - 95\% confidence interval for log-rank's p-value.
#' total.num.perms - number of permutations performed to estimate the p-value.
#' total.num.extreme.pvals - the number of permutations for which the log-rank statistic is (equally or) more
#' extreme than the log-rank stitistic for the data.
#' @examples
#' data('ColonSurvivalData')
#' logrank.heinze(data=colon.survival.data, max.num.perms=1000)
#' logrank.heinze(Survival + Death ~ cluster, data=colon.survival.data, max.num.perms=1000)
#' @export
logrank.heinze <- function(formula=NULL, data, max.num.perms=NULL,
                           num.perms.between.continue.tests=1000, significance.threshold=NULL,
                           confidence.relative.size=NULL, should.continue.func=NULL,
                           mc.cores=1, seed=NULL, verbose=TRUE) {

  if (all(is.null(max.num.perms) & is.null(significance.threshold) &
          is.null(confidence.relative.size) & is.null(should.continue.func))) {
    stop('No stopping criterion was set!')
  }

  if (!is.null(should.continue.func)) {
    if (!is.null(significance.threshold) | !is.null(confidence.relative.size)) {
      stop('Cannot set significance.threshold or confidence.relative.size if should.continue.func is supplied')
    }
  }

  stopifnot(mc.cores >= 1)
  if (!is.null(max.num.perms)) {
    stopifnot(max.num.perms >= 1)
  }
  if (!is.null(num.perms.between.continue.tests)) {
    stopifnot(num.perms.between.continue.tests >= 1)
  }
  if (!is.null(confidence.relative.size)) {
    stopifnot(confidence.relative.size >= 1)
  }

  if (is.null(formula)) {
    event.times = data[,1]
    is.censored = data[,2]
    groups = data[,3]
  } else {
    rhs = formula.tools::rhs.vars(formula)
    lhs = formula.tools::lhs.vars(formula)
    if (length(rhs) != 1 | length(lhs) != 2) {
      stop('Illegal formula format. Must be of the form: time + event.type ~ group.')
    }
    group.var = rhs[1]
    time.var = lhs[1]
    event.var = lhs[2]
    if (any(!(c(group.var, time.var, event.var) %in% colnames(data)))) {
      stop('Formula contains expressions which are not found in the supplied data')
    }
    event.times = data[,time.var]
    is.censored = data[,event.var]
    groups = data[,group.var]
  }

  surv.data = data.frame(Survival=event.times, Death=is.censored, cluster=groups)
  orig.pvalue = .inner.get.cond.asym.surv.given.data(surv.data)
  if (verbose) {
    print(sprintf('Asymptotic log-rank p-value: %s', orig.pvalue))
  }

  if (is.null(num.perms.between.continue.tests)) {
    num.perms.between.continue.tests = 1000
  }

  if (is.null(max.num.perms)) {
    num.remaining.perms = Inf
  } else {
    num.remaining.perms = max.num.perms
  }


  should.continue = TRUE
  total.num.perms = 0
  total.num.extreme.pvals = 0
  imp = imputeHeinze(surv.data)
  nsamples = nrow(surv.data)
  if (mc.cores == 1) {
    perm.func = lapply
    should.set.seed = FALSE
  } else {
    perm.func = function(...) {parallel::mclapply(..., mc.cores=mc.cores)}
    RNGkind("L'Ecuyer-CMRG")
    should.set.seed = TRUE
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  while (should.continue) {
    cur.num.perms = min(num.perms.between.continue.tests, num.remaining.perms)
    if (verbose) {
      print(sprintf('Will now perform %s permutations', cur.num.perms))
    }

    perm.pvals = as.numeric(perm.func(1:cur.num.perms, function(i) {
      perm = sample(1:nsamples, nsamples)
      perm.data = permuteHeinze(imp, perm)
      cur.pval = .inner.get.cond.asym.surv.given.data(perm.data)
      return(cur.pval)
    }))

    # advance random seed
    if (should.set.seed) {
      x <- .Random.seed
      for (i in seq_len(cur.num.perms)) {
        x = nextRNGStream(x)
        assign('.Random.seed', x, pos=.GlobalEnv)
      }
    }

    total.num.perms = total.num.perms + cur.num.perms
    total.num.extreme.pvals = total.num.extreme.pvals + sum(perm.pvals <= orig.pvalue)

    binom.ret = binom.test(total.num.extreme.pvals, total.num.perms)
    cur.pvalue = binom.ret$estimate
    cur.conf.int = binom.ret$conf.int

    if (verbose) {
      print(sprintf('Performed %s permutations so far, %s with equally or more statistic than original', total.num.perms, total.num.extreme.pvals))
      print(sprintf('Current p-value estimate is %s', cur.pvalue))
      print(paste0('Current 95% confidence interval for the p-value is [', cur.conf.int[1], ', ', cur.conf.int[2], ']'))
    }

    # check if we should continue
    num.remaining.perms = num.remaining.perms - cur.num.perms
    if (num.remaining.perms <= 0) {
      should.continue = FALSE
    } else if (!is.null(should.continue.func)) {
      should.continue = should.continue.func(total.num.perms, total.num.extreme.pvals, orig.pvalue)
    } else if (!is.null(significance.threshold) | !is.null(confidence.relative.size)) {
      is.threshold.in.interval = FALSE
      is.interval.small = TRUE
      if (!is.null(significance.threshold)) {
        is.threshold.in.interval = cur.conf.int[1] < significance.threshold & cur.conf.int[2] > significance.threshold
      }
      if (!is.null(confidence.relative.size)) {
        is.interval.small = (cur.pvalue / cur.conf.int[1] <= confidence.relative.size) &
          (cur.conf.int[2] / cur.pvalue <= confidence.relative.size)
      }
      should.continue = !is.interval.small | is.threshold.in.interval
    }
  }

  return(list(pvalue = unname(cur.pvalue), conf.int = cur.conf.int, total.num.perms=total.num.perms,
              total.num.extreme.pvals=total.num.extreme.pvals))
}

#' @title Heinze log-rank test (deprecated)
#' @name get.cond.perm.surv.given.data
#' @description Get the p-value for the Heinze log-rank test.
#' Permutations are performed until the 95\% confidence interval for log-rank's p-value does not contain
#' 0.05, and the size of the confidence interval is within 10\% of the estimated p-value.
#' @param surv.data A data.frame containing three columns:
#' Survival - the event time.
#' Death - binary value. 1 is the event time is observed, 0 if censored.
#' cluster - variable indicating the group of each sample.
#' @param max.num.perms If this number is set, max.num.perms are performed.
#' @param verbose True <=> info on the number of permutations should be printed.
#' @param mc.cores Number of cores to use to parallelize permutations.
#' @param seed Initialization of R's random seed, or NULL to not set the seed.
#' @return A list with the following entries:
#' pvalue - the estimation for log-rank's p-value.
#' conf.int - 95\% confidence interval for log-rank's p-value.
#' total.num.perms - number of permutations performed to estimate the p-value.
#' total.num.extreme.pvals - the number of permutations for which the log-rank statistic is (equally or) more
#' extreme than the log-rank stitistic for the data.
#' @examples
#' data('ColonSurvivalData')
#' get.cond.perm.surv.given.data(colon.survival.data)
#' @export
get.cond.perm.surv.given.data <- function(surv.data, max.num.perms=NULL, verbose=T, mc.cores=50, seed=42) {
  .Deprecated(new='', msg=paste0('This function is only maintained for backward compatibility. ',
                                 'Use logrank.heinze to perform the log-rank test.'))
  if (!is.null(seed)) {
    set.seed(seed)
  }

  orig.pvalue = .inner.get.cond.asym.surv.given.data(surv.data)

  if (is.null(max.num.perms)) {
    num.perms = round(min(max(10 / orig.pvalue, 1000), 1e5))
  } else {
    num.perms = max.num.perms
  }
  should.continue = T

  total.num.perms = 1
  total.num.extreme.pvals = 1

  imp = imputeHeinze(surv.data)
  nsamples = nrow(surv.data)

  while (should.continue) {
    if (verbose) {
      print('Another iteration in empirical survival calculation')
      print(num.perms)
    }
    perm.pvals = as.numeric(parallel::mclapply(1:num.perms, function(i) {
      perm = sample(1:nsamples, nsamples)
      perm.data = permuteHeinze(imp, perm)
      cur.pval = .inner.get.cond.asym.surv.given.data(perm.data)
      return(cur.pval)
    }, mc.cores=mc.cores))

    total.num.perms = total.num.perms + num.perms
    total.num.extreme.pvals = total.num.extreme.pvals + sum(perm.pvals <= orig.pvalue)

    binom.ret = binom.test(total.num.extreme.pvals, total.num.perms)
    cur.pvalue = binom.ret$estimate
    cur.conf.int = binom.ret$conf.int

    if (verbose) {
      print(c(total.num.extreme.pvals, total.num.perms))
      print(cur.pvalue)
      print(cur.conf.int)
    }

    sig.threshold = 0.05
    is.conf.small = ((cur.conf.int[2] - cur.pvalue) < min(cur.pvalue / 10, 0.01)) & ((cur.pvalue - cur.conf.int[1]) < min(cur.pvalue / 10, 0.01))
    is.threshold.in.conf = cur.conf.int[1] < sig.threshold & cur.conf.int[2] > sig.threshold
    if ((!is.null(max.num.perms)) | (is.conf.small & !is.threshold.in.conf) | (total.num.perms > 1e7)) {
      should.continue = F
    } else {
      num.perms = 1e4
    }
  }

  return(list(pvalue = cur.pvalue, conf.int = cur.conf.int, total.num.perms=total.num.perms,
              total.num.extreme.pvals=total.num.extreme.pvals))
}

# Code for the exact log-rank test conditioned on the observed follow-up from Heinze et. al,
# supporting 2 or more groups.
# Based on code from the R permGS package:
# https://github.com/cran/permGS

sampleFromKM <- function(n, fit, start=0, tmax=NULL, dv=1) {
  if(is.null(tmax)) tmax <- fit$time[length(fit$time)]
  vapply(runif(n, start, 1), function(v) {
    if(v > max(1-fit$surv)) c(tmax, 0)
    else c(min(fit$time[fit$surv <= 1-v]), dv)
  }, c(NA_real_, NA_integer_))
}

sampleFromCondKM <- function(U, fit, tmax=NULL, dv=1, f=NULL) {
  n <- length(U)
  if(is.null(f)) f <- approxfun(fit$time, fit$surv, method="constant", yleft=1, rule=2, f=0)
  sampleFromKM(n, fit, 1-f(U), tmax, dv)
}

imputeHeinze <- function(data) {
  time <- data[,1]
  status <- data[,2]

  groups <- data[,3]
  stopifnot(!is.factor(groups))
  ngroups = length(table(groups))
  stopifnot(all(sort(as.numeric(names(table(groups)))) == 1:ngroups))

  tmax <- max(time)

  surv.fit <- survival::survfit(survival::Surv(time, status) ~ 1)
  cen.fit.list = lapply(1:ngroups, function(g) {
    g.ind = data[,3] == g
    gdata = data[g.ind,]
    cur.time <- gdata[,1]
    cur.status <- gdata[,2]

    cen.fit <- survival::survfit(survival::Surv(cur.time, 1-cur.status) ~ 1)
    return(cen.fit)
  })

  surv.fit.inter <- approxfun(surv.fit$time, surv.fit$surv, method="constant", yleft=1, rule=2, f=0)
  cen.fit.inter.list <- lapply(cen.fit.list, function(x) approxfun(x$time, x$surv, method="constant", yleft=1, rule=2, f=0))
  ret = list(surv.fit=surv.fit, cen.fit.list=cen.fit.list,
             surv.fit.inter=surv.fit.inter, cen.fit.inter.list=cen.fit.inter.list, tmax=tmax, data=data)
  return(ret)
}

permuteHeinze <- function(imp, pp) {
  ## permute rows
  pdata <- imp$data[pp, ]

  t.orig <- imp$data[,1]
  delta.orig <- imp$data[,2]
  epsilon.orig <- 1 - delta.orig
  f.orig <- t.orig

  t.star <- pdata[,1]
  delta.star <- as.logical(pdata[,2])

  T = t.star
  C = t.orig

  # impute survival times
  groups = imp$data[,3]
  ngroups = length(table(groups))
  for (g in 1:ngroups) {
    group.cens <- !delta.star & groups == g
    if(any(group.cens)) {
      tmp <- sampleFromCondKM(T[group.cens], imp$surv.fit, imp$tmax, 1, imp$surv.fit.inter)
      T[group.cens] <- tmp[1,]
      delta.star[group.cens] <- tmp[2,]
    }
  }

  # impute censoring times
  for (g in 1:ngroups) {
    group.obs <- delta.orig & groups == g
    if(any(group.obs)) {
      C[group.obs] <- sampleFromCondKM(C[group.obs], imp$cen.fit.list[[g]], imp$tmax, 0, imp$cen.fit.inter.list[[g]])[1,]
    }
  }

  pY <- pmin(T, C)
  delta.star <- (T <= C) * delta.star

  ret = matrix(c(pY, delta.star, imp$data[,3]), ncol=3, nrow=length(pY))
  ret = as.data.frame(ret)
  colnames(ret) = c('Survival', 'Death', 'cluster')
  return(ret)
}


