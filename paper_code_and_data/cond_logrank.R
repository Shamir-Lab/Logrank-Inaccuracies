get.logrank.pvalue <- function(survdiff.res) {
  1 - pchisq(survdiff.res$chisq, length(survdiff.res$n) - 1)  
}

get.cond.asym.surv.given.data <- function(surv.data) {
  death.times = surv.data$Survival[surv.data$Death == 1]
  max.event.time = max(surv.data$Survival)
  if ((length(death.times) == 0) || (all(death.times == death.times[1]) & (death.times[1] == max.event.time))) {
    return(1)
  } else {
    surv.ret = survdiff(Surv(Survival, Death) ~ cluster, data=surv.data)
    return(get.logrank.pvalue(surv.ret))
  }
}

get.cond.perm.surv.given.data <- function(surv.data, max.num.perms=NULL, verbose=T) {
  set.seed(42)

  orig.pvalue = get.cond.asym.surv.given.data(surv.data)

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
    perm.pvals = as.numeric(mclapply(1:num.perms, function(i) {
        perm = sample(1:nsamples, nsamples)
        perm.data = permuteHeinze(imp, perm)
        cur.pval = get.cond.asym.surv.given.data(perm.data)
      return(cur.pval)
    }, mc.cores=50))
    
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

# Code for the exact log-rank test conditioned on the observed follow-up from Heinze et. al
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

    surv.fit <- survfit(Surv(time, status) ~ 1)
    cen.fit.list = lapply(1:ngroups, function(g) {
        g.ind = data[,3] == g
        gdata = data[g.ind,]
        cur.time <- gdata[,1]
        cur.status <- gdata[,2]

        cen.fit <- survfit(Surv(cur.time, 1-cur.status) ~ 1)
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

