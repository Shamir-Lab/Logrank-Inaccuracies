SUBTYPES.DATA = list(
  list(name='aml', only.primary=F, is.rna.seq=T, is.mirna.seq=T, display.name='AML'),
  list(name='breast', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='BIC'),
  list(name='colon', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='COAD'),
  list(name='gbm', only.primary=T, is.rna.seq=F, is.mirna.seq=F, display.name='GBM'),
  list(name='kidney', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='KIRC'),
  list(name='liver', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='LIHC'),
  list(name='lung', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='LUSC'),
  list(name='melanoma', only.primary=F, is.rna.seq=T, is.mirna.seq=T, display.name='SKCM'),
  list(name='ovarian', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='OV'),
  list(name='sarcoma', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='SARC'))

ALGORITHM.NAMES = c('kmeans', 'spectral', 'lracluster', 'pins', 'snf', 'mkl', 
                     'mcca', 'nmf', 'iCluster')
ALGORITHM.DISPLAY.NAMES = as.list(c('K-means', 'Spectral', 'LRAcluster', 'PINS', 
                           'SNF', 'rMKL-LPP', 'MCCA', 'MultiNMF', 'iClusterBayes'))
names(ALGORITHM.DISPLAY.NAMES) = ALGORITHM.NAMES

get.clustering.results.dir.path <- function() {
  return('')
}

get.empirical.pvalue.path <- function(subtype, algorithm) {
  emp.results.dir = ''
  res.path = file.path(emp.results.dir, 
                       sprintf('%s_%s', subtype, algorithm))
}

get.sim.results.dir <- function() {
  return('')
}

get.dataset.dir.path <- function() {
  return('')
}

get.clustering.results.file.path = function(subtype, alg.name) {
  results.path = get.clustering.results.dir.path()
  clustering.file.name = paste(subtype, alg.name, 'all', sep='_')
  survival.file.name = paste(clustering.file.name, 'surv', sep='_')
  return(list(clustering=file.path(results.path, clustering.file.name),
              surv=file.path(results.path, survival.file.name)))
}

get.figure.dir = function() {
  results.path = get.clustering.results.dir.path()
  return(file.path(results.path, 'logrank_letter_figures'))
}

get.subtype.survival.path <- function(subtype) {
  datasets.path = get.dataset.dir.path()
  survival.file.path = file.path(datasets.path, subtype, 'survival')
  return(survival.file.path)
}

get.surv.data <- function(groups, subtype, survival.file.path) {
  if (missing(survival.file.path)) {
    survival.file.path = get.subtype.survival.path(subtype)
  }
  survival.data = read.table(survival.file.path, header = TRUE)
  patient.names = names(groups)
  patient.names.in.file = as.character(survival.data[, 1])
  patient.names.in.file = toupper(gsub('-', '\\.', substring(patient.names.in.file, 1, 12)))

  stopifnot(all(patient.names %in% patient.names.in.file))
  
  indices = match(patient.names, patient.names.in.file)
  ordered.survival.data = survival.data[indices,]
  ordered.survival.data["cluster"] <- groups
  ordered.survival.data$Survival[is.na(ordered.survival.data$Survival)] = 0
  ordered.survival.data$Death[is.na(ordered.survival.data$Death)] = 0
  return(ordered.survival.data)
  
}

check.survival.coin <- function(groups, subtype) {
  return(check.survival(groups, subtype, is.coin=T))
}

check.survival <- function(groups, subtype, survival.file.path, is.coin=F) {
  if (missing(survival.file.path)) {
    survival.file.path = get.subtype.survival.path(subtype)
  }
  survival.data = read.table(survival.file.path, header = TRUE)
  patient.names = names(groups)
  patient.names.in.file = as.character(survival.data[, 1])
  patient.names.in.file = toupper(gsub('-', '\\.', substring(patient.names.in.file, 1, 12)))

  stopifnot(all(patient.names %in% patient.names.in.file))
  
  indices = match(patient.names, patient.names.in.file)
  ordered.survival.data = survival.data[indices,]
  ordered.survival.data["cluster"] <- groups
  ordered.survival.data$Survival[is.na(ordered.survival.data$Survival)] = 0
  ordered.survival.data$Death[is.na(ordered.survival.data$Death)] = 0
  if (is.coin) {
    ordered.survival.data$cluster <- as.factor(ordered.survival.data$cluster)
    return(coin::logrank_test(Surv(Survival, Death) ~ cluster, data=ordered.survival.data, distribution='asymptotic'))
  } else {
    return(survdiff(Surv(Survival, Death) ~ cluster, data=ordered.survival.data))
  }
  
}

get.cond.perm.surv <- function(clustering, subtype) {
  surv.data = get.surv.data(clustering, subtype)
  rownames(surv.data) = surv.data[,1]
  surv.data = surv.data[,-1]
  get.cond.perm.surv.given.data(surv.data)
}


################### Simulation code #####################

run.simulations.with.rate <- function(survival.rates, sim.name) {
  print('equal followup')
  
  print(run.simulations.spec.params(c(5, 5), c(0, 0), survival.rates, sim.ind1=sim.name, sim.ind2=1))
  print(run.simulations.spec.params(c(100, 100), c(0.04, 0.04), survival.rates, sim.ind1=sim.name, sim.ind2=2))

  print('unequal followup')
  print(run.simulations.spec.params(c(5, 5), c(0, 0.04), survival.rates, sim.ind1=sim.name, sim.ind2=3))
  print(run.simulations.spec.params(c(100, 100), c(0, 0.04), survival.rates, sim.ind1=sim.name, sim.ind2=4))
}

run.multi.group.simulations.with.rate <- function(survival.rates, sim.name) {
  print('equal followup')
  
  print(run.simulations.spec.params(c(100, 100, 100, 100), c(0.04, 0.04, 0.04, 0.04), survival.rates, sim.ind1=sim.name, sim.ind2=1))
  print(run.simulations.spec.params(c(20, 70, 130, 180), c(0.04, 0.04, 0.04, 0.04), survival.rates, sim.ind1=sim.name, sim.ind2=2))
  
  print('unequal followup')
  print(run.simulations.spec.params(c(100, 100, 100, 100), c(0, 0, 0.04, 0.04), survival.rates, sim.ind1=sim.name, sim.ind2=3))
  
  print(run.simulations.spec.params(c(20, 70, 130, 180), c(0.08, 0.04, 0.04, 0), survival.rates, sim.ind1=sim.name, sim.ind2=4))
  print(run.simulations.spec.params(c(20, 70, 130, 180), c(0, 0.04, 0.04, 0.08), survival.rates, sim.ind1=sim.name, sim.ind2=5))
}

run.simulations <- function() {
  run.simulations.with.rate(c(0.04, 0.04), 1)
  run.simulations.with.rate(c(0.04, 0.08), 2)
  run.simulations.with.rate(c(0.08, 0.04), 3)

  run.multi.group.simulations.with.rate(c(0.04, 0.04, 0.04, 0.04), 4)
  run.multi.group.simulations.with.rate(c(0.04, 0.04, 0.04, 0.08), 5)
  run.multi.group.simulations.with.rate(c(0.08, 0.04, 0.04, 0.04), 6)

  clusternomics.cluster.sizes = CLUSTERNOMICS.CLUSTER.SIZES
  num.clusters = length(clusternomics.cluster.sizes)
  print(run.simulations.spec.params(clusternomics.cluster.sizes, runif(num.clusters, 0, 0.04), rep(0.04, num.clusters), sim.ind1=7, sim.ind2=1))
  print(run.simulations.spec.params(clusternomics.cluster.sizes, runif(num.clusters, 0, 0.04), runif(num.clusters, 0.04, 0.08), sim.ind1=8, sim.ind2=1))
  
  mkl.lung.sizes = c(36, 40, 16, 43, 72, 134)
  num.clusters2 = length(mkl.lung.sizes)
  print(run.simulations.spec.params(mkl.lung.sizes, runif(num.clusters2, 0, 0.04), rep(0.04, num.clusters2), sim.ind1=7, sim.ind2=2))
  print(run.simulations.spec.params(mkl.lung.sizes, runif(num.clusters2, 0, 0.04), runif(num.clusters2, 0.04, 0.08), sim.ind1=8, sim.ind2=2))
}

# 1.3e4 permutations guarantee that for a binom variable with p=0.05, the 99% confidence interval will be < 0.01.
run.simulations.spec.params <- function(group.sizes, follow.up.rates, survival.rates, num.sims=1.3e4, sim.ind1=1, sim.ind2=1) {
        sim.dir = get.sim.results.dir()
        sim.path = file.path(sim.dir, paste0(sim.ind1, '_', sim.ind2))
        if (!file.exists(sim.path)) {
	  all.pvals = do.call(rbind, lapply(1:num.sims, function(i) {
                  set.seed(1e5 + i)
                print(paste('sim num', i))
	  	sim.data = simulate.data(group.sizes, follow.up.rates, survival.rates)
	  	cond.perm.pval = get.cond.perm.surv.given.data(sim.data, 1e3)$pvalue
                cond.asym.pval = get.cond.asym.surv.given.data(sim.data)

	  	return(c(cond.perm.pval, cond.asym.pval))
	  }))
          save(all.pvals, file=sim.path)
        }
        load(sim.path) 
        sig.threshold = 0.05
        return(colSums(all.pvals < sig.threshold) / num.sims)
}

simulate.data <- function(group.sizes, follow.up.rates, survival.rates, accrual.follow.up.range=c(12, 60)) {
  ngroups = length(group.sizes)
  sim.data = do.call(rbind, lapply(1:ngroups, function(i) {
    size = group.sizes[i]
        if (follow.up.rates[i] == 0) {
	  follow.up.times = rep(Inf, size)
        } else {
	  follow.up.times = rexp(size, follow.up.rates[i])
        }
	accrual.follow.up.times = runif(size, accrual.follow.up.range[1], accrual.follow.up.range[2])
	survival.times = rexp(size, survival.rates[i])
	obs.time = pmin(survival.times, follow.up.times, accrual.follow.up.times)
	cens = ifelse(survival.times < follow.up.times & survival.times < accrual.follow.up.times, 1, 0)
	ret = data.frame(Survival=obs.time, Death=cens, cluster=i)
	return(ret)
  }))
  return(sim.data)
}




###### Real data runs #######


get.all.logranks = function() {

  all.approx.logrank = matrix(NA, nrow=length(SUBTYPES.DATA), ncol=length(ALGORITHM.NAMES))
  all.empirical.logrank = matrix(NA, nrow=length(SUBTYPES.DATA), ncol=length(ALGORITHM.NAMES))
  upper.conf = matrix(NA, nrow=length(SUBTYPES.DATA), ncol=length(ALGORITHM.NAMES))
  lower.conf = matrix(NA, nrow=length(SUBTYPES.DATA), ncol=length(ALGORITHM.NAMES))
  
  results.path = get.clustering.results.dir.path()
  for (i in 1:length(SUBTYPES.DATA)) {
    for (j in 1:length(ALGORITHM.NAMES)) {
      subtype.name = SUBTYPES.DATA[[i]]$name
      algorithm.name = ALGORITHM.NAMES[[j]]
      results.ret = get.clustering.results.file.path(subtype.name, algorithm.name)
      # load solution
      clustering.file.name = results.ret$clustering
      
      survival.file.name = results.ret$surv
      load(survival.file.name)

      load(clustering.file.name)
      cur.survdiff = check.survival(clustering, subtype.name)
      all.approx.logrank[i, j] = get.logrank.pvalue(cur.survdiff)

      emp.res.path = get.empirical.pvalue.path(subtype.name, algorithm.name)
      if (!file.exists(emp.res.path)) {
        ret = get.cond.perm.surv(clustering, subtype.name)
        save(ret, file=emp.res.path)
      }
      load(emp.res.path)
      all.empirical.logrank[i, j] = ret$pvalue
      lower.conf[i, j] = ret$conf.int[1]
      upper.conf[i, j] = ret$conf.int[2]
    }
  }
  return(list(appr=all.approx.logrank, empirical=all.empirical.logrank,
              lower_conf=lower.conf, upper_conf=upper.conf))
}


plot.fig.a = function(all.logranks) {
  all.long.logranks = lapply(all.logranks, melt)
  num.out.of.conf = sum(all.long.logranks$lower_conf$value > all.long.logranks$appr$value | all.long.logranks$upper_conf$value < all.long.logranks$appr$value)
  print(paste0('approximations outside their confidence intervals: ', num.out.of.conf))
  num.appr.more.sig = sum(all.long.logranks$lower_conf$value > all.long.logranks$appr$value)
  print(paste0('approximations that are falsely more significant: ', num.appr.more.sig))
  
  num.significant = sum(all.long.logranks$appr$value <= 0.05)
  print(paste0('number of significant is ', num.significant))
  
  out.of.conf = sum((all.long.logranks$appr$value <= 0.05) & (all.long.logranks$lower_conf$value > all.long.logranks$appr$value | all.long.logranks$upper_conf$value < all.long.logranks$appr$value))
  num.out.of.conf = sum(out.of.conf)
  print(paste0('approximations outside their confidence intervals: ', num.out.of.conf))
  num.appr.more.sig = sum((all.long.logranks$appr$value <= 0.05) & (all.long.logranks$lower_conf$value > all.long.logranks$appr$value))
  print(paste0('approximations that are falsely more significant: ', num.appr.more.sig))
  
  falsely.called = all.long.logranks$appr$value <= 0.05 & all.long.logranks$empirical$value > 0.05 & out.of.conf
  num.falsely.called = sum(falsely.called)
  print(paste0('number of falsely called is ', num.falsely.called))
  print(all.long.logranks$appr[falsely.called,])
  print(all.long.logranks$empirical[falsely.called,])
  print(all.long.logranks$lower_conf[falsely.called,])
  print(all.long.logranks$upper_conf[falsely.called,])
  
  all.pvalues = c(-log10(all.long.logranks$appr$value), -log10(all.long.logranks$empirical$value))
  max.value = max(all.pvalues[is.finite(all.pvalues)])
  png(file.path(get.figure.dir(), 'fig_a.png'), width=1800, height=1800, res=300)
  dif = -log10(all.long.logranks$appr$value) + log10(all.long.logranks$empirical$value)
  print(paste0('number of red dots ', sum(dif > log10(2))))
  plot(-log10(all.long.logranks$appr$value), -log10(all.long.logranks$empirical$value), col=ifelse(dif > log10(2), 'red', 'black'),
       ylab='Empirical -log10(p-value)', xlab='Approximate -log10(p-value)', ylim=c(0, max.value), xlim=c(0, max.value), pch=19, cex.lab=1.3)
  abline(a=0, b=1)
  segments(-log10(all.long.logranks$appr$value), -log10(all.long.logranks$upper_conf$value), -log10(all.long.logranks$appr$value), -log10(all.long.logranks$lower_conf$value), col=ifelse(dif > log10(2), 'red', 'black'))
  dev.off()
}


get.empirical.pvalue.distrib = function(clustering, subtype, num.perms=1e6) {
  surv.data = get.surv.data(clustering, subtype)
  rownames(surv.data) = surv.data[,1]
  surv.data = surv.data[,-1]
  nsamples = nrow(surv.data)
  imp = imputeHeinze(surv.data)
  perm.pvals = as.numeric(mclapply(1:num.perms, function(i) {
      if (i %% 1e4 == 1) print(i)
      perm = sample(1:nsamples, nsamples)
      perm.data = permuteHeinze(imp, perm)
      surv.ret = survdiff(Surv(Survival, Death) ~ cluster, data=perm.data)
      cur.pvalue = get.logrank.pvalue(surv.ret)
      return(cur.pvalue)
    }, mc.cores=50))
  return(perm.pvals)
}

plot.kaplan.meier.curve <- function(surv.data, plot.path) {
  surv.ret = survfit(Surv(Survival, Death) ~ cluster, data=surv.data)
  png(plot.path, width=1800, height=1800, res=300)
  plot(surv.ret, ylab='Survival', xlab='Days', col=rainbow(length(table(surv.data$cluster))))
  dev.off()
}

plot.fig.b = function(subtype, alg.name) {

  subtype.file.paths = get.clustering.results.file.path(subtype, alg.name)
  surv.file.path = subtype.file.paths$surv
  cluster.file.path = subtype.file.paths$clustering
  load(cluster.file.path)
  load(surv.file.path)
  orig.pval = get.logrank.pvalue(check.survival(clustering, subtype))
  num.perms = 1e6
  
  surv.data = get.surv.data(clustering, subtype)
  plot.kaplan.meier.curve(surv.data, file.path(get.figure.dir(), 'mcca_km.png'))
  
  results.path = get.clustering.results.dir.path()
  distrib.path = file.path(results.path, 'perm_pvals')
  if (!file.exists(distrib.path)) {
    perm.pvalues = get.empirical.pvalue.distrib(clustering, subtype, num.perms)
    save(perm.pvalues, file=distrib.path)
  } else {
    load(distrib.path)
  }
  
  # plot the distribution of approx p.values  
  #plot(density(perm.pvalues))
  print(sum(perm.pvalues <= 0.05))
  print(sum(perm.pvalues <= 0.05) / num.perms)
  
  png(file.path(get.figure.dir(), 'fig_b.png'), width=1800, height=1800, res=300)
  histo<-hist(perm.pvalues)
  barplot(histo$counts/length(perm.pvalues), yaxt='n', )->bp 
  box()
  axis(2, at=seq(0, 1, by=0.01))
  axis(1,at=c(bp,24.6)-0.55,labels=histo$breaks[1:(length(histo$breaks))])
  abline(a=0.05, b=0, col='red')
  title(ylab="Relative Frequency",xlab="Approximate p-value", cex.lab=1.3)
  dev.off()
  
}

calc.empirical.significance.prob = function(clustering, subtype, orig.pvalue=0.05) {
  should.continue = T
  total.num.perms = 0
  total.num.sig.pvalue = 0
  
  while (should.continue) {
    num.perms = 1e5
    perm.pvalues = as.numeric(mclapply(1:num.perms, function(i) {
      if (i %% 1e5 == 1) print(i)
      cur.clustering = sample(clustering)
      names(cur.clustering) = names(clustering)
      cur.pvalue = get.logrank.pvalue(check.survival(cur.clustering, subtype))
      return(cur.pvalue)
    }, mc.cores=50))
    
    total.num.perms = total.num.perms + num.perms
    total.num.sig.pvalue = total.num.sig.pvalue + sum(perm.pvalues <= orig.pvalue)
    
    binom.ret = binom.test(total.num.sig.pvalue, total.num.perms)
    cur.pvalue = binom.ret$estimate
    cur.conf.int = binom.ret$conf.int
  
    sig.threshold = 0.05
    is.conf.small = ((cur.conf.int[2] - cur.pvalue) < min(cur.pvalue / 10, 0.01)) & ((cur.pvalue - cur.conf.int[1]) < min(cur.pvalue / 10, 0.01))
    is.threshold.in.conf = cur.conf.int[1] < sig.threshold & cur.conf.int[2] > sig.threshold
    if ((is.conf.small & !is.threshold.in.conf) | total.num.perms > 5e6) {
      should.continue = F
    }
  }
  
  return(binom.ret)
}


plot.fig.c = function(subtype) {
  
  load(get.clustering.results.file.path(subtype, 'mcca')$clustering)
  patient.names = names(clustering)
  diff.num.clusters.dir = file.path(get.clustering.results.dir.path(), 'different_num_clusters_exact')
  dir.create(diff.num.clusters.dir)
  
  max.num.clusters = 20
  all.binom.tests = as.list(rep(NA, max.num.clusters))
  for (num.clusters in seq(max.num.clusters, 2, by=-1)) {
    print(paste0('num.clusters is now ', num.clusters))
    #initial.solution = rep_len(1:num.clusters, length(clustering))
    initial.solution = do.call(c, lapply(1:(num.clusters - 1), function(x) rep(x, 10)))
    initial.solution = c(initial.solution, rep(num.clusters, length(clustering) - length(initial.solution)))
    names(initial.solution) = names(clustering)
    cur.file.path = file.path(diff.num.clusters.dir, num.clusters)
    if (!file.exists(cur.file.path)) {
      binom.test.ret = calc.empirical.significance.prob(initial.solution, subtype)
      save(binom.test.ret, file=cur.file.path)
    } else {
      load(cur.file.path)
    }
    
    all.binom.tests[[num.clusters]] = binom.test.ret
    print(binom.test.ret)
  }
  
  png(file.path(get.figure.dir(), 'fig_c.png'), width=1800, height=1800, res=300)
  estimates = sapply(all.binom.tests[2:max.num.clusters], function(x) x$estimate)
  max.deviation = max(sapply(all.binom.tests[2:max.num.clusters], function(x) x$conf.int[2] - x$conf.int[1]))
  print(paste0('max deviation is ', max.deviation))
  plot(2:max.num.clusters, estimates, ylim=c(0.04, 0.1*ceiling(max(estimates)/0.1)), ylab='Prob p-value <= 0.05', xlab='Number of clusters', pch=19, cex.lab=1.3)
  abline(a=0.05, b=0, col='red')
  dev.off()
}


create.letter.results = function() {
  par(mar=rep(0, 4))
  plot.fig.a(get.all.logranks())
  plot.fig.b('kidney', 'mcca')
  plot.fig.c('breast')
}


murine.analysis = function() {
  paper.groups = c(rep(1, 8), rep(2, 8))
  names(paper.groups) = paste0('MOUSE.', 1:16)
  surv.data = get.surv.data(paper.groups, 'mouse_test')
  surv.data$cluster = as.factor(surv.data$cluster)
  mid.ranks.pvalue = coin::pvalue(coin::logrank_test(Surv(Survival, Death) ~ cluster, data=surv.data, ties.method='mid-ranks', distribution='exact'))
  hothorn.pvalue = coin::pvalue(coin::logrank_test(Surv(Survival, Death) ~ cluster, data=surv.data, ties.method='Hothorn', distribution='exact'))
  ave.scores.pvalue = coin::pvalue(coin::logrank_test(Surv(Survival, Death) ~ cluster, data=surv.data, ties.method='average', distribution='exact'))
  asym.pvalue = get.logrank.pvalue(survdiff(Surv(Survival, Death) ~ cluster, data=surv.data))
  print(c(mid.ranks.pvalue, hothorn.pvalue, ave.scores.pvalue, asym.pvalue))
}

CLUSTERNOMICS.CLUSTER.SIZES = c(63, 58, 40, 31, 29, 25, 20, 19, 16, 10,  9,  9,  7,  5,  4,  3)

clusternomics.analysis = function() {
  
  subtype = 'breast_clusternomics'
  cluster.sizes = CLUSTERNOMICS.CLUSTER.SIZES
  survival.file.path = get.subtype.survival.path(subtype)
  ids = as.character(read.table(survival.file.path, header=T)[,1])

  clustering = c()
  for (i in seq_along(cluster.sizes)) {
    clustering = c(clustering, rep(i, cluster.sizes[i]))
  }
  names(clustering) = ids
  orig.pvalue = 0.0381
  
  set.seed(42)
  ret.cond = calc.empirical.significance.prob(clustering, subtype, orig.pvalue=orig.pvalue)
  print(ret.cond)
}
