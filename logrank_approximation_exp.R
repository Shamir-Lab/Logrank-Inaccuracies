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


get.all.logranks = function() {

  all.approx.logrank = matrix(NA, nrow=length(SUBTYPES.DATA), ncol=length(ALGORITHM.NAMES))
  all.empirical.logrank = matrix(NA, nrow=length(SUBTYPES.DATA), ncol=length(ALGORITHM.NAMES))
  upper.conf = matrix(NA, nrow=length(SUBTYPES.DATA), ncol=length(ALGORITHM.NAMES))
  lower.conf = matrix(NA, nrow=length(SUBTYPES.DATA), ncol=length(ALGORITHM.NAMES))
  
  results.path = get.clustering.results.dir.path()
  for (i in 1:length(SUBTYPES.DATA)) {
    for (j in 1:length(ALGORITHM.NAMES)) {
      results.ret = get.clustering.results.file.path(SUBTYPES.DATA[[i]]$name, ALGORITHM.NAMES[[j]])
      # load solution
      clustering.file.name = results.ret$clustering
      
      survival.file.name = results.ret$surv
      load(survival.file.name)
      all.empirical.logrank[i, j] = empirical.surv.ret$pvalue
      lower.conf[i, j] = empirical.surv.ret$conf.int[1]
      upper.conf[i, j] = empirical.surv.ret$conf.int[2]

      load(clustering.file.name)
      cur.survdiff = check.survival(clustering, SUBTYPES.DATA[[i]]$name)
      all.approx.logrank[i, j] = get.logrank.pvalue(cur.survdiff)
      
      #if (all.approx.logrank[i, j] == 0) {
      #  all.approx.logrank[i, j] = 1e-10
      #}
    }
  }
  return(list(appr=all.approx.logrank, empirical=all.empirical.logrank,  lower_conf=lower.conf, upper_conf=upper.conf))
}


plot.fig.a = function(all.logranks) {
  all.long.logranks = lapply(all.logranks, melt)
  # TODO: probably move to a different function
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
  perm.chisq = as.numeric(mclapply(1:num.perms, function(i) {
      if (i %% 1e4 == 1) print(i)
      cur.clustering = sample(clustering)
      names(cur.clustering) = names(clustering)
      cur.chisq = check.survival(cur.clustering, subtype)$chisq
      return(cur.chisq)
    }, mc.cores=50))
  return(perm.chisq)
}


plot.fig.b = function(subtype, alg.name) {

  subtype.file.paths = get.clustering.results.file.path(subtype, alg.name)
  surv.file.path = subtype.file.paths$surv
  cluster.file.path = subtype.file.paths$clustering
  load(cluster.file.path)
  load(surv.file.path)
  orig.chisq = check.survival(clustering, subtype)$chisq
  num.perms = 1e6
  
  results.path = get.clustering.results.dir.path()
  distrib.path = file.path(results.path, 'perm_chisq_values')
  if (!file.exists(distrib.path)) {
    perm.chisq = get.empirical.pvalue.distrib(clustering, subtype, num.perms)
    save(perm.chisq, file=distrib.path)
  } else {
    load(distrib.path)
  }
  
  
  perm.pvalues = (sapply(perm.chisq, function(chisq) {
    1 - pchisq(chisq, max(clustering) - 1)
  }))
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


calc.empirical.significance.prob = function(clustering, subtype) {
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
    total.num.sig.pvalue = total.num.sig.pvalue + sum(perm.pvalues <= 0.05)
    
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
  diff.num.clusters.dir = file.path(get.clustering.results.dir.path(), 'different_num_clusters')
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
  print(get.empirical.surv(paper.groups, 'mouse_analysis'))
  
}


clusternomics.analysis = function() {
  
  cluster.sizes = c(63, 58, 40, 31, 29, 25, 20, 19, 16, 10,  9,  9,  7,  5,  4,  3)

  clustering = c()
  for (i in seq_along(cluster.sizes)) {
    clustering = c(clustering, rep(i, cluster.sizes[i]))
  }
  names(clustering) = ids
  orig.pvalue = 0.0381
  subtype = 'breast_clusternomics'
    
  
  set.seed(42)
  # The initial number of permutations to run
  num.perms = round(min(max(10 / orig.pvalue, 1000), 1e6))
  should.continue = T
  
  total.num.perms = 0
  total.num.extreme.perms = 0
  
  while (should.continue) {
    print('Another iteration in empirical survival calculation')
    print(num.perms)
    perm.pval = as.numeric(mclapply(1:num.perms, function(i) {
      cur.clustering = sample(clustering)
      names(cur.clustering) = names(clustering)
      cur.pval = get.logrank.pvalue(check.survival(cur.clustering, subtype))
      return(cur.pval)
    }, mc.cores=1))
    
    total.num.perms = total.num.perms + num.perms
    total.num.extreme.perms = total.num.extreme.perms + sum(perm.pval <= orig.pvalue)
    
    binom.ret = binom.test(total.num.extreme.perms, total.num.perms)
    cur.pvalue = binom.ret$estimate
    cur.conf.int = binom.ret$conf.int
    
    print(c(total.num.extreme.perms, total.num.perms))
    print(cur.pvalue)
    print(cur.conf.int)
    
    sig.threshold = 0.05
    is.conf.small = ((cur.conf.int[2] - cur.pvalue) < min(cur.pvalue / 10, 0.01)) & ((cur.pvalue - cur.conf.int[1]) < min(cur.pvalue / 10, 0.01))
    is.threshold.in.conf = cur.conf.int[1] < sig.threshold & cur.conf.int[2] > sig.threshold
    if ((is.conf.small & !is.threshold.in.conf) | (total.num.perms > 2e7)) {
      should.continue = F
    } else {
      num.perms = 1e2
    }
  }
  
  print(list(pvalue = cur.pvalue, conf.int = cur.conf.int, total.num.perms=total.num.perms, 
                total.num.extreme.perms=total.num.extreme.perms))
  
}
