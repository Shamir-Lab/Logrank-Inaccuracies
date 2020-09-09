# Heinze's log-rank test
This repository contains an implementation of Heinze's conditional log-rank test.
Log-rank is a widely used test for differential survival. The widespread version of the test, implemented in R's survival package, is asymptotic. Its accuracy not only depends on the number of samples, but also on the number of events (deaths in the case of survival data) in each group. More details on inaccuracies of the log-rank test can be found in: https://doi.org/10.15252/msb.20188754.

This package is an R implementation of Heinze's version of the conditional log-rank test, which uses permutations of the data and does not rely on asymptotic assumptions. The implementation is based on code from the permGS R package, but supports more than two groups.

The package can be installed as follows:
```{r}
devtools::install_github("Shamir-Lab/Logrank-Inaccuracies/logrankHeinze")
```

The package includes one main function, logrank.heinze, and one dataset, containing survival data of colon cancer patients from TCGA (with an arbitrary partition to two groups). 

The basic usage is as follows:
```{r}
# loads colon.survival.data
data(ColonSurvivalData)
logrankHeinze::logrank.heinze(Survival + Death ~ cluster, data=colon.survival.data, max.num.perms = 1000)
# or alternatively:
logrankHeinze::logrank.heinze(data=colon.survival.data, max.num.perms = 1000)
```
The above example runs 1000 permutations to estimate the log-rank p-value. The package supports several termination conditions in addition to the number of iterations, such as the size of the 95% confidence interval for the p-value, and a significance threshold which should not be contained in the confidence interval. The package also supports parallelization, and reproducibility using seeds (even when using parallelization). All of these options and others are documented in the package.

# Reproducing the paper's results
Additionally, this repository contains code and data to reproduce the results presented in the "Inaccuracy of the logrank approximation in cancer data analysis" paper. These appear in the "paper_code_and_data" directory.

The implementation of the permutation-based method to compute the exact log-rank p-value is based on code from the permGS R package:
https://cran.r-project.org/web/packages/permGS/index.html.
It requires to install the R survival package, and to load the parallel package.

The data needed to reproduce the TCGA analyses are here: http://acgt.cs.tau.ac.il/multi_omic_benchmark/download.html.

The data needed to reproduce the clusternomics and murine model of endotoxemia are available in this repository.

# Contact
You can contact me at nimrodrap@gmail.com.

# Citation

    @article{Rappoport2019,
    abstract = {The log-rank test statistic is very broadly used in biology. Unfortunately, P-values  based on the popular chi-square approximation are often inaccurate and can be misleading.},
    author = {Rappoport, Nimrod and Shamir, Ron},
    doi = {10.15252/msb.20188754},
    issn = {1744-4292},
    journal = {Molecular Systems Biology},
    month = {aug},
    number = {8},
    publisher = {EMBO},
    title = {{Inaccuracy of the log‐rank approximation in cancer data analysis}},
    url = {https://onlinelibrary.wiley.com/doi/abs/10.15252/msb.20188754},
    volume = {15},
    year = {2019}
    }


