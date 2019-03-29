# Meta-analysis: status signalling hypothesis in house sparrows

This repository contains the R scripts use in the following study:

[Alfredo Sánchez-Tójar, Shinichi Nakagawa, Moises Sánchez-Fortún, Sukanya Ramani, Dominic Martin, Antje Girndt, Veronika Bókony, Bart Kempenaers, András Liker, David F. Westneat Terry Burke, Julia Schroeder. 2018. *Meta-analysis challenges a textbook example of status signalling and demonstrates publication bias*. eLife, 7:e37385. DOI: 10.7554/eLife.37385](https://doi.org/10.7554/eLife.37385)

More information and materials (e.g. data) available at the [OSF project](http://doi.org/10.17605/OSF.IO/CWKXB). For any further information, please contact: [Alfredo Sánchez-Tójar](https://scholar.google.co.uk/citations?hl=en&user=Sh-Rjq8AAAAJ&view_op=list_works&sortby=pubdate), email: alfredo.tojar@gmail.com

## Scripts:

`001_reanalysing_dominance_and_bib.R`: script to analyzed raw data obtained from previous publications and unpublished raw data. For each dataset, it estimates the correlation between bib size and dominance rank for each group of house sparrows studied. Dominance rank was inferred using the [randomized Elo-rating method](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2656.12776) implemented in the R package '[aniDom](https://cran.r-project.org/web/packages/aniDom/index.html)'.

`002_metadatasets.R`: script to estimate effect sizes and sampling variances, and overall, to build the datasets that are subsequently used for the meta-analyses and meta-regressions.

`003_metaanalyses.R`: script to run the actual multilevel meta-analyses and meta-regressions to test whether the status signalling hypothesis across house sparrow populations. The analyses are run using the R package '[MCMCglmm](https://cran.r-project.org/web/packages/MCMCglmm/index.html)'. In addition, it contains code for testing publication and time-lag bias, and for generating most of the plots presented in [our publication](https://doi.org/10.7554/eLife.37385).

`004_funnelplots.R`: script to create funnel plots based on the meta-analytic residuals extracted from the meta-analyses run using the previous script.

`005_poweranalysis.Rmd`: Rmarkdown script explaining the steps of the simple power analysis included in the discussion of [our publication](https://doi.org/10.7554/eLife.37385).

### Notes:

29th March 2019: I have recently realized that some of the code used to estimate the meta-analytic residuals of the multilevel meta-analyses is no longer supported by more recent versions of the R package '[MCMCglmm](https://cran.r-project.org/web/packages/MCMCglmm/index.html)'.

