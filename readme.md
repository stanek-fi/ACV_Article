Replication repository for the article "Optimal Out-of-Sample Forecast Evaluation Under Stationarity"

This repository contains the code necessary for replicating the article "Optimal Out-of-Sample Forecast Evaluation Under Stationarity." It is structured as follows:
- The folder "ACV" contains the implementation of the estimator/tests proposed in the article. The implementation may differ slightly from the up-to-date official version available at https://cran.r-project.org/web/packages/ACV/index.html, as it reflects the state at the time of writing of the article and does not include any subsequent improvements/bug-fixes.
- The folder "Scripts" contains individual computational experiments.
- The folder "Article" contains the LaTeX file.

To reproduce the results:
- Install CRAN dependencies.
- To access the M4 competition time series, install the M4comp2018 library from GitHub via:
    - install.packages("https://github.com/carlanetto/M4comp2018/releases/download/0.2.0/M4comp2018_0.2.0.tar.gz", repos = NULL)
- Run the following scripts:
    - Scripts\Illustration\Illustration_1_AR1.R
    - Scripts\Level\Level.R
    - Scripts\MCompetitions\M4Run.R
    - Scripts\MCompetitions\M4GenTables.R
    - Scripts\Power\Power.R
    - Scripts\VarianceRatio\VarianceRatio.R

The scripts are parallelized, so feel free to adjust the number of workers according to your needs.
Files containing the results of each experiment will be generated for the next step.

- Compile the article:
    - Article\ArticleJFore_Revision_2\ArticleJFore.tex
