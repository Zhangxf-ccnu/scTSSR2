README file for R package supporting the paper **"scTSSR2: an R package for imputing dropout events in single-cell RNA sequencing data using fast two-side self-representation"**.

The `scTSSR2` package has the following R-package dependencies: `SAVER`, `keras`, `tensorflow`. The dependent packages will be automatically installed along with `scTSSR2`. You can use the following commands to install `scTSSR2` from GitHub.

Installation
----------------------
**Step 1.** If the devtools package has been not installed, install the devtools package first. 

`install.packages("devtools")`

**Step 2.** Load the devtools package.

`library("devtools")`

**Step 3.** Install the scTSSR2 package from GitHub.

`install_github("Zhangxf-ccnu/scTSSR2")`

Useage
----------------------
Load the library scTSSR2 in R console, by running

`library(scTSSR2)`

Taking the baron dataset as an example, run the following code:

`data("baron")`

`baron_imputation_result <- scTSSR2(baron$count.samp)`

For detialed usages, please refer to "scTSSR2-manual.pdf".

Codes for reproducing the three downstream analyses (such as differential expression analysis, cell clustering analysis and pseudotime analysis) are available in [scTSSR-scTSSR2_experiments_codes](https://github.com/Zhangxf-ccnu/scTSSR-scTSSR2_experiments_codes) and Zenodo website with [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6423597.svg)](https://doi.org/10.5281/zenodo.6423597).

Contact
------------------------
Please do not hesitate to contact Miss Ke Jin ([kej13@mails.ccnu.edu.cn](kej13@mails.ccnu.edu.cn)) or Dr. Xiao-Fei Zhang ([zhangxf@mail.ccnu.edu.cn](zhangxf@mail.ccnu.edu.cn)) to seek any clarifications regarding any contents or operation of the archive.
