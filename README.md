README file for R package supporting the paper "scTSSR2: an R package for imputing dropout events in single-cell RNA sequencing data using fast two-side self-representation".

The scTSSR2 package has the following R-package dependencies: SAVER, keras, tensorflow. The dependent packages will be automatically installed along with scTSSR2. You can use the following commands to install scTSSR2 from GitHub.

# Step 1. If the devtools package has been not installed, install the devtools package first. Invoke R and then type

install.packages("devtools")

# Step 2. Load the devtools package.

library("devtools")

# Step 3. Install the scTSSR2 package from GitHub.

install_github("Zhangxf-ccnu/scTSSR2")

## Useage

Load the library scTSSR2 in R console, by running

library(scTSSR2)

Taking the baron dataset as an example, run the following code:

data("baron")

baron\_imputation_result <- scTSSR2(baron$count.samp)

For detialed usages, please refer to "scTSSR2-manual.pdf".

Please do not hesitate to contact Miss Ke Jin ([kej13@mails.ccnu.edu.cn](kej13@mails.ccnu.edu.cn)) or Dr. Xiao-Fei Zhang ([zhangxf@mail.ccnu.edu.cn](zhangxf@mail.ccnu.edu.cn)) to seek any clarifications regarding any contents or operation of the archive.