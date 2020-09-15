# novelty-cenozoic-microplankton

Code used to generate results, figures and tables for DOI: <Pending>.

The data used to generate this publication were obtained from the Neptune Sandbox. These are not provided alongside this repository, but are freely available from http://nsb-mfn-berlin.de/.

Two R scripts are provided in this repository, along with summarised data used for for the statistical models and analyses in the manuscript. These scripts makes extensive use of the code folding functionality in RStudio. Alt + O is the default shortcut on Windows and Linux versions of RStudio to collapse all folds.

**novel_comms_neptune_analysis.R** is a minimum working script that contains all of the data and analyses to produce the figures and tables found in the main text. It does not show how the raw data were transformed from the Neptune Sandbox downloads. Many of the supplementary analyses require these raw data files, so code to reproduce these analyses have not been included in this file. Finally, this script makes use of files in respository sub-folders, so please point your working directory at the start of this script to its folder path, and include all other directories in the repository.

**novel_comms_neptune.R** shows the complete data analysis pathway, as well as main and supplementary analyses, but will require Neptune Sandbox exports for calcareous nannoplankton, foraminifera, radiolarians and diatoms. This script requires respository sub-folders etc as per the analysis script above.

We have also provided variants of functions that run all main text analyses using Bayesian hierarchical models, which provided very similar results to the frequentist models we included in the published manuscript.

[![DOI](https://zenodo.org/badge/288867070.svg)](https://zenodo.org/badge/latestdoi/288867070)
