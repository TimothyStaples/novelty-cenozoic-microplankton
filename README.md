# novelty-cenozoic-microplankton

Code used to generate results, figures and tables for DOI: XXXXX.

The data used to generate this publication were obtained from the Neptune Sandbox. These are not provided alongside this repository, but are freely available from http://nsb-mfn-berlin.de/.

Two R scripts are provided in this repository, along with summarised data used for for the statistical models and analyses in the manuscript. These scripts makes extensive use of the code folding functionality in R-studio (Alt + O is the default shortcut on Window machines to collapse all folds).

**novel-comms-neptune-analysis.R** is a minimum working script that contains all of the data and analyses to produce the figures and tables found in the manuscript. It does not show how the raw data was transformed from the Neptune Sandbox downloads. It makes use of files in the "Data" and "Functions" sub-folder, so please point your working directory to the location of this script and include all other directories in the repository.

**novel-comms-neptune.R** shows the complete data analysis pathway, but will require Neptune Sandbox exports for calcareous nannoplankton, foraminifera, radiolarians and diatoms.
