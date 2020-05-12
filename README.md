# niwot_NEON_carabids
Prediction of carabid species abundance at Niwot Ridge.
End-of-semester (Spring 2020) presentation: [.md](https://github.com/EBIO6100Spring2020/niwot_NEON_carabids/blob/master/final_pres_done.md), [.Rmd](https://github.com/EBIO6100Spring2020/niwot_NEON_carabids/blob/master/final_pres_done.Rmd)

Subdirectories

data_raw - put raw data here and never modify it
data_derived - derived datasets (e.g. cleaned, sampled, reworked)
source - any custom functions sourced by scripts
output - temporary figures, tables, etc that you want to save
docs - documentation and resources
Scripts (.R, .Rmd) can live at the top level. Name scripts by project keyword, function, and perhaps include initials if you are working independently to start with (e.g. ticks_eda_bam.R). As a pipeline develops, include the sequence structure (e.g. carabid_01_download.R, carabid_02_clean.R, carabid_03_eda.R,...).

To set it up
Clone the repository to your computer using R studio as described in Happy Git with R 12.3. See the Git tutorials if you need a refresher.

Coding tips
Refer to files in scripts using relative rather than absolute paths.
Use relative paths so that code will work from any location on any computer. Don't use absolute paths in scripts such as

C:/user/jane/janescoolstuff/experiment2/data_raw/neon_ants.csv
This will break the script on a different user's computer. Instead, use relative paths, such as

data_raw/neon_ants.csv
Anyone can then run the code without needing to modify the file paths. This is especially important when collaborating via a repository.

Don't use setwd() in scripts.
This is not portable and will break the script on another person's computer. If you set up an RStudio project (as above), you will be in the correct working directory when you start RStudio by opening the project.

Don't save your R workspace
Start clean each time. RStudio setup: In Tools > Global options > General, set "save workspace" to "Never" and uncheck everything except "Automatically notify me of updates to RStudio". This ensures that all your work derives from code and provides a test of the code each time you work on the script.

Use the Issues tab in this repository to post questions about data, make TO DO lists etc.
