# HARSAT
Harmonized Regional Seas Assessment Tool


Requirements:
-	R programming language version 4.2.1
-	RStudio Integrated Development Environment version 2022.07.1 Build 554
-	Additional R packages to those which come with the standard RStudio installation. They have to be installed, either using the RStudio GUI or by the command install.packages e.g. install.packages(“lme4”).
o	rmvtnorm
o	tidyverse
o	lme4
o	numDeriv
o	optimx
o	pbapply
o	rmarkdown
o	htmlwidgets
o	DT
File -> New Project in RStudio to create a new project. File -> New Project -> Version Control -> GIT -> https://github.com/osparcomm/HARSAT to clone from repository.
All datasets can be run with the same directory structure (although some of the directories will be redundant). The following directories should be manually created:
- data
- functions
- information
- output
     - summary
- RData
    - biota assessment components
    - sediment assessment components
    - water assessment components
The oddities directory will be created automatically when running the code.

![Alt text](images/fig1.jpg?raw=true "Required Project File Structure for OSPAR dataset.")
Figure 1. Required Project File Structure for OSPAR dataset.

![Alt text](images/fig2.jpg?raw=true "Required Project File Structure for HELCOM dataset.")
Figure 2. Required Project File Structure for HELCOM dataset.

All R functions from the GitHub code repository should be moved to ‘functions’, and csv files with the reference tables should go to ‘information’. 
It is recommended that the R files with examples, and the data and output directories should be in the same (project) directory e.g.
HARSAT\data
HARSAT\output\summary
HARSAT\OSPAR 2022.r
although one can adjust that with the path argument to ctsm_read_data and ctsm.summary.table.
The ‘functions’ and ‘information’ directories can be anywhere (although it is recommended that they got to the same directory), with the function_path variable (at the top of OSPAR 2022.r) pointing to the former.
Data required for OSPAR:
 
![Alt text](images/fig3.jpg?raw=true "OSPAR data.")
Figure 3. OSPAR data.

Data required for HELCOM:
 
![Alt text](images/fig4.jpg?raw=true "HELCOM data.") 
Figure 4. HELCOM data.

HELCOM dataset requires additional R files to work, and different reference tables for ‘assessment criteria biota.csv’ and ‘assessment criteria sediment.csv’.
HELCOM also requires different information functions from the repository with a different version that has HELCOM specific functions (‘information functions v2_68.r’).

![Alt text](images/fig5.jpg?raw=true "Additional R files required for HELCOM dataset.") 
Figure 5. Additional R files required for HELCOM dataset.

