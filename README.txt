The R code and results are divided into two directories, called `simulationstudy' and `asymptoticstudy'.  So, this README.txt file is divided into two parts.  The first part concerns the simulationstudy directory; the second part concerns the asymptoticstudy directory.  In addition to these two directories, you will see the file `copulaplots.R'.  The code in that file reproduces Figure 2 in the article.


simulationstudy directory
-------------------------
-------------------------

Reproducing Tables 1-3 and Table A7
=================================================================

The full set of results shown in Tables 1-3 and Table A7 was generated using a high performance cluster.  The full set of analyses could take several days to run on a standard computer.

To reproduce the results, run the R code in the file `simulation_master.R'.
This calls the files `simulation_analysis.R' and `simulation_saveresults.R' and `simulation_createtables.R'.
The file `simulation_analysis.R' calls the file `simulation_data_generate.R'.
Files containing intermediate results will be placed in a subdirectory called `intermediate_results'.
The files containing final results will be placed in a subdirectory called `results'.


Quickly testing the code
========================

Running the code in the `simulation_master.R' file could take several days.

In order more quickly to check the code is working as it should, you can reduce the number of simulations to 20 and the number of bootstraps to 200.
To do this, change "nsim<-1000" on line 23 of `simulation_master.R' to "nsim<-20" and change "nboot<-1000" on line 26 to "nboot<-200".
The results we obtained when we did this can be found in the subdirectory `results_fewsimulations' (with intermediate results in the subdirectory `intermediate_results_fewsimulations').


sessionInfo() for Tables 1-3 and Table A7
=========================================

When the code for reproducing Tables 1-3 and Table A7 was run on a Dell laptop, the information from sessionInfo() was as follows:

R version 4.4.1 (2024-06-14)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 20.04.6 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0 
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0

locale:
 [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
 [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8    LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       

time zone: Europe/London
tzcode source: system (glibc)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] dplyr_1.1.4    sandwich_3.0-2 survival_3.7-0

loaded via a namespace (and not attached):
 [1] utf8_1.2.3       R6_2.5.1         tidyselect_1.2.0 Matrix_1.6-5     lattice_0.22-5   magrittr_2.0.3  
 [7] splines_4.4.1    glue_1.6.2       zoo_1.8-12       tibble_3.2.1     pkgconfig_2.0.3  generics_0.1.3  
[13] lifecycle_1.0.3  cli_3.6.2        fansi_1.0.4      grid_4.4.1       vctrs_0.6.5      compiler_4.4.1  
[19] tools_4.4.1      pillar_1.9.0     rlang_1.1.0 



asymptoticstudy directory
-------------------------
-------------------------

Reproducing Tables A1 - A3
==========================

To reproduce Tables A1-A3, run the R code in the file `simulateMSM_runall.R'.
This will recreate the files `Table_A1_samplesize1e+06.txt', `Table_A2_samplesize1e+06.txt', and `Table_A3_samplesize1e+06.txt' provided in the `results_previous' subdirectory.
The new results files will be created in the `results' subdirectory.
The code provided uses a sample size of n=1000000 and will take several hours to run.

Quickly testing the code
------------------------

For purpose of checking that the code is working as it should, you can reduce the running time by using a smaller sample size.
To do this, edit line 3 of `simulateMSM_runall.R' to change "n <- 1e6" to "n <- 50000".
This will recreate the files `Table_A1_samplesize50000.txt', `Table_A2_samplesize50000.txt', and `Table_A3_samplesize50000.txt' provided in the `results_previous' subdirectory.
The new results files will be created in the `results' subdirectory.


Reproducing Tables A5 and A6
============================

To reproduce Tables A5 and A6, run the R code in the file `simulateMSM_sensitivity_runall.R'.
This will recreate the files `Table_A5_samplesize1e+05.txt' and `Table_A6_samplesize1e+05.txt' provided in the `results_previous' subdirectory.
However, this uses a sample size of n=100000 and will take several hours to run.
The new results files will be created in the `results' subdirectory.


Quickly testing the code
------------------------

For purpose of checking that the code is working as it should, you can reduce the running time by using a smaller sample size.
To do this, edit line 3 of `simulateMSM_sensitivity_runall.R' to change "n <- 100000" to "n <- 5000".
Running the resulting code will recreate the files `Table_A5_samplesize5000' and `Table_A6_samplesize5000' provided in the `results_previous' subdirectory
The new results files will be created in the `results' subdirectory.


sessionInfo() for Tables A1-A3, A5 and A6
=========================================

When the code for reproducing Tables A1-A3, A5 and A6 was run on a Dell laptop, the information from sessionInfo() was as follows:

R version 4.4.0 (2024-04-24)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 20.04.6 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0 
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0

locale:
 [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
 [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
 [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       

time zone: Europe/London
tzcode source: system (glibc)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] sandwich_3.0-2

loaded via a namespace (and not attached):
[1] zoo_1.8-12     compiler_4.4.0 tools_4.4.0    grid_4.4.0     lattice_0.22-5
