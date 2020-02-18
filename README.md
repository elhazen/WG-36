# WG-36

This repository includes three test scripts and a test dataset from the California Current to be tested by Yokohama PICES meeting. Please try all the code on the test dataset (and your dataset if you can) before the meeting, and bring as many indicators that align with WG28's summaries that you can to the workshop. Specifically, this includes:

1) Test data - "coastwide data for reference points.csv"
2) Dynamic Factor Analysis code - "DFA code v2.0_ELH.R" that will create plots and create a clean dataset that will be used in the next two scripts.
3) Single factor Generalized Additive Model code "Single_Driver_ResponseGAM_v2.R" - GAM and inflections to identify non-linear thersholds.
4) Gradient Forest Analysis code - "gradientForestAnalysis.R" - a multi-factor regression approach to identify non-linear thresholds. There is also a function that can be used in your own code-writing that Scott graciously put together "gradForestFunctions.R"
5) INDperform package - "IndicatorPerformancePackage.R" - A method for testing redundancy and utility of indicators presented by Saskia Otto at the ECCWO W11.


Updates:

5/29/18 - adding Samhouri et al. 2017 EcoSphere code for DFA time series trends - https://people.ucsc.edu/~elhazen/hazen/Publications_files/Samhouri_et_al-2017-Ecosphere.pdf

6/3/18 - added additional code from Large et al. 2015 - MEPS and Saskio et al. 2018 - Ecol Indicators and updated this ReadMe.

2/18/20 - Saskia: run of INDperform analysis to certain extent -> needs further discussion on a) dealing with temporal autocorrelation (difficult to identify as well as handle since time series is so short) and b) choice of threshold variable and interaction testing (all combinations lead to 1890 models to test for) --> see also my comments in the R script.
