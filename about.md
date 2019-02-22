An IC 50 Calculator
=======================

IC 50 can be automatic calculated from laboratory data.

The calculator needs packaes "shiny" and "drc". Make sure you have these packages installed.

The data need to be in a table as .csv file.

1. Go to Upload file and browse to your data table.
2. Check options to make sure your data table look correct. THis should look as follow.
  1. first column should comprise the concentration
  2. Second column should contain the Inhibition values.
  3. Headers should be conz and value for the first and second column, respectively.
3. Go to generate IC50 tab
Select which model you want to keep (default = 1, which should be the best possible)
Choose a level of confidence for generating the confidence interval
click on "Calculate IC50".

The selected model will be plotted and the raw data are shown as blue stars.
Changes can be made and the output will be updated by clicking on "Calculate IC50".

If the dataset seem to have correctly been uploaded but no model can be fit, an output will be printed asking to check the values of your data.


2019 Dominic Rittler, Nelson Marreros

