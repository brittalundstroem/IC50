An IC 50 Calculator
=======================

The calculator needs packages "shiny" and "drc" to work properly.

The data need to be imported as a .csv file. The concentration has to be in the first column and the inhibition values in the second. Make sure your table follows this structure before you upload the file.

1. Go to `Upload File` and browse to upload your dataset.
2. Check the printed table to ensure your dataset has the right structure. It is possible to use different column names, but the *order of the columns* should remain the same.
    1. If your dataset has no header, remove the checkmark from the `Heather` checkbox.
    2. If the cell separator is different than a semicolon (;), choose another option.
3. Go to the `Generate IC50` tab
    1. Select which model you want to fit: 1 means the best model according to the Akaike information criterion.
    2. Choose an alpha-level of error for generating the confidence interval, default is 95%
4. click on `(Re-)Calculate IC50`.

The fitted line and data will be plotted. In addition, the raw data are plotted as blue stars.


It is strongly recommended to compare several models, since the AIC is no more than an automated selection process. After changes are made, the plot and ouptut can be updated by clicking again on `(re-)calculate IC50`.

If the dataset seems to be uploaded but no model can be fit, an output will be printed asking to check the values of your data.




2019 Dominic Rittler, Nelson Marreros

