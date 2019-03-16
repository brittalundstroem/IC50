An IC 50 Calculator
=======================

Calculate IC50 and confidence interval.

The data need to be imported as a .csv file. The first column should contain concentration values and the second column should contain the inhibition values. Make sure your table follows this structure before you upload the file.

1. Go to `Upload File` and browse to upload your dataset.
2. Check the printed table to ensure your dataset has the right structure. It is possible to use different column names, but the *order of the columns* should remain the same.
    1. If your dataset has no header, remove the checkmark from the `Heather` checkbox.
    2. If the cell separator is different than a semicolon (;), choose another option.
3. Go to the `Generate IC50` tab
    1. click on `Check the dataset`.
      1. If a model can be fitted a list of models is printed below the button. First model on the list is the best according to the Akaike information criterion. The second model on the list is the second best and so on.
      2. If no model can be fitted, a message asking to check the dat	aset is printed.
    1. Select which model you want to fit from the selection box.
    2. Choose an alpha-error level for generating the confidence interval, default is 95%.
4. click on `(Re-)Calculate IC50`.

The model will be plotted with the raw data.
A summary of the model will be printed below the plot.

I strongly recommended to compare several models, since the AIC is no more than an automated selection process. After changes are made, the plot and ouptut can be updated by clicking again on `(re-)calculate IC50`.

The visible plot and model output can be downloaded in a .zip archive if desired.



2019 Dominic Rittler, Nelson Marreros
