# testt
Student's t test for unpaired or paired samples.<br/>
This file is applicable for equal or unequal sample sizes; for paired or 
unpaired samples. When the test is unpaired, the Fisher-Snedecor F-test is
performed to assess the equality of variance. If variances are not equal,
Satterthwaite's approximate t test is performed.

Syntax: 	TESTT(X1,X2,TST,ALPHA,TAIL)
     
    Inputs:
          X1 and X2 - data vectors (mandatory). 
          TST - unpaired (0) or paired (1) test (default = 0).
          ALPHA - significance level (default = 0.05).
          TAIL - 1-tailed test (1) or 2-tailed test (2). (default = 1).
    Outputs:
          - t value.
          - degrees of freedom.
          - Confidence interval of means difference (for paired test)
          - Critical value
          - p-value

     Example: 

          X1=[77 79 79 80 80 81 81 81 81 82 82 82 82 83 83 84 84 84 84 85 ...
          85 86 86 87 87];

          X2=[82 82 83 84 84 85 85 86 86 86 86 86 86 86 86 86 87 87 87 88 ...
          88 88 89 90 90];

          Calling on Matlab the function: testt

          Answer is:

FISHER-SNEDECOR F-TEST FOR EQUALITY OF VARIANCES
 
      F       DF_numerator    DF_denominator    p_value
    ______    ____________    ______________    _______

    1.5379    24              24                0.29861

--------------------------------------------------------------------------------
Variances are equal
--------------------------------------------------------------------------------
 
STUDENT'S T-TEST FOR UNPAIRED SAMPLES
 
      t       DF    tail     p_value      Power 
    ______    __    ____    _________    _______

    5.2411    48    1       1.765e-06    0.99958

STATS=TESTT(...) returns a structure with all test(s) statistics

          Created by Giuseppe Cardillo
          giuseppe.cardillo-edta@poste.it

To cite this file, this would be an appropriate format:
Cardillo G. (2006). Student t-Test for unpaired or paired samples.
http://www.mathworks.com/matlabcentral/fileexchange/12699
