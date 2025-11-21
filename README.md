[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=dnafinder/testt)

🌐 Overview
testt.m implements Student’s t-test for two samples, supporting both unpaired and paired designs. For unpaired data, it automatically checks the equality of variances via a Fisher-Snedecor F-test and chooses between the classical pooled-variance t-test and the Satterthwaite (Welch) approximate t-test when variances differ. For paired data, it operates on the differences between paired observations, providing descriptive statistics, a confidence interval for the mean difference, the t statistic, the p-value and an approximate power estimate. The 2025 refactoring preserves the original statistical logic while improving robustness, clarity and documentation.

⭐ Features
- Unpaired two-sample t-test with:
  - Fisher-Snedecor F-test for equality of variances
  - Pooled-variance t-test when variances are equal
  - Satterthwaite (Welch) approximate t-test when variances differ
- Paired t-test on the differences:
  - Mean difference, standard deviation and standard error
  - (1 − alpha) two-sided confidence interval for the mean difference
- Flexible specification of:
  - Test type: unpaired (0) or paired (1)
  - Significance level alpha
  - One-tailed (1) or two-tailed (2) test
- Printed summary tables for:
  - F-test (unpaired case)
  - Paired differences (paired case)
  - Final t-test result with t, df, tail, alpha, p-value and power
- Returns a structure STATS with core test statistics for further processing.

🛠️ Installation
1. Clone or download the repository:
   https://github.com/dnafinder/testt
2. Add the folder containing testt.m to your MATLAB path:
   addpath('path_to_testt_folder')
3. You can also open and run testt.m directly in MATLAB Online using the badge at the top of this file.

▶️ Usage
Unpaired example:
    X1 = [77 79 79 80 80 81 81 81 81 82 82 82 82 83 83 84 84 84 84 85 ...
           85 86 86 87 87];
    X2 = [82 82 83 84 84 85 85 86 86 86 86 86 86 86 86 86 87 87 87 88 ...
           88 88 89 90 90];
    STATS = testt(X1, X2);              % default: unpaired, alpha=0.05, tail=1

Paired example:
    x1 = [60 70 40 41 40 40 45 48 30 50];
    x2 = [61 71 38 39 38 33 55 56 38 68];
    STATS = testt(x1, x2, 1, 0.05, 2);  % paired, alpha=0.05, two-tailed

The function prints relevant tables and a final summary line with t, df, tail, alpha, p-value and power.

🔣 Inputs
- x1
  Type: numeric vector (row or column).
  Description: first sample of observations.

- x2
  Type: numeric vector (row or column).
  Description: second sample of observations. For a paired test, x1 and x2 must have the same length.

- TST (optional)
  Type: scalar (0 or 1).
  Description:
    0 = unpaired test (default)
    1 = paired test

- ALPHA (optional)
  Type: scalar in (0,1).
  Description: significance level of the test. Default is 0.05. The confidence interval for the paired mean difference is always computed as (1 − ALPHA) two-sided.

- TAIL (optional)
  Type: scalar (1 or 2).
  Description:
    1 = one-tailed test
    2 = two-tailed test
  Default is 1.

📤 Outputs
The function returns a struct STATS when called with an output argument:

- STATS.tvalue
  The t statistic of the test.

- STATS.tdf
  The degrees of freedom used for the test:
  - For unpaired tests with equal variances: n1 + n2 − 2.
  - For unpaired tests with unequal variances: Satterthwaite’s approximate df.
  - For paired tests: n − 1, where n is the number of pairs.

- STATS.ttail
  The number of tails specified for the test (1 or 2).

- STATS.tpvalue
  The p-value associated with the t statistic, based on TAIL and ALPHA.

- STATS.power
  An approximate power estimate for the test, computed from the observed t value and the chosen tails.

For unpaired tests (TST = 0), additional fields are present:

- STATS.Fvalue
  F statistic for the test of equality of variances.

- STATS.DFn
  Numerator degrees of freedom of the F-test.

- STATS.DFd
  Denominator degrees of freedom of the F-test.

- STATS.FPvalue
  p-value of the F-test.

📘 Interpretation
- Unpaired case:
  - The F-test indicates whether the variances can be assumed equal.
  - If variances are equal, the pooled-variance t-test is used.
  - If variances differ, the Satterthwaite (Welch) approximation provides a robust alternative.
  - The sign and magnitude of the mean difference (in the internal code) indicate which group tends to have larger values, while the p-value indicates whether that difference is statistically significant.

- Paired case:
  - The test works on the differences d = x1 − x2.
  - The reported mean difference, standard deviation and standard error summarize how large and how variable the paired differences are.
  - The (1 − ALPHA) two-sided confidence interval for the mean difference provides a range of plausible values for the true mean difference.
  - If the CI does not include zero, and the p-value is below ALPHA, there is evidence of a systematic difference between conditions.

📝 Notes
- The refactored version (2025) keeps the original core formulas and decision rules while:
  - allowing both row and column vector inputs,
  - improving input checking and error messages,
  - clarifying the paired-case confidence interval formula.
- For paired data, the confidence interval is always interpreted as a two-sided (1 − ALPHA) interval, regardless of the one- or two-tailed choice of the test.
- The power calculations are heuristic and based on the observed t statistic and the chosen tails; they are not a substitute for a full a priori power analysis.

📚 Citation
If you use this function in scientific work, please cite:

Cardillo G. (2006). Student t-Test for unpaired or paired samples.  
GitHub repository: https://github.com/dnafinder/testt

👤 Author
Author : Giuseppe Cardillo  
Email  : giuseppe.cardillo.75@gmail.com  
GitHub : https://github.com/dnafinder  

⚖️ License
This code is distributed under the MIT License. You are free to use, modify and redistribute it, provided that proper credit is given to the original author and the source repository.
