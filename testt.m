function STATS = testt(x1, x2, varargin)
%TESTT Student's t-test for unpaired or paired samples.
%
%   TESTT(X1, X2) performs a Student's t-test on the data vectors X1 and X2.
%   It supports:
%       - unpaired two-sample t-test with equal or unequal sample sizes
%       - paired t-test when measurements are paired by subject
%
%   When the test is unpaired, a Fisher-Snedecor F-test is performed to
%   assess the equality of variances. If variances are unequal, the
%   Satterthwaite (Welch) approximate t-test is used. If variances are
%   equal, the classical pooled-variance t-test is applied.
%
%   TESTT(X1, X2, TST, ALPHA, TAIL) allows full control:
%       TST   : 0 = unpaired (default), 1 = paired
%       ALPHA : significance level of the test (default = 0.05)
%       TAIL  : 1 = one-tailed test, 2 = two-tailed test (default = 1)
%
%   For the paired test (TST = 1), the function reports descriptive
%   statistics on the differences (x1 - x2), including mean, standard
%   deviation, standard error and the (1 - ALPHA) two-sided confidence
%   interval of the mean difference.
%
%   STATS = TESTT(...) returns a structure with test statistics.
%
%   ------------------------------------------------------------------
%   Syntax:
%       testt(x1, x2)
%       testt(x1, x2, TST)
%       testt(x1, x2, TST, ALPHA)
%       testt(x1, x2, TST, ALPHA, TAIL)
%       STATS = testt(...)
%
%   Inputs:
%       x1      - Numeric vector (row or column), first sample.
%       x2      - Numeric vector (row or column), second sample.
%       TST     - (optional) test type:
%                   0 = unpaired (default)
%                   1 = paired
%       ALPHA   - (optional) significance level of the test:
%                   0 < ALPHA < 1, default = 0.05.
%                 The confidence interval for paired samples is always
%                 computed as (1 - ALPHA) two-sided (e.g. 95% CI).
%       TAIL    - (optional) number of tails for the test:
%                   1 = one-tailed
%                   2 = two-tailed (default = 1)
%
%   Outputs:
%       STATS   - Structure with fields:
%                   STATS.tvalue  : t statistic
%                   STATS.tdf     : degrees of freedom
%                   STATS.ttail   : number of tails (1 or 2)
%                   STATS.tpvalue : p-value
%                   STATS.power   : approximate power of the test
%
%                 For unpaired tests (TST = 0), the structure also
%                 contains:
%                   STATS.Fvalue  : F statistic for equality of variances
%                   STATS.DFn     : numerator degrees of freedom (F-test)
%                   STATS.DFd     : denominator degrees of freedom (F-test)
%                   STATS.FPvalue : p-value of the F-test
%
%   Example:
%       % Unpaired test
%       X1 = [77 79 79 80 80 81 81 81 81 82 82 82 82 83 83 84 84 84 84 85 ...
%              85 86 86 87 87];
%       X2 = [82 82 83 84 84 85 85 86 86 86 86 86 86 86 86 86 87 87 87 88 ...
%              88 88 89 90 90];
%       STATS = testt(X1, X2, 0, 0.05, 1);
%
%       % Paired test
%       x1 = [60 70 40 41 40 40 45 48 30 50];
%       x2 = [61 71 38 39 38 33 55 56 38 68];
%       STATS = testt(x1, x2, 1, 0.05, 1);
%
%   ------------------------------------------------------------------
%   Metadata:
%       Author : Giuseppe Cardillo
%       Email  : giuseppe.cardillo.75@gmail.com
%       GitHub : https://github.com/dnafinder
%       Created: 2006-01-01
%       Updated: 2025-11-21
%       Version: 2.0.0
%
%   Citation:
%       Cardillo G. (2006). Student t-Test for unpaired or paired samples.
%       GitHub: https://github.com/dnafinder/testt
%
%   License:
%       This code is distributed under the MIT License.
%   ------------------------------------------------------------------

%% Input parsing
p = inputParser;
p.FunctionName = 'testt';

addRequired(p, 'x1', @(x) validateattributes( ...
    x, {'numeric'}, {'vector','real','finite','nonnan','nonempty'}));

addRequired(p, 'x2', @(x) validateattributes( ...
    x, {'numeric'}, {'vector','real','finite','nonnan','nonempty'}));

addOptional(p, 'tst',   0, @(x) isnumeric(x) && isscalar(x) && ismember(x, [0 1]));
addOptional(p, 'alpha', 0.05, @(x) validateattributes( ...
    x, {'numeric'}, {'scalar','real','finite','nonnan','>',0,'<',1}));
addOptional(p, 'tail',  1, @(x) isnumeric(x) && isscalar(x) && ismember(x, [1 2]));

parse(p, x1, x2, varargin{:});
tst   = p.Results.tst;
alpha = p.Results.alpha;
tail  = p.Results.tail;
clear p

% Ensure column vectors
x1 = x1(:);
x2 = x2(:);

if tst == 1
    assert(numel(x1) == numel(x2), ...
        'testt:SizeMismatch', ...
        'For a paired test, x1 and x2 must have the same length.');
end

tr = repmat('-', 1, 80);

%% Main switch: unpaired vs paired
switch tst
    case 0  % unpaired test
        n = [numel(x1), numel(x2)];    % sample sizes
        m = [mean(x1),  mean(x2)];     % sample means
        v = [var(x1),   var(x2)];      % sample variances
        
        % Summary of groups
        disp(table([n(1); m(1); v(1)], [n(2); m(2); v(2)], ...
            'VariableNames', {'Group_1','Group_2'}, ...
            'RowNames', {'Numerosity','Mean','Variance'}));
        
        % Fisher-Snedecor F-test for equality of variances
        % Ensure F >= 1 by swapping groups if needed
        if v(2) > v(1)
            v = fliplr(v);
            m = fliplr(m);
            n = fliplr(n);
        end
        
        F  = v(1) / v(2);       % variance ratio
        DF = n - 1;             % degrees of freedom for numerator/denominator
        pF = fcdf(F, DF(1), DF(2));
        pF = 2 * min(pF, 1 - pF);   % two-sided p-value for F
        
        if nargout
            STATS.Fvalue  = F;
            STATS.DFn     = DF(1);
            STATS.DFd     = DF(2);
            STATS.FPvalue = pF;
        end
        
        % Display F-test results
        disp('FISHER-SNEDECOR F-TEST FOR EQUALITY OF VARIANCES')
        disp(' ')
        disp(array2table([F DF(1) DF(2) pF], ...
            'VariableNames', {'F','DF_numerator','DF_denominator','p_value'}));
        disp(tr)
        
        if pF < alpha
            % Unequal variances (Behrens-Welch / Satterthwaite case)
            fprintf('Variances are different: Behrens-Welch problem\n');
            disp(tr); disp(' ');
            
            % Satterthwaite's approximate t-test
            a     = v ./ n;
            b     = sum(a);
            denom = sqrt(b);
            gl    = b^2 / sum(a.^2 ./ (n - 1));   % effective degrees of freedom
            disp('SATTERTHWAITE''S APPROXIMATE T-TEST FOR UNPAIRED SAMPLES')
            disp(' ')
        else
            % Equal variances (classical pooled-variance t-test)
            fprintf('Variances are equal\n');
            disp(tr); disp(' ');
            
            gl    = sum(n) - 2;                        % degrees of freedom
            s     = sum((n - 1) .* v) / (sum(n) - 2);  % pooled variance
            denom = sqrt(sum(s ./ n));                 % SE of mean difference
            disp('STUDENT''S T-TEST FOR UNPAIRED SAMPLES')
            disp(' ')
        end
        
        dm = diff(m);    % difference of means (m2 - m1)
        clear n m v a b s  % clean up
        
    case 1  % paired test
        disp('STUDENT''S T-TEST FOR PAIRED SAMPLES')
        disp(' ')
        disp(tr)
        
        n  = numel(x1);         % number of pairs
        gl = n - 1;             % degrees of freedom
        d  = x1 - x2;           % differences
        dm = mean(d);           % mean difference
        sd = std(d, 0);         % sample standard deviation (denominator n-1)
        
        % Standard error of the mean difference
        denom = sd / sqrt(n);
        
        % Critical value for (1-ALPHA) two-sided CI (e.g. 95% CI if alpha=0.05)
        vc = tinv(1 - alpha/2, gl);
        
        % Confidence interval for mean difference
        ci = dm + [-1 1] .* vc .* denom;
        
        % Table with N, mean, SD, SE, CI
        TBL = table(n, dm, sd, denom, ci, ...
            'VariableNames', {'N_pairs','Mean_difference','Std_dev','Std_error','Conf_interval'});
        disp(TBL)
        
        clear n d  % clean up
end

%% t statistic and p-value
t = abs(dm) / denom;          % t value
p = (1 - tcdf(t, gl)) * tail; % p-value according to tail

%% Power estimation
switch tail
    case 1  % one-tailed
        if p >= alpha
            tp    = tinv(1 - alpha, gl) - t;
            Power = 1 - tcdf(tp, gl);
        else
            tp    = t - tinv(1 - alpha, gl);
            Power = tcdf(tp, gl);
        end
    case 2  % two-tailed
        if p >= alpha
            tp1   = tinv(1 - alpha/2, gl) - t;
            tp2   = t + tinv(1 - alpha/2, gl);
            Power = 2 - tcdf(tp1, gl) - tcdf(tp2, gl);
        else
            tb1   = t - tinv(1 - alpha/2, gl);
            tb2   = t + tinv(1 - alpha/2, gl);
            Power = 1 - (tcdf(tb2, gl) - tcdf(tb1, gl));
        end
end

% Display final t-test summary
disp(array2table([t gl tail alpha p Power], ...
    'VariableNames', {'t','DF','tail','alpha','p_value','Power'}));

%% Output struct
if nargout
    STATS.tvalue  = t;
    STATS.tdf     = gl;
    STATS.ttail   = tail;
    STATS.tpvalue = p;
    STATS.power   = Power;
end

end
