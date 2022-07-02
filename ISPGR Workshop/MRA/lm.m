function [betas, pred, resid] = lm(x, y, order)
% [BETAS, PRED, RESID] = LM(X, Y, ORDER)
%
% Simple function for computing a linear model that returns coefficients,
% predicted values, and residuls.
%   Author: Aaron D. Likens (2022)
% Reference:
% Cohen, Cohen, West, & Aiken (2004).

    if nargin < 3
        order = 1;
    end
        
    betas = polyfit(x, y, order);
    pred = polyval(betas, x);
    resid = y - pred;
    
end
