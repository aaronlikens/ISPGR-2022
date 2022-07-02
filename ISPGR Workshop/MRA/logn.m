function logx_n = logn(x, n)
% LOGX_N = LOGN(X, N)
%
%  Change of Logarthimc Bases
%
% log(x, n) takes the logarithm of an abitrary value n by applying the
% change of base formula to a natural logarithm. This is also a helper
% function in dfa and mfdfa.
%
% Input parameter:
% x is a real valued vector of values
%
% n is a real number >= 0
%
% Output parameter:
% logx_n 
%
%   Example:
%       x = abs(rand(100,1));
%       n = 1.1;
%       logx_n = logn(x, n);
%   Author: Aaron D. Likens (2022)

% function for change of base formula.
logx_n = log(x)/log(n);
    
end