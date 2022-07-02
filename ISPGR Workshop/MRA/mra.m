function [s, fsx, fsy, fsxy, betaxy] = mra(x, y, order, mins, maxs, ...
    scale_ratio, overlap)
% [S, FSX, FSY, FSXY, BETAXY] = MRA(X, Y, ORDER, MINS, MAXS,...
%   SCALE_RATIO, OVERLAP)
% 
%   Multiscale (Fractal) Regression Analysis
%
%   mra(x, y, order, mins, maxs, scale_ratio, overlap) calcuations
%   multiscale regression analysis for a bivariate time series using either
%   overlapping or, more typically, non-overlapping windows.
%   (see references).
%
%   Input parameters:
%   x is a real valued time series, of equal length to y. x is taken to be
%   the independent variable.
%
%   y is a real valued time series, of equal length to x. y is taken to be
%   the dependent variable.
%
%   order is the polynomial order of the detrending function used within
%   local windows.
%
%   mins is an integer indicating the minimum time scale over which to
%   estimate the fluctuation function
%
%   maxs is an integer indicating the maximum time scale over which to 
%   estimate the fluctuation function
%
%   scale_ratio is real number > 1 that determines the logarithmic distance 
%   between neighboring scales
%
%   overlap is a logical indicating whether overlapping or non-overlapping
%   windows should be used in calculating scale-wise fluctuations
%
%   Output Parameters:
%   s is a vector containing the scales used in computing the fluctuation
%   function
%
%   fsx is a vector containing the fluctuation function for x estimated
%   over s
%
%   fsy is a vector containing the fluctuation function for x estimated
%   over s
% 
%   fsxy is a vector containing the detrended covariance function between x
%   and y estimaated at each scale s
%
%   betaxy is a vector containing the scale-wise regression coefficients 
%   where x predicts y at each scale, s
%   
%   Example:
%       t = 10000;
%       x = fgn_sim(t, .8)';
%       y = .7*x + randn(t, 1);
%       order = 1;
%       mins = 16;
%       maxs = floor(n/4);
%       scale_ratio = 2;
%       overlap = false;
%       [s, fsx, fsy, fsxy, betaxy] = mra(x, y, order, mins, maxs, ...
%           scale_ratio, overlap)
%       plot(s, betaxy)
%       xlabel('s'); ylabel('rho(x,y)');
%   Author: Aaron D. Likens (2022)
%
%   References:
%   Kristoufek, L. (2015). Detrended fluctuation analysis as a regression
%   framework: Estimating dependence at different scales. Physical Review 
%   E, 91(2), 022802.
%
%   Likens, A. D., Amazeen, P. G., West, S. G., & Gibbons, C. T. (2019).
%   Statistical properties of Multiscale Regression Analysis: Simulation 
%   and application to human postural control. Physica A: Statistical 
%   Mechanics and its Applications, 532, 121580.

n = length(x);

% create a scale vector that is evenly spaced in the log domain
s=0;
counter=0;
while s <= maxs
    counter=counter+1;
    s(counter,1)=ceil(scale_ratio^counter);
end
s = s(s >= mins & s <= maxs);
s = sort(unique(s));

% allocate output vector
fsx     = zeros(length(s), 1);
fsy     = zeros(length(s), 1);
fsxy    = zeros(length(s), 1);

% create the profile
X = cumsum(x-mean(x));
X_rev = flipud(X);

Y = cumsum(y-mean(y));
Y_rev = flipud(Y);

% choose between overlapping and non-overlapping DFA
if overlap == true
    for i = 1:length(s)
        s_i = s(i);
        ind = (1:s_i)';
        
        %determine number of windows for a scale, s
        nwins = n - s_i;
        
        % allocate vectors for variances and covariance
        f2sx = zeros(nwins, 1);
        f2sy = zeros(nwins, 1);
        f2xy = zeros(nwins, 1);

        for j = 1:nwins
            % obtain variances for x
            [~,~,resx] = lm(ind, X(ind), order);
            f2sx(j) = sum(resx.^2)/(s_i - 1);
            
            % obtain variances for y
            [~,~,resy] = lm(ind, Y(ind), order);
            f2sy(j) = sum(resy.^2)/(s_i - 1);
            
            % obtain covariances between X and Y
            covxy = resx.*resy;
            f2xy(j) = sum(covxy)/(s_i - 1);
            
            % increment window
            ind = ind + 1;
        end
        
        fsx(i) = sum(f2sx)/nwins;
        fsy(i) = sum(f2sy)/nwins;
        fsxy(i) = sum(f2xy)/nwins;
    end
    
    % compute the scale-wise correlation coefficient
    betaxy = fsxy./fsx;
    fsx = sqrt(fsx);
    fsy = sqrt(fsy);
    
else
    
    % loop through each scale, s
    for i = 1:length(s)
        
        % loop through each window at each scale s
        % also reverse x to perform procedure in each direction
        s_i = s(i);
        ind = (1:s_i)';
        nwins = floor(n/s_i);
        f2sx = zeros(nwins, 1);
        f2sx_rev = zeros(nwins,1);
        f2sy = zeros(nwins, 1);
        f2sy_rev = zeros(nwins,1);
        f2xy = zeros(nwins,1);
        f2xy_rev = zeros(nwins,1);
        for j = 1:nwins
            
            % obtain variances for x
            [~,~,resx] = lm(ind, X(ind), order);
            f2sx(j) = sum(resx.^2)/(s_i - 1);
            [~,~,resx_rev] = lm(ind, X_rev(ind), order);
            f2sx_rev(j) = sum(resx_rev.^2)/(s_i - 1);
            
            % obtain variances for y
            [~,~,resy] = lm(ind, Y(ind), order);
            f2sy(j) = sum(resy.^2)/(s_i - 1);
            [~,~,resy_rev] = lm(ind, Y_rev(ind), order);
            f2sy_rev(j) = sum(resy_rev.^2)/(s_i - 1);
            
            
            % obtain covariances between X and Y
            covxy = resx.*resy;
            f2xy(j) = sum(covxy)/(s_i - 1);
            covxy_rev = resx_rev.*resy_rev;
            f2xy_rev(j) = sum(covxy_rev)/(s_i - 1);
            
            % increment window
            ind = ind + s_i;
        end
        fsx(i) = sum([f2sx; f2sx_rev])/(2*nwins);
        fsy(i) = sum([f2sy; f2sy_rev])/(2*nwins);
        fsxy(i) = sum([f2xy; f2xy_rev])/(2*nwins);
        
    end
    
    % compute the scale-wise correlation coefficient
    betaxy = fsxy./fsx;
    fsx = sqrt(fsx);
    fsy = sqrt(fsy);
    
end

end