% fitlinerrors.m
% Fit a line to data points
% Accepts x, y and sigma (y-error) as arguments
% Returns a (intercept), b (slope) 
% and errors (variance = sigma^2) for both
% Also return r^2 value
% (coefficient of determination)

% From Numerical Recipes in C
% 15.2 Fitting Data to a Straight Line

function [a,b,var_a,var_b,chi2,Q] = linfiterrors(x,y,sigma)

% Basic statistical quantities
S  = sum(1./sigma.^2);
Sx = sum(x./sigma.^2);
Sy = sum(y./sigma.^2);

% Rescaled model parameters
t = (x-Sx/S)./sigma;
Stt = sum(t.^2);

% Fit solution
b = 1/Stt*sum(t.*y./sigma);
a = ((Sy-Sx*b)/S);

% Fit errors
var_a = 1./S*(1+Sx.^2./(S.*Stt));
var_b = 1./Stt;

% Chi Squared
chi2 = sum(((y-a-b*x)./sigma).^2);

% Goodness-of-fit
N = length(x);
Q = gammainc((N-2)/2,chi2/2);

end