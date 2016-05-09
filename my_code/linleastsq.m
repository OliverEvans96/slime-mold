% Oliver Evans
% 5-4-2016

% Calculate the slope (m) and intercept (b)
% of the line of best fit for a set of data
% based a 2D least squares approximation.

% Also calculate error for m and b
% based on error in data


function [m,b] = linleastsq(x,y)
    % Construct coefficient matrix
    A = zeros(2);
    A(1,1) = sum(x.^2);
    A(1,2) = sum(x);
    A(2,1) = sum(x);
    A(2,2) = length(x);
    
    % Right hand side
    b = zeros(2,1);
    b(1) = sum(x.*y);
    b(2) = sum(y);
    
    % Solve matrix equation
    sol = A\b;
    
    % Return results
    m = sol(1);
    b = sol(2);
end
