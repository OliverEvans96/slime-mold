% fpivot.m 
% JLS April 17, 2002, revised April 2016
% A function that performs a pivot move on a chain in two dimensions
% Input are the pivot site of the chain (npivot) and the type of move (rpivot)
% If the move is successful, the function returns the coordinates of the beads 
% of the chain after the move, if the move is not successful, it returns the 
% original coordinates

function [xnew,ynew,rflag,overlap] = fpivot(N,beta,current_overlap,npivot,rpivot,x,y)

%global N % chain length 

tiny = 1.0e-4; % small number for distance comparison

xtemp = x; % initialize the trial x-coordinates to the original coordinates
ytemp = y; % initialize the trial y-coordinates to the original coordinates 
distsq = 0; % initialize variable for distance squared between beads

if rpivot == 1 % rotation through 90 degree counterclockwise
    if npivot < N/2  % pivot site on first half of chain
        for kk = 1:npivot-1  
            xtemp(kk) = x(npivot)+y(kk)-y(npivot);
            ytemp(kk) = y(npivot)-x(kk)+x(npivot);
        end
    else
        for kk = npivot+1:N
            xtemp(kk) = x(npivot)+y(kk)-y(npivot);
            ytemp(kk) = y(npivot)-x(kk)+x(npivot);
        end
    end
elseif rpivot == 2 % rotation through 90 degree clockwise
    if npivot < N/2  
        for kk = 1:npivot-1
            xtemp(kk) = x(npivot)-y(kk)+y(npivot);
            ytemp(kk) = y(npivot)+x(kk)-x(npivot);
        end
    else
        for kk = npivot+1:N
            xtemp(kk) = x(npivot)-y(kk)+y(npivot);
            ytemp(kk) = y(npivot)+x(kk)-x(npivot);
        end
    end
else % rotation through 180 degree
    if npivot < N/2  
        for kk = 1:npivot-1
            xtemp(kk) = 2*x(npivot)-x(kk);
            ytemp(kk) = 2*y(npivot)-y(kk);
        end
    else
        for kk = npivot+1:N
            xtemp(kk) = 2*x(npivot)-x(kk);
            ytemp(kk) = 2*y(npivot)-y(kk);
        end
    end
end

% check for overlap by calculating the distance between 
% all pairs of non-bonded beads
overlap = 0; % overlap counter, add 1 for each overlap
% initialize the counter for the outer loop
for k=1:N-1
    % initialize the counter for the inner loop
    for m=(k+1):N
        % calculate the distance between beads
        distsq = ( xtemp(m)-xtemp(k) )^2 + ( ytemp(m)-ytemp(k) )^2;
        if distsq < tiny    % check if two beads overlap
            overlap = overlap + 1;
        end
    end
end   

rflag = 1; % rejection flag, initialize to 1 = rejected  

% Difference in number of overlaps between current and proposed states
del_overlap = overlap - current_overlap;

% Threshold for accepting new state
if (del_overlap <= 0)
    prob_threshold = 1;
else
    prob_threshold = exp(-beta*del_overlap);
end

% Generate random number
r = rand;

if (r < prob_threshold)  % move was successful, return new coordinates
    xnew = xtemp;
    ynew = ytemp;
    rflag = 0; % move was not rejected
%     disp('accept')
else % move was unsuccessful, return old coordinates, reset overlap counter
    xnew = x;
    ynew = y;
    overlap = 0;
%     disp('reject')
end

% fprintf('\n\n\n')

end
