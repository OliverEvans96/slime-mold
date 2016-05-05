% fpivot.m 
% JLS April 17, 2002, revised April 2016
% A function that performs a pivot move on a chain in two dimensions
% Input are the pivot site of the chain (npivot) and the type of move (rpivot)
% If the move is successful, the function returns the coordinates of the beads 
% of the chain after the move, if the move is not successful, it returns the 
% original coordinates

function [xnew,ynew,rflag] = fpivot(npivot,rpivot,x,y)

global N % chain length 

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
overlap = 0; % overlap flag, set to 1 when beads overlap
k = 1; % initialize the counter for the outer loop
while( k < N && overlap==0 )
    m=k+1; % initialize the counter for the inner loop
    while( m <= N && overlap==0 )
        % calculate the distance between beads
        distsq = ( xtemp(m)-xtemp(k) )^2 + ( ytemp(m)-ytemp(k) )^2;
        if distsq < tiny    % check if two beads overlap
            overlap = 1;
        end             
        m = m + 1; % augment the counter for the inner loop
    end
    k = k + 1;  % augment the counter for the outer loop
end   

rflag = 1; % rejection flag, initialize to 1 = rejected  
if overlap == 0   % move was successful, return new coordinates
    xnew = xtemp;
    ynew = ytemp;
    rflag = 0; % move was not rejected
else % move was unsuccessful, return old coordinates
    xnew = x;
    ynew = y;
end