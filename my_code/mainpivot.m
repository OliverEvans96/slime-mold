% JLS, April 17, 2001, revised 2002, 2016
% mainpivot.m
% This program uses the pivot algorithm to generate a self-avoiding random walk
% from an initially straight chain conformation.
% mainpivot calls the function fpivot which attempts a pivot move and returns
% the new coordinates, if the move was successful (rflag=0), and 
% the original coordinates if the move was not successful (rflag=1). 

function [accept_rate,avgResq,stdResq] = mainpivot(N,beta,MCsteps,draw)

%clear all; % clear all variables
%close all; % close all figures

%global N % make variable globally visible

% set flags for animation and evaluation
illustration = 0; % illustration flag, set to zero to switch off the display of
% intermediate configurations
evaluation = 1; % evaluation flag, set to one to perform block averaging 

%N = 21; % length of chain
%MCsteps = 10000; % total number of Monte Carlo steps 
MCequilib = MCsteps/10; % number of MC steps for equilibration
accept_rate = 0; % acceptance rate, ratio of attempted to successful moves
current_overlap = 0; % number of overlaps this time step
naccept = 0; % counter for accepted steps
Resq = zeros(1,MCsteps); % end-to-end vector squared
x = zeros(1,N); % x coordinates of chain sites
y = zeros(1,N); % y coordinates of chain sites
xtemp = zeros(1,N); %temporary x-value
ytemp = zeros(1,N); %temporary y-value
npivot = 0; % variable for the pivot site on the chain
rpivot = 0; % variable for the type of pivot move (see fpivot.m)
nblock = 10; % number of blocks for evaluation, Mcsteps/nblock has to be an integer
mblock = MCsteps/nblock; % number of MC steps in a block
bini = 0; % index for first member of the block
bfin = 0; % index for last member of the block
bResq = zeros(1,nblock); % average value per block
avgResq = 0;  % average value of the end-to-end distance
stdResq = 0;  % standard deviation of the end to end distance

seed = sum(1000*clock); % generate a seed from the clock time
% seed = 7; % use a fixed seed for debugging
rand('state',seed); % seed Matlab's random number generator

%create the initial configuration
for i=1:N
    x(i)=i;
    y(i)=0;
end

% Whether to draw initial & final configuration
if (draw == 1)
    figure(11); clf
    plot(x,y,'ro-','MarkerFaceColor','r','LineWidth',2)
    xlabel('x')
    ylabel('y')
    title(['Initial configuration, N = ',num2str(N),', n = ',num2str(N-1)])
    axis ([-1 N+1 -N/2 N/2])  
    axis equal
end

% if (illustration == 1)
%     figure(11); clf
%end
rflag = 0; % rejection flag, 0 for a successful move, 1 otherwise
for m = 1:MCsteps+MCequilib   % loop over Monte Carlo steps
    for n = 1:N     % loop for a single Monte Carlo step
        %fprintf('m=%d\nn=%d\n\n',m,n)
        npivot = 1 + ceil(rand*(N-2)); % pick a pivot site
        rpivot = ceil(rand*3); % pick a pivot move
        % the function fpivot to performs a trial move and 
        % overwrites the x,y coordinates with the new coordinates if successful
        [x,y,rflag,overlap] = fpivot(N,beta,current_overlap,npivot,rpivot,x,y);
        
        % Update number of overlaps
        current_overlap = overlap;
        
        % display the current configuration if the illustration flag is set to 1
        if (rflag == 0)  
            naccept = naccept+1;  % augment the counter of accepted steps 
        end
%         if illustration == 1
%             plot(x,y,'ro-');   
%             axis equal
%             axis ([-1 N+1 -N/2 N/2])
%             pause(.1)
%         end
    end
    if m == MCequilib
        accept_rate = naccept/(MCequilib*N);  % display the acceptance rate
        naccept = 0; % reset the counter to zero
    elseif m > MCequilib
        % once the system is equilibrated, calculate the end-to-end distance 
        % after each Monte Carlo step 
        Resq(m-MCequilib) = (x(N)-x(1))^2+(y(N)-y(1))^2; % end-to-end vector squared
    end
end

if evaluation == 1
    accept_rate = naccept/(N*MCsteps); % calculate the acceptance rate
    % perform block averaging to determine the average value of the end-to-end
    % distance and its error
    for nb=1:nblock
        bini=(nb-1)*mblock+1; % index for first member of the block
        bfin=nb*mblock; % index for last member of the block
        bResq(nb)=mean(Resq(bini:bfin));
    end
    avgResq = mean(bResq);
    stdResq = std(bResq);

    % display results for the end-to-end distance
%     figure(2); clf
%     plot(1:MCsteps,Resq,'k-')
%     xlabel('MC steps')
%     ylabel('Square end-to-end distance')
%     hold on
%     plot(mblock/2:mblock:MCsteps,bResq,'rs-','MarkerFaceColor','r','LineWidth',2)
%     title({['n = ',num2str(N-1),', seed = ',num2str(seed),...
%         ', MCsteps = ',num2str(MCsteps,5),', MCequilib = ',num2str(MCequilib,5)],...
%         ['acceptance rate = ',num2str(100*accept_rate,3),'% ,  R_e^2 = ',num2str(avgResq,4),...
%         ', \sigma_{Re^2} = ',num2str(stdResq,2)]});
%     legend('raw data','block averages')
%     hold off
end

if(draw)
    figure(11); 
    title(['Initial and final configuration, N = ',num2str(N),', n = ',num2str(N-1)])
    hold on
    plot(x,y,'bd-','MarkerFaceColor','b','LineWidth',2) % display the final configuration
    axis ([-1 N+1 -N/2 N/2])  
    axis equal
    legend('initial','final')
end

end
