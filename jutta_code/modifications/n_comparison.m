% Oliver Evans
% Thermal Physics
% 5-4-2016
% Run mainpivot for a variety of N values
% Plot accept_rate, avgResq, stdResq as functions of N

% Clear all variables
clear all
% Close all figures
% close all

% A list of some colors
colors = ['b','g','r','m','c'];
n_colors = length(colors);

% Energy penalty for overlap
penalty_vals = linspace(0.5,1,2);
beta_vals = -log(penalty_vals(1:1)); % SKIP BETA=INF: ALREADY DONE

% Save plot handles from log plot to remove some labels
plot_handles = zeros(1,length(beta_vals));
log_plot_handles = zeros(1,length(beta_vals));

% Set simulation parameters
k_max = 1;
N_vals = 11:10:k_max*10+1;
N_vals = 31:31;
MCsteps = 2000;

seed = sum(1000*clock); % generate a seed from the clock time
% seed = 7; % use a fixed seed for debugging
rand('state',seed); % seed Matlab's random number generator

% Only run one long chain
% beta_vals = 1.61;
% N_vals = 121;

beta_num = 0;
figure(1); clf;
figure(2); clf;
figure(3); clf;

for beta = beta_vals
    fprintf('\n---------')
    % Counter while looping through energy penalties
    beta_num = beta_num + 1;
    color_num = mod(beta_num-1,n_colors)+1;
    fprintf('\n')
    fprintf('beta = %.2f\n',beta)
    fprintf('beta_num = %d\n',beta_num)

    % Allocate results arrays
    accept_rate = zeros(1,length(N_vals));
    avgResq = zeros(1,length(N_vals));
    stdResq = zeros(1,length(N_vals));

    k=1; % counter
    for N = N_vals
        fprintf('\nN = %d\n',N)
        % Run simulation
        [accept_rate(k),avgResq(k),stdResq(k)] = mainpivot(N,beta,MCsteps,0,seed);

%         fprintf('k = %3d: ',k)
%         fprintf('%5.2f ',[accept_rate(k),avgResq(k),stdResq(k)])
%         fprintf('\n')

        % Increment counter
        k = k + 1;
    end

    % Calculate line of best fit for double log
    [m_line,b_line] = linleastsq(log(N_vals-1),log(avgResq));
    x_line = [log(N_vals(1)-1),log(N_vals(end)-1)];
    y_line = m_line * x_line + b_line;
    fprintf('m = %.2f \nb = %.2f\n',m_line,b_line)

    % Create figure
    figure(1)

    % Plot results
    plot_handles(beta_num) = errorbar(N_vals-1,avgResq,stdResq,'o-',...
        'color',colors(color_num),...
        'DisplayName',sprintf('\\beta=%.2f: m=%.2f',beta,m_line));
    hold on
    title('Regular plot')
    ylabel('\langle R_e^2 \rangle')
    xlabel('n')
    drawnow

    % Log plot
    figure(2)
    hold on
    log_plot_handles(beta_num) = plot(log(N_vals-1),log(avgResq),'o',...
        'color',colors(color_num),...
        'DisplayName',sprintf('\\beta=%.2f: m=%.2f',beta,m_line));
    plot(x_line,y_line,'--','color',colors(color_num))
    title('Log Plot')
    ylabel('ln \langle R_e^2 \rangle')
    xlabel('ln n')
    drawnow
    
    % Slope plot
    figure(3)
    hold on
    plot(exp(-beta),m_line,'o-','color',colors(color_num))
    title('Slope Plot')
    xlabel('e^{-\beta}')
    ylabel('slope (m\_line)')
    drawnow
end

% Show legends
figure(1);
legend(plot_handles,'Location','northwest')

figure(2);
log_legend = legend(log_plot_handles,'Location','northwest');

figure(1);
legend(plot_handles,'Location','northwest')

% Save plots
saveas(1,'regular.png')
saveas(2,'log.png')
saveas(3,'slope.png')

% Save data
filename = sprintf('n_comparison.mat',N,beta);
fprintf('dir: %s\n',pwd)
fprintf('fname: n_comparison.mat\n',N,beta);
save(filename)
