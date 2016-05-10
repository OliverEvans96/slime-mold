% Oliver Evans
% Thermal Physics
% 5-4-2016
% Load data generated by n_comparison
% Plot accept_rate, avgResq, stdResq as functions of N

% Clear all variables
clear all
% Close all figures
% close all

% Save data
filename = sprintf('../data/n_comparison.mat');
fprintf('n_comparison.m\n')
fprintf('dir: %s\n',pwd)
fprintf('fname: ../data/n_comparison.mat\n');
load(filename)

% A list of some colors
%colors = ['b','g','r','m','c'];
%n_colors = length(colors);

% Energy penalty for overlap
%penalty_vals = linspace(0,2,10);
%beta_vals = -log(penalty_vals.^2);
%beta_vals = penalty_vals.^2;

% Set simulation parameters
%k_max = 4;
%N_vals = (11:10:k_max*10+1)'; % Evenly spaced in linear scale
%N_vals = (2.^(4:7))'      % Evenly spaced in log scale
%MCsteps = 10000

%seed = sum(1000*clock); % generate a seed from the clock time
%seed = 7; % use a fixed seed for debugging
%rand('state',seed); % seed Matlab's random number generator

% Only run one long chain
% beta_vals = 1.61;
% N_vals = 121;

beta_num = 0;
figure(1); clf;
hold on
figure(2); clf;
hold on
figure(3); clf;
hold on
figure(4); clf;
hold on
figure(5); clf;
hold on

% % Allocate results arrays
% accept_rate = zeros(length(N_vals),length(beta_vals));
% avgResq = zeros(length(N_vals),length(beta_vals));
% stdResq = zeros(length(N_vals),length(beta_vals));
% avgoverlap = zeros(length(N_vals),length(beta_vals));
% stdoverlap = zeros(length(N_vals),length(beta_vals));
% avgCV = zeros(length(N_vals),length(beta_vals));
% stdCV = zeros(length(N_vals),length(beta_vals));
% logerrResq = zeros(length(N_vals),length(beta_vals));
% 
% % Allocate line fit arrays
% x_line = zeros(2,length(beta_vals));
% y_line = zeros(2,length(beta_vals));
% b_line = zeros(1,length(beta_vals));
% m_line = zeros(1,length(beta_vals));
% b_var = zeros(1,length(beta_vals));
% m_var = zeros(1,length(beta_vals));
% chi2 = zeros(1,length(beta_vals));
% Q = zeros(1,length(beta_vals));
% 
% % Save plot handles
% Resq_plot_handles =  zeros(1,length(beta_vals));
% log_plot_handles = zeros(1,length(beta_vals));
% slope_plot_handles = zeros(1,length(beta_vals));
% slope_error_bars = zeros(3,length(beta_vals));
% CV_plot_handles =  zeros(1,length(length(N_vals)));
% overlap_plot_handles =  zeros(1,length(length(N_vals)));

% m_std from m_var
m_std = sqrt(m_var);

% Loop through beta values (beta = 1/kT)
for beta = beta_vals
    fprintf('\n---------')
    % Counter while looping through energy penalties
    beta_num = beta_num + 1;
    color_num = mod(beta_num-1,n_colors)+1;
    fprintf('\n')
    fprintf('beta = %.2f\n',beta)
    fprintf('beta_num = %d\n',beta_num)
    fprintf('m = %.2f \nb = %.2f\n',m_line(beta_num),b_line(beta_num))

	% Resq Plot results
	figure(1)
	Resq_plot_handles(beta_num) = errorbar(N_vals-1,avgResq(:,beta_num),...
		stdResq(:,beta_num),'o-','color',colors(color_num),...
		'DisplayName',sprintf('\\beta=%.2f',beta));
	title(['R_e^2 plot','   (seed=',num2str(seed),')'])
	ylabel('\langle R_e^2 \rangle')
	xlabel('n')
	drawnow

	% Log plot
	figure(2)
	log_plot_handles(beta_num) = errorbar(log(N_vals-1),log(avgResq(:,beta_num)),...
		logerrResq(:,beta_num),'o','color',colors(color_num),...
		'DisplayName',sprintf('\\beta=%.2f: m=%.2f',beta,m_line(beta_num)));
	plot(x_line(:,beta_num),y_line(:,beta_num),'--','color',colors(color_num))
	title(['Log Plot','   (seed=',num2str(seed),')'])
	ylabel('ln \langle R_e^2 \rangle')
	xlabel('ln n')
	drawnow

	% Slope plot
	figure(3)
    % Slope plotting variables
    spx = beta;
%     spy = log(m_line(beta_num));
%     spe = m_std(beta_num)/m_line(beta_num);
    spy = m_line(beta_num);
    spe = m_std(beta_num);
    % Axis labels
	xlabel('\beta')
	ylabel('slope (m\_line)')
    
    % Draw plot w/ error bars
	%slope_plot_handles(beta_num) = errorbar(spx,spy,spe,...
	%	'o-','color',colors(color_num),...
	%	'DisplayName',sprintf('Q=%.2f',Q(beta_num)));
    title(['Slope Plot','   (seed=',num2str(seed),')'])

	slope_error_bars(:,beta_num) = terrorbar(spx,spy,spe,0.2);
	slope_plot_handles(beta_num) = plot(spx,spy,...
		'o-','color',colors(color_num),...
		'DisplayName',sprintf('Q=%.2f',Q(beta_num)));
	set(slope_error_bars(:,beta_num),'Color',colors(color_num))
	drawnow

	% Overlap Plot results
	figure(5)
	overlap_plot_handles(beta_num) = errorbar(N_vals-1,avgoverlap(:,beta_num)./(N-1),...
		stdoverlap(:,beta_num)/(N-1),'o-','color',colors(color_num),...
		'DisplayName',sprintf('\\beta=%.2f',beta));
	title(['Overlap plot','   (seed=',num2str(seed),')'])
	xlabel('n')
	ylabel('Number of overlaps')
	drawnow
end

% CV Plot results
figure(4)
k=1;
color_num=1;
for N = N_vals'
 	CV_plot_handles(k) = errorbar(beta_vals,avgCV(k,:)/(N-1),...
		stdCV(k,:)/(N-1),'o-','color',colors(color_num),...
		'DisplayName',sprintf('N=%d',N));
	title(['CV plot','   (seed=',num2str(seed),')'])
	xlabel('\beta')
	ylabel('$\displaystyle \langle \frac{C_V}{n} \rangle$','Interpreter','latex')
	%ylabel('C_V/n')
	drawnow
	k = k + 1;
    color_num = mod(k-1,n_colors)+1;
end

% Linear interpolation on slope plot
figure(3)
plot(beta_vals,m_line,'k--','DisplayName','lin. interp.')

% Add legends
legend(Resq_plot_handles,'Location','northwest')
legend(log_plot_handles,'Location','northwest')
legend(slope_plot_handles,'Location','southeast')
legend(CV_plot_handles,'Location','northwest')
legend(overlap_plot_handles,'Location','northwest')

% Save plots
saveas(1,'resq.png')
saveas(2,'log.png')
saveas(3,'slope.png')
saveas(4,'CV.png')
saveas(4,'overlap.png')

% % Save workspace
% save(filename)
