%testlinfiterrors.m

x = 1:10;
y = 2*x+3 + 2*rand(1,10)-1;
sigma = 0*x + 1;

[a,b,var_a,var_b,chi2] = linfiterrors(x,y,sigma);

figure(1); clf
hold on
errorbar(x,y,sigma,'bo')
plot(x,b*x+a,'r-');