addpath(genpath('./'))

x = linspace(0,1,10);
y = sin(x*3) + randn(1,length(x))*0.1;

figure(1)
plot(x,sin(x*3))
title('Noisy samples from function sin(x*3)')

%__________________________________________________________________________
%Example 1: 
%Pick a length scale of 0.01, and variance of 1 (noise 0.1)
hyp.cov = log([0.01,1]);
hyp.lik = log(0.1);

%Now fit a Gaussian process to the data
xprime = linspace(0,1,100);
[ymu ys2 fmu fs2   ] = gp(hyp, 'infExact', [], 'covSEiso', 'likGauss', x', y', xprime');

[L(1,1) dnlZ   ] = gp(hyp, 'infExact', [], 'covSEiso', 'likGauss', x', y');

figure(2)
plot(x,sin(x*3),'o'),hold on,plot(xprime,fmu,'r-')


%__________________________________________________________________________
%Example 2:
%Pick a length scale of 3, and variance of 1 (noise 0.1)
hyp.cov = log([3,1]);
hyp.lik = log(0.1);
[ymu ys2 fmu fs2   ] = gp(hyp, 'infExact', [], 'covSEiso', 'likGauss', x', y', xprime');
figure(3)
plot(x,sin(x*3),'o'),hold on,plot(xprime,fmu,'r-')
[L(2,1) dnlZ   ] = gp(hyp, 'infExact', [], 'covSEiso', 'likGauss', x', y');
%Now fit a Gaussian process

%Example 3
%__________________________________________________________________________
%Now choose the parameter by optimising the marginal likelihood with
%respect to the hyperparameters
MLhyp  = minimise(hyp, @gp, -100, 'infExact', [], 'covSEiso', 'likGauss', x', y');

[L(3,1) dnlZ   ] = gp(MLhyp, 'infExact', [], 'covSEiso', 'likGauss', x', y');
%Now fit a Gaussian process
[ymu ys2 fmu fs2   ] = gp(MLhyp, 'infExact', [], 'covSEiso', 'likGauss', x', y', xprime');
figure(4)
plot(x,sin(x*3),'o'),hold on,plot(xprime,fmu,'r-')


rmpath(genpath('./'))