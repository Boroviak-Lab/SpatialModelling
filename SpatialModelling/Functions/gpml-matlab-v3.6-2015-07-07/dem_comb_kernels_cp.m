
%Generate 3 lots of test data. These are GP with sum of SE covariance
%function, product of SE covariance functions and a single SE covariance
%function.
x1   = linspace(0,60,300);
K1   = (covSEiso(log([.5;5]),x1') + covSEiso(log([5;5]),x1') ) + covNoise(log(0.5),x1');
K2   = (covSEiso(log([.5;5]),x1').*covSEiso(log([5;5]),x1')) + covNoise(log(0.5),x1');
K3   = covSEiso(log([5;5]),x1') + covNoise(log(0.5),x1');
y{1} = real(gsamp(zeros(1,300),K1,10)); %1 = sum
y{2} = real(gsamp(zeros(1,300),K2,10)); %2 = prod
y{3} = real(gsamp(zeros(1,300),K3,10)); %3 = ind
for j = 1:3
for i = 1:10
y{j}(i,:) = (y{j}(i,:)-mean(y{j}(i,:)))./std((y{j}(i,:)-mean(y{j}(i,:))));
end
end

%Plot these functions to see that they're fairly different by eye
subplot(3,2,1);plot(real(gsamp(zeros(1,300),K3,10))')
subplot(3,2,3);plot(real(gsamp(zeros(1,300),K3,10))')
subplot(3,2,5);plot(real(gsamp(zeros(1,300),K3,10))')

Modell = {'Sum','Prod','Ind'};
plotind = [2,4,6];

%Single noise model
likfunc = @likGauss; sn = 0.1; hyp.lik = log(sn);

%Loop over the three datasets and try to fit one of 3 models to each one
for inds   = 1:3
H1 = zeros(5,10);
H2 = zeros(5,10);
H3 = zeros(3,10);

for i = 1:10

%Model is sum of base kernels   
covfunc = {@covSum, {@covSEiso, @covSEiso}};
hyp2.cov = log([0.5; 1; 4; 2]); 
hyp2.lik = log(0.1);
hyp2 = minimize(hyp2, @gp, -100, @infExact, [], covfunc, likfunc, x1', y{inds}(i,:)');
H1(:,i) = [hyp2.cov;hyp2.lik];
nlml(1,i) = -2*gp(hyp2, @infExact, [], covfunc, likfunc, x1', y{inds}(i,:)') + (length(hyp2.cov)+1)*log(length(x1));

%Product of base kernels
covfunc = {@covProd, {@covSEiso, @covSEiso}};
hyp2.cov = log([0.5; 1; 4; 2]); 
hyp2.lik = log(0.1);
hyp2 = minimize(hyp2, @gp, -100, @infExact, [], covfunc, likfunc, x1', y{inds}(i,:)');
H2(:,i) = [hyp2.cov;hyp2.lik];
nlml(2,i) = -2*gp(hyp2, @infExact, [], covfunc, likfunc, x1', y{inds}(i,:)') + (length(hyp2.cov)+1)*log(length(x1));

%Single kernel
covfunc = @covSEiso; 
hyp2.cov = log([5; 1]); 
hyp2.lik = log(0.1);
hyp2 = minimize(hyp2, @gp, -100, @infExact, [], covfunc, likfunc, x1', y{inds}(i,:)');
H3(:,i) = [hyp2.cov;hyp2.lik];
nlml(3,i) = -2*gp(hyp2, @infExact, [], covfunc, likfunc, x1', y{inds}(i,:)') + (length(hyp2.cov)+1)*log(length(x1));
end

%Now plot the inferred models
subplot(3,2,plotind(inds)); plot(mean(nlml,2),'o'), title(['Best model is ' Modell{(find(sum(nlml,2)==max(sum(nlml,2))))} ', True model is ' Modell{inds}])

end