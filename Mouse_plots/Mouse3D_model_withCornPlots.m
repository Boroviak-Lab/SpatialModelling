h4=figure;

addpath(genpath('/Users/christopherpenfold/Desktop/Code/BranchingGPs/code/gpml-matlab-v3.6-2015-07-07'))
%Load in the mouse legends (and any other legenes if we have more datasets)
Leg = importdata('/Users/christopherpenfold/Desktop/Thorsten/CombinedModelling/mouse/mouseKey.csv');
LegSC = importdata('/Users/christopherpenfold/Desktop/Thorsten/CombinedModelling/mouse/mousescKey.csv');
LegSp = importdata('/Users/christopherpenfold/Desktop/Thorsten/CombinedModelling/mouse/mousespatialKey.csv');

Dexp = importdata('/Users/christopherpenfold/Desktop/Thorsten/CombinedModelling/mouse/featurecountsMouseSpatialAll.csv')
Dexp2 = importdata('/Users/christopherpenfold/Desktop/Thorsten/CombinedModelling/mouse/featurecountsMouseSCAll.csv')

Ind = find(LegSp.data(:,4)==1);

genes = Dexp.textdata(2:end,1);
gene2test2 = {'WNT3','EOMES','FOXA2','OTX2','POU5F1','CER1','SOX2','MIXL1'}
gene2test2 = {'OTX2','EOMES','FOXA2','MIXL1','POU5F1','CER1','SOX2'}

for i = 1:length(gene2test2)
    inds(i,1) = find(strcmp(genes, gene2test2{i} ))
end

Dexp1 = Dexp.data(:,Ind);
Dexp1 = Dexp1(inds,:);

Xtrain = [LegSp.data(:,5),LegSp.data(:,6),1.5*LegSp.data(:,7)];
Xtrain = Xtrain(Ind,:);
Xtrain = Xtrain + randn(size(Xtrain,1),3)*0.1;

fileids = LegSp.textdata(2:end,1);
fileids = fileids(Ind);
cellids = LegSp.textdata(2:end,5);
cellids = cellids(Ind);

ucellids = unique(cellids);
Dcorn = zeros(size(Dexp1,1),length(ucellids));
for i = 1:length(ucellids)
    ind1 = find(strcmp(cellids,ucellids{i}));    
    Dcorn(:,i) = mean(Dexp1(:,ind1),2);    
end
    
%Load and plot the cornplots
%load('MyColormaps.mat','mycmapnew');
for i=1:length(ucellids)
    tmp=ucellids{i};
    pos=strfind(tmp,'.');
    nsample(i)={tmp(1:pos(1)-1)};
    stage(i)={tmp(pos(1)+1:end)};
end

sv(1).sample=nsample;
sv(1).value=log10(Dcorn+1);

%mycmapnew = colormap(default)
addpath(genpath('/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/RCode/CornPlot/'))
 
 for i = 1:9
 subplot(9,5,(i-1)*5+1);
 %E70EpiCornPlot(sv(1).value(i,:),sv(1).sample,'E7.0',mycmapnew,0,max(sv(1).value(i,:)))
 %caxis([0,max(sv(1).value(i,:))] )
 
 E70EpiCornPlot(sv(1).value(i,:),sv(1).sample,'E7.0',1,0,4)
 caxis([0,4] )
 
 title(gene2test2{i})
 #pause
 clf
 end


%Generate interpolation test locations
% 
% [x,y,z] = sph2cart(2*pi*rand(2000,1), asin(2*rand(2000,1)-1),(rand(2000,1).^(1/3)));
% z = z*2.5;
% x = x(find(z<0));
% y = y(find(z<0));
% z = z(find(z<0));
% r22 = x.^2 + y.^2 + z.^2;
% 
% x = x(find(r22>0.5^2));
% y = y(find(r22>0.5^2));
% z =  z(find(r22>0.5^2));
% subplot(1,3,1); plot3(x,y,z,'.'),xlim([-3 3]),ylim([-3 3]),zlim([-3 3])
% 
[x,y,z] = sph2cart(2*pi*rand(200000,1), asin(2*rand(200000,1)-1),(rand(200000,1).^(1/3)));
z = z;
x = 1.2*x(find(z<0));
y = 1.2*y(find(z<0));
z = 1.2*z(find(z<0));
r22 = x.^2 + y.^2 + z.^2;

x = x(find(r22>0.75^2));
y = y(find(r22>0.75^2));
z =  2.5*z(find(r22>0.75^2));
%subplot(1,3,2); 
%plot3(x,y,z,'.'),xlim([-3 3]),ylim([-3 3]),zlim([-3 3])
%hold on
%plot3(LegSp.data(:,5),LegSp.data(:,6),LegSp.data(:,7),'rs')

[xsph_test,ysph_test,zsph_test] = cart2sph(Xtrain(:,1), Xtrain(:,2),Xtrain(:,3));
[xsph,ysph,zsph] = cart2sph(x,y,z);
Xtest = [x,y,z];

Xtrain2 = [xsph_test,ysph_test,zsph_test];
Xtest2 = [xsph,ysph,zsph];

    
Ytrain = log10(Dexp1+1);
for i = 1:size(Ytrain,1)
    Mnorm(i,1) = mean(Ytrain(i,:));
    Stdnnorm(i,1) = std((Ytrain(i,:)-mean(Ytrain(i,:))));
    Ytrain(i,:) = (Ytrain(i,:)-mean(Ytrain(i,:)))./std((Ytrain(i,:)-mean(Ytrain(i,:))) ); 
end

%Let us first set the priors
pgls = {@priorGauss,0,1}; 
pgvar = {@priorGauss,0,1}; 
pglik = {@priorGauss,log(1/3),1}; 
meanc = {@meanConst}; cov = {@covSEard}; lik = {@likGauss};

for i = 1:size(Ytrain,1)

prior.cov = {pgls;pgls;pgls;pgls};  % Gaussian prior for first, clamp second par
prior.mean  = {[]}; % box prior for first, nothing for second par
prior.lik = {pglik};

im = {@infPrior,@infExact,prior};
par = {meanc,cov,lik,Xtrain(:,1:3), Ytrain(i,:)'}; 
mfun = @minimize; % input for GP function

hyp.cov  = [std(Ytrain(i,:)); log(0.1); log(0.1); log(0.1)];
hyp.lik  = [std(Ytrain(i,:))/3];    
hyp.mean = mean(Ytrain(i,:));       

im = {@infPrior,@infExact,prior};                % inference method
hyp2 = feval(mfun, hyp, @gp, -1000, im, par{:});         % optimise

%hyp2 = minimize(hyp, @gp, -1000, 'infExact', @meanConst, @covSEard, 'likGauss', Xtrain(:,1:3), Ytrain(i,:)');
[M S] = gp(hyp2, 'infExact', @meanConst, @covSEard, 'likGauss',  Xtrain(:,1:3), Ytrain(i,:)', Xtest(:,1:3));

m1{i} = M;
s1{i} = S;
H3{i} = hyp2;
end
    
jit = randn(size(Xtrain,1),size(Xtrain,2))*0.05;

for i = 1:size(Ytrain,1)
subplot(9,5,(i-1)*5+2);
scatter3(Xtrain(:,1)+jit(:,1),Xtrain(:,2)+jit(:,2),Xtrain(:,3)+jit(:,3),25,Ytrain(i,:)*Stdnnorm(i,1) + Mnorm(i,1),'fill') 
%set(gca,'XTick',[],'YTick',[],'ZTick',[])
set(gca,'visible','off')
%caxis([0,max(sv(1).value(i,:))] )
caxis([0,4] )
view([-158 29])

subplot(9,5,(i-1)*5+3);
scatter3(Xtest(:,1),Xtest(:,2),Xtest(:,3),25,m1{i}*Stdnnorm(i,1) + Mnorm(i,1),'fill') 
%set(gca,'XTick',[],'YTick',[],'ZTick',[])
set(gca,'visible','off')
caxis([0,4] )
view([-1 1 1])
view([-158 29])
%view([-180 0])
%view([-180 20])
%view([-180 -40])

subplot(9,5,(i-1)*5+4);
scatter3(Xtest(:,1),Xtest(:,2),Xtest(:,3),25,m1{i}*Stdnnorm(i,1) + Mnorm(i,1),'fill') 
%set(gca,'XTick',[],'YTick',[],'ZTick',[])
set(gca,'visible','off')
caxis([0,4] )
view([-158 0])

subplot(9,5,(i-1)*5+5);
scatter3(Xtest(:,1),Xtest(:,2),Xtest(:,3),25,m1{i}*Stdnnorm(i,1) + Mnorm(i,1),'fill') 
%set(gca,'XTick',[],'YTick',[],'ZTick',[])
set(gca,'visible','off')
caxis([0,4] )
view([-180 -40])


end

set(h4,'PaperSize',[12 4*12 ]);
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6*2 12*2];
%set(h4,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

print('./TestSC','-dpdf','-r0','-painters')

return
