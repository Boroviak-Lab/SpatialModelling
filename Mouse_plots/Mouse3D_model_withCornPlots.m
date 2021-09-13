h4=figure;

addpath(genpath('/Users/christopherpenfold/Desktop/Code/BranchingGPs/code/gpml-matlab-v3.6-2015-07-07'))
%Load in the mouse legends (and any other legenes if we have more datasets)
Leg = importdata('/Users/christopherpenfold/Desktop/Thorsten/CombinedModelling/mouse/mouseKey.csv');
LegSC = importdata('/Users/christopherpenfold/Desktop/Thorsten/CombinedModelling/mouse/mousescKey.csv');
LegSp = importdata('/Users/christopherpenfold/Desktop/Thorsten/CombinedModelling/mouse/mousespatialKey.csv');

%Load in the normalised data (and other data we may want to plot)
%Dexp = importdata('/Users/christopherpenfold/Desktop/Thorsten/FINAL/Mouse2_Modelling20k/NormData.csv')

Dexp = importdata('/Users/christopherpenfold/Desktop/Thorsten/CombinedModelling/mouse/featurecountsMouseSpatialAll.csv')
Dexp2 = importdata('/Users/christopherpenfold/Desktop/Thorsten/CombinedModelling/mouse/featurecountsMouseSCAll.csv')

%trainingset = LegSp.textdata(setdiff(1:length(LegSp.textdata(1:end,3)),strmatch('sc',LegSp.textdata(1:end,3))),1);
%trainingset = trainingset(2:end);

%headers = Dexp.textdata(1,2:end);
%headers = strrep(headers,'"','');


%view([-180 0])
%view([-180 20])
%view([-180 -40])

Ind = find(LegSp.data(:,4)==1);

%for i = 1:length(trainingset)
%   Ind(i,1) = strmatch(trainingset{i},headers);
%end

genes = Dexp.textdata(2:end,1);
%gene2test2 = {'OTX2','HESX1','FOXA2','NOTO','TBXT','MIXL1','EOMES','WNT5A','TDGF1'}

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
 pause
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


% hyp = [ log(ell_1)
%         log(ell_2)
%          .
%         log(ell_D)
%         log(sf) ]


%Let us first set the priors
pgls = {@priorGauss,0,1}; 
pgvar = {@priorGauss,0,1}; 
pglik = {@priorGauss,log(1/3),1}; 

%par = {mean,cov,lik,x,y}; mfun = @minimize; % input for GP function

% a) plain marginal likelihood optimisation (maximum likelihood)
%im = @infExact;                                  % inference method
%hyp_plain = feval(mfun, hyp, @gp, -10, im, par{:});      % optimise

% b) regularised optimisation (maximum a posteriori) with 1d priors
%prior.mean = {pg;pc};  % Gaussian prior for first, clamp second par
%prior.cov  = {p1;[]}; % box prior for first, nothing for second par
%im = {@infPrior,@infExact,prior};                % inference method
%hyp_p1 = feval(mfun, hyp, @gp, -10, im, par{:});         % optimise

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

%hyp3 = minimize(hyp, @gp, -1000, 'infExact', @meanConst, @covSEard, 'likGauss', Xtrain2(:,1:3), Ytrain(i,:)');
%[M1 S1] = gp(hyp3, 'infExact', @meanConst, @covSEard, 'likGauss',  Xtrain2(:,1:3), Ytrain(i,:)', Xtest2(:,1:3));
%
%m1{i} = M1;
%s1{i} = S1;
%H3{i} = hyp3;

%hyp.cov  = [var(Ytrain(i,:)); log(0.1); log(0.1); log(0.1); log(0.1); log(0.1); log(0.1)]';
%hyp.lik  = [var(Ytrain(i,:))/3];    
%hyp.mean = mean(Ytrain(i,:));       
%hyp4 = minimize(hyp, @gp, -100, 'infExact', @meanConst, @covSEard, 'likGauss', [Xtrain(:,1:3),Xtrain2(:,1:3)], Ytrain(i,:)');

%[M2 S2] = gp(hyp4, 'infExact', @meanConst, @covSEard, 'likGauss',  [Xtrain(:,1:3),Xtrain2(:,1:3)], Ytrain(i,:)', [Xtest(:,1:3),Xtest2(:,1:3)] );

%m2{i} = M2;
%s2{i} = S2;
%H4{i} = hyp4;

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


print('~/Desktop/TestSC','-dpdf','-r0','-painters')

return
set(gca,'position',[0.75 0.25 0.14 0.54]); 


h4=figure;

Dsc = importdata('/Users/christopherpenfold/Desktop/Thorsten/CombinedModelling/mouse/featurecountsMouseSCAll.csv');
cells = Dsc.textdata(1,2:end);
genes = Dsc.textdata(2:end,1);

%gene2test2 = {'OTX2','HESX1','FOXA2','NOTO','TBXT','MIXL1','EOMES','WNT5A','TDGF1'}
for i = 1:length(gene2test2)
    INSC(i,1) = find(strcmp(genes,gene2test2{i})==1)
end

%Now load the mouse embryo data

DPCA = importdata('/Users/christopherpenfold/Desktop/Thorsten/FINAL/Mouse_Modelling20k_4k_all4/DimRed/EmbeddingsPCA.csv');
PCAdist = DPCA.data(:,1:2);
Delta = dist(PCAdist');
PCAcells = DPCA.textdata(2:end,1);

isinsc = contains(PCAcells,'X3289STDY');

%contains(PCAcells,'SRR1803')
YSC = zeros(length(gene2test2),length(fileids));
for i = 1:length(fileids)
    %Where the cell matches in the PCA
    inds = find(strcmp(fileids{i},PCAcells)==1) 
    distances = Delta(inds,:);
    distances(1,inds) = Inf;
    distance = distances./isinsc';
    inds2 = find(distance==min(distance));
    ClosestNe = PCAcells{inds2}(2:end); %find(distance==min(distance))  
    
    indinSC = find(strcmp(cells,ClosestNe)==1);

    YSC(:,i) = Dsc.data(INSC,indinSC)
    
end


YSCtrain = log10(YSC+1);
for i = 1:size(YSCtrain,1)
    Mnorm1(i,1) = mean(YSCtrain(i,:));
    Stdnnorm1(i,1) = std((YSCtrain(i,:)-mean(YSCtrain(i,:))));
    YSCtrain(i,:) = (YSCtrain(i,:)-mean(YSCtrain(i,:)))./std((YSCtrain(i,:)-mean(YSCtrain(i,:))) ); 
end


for i = 1:size(YSCtrain,1)
hyp.cov  = [var(YSCtrain(i,:)); log(0.1); log(0.1); log(0.1)];
hyp.lik  = [var(YSCtrain(i,:))/3];    
hyp.mean = mean(YSCtrain(i,:));       
hyp5{i} = minimize(hyp, @gp, -1000, 'infExact', @meanConst, @covSEard, 'likGauss', Xtrain(:,1:3), YSCtrain(i,:)');
[m5{i} s5{i}] = gp(hyp5{i}, 'infExact', @meanConst, @covSEard, 'likGauss',  Xtrain(:,1:3), YSCtrain(i,:)', Xtest(:,1:3));
end


for i = [1:size(YSCtrain,1)]
subplot(9,5,(i-1)*5+2);
scatter3(Xtrain(:,1)+jit(:,1),Xtrain(:,2)+jit(:,2),Xtrain(:,3)+jit(:,3),25,YSCtrain(i,:)*Stdnnorm1(i,1) + Mnorm1(i,1),'fill') 
%set(gca,'XTick',[],'YTick',[],'ZTick',[])
set(gca,'visible','off')
xlim([-3 3])
ylim([-3 3])
zlim([-3 3])
%caxis([0,max(sv(1).value(i,:))] )
caxis([0,4] )
view([-1 1 1])
end


for i = [1:size(YSCtrain,1)]

subplot(9,5,(i-1)*5+5);
scatter3(Xtest(:,1),Xtest(:,2),Xtest(:,3),25,m5{i}*Stdnnorm1(i,1) + Mnorm1(i,1),'fill') 
%set(gca,'XTick',[],'YTick',[],'ZTick',[])
set(gca,'visible','off')
%caxis([0,4] )
xlim([-3 3])
ylim([-3 3])
zlim([-3 3])
view([-1 1 1])
end


%set(h4,'PaperSize',[12 4*12 ]);

%set(h4,'Units','Inches');
%pos = get(h4,'Position');

%set(h4,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6*2 12*2];
%set(h4,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])


print('~/Desktop/TestSC','-dpdf','-r0','-painters')




print('~/Desktop/TestSC','-svg','-painters')

return

set(gca,'position',[0.75 0.25 0.14 0.54]); 
















find(strcmp(cells,ClosestNe)==1)
P = linspace(0,2*pi,200);

xx = linspace(0,12*pi,200);
x = ones(1,200);
y = sin(xx)*2 -1;
z = ones(1,200);
% 
% for i = 1:10
% X = x(i); Y = y(i); Z = z(i);
%     subplot(2,2,1);
%     scatter3(Xtrain(:,1)+jit(:,1),Xtrain(:,2)+jit(:,2),Xtrain(:,3)+jit(:,3),25,Ytrain,'fill') 
%     set(gca,'XTick',[],'YTick',[],'ZTick',[])
%     view([X Y Z ])
%     subplot(2,2,2);
%     scatter3(Xtest(:,1),Xtest(:,2),Xtest(:,3),25,m,'fill')
%     set(gca,'XTick',[],'YTick',[],'ZTick',[])    
%     view([X Y Z ])    
%     subplot(2,2,3);
%     scatter3(Xtrain(:,1)+jit(:,1),Xtrain(:,2)+jit(:,2),Xtrain(:,3)+jit(:,3),25,Ytrain2,'fill') 
%     set(gca,'XTick',[],'YTick',[],'ZTick',[])
%     view([X Y Z ])
%     subplot(2,2,4);
%     scatter3(Xtest(:,1),Xtest(:,2),Xtest(:,3),25,m2,'fill') 
%     set(gca,'XTick',[],'YTick',[],'ZTick',[])    
%     view([X Y Z ])    
% 
%     pause
%     clf
% end

h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated.gif';
set(gcf,'color','w');

for n = 1:length(P)
    
    X = x(n); Y = y(n); Z = z(n);
    
    % Draw plot for y = x.^n
    subplot(2,2,1);
    scatter3(Xtrain(:,1)+jit(:,1),Xtrain(:,2)+jit(:,2),Xtrain(:,3)+jit(:,3),25,Ytrain,'fill') 
        set(gca,'XTick',[],'YTick',[],'ZTick',[])

    view([X Y Z ])
    grid off
    subplot(2,2,2);
    scatter3(Xtest(:,1),Xtest(:,2),Xtest(:,3),25,m,'fill') 
        set(gca,'XTick',[],'YTick',[],'ZTick',[])

    view([X Y Z ])  
    grid off
    subplot(2,2,3);
    scatter3(Xtrain(:,1)+jit(:,1),Xtrain(:,2)+jit(:,2),Xtrain(:,3)+jit(:,3),25,Ytrain2,'fill') 
        set(gca,'XTick',[],'YTick',[],'ZTick',[])

    view([X Y Z ])
    grid off
    subplot(2,2,4);
    scatter3(Xtest(:,1),Xtest(:,2),Xtest(:,3),25,m2,'fill') 
        set(gca,'XTick',[],'YTick',[],'ZTick',[])

    view([X Y Z ])    
    grid off
    drawnow 
      % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if n == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime',0.1); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime',0.1); 
      end 
  end


%hold on
%scatter3(Xtest(:,1),Xtest(:,2),Xtest(:,3),'rs','fill') 
%set(gca,'xtick',[],'ytick',[],'ztick',[])
%caxis manual; caxis([bottom top]);
%view([1 1 1 ])
%legend({'SHOT','Embryo'})
%print(['./Plots/' gene2test2{i} '_mouse.pdf'],'-dpdf')
%clf

 

%gene2test = 'SOX2'

%Load in the E15 and the key
%D1 = importdata('/Users/christopherpenfold/Desktop/Thorsten/shots_new.csv')
%D2 = importdata('/Users/christopherpenfold/Desktop/Thorsten/embryo_new.csv')
%D3 = importdata('/Users/christopherpenfold/Desktop/Thorsten/CombinedModelling/marmoset/KeyMATLAB.csv')

%Get look up the cell IDs. Need to do some jugglinig IDs around
%locations = D1.textdata(2:end,1)
%%for i = 1:length(locations)
%    try
%    idx(i,1) = find(strcmp(D3.textdata(:,4) , ['P1-E15B-',strrep(locations{i},'_','-')]));   
%    catch
%    idx(i,1) = find(strcmp(D3.textdata(:,4) , ['P2-E15B-',strrep(locations{i},'_','-')]));    
%    end    
%end
%CellIDs = D3.textdata(idx,1);
%D3subs = D3.data(idx,:);

%celltype = D3.datsa

%Load in expression data
%Dexp = importdata('/Users/christopherpenfold/Desktop/Thorsten/FINAL/Marmoset_Cyno_Human/NormData.csv');
%genes = Dexp.textdata(2:end,1);
%headers = Dexp.textdata(1,2:end);
%headers = strrep(headers,'"','');
%for i = 1:length(CellIDs)
%    try
%    indx(i,1) = find(strcmp(headers , CellIDs{i} ));
%    catch
%    end
%end
%
%
%CellIDs = CellIDs(find(indx~=0));
%XYZ = D1.data(find(indx~=0),:);
%
%
%Dexp1 = Dexp.data(:,indx(find(indx~=0)));
%
%Xtrain = XYZ;
%Xtrain(:,1:3) = Xtrain(:,1:3) / 400;
%Xtest = D2(:,[5:7,4]);
%Xtest(:,1:3) =Xtest(:,1:3)/400;
%Ytrain = Dexp1(find(strcmp(genes,gene2test )),:);


%figure(1)
%%scatter3(Xtrain(:,1),Xtrain(:,2),Xtrain(:,3),'bo','fill') 
%hold on
%scatter3(Xtest(:,1),Xtest(:,2),Xtest(:,3),'rs','fill') 
%set(gca,'xtick',[],'ytick',[],'ztick',[])
%caxis manual; caxis([bottom top]);
%view([1 1 1 ])
%legend({'SHOT','Embryo'})
%print(['./Plots/XYZ.pdf'],'-dpdf')
%clf

hyp.cov  = [var(Ytrain); log(0.1); log(0.1); log(0.1)];
hyp.lik  = [var(Ytrain)/3];    
hyp.mean = mean(Ytrain);       
hyp2 = minimize(hyp, @gp, -100, 'infExact', @meanConst, @covSEard, 'likGauss', Xtrain(:,1:3), Ytrain');
[m s2] = gp(hyp2, 'infExact', @meanConst, @covSEard, 'likGauss',  Xtrain(:,1:3), Ytrain', Xtest(:,1:3));
    
%figure(1)
%subplot(4,4,[1,2,5,6])
%scatter3(Xtrain(:,1),Xtrain(:,2),Xtrain(:,3),25,Ytrain,'fill') 
%title(gene2test)
%view([1 1 1 ])
%subplot(4,4,[9,10,13,14])
%scatter3(Xtest(:,1),Xtest(:,2),Xtest(:,3),25,m,'fill') 
%title('Marmoset')
%view([1 1 1 ])

DUMAP = importdata('/Users/christopherpenfold/Desktop/Thorsten/FINAL/Mouse2_Modelling20k/DimRed/Embeddings.csv')
Distancematrix = dist(DUMAP.data'); 



headers2 = Dexp2.textdata(1,2:end);
umapid = DUMAP.textdata(2:end,1);
for i = 1:length(trainingset)
    
    %Get the closest cynomologous data point and replace 
    umindx(i,1) = find(strcmp(umapid , trainingset{i} ));
    
    
    %umindx(i,1) = find(strcmp(umapid , CellIDs{i} ));
    coord1 = Distancematrix(umindx(i,1),:);    
    
    minD = min(coord1(find(contains(umapid , '3289STDY'))));     
    
    
    closestInd = umapid{find(coord1==minD)}
    
    %Something weird happening ... seurat addig X to header row for UMAP
    %output
    
    Clindx(i,1) = find(strcmp(headers2 , closestInd(2:end) ));
    
    Type1Store{i,1} = closestInd(2:end);
    
    %Get the closes human data point from a mix of the datasets
    %minD = min(coord1(find(contains(umapid , 'qc.'))));     
    %minD2 = min(coord1(find(contains(umapid , 'Aligned.sorted')~=1 & contains(umapid , 'qc.')~=1) ))
    
    %minD2 = min(coord1(find(contains(umapid , 'SLX_'))));    
    %closestInd2 = DUMAP.textdata{find(coord1==minD2)}
    %Clindx2(i,1) = find(strcmp(headers , closestInd2 ));
    
    %Type2Store{i,1} = closestInd2;
    
end



%Dexp2 = Dexp.data(:,Clindx);
Dexp3 = Dexp2.data(:,Clindx);


%XYZ = D1.data(find(indx~=0),:);
%Xtrain = XYZ(:,[1:4]) 
%Xtrain(:,1:3) = Xtrain(:,1:3) / 400;
%Xtest = D2(:,[5:7,4]) 
%Xtest(:,1:3) = Xtest(:,1:3) / 400;
genes = Dexp2.textdata(2:end,1);
Ytrain3 = Dexp3(find(strcmp(genes,gene2test )),:);
Ytrain3 = (Ytrain3-mean(Ytrain3))./std((Ytrain3-mean(Ytrain3)) ); 


hyp.cov  = [var(Ytrain3); log(0.1); log(0.1); log(0.1)];
hyp.lik  = [var(Ytrain3)/3];    
hyp.mean = mean(Ytrain3);       
hyp3 = minimize(hyp, @gp, -100, 'infExact', @meanConst, @covSEard, 'likGauss', Xtrain(:,1:3), Ytrain3');
[m2 s2] = gp(hyp3, 'infExact', @meanConst, @covSEard, 'likGauss',  Xtrain(:,1:3), Ytrain3', Xtest(:,1:3));


subplot(4,4,[1,2,5,6])
scatter3(Xtest(:,1),Xtest(:,2),Xtest(:,3),25,m,'fill') 
title('Spatial trascriptomics')
view([-1 -1 1 ])

subplot(4,4,[3,4,7,8])
scatter3(Xtest(:,1),Xtest(:,2),Xtest(:,3),25,m,'fill') 
title('Spatial trascriptomics')
view([1 1 1 ])


subplot(4,4,[9,10,13,14])
scatter3(Xtest(:,1),Xtest(:,2),Xtest(:,3),25,m2,'fill') 
title('Single cell')
view([-1 -1 1 ])

subplot(4,4,[11,12,15,16])
scatter3(Xtest(:,1),Xtest(:,2),Xtest(:,3),25,m2,'fill') 
title('Single cell')
view([1 1 1 ])


return

%scatter3(XYZ(:,1),XYZ(:,2),XYZ(:,3),15,Dexp(find(strcmp(genes,'SOX17' )),:),'filled')

%Load in the normalised expression data from Seurat

%plot3(D2(:,5),D2(:,6),D2(:,7),'o')
%hold on
%plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),'ro')
M1 = zeros(2000,size(Xtrain,1)+size(Xtest,1));
M2 = zeros(2000,size(Xtrain,1)+size(Xtest,1));
M2B = zeros(2000,size(Xtrain,1)+size(Xtest,1));
vM1 = zeros(2000,size(Xtrain,1)+size(Xtest,1));
vM2 = zeros(2000,size(Xtrain,1)+size(Xtest,1));
vM2B = zeros(2000,size(Xtrain,1)+size(Xtest,1));


XY = repmat([Xtest(:,1),Xtest(:,2)],10,1)  ;
Z = repelem(linspace(min(Xtest(:,3)),max(Xtest(:,3)),10),size(Xtest,1))';
idZ = repmat(Xtest(:,4),10,1)  ;


M3 = zeros(2000,size(XY,1));
M4 = zeros(2000,size(XY,1));
M4B = zeros(2000,size(XY,1));
vM3 = zeros(2000,size(XY,1));
vM4 = zeros(2000,size(XY,1));
vM4B = zeros(2000,size(XY,1));

for i = 1:2000

Ytrain = Dexp1(i,:);    
hyp.cov  = [var(Ytrain); log(0.1); log(0.1); log(0.1)];
hyp.lik  = [var(Ytrain)/3];    
hyp.mean = mean(Ytrain);       
hyp2 = minimize(hyp, @gp, -100, 'infExact', @meanConst, @covSEard, 'likGauss', Xtrain(:,1:3), Ytrain');
[m s] = gp(hyp2, 'infExact', @meanConst, @covSEard, 'likGauss',  Xtrain(:,1:3), Ytrain', [Xtrain(:,1:3);Xtest(:,1:3)]);
[m_ s_] = gp(hyp2, 'infExact', @meanConst, @covSEard, 'likGauss',  Xtrain(:,1:3), Ytrain', [XY,Z]);

Ytrain2 = Dexp2(i,:);    
hyp.cov  = [var(Ytrain2); log(0.1); log(0.1); log(0.1)];
hyp.lik  = [var(Ytrain2)/3];    
hyp.mean = mean(Ytrain2);       
hyp3 = minimize(hyp, @gp, -100, 'infExact', @meanConst, @covSEard, 'likGauss', Xtrain(:,1:3), Ytrain2');
[m2 s2] = gp(hyp3, 'infExact', @meanConst, @covSEard, 'likGauss',  Xtrain(:,1:3), Ytrain2', [Xtrain(:,1:3);Xtest(:,1:3)]);
[m2_ s2_] = gp(hyp2, 'infExact', @meanConst, @covSEard, 'likGauss',  Xtrain(:,1:3), Ytrain2', [XY,Z]);

Ytrain3 = Dexp3(i,:);    
hyp.cov  = [var(Ytrain2); log(0.1); log(0.1); log(0.1)];
hyp.lik  = [var(Ytrain2)/3];    
hyp.mean = mean(Ytrain2);       
hyp3 = minimize(hyp, @gp, -100, 'infExact', @meanConst, @covSEard, 'likGauss', Xtrain(:,1:3), Ytrain3');
[m3 s3] = gp(hyp3, 'infExact', @meanConst, @covSEard, 'likGauss',  Xtrain(:,1:3), Ytrain3', [Xtrain(:,1:3);Xtest(:,1:3)]);
[m3_ s3_] = gp(hyp2, 'infExact', @meanConst, @covSEard, 'likGauss',  Xtrain(:,1:3), Ytrain3', [XY,Z]);

M1(i,:) = m;
M2(i,:) = m2;
M2B(i,:) = m3;
M3(i,:) = m_;
M4(i,:) = m2_;
M4B(i,:) = m3_;

vM1(i,:) = s;
vM2(i,:) = s2;
vM2B(i,:) = s3;
vM3(i,:) = s_;
vM4(i,:) = s2_;
vM4B(i,:) = s3_;

disp([num2str(i)])

Delta(1,i) = sum( (Ytrain' - m(1:length(Ytrain))).^2 );
Delta(2,i) = sum( (Ytrain2' - m2(1:length(Ytrain))).^2 );
Delta(3,i) = sum( (Ytrain3' - m3(1:length(Ytrain))).^2 );

Delta(4,i) = sum( (m_ - m2_).^2 );
Delta(5,i) = sum( (m_ - m3_).^2 );
Delta(6,i) = sum( (m2_ - m3_).^2 );

end

figure(1)
subplot(3,1,1); bar( sqrt(Delta(1,:)) / length(Xtrain));
subplot(3,1,2); bar( sqrt(Delta(2,:)) / length(Xtrain));
subplot(3,1,3); bar( sqrt(Delta(3,:)) / length(Xtrain));
print(['./Plots/MSE.pdf'],'-dpdf')

figure(2)
subplot(3,1,1); bar( sqrt(Delta(4,:)) / length(Xtrain));
subplot(3,1,2); bar( sqrt(Delta(5,:)) / length(Xtrain));
subplot(3,1,3); bar( sqrt(Delta(6,:)) / length(Xtrain));
print(['./Plots/MSE_vs_marmoset.pdf'],'-dpdf')



Xtests = [Xtrain;Xtest]
for i = 1:100
    
bottom = min(min([Dexp1(i,:);Dexp2(i,:);Dexp3(i,:)]));
top = max(max([Dexp1(i,:);Dexp2(i,:);Dexp3(i,:)]));
 
figure(1)
subplot(3,3,1)
scatter3(Xtrain(:,1),Xtrain(:,2),Xtrain(:,3),25,Dexp1(i,:),'fill') 
set(gca,'xtick',[],'ytick',[],'ztick',[])
title(['Marmoset ' genes{i}])
caxis manual; caxis([bottom top]);
view([1 1 1 ])
subplot(3,3,2)
scatter3(Xtrain(:,1),Xtrain(:,2),Xtrain(:,3),25,Dexp2(i,:),'fill') 
set(gca,'xtick',[],'ytick',[],'ztick',[])
title(['Cynomologous ' genes{i}])
caxis manual; caxis([bottom top]);
view([1 1 1 ])
subplot(3,3,3)
scatter3(Xtrain(:,1),Xtrain(:,2),Xtrain(:,3),25,Dexp3(i,:),'fill') 
set(gca,'xtick',[],'ytick',[],'ztick',[])
title(['Porcine ' genes{i}])
caxis manual; caxis([bottom top]);
view([1 1 1 ])

subplot(3,3,4)
scatter3(Xtests(:,1),Xtests(:,2),Xtests(:,3),25,M1(i,:),'fill') 
set(gca,'xtick',[],'ytick',[],'ztick',[])
title('Marmoset')
caxis manual; caxis([bottom top]);
view([1 1 1 ])

subplot(3,3,5)
scatter3(Xtests(:,1),Xtests(:,2),Xtests(:,3),25,M2(i,:),'fill') 
set(gca,'xtick',[],'ytick',[],'ztick',[])
title('Cynomologous')
caxis manual; caxis([bottom top]);
view([1 1 1 ])

subplot(3,3,6)
scatter3(Xtests(:,1),Xtests(:,2),Xtests(:,3),25,M2B(i,:),'fill') 
set(gca,'xtick',[],'ytick',[],'ztick',[])
title('Porcine')
caxis manual; caxis([bottom top]);
view([1 1 1 ])

bottom = min(min([vM1(i,:),vM2(i,:),vM2B(i,:)]));
top = max(max([vM1(i,:),vM2(i,:),vM2B(i,:)]));

subplot(3,3,7)
scatter3(Xtests(:,1),Xtests(:,2),Xtests(:,3),25,vM1(i,:),'fill') 
set(gca,'xtick',[],'ytick',[],'ztick',[])
title('Marmoset')
caxis manual; caxis([bottom top]);
view([1 1 1 ])

subplot(3,3,8)
scatter3(Xtests(:,1),Xtests(:,2),Xtests(:,3),25,vM2(i,:),'fill') 
set(gca,'xtick',[],'ytick',[],'ztick',[])
title('Cynomologous')
caxis manual; caxis([bottom top]);
view([1 1 1 ])

subplot(3,3,9)
scatter3(Xtests(:,1),Xtests(:,2),Xtests(:,3),25,vM2B(i,:),'fill') 
set(gca,'xtick',[],'ytick',[],'ztick',[])
title('Porcine')
caxis manual; caxis([bottom top]);
view([1 1 1 ])

print(['./Plots/' genes{i} '.pdf'],'-dpdf')
clf
end


for i = 1:100

bottom = min(min([Dexp1(i,:);Dexp2(i,:);Dexp3(i,:)]));
top = max(max([Dexp1(i,:);Dexp2(i,:);Dexp3(i,:)]));

figure(1)
subplot(1,3,1)
scatter3(XY(:,1),XY(:,2),Z(:,1),25,M3(i,:),'fill') 
title(['Marmoset ' genes{i}])
caxis manual; caxis([bottom top]);
view([1 1 1 ])
subplot(1,3,2)
scatter3(XY(:,1),XY(:,2),Z(:,1),25,M4(i,:),'fill') 
title(['Cynomologous ' genes{i}])
caxis manual; caxis([bottom top]);
view([1 1 1 ])
subplot(1,3,3)
scatter3(XY(:,1),XY(:,2),Z(:,1),25,M4B(i,:),'fill') 
title(['Porcine ' genes{i}])
caxis manual; caxis([bottom top]);
view([1 1 1 ])

print(['./Plots/' genes{i} '_section.pdf'],'-dpdf')
clf
end

steps = 1:3:28;
uZ = unique(Z);

figure(1)
scatter3(XY(:,1),XY(:,2),Z(:,1),'bo','fill') 
set(gca,'xtick',[],'ytick',[],'ztick',[])
caxis manual; caxis([bottom top]);
view([1 1 1 ])
%fig = gcf;
%fig.PaperUnits = 'inches';
%fig.PaperPosition = [0 0 3 10];
print(['./Plots/MAXProj.pdf'],'-dpdf')
clf


for i = 1:100

bottom = min(min([Dexp1(i,:);Dexp2(i,:);Dexp3(i,:)]));
top = max(max([Dexp1(i,:);Dexp2(i,:);Dexp3(i,:)]));

for j = 1:10
inds = find(Z==uZ(j));
subplot(10,3,steps(j)) %[1,2,5,6,9,10,13,14])
scatter(XY(inds,1),XY(inds,2),25,M3(i,inds),'fill') 
set(gca,'xtick',[],'ytick',[]) %title(['Marmoset ' genes{i}])
caxis manual; caxis([bottom top]);
subplot(10,3,steps(j)+1) %[3,4,7,8,11,12,15,16])
scatter(XY(inds,1),XY(inds,2),25,M4(i,inds),'fill') 
%title(['Cynomologous ' genes{i}])
set(gca,'xtick',[],'ytick',[]) 
caxis manual; caxis([bottom top]);
subplot(10,3,steps(j)+2)%[3,4,7,8,11,12,15,16])
scatter(XY(inds,1),XY(inds,2),25,M4B(i,inds),'fill') 
set(gca,'xtick',[],'ytick',[]) 
%title(['Porcine ' genes{i}])
caxis manual; caxis([bottom top]);
end
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 3 10];
print(['./Plots/' genes{i} '_section2.pdf'],'-dpdf')
clf
end


%COEFF = pca(M1(:,85:end));
%
%plot(COEFF(D2(:,4)==0,1),COEFF(D2(:,4)==0,2),'ro','MarkerFaceColor','r')
%hold on
%plot(COEFF(D2(:,4)==1,1),COEFF(D2(:,4)==1,2),'bo','MarkerFaceColor','b')
%plot(COEFF(D2(:,4)==2,1),COEFF(D2(:,4)==2,2),'go','MarkerFaceColor','g')
%plot(COEFF(D2(:,4)==3,1),COEFF(D2(:,4)==3,2),'ko','MarkerFaceColor','k')
%plot(COEFF(D2(:,4)==4,1),COEFF(D2(:,4)==4,2),'yo','MarkerFaceColor','y')


for i = 1:2000
    
bottom = min(min([Dexp1(i,:);Dexp2(i,:);Dexp3(i,:)]));
top = max(max([Dexp1(i,:);Dexp2(i,:);Dexp3(i,:)]));
 
figure(1)
subplot(6,3,1)
scatter3(Xtrain(:,1),Xtrain(:,2),Xtrain(:,3),25,Dexp1(i,:),'fill') 
set(gca,'xtick',[],'ytick',[],'ztick',[])
title(['Marmoset ' genes{i}])
caxis manual; caxis([bottom top]);
view([1 1 1 ])
subplot(6,3,2)
scatter3(Xtrain(:,1),Xtrain(:,2),Xtrain(:,3),25,Dexp2(i,:),'fill') 
set(gca,'xtick',[],'ytick',[],'ztick',[])
title(['Cynomologous ' genes{i}])
caxis manual; caxis([bottom top]);
view([1 1 1 ])
subplot(6,3,3)
scatter3(Xtrain(:,1),Xtrain(:,2),Xtrain(:,3),25,Dexp3(i,:),'fill') 
set(gca,'xtick',[],'ytick',[],'ztick',[])
title(['Human ' genes{i}])
caxis manual; caxis([bottom top]);
view([1 1 1 ])

indsub = find(Xtests(:,4)==0);
subplot(6,3,4)
scatter3(Xtests(indsub,1),Xtests(indsub,2),Xtests(indsub,3),25,M1(i,indsub),'fill') 
set(gca,'xtick',[],'ytick',[],'ztick',[])
caxis manual; caxis([bottom top]);
view([1 1 1 ])
subplot(6,3,5)
scatter3(Xtests(indsub,1),Xtests(indsub,2),Xtests(indsub,3),25,M2(i,indsub),'fill') 
set(gca,'xtick',[],'ytick',[],'ztick',[])
caxis manual; caxis([bottom top]);
view([1 1 1 ])
subplot(6,3,6)
scatter3(Xtests(indsub,1),Xtests(indsub,2),Xtests(indsub,3),25,M2B(i,indsub),'fill') 
set(gca,'xtick',[],'ytick',[],'ztick',[])
caxis manual; caxis([bottom top]);
view([1 1 1 ])

indsub = find(Xtests(:,4)==1);
subplot(6,3,7)
scatter3(Xtests(indsub,1),Xtests(indsub,2),Xtests(indsub,3),25,M1(i,indsub),'fill') 
set(gca,'xtick',[],'ytick',[],'ztick',[])
caxis manual; caxis([bottom top]);
view([1 1 1 ])
subplot(6,3,8)
scatter3(Xtests(indsub,1),Xtests(indsub,2),Xtests(indsub,3),25,M2(i,indsub),'fill') 
set(gca,'xtick',[],'ytick',[],'ztick',[])
caxis manual; caxis([bottom top]);
view([1 1 1 ])
subplot(6,3,9)
scatter3(Xtests(indsub,1),Xtests(indsub,2),Xtests(indsub,3),25,M2B(i,indsub),'fill') 
set(gca,'xtick',[],'ytick',[],'ztick',[])
caxis manual; caxis([bottom top]);
view([1 1 1 ])

indsub = find(Xtests(:,4)==2);
subplot(6,3,10)
scatter3(Xtests(indsub,1),Xtests(indsub,2),Xtests(indsub,3),25,M1(i,indsub),'fill') 
set(gca,'xtick',[],'ytick',[],'ztick',[])
caxis manual; caxis([bottom top]);
view([1 1 1 ])
subplot(6,3,11)
scatter3(Xtests(indsub,1),Xtests(indsub,2),Xtests(indsub,3),25,M2(i,indsub),'fill') 
set(gca,'xtick',[],'ytick',[],'ztick',[])
caxis manual; caxis([bottom top]);
view([1 1 1 ])
subplot(6,3,12)
scatter3(Xtests(indsub,1),Xtests(indsub,2),Xtests(indsub,3),25,M2B(i,indsub),'fill') 
set(gca,'xtick',[],'ytick',[],'ztick',[])
caxis manual; caxis([bottom top]);
view([1 1 1 ])

indsub = find(Xtests(:,4)==3);
subplot(6,3,13)
scatter3(Xtests(indsub,1),Xtests(indsub,2),Xtests(indsub,3),25,M1(i,indsub),'fill') 
set(gca,'xtick',[],'ytick',[],'ztick',[])
caxis manual; caxis([bottom top]);
view([1 1 1 ])
subplot(6,3,14)
scatter3(Xtests(indsub,1),Xtests(indsub,2),Xtests(indsub,3),25,M2(i,indsub),'fill') 
set(gca,'xtick',[],'ytick',[],'ztick',[])
caxis manual; caxis([bottom top]);
view([1 1 1 ])
subplot(6,3,15)
scatter3(Xtests(indsub,1),Xtests(indsub,2),Xtests(indsub,3),25,M2B(i,indsub),'fill') 
set(gca,'xtick',[],'ytick',[],'ztick',[])
caxis manual; caxis([bottom top]);
view([1 1 1 ])

indsub = find(Xtests(:,4)==4);
subplot(6,3,16)
scatter3(Xtests(indsub,1),Xtests(indsub,2),Xtests(indsub,3),25,M1(i,indsub),'fill') 
set(gca,'xtick',[],'ytick',[],'ztick',[])
caxis manual; caxis([bottom top]);
view([1 1 1 ])
subplot(6,3,17)
scatter3(Xtests(indsub,1),Xtests(indsub,2),Xtests(indsub,3),25,M2(i,indsub),'fill') 
set(gca,'xtick',[],'ytick',[],'ztick',[])
caxis manual; caxis([bottom top]);
view([1 1 1 ])
subplot(6,3,18)
scatter3(Xtests(indsub,1),Xtests(indsub,2),Xtests(indsub,3),25,M2B(i,indsub),'fill') 
set(gca,'xtick',[],'ytick',[],'ztick',[])
caxis manual; caxis([bottom top]);
view([1 1 1 ])

%bottom = min(min([vM1(i,:),vM2(i,:),vM2B(i,:)]));
%#top = max(max([vM1(i,:),vM2(i,:),vM2B(i,:)]));
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6 10];
print(['./Plots/' genes{i} '_seperate.pdf'],'-dpdf')
clf
end



set(h4,'PaperSize',[12 4*12 ]);
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 6*2 12*2];
%set(h4,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])


print('~/Desktop/TestSC','-dpdf','-r0')
print('~/Desktop/TestSC','-svg','-r0')



%hold on
%plot(COEFF(1:84,1),COEFF(1:84,2),'ro')
%COEFF = pca([Dexp2,M3])
