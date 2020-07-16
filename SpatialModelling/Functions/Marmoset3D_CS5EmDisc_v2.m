function [m_1,s_1,Output] = Marmoset3D_CS5EmDisc_v2(Dexp1,Output,Model,Arc,Arc2,gene2test)

LabBin = Output.LabBin;
genes = Output.genes;
ind0 = Output.ind0;
ind1 = Output.ind1;
ind2 = Output.ind2;
ind3 = Output.ind3;
ind4 = Output.ind4;
ind5 = Output.ind5;

Xtrain = Output.Xtrain;
Xtest = Output.Xtest;


pg1 = {@priorGauss,0,2};  
pg2 = {@priorGauss,0,2};  
pg3 = {@priorGauss,log(0.5),2};  
pc = {@priorClamped};

%Prior for base model
prior1.mean = {[]};  
prior1.cov  = {pg1;pg2;pg2;pg2};
prior1.lik  = {pg3};

%Prior for long lenghth-scale model
prior2.mean = {[]};  
prior2.cov  = {pg1;pc;pc;pc};
prior2.lik  = {pg3};

%Prior for exteded model
prior3.mean = {[]}; 
prior3.cov  = {pg1;pg2;pg2;pg2;[];[];[];[]}; 
prior3.lik  = {pg3};


im1 = {@infPrior,@infExact,prior1};                % inference method
im2 = {@infPrior,@infExact,prior2};                % inference method
im3 = {@infPrior,@infExact,prior3};                % inference method

gind = find(strcmp(genes,gene2test)==1);

Ytrain = Dexp1(gind,:);


if strcmp(Model,'Base')==1
    
    
Ys = nthroot((reshape(Output.P(find(LabBin==1),:),size(Output.P(find(LabBin==1),:),1)*size(Output.P(find(LabBin==1),:),2),1) ), 1)+0.00001   ;
    
pg1 = {@priorGauss,0,1};  
pg2 = {@priorGauss,0,1};  
pg3 = {@priorGauss,log(0.5),1};  
pc = {@priorClamped};

%Prior for base model
prior1.mean = {[]};  
prior1.cov  = {[];[];[]}
%prior1.cov  = {pg1;pg2;pg3};
prior1.lik  = {[]};

%Prior for long lenghth-scale model
prior2.mean = {[]};  
prior2.cov  = {pg1;pc};
prior2.lik  = {[]};  %log(4.7);

%{@likBeta,'logit'};     phi = 4.7; hyp.lik = log(phi);


im1 = {@infPrior,@infLaplace,prior1};                % inference method
im2 = {@infPrior,@infLaplace,prior2};                % inference method


x=repmat(Output.XYZP(find(LabBin==1),1),size(Output.P,2),1);
y=repmat(Output.XYZP(find(LabBin==1),2),size(Output.P,2),1);
z=repmat(Output.XYZP(find(LabBin==1),3),size(Output.P,2),1);    
    
%Initialise first guess of hyperparameters
hyp1.cov  = [var(Ys); log(1); log(0.1)];
hyp1.lik  = log(1);  
hyp1.mean = 0;     

hyp1B.cov  = [var(Ys); log(40)];
hyp1B.lik  = log(4.7);   
hyp1B.mean = mean(Ys);       

lik = {'likBeta','expexp'}

%par  = {@meanConst,'covSEiso',lik,[x,y,z], reshape(Output.P,size(Output.P,1)*size(Output.P,2),1)};
%keyboard
%Optimise and get MLs/predictions

hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, @meanConst,{@covSum,{@covSEiso,@covNoise}},lik,[x,y,z],Ys);         % optimise

%hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
%hyp3 = feval(@minimize, hyp1B, @gp, -1000, im2, par{:});         % optimise
[m_1 s_1] = gp(hyp2, @infLaplace, @meanConst, {@covSum,{@covSEiso,@covNoise}}, {@likBeta,'expexp'},  [x,y,z], Ys, [Xtest(ind1,1:3)]);

%[m_1 s_1] = gp(hyp2, @infLaplace, @meanConst, {@covSum,{@covSEiso,@covNoise}}, {@likBeta,'expexp'},  [x,y,z], Ys, [x,y,z]);

%keyboard
%

%[L1 ] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(:,1:3), Ytrain');
%[L2 ] = gp(hyp3, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(:,1:3), Ytrain');

Output.Ytrain = Ytrain;
Output.Xtrain = Xtrain;
Output.m_1 = m_1;
%Output.L1 = L1;
%Output.L2 = L2;
Output.hyp2 = hyp1;
Output.hyp2 = hyp2;
%Output.hyp3 = hyp3;

elseif strcmp(Model,'Base2')==1


pg1 = {@priorGauss,0,1};  
pg2 = {@priorGauss,0,1};  
pg3 = {@priorGauss,log(0.5),1};  
pc = {@priorClamped};

%Prior for base model
prior1.mean = {[]};  
prior1.cov  = {pg1;pg2};
prior1.lik  = {pg3};

%Prior for long lenghth-scale model
prior2.mean = {[]};  
prior2.cov  = {pc;pg1};
prior2.lik  = {pg3};


im1 = {@infPrior,@infExact,prior1};                % inference method
im2 = {@infPrior,@infExact,prior2};                % inference method


        
    HYP2 = cell(1,6);
    HYP3 = cell(1,6);
    
hyp1.cov  = [log(2) ; var(Ytrain)];
hyp1.lik  = [var(Ytrain)/2];    
hyp1.mean = mean(Ytrain);       

hyp1B.cov  = [40; var(Ytrain)];
hyp1B.lik  = [var(Ytrain)/2];    
hyp1B.mean = mean(Ytrain); 

par = {@meanConst,@covSEiso,'likGauss',Xtrain(:,1:3), Ytrain'};
%par1 = {@meanConst,@covSEard,'likGauss',Xtrain(:,1:3), Ytrain',hyp1};
%par2 = {@meanConst,@covSEard,'likGauss',Xtrain(:,1:3), Ytrain',hyp3};
    
indT = find(LabBin==0);
par = {@meanConst,@covSEiso,'likGauss',Xtrain(indT,1:3), Ytrain(indT)'};
hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
hyp3 = feval(@minimize, hyp1B, @gp, -1000, im2, par{:});         % optimise

[m0 s0] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', [Xtest([ind0],1:3)]);


[L1_0 ] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');
[L2_0 ] = gp(hyp3, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');

    
    HYP2{1,1} = hyp2;
    HYP3{1,1} = hyp3;
    
    
indT = find(LabBin==1);
par = {@meanConst,@covSEiso,'likGauss',Xtrain(indT,1:3), Ytrain(indT)'};
hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
hyp3 = feval(@minimize, hyp1B, @gp, -1000, im2, par{:});         % optimise

[m1 s1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', [Xtest([ind1],1:3)]);

[L1_1 ] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');
[L2_1 ] = gp(hyp3, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');

Vind = var(Ytrain(find(LabBin ==1)));

[marc1 sarc1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', Arc);
[marc2 sarc2] = gp(hyp3, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', Arc);


    
    HYP2{1,2} = hyp2;
    HYP3{1,2} = hyp3;    

indT = find(LabBin==2);
par = {@meanConst,@covSEiso,'likGauss',Xtrain(indT,1:3), Ytrain(indT)'};
hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
hyp3 = feval(@minimize, hyp1B, @gp, -1000, im2, par{:});         % optimise
[m2 s2] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', [Xtest([ind2],1:3)]);

    HYP2{1,3} = hyp2;
    HYP3{1,3} = hyp3;


[L1_2 ] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');
[L2_2 ] = gp(hyp3, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');




indT = find(LabBin==3);
par = {@meanConst,@covSEiso,'likGauss',Xtrain(indT,1:3), Ytrain(indT)'};
hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
hyp3 = feval(@minimize, hyp1B, @gp, -1000, im2, par{:});         % optimise

[m3 s3] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', [Xtest([ind3],1:3)]);

    HYP2{1,4} = hyp2;
    HYP3{1,4} = hyp3;


[L1_3 ] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');
[L2_3 ] = gp(hyp3, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');


indT = find(LabBin==4);
par = {@meanConst,@covSEiso,'likGauss',Xtrain(indT,1:3), Ytrain(indT)'};
hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
hyp3 = feval(@minimize, hyp1B, @gp, -1000, im2, par{:});         % optimise

[m4 s4] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', [Xtest([ind4],1:3)]);

    HYP2{1,5} = hyp2;
    HYP3{1,5} = hyp3;


    Vind2 = var(Ytrain(find(LabBin ==4)));

[marc3 sarc3] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', Arc2);
[marc4 sarc4] = gp(hyp3, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', Arc2);

    
[L1_4 ] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');
[L2_4 ] = gp(hyp3, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');



indT = find(LabBin==5);
par = {@meanConst,@covSEiso,'likGauss',Xtrain(indT,1:3), Ytrain(indT)'};
hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
hyp3 = feval(@minimize, hyp1B, @gp, -1000, im2, par{:});         % optimise

    HYP2{1,6} = hyp2;
    HYP3{1,6} = hyp3;


[m5 s5] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', [Xtest([ind5],1:3)]);
    
[L1_5 ] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');
[L2_5 ] = gp(hyp3, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');

    
m_1 = [m0;m1;m2;m3;m4;m5];
s_1 = [s0;s1;s2;s3;s4;s5];

Output.m_1 = m_1;
Output.s_1 = s_1;
Output.Ytrain = Ytrain;
Output.Xtrain = Xtrain;
Output.L1 = [L1_0,L1_1,L1_2,L1_3,L1_4,L1_5];
Output.L2 = [L2_0,L2_1,L2_2,L2_3,L2_4,L2_5];
%Output.hyp2 = hyp1;
Output.hyp2 = HYP2;
Output.hyp3 = HYP3;
    Output.marc1 = flipdim(marc1,1);
Output.marc2 = flipdim(marc2,1);
   Output.marc3 = flipdim(marc3,1);
Output.marc4 = flipdim(marc4,1);

Output.sarc1 = flipdim(sarc1,1);
Output.sarc2 = flipdim(sarc2,1);
Output.sarc3 = flipdim(sarc3,1);
Output.sarc4 = flipdim(sarc4,1);

Output.Vind = Vind;
Output.Vind2 = Vind2;

elseif strcmp(Model,'Independent')==1

        
    HYP2 = cell(1,4);
    HYP3 = cell(1,4);
    
hyp1.cov  = [var(Ytrain); log(2); log(2); log(2)];
hyp1.lik  = [var(Ytrain)/2];    
hyp1.mean = mean(Ytrain);     

hyp1B.cov  = [var(Ytrain); log(40); log(40); log(40)];
hyp1B.lik  = [var(Ytrain)/2];    
hyp1B.mean = mean(Ytrain);  

par = {@meanConst,@covSEard,'likGauss',Xtrain(:,1:3), Ytrain'};
%par1 = {@meanConst,@covSEard,'likGauss',Xtrain(:,1:3), Ytrain',hyp1};
%par2 = {@meanConst,@covSEard,'likGauss',Xtrain(:,1:3), Ytrain',hyp3};
    
indT = find(LabBin==0 | LabBin==1);
par = {@meanConst,@covSEard,'likGauss',Xtrain(indT,1:3), Ytrain(indT)'};
hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
hyp3 = feval(@minimize, hyp1B, @gp, -1000, im2, par{:});         % optimise

[m0 s0] = gp(hyp2, 'infExact', @meanConst, @covSEard, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', [Xtest([ind0],1:3)]);
[m1 s1] = gp(hyp2, 'infExact', @meanConst, @covSEard, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', [Xtest([ind1],1:3)]);


[L1_0 ] = gp(hyp2, 'infExact', @meanConst, @covSEard, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');
[L2_0 ] = gp(hyp3, 'infExact', @meanConst, @covSEard, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');

    
    HYP2{1,1} = hyp2;
    HYP3{1,1} = hyp3;

indT = find(LabBin==2);
par = {@meanConst,@covSEard,'likGauss',Xtrain(indT,1:3), Ytrain(indT)'};
hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
hyp3 = feval(@minimize, hyp1B, @gp, -1000, im2, par{:});         % optimise
[m2 s2] = gp(hyp2, 'infExact', @meanConst, @covSEard, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', [Xtest([ind2],1:3)]);

    HYP2{1,2} = hyp2;
    HYP3{1,2} = hyp3;


[L1_2 ] = gp(hyp2, 'infExact', @meanConst, @covSEard, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');
[L2_2 ] = gp(hyp3, 'infExact', @meanConst, @covSEard, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');




indT = find(LabBin==3 | LabBin==4);
par = {@meanConst,@covSEard,'likGauss',Xtrain(indT,1:3), Ytrain(indT)'};
hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
hyp3 = feval(@minimize, hyp1B, @gp, -1000, im2, par{:});         % optimise

[m3 s3] = gp(hyp2, 'infExact', @meanConst, @covSEard, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', [Xtest([ind3],1:3)]);
[m4 s4] = gp(hyp2, 'infExact', @meanConst, @covSEard, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', [Xtest([ind4],1:3)]);

    HYP2{1,3} = hyp2;
    HYP3{1,3} = hyp3;


[L1_3 ] = gp(hyp2, 'infExact', @meanConst, @covSEard, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');
[L2_3 ] = gp(hyp3, 'infExact', @meanConst, @covSEard, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');


indT = find(LabBin==5);
par = {@meanConst,@covSEard,'likGauss',Xtrain(indT,1:3), Ytrain(indT)'};
hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
hyp3 = feval(@minimize, hyp1B, @gp, -1000, im2, par{:});         % optimise

    HYP2{1,4} = hyp2;
    HYP3{1,4} = hyp3;


[m5 s5] = gp(hyp2, 'infExact', @meanConst, @covSEard, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', [Xtest([ind5],1:3)]);
    
[L1_5 ] = gp(hyp2, 'infExact', @meanConst, @covSEard, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');
[L2_5 ] = gp(hyp3, 'infExact', @meanConst, @covSEard, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');

    
m_1 = [m0;m1;m2;m3;m4;m5];
s_1 = [s0;s1;s2;s3;s4;s5];

Output.m_1 = m_1;
Output.s_1 = s_1;
Output.Ytrain = Ytrain;
Output.Xtrain = Xtrain;
Output.L1 = [L1_0,L1_2,L1_3,L1_5];
Output.L2 = [L2_0,L2_2,L2_3,L2_5];
Output.hyp2 = hyp1;
Output.hyp2 = HYP2;
Output.hyp3 = HYP3;
    
    
elseif strcmp(Model,'Independent2')==1
 %Prior for exteded model
prior3.mean = {[]}; 
prior3.cov  = {pg1;pg2;pg2;pg2;[];[];[];[];[];[]}; 
prior3.lik  = {pg3}


im1 = {@infPrior,@infExact,prior3};                % inference method
    
    
xtaux = zeros(size(LabBin,1),6);
xtaux(find(LabBin==0 ),1) = ones;
xtaux(find(LabBin==1 ),2) = ones;
xtaux(find(LabBin==2),3) = ones;
xtaux(find(LabBin==3),4) = ones;
xtaux(find(LabBin==4),5) = ones;
xtaux(find(LabBin==5),6) = ones;

xpaux = zeros(size(Xtest,1),6);
xpaux([ind0],1) = ones;
xpaux([ind1],2) = ones;
xpaux([ind2],3) = ones;
xpaux([ind3],4) = ones;
xpaux([ind4],5) = ones;
xpaux([ind5],6) = ones;

hyp1.cov  = [var(Ytrain); log(2); log(2); log(2); log(2); log(2); log(2); log(2); log(2); log(2)];
hyp1.lik  = [var(Ytrain)/2];    
hyp1.mean = mean(Ytrain);       

par = {@meanConst,@covSEard,'likGauss',[Xtrain(:,1:3),xtaux], Ytrain'};

hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
[m_1 s_1] = gp(hyp2, 'infExact', @meanConst, @covSEard, 'likGauss',  [Xtrain(:,1:3),xtaux], Ytrain', [Xtest(:,1:3),xpaux]);
            
else
    
    disp('Model not defined ...')
end
    
