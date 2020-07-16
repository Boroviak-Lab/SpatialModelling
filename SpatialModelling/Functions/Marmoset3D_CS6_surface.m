function [Output] = Marmoset3D_CS6_surface(Dexp1,Output,Model,gene2test)

LabBin = Output.LabBin;
genes = Output.genes;
ind0 = Output.ind0;
ind1 = Output.ind1;
ind2 = Output.ind2;
ind3 = Output.ind3;
ind4 = Output.ind4;
ind5 = Output.ind5;
ind6 = Output.ind6;

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

    
    
%Initialise first guess of hyperparameters
hyp1.cov  = [var(Ytrain); log(2)];
hyp1.lik  = [var(Ytrain)/2];    
hyp1.mean = mean(Ytrain);       

par  = {@meanConst,@covSEiso,'likGauss',Xtrain(:,1:3), Ytrain'};

%Optimise and get MLs/predictions
hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
hyp3 = feval(@minimize, hyp1, @gp, -1000, im2, par{:});         % optimise

[m_1 s_1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(:,1:3), Ytrain', [Xtest(:,1:3)]);
    
[L1 ] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(:,1:3), Ytrain');
[L2 ] = gp(hyp3, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(:,1:3), Ytrain');

Output.Ytrain = Ytrain;
Output.Xtrain = Xtrain;
Output.m_1 = m_1;
Output.L1 = L1;
Output.L2 = L2;
Output.hyp2 = hyp1;
Output.hyp2 = hyp2;
Output.hyp3 = hyp3;

elseif strcmp(Model,'Base2')==1
    
OBJ1=importdata('/Volumes/GoogleDrive/My Drive/Marmoset_shared_folder/3D_Modellinng/Blender meshes/Cross sections/CS6_section_cs2.mat');

%Xex2 = OBJ1.vertices(OBJ1.objects(12).data.vertices,:);
%Xem2 = OBJ1.vertices(OBJ1.objects(8).data.vertices,:);
%Xtroph2 = OBJ1.vertices(OBJ1.objects(20).data.vertices,:);
%Xam2 = OBJ1.vertices(OBJ1.objects(4).data.vertices,:);
%Xsys2 = OBJ1.vertices(OBJ1.objects(16).data.vertices,:);
%Xve2 = OBJ1.vertices(OBJ1.objects(24).data.vertices,:);
%D2 = [Xam2,zeros(size(Xam2,1),1) ; Xem2,ones(size(Xem2,1),1) ; Xsys2,2*ones(size(Xsys2,1),1) ; Xex2,3*ones(size(Xex2,1),1) ; Xve2,4*ones(size(Xve2,1),1) ; Xtroph2,5*ones(size(Xtroph2,1),1) ];
    
    
   
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

    HYP2 = cell(1,7);
    HYP3 = cell(1,7);    

hyp1.cov  = [log(2); var(Ytrain)];
hyp1.lik  = [var(Ytrain)/2];    
hyp1.mean = mean(Ytrain);       

hyp1B.cov  = [(40) ; var(Ytrain)];
hyp1B.lik  = [var(Ytrain)/2];    
hyp1B.mean = mean(Ytrain);       

par = {@meanConst,@covSEiso,'likGauss',Xtrain(:,1:3), Ytrain'};

indT = find(LabBin==0);
par = {@meanConst,@covSEiso,'likGauss',Xtrain(indT,1:3), Ytrain(indT)'};
hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
hyp3 = feval(@minimize, hyp1B, @gp, -1000, im2, par{:});         % optimise
    HYP2{1,1} = hyp2;
    HYP3{1,1} = hyp3;
[m0 s0] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', OBJ1.vertices/400);
[L1_0 ] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');
[L2_0 ] = gp(hyp3, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');

    
indT = find(LabBin==1);
par = {@meanConst,@covSEiso,'likGauss',Xtrain(indT,1:3), Ytrain(indT)'};
hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
hyp3 = feval(@minimize, hyp1B, @gp, -1000, im2, par{:});         % optimise
    HYP2{1,2} = hyp2;
    HYP3{1,2} = hyp3;
[m1 s1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', OBJ1.vertices/400);
[L1_1 ] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');
[L2_1 ] = gp(hyp3, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');



indT = find(LabBin==2);
par = {@meanConst,@covSEiso,'likGauss',Xtrain(indT,1:3), Ytrain(indT)'};
hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
hyp3 = feval(@minimize, hyp1B, @gp, -1000, im2, par{:});         % optimise
[m2 s2] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', OBJ1.vertices/400);
    HYP2{1,3} = hyp2;
    HYP3{1,3} = hyp3;    
[L1_2 ] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');
[L2_2 ] = gp(hyp3, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');


%keyboard
indT = find(LabBin==3);
par = {@meanConst,@covSEiso,'likGauss',Xtrain(indT,1:3), Ytrain(indT)'};
hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
hyp3 = feval(@minimize, hyp1B, @gp, -1000, im2, par{:});         % optimise
    HYP2{1,4} = hyp2;
    HYP3{1,4} = hyp3;
Vind = var(Ytrain(find(LabBin==1 )));
[m3 s3] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', OBJ1.vertices/400);
[L1_3 ] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');
[L2_3 ] = gp(hyp3, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');


indT = find(LabBin==4);
par = {@meanConst,@covSEiso,'likGauss',Xtrain(indT,1:3), Ytrain(indT)'};
hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
hyp3 = feval(@minimize, hyp1B, @gp, -1000, im2, par{:});         % optimise
[m4 s4] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', OBJ1.vertices/400);
    HYP2{1,5} = hyp2;
    HYP3{1,5} = hyp3;    
[L1_4 ] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');
[L2_4 ] = gp(hyp3, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');


indT = find(LabBin==5);
par = {@meanConst,@covSEiso,'likGauss',Xtrain(indT,1:3), Ytrain(indT)'};
hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
hyp3 = feval(@minimize, hyp1B, @gp, -1000, im2, par{:});         % optimise
[m5 s5] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', OBJ1.vertices/400);
    HYP2{1,6} = hyp2;
    HYP3{1,6} = hyp3;
[L1_5 ] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');
[L2_5 ] = gp(hyp3, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');
    


indT = find(LabBin==6);
par = {@meanConst,@covSEiso,'likGauss',Xtrain(indT,1:3), Ytrain(indT)'};
hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
hyp3 = feval(@minimize, hyp1B, @gp, -1000, im2, par{:});         % optimise
[m6 s6] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', OBJ1.vertices/400);
[L1_6 ] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');
[L2_6 ] = gp(hyp3, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');
    HYP2{1,7} = hyp2;
    HYP3{1,7} = hyp3;
 
Output.m_0 = m0;
Output.m_1 = m1;
Output.m_2 = m2;
Output.m_3 = m3;
Output.m_4 = m4;
Output.m_5 = m5;
Output.m_6 = m6;

%Output.m_1 = m_1;
%Output.s_1 = s_1;
Output.Ytrain = Ytrain;
Output.Xtrain = Xtrain;
Output.L1 = [L1_0,L1_1,L1_2,L1_3,L1_4,L1_5,L1_6];
Output.L2 = [L2_0,L2_1,L2_2,L2_3,L2_4,L2_5,L2_6];
Output.hyp2 = hyp1;
Output.hyp2 = HYP2;
Output.hyp3 = HYP3;
%Output.marc1 = marc1;
%Output.marc2 = marc2;
%Output.sarc1 = sarc1;
%Output.sarc2 = sarc2;
Output.Vind = Vind;

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
    %keyboard
indT = find(LabBin==0);
par = {@meanConst,@covSEard,'likGauss',Xtrain(indT,1:3), Ytrain(indT)'};
hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
hyp3 = feval(@minimize, hyp1B, @gp, -1000, im2, par{:});         % optimise
    HYP2{1,1} = hyp2;
    HYP3{1,1} = hyp3;

[m0 s0] = gp(hyp2, 'infExact', @meanConst, @covSEard, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', [Xtest([ind0],1:3)]);

[L1_0 ] = gp(hyp2, 'infExact', @meanConst, @covSEard, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');
[L2_0 ] = gp(hyp3, 'infExact', @meanConst, @covSEard, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');



indT = find(LabBin==2 | LabBin==4);
par = {@meanConst,@covSEard,'likGauss',Xtrain(indT,1:3), Ytrain(indT)'};
hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
hyp3 = feval(@minimize, hyp1B, @gp, -1000, im2, par{:});         % optimise

[m2 s2] = gp(hyp2, 'infExact', @meanConst, @covSEard, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', [Xtest([ind2],1:3)]);
[m4 s4] = gp(hyp2, 'infExact', @meanConst, @covSEard, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', [Xtest([ind4],1:3)]);

    HYP2{1,2} = hyp2;
    HYP3{1,2} = hyp3;


[L1_2 ] = gp(hyp2, 'infExact', @meanConst, @covSEard, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');
[L2_2 ] = gp(hyp3, 'infExact', @meanConst, @covSEard, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');


indT = find(LabBin==1 | LabBin==3 | LabBin ==6);
par = {@meanConst,@covSEard,'likGauss',Xtrain(indT,1:3), Ytrain(indT)'};
hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
hyp3 = feval(@minimize, hyp1B, @gp, -1000, im2, par{:});         % optimise

    HYP2{1,3} = hyp2;
    HYP3{1,3} = hyp3;


[m1 s1] = gp(hyp2, 'infExact', @meanConst, @covSEard, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', [Xtest([ind1],1:3)]);
[m3 s3] = gp(hyp2, 'infExact', @meanConst, @covSEard, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', [Xtest([ind3],1:3)]);
[m6 s6] = gp(hyp2, 'infExact', @meanConst, @covSEard, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', [Xtest([ind6],1:3)]);

[L1_3 ] = gp(hyp2, 'infExact', @meanConst, @covSEard, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');
[L2_3 ] = gp(hyp3, 'infExact', @meanConst, @covSEard, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');

keyboard
indT = find(LabBin==5);
par = {@meanConst,@covSEard,'likGauss',Xtrain(indT,1:3), Ytrain(indT)'};
hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
hyp3 = feval(@minimize, hyp1B, @gp, -1000, im2, par{:});         % optimise

[m5 s5] = gp(hyp2, 'infExact', @meanConst, @covSEard, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', [Xtest([ind5],1:3)]);

    HYP2{1,4} = hyp2;
    HYP3{1,4} = hyp3;


[L1_5 ] = gp(hyp2, 'infExact', @meanConst, @covSEard, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');
[L2_5 ] = gp(hyp3, 'infExact', @meanConst, @covSEard, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');
    
    
m_1 = [m0;m1;m2;m3;m4;m5;m6];
s_1 = [s0;s1;s2;s3;s4;s5;s6];

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
prior3.cov  = {pg1;pg2;pg2;pg2;[];[];[];[];[];[];[]}; 
prior3.lik  = {pg3}


im1 = {@infPrior,@infExact,prior3};                % inference method
        
    
xtaux = zeros(size(LabBin,1),7);
xtaux(find(LabBin==0 ),1) = ones;
xtaux(find(LabBin==1 ),2) = ones;
xtaux(find(LabBin==2),3) = ones;
xtaux(find(LabBin==3),4) = ones;
xtaux(find(LabBin==4),5) = ones;
xtaux(find(LabBin==5),6) = ones;
xtaux(find(LabBin==6),7) = ones;


xpaux = zeros(size(Xtest,1),7);
xpaux([ind0],1) = ones;
xpaux([ind1],2) = ones;
xpaux([ind2],3) = ones;
xpaux([ind3],4) = ones;
xpaux([ind4],5) = ones;
xpaux([ind5],6) = ones;
xpaux([ind5],7) = ones;

hyp1.cov  = [var(Ytrain); log(2); log(2); log(2); log(2); log(2); log(2); log(2); log(2); log(2); log(2)];
hyp1.lik  = [var(Ytrain)/2];    
hyp1.mean = mean(Ytrain);       

par = {@meanConst,@covSEard,'likGauss',[Xtrain(:,1:3),xtaux], Ytrain'};

hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
[m_1 s_1] = gp(hyp2, 'infExact', @meanConst, @covSEard, 'likGauss',  [Xtrain(:,1:3),xtaux], Ytrain', [Xtest(:,1:3),xpaux]);
            
else
    
    disp('Model not defined ...')
end
    
