function [Output] = MarmosetGP_CS7_v3(D,Output,gene2test);

%Scale the input data
Xtrain = Output.cleanX / Output.scalefactor;

gind = find(strcmp(D.textdata(2:end,1),gene2test)==1);
Ytrain = D.data(gind,Output.cleanindex);

pg1 = {@priorGauss,0,1};  
pg2 = {@priorGauss,1,1};  
pg3 = {@priorGauss,log(0.5),1};  
pc = {@priorClamped};

%Prior for base model
prior1.mean = {[]};  
prior1.cov  = {pg2;pg1};
prior1.lik  = {pg3};

%Prior for long lenghth-scale model
prior2.mean = {[]};  
prior2.cov  = {pc;pg1};
prior2.lik  = {pg3};

im1 = {@infPrior,@infExact,prior1};
im2 = {@infPrior,@infExact,prior2};

%Hyparameters: Am, EmDisc+Stalk, PGC, VE, SYS, ExMes, Tb    
HYP2 = cell(1,8);
HYP3 = cell(1,8);

%Initalise a high noiose std and long length scale
hyp1.cov  = [ (4) ; std(Ytrain) ];
hyp1.lik  = [std(Ytrain)/2];    
hyp1.mean = mean(Ytrain);       

hyp1B.cov  = [(40) ; var(Ytrain)];
hyp1B.lik  = [var(Ytrain)/2];    
hyp1B.mean = mean(Ytrain); 

%par = {@meanConst,@covSEiso,'likGauss',Xtrain(:,1:3), Ytrain'};
%par1 = {@meanConst,@covSEard,'likGauss',Xtrain(:,1:3), Ytrain',hyp1};
%par2 = {@meanConst,@covSEard,'likGauss',Xtrain(:,1:3), Ytrain',hyp3};
 
%First is all
par = {@meanConst,@covSEiso,'likGauss',Xtrain(:,1:3), Ytrain'};
hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
hyp3 = feval(@minimize, hyp1B, @gp, -1000, im2, par{:});         % optimise
HYP2{1,9} = hyp2;
HYP3{1,9} = hyp3;
%[m0 s0] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', OBJ1.vertices/400);
[L1_9 ] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(:,1:3), Ytrain');
[L2_9 ] = gp(hyp3, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(:,1:3), Ytrain');
    
%First is amnion
indT = find(Output.cleanAnotaton=="Am_CS7");
par = {@meanConst,@covSEiso,'likGauss',Xtrain(indT,1:3), Ytrain(indT)'};
hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
hyp3 = feval(@minimize, hyp1B, @gp, -1000, im2, par{:});         % optimise
HYP2{1,1} = hyp2;
HYP3{1,1} = hyp3;
%[m0 s0] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', OBJ1.vertices/400);
[L1_1 ] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');
[L2_1 ] = gp(hyp3, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');


%Second is EmDisk
indT = find(Output.cleanAnotaton=="EmDisc_CS7");
par = {@meanConst,@covSEiso,'likGauss',Xtrain(indT,1:3), Ytrain(indT)'};
hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
hyp3 = feval(@minimize, hyp1B, @gp, -1000, im2, par{:});         % optimise
HYP2{1,2} = hyp2;
HYP3{1,2} = hyp3;
%[m1 s1] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', OBJ1.vertices/400);
%[m3 s3] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', OBJ1.vertices/400);
[L1_2 ] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');
[L2_2 ] = gp(hyp3, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');

L1_3 = NaN;
L2_3 = NaN;
    
L1_4 = NaN;
L2_4 = NaN;

%Next we have SYS
indT = find(Output.cleanAnotaton=="SYS_CS7");
par = {@meanConst,@covSEiso,'likGauss',Xtrain(indT,1:3), Ytrain(indT)'};
hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
hyp3 = feval(@minimize, hyp1B, @gp, -1000, im2, par{:});         % optimise
HYP2{1,5} = hyp2;
HYP3{1,5} = hyp3;
%[m3 s3] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', OBJ1.vertices/400);
[L1_5 ] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');
[L2_5 ] = gp(hyp3, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');

%ExMes
indT = find(Output.cleanAnotaton=="ExMes_CS7");
par = {@meanConst,@covSEiso,'likGauss',Xtrain(indT,1:3), Ytrain(indT)'};
hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
hyp3 = feval(@minimize, hyp1B, @gp, -1000, im2, par{:});         % optimise
HYP2{1,6} = hyp2;
HYP3{1,6} = hyp3;
%[m5 s5] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', OBJ1.vertices/400);
[L1_6 ] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');
[L2_6 ] = gp(hyp3, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');

%ExMes
indT = find(Output.cleanAnotaton=="Stalk_CS7"  | Output.cleanAnotaton=="EmDisc_stalk_CS7");
par = {@meanConst,@covSEiso,'likGauss',Xtrain(indT,1:3), Ytrain(indT)'};
hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
hyp3 = feval(@minimize, hyp1B, @gp, -1000, im2, par{:});         % optimise
HYP2{1,8} = hyp2;
HYP3{1,8} = hyp3;
%[m5 s5] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', OBJ1.vertices/400);
[L1_8 ] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');
[L2_8 ] = gp(hyp3, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');

    
%ExMes
indT = find(Output.cleanAnotaton=="Tb_CS7");
par = {@meanConst,@covSEiso,'likGauss',Xtrain(indT,1:3), Ytrain(indT)'};
hyp2 = feval(@minimize, hyp1, @gp, -1000, im1, par{:});         % optimise
hyp3 = feval(@minimize, hyp1B, @gp, -1000, im2, par{:});         % optimise
HYP2{1,7} = hyp2;
HYP3{1,7} = hyp3;
%[m5 s5] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)', OBJ1.vertices/400);
[L1_7 ] = gp(hyp2, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');
[L2_7 ] = gp(hyp3, 'infExact', @meanConst, @covSEiso, 'likGauss',  Xtrain(indT,1:3), Ytrain(indT)');
 
%Output.m_0 = m0;
%Output.m_1 = m1;
%Output.m_2 = m2;
%Output.m_3 = m3;
%Output.m_5 = m5;

%Output.Ytrain = Ytrain;
%Output.Xtrain = Xtrain;
Output.L1 = [L1_1,L1_2,L1_3,L1_4,L1_5,L1_6,L1_7,L1_8,L1_9];
Output.L2 = [L2_1,L2_2,L2_3,L2_4,L2_5,L2_6,L2_7,L2_8,L2_9];
Output.HYP1 = HYP2;
Output.HYP2 = HYP3;

Output.im1 = im1;
Output.im2 = im2;
Output.im1 = prior1;
Output.im2 = prior2;

Output.Ytrain = Ytrain;

%Output.hyp2 = hyp1;
%Output.hyp2 = hyp2;
%Output.hyp3 = hyp3;

